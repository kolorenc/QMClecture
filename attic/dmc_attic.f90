

  subroutine dmc_init_walkers(E_vmc)
    ! {{{ VMC sampling of trial wavefunction, we get energy as well
    real(dp), intent(out) :: E_vmc
    real(dp), dimension(dmc_dim) :: r
    type(t_walker) :: wlkr
    real(dp) :: weight, psi, psin, EL
    real(dp) :: E2
    integer :: i, error
    integer, parameter :: t_therm=1000  ! thermalization time
    NW=NWopt
    ! starting point
    call ran1(r)
    wlkr%r=r
    call psiT2(wlkr)
    psi=wlkr%psiTsq
    do i=1, t_therm     ! thermalization (to forget the starting point)
       call metropolis_step()
    end do
    E2=0.0_dp
    do i=NW, 1, -1   ! harvest
       call metropolis_step()
       walker(i)%r=r
       walker(i)%age=0
       call EL_drift(walker(i))
       E2=E2+real(walker(i)%EL,dp)
    end do
    E_vmc=real(E2/NW,dp)
  contains
    subroutine metropolis_step()
      ! {{{ one step of Metropolis algorithm
      real(dp), dimension(dmc_dim) :: vrnd, rn
      real(dp) :: rnd
      call gasdev(vrnd)
      rn=r+vrnd
      wlkr%r=rn
      call psiT2(wlkr)
      psin=wlkr%psiTsq
      weight=min(1.0_dp,psin/psi)
      call ran1(rnd)
      if ( rnd <= weight ) then
         r=rn
         psi=psin
      end if
      ! }}}
    end subroutine metropolis_step
    ! }}}
  end subroutine dmc_init_walkers

  subroutine add_walker(wlkr)
    ! {{{ adds walker 'wlkr' after the last element of array
    type(t_walker), intent(in) :: wlkr
    NW=NW+1
    if ( NW>NWmax ) then
       print *, "Too many walkers."
       stop
    end if
    walker(NW)=wlkr
    ! }}}
  end subroutine add_walker

  subroutine kill_walker(i)
    ! {{{ removes walker
    integer, intent(in) :: i
    if ( i==NW ) then
       NW=NW-1
    else
       walker(i)=walker(NW)
       NW=NW-1
    end if
    if ( NW < 1 ) then
       print *, "All walkers have gone ;-("
       stop
    end if
    ! }}}
  end subroutine kill_walker

  subroutine get_walker(i,wlkr)
    ! {{{ .
    integer, intent(in) :: i
    type(t_walker), intent(out) :: wlkr
    wlkr=walker(i)
    ! }}}
  end subroutine get_walker
  
  function Gd_ratio(wlkr_from,wlkr_to)
    ! {{{ ratio of Green's functions for imposing detailed balance
    real(dp) :: Gd_ratio
    type(t_walker), intent(in) :: wlkr_from, wlkr_to
    real(dp) :: forw, backw
    forw=sum((wlkr_to%r-wlkr_from%r-tau*wlkr_from%vD)**2)
    backw=sum((wlkr_from%r-wlkr_to%r-tau*wlkr_to%vD)**2)
    Gd_ratio=exp((forw-backw)/2.0_dp/tau)
    ! }}}
  end function Gd_ratio

  subroutine dmc_step(ET,energy,NW_out)
    ! {{{ time-step evolution of walker population, it seems to be (much)
    !  better (smaller fluctuations) to collect energy rather than ET that
    !  contains a correction for constant number of walkers
    real(dp), intent(in) :: ET
    real(dp), intent(out) :: energy
    integer, intent(out) :: NW_out
    real(dp) :: weight, rnd, balance, vDlsq, vDlsq_new, vDbarlsq, vDbarlsq_new
    real(dp) :: ET2
    real(dp), dimension(dmc_dim) :: vrnd, vD
    type(t_walker) :: wlkr
    integer :: num, i, j, start
    ET2=0.0_dp
    start=NW
    do i=start, 1, -1
       call gasdev_v(vrnd)
       ! Umrigar et al. trick
       !vDlsq=sum(walker(i)%vD**2)
       !vD=walker(i)%vD*(sqrt(1.0_dp+2.0_dp*vDlsq*tau)-1.0_dp)/vDlsq/tau
       !wlkr%r=walker(i)%r+sqrtau*vrnd+vD*tau
       wlkr%r=walker(i)%r+sqrtau*vrnd+walker(i)%vD*tau
       call EL_drift(wlkr)
       weight=exp(-tau/2.0_dp*(wlkr%EL+walker(i)%EL-2.0_dp*ET))
       !vDbarlsq=sum(vD**2)
       !vDlsq_new=sum(wlkr%vD**2)
       !vD=wlkr%vD*(sqrt(1.0_dp+2.0_dp*vDlsq_new*tau)-1.0_dp)/vDlsq_new/tau
       !vDbarlsq_new=sum(vD**2)
       !weight=exp(-tau/2.0_dp*((wlkr%EL-ET)*sqrt(vDbarlsq_new/vDlsq_new) &
       !     & +(walker(i)%EL-ET)*sqrt(vDbarlsq/vDlsq)))
       ! accept step only if walker didn't cross a node, otherwise leave
       ! the old positions
       if ( walker(i)%sgn == wlkr%sgn ) then
          ! additional detailed balance
          balance=min(1.0_dp,Gd_ratio(walker(i),wlkr) &
               & *wlkr%psiTsq/walker(i)%psiTsq)
          if ( walker(i)%age > age_limit ) then
             if ( walker(i)%age == age_limit+1 ) too_old=too_old+1
             balance=balance*exp(1.1_dp*(walker(i)%age-age_limit))
          end if
          call ran1(rnd)
          if ( rnd <= balance ) then
             walker(i)=wlkr
             db_accept=db_accept+1
          else
             walker(i)%age=walker(i)%age+1
             db_reject=db_reject+1
          end if
       end if
       ! birth/death
       call ran1(rnd)
!!$       num=min(floor(weight+rnd),2)
!!$       select case(num)
!!$       case(0)
!!$          call kill_walker(i)
!!$       case(1)
!!$          ET2=ET2+real(walker(i)%EL,dp)     ! sums to be done in double prec.
!!$       case(2)
!!$          call add_walker(walker(i))
!!$          ET2=ET2+2.0_dp*real(walker(i)%EL,dp)
!!$       end select
       num=floor(weight+rnd)
       select case(num)
       case(0)
          call kill_walker(i)
       case(1)
          ET2=ET2+real(walker(i)%EL,dp)     ! sums to be done in double prec.
       case default
          do j=2, num
             call add_walker(walker(i))
             ET2=ET2+2.0_dp*real(walker(i)%EL,dp)
          end do
       end select
    end do
    energy=real(ET2/NW,dp)
    NW_out=NW
    !ET=energy+(1.0_dp-NW*1.0_dp/NWopt)/tau
    !ET=energy-log(NW*1.0_dp/NWopt)
    ! }}}
  end subroutine dmc_step

  subroutine dmc_get_db_stats(db_acc,db_rej,old)
    ! {{{ get statistics about detailed balance accept/reject step, we use
    !  real numbers, since MPI does not support long integers (?) 
    real(dp) :: db_acc, db_rej
    integer, intent(out) :: old
    db_acc=real(db_accept,dp)
    db_rej=real(db_reject,dp)
    old=too_old
    ! }}}
  end subroutine dmc_get_db_stats
