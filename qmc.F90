! You need to supply 2 subroutines:
! EL_drift(sys,wlkr) and psiT2(sys,wlkr) in module dmc_input

! The derived type t_qmc_data contains allocatable arrays which is not
! possible in standard Fortran95. It is defined as TR15581 extension to F95
! (which is a part of F03).

#define VMC_DIFFUSION_DRIFT 1
#define ET_RUNNING_AVERAGE 1

module qmc
  use types_const, only: dp, i8b
  use rnd
  use qmc_input
  implicit none
  private
  
  public :: t_qmc_data
  public :: init_qmc_data, init_qmc
  public :: vmc_run, dmc_run

  real(dp), parameter :: q=2.0_dp    ! NWmax/NWopt
  integer, parameter :: age_limit=30 ! older walkers are "encouraged" to move

  type t_qmc_data
     real(dp) :: tau, sqrtau         ! time step

     ! size of the walker population
     integer :: NWopt                ! desired number of walkers
     integer :: NWmax                ! maximal number of walkers
     integer :: NW                   ! instant number of walkers

     real(dp) :: ET                  ! DMC trial energy

     integer :: Nsteps               ! number of MC steps in a run

     ! number of single walker moves can get really large, for monitoring
     ! the detailed balance accept/reject we thus need long (8-byte) integers
     integer(i8b) :: accepted_moves
     integer(i8b) :: total_moves     ! for Metropolis step monitoring
     integer :: too_old              ! number trapped walkers (DMC)

     ! walker population
     type(t_walker), dimension(:), allocatable :: walker

     ! state of the random number generator for each thread
     integer :: Nthreads
     type(rnd_state_vector), dimension(:), allocatable :: rnd_state

     ! total energy (average over population, full history over steps)
     real(dp), dimension(:), allocatable :: Etot

     ! total energy squared (average over population, full history over steps);
     ! needed to get the variance
     real(dp), dimension(:), allocatable :: Etot2

     ! growth estimator of the total energy
     real(dp), dimension(:), allocatable :: EtotG

     ! trace of the population size
     integer, dimension(:), allocatable :: NWtr
  end type t_qmc_data

contains

  subroutine init_qmc_data(qmc_data,sys,NW_opt,N_threads)
    ! {{{ initialization of the QMC process (common stuf for VMC and DMC),
    !     the storage is unnecessarily large for VMC but it should not matter
    !     much considering our goals
    type(t_qmc_data), intent(inout) :: qmc_data
    type(t_sys), intent(in) :: sys
    integer, intent(in) :: NW_opt, N_threads
    integer :: error

    ! allocate walker storage, make extra space for DMC
    qmc_data%NWopt=NW_opt
    qmc_data%NWmax=int(q*NW_opt)
    if ( allocated(qmc_data%walker) ) deallocate(qmc_data%walker)
    allocate( qmc_data%walker(qmc_data%NWmax), stat=error )
    if ( error /= 0 ) then
       print *, "qmc_init: cannot allocate array of walkers."
       stop
    end if

    ! initialize random numbers
    qmc_data%Nthreads=N_threads
    if ( allocated(qmc_data%rnd_state) ) deallocate(qmc_data%rnd_state)
    allocate( qmc_data%rnd_state(N_threads), stat=error )
    if ( error /= 0 ) then
       print *, "qmc_init: cannot allocate storage for rnd state."
       stop
    end if
    call rnd_init(qmc_data%rnd_state)

    ! initialize the walker population
    call randomize_walker_population(qmc_data,sys)
    ! }}}
  end subroutine init_qmc_data

  subroutine randomize_walker_population(qmc_data,sys)
    ! {{{ initial random population of walkers, this can be single thread
    type(t_qmc_data), intent(inout) :: qmc_data
    type(t_sys), intent(in) :: sys
    type(rnd_state_vector) :: rnd_state
    real(dp) :: x
    integer :: i, iw, NW
    rnd_state=qmc_data%rnd_state(1)
    NW=qmc_data%NWopt
    do iw=1, NW
       do i=1, sys_dim
          call ran1(rnd_state,x)
          qmc_data%walker(iw)%r(i)=x-0.5_dp
       end do
       call EL_drift(sys,qmc_data%walker(i))
    end do
    qmc_data%NW=NW
    qmc_data%ET=sum(qmc_data%walker%EL)/qmc_data%NWopt
    qmc_data%rnd_state(1)=rnd_state
    ! }}}
  end subroutine randomize_walker_population

  subroutine init_qmc(qmc_data,t_step,N_steps)
    ! {{{ initialization of a VMC or DMC run
    type(t_qmc_data), intent(inout) :: qmc_data
    real(dp), intent(in) :: t_step
    integer, intent(in) :: N_steps
    integer :: error

    ! parameters of the walk
    qmc_data%tau=t_step
    qmc_data%sqrtau=sqrt(t_step)
    qmc_data%Nsteps=N_steps

    ! reset counters
    qmc_data%accepted_moves=0
    qmc_data%total_moves=0
    qmc_data%too_old=0

    ! reset walkers (but we want to keep positions etc.)
    qmc_data%walker%age=0
    qmc_data%walker%wt=1.0_dp

    ! arrays for measuring observables
    if ( allocated(qmc_data%Etot) ) deallocate(qmc_data%Etot)
    if ( allocated(qmc_data%Etot2) ) deallocate(qmc_data%Etot2)
    if ( allocated(qmc_data%EtotG) ) deallocate(qmc_data%EtotG)
    if ( allocated(qmc_data%NWtr) ) deallocate(qmc_data%NWtr)
    allocate( qmc_data%Etot(N_steps), qmc_data%Etot2(N_steps), stat=error )
    if ( error /= 0 ) then
       print *, "qmc_init: cannot allocate arrays of energies."
       stop
    end if
    ! }}}
  end subroutine init_qmc


  ! ==========================================================================
  ! variational Monte Carlo
  ! ==========================================================================

  function vmc_walker_move(vmc_data,sys,rnd_state,iwlkr) result(accept)
    ! {{{ VMC step on a single walker, all electrons are moved at once;
    !     we want this to be thread safe
    type(t_qmc_data), intent(inout) :: vmc_data
    type(t_sys), intent(in) :: sys
    type(rnd_state_vector), intent(inout) :: rnd_state
    integer, intent(in) :: iwlkr
    integer :: accept
    real(dp), dimension(sys_dim) :: vrnd, vDi, vDf, ri, rf
    type(t_walker) :: wlkr
    real(dp) :: r, q, y
    integer :: i

#ifndef VMC_DIFFUSION_DRIFT

    ! Metropolis step [MRRT&T, J. Chem. Phys. 21, 1087 (1953)];
    ! not the mot efficient sampling with respect to the number of steps but
    ! each of the step is fast;
    ! the first quantum Monte Carlo of this type was
    ! [McMillan, Phys. Rev. 138, A442 (1965)]
    do i=1, sys_dim
       call ran1(rnd_state,vrnd(i))
    end do
    wlkr%r = vmc_data%walker(iwlkr)%r + vmc_data%tau*(vrnd-0.5_dp)
    call psiT2(sys,wlkr)
    r = wlkr%psiTsq / vmc_data%walker(iwlkr)%psiTsq

#else

    ! diffusion-drift move motivated by importance-sampled DMC,
    ! one of the earliest references is [Rossky, Doll & Friedman,
    ! J. Chem. Phys. 69, 4628 (1978)] but similar "smart" sampling is
    ! described also in [Ceperley, Chester & Kalos, Phys. Rev. B 16,
    ! 3081 (1977)]
    ! REMEMBER: the proposal probablity is assymetric and hence the ratio
    ! for accept-reject is more complicated than in the pure Metropolis
    ! [Hastings, Biometrika 57, 97 (1970)]
    do i=1, sys_dim
       call gasdev(rnd_state,vrnd(i))
    end do
    vDi=vmc_data%walker(iwlkr)%vD
    ri = vmc_data%walker(iwlkr)%r
    rf = ri + vmc_data%sqrtau*vrnd + vDi*vmc_data%tau
    wlkr%r = rf
    call EL_drift(sys,wlkr)
    vDf=wlkr%vD
    q = 0.5_dp*vmc_data%tau*(sum(vDi*vDi)-sum(vDf*vDf)) &
         + sum( (ri-rf)*(vDi+vDf) )
    r = wlkr%psiTsq / vmc_data%walker(iwlkr)%psiTsq * exp(q)

#endif

    ! accept-reject
    accept=0
    call ran1(rnd_state,y)   ! y is always less than 1
    if ( y < r ) then
#ifndef VMC_DIFFUSION_DRIFT
       call EL_drift(sys,wlkr)   ! for diffusion-drift this is already done
#endif
       vmc_data%walker(iwlkr)=wlkr
       vmc_data%walker(iwlkr)%age=0
       accept=1
    else
       vmc_data%walker(iwlkr)%age = vmc_data%walker(iwlkr)%age + 1
       accept=0
    end if
    ! }}}
  end function vmc_walker_move

  function vmc_step_thread(vmc_data,sys,threadID,iwMin,iwMax) result(accept)
    ! {{{ VMC moves on a chunk of the total population
    type(t_qmc_data), intent(inout) :: vmc_data
    type(t_sys), intent(in) :: sys
    integer, intent(in) :: threadID, iwMin, iwMax
    integer :: accept
    type(rnd_state_vector) :: rnd_state
    integer :: iw
    rnd_state=vmc_data%rnd_state(threadID)
    accept=0
    do iw=iwMin, iwMax
       accept=accept+vmc_walker_move(vmc_data,sys,rnd_state,iw)
    end do
    vmc_data%rnd_state(threadID)=rnd_state
    ! }}}
  end function vmc_step_thread

  subroutine vmc_run(vmc_data,sys)
    ! {{{ the whole VMC run
    type(t_qmc_data), intent(inout) :: vmc_data
    type(t_sys), intent(in) :: sys
    integer, dimension(:), allocatable :: limits, accept
    real(dp), dimension(:), allocatable :: Etot, Etot2
    integer :: j, k, chunk, Nthreads, NW, Nsteps, istep, error

    Nthreads=vmc_data%Nthreads
    NW=vmc_data%NW
    Nsteps=vmc_data%Nsteps

    allocate( limits(Nthreads+1), accept(Nthreads), Etot(Nthreads), &
          Etot2(Nthreads), stat=error )
    if ( error /= 0 ) then
       print *, "vmc_run: cannot allocate arrays."
       stop
    end if

    ! equal workload for each thread
    chunk = int(NW/Nthreads)
    do k=2, Nthreads
       limits(k)=(k-1)*chunk
    end do
    limits(1)=0
    limits(Nthreads+1)=NW

    ! one VMC step is one move of each walker in the population
    do istep=1, Nsteps

       Etot=0.0_dp
       Etot2=0.0_dp
       !schedule(dynamic) schedule(dynamic,10) schedule(static)
       !$omp parallel do default(shared) private(k) schedule(guided,1)
       do k=1, Nthreads
          accept(k)=vmc_step_thread(vmc_data,sys,k,limits(k)+1,limits(k+1))
          do j=limits(k)+1, limits(k+1)
             Etot(k)=Etot(k)+vmc_data%walker(j)%EL
             Etot2(k)=Etot2(k)+vmc_data%walker(j)%EL**2
          end do
       end do
       !$omp end parallel do

       vmc_data%total_moves    = vmc_data%total_moves + NW
       vmc_data%accepted_moves = vmc_data%accepted_moves + sum(accept)
       vmc_data%Etot(istep)=sum(Etot)/NW
       vmc_data%Etot2(istep)=sum(Etot2)/NW

    end do
       
    deallocate(limits,accept,Etot,Etot2)
    ! }}}
  end subroutine vmc_run


  ! ==========================================================================
  ! diffusion Monte Carlo
  ! ==========================================================================

  function dmc_walker_move(dmc_data,sys,rnd_state,iwlkr) result(accept)
    ! {{{ (pure) DMC step on a single walker using an all-electrons move;
    !     NB: this is the same as VMC diffusion-drift move, extra is only
    !     the weight and the fixed-node constraint
    type(t_qmc_data), intent(inout) :: dmc_data
    type(t_sys), intent(in) :: sys
    type(rnd_state_vector), intent(inout) :: rnd_state
    integer, intent(in) :: iwlkr
    integer :: accept
    real(dp), dimension(sys_dim) :: vrnd, vDi, vDf, ri, rf
    real(dp) :: ELi, r, q, y
    type(t_walker) :: wlkr
    integer :: i, age
    logical :: acc

    ! diffusion-drift
    do i=1, sys_dim
       call gasdev(rnd_state,vrnd(i))
    end do
    vDi=dmc_data%walker(iwlkr)%vD
    ri = dmc_data%walker(iwlkr)%r
    rf = ri + dmc_data%sqrtau*vrnd + vDi*dmc_data%tau
    wlkr%r = rf
    call EL_drift(sys,wlkr)

    ! calculation of weight
    ELi=dmc_data%walker(iwlkr)%EL
    wlkr%wt = dmc_data%walker(iwlkr)%wt * &
         exp(-dmc_data%tau*0.5_dp*(wlkr%EL+ELi-2.0_dp*dmc_data%ET))

    ! detailed balance
    acc=.false.
    if ( dmc_data%walker(iwlkr)%sgn == wlkr%sgn ) then
       ! additional detailed balance
       vDf=wlkr%vD
       q = 0.5_dp*dmc_data%tau*(sum(vDi*vDi)-sum(vDf*vDf)) &
            + sum( (ri-rf)*(vDi+vDf) )
       r = wlkr%psiTsq / dmc_data%walker(iwlkr)%psiTsq * exp(q)
       ! if the walker is stuck for a long time, encourage it to move; this
       ! hack should be needed only for poor trial functions
       age=dmc_data%walker(iwlkr)%age
       if ( age > age_limit ) then
          if ( age == age_limit+1 ) then
             dmc_data%too_old = dmc_data%too_old + 1
          end if
          r=r*exp(1.1_dp*(age-age_limit))
       end if
       call ran1(rnd_state,y)
       if ( y < r ) then
          dmc_data%walker(iwlkr)=wlkr
          acc=.true.
       else
          acc=.false.
       end if
    end if

    ! this extra condition is needed to correctly increment age (but is this
    ! really the correct implementation?)
    if ( acc ) then
       accept=1
    else
       accept=0
       dmc_data%walker(iwlkr)%age = dmc_data%walker(iwlkr)%age + 1
    end if
    ! }}}
  end function dmc_walker_move

  function dmc_step_thread(dmc_data,sys,threadID,iwMin,iwMax) result(accept)
    ! {{{ DMC moves on a chunk of the total population
    type(t_qmc_data), intent(inout) :: dmc_data
    type(t_sys), intent(in) :: sys
    integer, intent(in) :: threadID, iwMin, iwMax
    integer :: accept
    type(rnd_state_vector) :: rnd_state
    integer :: iw
    rnd_state=dmc_data%rnd_state(threadID)
    accept=0
    do iw=iwMin, iwMax
       accept=accept+dmc_walker_move(dmc_data,sys,rnd_state,iw)
    end do
    dmc_data%rnd_state(threadID)=rnd_state
    ! }}}
  end function dmc_step_thread

  subroutine dmc_run(dmc_data,sys)
    ! {{{ the whole DMC run
    type(t_qmc_data), intent(inout) :: dmc_data
    type(t_sys), intent(in) :: sys
    type(rnd_state_vector) :: rnd_state
    integer, dimension(:), allocatable :: limits, accept
    integer :: j, k, chunk, Nthreads, NW, Nsteps, istep, error, kstart
    integer :: num
    real(dp) :: Etot, Etot2, y, ET

    Nthreads=dmc_data%Nthreads
    NW=dmc_data%NW
    Nsteps=dmc_data%Nsteps

    allocate( limits(Nthreads+1), accept(Nthreads), &
         dmc_data%NWtr(Nsteps), stat=error )
    !allocate( dmc_data%EtotG(Nsteps) )
    if ( error /= 0 ) then
       print *, "dmc_run: cannot allocate arrays."
       stop
    end if

    ! initial trial energy
    dmc_data%ET=sum(dmc_data%walker%EL)/NW
    ET=0.0_dp

    ! one DMC step is one move of each walker in the population
    do istep=1, Nsteps

       ! equal workload for each thread; chunk changes during the simulation
       ! because population size changes
       chunk = int(NW/Nthreads)
       do k=2, Nthreads
          limits(k)=(k-1)*chunk
       end do
       limits(1)=0
       limits(Nthreads+1)=NW

       ! diffusion-drift moves are done in parallel
       !schedule(dynamic) schedule(dynamic,10) schedule(static)
       !$omp parallel do default(shared) private(k) schedule(guided,1)
       do k=1, Nthreads
          accept(k)=dmc_step_thread(dmc_data,sys,k,limits(k)+1,limits(k+1))
       end do
       !$omp end parallel do

       dmc_data%total_moves    = dmc_data%total_moves + NW
       dmc_data%accepted_moves = dmc_data%accepted_moves + sum(accept)

       ! birth-death moves are done in a single thread, otherwise we would
       ! need an extra walker array for each CPU core
       rnd_state=dmc_data%rnd_state(1)
       kstart=NW
       Etot=0.0_dp
       Etot2=0.0_dp
       do k=kstart, 1, -1
          call ran1(rnd_state,y)
          num=floor(dmc_data%walker(k)%wt+y)
          dmc_data%walker(k)%wt=1.0_dp  ! reset the weight
          select case(num)
          case(0)
             ! kill this walker
             if ( k==NW ) then
                NW=NW-1
             else
                dmc_data%walker(k)=dmc_data%walker(NW)
                NW=NW-1
             end if
             if ( NW < 1 ) then
                print *, "All walkers are gone ;-("
                stop
             end if
          case(1)
             ! one walker continues, do nothing here
          case default
             ! walker will be spawned
              do j=2, num
                ! add a copy (or copies) of the walker
                NW=NW+1
                if ( NW > dmc_data%NWmax ) then
                   print *, "Too many walkers."
                   stop
                end if
                dmc_data%walker(NW)=dmc_data%walker(k)
             end do
          end select
          Etot = Etot + num*dmc_data%walker(k)%EL
          Etot2= Etot2+ num*dmc_data%walker(k)%EL**2
       end do
       dmc_data%rnd_state(1)=rnd_state

       Etot=Etot/NW
       dmc_data%Etot(istep)=Etot
       dmc_data%Etot2(istep)=Etot2/NW

#ifdef ET_RUNNING_AVERAGE
       ! trial energy as the best estimate we have so far (this looks like it
       ! can be biased by initial projection)
       ET=ET+Etot          !sum(dmc_data%Etot(1:istep))/istep
       dmc_data%ET=ET/istep-log(NW*1.0_dp/dmc_data%NWopt)
#else
       ! trial energy as an average over the population at each time step (this
       ! should fluctuate more but these fluctuations are possibly overshadowed
       ! by the population-control term enyway)
       dmc_data%ET=Etot-log(NW*1.0_dp/dmc_data%NWopt)
#endif

       !dmc_data%EtotG(istep)=dmc_data%ET
       dmc_data%NWtr(istep)=NW

    end do ! istep

    deallocate(limits,accept)
    ! }}}
  end subroutine dmc_run

end module qmc


! Local variables:
! folded-file: t
! End:
