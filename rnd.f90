! ============================================================================
! Adaptation of some stuff from Numerical Recipes in Fortran 90, chapter B7,
! to be at least somewhat OpenMP friendly. Original intentions were quite
! different from the end result. I wanted a module-global array of state
! vectors and ran1 to call omp_get_thread_num() to figure out which element to
! use. But accessing such array shared between all threads appears slow and
! performance scaling with the number of processes is strange (but maybe I am
! just missing something). The function omp_get_thread_num() is not all that
! fast either to be called at every rnd request.
!
! I can get a consistent performance only when I set a thread-local scalar
! rnd_state_loc at the beginning of a thread and use that as an argument to
! ran1. It is silly. It should be possible to hide all these state vectors
! from the end user somehow.
!
! BTW, the pointer gymnastics from NR is quite slow, but nothing of that
! survived here anyway.
!
! NB: gfortran needs -fno-strict-overflow, otherwise 'arith assump 3 fails'
! in rnd_init(); gfortran's -pedantic will discover range overflows
!
! Some sites to check (Mersenne Twister)
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!   http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html
!   http://www.cs.hmc.edu/~geoff/mtwist.html
! ============================================================================

module rnd
  use types_const
  implicit none
  private

  public :: rnd_init
  public :: ran1
  public :: gasdev
  public :: t_rnd_state

  integer, parameter :: k4b=selected_int_kind(9)
  integer(k4b), parameter :: hg=huge(1_k4b), hgm=-hg, hgng=hgm-1
  integer(k4b), parameter :: seq=0
  real(dp), save :: amm

  type t_rnd_state
     integer(k4b) ::  iran, jran, kran, nran, mran, ranv
     real(dp) :: g
     logical :: gaus_stored
  end type t_rnd_state

contains

  subroutine rnd_init(rnd_state)
    ! {{{ initialize an array of generator state vectors, one such vector
    !     for each thread
    implicit none
    type(t_rnd_state), dimension(:), intent(out) :: rnd_state
    integer(k4b) :: j, hgt, length
    integer(k4b), dimension(:,:), allocatable :: ranseeds

    hgt=hg
    if (hg /= 2147483647) call error_msg('rnd_init: arith assump 1 fails')
    if (hgng >= 0)        call error_msg('rnd_init: arith assump 2 fails')
    if (hgt+1 /= hgng)    call error_msg('rnd_init: arith assump 3 fails')
    if (not(hg) >= 0)     call error_msg('rnd_init: arith assump 4 fails')
    if (not(hgng) < 0)    call error_msg('rnd_init: arith assump 5 fails')
    if (hg+hgng >= 0)     call error_msg('rnd_init: arith assump 6 fails')
    if (not(-1_k4b) < 0)  call error_msg('rnd_init: arith assump 7 fails')
    if (not(0_k4b) >= 0)  call error_msg('rnd_init: arith assump 8 fails')
    if (not(1_k4b) >= 0)  call error_msg('rnd_init: arith assump 9 fails')

    length=size(rnd_state)
    allocate( ranseeds(length,5) )

    ! stupid Fortran cannot have nearest() in a parameter initialization
    ! expression
    amm=nearest(1.0_dp,-1.0_dp)/hgng
    if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
         call error_msg('rnd_init: arth assump 10 fails')
    
    ranseeds(:,1)=seq
    ranseeds(:,2:5)=spread(arth(1,1,size(ranseeds(:,1))),2,4)
    do j=1,4
       call ran_hash(ranseeds(:,j),ranseeds(:,j+1))
    end do
    where (ranseeds(:,1:3) < 0)  ranseeds(:,1:3)=not(ranseeds(:,1:3))
    where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1

    rnd_state%iran = ranseeds(:,1)
    rnd_state%jran = ranseeds(:,2)
    rnd_state%kran = ranseeds(:,3)
    rnd_state%mran = ranseeds(:,4)
    rnd_state%nran = ranseeds(:,5)
    rnd_state%ranv = rnd_state%nran

    rnd_state%gaus_stored=.false.

    deallocate(ranseeds)
    ! }}}
  end subroutine rnd_init

  subroutine ran_hash(il,ir)
    ! {{{ 
    implicit none
    integer(k4b), dimension(:), intent(inout) :: il,ir
    integer(k4b), dimension(size(il)) :: is
    integer(k4b) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
    ! }}}
  end subroutine ran_hash
  
  pure subroutine ran1(rnd_state,harvest)
    ! {{{ lagged Fibonacci generator combined with two Marsaglia shift
    !     sequences. On output, returns as harvest a uniform random deviate
    !     between 0.0 and 1.0 (exclusive of the endpoint values). The period
    !     of this generator is about 8.5e37. Validity of the integer model
    !     assumed by this generator is tested in rnd_init().
    !     [NR F90, chapter B7, page 1149]
    implicit none
    type(t_rnd_state), intent(inout) :: rnd_state
    real(dp), intent(out) :: harvest
    rnd_state%ranv=rnd_state%iran-rnd_state%kran
    if (rnd_state%ranv < 0) rnd_state%ranv=rnd_state%ranv+2147483579_k4b
    rnd_state%iran=rnd_state%jran
    rnd_state%jran=rnd_state%kran
    rnd_state%kran=rnd_state%ranv
    rnd_state%nran=ieor(rnd_state%nran,ishft(rnd_state%nran,13))
    rnd_state%nran=ieor(rnd_state%nran,ishft(rnd_state%nran,-17))
    rnd_state%nran=ieor(rnd_state%nran,ishft(rnd_state%nran,5))
    if (rnd_state%nran == 1) rnd_state%nran=270369_k4b
    rnd_state%mran=ieor(rnd_state%mran,ishft(rnd_state%mran,5))
    rnd_state%mran=ieor(rnd_state%mran,ishft(rnd_state%mran,-13))
    rnd_state%mran=ieor(rnd_state%mran,ishft(rnd_state%mran,6))
    rnd_state%ranv=ieor(rnd_state%nran,rnd_state%ranv)+rnd_state%mran
    harvest=amm*merge(rnd_state%ranv,not(rnd_state%ranv), rnd_state%ranv<0 )
    ! }}}
  end subroutine ran1

  subroutine gasdev(rnd_state,harvest)
    ! {{{ Returns in harvest a normally distributed deviate with zero mean
    !     and unit variance, using ran1 as the source of uniform deviates.
    type(t_rnd_state), intent(inout) :: rnd_state
    real(dp), intent(out) :: harvest
    real(dp) :: rsq,v1,v2
    if (rnd_state%gaus_stored) then
       harvest=rnd_state%g
       rnd_state%gaus_stored=.false.
    else
       do
          call ran1(rnd_state,v1)
          call ran1(rnd_state,v2)
          v1=2.0_dp*v1-1.0_dp
          v2=2.0_dp*v2-1.0_dp
          rsq=v1**2+v2**2
          if (rsq > 0.0_dp .and. rsq < 1.0_dp) exit
       end do
       rsq=sqrt(-2.0_dp*log(rsq)/rsq)
       harvest=v1*rsq
       rnd_state%g=v2*rsq
       rnd_state%gaus_stored=.true.
    end if
    ! }}}
  end subroutine gasdev

  subroutine error_msg(string)
    ! {{{ 
    character(len=*), intent(in) :: string
    write (*,*) string
    stop
    ! }}}
  end subroutine error_msg

  function arth(first,increment,n)
    ! {{{ returns arithmetic progression as an array
    integer(i4b), parameter :: NPAR_ARTH=16,NPAR2_ARTH=8
    integer(i4b), intent(IN) :: first,increment,n
    integer(i4b), dimension(n) :: arth
    integer(i4b) :: k,k2,temp
    if (n > 0) arth(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth(k)=arth(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth(k)=arth(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
    ! }}}
  end function arth
  
end module rnd

! Local variables:
! folded-file: t
! End:
