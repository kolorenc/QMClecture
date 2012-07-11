! You need to supply 2 subroutines: EL_drift(wlkr) and psiT2(wlkr) in module
! dmc_input

! The derived type t_qmc_data contains allocatable arrays which is not
! possible in standard Fortran95. It is defined as TR15581 extension to F95
! (which is a part of F03).

#define VMC_DIFFUSION_DRIFT 1

module qmc
  use types_const, only: dp, i8b
  use rnd
  use qmc_input
  implicit none
  private
  
  public :: t_qmc_data
  public :: init_qmc_data, init_qmc
  public :: vmc_run

  real(dp), parameter :: q=2.0_dp    ! NWmax/NWopt
  integer, parameter :: age_limit=30 ! older walkers are "encouraged" to move

  type t_qmc_data
     real(dp) :: tau, sqrtau         ! time step

     ! size of the walker population
     integer :: NWopt                ! desired number of walkers
     integer :: NWmax                ! maximal number of walkers
     integer :: NW                   ! instant number of walkers

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

     ! trace of the total energy over the MC process (average over population)
     real(dp), dimension(:), allocatable :: Etot
  end type t_qmc_data

contains

  subroutine init_qmc_data(qmc_data,NW_opt,N_threads)
    ! {{{ initialization of the QMC process (common stuf for VMC and DMC),
    !     the storage is unnecessarily large for VMC but it should not matter
    !     much considering our goals
    type(t_qmc_data), intent(inout) :: qmc_data
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
    call randomize_walker_population(qmc_data)
    ! }}}
  end subroutine init_qmc_data

  subroutine randomize_walker_population(qmc_data)
    ! {{{ initial random population of walkers, this can be single thread
    type(t_qmc_data), intent(inout) :: qmc_data
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
       qmc_data%walker(iw)%age=0
       call EL_drift(qmc_data%walker(i))
    end do
    qmc_data%NW=NW
    qmc_data%rnd_state(1)=rnd_state
    ! }}}
  end subroutine randomize_walker_population

  subroutine init_qmc(qmc_data,t_step,N_steps)
    ! {{{ initialization of a VMC or DMC run
    type(t_qmc_data), intent(inout) :: qmc_data
    real(dp), intent(in) :: t_step
    integer, intent(in) :: N_steps
    integer :: error
    qmc_data%tau=t_step
    qmc_data%sqrtau=sqrt(t_step)
    qmc_data%accepted_moves=0
    qmc_data%total_moves=0
    qmc_data%walker%age=0
    qmc_data%too_old=0
    if ( allocated(qmc_data%Etot) ) deallocate(qmc_data%Etot)
    allocate( qmc_data%Etot(N_steps), stat=error )
    if ( error /= 0 ) then
       print *, "qmc_init: cannot allocate array of energies."
       stop
    end if
    qmc_data%Nsteps=N_steps
    ! }}}
  end subroutine init_qmc

  ! ==========================================================================
  ! variational Monte Carlo
  ! ==========================================================================

  function vmc_walker_move(vmc_data,rnd_state,iwlkr) result(accept)
    ! {{{ VMC step on a single walker, all electrons are moved at once;
    !     we want this to be thread safe
    type(t_qmc_data), intent(inout) :: vmc_data
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
    call psiT2(wlkr)
    r = wlkr%psiTsq / vmc_data%walker(iwlkr)%psiTsq

#else

    ! diffusion-drift move motivated by importance-sampled DMC
    ! one of the earliest references is [Rossky, Doll & Friedman,
    ! J. Chem. Phys. 69, 4628 (1978)]
    ! REMEMBER: the proposal probablity is assymetric and hence the ratio
    ! for accept-reject is more complicated than in pure Metropolis
    ! [Hastings, Biometrika 57, 97 (1970)]
    do i=1, sys_dim
       call gasdev(rnd_state,vrnd(i))
    end do
    vDi=vmc_data%walker(iwlkr)%vD
    ri = vmc_data%walker(iwlkr)%r
    rf = ri + vmc_data%sqrtau*vrnd + vDi*vmc_data%tau
    wlkr%r = rf
    call EL_drift(wlkr)
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
       call EL_drift(wlkr)   ! for diffusion-drift this is already done
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

  function vmc_step_thread(vmc_data,threadID,iwMin,iwMax) result(accept)
    ! {{{ VMC moves on a chunk of the total population
    type(t_qmc_data), intent(inout) :: vmc_data
    integer, intent(in) :: threadID, iwMin, iwMax
    integer :: accept
    type(rnd_state_vector) :: rnd_state
    integer :: iw
    rnd_state=vmc_data%rnd_state(threadID)
    accept=0
    do iw=iwMin, iwMax
       accept=accept+VMC_walker_move(vmc_data,rnd_state,iw)
    end do
    vmc_data%rnd_state(threadID)=rnd_state
    ! }}}
  end function vmc_step_thread

  subroutine vmc_run(vmc_data)
    ! {{{ one move for each walker in the population
    type(t_qmc_data), intent(inout) :: vmc_data
    integer, dimension(:), allocatable :: limits, accept
    real(dp), dimension(:), allocatable :: Etot
    integer :: j, k, chunk, Nthreads, NW, Nsteps, istep, error

    Nthreads=vmc_data%Nthreads
    NW=vmc_data%NW
    Nsteps=vmc_data%Nsteps

    allocate( limits(Nthreads+1), accept(Nthreads), Etot(Nthreads), &
         stat=error )
    if ( error /= 0 ) then
       print *, "vmc_run: cannot allocate arrays."
       stop
    end if
    chunk = int(NW/Nthreads)
    do k=2, Nthreads
       limits(k)=(k-1)*chunk
    end do
    limits(1)=0
    limits(Nthreads+1)=NW

    do istep=1, Nsteps

       Etot=0.0_dp
       !schedule(dynamic) schedule(dynamic,10) schedule(static)
       !$omp parallel do default(shared) private(k) schedule(guided,1)
       do k=1, Nthreads
          accept(k)=VMC_step_thread(vmc_data,k,limits(k)+1,limits(k+1))
          do j=limits(k)+1, limits(k+1)
             Etot(k)=Etot(k)+vmc_data%walker(j)%EL
          end do
       end do
       !$omp end parallel do

       vmc_data%total_moves    = vmc_data%total_moves + NW
       vmc_data%accepted_moves = vmc_data%accepted_moves + sum(accept)
       vmc_data%Etot(istep)=sum(Etot)/NW

    end do
       
    deallocate(limits,accept,Etot)
    ! }}}
  end subroutine vmc_run


  ! ==========================================================================
  ! diffusion Monte Carlo
  ! ==========================================================================

end module qmc


! Local variables:
! folded-file: t
! End:
