program vmc_optimize
  use nrtype
  use nr_random
  use dmc_input
  implicit none

  include 'mpif.h'

  integer, parameter :: NW_tot=1000000
  integer, parameter :: t_therm=1000
  real(sp) :: EL, E_vmc, E_vmc_loc, psi, psin, weight
  real(dp) :: EL2
  real(sp), dimension(dmc_dim) :: r
  type(t_walker) :: wlkr
  integer :: i, NW, k, kmax
  real(sp), dimension(:), allocatable :: rnd_init
  integer, dimension(:), allocatable :: rnd_seq
  integer :: rnd_seq_local

  integer :: ierr, rank, comm_size

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_size,ierr)

  if ( rank==0 ) then
     allocate(rnd_init(0:comm_size-1), rnd_seq(0:comm_size-1), stat=ierr)
     if ( ierr /=0 ) then
        print *, "cannot allocate arrays."
        stop
     end if
     call ran1(rnd_init)
     rnd_seq=int(1.0e6_sp*rnd_init)
     NW=NW_tot/comm_size
     open(unit=101,file="optimize.dat",action="write")
  end if
  
  call MPI_SCATTER(rnd_seq,1,MPI_INTEGER,rnd_seq_local,1, &
       & MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call ran_seed(sequence=rnd_seq_local)
  
  call MPI_BCAST(NW,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  kmax=20
  do k=0, kmax
     Z1=1.0_sp+k*1.0_sp/kmax

     call ran1(r)
     r=r-(/ .5_sp, .5_sp, .5_sp /)
     wlkr%r=r
     call psiT2_EL(wlkr)
     psi=wlkr%psiTsq
     do i=1, t_therm     ! thermalization (to forget the starting point)
        call metropolis_step()
     end do
     EL2=0.0_dp
     do i=NW, 1, -1   ! harvest
        call metropolis_step()
        EL2=EL2+EL
     end do
     E_vmc_loc=real(EL2,sp)
     
     call MPI_REDUCE(E_vmc_loc,E_vmc,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     
     if ( rank==0 ) then
        E_vmc=E_vmc/NW/comm_size
        print *, Z1, E_vmc
        write(unit=101,fmt=*) Z1, E_vmc
     end if
     
  end do

  if ( rank==0 ) close(unit=101)

  call MPI_FINALIZE(ierr)
  
contains

  subroutine metropolis_step()
    ! {{{ one step of Metropolis algorithm
    real(sp), dimension(dmc_dim) :: vrnd, rn
    real(sp) :: rnd
    call gasdev(vrnd)
    rn=r+vrnd
    wlkr%r=rn
    call psiT2_EL(wlkr)
    psin=wlkr%psiTsq
    weight=min(1.0_sp,psin/psi)
    call ran1(rnd)
    if ( rnd <= weight ) then
       r=rn
       psi=psin
       EL=wlkr%EL
    end if
    ! }}}
  end subroutine metropolis_step
  

end program vmc_optimize

! Local variables:
! folded-file: t
