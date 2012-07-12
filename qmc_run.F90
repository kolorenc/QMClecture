! Histogram in gnuplot:
! bin(x,width)=width*floor(x/width) + width/2.0
! plot 'Evmc_trace.dat' using (bin($1,.001)):(1.0) smooth freq with boxes

program qmc_run
#ifdef _OPENMP
  use omp_lib, only: omp_get_max_threads
#endif
  use types_const, only: dp
  use qmc
  use reblocking
  implicit none

  type(t_qmc_data) :: qmc_data
  integer, parameter :: NW=2**11
  integer, parameter :: Nsteps=2**16
  integer :: Nthreads
  !real(dp), parameter :: timestep=0.8_dp ! for pure Metropolis VMC
  real(dp), parameter :: VMCtstep=0.1_dp
  real(dp), parameter :: DMCtstep=0.01_dp

#ifdef _OPENMP
  ! the return value of omp_get_max_threads() is controlled by
  !   export OMP_NUM_THREADS=n
  Nthreads=omp_get_max_threads()
#else
  Nthreads=1
#endif
  print *, "number of threads:", Nthreads
  call init_qmc_data(qmc_data,NW,Nthreads)

  ! thermalization
  call init_qmc(qmc_data,VMCtstep,2048)
  call vmc_run(qmc_data)

  ! actual data harvest
  call init_qmc(qmc_data,VMCtstep,Nsteps)
  call vmc_run(qmc_data)

  call process_results("Evmc_trace.dat","Evmc_reblock.dat","E2vmc_reblock.dat")

  ! thermalization
  call init_qmc(qmc_data,DMCtstep,2048)
  call dmc_run(qmc_data)

  ! actual data harvest
  call init_qmc(qmc_data,DMCtstep,Nsteps)
  call dmc_run(qmc_data)

  call process_results("Edmc_trace.dat","Edmc_reblock.dat","E2dmc_reblock.dat")

contains

  subroutine process_results(file_Etrace,file_Ereblock,file_E2reblock)
    ! {{{ calculate errorbars and statistics of the MC process
    character(len=*), intent(in) :: file_Etrace, file_Ereblock, file_E2reblock
    real(dp) :: En, errEn, En2, errEn2, Var, errVar, corrlen, ageAV
    integer :: i, ageM, finalNblock

    open(unit=21,file=file_Etrace)
    write(unit=21,fmt="(e14.7)") qmc_data%Etot
    close(unit=21)

    if ( allocated(qmc_data%NWtr) ) then
       open(unit=21,file="NWdmc_trace.dat")
       write(unit=21,fmt="(i8)") qmc_data%NWtr
       close(unit=21)
    end if

    print *
    print *, "acceptance ratio:", &
         qmc_data%accepted_moves / (1.0_dp*qmc_data%total_moves)
    ageM=qmc_data%walker(1)%age
    ageAV=ageM
    do i=2, NW
       if ( qmc_data%walker(i)%age > ageM ) ageM=qmc_data%walker(i)%age
       ageAV=ageAV+qmc_data%walker(i)%age
    end do
    ageAV=ageAV/NW
    print *, "maximal age:", ageM
    print *, "average age:", ageAV

    En=sum(qmc_data%Etot)/Nsteps
    errEn=errorbar(qmc_data%Etot,En,file_Ereblock,corrlen,finalNBlock)
    print *, "---"
    print *, "total energy:", En, "+/-", errEn
    print *, "correlation length for E (steps):", corrlen
    print *, "final number of block in the statistics:", finalNblock

    En2=sum(qmc_data%Etot2)/Nsteps
    errEn2=errorbar(qmc_data%Etot2,En2,file_E2reblock,corrlen,finalNBlock)
    Var=En2-En**2
    errVar=sqrt(errEn2**2+(2*En*errEn)**2)
    print *, "---"
    print *, "variance:", Var, "+/-", errVar
    print *, "correlation length for E^2 (steps):", corrlen
    print *, "final number of block in the statistics:", finalNblock

    if ( allocated(qmc_data%EtotG) ) then
       En=sum(qmc_data%EtotG)/Nsteps
       errEn=errorbar(qmc_data%EtotG,En,file_Ereblock,corrlen,finalNBlock)
       print *, "---"
       print *, "total energy (growth):", En, "+/-", errEn
       print *, "correlation length:", corrlen
       print *, "final number of block in the statistics:", finalNblock
    end if
    ! }}}
  end subroutine process_results

end program qmc_run


! Local variables:
! folded-file: t
! End:
