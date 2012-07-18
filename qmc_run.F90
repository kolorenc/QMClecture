! Histogram in gnuplot:
! bin(x,width)=width*floor(x/width) + width/2.0
! plot 'Evmc_trace.dat' using (bin($1,.001)):(1.0) smooth freq with boxes

!f(a,b,c,x)=a*x**2+b*x+c
!fit f(a2,b2,c2,x) 'ET_comparison.dat' using 1:4:5 via a2,b2,c2
!plot [0:0.011]'ET_comparison.dat' using 1:2:3 w errorbars, f(a1,b1,c1,x) lt 1, 'ET_comparison.dat' using 1:4:5 w errorbars, -2.903724

program qmc_run
#ifdef _OPENMP
  use omp_lib, only: omp_get_max_threads
#endif
  use types_const, only: dp
  use qmc_input
  use qmc
  use reblocking
  implicit none

  type(t_sys) :: sys
  type(t_qmc_data) :: qmc_data
  integer, parameter :: NW=2**11
  integer, parameter :: Nsteps=2**16
  integer :: Nthreads
  !real(dp), parameter :: timestep=0.8_dp ! for pure Metropolis VMC
  real(dp), parameter :: VMCtstep=0.1_dp
  real(dp), parameter :: DMCtstep=0.01_dp
  integer, parameter :: extra_unit=22

  integer :: i

  ! watch out for the DMC acceptance ratio at larger Znuc - smaller timestep is
  ! needed in the high-density region (or smarter moves such as Umrigar,
  ! Nightingale & Runge
  call init_sys(sys,2.0_dp)

#ifdef _OPENMP
  ! the return value of omp_get_max_threads() is controlled by
  !   export OMP_NUM_THREADS=n
  Nthreads=omp_get_max_threads()
#else
  Nthreads=1
#endif
  print *, "number of threads:", Nthreads
  call init_qmc_data(qmc_data,sys,NW,Nthreads)

  ! thermalization
  call init_qmc(qmc_data,VMCtstep,2048)
  call vmc_run(qmc_data,sys)

  ! actual data harvest
  call init_qmc(qmc_data,VMCtstep,Nsteps)
  call vmc_run(qmc_data,sys)
  call show_results(qmc_data,"Evmc_reblock.dat","E2vmc_reblock.dat", &
       "Evmc_trace.dat")

  ! thermalization
  call init_qmc(qmc_data,DMCtstep,2048)
  call dmc_run(qmc_data,sys)

  open(unit=extra_unit,file="tstep_ET_instant.dat")
  write(unit=extra_unit,fmt=*) "# time step   energy   errorbar  acc"
  do i=0, 9
     ! actual data harvest
     call init_qmc(qmc_data,DMCtstep*(10.0_dp-i)/10.0_dp, &
          int(Nsteps*10.0_dp/(10-i)))
     call dmc_run(qmc_data,sys)
     call show_results(qmc_data,"Edmc_reblock.dat","E2dmc_reblock.dat", &
          "Edmc_trace.dat","NWdmc_trace.dat",extra_unit)
  end do
  close(unit=extra_unit)

  ! actual data harvest
!!$  call init_qmc(qmc_data,DMCtstep/2,Nsteps*2)
!!$  call dmc_run(qmc_data,sys)
!!$  call process_results(qmc_data,"Edmc_trace.dat","Edmc_reblock.dat", &
!!$       "E2dmc_reblock.dat")

  ! actual data harvest
!!$  call init_qmc(qmc_data,DMCtstep/4,Nsteps*4)
!!$  call dmc_run(qmc_data,sys)
!!$  call process_results(qmc_data,"Edmc_trace.dat","Edmc_reblock.dat", &
!!$       "E2dmc_reblock.dat")

contains

  subroutine show_results(qmc_data,file_Ereblock,file_E2reblock, &
       file_Etrace,file_NWtrace,extra_unit)
    ! {{{ calculate and write measured quantities, their errorbars and
    !     statistics of the MC process
    !     argument qmc_data is inout due to the implementation of errorbar()
    type(t_qmc_data), intent(inout) :: qmc_data
    character(len=*), intent(in) :: file_Ereblock, file_E2reblock
    character(len=*), optional, intent(in) :: file_Etrace, file_NWtrace
    integer, optional, intent(in) :: extra_unit
    real(dp) :: En, errEn, En2, errEn2, Var, errVar, corrlen, ageAV
    integer :: i, ageM, finalNblock, NW

    NW=qmc_data%NW

    if ( present(file_Etrace) ) then
       open(unit=21,file=file_Etrace)
       write(unit=21,fmt="(e14.7)") qmc_data%Etot
       close(unit=21)
    end if

    if ( present(file_NWtrace) .and. allocated(qmc_data%NWtr) ) then
       open(unit=21,file=file_NWtrace)
       write(unit=21,fmt="(i8)") qmc_data%NWtr
       close(unit=21)
    end if

    print *
    print *, trim(qmc_data%method_id)
    print *, "population size:", qmc_data%NWopt
    print *, "number of steps:", qmc_data%Nsteps
    print *, "time step:", qmc_data%tau
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

    En=sum(qmc_data%Etot)/qmc_data%Nsteps
    errEn=errorbar(qmc_data%Etot,En,file_Ereblock,corrlen,finalNBlock)
    print *, "---"
    print *, "total energy:", En, "+/-", errEn
    print *, "correlation length for E (steps):", corrlen
    print *, "final number of blocks in the statistics:", finalNblock

    En2=sum(qmc_data%Etot2)/qmc_data%Nsteps
    errEn2=errorbar(qmc_data%Etot2,En2,file_E2reblock,corrlen,finalNBlock)
    Var=En2-En**2
    errVar=sqrt(errEn2**2+(2*En*errEn)**2)
    print *, "---"
    print *, "variance:", Var, "+/-", errVar
    print *, "correlation length for E^2 (steps):", corrlen
    print *, "final number of block in the statistics:", finalNblock

    if ( allocated(qmc_data%EtotG) ) then
       En=sum(qmc_data%EtotG)/qmc_data%Nsteps
       errEn=errorbar(qmc_data%EtotG,En,file_Ereblock,corrlen,finalNBlock)
       print *, "---"
       print *, "total energy (growth):", En, "+/-", errEn
       print *, "correlation length:", corrlen
       print *, "final number of blocks in the statistics:", finalNblock
    end if

    if ( present(extra_unit) ) then
       write(unit=22,fmt='(4e18.10)') qmc_data%tau, En, errEn, &
            qmc_data%accepted_moves / (1.0_dp*qmc_data%total_moves)
    end if

    ! }}}
  end subroutine show_results

end program qmc_run


! Local variables:
! folded-file: t
! End:
