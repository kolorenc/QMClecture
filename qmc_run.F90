! ==========================================================================
!
! Basic VMC/DMC code for educational purposes.
! This file: top-level driver, compiles to the executable
!
! Copyright (C) 2012 Jindrich Kolorenc
!
! The software is released under MIT/X11 license.
!
! ==========================================================================

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
  use types_const, only: dp, missing
  use cfparser
  use qmc_input
  use qmc
  use reblocking
  implicit none

  type t_run_param
     integer :: NW=2**11
     integer :: VMCsteps=2**16
     integer :: VMCtherm=2**11
     real(dp) :: VMCtstep=0.1_dp
     integer :: DMCsteps=2**16
     integer :: DMCtherm=2**11
     real(dp) :: DMCtstep=0.01_dp
     integer :: DMCruns=1
  end type t_run_param

  integer, parameter :: extra_unit=22

  type(t_sys) :: sys
  type(t_qmc_data) :: qmc_data
  type(t_run_param) :: run_param
  type(t_words) :: words
  real(dp) :: f
  integer :: i, Nthreads

  ! watch out for the DMC acceptance ratio at larger Znuc - a smaller timestep
  ! is needed in the high-density region, or smarter moves such as [Umrigar,
  ! Nightingale & Runge, J. Chem. Phys. 99, 2865 (1993)]
  call readConf(words,"qmc_run.cfg")
  call init_run_param(run_param,words)
  call init_sys(sys,words)
  call clear(words)

#ifdef _OPENMP
  ! the return value of omp_get_max_threads() is controlled by
  !   export OMP_NUM_THREADS=n
  Nthreads=omp_get_max_threads()
#else
  Nthreads=1
#endif
  print *, "number of threads:", Nthreads
  call init_qmc_data(qmc_data,sys,run_param%NW,Nthreads)

  ! thermalization
  call init_qmc(qmc_data,run_param%VMCtstep,run_param%VMCtherm)
  call vmc_run(qmc_data,sys)

  ! data harvest
  call init_qmc(qmc_data,run_param%VMCtstep,run_param%VMCsteps)
  call vmc_run(qmc_data,sys)
  call show_results(qmc_data,"Evmc_reblock.dat","E2vmc_reblock.dat", &
       "Evmc_trace.dat")

  ! thermalization
  call init_qmc(qmc_data,run_param%DMCtstep,run_param%DMCtherm)
  call dmc_run(qmc_data,sys)

  ! data harvest
  !call init_qmc(qmc_data,run_param%DMCtstep,run_param%DMCsteps)
  !call dmc_run(qmc_data,sys)
  !call show_results(qmc_data,"Edmc_reblock.dat","E2dmc_reblock.dat", &
  !     "Edmc_trace.dat","NWdmc_trace.dat")

  ! data harvest (with optional time-step extrapolation)
  open(unit=extra_unit,file="tstep_extrap.dat")
  write(unit=extra_unit,fmt=*) "# time step   energy   errorbar   acceptance"
  do i=0, run_param%DMCruns-1
     f=(run_param%DMCruns-i)/real(run_param%DMCruns,dp)
     call init_qmc(qmc_data,run_param%DMCtstep*f, &
          int(run_param%DMCsteps/f))
     call dmc_run(qmc_data,sys)
     call show_results(qmc_data,"Edmc_reblock.dat","E2dmc_reblock.dat", &
          "Edmc_trace.dat","NWdmc_trace.dat",extra_unit)
  end do
  close(unit=extra_unit)

contains

  subroutine init_run_param(run_param,words)
    ! {{{ take run parameters from the list of words
    type(t_run_param), intent(out) :: run_param
    type(t_words), intent(in) :: words
    character(len=8) :: var

    var="NW"
    if ( .not.readvalue(words,run_param%NW,var) ) call missing(var)

    var="VMCtherm"
    if ( .not.readvalue(words,run_param%VMCtherm,var) ) call missing(var)
    var="VMCsteps"
    if ( .not.readvalue(words,run_param%VMCsteps,var) ) call missing(var)
    var="VMCtstep"
    if ( .not.readvalue(words,run_param%VMCtstep,var) ) call missing(var)

    var="DMCtherm"
    if ( .not.readvalue(words,run_param%DMCtherm,var) ) call missing(var)
    var="DMCsteps"
    if ( .not.readvalue(words,run_param%DMCsteps,var) ) call missing(var)
    var="DMCtstep"
    if ( .not.readvalue(words,run_param%DMCtstep,var) ) call missing(var)
    var="DMCruns"
    if ( .not.readvalue(words,run_param%DMCruns,var) ) call missing(var)

    run_param%NW=2**run_param%NW
    run_param%VMCsteps=2**run_param%VMCsteps
    run_param%VMCtherm=2**run_param%VMCtherm
    run_param%DMCsteps=2**run_param%DMCsteps
    run_param%DMCtherm=2**run_param%DMCtherm

    ! }}}
  end subroutine init_run_param

  subroutine show_results(qmc_data,file_Ereblock,file_E2reblock, &
       file_Etrace,file_NWtrace,extra_unit)
    ! {{{ calculate and write the measured quantities, their errorbars and
    !     statistics of the Monte Carlo process;
    !     argument qmc_data is 'inout' due to the implementation of errorbar()
    type(t_qmc_data), intent(inout) :: qmc_data
    character(len=*), intent(in) :: file_Ereblock, file_E2reblock
    character(len=*), optional, intent(in) :: file_Etrace, file_NWtrace
    integer, optional, intent(in) :: extra_unit
    real(dp) :: En, errEn, En2, errEn2, Var, errVar, corrlen, corrlenEn, ageAV
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
    errEn=errorbar(qmc_data%Etot,En,file_Ereblock,corrlenEn,finalNBlock)
    print *, "---"
    print *, "total energy:", En, "+/-", errEn
    print *, "correlation length for E (steps):", corrlenEn
    print *, "final number of blocks in the statistics:", finalNblock

    En2=sum(qmc_data%Etot2)/qmc_data%Nsteps
    errEn2=errorbar(qmc_data%Etot2,En2,file_E2reblock,corrlen,finalNBlock)
    Var=En2-En**2
    errVar=sqrt(errEn2**2+(2*En*errEn)**2)
    print *, "---"
    print *, "variance:", Var, "+/-", errVar
    print *, "correlation length for E^2 (steps):", corrlen
    print *, "final number of blocks in the statistics:", finalNblock

    if ( allocated(qmc_data%EtotG) ) then
       En=sum(qmc_data%EtotG)/qmc_data%Nsteps
       errEn=errorbar(qmc_data%EtotG,En,file_Ereblock,corrlen,finalNBlock)
       print *, "---"
       print *, "total energy (growth):", En, "+/-", errEn
       print *, "correlation length:", corrlen
       print *, "final number of blocks in the statistics:", finalNblock
    end if

    if ( present(extra_unit) ) then
       write(unit=extra_unit,fmt='(4e18.10)') qmc_data%tau, En, errEn, &
            qmc_data%accepted_moves / (1.0_dp*qmc_data%total_moves)
    end if

    ! }}}
  end subroutine show_results

end program qmc_run


! Local variables:
! folded-file: t
! End:
