! Histogram in gnuplot:
! bin(x,width)=width*floor(x/width) + width/2.0
! plot 'Evmc_trace.dat' using (bin($1,.001)):(1.0) smooth freq with boxes

program vmc
  use omp_lib, only: omp_get_max_threads
  use types_const, only: dp
  use qmc
  use reblocking
  implicit none

  type(t_qmc_data) :: vmc_data
  integer, parameter :: NW=2**11
  integer, parameter :: Nsteps=2**16
  integer :: Nthreads
  !real(dp), parameter :: timestep=0.8_dp ! for pure Metropolis
  real(dp), parameter :: timestep=0.1_dp
  real(dp) :: En, errEn, corrlen, ageAV
  integer :: i, ageM, finalNblock

  ! the return value of omp_get_max_threads() is controlled by
  !   export OMP_NUM_THREADS=n
  Nthreads=omp_get_max_threads()
  print *, "number of threads:", Nthreads
  
  call init_qmc_data(vmc_data,NW,Nthreads)

  ! thermalization
  call init_qmc(vmc_data,timestep,2048)
  call vmc_run(vmc_data)

  ! actual data harvest
  call init_qmc(vmc_data,timestep,Nsteps)
  call vmc_run(vmc_data)

  open(unit=21,file="Evmc_trace.dat")
  write(unit=21,fmt="(e14.7)") vmc_data%Etot
  close(unit=21)

  En=sum(vmc_data%Etot)/Nsteps
  errEn=errorbar(vmc_data%Etot,En,"Evmc_reblock.dat",corrlen,finalNBlock)
  print *, "total energy:", En, "+/-", errEn
  print *, "acceptance ratio:", &
       vmc_data%accepted_moves / (1.0_dp*vmc_data%total_moves)
  print *, "correlation length:", corrlen
  print *, "final number of block in the statistics:", finalNblock

  ageM=vmc_data%walker(1)%age
  ageAV=ageM
  do i=2, NW
     if ( vmc_data%walker(i)%age > ageM ) ageM=vmc_data%walker(i)%age
     ageAV=ageAV+vmc_data%walker(i)%age
  end do
  ageAV=ageAV/NW
  print *, "maximal age:", ageM
  print *, "average age:", ageAV

end program vmc
