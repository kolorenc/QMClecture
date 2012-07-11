! Single node (non-MPI) implementation of DMC using subroutines from modules
! dmc and amc_input. Mainly for benchmarking the MPI version.
program dmc_run
  use nrtype
  use dmc
  use dmc_input
  use nr_random
  implicit none
  
  integer :: NW_opt, t_therm, t_harvest, t_block, nblock
  real(sp) :: tau
  
  integer, parameter :: dens_res=200
  integer, dimension(dens_res,dens_res) :: dens
  real(sp), parameter :: xMin=-1.2_sp, xMax=1.2_sp
  real(sp), parameter :: zMin=-1.0_sp, zMax=2.4_sp
  real(sp), parameter :: yMin=-0.05_sp, yMax=0.05_sp
  integer :: i, j
  real(sp) :: x, z

  real(sp), dimension(:), allocatable :: Earr, Earr_block
  real(sp) :: ET, energy, db_accept, db_reject
  real(sp) :: E_dmc, dE_dmc, E_vmc, db_ratio
  integer :: too_old, NW, t, ierr

  open(unit=101,file="dmc.cfg",action="read")
  read(unit=101,fmt=*) NW_opt, tau, t_therm, nblock, t_block
  close(unit=101)
  t_harvest=nblock*t_block
  
  open(unit=101,file="population.log",action="write")

  write(unit=101,fmt=*) "# -------------------"
  write(unit=101,fmt=*) "# number of walkers: ", NW_opt
  write(unit=101,fmt=*) "# time step:         ", tau
  write(unit=101,fmt=*) "# projection time:   ", t_therm, "steps"
  write(unit=101,fmt=*) "# number of blocks:  ", nblock
  write(unit=101,fmt=*) "# time per block:    ", t_block, "steps"
  write(unit=101,fmt=*) "# -------------------"

  allocate( Earr(t_harvest), Earr_block(nblock), stat=ierr )
  if ( ierr /=0 ) then
     print *, "dmc_run: cannot allocate arrays."
     stop
  end if

  call dmc_init(NW_opt,tau)
  call dmc_init_walkers(E_vmc)
  ET=E_vmc
  ! termalize first
  do t=1, t_therm
     call dmc_step(ET,energy,NW)
     write(unit=101,fmt="(2i10,e18.6)") t, NW, energy
     ET=energy-log(NW*1.0_sp/NW_opt)
  end do
  ! collect averages
  do i=1, nblock
     do t=(i-1)*t_block+1, i*t_block
        call dmc_step(ET,energy,NW)
        Earr(t)=energy
        write(unit=101,fmt="(2i10,e18.6)") t+t_therm, NW, energy
        !call dens_add()
        ET=energy-log(NW*1.0_sp/NW_opt)
     end do
     Earr_block(i)=real(sum(real(Earr((i-1)*t_block+1:i*t_block),dp)) &
          & /t_block,sp)
     dE_dmc=real(sqrt(sum(real((Earr((i-1)*t_block+1:i*t_block) &
          & -Earr_block(i))**2,dp))/(t_block-1)),sp)
     print *, "block:", Earr_block(i)," +/-", dE_dmc
     write(unit=101,fmt=*) "# block:", Earr_block(i)," +/-", dE_dmc
  end do
  ! calculate results
  E_dmc=real(sum(real(Earr_block,dp))/nblock,sp) ! sums to be done in double
  dE_dmc=real(sqrt(sum(real((Earr_block-E_dmc)**2,dp))/(nblock-1)),sp)
  call dmc_get_db_stats(db_accept,db_reject,too_old)
  db_ratio=(db_accept-db_reject)*100.0_sp/(db_accept+db_reject)

  print *, "VMC:", E_vmc
  print *, "DMC:", E_dmc," +/-", dE_dmc
  print *, "Detailed balance accepted", db_ratio,"% moves."
  print *, "Encountered", too_old, "too old (trapped) walkers."
  write(unit=101,fmt=*) "# -------------------"
  write(unit=101,fmt=*) "# VMC:               ", E_vmc
  write(unit=101,fmt=*) "# DMC:               ", E_dmc, &
       & " +/-", dE_dmc
  write(unit=101,fmt=*) "# Detailed balance accepted", db_ratio,"% moves."
  write(unit=101,fmt=*) "# Encountered", too_old, &
       & "too old (trapped) walkers."
  write(unit=101,fmt=*) "# -------------------"

  close(unit=101)

  stop

  open(unit=101,file="density_gnuplot.dat",action="write")
  do i=1, dens_res
     x=xMin+(i-0.5_sp)*(xMax-xMin)/dens_res
     do j=1, dens_res
        z=zMin+(j-0.5_sp)*(zMax-zMin)/dens_res
        write(unit=101,fmt=*) x, z, dens(i,j)
     end do
     write(unit=101,fmt=*) ' '
  end do
  close(unit=101)
  open(unit=101,file="density_mathematica.dat",action="write")
  do i=1, dens_res
     write(unit=101,fmt="(500i)") dens(i,:)
  end do
  close(unit=101)


contains

  subroutine dens_add()
    ! {{{ .
    type(t_walker) :: walker
    integer :: ix, iz
    do i=1, NW
       call get_walker(i,walker)
       do j=1, dmc_dim, 3
          if ( (walker%r(j)>xMin).and.(walker%r(j)<xMax) &
               & .and.(walker%r(j+1)>yMin).and.(walker%r(j+1)<yMax) &
               & .and.(walker%r(j+2)>zMin).and.(walker%r(j+2)<zMax) ) then
             ix=int((walker%r(j)-xMin)/(xMax-xMin)*dens_res)+1
             iz=int((walker%r(j+2)-zMin)/(zMax-zMin)*dens_res)+1
             dens(ix,iz)=dens(ix,iz)+1
          end if
       end do
    end do
    open(unit=102,file="dens.dat",action="write",form="unformatted")
    write(unit=102) dens
    close(unit=102)
    ! }}}
  end subroutine dens_add

end program dmc_run


! Local variables:
! folded-file: t
! End:
