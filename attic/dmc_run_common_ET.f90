! Here we use common trial energy ET for all nodes. This way it is more fun,
! definitelly. Distribution of walkers among nodes needs slight adjustments
! during run, otherways the populations of walkers on different nodes tend
! to differ considerably, which makes computation rather inefficient. 
program dmc_run
  use nrtype
  use dmc
  use dmc_input
  use nr_random
  !use random
  implicit none

  include 'mpif.h'

  real(sp) :: E_dmc, E_dmc2, dE_dmc, E_vmc, dE_vmc
  real(sp), dimension(:), allocatable :: Earr, ETarr

  integer, dimension(:), allocatable :: NW_nodes, NW_opt_nodes, dNW_nodes
  real(sp) :: ET, energy
  integer:: NW_opt, t_therm, t_harvest
  real(sp) :: tau
  integer :: NW_local, NW_opt_local, NW, t, i
  integer :: too_old_local, too_old
  real(sp) :: db_accept_local, db_accept, db_reject_local, db_reject
  real(sp) :: db_ratio

  integer :: ierr, rank, comm_size

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,comm_size,ierr)

  if ( rank == 0 ) then

     open(unit=101,file="dmc.cfg",action="read")
     read(unit=101,fmt=*) NW_opt, tau, t_therm, t_harvest
     close(unit=101)

     open(unit=101,file="population.log",action="write")
     open(unit=102,file="transfer.log",action="write")

     print *, "-------------------"
     print *, "number of walkers: ", NW_opt
     print *, "time step:         ", tau
     print *, "projection time:   ", t_therm, "steps"
     print *, "harvest time:      ", t_harvest, "steps"
     print *, "-------------------"
     write(unit=101,fmt=*) "# -------------------"
     write(unit=101,fmt=*) "# number of walkers: ", NW_opt
     write(unit=101,fmt=*) "# time step:         ", tau
     write(unit=101,fmt=*) "# projection time:   ", t_therm, "steps"
     write(unit=101,fmt=*) "# harvest time:      ", t_harvest, "steps"
     write(unit=101,fmt=*) "# -------------------"

     allocate( Earr(t_harvest), ETarr(t_harvest), NW_nodes(0:comm_size-1), &
          & NW_opt_nodes(0:comm_size-1), &
          & dNW_nodes(0:comm_size-1), stat=ierr )
     if ( ierr /= 0 ) then
        print *, "dmc_run: cannot allocate arrays."
        stop
     end if
  end if

  call MPI_BCAST(tau,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t_therm,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(t_harvest,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  call random_init()
  !call workdistr_bogo()
  call workdistr_etalon()

  ! we can start the calculation now
  NW_local=NW_opt_local
  call dmc_init(NW_local,tau)

  Z1=2.0_dp
  call dmc_init_walkers(E_vmc)
  energy=E_vmc*NW_local
  call MPI_REDUCE(energy,ET,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  if ( rank == 0 ) then
     NW=sum(NW_opt_nodes)
     ET=ET/NW
     E_vmc=ET
  end if
  call MPI_BCAST(ET,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

  do t=1, t_therm
     call step(.false.)
  end do
  !dens_loc=0
  do t=1, t_harvest
     call step(.true.)
  end do

  call dmc_get_db_stats(db_accept_local,db_reject_local,too_old_local)
  call MPI_REDUCE(db_accept_local,db_accept, &
       & 1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(db_reject_local,db_reject, &
       & 1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  call MPI_REDUCE(too_old_local,too_old, &
       & 1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

  if ( rank==0 ) then
     E_dmc=real(sum(real(Earr,dp))/t_harvest,sp) ! sums to be done in double
     E_dmc2=real(sum(real(ETarr,dp))/t_harvest,sp)
     dE_dmc=real(sqrt(sum(real((Earr-E_dmc)**2,dp))/(t_harvest-1)),sp)
     db_ratio=db_accept*100.0_sp/(db_accept+db_reject)

     print *, "VMC:               ", E_vmc
     print *, "DMC (potential):   ", E_dmc," +/-", dE_dmc
     print *, "DMC (trial energy):", E_dmc2," +/-", dE_dmc
     print *, "Detailed balance accepted", db_ratio,"% moves."
     print *, "Encountered", too_old, "too old (trapped) walkers."
     print *, "-------------------"
     write(unit=101,fmt=*) "# VMC:               ", E_vmc
     write(unit=101,fmt=*) "# DMC (potential):   ", E_dmc, &
          & " +/-", dE_dmc
     write(unit=101,fmt=*) "# DMC (trial energy):", E_dmc2, &
          & " +/-", dE_dmc
     write(unit=101,fmt=*) "# Detailed balance accepted", db_ratio,"% moves."
     write(unit=101,fmt=*) "# Encountered", too_old, &
          & "too old (trapped) walkers."
     write(unit=101,fmt=*) "# -------------------"

     close(unit=101)
     close(unit=102)
     
  end if

  call dump_walkers()


contains

  subroutine random_init()
    ! {{{ initialization/desynchronization of random numbers
    real(sp), dimension(:), allocatable :: rnd_init
    integer, dimension(:), allocatable :: rnd_seq
    integer :: rnd_seq_local
    if ( rank == 0 ) then
       allocate(rnd_init(0:comm_size-1), rnd_seq(0:comm_size-1), &
            & stat=ierr)
       if ( ierr /= 0 ) then
          print *, "dmc_run: cannot allocate random seed arrays."
          stop
       end if
       call ran1(rnd_init)
       rnd_seq=int(1.0e8_sp*rnd_init) ! max of 4-byte integer is cca 2e9
    end if
    call MPI_SCATTER(rnd_seq,1,MPI_INTEGER,rnd_seq_local,1, &
         & MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call ran_seed(sequence=rnd_seq_local)
    if ( rank == 0 ) deallocate(rnd_init,rnd_seq)
    ! }}}
  end subroutine random_init

  subroutine workdistr_bogo()
    ! {{{ optimal work distribution (delivers NW_opt_local to nodes)
    real(sp), dimension(:), allocatable :: bogo_nodes
    real(sp) :: bogo, bogo_local
    character(255) :: line
    integer :: i
    if ( rank == 0 ) then
       allocate(bogo_nodes(0:comm_size-1), stat=ierr)
       if ( ierr /= 0 ) then
          print *, "dmc_run: cannot allocate bogoMIPS arrays."
          stop
       end if
    end if
    open(unit=103,file="/home/scratch/kolorenc/qmc_mpi/cpuinfo", &
         & action="read",form="formatted")
    do i=1, 100
       read(unit=103,fmt="(A)") line
       if ( index(line,"ogo")>0 ) exit
    end do
    backspace(unit=103)
    read(unit=103,fmt=*) line, bogo_local
    close(unit=103)
    call MPI_GATHER(bogo_local,1,MPI_REAL,bogo_nodes,1,MPI_REAL,0, &
         & MPI_COMM_WORLD,ierr)
    if ( rank == 0 ) then
       bogo=sum(bogo_nodes)
       print *, "Total bogoMIPS:    ", bogo," on", comm_size, "processors."
       print *, "-------------------"
       write(unit=101,fmt=*) "# Total bogoMIPS:    ", bogo," on", &
            & comm_size, "processors."
       write(unit=101,fmt=*) "# -------------------"
       NW_opt_nodes=int(NW_opt*bogo_nodes/bogo)
    end if
    call MPI_SCATTER(NW_opt_nodes,1,MPI_INTEGER,NW_opt_local,1, &
         & MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( rank == 0 ) deallocate(bogo_nodes)
    ! }}}
  end subroutine workdistr_bogo

  subroutine workdistr_etalon()
    ! {{{ optimal work distribution (delivers NW_opt_local to nodes)
    real(sp), dimension(:), allocatable :: etalon_nodes
    real(sp) :: total_power, etalon_local
    if ( rank == 0 ) then
       allocate(etalon_nodes(0:comm_size-1), stat=ierr)
       if ( ierr /= 0 ) then
          print *, "dmc_run: cannot allocate array for speed etalon."
          stop
       end if
    end if
    open(unit=103,file="/home/scratch/kolorenc/qmc_mpi/speed_etalon.txt", &
         & action="read",form="formatted")
    read(unit=103,fmt=*) etalon_local
    close(unit=103)
    etalon_local=1.0_sp/etalon_local
    call MPI_GATHER(etalon_local,1,MPI_REAL,etalon_nodes,1,MPI_REAL,0, &
         & MPI_COMM_WORLD,ierr)
    if ( rank == 0 ) then
       total_power=sum(etalon_nodes)
       NW_opt_nodes=int(NW_opt*etalon_nodes/total_power)
       print *, NW_opt_nodes
    end if
    call MPI_SCATTER(NW_opt_nodes,1,MPI_INTEGER,NW_opt_local,1, &
         & MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if ( rank == 0 ) deallocate(etalon_nodes)
    ! }}}
  end subroutine workdistr_etalon

  subroutine step(harvest)
    ! {{{ .
    logical, intent(in) :: harvest
    integer, dimension(MPI_STATUS_SIZE) :: status
    type(t_walker) :: walker
    integer :: dNW_min_rank, dNW_max_rank
    integer :: NW2move, i, dNW_min, dNW_max, NW_max
    !
    call dmc_step(ET,energy,NW_local)
    energy=energy*NW_local
    call MPI_REDUCE(energy,ET,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_GATHER(NW_local,1,MPI_INTEGER,NW_nodes,1,MPI_INTEGER,0, &
         & MPI_COMM_WORLD,ierr)
    if ( rank == 0 ) then
       NW=sum(NW_nodes)
       energy=ET/NW
       ET=energy-log(NW*1.0_sp/NW_opt)*0.025_sp
       if ( harvest ) then
          Earr(t)=energy
          ETarr(t)=ET
          write(unit=101,fmt="(2i10,e18.6,15i10)") t_therm+t, NW, &
               & energy, NW_nodes
          !call dens_add()
       else
          write(unit=101,fmt="(2i10,e18.6,15i10)") t, NW, energy, NW_nodes
       end if
       dNW_nodes=NW_nodes-NW_opt_nodes
       dNW_min_rank=0
       dNW_max_rank=0
       dNW_min=dNW_nodes(0)
       dNW_max=dNW_nodes(0)
       do i=1, comm_size-1
          if ( dNW_nodes(i) < dNW_min ) then
             dNW_min_rank=i
             dNW_min=dNW_nodes(i)
          end if
          if ( dNW_nodes(i) > dNW_max ) then
             dNW_max_rank=i
             dNW_max=dNW_nodes(i)
          end if
       end do
       NW2move=(dNW_max-dNW_min)/50
       NW_max=NW_nodes(dNW_max_rank)
       write(unit=102,fmt=*) NW2move, "walkers", &
            & dNW_max_rank, "-->", dNW_min_rank
    end if

    call MPI_BCAST(ET,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(NW2move,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

    if ( NW2move > 0 ) then

       call MPI_BCAST(dNW_min_rank,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(dNW_max_rank,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(NW_max,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       do i=NW_max, NW_max-NW2move,-1
          if ( rank == dNW_max_rank ) then
             call get_walker(i,walker)
             call kill_walker(i)
             call MPI_SEND(walker%r,dmc_dim,MPI_REAL,dNW_min_rank,1001, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%vD,dmc_dim,MPI_REAL,dNW_min_rank,1002, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%EL,1,MPI_REAL,dNW_min_rank,1003, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%psiTsq,1,MPI_REAL,dNW_min_rank,1004, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%sgn,1,MPI_INTEGER,dNW_min_rank,1005, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%age,1,MPI_INTEGER,dNW_min_rank,1006, &
                  & MPI_COMM_WORLD,ierr)
          else if (rank == dNW_min_rank ) then
             call MPI_RECV(walker%r,dmc_dim,MPI_REAL,dNW_max_rank,1001, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%vD,dmc_dim,MPI_REAL,dNW_max_rank,1002, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%EL,1,MPI_REAL,dNW_max_rank,1003, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%psiTsq,1,MPI_REAL,dNW_max_rank,1004, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%sgn,1,MPI_INTEGER,dNW_max_rank,1005, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%age,1,MPI_INTEGER,dNW_max_rank,1006, &
                  & MPI_COMM_WORLD,status,ierr)
             call add_walker(walker)
          end if
       end do

    end if
    ! }}}
  end subroutine step

  subroutine dump_walkers()
    ! {{{ dump last walker positions
    type(t_walker) :: walker
    integer :: i, j, NW_j
    integer, dimension(MPI_STATUS_SIZE) :: status
    real(dp) :: EL_sum2, EL_sum_local
    real(sp) :: EL_sum
    real(sp), dimension(comm_size) :: EL_sum_nodes
    
    ! just to check, whether we do the MPI send/receive correctly so that
    ! we get all the walkers fine
    EL_sum_local=0.0_dp
    do i=1, NW_local
       call get_walker(i,walker)
       EL_sum_local=EL_sum_local+real(walker%EL,dp)
    end do
    call MPI_GATHER(real(EL_sum_local,sp),1,MPI_REAL,EL_sum_nodes,1,&
         & MPI_REAL,0, MPI_COMM_WORLD,ierr)
    EL_sum=sum(EL_sum_nodes)


    call MPI_GATHER(NW_local,1,MPI_INTEGER,NW_nodes,1,MPI_INTEGER,0, &
         & MPI_COMM_WORLD,ierr)

    if ( rank == 0 ) then
       EL_sum2=0.0_dp
       NW=sum(NW_nodes)
       open(unit=101,file="walker_positions.dat", &
            & action="write",form="unformatted")
       write(unit=101) NW
       do i=1, NW_nodes(0)
          call get_walker(i,walker)
          write(unit=101) walker%r
          EL_sum2=EL_sum2+real(walker%EL,dp)
       end do
    end if

    do j=1, comm_size-1
       if ( rank == 0 ) then
          !print *, j
          NW_j=NW_nodes(j)
       end if
       call MPI_BCAST(NW_j,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       do i=1, NW_j
          if ( rank == j ) then
             call get_walker(i,walker)
             call MPI_SEND(walker%r,dmc_dim,MPI_REAL,0,1001, &
                  & MPI_COMM_WORLD,ierr)
             call MPI_SEND(walker%EL,1,MPI_REAL,0,1003, &
                  & MPI_COMM_WORLD,ierr)
          else if ( rank == 0 ) then
             call MPI_RECV(walker%r,dmc_dim,MPI_REAL,j,1001, &
                  & MPI_COMM_WORLD,status,ierr)
             call MPI_RECV(walker%EL,1,MPI_REAL,j,1003, &
                  & MPI_COMM_WORLD,status,ierr)
             write(unit=101) walker%r
             EL_sum2=EL_sum2+real(walker%EL,dp)
          end if
          !call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       end do
    end do

    if ( rank == 0 ) then
       print *, EL_sum, EL_sum2
       close(unit=101)
    end if

    ! }}}
  end subroutine dump_walkers

end program dmc_run


! Local variables:
! folded-file: t
