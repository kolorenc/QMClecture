! problem specific subroutines for He atom (orthohelium=parallel spins)
module qmc_input
  use types_const, only: dp, missing
  implicit none
  private

  public :: sys_dim
  public :: t_walker, t_sys
  public :: init_sys
  public :: EL_drift, psiT2

  integer, parameter :: sys_dim=6       ! dimension of walker

  type t_walker
     real(dp), dimension(sys_dim) :: r  ! position
     real(dp), dimension(sys_dim) :: vD ! drift velocity
     real(dp) :: EL                     ! local energy
     real(dp) :: psiTsq                 ! square of trial wavefunction
     real(dp) :: wt=1.0_dp              ! DMC weight
     integer :: sgn=1                   ! sign of trial wavefunction at r
     integer :: age=0                   ! how many times rejected by det. bal.
  end type t_walker

  type t_sys
     character(11) :: id="orthohelium"  ! identification of the system
     real(dp) :: Znuc=2.0_dp            ! nuclear charge
     ! effective charges optimized in VMC
     !real(dp) :: Z1=1.98_dp            ! effective charge for 1s wave
     !real(dp) :: Z2=1.55_dp            ! effective charge for 2s wave
     ! from analytics in Mathematica, corresponding energy is -2.166639875
     !real(dp) :: Z1=1.99364_dp
     !real(dp) :: Z2=1.55094_dp
     ! our trial wave function has the same nodes for all charges; it is thus
     ! actually better to take unscreened charges because in that case there
     ! are no 1/r divergences in EL due to core (remains only that due to e-e
     ! interaction)
     real(dp) :: Z1=2.0_dp
     real(dp) :: Z2=2.0_dp
  end type t_sys

contains

  subroutine init_sys(sys,words)
    ! {{{ initialize the system and wave-function parameters
    use cfparser, only: t_words, readsection, readvalue, clear
    type(t_sys), intent(inout) :: sys
    type(t_words), intent(in) :: words
    type(t_words) :: syswords

    if ( .not.readsection(words,syswords,"system") ) then
       print *, "Missing section 'system' in the input file."
       stop
    end if

    if ( .not.readvalue(syswords,sys%Znuc,"Znuc") ) call missing("Znuc")
    if ( .not.readvalue(syswords,sys%Z1,"Z1") ) sys%Z1 = sys%Znuc
    if ( .not.readvalue(syswords,sys%Z2,"Z2") ) sys%Z2 = sys%Znuc
    call clear(syswords)

    write(unit=*,fmt='(1x,a)') "system:"
    write(unit=*,fmt='(3x,a)') sys%id
    write(unit=*,fmt='(3x,a,f8.4)') "nuclear charge [Znuc]:", sys%Znuc
    write(unit=*,fmt=*)
    ! }}}
  end subroutine init_sys

  pure subroutine EL_drift(sys,wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l, psiA, psiB, psi, s12, s21
    real(dp), dimension(3) :: r1, r2, r10, r20
    real(dp) :: Z2half
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12l=sqrt(sum((r1-r2)**2))
    r10=r1/r1l
    r20=r2/r2l
    Z2half=sys%Z2/2.0_dp
    s12=exp(-sys%Z1*r1l-Z2half*r2l)
    s21=exp(-sys%Z1*r2l-Z2half*r1l)
    psiA=(1.0_dp-Z2half*r2l)*s12
    psiB=(1.0_dp-Z2half*r1l)*s21
    psi=psiA-psiB
    wlkr%sgn=int(sign(1.0_dp,psi))
    wlkr%EL=-(sys%Z1**2+Z2half**2)/2.0_dp+1.0_dp/r12l &
         +( ( (sys%Z1-sys%Znuc)/r1l+(sys%Z2-sys%Znuc)/r2l )*psiA &
         -( (sys%Z1-sys%Znuc)/r2l+(sys%Z2-sys%Znuc)/r1l )*psiB )/psi
    wlkr%psiTsq=psi**2
    wlkr%vD(1:3)=( -sys%Z1*psiA + Z2half*psiB + Z2half*s21 )*r10
    wlkr%vD(4:6)=(  sys%Z1*psiB - Z2half*psiA - Z2half*s12 )*r20
    wlkr%vD=wlkr%vD/psi
    ! }}}
  end subroutine EL_drift

  pure subroutine psiT2(sys,wlkr)
    ! {{{ trial wave-function squared and local energy
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, psiA, psiB, psi
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    psiA=(1.0_dp-sys%Z2*r2l/2.0_dp)*exp(-sys%Z1*r1l-sys%Z2*r2l/2.0_dp)
    psiB=(1.0_dp-sys%Z2*r1l/2.0_dp)*exp(-sys%Z1*r2l-sys%Z2*r1l/2.0_dp)
    psi=psiA-psiB
    wlkr%psiTsq=psi**2
    ! }}}
  end subroutine psiT2

end module qmc_input


! Local variables:
! folded-file: t
! End:
