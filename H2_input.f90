! problem specific subroutines for H2 molecule

module qmc_input
  use types_const, only: dp, missing
  implicit none
  private

  public :: sys_dim
  public :: t_walker, t_sys
  public :: init_sys
  public :: EL_drift, psiT2

  ! logically, this should go to t_sys, but I really want to keep static arrays
  ! in t_walker, which enforces this to be a global constant
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
     character(11) :: id="H2 molecule"  ! identification of the system
     real(dp) :: R=1.4_dp               ! proton-proton distance
  end type t_sys

  interface EL_drift
     module procedure EL_drift_HeitlerLondon
  end interface EL_drift

  interface psiT2
     module procedure psiT2_HeitlerLondon
  end interface psiT2

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

    if ( .not.readvalue(syswords,sys%R,"R") ) call missing("R")
    call clear(syswords)

    write(unit=*,fmt='(1x,a)') "system:"
    write(unit=*,fmt='(3x,a)') sys%id
    write(unit=*,fmt='(3x,a,f8.4)') "proton-proton distance [R]:", sys%R
    write(unit=*,fmt=*)
    ! }}}
  end subroutine init_sys

  pure subroutine EL_drift_HeitlerLondon(sys,wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp), dimension(3) :: r1, r2
    real(dp) :: rho1, r1a, r1b, rho2, r2a, r2b, r12
    real(dp) :: phiA1B2, phiA2B1, psi
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    rho1=r1(1)**2+r1(2)**2
    r1a=sqrt(rho1+(r1(3)-sys%R)**2)
    r1b=sqrt(rho1+r1(3)**2)
    rho2=r2(1)**2+r2(2)**2
    r2a=sqrt(rho2+(r2(3)-sys%R)**2)
    r2b=sqrt(rho2+r2(3)**2)
    r12=sqrt(sum((r1-r2)**2))
    phiA1B2=exp(-r1a-r2b)
    phiA2B1=exp(-r2a-r1b)
    psi=phiA1B2+phiA2B1
    wlkr%psiTsq=psi**2
    wlkr%EL=-1.0_dp-((1.0_dp/r1b+1.0_dp/r2a)*phiA1B2 &
         +(1.0_dp/r1a+1.0_dp/r2b)*phiA2B1)/psi+1.0_dp/r12+1.0_dp/sys%R
    wlkr%vD(1:2)=-phiA1B2*r1(1:2)/r1a-phiA2B1*r1(1:2)/r1b
    wlkr%vD(3)=-phiA1B2*(r1(3)-sys%R)/r1a-phiA2B1*r1(3)/r1b
    wlkr%vD(4:5)=-phiA1B2*r2(1:2)/r2b-phiA2B1*r2(1:2)/r2a
    wlkr%vD(6)=-phiA1B2*r2(3)/r2b-phiA2B1*(r2(3)-sys%R)/r2a
    wlkr%vD=wlkr%vD/psi
    wlkr%sgn=1
    ! }}}
  end subroutine EL_drift_HeitlerLondon

  pure subroutine psiT2_HeitlerLondon(sys,wlkr)
    ! {{{ trial wave-function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp), dimension(3) :: r1, r2
    real(dp) :: rho1, r1a, r1b, rho2, r2a, r2b, r12
    real(dp) :: phiA1B2, phiA2B1, psi
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    rho1=r1(1)**2+r1(2)**2
    r1a=sqrt(rho1+(r1(3)-sys%R)**2)
    r1b=sqrt(rho1+r1(3)**2)
    rho2=r2(1)**2+r2(2)**2
    r2a=sqrt(rho2+(r2(3)-sys%R)**2)
    r2b=sqrt(rho2+r2(3)**2)
    r12=sqrt(sum((r1-r2)**2))
    phiA1B2=exp(-r1a-r2b)
    phiA2B1=exp(-r2a-r1b)
    psi=phiA1B2+phiA2B1
    wlkr%psiTsq=psi**2
    wlkr%sgn=1
    ! }}}
  end subroutine psiT2_HeitlerLondon

end module qmc_input


! Local variables:
! folded-file: t
! End:
