! ==========================================================================
!
! Basic VMC/DMC code for educational purposes.
! This file: system-specific subroutines for He atom (parahelium = spin
! singlet), it is good for any helium-like ion, for instance H^-
!
! Copyright (C) 2012 Jindrich Kolorenc
!
! The software is released under MIT/X11 license.
!
! ==========================================================================

! Trial wave function:
! Hylleraas-type ansatz with a single variational parameter 'a'. The wave
! function fulfills electron-ion as well as electron-electron cusp conditions.
!
! \begin{equation}
! \Psi_\text{T}(\mathbf r_1,\mathbf r_2)
! =\biggl(1+\frac12 r_{12}\biggr)\Bigl[1+a\bigl(r_1^2+r_2^2\bigr)\Bigr]
!   \,\rme^{-Z(r_1+r_2)}
! \end{equation}
!
! The total energy can be calculated analytically. The optimal value of the
! variational parameter is a=0.0878628 and gives the energy -2.89537.

! some of the variational energies we know from other sources
!  "no cusp", Z1=2            -2.75
!  "no cusp", Z1=2-5/16       -2.84766
!  "all cusps", a=0.0878628   -2.89537

! Hartree-Fock limit: -2.8617
! [Physical chemistry: a molecular approach By Donald Allan McQuarrie,
! John Douglas Simon, p. 283]

! Exact energy:   -2.903724
! [C. Schwartz, arXiv:physics/0208004v1, published as Int. J. Mod. Phys. E 15,
! 877 (2006)]

module qmc_input
  use types_const, only: dp, finish_line
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
     ! identification of the system
     character(48) :: id="parahelium (wave function with effective charge)"

     ! nuclear charge
     real(dp) :: Znuc=2.0_dp

     ! parameters in the trial wave function
     real(dp) :: Z1=27.0_dp/16.0_dp     ! effective charge in the s wave
     real(dp) :: b=0.5_dp               ! denominator in the ee-cusp term
  end type t_sys

  interface EL_drift
     module procedure EL_drift_Zstar
     ! module procedure EL_drift_Zstar_ee_cusp
  end interface EL_drift

  interface psiT2
     module procedure psiT2_Zstar
     ! module procedure psiT2_Zstar_ee_cusp
  end interface psiT2

contains

  subroutine init_sys(sys,words)
    ! {{{ initialize the system and wave-function parameters
    use cfparser, only: t_words, readsection, readvalue, clear
    type(t_sys), intent(inout) :: sys
    type(t_words), intent(in) :: words
    type(t_words) :: syswords
    logical :: Znuc_given, Z1_given

    ! read the parameters of the system and the trial wave function from the
    ! configuration file
    Znuc_given=.false.
    Z1_given=.false.
    if ( readsection(words,syswords,"system") ) then
       Znuc_given=readvalue(syswords,sys%Znuc,"Znuc")
       Z1_given=readvalue(syswords,sys%Z1,"Z1")
    end if
    call clear(syswords)

    ! set helium (Z=2) if the nuclear charge is not given
    if ( .not.Znuc_given ) sys%Znuc=2.0_dp

    ! if the variational parameter a is not given then use the optimal value
    if ( .not.Z1_given ) then
       sys%Z1=sys%Znuc-5.0_dp/16.0_dp
    end if

    ! output the settings to STDOUT
    write(unit=*,fmt='(1x,a)') "system and trial wave function:"
    write(unit=*,fmt='(3x,a)') sys%id

    write(unit=*,fmt='(3x,a,f12.8)',advance="no") &
            "nuclear charge                 [Znuc]:", sys%Znuc
    call finish_line(Znuc_given,"(default)")

    write(unit=*,fmt='(3x,a,f12.8)',advance="no") &
            "charge in the wave function    [Z1]  :", sys%Z1
    call finish_line(Z1_given,"(the optimal value)")

    write(unit=*,fmt=*)
    ! }}}
  end subroutine init_sys

  ! ===========================================================================
  ! simple trial wave function, just a product of two s orbitals with effective
  ! charge; all cusps are broken when the effective charge deviates from 2;
  ! the earliest reference I could find is [Kellner, Z. Phys. 44, 91 (1927)]
  ! ===========================================================================

  pure subroutine EL_drift_Zstar(sys,wlkr)
    ! {{{ local energy, drift velocity, wave function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12l=sqrt(sum((r1-r2)**2))
    wlkr%EL=-sys%Z1**2-(sys%Znuc-sys%Z1)*(1.0_dp/r1l+1.0_dp/r2l)+1.0_dp/r12l
    wlkr%psiTsq=exp(-2.0_dp*sys%Z1*(r1l+r2l))
    wlkr%vD(1:3)=-sys%Z1*r1/r1l
    wlkr%vD(4:6)=-sys%Z1*r2/r2l
    wlkr%sgn=1
    ! }}}
  end subroutine EL_drift_Zstar

  pure subroutine psiT2_Zstar(sys,wlkr)
    ! {{{ trial wave function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    wlkr%psiTsq=exp(-2.0_dp*sys%Z1*(r1l+r2l))
    ! }}}
  end subroutine psiT2_Zstar


  ! ===========================================================================
  ! the above function multiplied by an e-e cusp-fixing Jastrow factor,
  ! see for instance [Reynolds, Ceperley, Alder & Lester, J. Chem. Phys. 77,
  ! 5593 (1982)] but this functional form is surely older
  ! ===========================================================================

  pure subroutine EL_drift_Zstar_ee_cusp(sys,wlkr)
    ! {{{ local energy, drift velocity, wave function squared, implemented
    !     electron-electron cusp
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l, denom
    real(dp), dimension(3) :: r1, r2, r12, grad1_Jw, grad1_psi0, grad2_psi0
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12=r1-r2
    r12l=sqrt(sum(r12**2))
    denom=1.0_dp+sys%b*r12l
    grad1_Jw=r12/r12l/2.0_dp/denom**2
    grad1_psi0=-sys%Z1*r1/r1l
    grad2_psi0=-sys%Z1*r2/r2l
    wlkr%EL=-sys%Z1**2-(sys%Znuc-sys%Z1)*(1.0_dp/r1l+1.0_dp/r2l) &
         & +sys%b/denom**3+sys%b*(2.0_dp+sys%b*r12l)/denom**2 &
         & -sum(grad1_Jw**2)-sum(grad1_psi0*grad1_Jw)+sum(grad2_psi0*grad1_Jw)
    wlkr%psiTsq=exp(-2.0_dp*sys%Z1*(r1l+r2l)+r12l/denom)
    wlkr%vD(1:3)=grad1_psi0+grad1_Jw
    wlkr%vD(4:6)=grad2_psi0-grad1_Jw
    wlkr%sgn=1
    ! }}}
  end subroutine EL_drift_Zstar_ee_cusp

  pure subroutine psiT2_Zstar_ee_cusp(sys,wlkr)
    ! {{{ trial wave function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l, denom
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12l=sqrt(sum((r1-r2)**2))
    denom=1.0_dp+sys%b*r12l
    wlkr%psiTsq=exp(-2.0_dp*sys%Z1*(r1l+r2l)+r12l/denom)
    ! }}}
  end subroutine psiT2_Zstar_ee_cusp

end module qmc_input


! Local variables:
! folded-file: t
! End:
