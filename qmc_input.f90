! problem specific subroutines for He atom (parahelium=opposite spins), good
! also for any helium-like ion, for instance H^-

! some of the VMC energies we know from other sources
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
  use types_const, only: dp
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
     real(dp) :: Znuc=2.0_dp            ! nuclear charge

     ! parameter in the Hylleraas ansatz
     ! a=4*(-123371.0_dp+sqrt(32098719641.0_dp))/2539875.0_dp for Znuc=2
     ! a=2*(-685.0_dp+3*sqrt(153862015.0_dp))/800167.0_dp for Znuc=1
     real(dp) :: a=0.08786283656089399052453400224163986218292_dp

     ! parameters in the Slater-Jastrow ansatz
     real(dp) :: Z1=27.0_dp/16.0_dp     ! effective charge in the s wave
     real(dp) :: b=0.5_dp               ! denominator in the ee-cusp term
  end type t_sys

  interface EL_drift
     module procedure EL_drift_all_cusps
  end interface EL_drift

  interface psiT2
     module procedure psiT2_all_cusps
  end interface psiT2

contains

  subroutine init_sys(sys,Znuc)
    ! {{{ initialize the system and wave-function parameters
    type(t_sys), intent(inout) :: sys
    real(dp), intent(in) :: Znuc
    sys%Znuc=Znuc
    sys%Z1=Znuc-5.0_dp/16.0_dp
    sys%a=(2*(Znuc**2*(18105 -                &
              2*Znuc*(-18912 +                &
                 Znuc*(-1848 + Znuc*(16329 + 8*Znuc*(1301 + 256*Znuc))) &
                 )) + sqrt(Znuc**4*           &
             (128388465 +                     &
               2*Znuc*(327330432 +            &
                  Znuc*(667667520 +           &
                     Znuc*(730902672 +        &
                        Znuc*                 &
                         (733280304 +         &
                          Znuc*               &
                          (1064874760 +       &
                          Znuc*               &
                          (1264584679 +       &
                          16*Znuc*            &
                          (56671967 +         &
                          8*Znuc*             &
                          (2950393 + 512*Znuc*(1309 + 128*Znuc))) &
                          ))))))))))/         &
       (3.*(-98910 + Znuc*(-43776 +           &
              Znuc*(304160 + Znuc*(413413 + 4096*Znuc*(47 + 8*Znuc))))) &
         )
    ! }}}
  end subroutine init_sys

  ! ===========================================================================
  ! simple trial wave function, just a product of two s orbitals with effective
  ! charge; all cusps are broken when the effective charge deviates from 2;
  ! the earliest reference I could find is [Kellner, Z. Phys. 44, 91 (1927)]
  ! ===========================================================================

  subroutine EL_drift_no_cusp(sys,wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
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
  end subroutine EL_drift_no_cusp

  subroutine psiT2_no_cusp(sys,wlkr)
    ! {{{ trial wave-function squared
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
  end subroutine psiT2_no_cusp


  ! ===========================================================================
  ! the above function multiplied by an e-e cusp-fixing Jastrow factor,
  ! see for instance [Reynolds, Ceperley, Alder & Lester, J. Chem. Phys. 77,
  ! 5593 (1982)] but this functional form is surely older
  ! ===========================================================================

  subroutine EL_drift_ee_cusp(sys,wlkr)
    ! {{{ local energy, drift velocity, wave-function squared, implemented
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
  end subroutine EL_drift_ee_cusp

  subroutine psiT2_ee_cusp(sys,wlkr)
    ! {{{ trial wave-function squared
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
  end subroutine psiT2_ee_cusp


  ! ===========================================================================
  ! wave function of the Hylleraas form, cusps fullfiled; there is only a
  ! single variational parameter and hence the energy is no miracle, but it is
  ! possible to analytically optimize this single parameter from start to
  ! finish (although it is tedious for sure);
  ! see [Hylleraas, Z. Phys. 54, 347 (1929)] but he does not care about cusps
  ! (English translation is in the Hettema's book)
  ! ===========================================================================

  subroutine EL_drift_all_cusps(sys,wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp), dimension(3) :: r1, r2, r12
    real(dp) :: u, v, w, u2, v2, w2, v2_w2, vw, twouvw, vDu, vDv, vDw, Xa, u_2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    v2=sum(r1**2)
    w2=sum(r2**2)
    r12=r1-r2
    u2=sum(r12**2)
    u=sqrt(u2)
    v=sqrt(v2)
    w=sqrt(w2)
    v2_w2=v2+w2
    vw=v*w
    twouvw=2.0_dp*u*vw
    Xa=1.0_dp + sys%a*v2_w2
    u_2=u+2.0_dp
    wlkr%PsiTsq=(u_2**2*Xa**2)/(4.0_dp*exp(2.0_dp*(v + w)*sys%Znuc))
    wlkr%EL=(twouvw*(Xa + sys%a*(-12.0_dp - 8.0_dp*u)) + &
         (v + w)*(4.0_dp*sys%a*twouvw - (v - w)**2*Xa + &
            u2*(Xa + sys%a*4.0_dp*vw))*sys%Znuc)/(u_2*twouvw*Xa) - &
         sys%Znuc**2
    vDu=1.0_dp/(u_2*u)
    vDv=2.0_dp*sys%a/Xa
    vDw=vDv-sys%Znuc/w
    vDv=vDv-sys%Znuc/v
    wlkr%vD(1:3)=vDu*r12
    wlkr%vD(4:6)=vDw*r2-wlkr%vD(1:3)
    wlkr%vD(1:3)=vDv*r1+wlkr%vD(1:3)
    ! }}}
  end subroutine EL_drift_all_cusps

  subroutine psiT2_all_cusps(sys,wlkr)
    ! {{{ trial wave-function squared
    type(t_sys), intent(in) :: sys
    type(t_walker), intent(inout) :: wlkr
    real(dp), dimension(3) :: r1, r2, r12
    real(dp) :: u, v, w, u2, v2, w2, v2_w2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    v2=sum(r1**2)
    w2=sum(r2**2)
    r12=r1-r2
    u2=sum(r12**2)
    u=sqrt(u2)
    v=sqrt(v2)
    w=sqrt(w2)
    v2_w2=v2+w2
    wlkr%PsiTsq=((2.0_dp + u)**2*(1.0_dp + sys%a*v2_w2)**2)/(4.0_dp*exp(2.0_dp*(v + w)*sys%Znuc))
    ! }}}
  end subroutine PsiT2_all_cusps

end module qmc_input


! Local variables:
! folded-file: t
! End:
