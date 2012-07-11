! problem specific subroutines for He atom (parahelium=opposite spins), good
! also for any helium-like ion, for instance H^-

! some of the VMC energies we know from other sources
!  "no cusp", Z1=2            -2.75
!  "no cusp", Z1=2-5/16       -2.84766
!  "all cusps", a=0.0878628   -2.89537

module qmc_input
  use types_const, only: dp
  implicit none
  private

  public :: sys_dim
  public :: t_walker
  public :: EL_drift, psiT2
  
  integer, parameter :: sys_dim=6       ! dimension of walker

  type t_walker
     real(dp), dimension(sys_dim) :: r  ! position
     real(dp), dimension(sys_dim) :: vD ! drift velocity
     real(dp) :: EL                     ! local energy  
     real(dp) :: psiTsq                 ! square of trial wavefunction
     integer :: sgn                     ! sign of trial wavefunction at r
     integer :: age                     ! how many times rejected by det. bal.
  end type t_walker

  real(dp), parameter :: ZHe=2.0_dp
  real(dp), parameter :: Z1=ZHe-5.0_dp/16.0_dp  ! effective charge for 1s wave
  real(dp), parameter :: b_cusp=0.5_dp
  real(dp), parameter :: a=4*(-123371.0_dp+sqrt(32098719641.0_dp))/2539875.0_dp

  interface EL_drift
     module procedure EL_drift_all_cusps
  end interface EL_drift

  interface psiT2
     module procedure psiT2_all_cusps
  end interface psiT2

contains

  ! ===========================================================================
  ! simple trial wave function, just a product of two s orbitals with effective
  ! charge; all cusps are broken when the effective charge deviates from 2
  ! ===========================================================================

  subroutine EL_drift_no_cusp(wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12l=sqrt(sum((r1-r2)**2))
    wlkr%EL=-Z1**2-(ZHe-Z1)*(1.0_dp/r1l+1.0_dp/r2l)+1.0_dp/r12l
    wlkr%psiTsq=exp(-2.0_dp*Z1*(r1l+r2l))
    wlkr%vD(1:3)=-Z1*r1/r1l
    wlkr%vD(4:6)=-Z1*r2/r2l
    wlkr%sgn=1
    ! }}}
  end subroutine EL_drift_no_cusp

  subroutine psiT2_no_cusp(wlkr)
    ! {{{ trial wave-function squared
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    wlkr%psiTsq=exp(-2.0_dp*Z1*(r1l+r2l))
    ! }}}
  end subroutine psiT2_no_cusp


  ! ===========================================================================
  ! the above function multiplied by a cusp-fixing Jastrow factor
  ! ===========================================================================

  subroutine EL_drift_ee_cusp(wlkr)
    ! {{{ local energy, drift velocity, wave-function squared, implemented
    !     electron-electron cusp
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l, denom
    real(dp), dimension(3) :: r1, r2, r12, grad1_Jw, grad1_psi0, grad2_psi0
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12=r1-r2
    r12l=sqrt(sum(r12**2))
    denom=1.0_dp+b_cusp*r12l
    grad1_Jw=r12/r12l/2.0_dp/denom**2
    grad1_psi0=-Z1*r1/r1l
    grad2_psi0=-Z1*r2/r2l
    wlkr%EL=-Z1**2-(ZHe-Z1)*(1.0_dp/r1l+1.0_dp/r2l) &
         & +b_cusp/denom**3+b_cusp*(2.0_dp+b_cusp*r12l)/denom**2 &
         & -sum(grad1_Jw**2)-sum(grad1_psi0*grad1_Jw)+sum(grad2_psi0*grad1_Jw)
    wlkr%psiTsq=exp(-2.0_dp*Z1*(r1l+r2l)+r12l/denom)
    wlkr%vD(1:3)=grad1_psi0+grad1_Jw
    wlkr%vD(4:6)=grad2_psi0-grad1_Jw
    wlkr%sgn=1
    ! }}}
  end subroutine EL_drift_ee_cusp

  subroutine psiT2_ee_cusp(wlkr)
    ! {{{ trial wave-function squared
    type(t_walker), intent(inout) :: wlkr
    real(dp) :: r1l, r2l, r12l, denom
    real(dp), dimension(3) :: r1, r2
    r1=wlkr%r(1:3)
    r2=wlkr%r(4:6)
    r1l=sqrt(sum(r1**2))
    r2l=sqrt(sum(r2**2))
    r12l=sqrt(sum((r1-r2)**2))
    denom=1.0_dp+b_cusp*r12l
    wlkr%psiTsq=exp(-2.0_dp*Z1*(r1l+r2l)+r12l/denom)
    ! }}}
  end subroutine psiT2_ee_cusp


  ! ===========================================================================
  ! wave function of the Hylleraas form, cusps fullfiled; there is only a
  ! single variational parameter and hence the energy is no miracle, but it is
  ! possible to analytically optimize this single parameter from start to
  ! finish (although it is tedious for sure)
  ! ===========================================================================

  subroutine EL_drift_all_cusps(wlkr)
    ! {{{ local energy, drift velocity, wave-function squared
    type(t_walker), intent(inout) :: wlkr
    real(dp), dimension(3) :: r1, r2, r12
    real(dp) :: u, v, w, u2, v2, w2, v2_w2, vw, vDu, vDv, vDw
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
    wlkr%PsiTsq=((2 + u)**2*(1 + a*v2_w2)**2)/(4*exp(2*(v + w)*ZHe))
    wlkr%EL=(2*u*vw*(1 + a*(-12 - 8*u + v2_w2)) + &
         (v + w)*(8*a*u*vw - (v - w)**2*(1 + a*v2_w2) + &
            u**2*(1 + a*(v2_w2 + 4*v*w)))*ZHe - &
         2*u*(2 + u)*vw*(1 + a*v2_w2)*ZHe**2)/ &
       (2*u*(2 + u)*vw*(1 + a*v2_w2))
    vDu=1/(2+u)
    vDv=2*a*v/(1+a*v2_w2)-ZHe
    vDw=2*a*w/(1+a*v2_w2)-ZHe
    wlkr%vD(1:3)=vDv*r1/v+vDu*r12/u
    wlkr%vD(4:6)=vDw*r2/w-vDu*r12/u
    ! }}}
  end subroutine EL_drift_all_cusps

  subroutine psiT2_all_cusps(wlkr)
    ! {{{ trial wave-function squared
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
    wlkr%PsiTsq=((2 + u)**2*(1 + a*v2_w2)**2)/(4*exp(2*(v + w)*ZHe))
    ! }}}
  end subroutine PsiT2_all_cusps

end module qmc_input


! Local variables:
! folded-file: t
! End:
