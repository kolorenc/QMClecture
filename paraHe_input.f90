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

! Hartree-Fock limit: -2.8617
! [Donald Allan McQuarrie, John Douglas Simon: "Physical chemistry: a molecular
! approach", p. 283]

! Exact energy: -2.903724
! [C. Schwartz, Int. J. Mod. Phys. E 15, 877 (2006)] also available as
! [arXiv:physics/0208004v1]

module qmc_input
  use types_const, only: dp, finish_line
  implicit none
  private

  public :: sys_dim
  public :: t_walker, t_sys
  public :: init_sys
  public :: EL_drift, psiT2

  ! logically, this should go to t_sys, but I really want to keep static arrays
  ! in t_walker, which enforces 'sys_dim' to be a global constant
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
     character(36) :: id="parahelium (Hylleraas wave function)"

     ! nuclear charge
     real(dp) :: Znuc=2.0_dp

     ! parameter in the Hylleraas ansatz
     ! a=4*(-123371.0_dp+sqrt(32098719641.0_dp))/2539875.0_dp for Znuc=2
     ! a=2*(-685.0_dp+3*sqrt(153862015.0_dp))/800167.0_dp for Znuc=1
     real(dp) :: a=0.08786283656089399052453400224163986218292_dp !Znuc=2
  end type t_sys

contains

  subroutine init_sys(sys,words)
    ! {{{ initialize the system and wave-function parameters
    use cfparser, only: t_words, readsection, readvalue, clear
    type(t_sys), intent(inout) :: sys
    type(t_words), intent(in) :: words
    type(t_words) :: syswords
    logical :: Znuc_given, a_given
    real(dp) :: Znuc, Etot

    ! read the parameters of the system and the trial wave function from the
    ! configuration file
    Znuc_given=.false.
    a_given=.false.
    if ( readsection(words,syswords,"system") ) then
       Znuc_given=readvalue(syswords,sys%Znuc,"Znuc")
       a_given=readvalue(syswords,sys%a,"a")
    end if
    call clear(syswords)

    ! set helium (Z=2) if the nuclear charge is not given
    if ( .not.Znuc_given ) sys%Znuc=2.0_dp

    ! if the variational parameter a is not given then use the optimal value
    Znuc=sys%Znuc
    if ( .not.a_given ) then
       !sys%Z1=Znuc-5.0_dp/16.0_dp
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
          (3.0_dp*(-98910 + Znuc*(-43776 +           &
                    Znuc*(304160 + Znuc*(413413 + 4096*Znuc*(47 + 8*Znuc))))) &
          )
    end if

    ! analytically calculated total energy
    Etot=(Znuc*(-8*sys%a*Znuc**2*(-567 + 32*Znuc*(27 + 8*Znuc*(7 + 3*Znuc))) - &
            8*Znuc**4*(-35 + 4*Znuc*(4 + Znuc*(25 + 16*Znuc))) -               &
            3*sys%a**2*(-10065 + 4*Znuc*(4992 + Znuc*(6903 + 2176*Znuc)))))/   &
         (4.0_dp*(8*Znuc**4*(24 + Znuc*(35 + 16*Znuc)) +     &
            24*sys%a*Znuc**2*(168 + Znuc*(189 + 64*Znuc)) +  &
            9*sys%a**2*(3680 + Znuc*(3355 + 896*Znuc))))

    ! output the settings to STDOUT
    write(unit=*,fmt='(1x,a)') "system and trial wave function:"
    write(unit=*,fmt='(3x,a)') sys%id

    write(unit=*,fmt='(3x,a,f12.8)',advance="no") &
            "nuclear charge                 [Znuc]:", sys%Znuc
    call finish_line(Znuc_given,"(default)")

    write(unit=*,fmt='(3x,a,f12.8)',advance="no") &
            "parameter in the wave function [a]   :", sys%a
    call finish_line(a_given,"(the optimal value)")

    write(unit=*,fmt='(3x,a,f12.8)') &
            "analytically computed total energy   :", Etot

    write(unit=*,fmt=*)
    ! }}}
  end subroutine init_sys


  ! ===========================================================================
  ! wave function of the Hylleraas form, cusps fullfiled; there is only a
  ! single variational parameter and hence the energy is no miracle, but it is
  ! possible to analytically optimize this single parameter from start to
  ! finish (although it is tedious for sure);
  ! see [Hylleraas, Z. Phys. 54, 347 (1929)] but he does not care about cusps
  ! (English translation is in the Hettema's book)
  ! ===========================================================================

  pure subroutine EL_drift(sys,wlkr)
    ! {{{ local energy, drift velocity, wave function squared
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
  end subroutine EL_drift

  pure subroutine psiT2(sys,wlkr)
    ! {{{ trial wave function squared
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
  end subroutine PsiT2

end module qmc_input


! Local variables:
! folded-file: t
! End:
