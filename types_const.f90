module types_const
  implicit none

  integer, parameter :: sp=kind(1.0)
  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: spc=kind((1.0,1.0))
  integer, parameter :: dpc=kind((1.0d0,1.0d0))

  integer, parameter :: i2b=selected_int_kind(4)  ! 2 byte int < 2**15
  integer, parameter :: i4b=selected_int_kind(9)  ! 4 byte int < 2**31
  integer, parameter :: i8b=selected_int_kind(10) ! 8 byte int < 2**63

  integer, parameter :: lgt=kind(.true.)

  complex(dpc), parameter :: ii=(0.0_dp,1.0_dp)

  real(dp), parameter :: pi_d=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: twopi_d=6.283185307179586476925286766559005768394_dp
  real(dp), parameter :: sqrt2_d=1.414213562373095048801688724209698078570_dp
  real(dp), parameter :: sqrt3_d=1.732050807568877293527446341505872366943_dp
  real(dp), parameter :: sqrt5_d=2.236067977499789696409173668731276235441_dp
  real(sp), parameter :: pi=3.141592653589793238462643383279502884197_sp
  real(sp), parameter :: twopi=6.283185307179586476925286766559005768394_sp
  real(sp), parameter :: sqrt2=1.414213562373095048801688724209698078570_sp
  real(dp), parameter :: sqrt3=1.732050807568877293527446341505872366943_sp
  real(dp), parameter :: sqrt5=2.236067977499789696409173668731276235441_sp

  type one_vec
     real(dp), dimension(:), allocatable :: vec
  end type one_vec

  type one_mtrx
     real(dp), dimension(:,:), allocatable :: mtrx
  end type one_mtrx

end module types_const
