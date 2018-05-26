! ==========================================================================
!
! Basic VMC/DMC code for educational purposes.
! This file: estimation of error bars from correlated data.
!
! Copyright (C) 2012 Jindrich Kolorenc
!
! The software is released under MIT/X11 license.
!
! ==========================================================================

! The error bars are estimated with the blocking method [Flyvbjerg &
! Petersen, J. Chem. Phys. 91, 461-6 (1989)].
! The subroutines assume that we start with 2^N data points (the reblocking
! is done by pairing the points/blocks).

module reblocking
  use types_const, only: dp
  implicit none
  private

  public :: errorbar

contains

  subroutine reblock(A)
    ! {{{ reblocking by pairing
    real(dp), dimension(:), intent(inout) :: A
    integer :: i
    do i=1, size(A), 2
       A(i/2+1)=(A(i)+A(i+1))*0.5_dp
    end do
    ! }}}
  end subroutine reblock

  function varmean(A,mean) result(s)
    ! {{{ variance of the mean
    real(dp), dimension(:), intent(in) :: A
    real(dp), intent(in) :: mean
    real(dp) :: s
    integer :: dim
    dim=size(A)
    s=sum((A-mean)**2)/(dim-1)/dim
    ! }}}
  end function varmean

  function errorbar(A,mean,filename,corrlen,finalNblock) result(s)
    ! {{{ NB: on output, the array A is modified, since reblocking is
    !     called directly on it
    real(dp), dimension(:), intent(inout) :: A
    real(dp), intent(in) :: mean
    character(*), intent(in), optional :: filename
    real(dp), intent(out), optional :: corrlen
    integer, intent(out), optional :: finalNblock
    real(dp) :: s0, s, s1
    integer :: dim
    logical :: writeout
    writeout=present(filename)
    if ( writeout ) then
       open(unit=99,file=filename,action="write")
       write(unit=99,fmt=*) "# number of blocks    variance    errorbar"
    end if
    dim=size(A)
    s0=varmean(A,mean)
    s=s0
    do while ( dim > 2 )
       call reblock(A(1:dim))
       dim=dim/2
       s1=varmean(A(1:dim),mean)
       if ( writeout ) write(unit=99,fmt=*) dim, s1, sqrt(s1)
       if ( s1 <= s ) exit
       s=s1
    end do
    if ( writeout ) close(unit=99)
    if ( present(corrlen) ) corrlen=s/s0
    if ( present(finalNBlock) ) finalNblock=dim
    s=sqrt(s)
    ! }}}
  end function errorbar

end module reblocking


! Local variables:
! folded-file: t
! End:
