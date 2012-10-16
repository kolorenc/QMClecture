! Random number generators from Numerical Recipes, minimum compilation
! of needed subroutines
module nr_random
  use nrtype
  implicit none

  integer, parameter :: K4B=selected_int_kind(9)
  integer(K4B), parameter :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
  integer(K4B), save :: lenran=0, seq=0
  integer(K4B), save :: iran0,jran0,kran0,nran0,mran0,rans
  integer(K4B), dimension(:,:), pointer, save :: ranseeds
  integer(K4B), dimension(:), pointer, save :: iran,jran,kran, &
       nran,mran,ranv
  real(SP), save :: amm

  interface ran_hash
     module procedure ran_hash_s, ran_hash_v
  end interface

!!$  interface array_copy
!!$     module procedure array_copy_r, array_copy_d, array_copy_i
!!$  end interface

  interface ran1
     module procedure ran1_s, ran1_v
  end interface

  interface gasdev
     module procedure gasdev_s, gasdev_v
  end interface

contains

  
  subroutine ran1_s(harvest)
    ! {{{ .
    implicit none
    real(SP), intent(OUT) :: harvest
    if (lenran < 1) call ran_init(1)
    rans=iran0-kran0
    if (rans < 0) rans=rans+2147483579_k4b
    iran0=jran0
    jran0=kran0
    kran0=rans
    nran0=ieor(nran0,ishft(nran0,13))
    nran0=ieor(nran0,ishft(nran0,-17))
    nran0=ieor(nran0,ishft(nran0,5))
    if (nran0 == 1) nran0=270369_k4b
    mran0=ieor(mran0,ishft(mran0,5))
    mran0=ieor(mran0,ishft(mran0,-13))
    mran0=ieor(mran0,ishft(mran0,6))
    rans=ieor(nran0,rans)+mran0
    harvest=amm*merge(rans,not(rans), rans<0 )
    ! }}}
  end subroutine ran1_s

  subroutine ran1_v(harvest)
    ! {{{ .
    implicit none
    real(SP), dimension(:), intent(OUT) :: harvest
    integer(K4B) :: n
    n=size(harvest)
    if (lenran < n+1) call ran_init(n+1)
    ranv(1:n)=iran(1:n)-kran(1:n)
    where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
    iran(1:n)=jran(1:n)
    jran(1:n)=kran(1:n)
    kran(1:n)=ranv(1:n)
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
    nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
    where (nran(1:n) == 1) nran(1:n)=270369_k4b
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
    mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
    ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
    harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
    ! }}}
  end subroutine ran1_v

  subroutine gasdev_s(harvest)
    ! {{{ Returns in harvest a normally distributed deviate with zero mean
    !   and unit variance, using ran1 as the source of uniform deviates.
    real(sp), intent(out) :: harvest
    real(sp) :: rsq,v1,v2
    real(sp), save :: g
    logical, save :: gaus_stored=.false.
    if (gaus_stored) then
       harvest=g
       gaus_stored=.false.
    else
       do
          call ran1(v1)
          call ran1(v2)
          !call random_number(v1)
          !call random_number(v2)
          v1=2.0_sp*v1-1.0_sp
          v2=2.0_sp*v2-1.0_sp
          rsq=v1**2+v2**2
          if (rsq > 0.0_sp .and. rsq < 1.0_sp) exit
       end do
       rsq=sqrt(-2.0_sp*log(rsq)/rsq)
       harvest=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
    end if
    ! }}}
  end subroutine gasdev_s

  subroutine gasdev_v(harvest)
    ! {{{ vector version of gasdev
    use nrutil, only : array_copy
    real(SP), dimension(:), intent(OUT) :: harvest
    real(SP), dimension(size(harvest)) :: rsq,v1,v2
    real(SP), allocatable, dimension(:), save :: g
    integer(I4B) :: n,ng,nn,m
    integer(I4B), save :: last_allocated=0
    logical, save :: gaus_stored=.false.
    logical, dimension(size(harvest)) :: mask
    n=size(harvest)
    if (n /= last_allocated) then
       if (last_allocated /= 0) deallocate(g)
       allocate(g(n))
       last_allocated=n
       gaus_stored=.false.
    end if
    if (gaus_stored) then
       harvest=g
       gaus_stored=.false.
    else
       ng=1
       do
          if (ng > n) exit
          call ran1(v1(ng:n))
          call ran1(v2(ng:n))
          !call random_number(v1(ng:n))
          !call random_number(v2(ng:n))
          v1(ng:n)=2.0_sp*v1(ng:n)-1.0_sp
          v2(ng:n)=2.0_sp*v2(ng:n)-1.0_sp
          rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
          mask(ng:n)=(rsq(ng:n)>0.0 .and. rsq(ng:n)<1.0)
          call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
          v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
          rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
          ng=ng+nn
       end do
       rsq=sqrt(-2.0_sp*log(rsq)/rsq)
       harvest=v1*rsq
       g=v2*rsq
       gaus_stored=.true.
    end if
    ! }}}
  end subroutine gasdev_v

  ! {{{ module ran_state

  subroutine ran_init(length)
    ! {{{ .
    use nrtype; use nrutil, only : arth,nrerror,reallocate
    implicit none
    integer(K4B), intent(IN) :: length
    integer(K4B) :: new,j,hgt
    if (length < lenran) return
    hgt=hg
    if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
    if (hgng >= 0)        call nrerror('ran_init: arith assump 2 fails')
    if (hgt+1 /= hgng)    call nrerror('ran_init: arith assump 3 fails')
    if (not(hg) >= 0)     call nrerror('ran_init: arith assump 4 fails')
    if (not(hgng) < 0)    call nrerror('ran_init: arith assump 5 fails')
    if (hg+hgng >= 0)     call nrerror('ran_init: arith assump 6 fails')
    if (not(-1_k4b) < 0)  call nrerror('ran_init: arith assump 7 fails')
    if (not(0_k4b) >= 0)  call nrerror('ran_init: arith assump 8 fails')
    if (not(1_k4b) >= 0)  call nrerror('ran_init: arith assump 9 fails')
    if (lenran > 0) then
       ranseeds=>reallocate(ranseeds,length,5)
       ranv=>reallocate(ranv,length-1)
       new=lenran+1
    else
       allocate(ranseeds(length,5))
       allocate(ranv(length-1))
       new=1
       amm=nearest(1.0_sp,-1.0_sp)/hgng
       if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            call nrerror('ran_init: arth assump 10 fails')
    end if
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
    do j=1,4
       call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    end do
    where (ranseeds(new:,1:3) < 0) &
         ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
    where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1
    if (new == 1) then
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
       rans=nran0
    end if
    if (length > 1) then
       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran
    end if
    lenran=length
    ! }}}
  end subroutine ran_init

  subroutine ran_deallocate
    ! {{{ .
    if (lenran > 0) then
       deallocate(ranseeds,ranv)
       nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
       lenran = 0
    end if
    ! }}}
  end subroutine ran_deallocate

  subroutine ran_seed(sequence,size,put,get)
    ! {{{ .
    implicit none
    integer, optional, intent(IN) :: sequence
    integer, optional, intent(OUT) :: size
    integer, dimension(:), optional, intent(IN) :: put
    integer, dimension(:), optional, intent(OUT) :: get
    if (present(size)) then
       size=5*lenran
    else if (present(put)) then
       if (lenran == 0) return
       ranseeds=reshape(put,shape(ranseeds))
       where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
          where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
          iran0=ranseeds(1,1)
          jran0=ranseeds(1,2)
          kran0=ranseeds(1,3)
          mran0=ranseeds(1,4)
          nran0=ranseeds(1,5)
    else if (present(get)) then
          if (lenran == 0) return
          ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
          get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
          call ran_deallocate
          seq=sequence
    end if
    ! }}}
  end subroutine ran_seed

  subroutine ran_hash_s(il,ir)
    ! {{{ .
    implicit none
    integer(K4B), intent(INOUT) :: il,ir
    integer(K4B) :: is,j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
    ! }}}
  end subroutine ran_hash_s

  subroutine ran_hash_v(il,ir)
    ! {{{ .
    implicit none
    integer(K4B), dimension(:), intent(INOUT) :: il,ir
    integer(K4B), dimension(size(il)) :: is
    integer(K4B) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
    ! }}}
  end subroutine ran_hash_v

  ! }}}
  
!!$  ! {{{ array_copy definitions
!!$
!!$  subroutine array_copy_r(src,dest,n_copied,n_not_copied)
!!$    real(SP), dimension(:), intent(IN) :: src
!!$    real(SP), dimension(:), intent(OUT) :: dest
!!$    integer(I4B), intent(OUT) :: n_copied, n_not_copied
!!$    n_copied=min(size(src),size(dest))
!!$    n_not_copied=size(src)-n_copied
!!$    dest(1:n_copied)=src(1:n_copied)
!!$  end subroutine array_copy_r
!!$
!!$  subroutine array_copy_d(src,dest,n_copied,n_not_copied)
!!$    real(DP), dimension(:), intent(IN) :: src
!!$    real(DP), dimension(:), intent(OUT) :: dest
!!$    integer(I4B), intent(OUT) :: n_copied, n_not_copied
!!$    n_copied=min(size(src),size(dest))
!!$    n_not_copied=size(src)-n_copied
!!$    dest(1:n_copied)=src(1:n_copied)
!!$  end subroutine array_copy_d
!!$
!!$  subroutine array_copy_i(src,dest,n_copied,n_not_copied)
!!$    integer(I4B), dimension(:), intent(IN) :: src
!!$    integer(I4B), dimension(:), intent(OUT) :: dest
!!$    integer(I4B), intent(OUT) :: n_copied, n_not_copied
!!$    n_copied=min(size(src),size(dest))
!!$    n_not_copied=size(src)-n_copied
!!$    dest(1:n_copied)=src(1:n_copied)
!!$  end subroutine array_copy_i
!!$
!!$  ! }}}

end module nr_random

! Local variables:
! folded-file: t
