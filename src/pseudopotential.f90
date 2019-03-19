module pseudopotential
  implicit none
  type t_pseudo
     integer::n
     integer :: npot,npotu
     double precision :: a,b,zval
     double precision,allocatable::r(:)
     double precision,allocatable::pot(:,:)
  end type t_pseudo
  type(t_pseudo) :: pp
contains
    ! --------------------------------------------------------------------------------------
  !
  !             read_pp()
  !
  ! read pseudpopotential file
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_pp(pp)
    type(t_pseudo)::pp
    integer ::i,j
    open(unit=1,file='/home/bulou/ownCloud/src/Octopus/octopus-7.1/share/pseudopotentials/PSF/H.psf',&
         action="read",form='formatted',status='old')
    read(1, *) 
    read(1, *)
    read(1, *) 
    read(1,*) pp%npot, pp%npotu, pp%n, pp%b, pp%a, pp%zval
    print *, pp%npot, pp%npotu, pp%n, pp%b, pp%a, pp%zval
    ! add extra point for the zero
    pp%n=pp%n+1
    print *,pp%n
    allocate(pp%r(pp%n))
    allocate(pp%pot(pp%npot,pp%n))
    pp%r(1)=0.0
    read(1, *)
    read(1,*) (pp%r(i),i=2,pp%n)
    do i=1,pp%npot
       read(1,*)
       pp%pot(i,1)=0.0
       read(1,*) (pp%pot(i,j),j=2,pp%n)
    end do
    ! do i=2,pp%n
    !    read(1,*) pp%r(i)
    !    print *,i,pp%r(i)
    !    j=j+1
    !    if(j.eq.4) then
    !       read(1,*)
    !       j=0
    !    end if
    ! end do

    close(1)

    do i=1,pp%npot
       pp%pot(i,2:pp%n)=pp%pot(i,2:pp%n)/pp%r(2:pp%n)
    end do
    open(unit=1,file='pot.dat',form='formatted',status='unknown')
    do i=1,pp%n
       write(1,*) pp%r(i),(pp%pot(j,i),j=1,pp%npot)
    end do
    close(1)



  end subroutine read_pp

end module pseudopotential
