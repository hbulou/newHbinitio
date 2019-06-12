module pseudopotential
  use global
  implicit none
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
    open(unit=1,file=trim(pp%file),action="read",form='formatted',status='old')
    read(1, *) 
    read(1, *)
    read(1, *) 
    read(1,*) pp%lmax, pp%npotu, pp%n, pp%b, pp%a, pp%zval
    print *, 'PSEUDO> lmax= ',pp%lmax, pp%npotu, pp%n, pp%b, pp%a, 'zval= ',pp%zval
    !
    ! read the grid
    !
    allocate(pp%r(0:pp%n))
    read(1,*)
    read(1,*) (pp%r(i),i=1,pp%n)
    !
    ! read the pseudopotentials
    !
    allocate(pp%pot(0:pp%lmax-1,0:pp%n))
    do i=0,pp%lmax-1
       print *,"PSEUDO > l=",i
       read(1,*) ;        read(1,*)
       read(1,*) (pp%pot(i,j),j=1,pp%n)
    end do
    allocate(pp%zv(0:pp%n))
    read(1,*) ;     read(1,*) (pp%zv(j),j=1,pp%n)
    read(1,*) ;     read(1,*) (pp%zv(j),j=1,pp%n)


    close(1)

    do i=0,pp%lmax-1
       pp%pot(i,:)=pp%pot(i,:)/pp%r(:)
    end do
    open(unit=1,file='pot.dat',form='formatted',status='unknown')
    do i=1,pp%n
       write(1,*) pp%r(i),(pp%pot(j,i),j=0,pp%lmax-1),-2*pp%zval/pp%r(i)
    end do
    close(1)
     


  end subroutine read_pp

end module pseudopotential
