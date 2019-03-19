module conjugate_gradient_mod
  use global
  use mesh_mod
  use IO
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !             Conjugate_gradient()
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Conjugate_gradient(b,x,n,h)
    implicit none
    double precision::b(:),x(:)
    integer::n
    double precision::h

    double precision,allocatable::g(:)
    double precision,allocatable::d(:)
    double precision,allocatable::Ad(:)
    double precision,allocatable::Ax(:)
    double precision::alpha,beta,dAd,f,fold,cvg
    integer::i,iloop,iloopmax
    double precision, external :: ddot
    allocate(g(n))
    allocate(d(n))
    allocate(Ad(n))
    allocate(Ax(n))


    
    
    g(1)=(x(2)-2*x(1))/h**2+b(1)
    do i=2,n-1
       g(i)=(x(i+1)-2*x(i)+x(i-1))/h**2+b(i)
    end do
    g(n)=(x(n-1)-2*x(n))/h**2+b(n)


    call dcopy(n,g,1,d,1) ! g->d

!    alpha=-1.0 ;    call dscal(n,alpha,d,1) ! d-> -d
!    print *,ddot(n,d,1,d,1)
!    stop

    cvg=1.0e-8
    fold=0.0
    f=2*cvg
    iloopmax=10000
    iloop=1
    do while(.not.((iloop.le.iloopmax).and.(abs(f-fold).lt.cvg))) 
!       do iloop=1,iloopmax

       Ad(1)=(d(2)-2*d(1))/h**2
       do i=2,n-1
          Ad(i)=(d(i+1)-2*d(i)+d(i-1))/h**2
       end do
       Ad(n)=(d(n-1)-2*d(n))/h**2

       dAd=ddot(n,d,1,Ad,1)
       alpha=-ddot(n,g,1,d,1)/dAd
       
       call daxpy(n,alpha,d,1,x,1)  ! x+alpha*d -> x
       
       call daxpy(n,alpha,Ad,1,g,1) ! g+alpha*Ad -> g 

       fold=f
       f=ddot(n,b,1,x,1)
       Ax(1)=(x(2)-2*x(1))/h**2
       do i=2,n-1
          Ax(i)=(x(i+1)-2*x(i)+x(i-1))/h**2
       end do
       Ax(n)=(x(n-1)-2*x(n))/h**2
       f=f+0.5*ddot(n,x,1,Ax,1)
!       print *,iloop,f,abs(f-fold),ddot(n,g,1,g,1)

       beta=ddot(n,g,1,Ad,1)/dAd
       call dscal(n,beta,d,1) ! beta*d-> d
       alpha=-1.0
       call daxpy(n,alpha,g,1,d,1) ! -g+d -> d 
       iloop=iloop+1
       if(iloop.eq.iloopmax) then
          print *,"ERROR in Conjugate_gradient"
          stop
       end if
    end do
    
    deallocate(g)
    deallocate(d)
    deallocate(Ad)
    deallocate(Ax)
  end subroutine Conjugate_gradient
  ! --------------------------------------------------------------------------------------
  !
  !             Conjugate_gradient_3D()
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Conjugate_gradient_3D(b,x,n,h,m)
    implicit none
    double precision::b(:),x(:)
    integer::n
    double precision::h,hsqr
    type(t_mesh)::m

    double precision,allocatable::g(:)
    double precision,allocatable::d(:)
    double precision,allocatable::Ad(:)
    double precision,allocatable::Ax(:)
    double precision::alpha,beta,dAd,f,fold,cvg
    integer::i,iloop,iloopmax,l
    character (len=1024) :: filecube

    double precision, external :: ddot
    allocate(g(n))
    allocate(d(n))
    allocate(Ad(n))
    allocate(Ax(n))
    hsqr=h**2
    
    do i=1,n
       g(i)=-2*m%dim*x(i)/hsqr
       do l=1,m%node(i)%n_neighbors
          g(i)=g(i)+x(m%node(i)%list_neighbors(l))/hsqr
       end do
       g(i)=g(i)+b(i)
    end do

    call dcopy(n,g,1,d,1) ! g->d

!    alpha=-1.0 ;    call dscal(n,alpha,d,1) ! d-> -d
!    print *,ddot(n,d,1,d,1)
!    stop

    cvg=1.0e-8
    fold=0.0
    f=2*cvg
    iloopmax=10000
    iloop=1
    do while(.not.((iloop.le.iloopmax).and.(abs(f-fold).lt.cvg))) 
!       do iloop=1,iloopmax

       do i=1,n
          Ad(i)=-2*m%dim*d(i)/hsqr
          do l=1,m%node(i)%n_neighbors
             Ad(i)=Ad(i)+d(m%node(i)%list_neighbors(l))/hsqr
          end do
        end do

       dAd=ddot(n,d,1,Ad,1)
       alpha=-ddot(n,g,1,d,1)/dAd
       
       call daxpy(n,alpha,d,1,x,1)  ! x+alpha*d -> x
       call daxpy(n,alpha,Ad,1,g,1) ! g+alpha*Ad -> g 

       fold=f
       f=ddot(n,b,1,x,1)
       !print *,f-fold
       Ax(1)=(x(2)-2*x(1))/h**2
       do i=2,n-1
          Ax(i)=(x(i+1)-2*x(i)+x(i-1))/h**2
       end do
       Ax(n)=(x(n-1)-2*x(n))/h**2
       f=f+0.5*ddot(n,x,1,Ax,1)
!       print *,iloop,f,abs(f-fold),ddot(n,g,1,g,1)

       beta=ddot(n,g,1,Ad,1)/dAd
       call dscal(n,beta,d,1) ! beta*d-> d
       alpha=-1.0
       call daxpy(n,alpha,g,1,d,1) ! -g+d -> d 
       iloop=iloop+1
       if(iloop.eq.iloopmax) then
          print *,"ERROR in COnjugtae_gradient"
          stop
       end if
    end do

!    write(filecube,'(a)') 'hartee.cube'
!    call save_cube_3D(x,filecube,m)
    
    deallocate(g)
    deallocate(d)
    deallocate(Ad)
    deallocate(Ax)
  end subroutine Conjugate_gradient_3D
end module conjugate_gradient_mod
