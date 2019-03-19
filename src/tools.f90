module tools
  use global
  use mesh_mod
  use param_mod
 
  implicit none
  !------------------------------------------
  type t_GramSchmidt
     integer :: nindep
     integer :: ndep ! number of linear dependencies discovered
  end type t_GramSchmidt
contains
  ! --------------------------------------------------------------------------------------
  !
  !              GramSchmidt()
  !
  ! --------------------------------------------------------------------------------------
  ! 
  ! Ref.: D. G. Clayton "Gram-Schmidt Orthogonalization", J. Roy. Stat. Soc. C 20, 335 (1971)
  !
  subroutine GramSchmidt(Vout,Vin,nvec,m,GS)
    implicit none
    integer ::  nvec
    double precision,allocatable :: Vin(:,:),Vout(:,:)
    type(t_mesh)::m
    type(t_GramSchmidt) :: GS
    
    integer :: i,k,i0
    double precision, parameter :: ETA=1.0e-6
    double precision,allocatable :: a(:)
    double precision :: normloc
    double precision, external :: ddot

    allocate(a(nvec))
    print *,'GS> ----------------------'
    print *,"GS> Gram-Schmidt algorithm"
    print *,'GS> ----------------------'
    print *,'GS> ',nvec,' vectors to orthogonalize'
    print *,'GS> ',GS%nindep,' vectors are already orthogonalized'
    Vout(:,1:GS%nindep)=Vin(:,1:GS%nindep)
    GS%ndep=0 
    i0=GS%nindep
    do i=i0+1,nvec
       Vout(:,GS%nindep+1)=Vin(:,i)
       do k=1,GS%nindep
          ! We compute the projection of Vini(:,i) on V(:,1-nindep)
          a(k)=ddot(m%nactive,Vout(:,k),1,Vin(:,i),1)
          !print *,'GS > ',i,k,a(k)
          ! then we remove V(:,k) from V(:,nindep+1)
          call daxpy(m%nactive,-a(k),Vout(:,k),1,Vout(:,GS%nindep+1),1)
       end do
       ! now wre compute the norm of V(:,nindep+1)
       normloc=sqrt(ddot(m%nactive,Vout(:,GS%nindep+1),1,Vout(:,GS%nindep+1),1))
       !print *,'GS > norm(',i,')=',norm
       if (normloc.le.ETA) then
          GS%ndep=GS%ndep+1 ! V(:,nindep+1) is not linearly inependent
       else
          normloc=1.0/normloc
          call dscal(m%nactive,normloc,Vout(:,GS%nindep+1),1)
          !              do k=2,icur
          !                 print *,'<U',k-1,'|U',i,'>=',ddot(N,V(:,k-1),1,V(:,icur),1)
          !              end do
          GS%nindep=GS%nindep+1
       end if
    end do
    print *,'GS> ',GS%ndep,' vectors linearly dependant'
    print *,'GS> ',GS%nindep,' vectors linearly independant'
    print *,'GS> Size of the basis from ',nvec,' to ',GS%nindep
    !call check_ortho(Vout,nvec,m%nactive)
    !stop
    nvec=GS%nindep
    deallocate(a)
  end subroutine GramSchmidt
    ! --------------------------------------------------------------------------------------
  !
  !              simpson()
  !
  ! --------------------------------------------------------------------------------------
  function  simpson(m,f)
    implicit none
    type (t_mesh)::m
    double precision::simpson
    double precision :: f(:)
    integer::i
    simpson=0.0
    do i=1,m%nactive-2,2
       simpson=simpson+f(i)*f(i)+4*f(i+1)*f(i+1)+f(i+2)*f(i+2)
    end do
    simpson=m%dv*simpson/3.0
  end function simpson
  ! --------------------------------------------------------------------------------------
  !
  !              trapz()
  !
  ! --------------------------------------------------------------------------------------
  function  trapz(m,f)
    implicit none
    type (t_mesh)::m
    double precision::trapz
    double precision :: f(:)
    integer::i
     trapz=0.0
    do i=1,m%nactive-1
        trapz=trapz+f(i)*f(i)+f(i+1)*f(i+1)
    end do
    trapz=0.5*m%dv*trapz
  end function trapz
  ! --------------------------------------------------------------------------------------
  !
  !              norm()
  !
  ! --------------------------------------------------------------------------------------
  subroutine  norm(m,evec) 
    implicit none
    double precision :: evec(:),normloc
    double precision, external :: ddot
    type(t_mesh)::m
    if(m%dim.eq.3) then
       normloc=1.0/sqrt(m%dv*ddot(m%nactive,evec(:),1,evec(:),1))
    else if(m%dim.eq.2) then
       normloc=1.0/sqrt(m%dv*ddot(m%nactive,evec(:),1,evec(:),1))
    else    if(m%dim.eq.1) then
       !       normloc=1.0/sqrt(trapz(m,evec))
       !       print *,normloc,sqrt(trapz(m,evec))
       normloc=1.0/sqrt(simpson(m,evec))
    else
       print *,' STOP in norm(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
    call dscal(m%nactive,normloc,evec(:),1)
  end subroutine norm
  ! --------------------------------------------------------------------------------------
  !
  !              DIAGONALIZATION()
  !
  ! INPUT => H: Symmetric NxN matrix
  ! OUTPUT => H: eigenvectors (N)
  !                     S: eigenvalues (N)
  ! --------------------------------------------------------------------------------------
  subroutine diagonalization(S,H,N)
    implicit none
    integer :: N
    double precision :: H(:,:),S(:)

    integer :: lwork,info
    integer :: lwmax
    double precision,allocatable::work(:)
    parameter(lwmax=100000)
    allocate(work(lwmax))
    lwork=-1
    call  dsyev('vectors','Upper',N,H,N,S,work,lwork,info)
    lwork=min(lwmax,int(work(1)))
    if (lwork.ge.lwmax) then
       write(*,*) 'Diagonalization > WARNING info = ',info
       write(*,*) 'Diagonalization > WARNING lwork=',lwork
       write(*,*) 'Diagonalization > WARNING size of work(1)',int(work(1))
       stop
    end if
    call  dsyev('vectors','Upper',N,H,N,S,work,lwork,info)
    if(info.gt.0) then
       write(*,*) "Diagonalization > WARNING The algorithm computing failed"
       stop
    end if
    deallocate(work)
  end subroutine diagonalization
  ! --------------------------------------------------------------------------------------
  !
  !              Check_ortho()
  !
  ! INPUT => P(N,nvec) a Nxnvec matrix which contains a set of vectors (column)
  ! --------------------------------------------------------------------------------------
  subroutine check_ortho(P,nvec,N)
    implicit none
    integer :: nvec, N
    double precision :: P(:,:)
    double precision, parameter :: ETA=1.0e-6
    integer :: i,j,nfail
    double precision :: pscal
    double precision, external :: ddot
    nfail=-nvec
    print *,"--- check_ortho() ---"      
    do i=1,nvec
       do j=1,nvec
          pscal=ddot(N,P(:,i),1,P(:,j),1)
          write(*,'(F10.2)',advance='no') pscal
          if (pscal.gt.ETA) nfail=nfail+1
       end do
       write(*,*)
    end do
    if (nfail.gt.0) then
       print *,'### ERROR in check_ortho() ###',nfail,' fail(s)'
       stop
    end if
  end subroutine check_ortho
  ! --------------------------------------------------------------------------------------
  !
  !              compute_density()
  !
  ! --------------------------------------------------------------------------------------
       subroutine compute_density(molecule,param,V,mesh)
         implicit none
         type(t_molecule)::molecule
         type(t_param)::param
         type(t_mesh)::mesh
         double precision::V(:,:)
         integer::i,j
         double precision::charge
         double precision,allocatable::b(:)

         do j=1,param%noccstate
            call dcopy(molecule%mesh%nactive,V(:,j),1,molecule%wf%wfc(:,j),1)
            call norm(mesh,molecule%wf%wfc(:,j))
         end do



         molecule%rho=0.0
         do j=1,param%noccstate
            do i=1,mesh%nactive
               molecule%rho(i)=molecule%rho(i)-molecule%wf%wfc(i,j)**2
            end do
         end do
         charge=mesh%dv*sum(molecule%rho)
         print *,"Compute density > Charge ",charge
         print *,"Compute density > mesh%dv=",mesh%dv
         print *,"Compute density > mesh%dx=",mesh%dx
         print *,"Compute density > mesh%Nx=",mesh%Nx
         print *,"Compute density > mesh%Ntot=",mesh%Ntot
         print *,"Compute density > mesh%Nactive=",mesh%nactive
         print *,"Compute density > mesh%Nunactive=",mesh%nunactive
         do i=1,mesh%nbound
            mesh%bound(i)%val=charge/mesh%bound(i)%d
            !print *,mesh%bound(i)%q,mesh%bound(i)%val
         end do
         

    allocate(b(mesh%nactive))
    do i=1,mesh%nactive
       b(i)=molecule%rho(i)
       do j=1,mesh%n_bound(i)
          b(i)=b(i)+mesh%bound(mesh%list_bound(i,j))%val/mesh%dx**2
       end do
!       print *,b(i)
    end do
!    b(mesh%nactive-1)=b(mesh%nactive-1)-U(mesh%nactive)/mesh%dx**2
!    call Conjugate_gradient_3D(-b,pot%hartree,mesh%nactive,mesh%dx,mesh)    
  end subroutine compute_density
  ! --------------------------------------------------------------------------------------
  !
  !             interpolate()
  !
  ! --------------------------------------------------------------------------------------
  function interpolate(x,y,z,mesh,molecule,l)
    implicit none
    double precision::x,y,z
    integer :: i,j,k,n,i0,j0,k0,i1,j1,k1
    double precision::x0,y0,z0,x1,y1,z1
    type(t_mesh)::mesh
    type(t_molecule)::molecule
    integer,parameter::ndim=8
    double precision::A(ndim,ndim),C(ndim)
    integer::info,ipiv(ndim),l
    double precision::interpolate
    
    i0=floor(x/mesh%dx)
    j0=floor(y/mesh%dy)
    k0=floor(z/mesh%dz)

    i1=i0+1
    j1=j0+1
    k1=k0+1
    x0=i0*mesh%dx
    y0=j0*mesh%dy
    z0=k0*mesh%dz
    x1=i1*mesh%dx
    y1=j1*mesh%dy
    z1=k1*mesh%dz

    i=1;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y0 ; A(i,4)=z0 ; A(i,5) = x0*y0 ; A(i,6)=x0*z0 ; A(i,7) = y0*z0 ; A(i,8) = x0*y0*z0
    i=2;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y0 ; A(i,4)=z0 ; A(i,5) = x1*y0 ; A(i,6)=x1*z0 ; A(i,7) = y0*z0 ; A(i,8) = x1*y0*z0
    i=3;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y1 ; A(i,4)=z0 ; A(i,5) = x0*y1 ; A(i,6)=x0*z0 ; A(i,7) = y1*z0 ; A(i,8) = x0*y1*z0
    i=4;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y1 ; A(i,4)=z0 ; A(i,5) = x1*y1 ; A(i,6)=x1*z0 ; A(i,7) = y1*z0 ; A(i,8) = x1*y1*z0
    i=5;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y0 ; A(i,4)=z1 ; A(i,5) = x0*y0 ; A(i,6)=x0*z1 ; A(i,7) = y0*z1 ; A(i,8) = x0*y0*z1
    i=6;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y0 ; A(i,4)=z1 ; A(i,5) = x1*y0 ; A(i,6)=x1*z1 ; A(i,7) = y0*z1 ; A(i,8) = x1*y0*z1
    i=7;     A(i,1)=1.0 ; A(i,2)=x0 ; A(i,3) = y1 ; A(i,4)=z1 ; A(i,5) = x0*y1 ; A(i,6)=x0*z1 ; A(i,7) = y1*z1 ; A(i,8) = x0*y1*z1
    i=8;     A(i,1)=1.0 ; A(i,2)=x1 ; A(i,3) = y1 ; A(i,4)=z1 ; A(i,5) = x1*y1 ; A(i,6)=x1*z1 ; A(i,7) = y1*z1 ; A(i,8) = x1*y1*z1

    C=0.0
    i=i0 ; j=j0 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(1)=molecule%wf%wfc(n,l)
    i=i1 ; j=j0 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(2)=molecule%wf%wfc(n,l)
    i=i0 ; j=j1 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(3)=molecule%wf%wfc(n,l)
    i=i1 ; j=j1 ; k=k0 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(4)=molecule%wf%wfc(n,l)
    i=i0 ; j=j0 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(5)=molecule%wf%wfc(n,l)
    i=i1 ; j=j0 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(6)=molecule%wf%wfc(n,l)
    i=i0 ; j=j1 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(7)=molecule%wf%wfc(n,l)
    i=i1 ; j=j1 ; k=k1 ;     n=j+(i-1)*mesh%Ny+(k-1)*mesh%Ny*mesh%Nx;    if(n.le.mesh%nactive) C(8)=molecule%wf%wfc(n,l)

    ! do i=1,8
    !    print *,(A(i,j),j=1,8)
    ! end do
    ! do i=1,8
    !    print *,C(i)
    ! end do


!    CALL DGETRS('N' ,  ndim ,  ndim, A , ndim , IPIV, C , ndim , INFO)
    CALL DGESV( ndim , 1, A , ndim , IPIV, C , ndim , INFO)
    ! print *,info
    ! do i=1,8
    !    print *,'C(',i,')=',C(i)
    ! end do
    interpolate=C(1)+C(2)*x+C(3)*y+C(4)*z+C(5)*x*y+C(6)*x*z+C(7)*y*z+C(8)*x*y*z

  end function interpolate

end module tools
