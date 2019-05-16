module mesh_mod
  use global
  implicit none
contains
  ! -----------------------------------------------
  !
  !       new_mesh(m)
  !
  ! -----------------------------------------------
  subroutine new_mesh(mesh)
    implicit none
    type(t_mesh)::mesh
    integer::i,j
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !                3D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(mesh%dim.eq.3) then
       mesh%Ny=mesh%Nx
       mesh%Nz=mesh%Nx
       mesh%dx=mesh%box%width/(mesh%Nx+1)
       mesh%dy=mesh%box%width/(mesh%Ny+1)
       mesh%dz=mesh%box%width/(mesh%Nz+1)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !                2D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if(mesh%dim.eq.2) then
       mesh%Ny=mesh%Nx
       mesh%Nz=1
       mesh%dx=mesh%box%width/(mesh%Nx+1)
       mesh%dy=mesh%box%width/(mesh%Ny+1)
       mesh%dz=1.0
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !                1D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else   if(mesh%dim.eq.1) then
       mesh%Ny=1
       mesh%Nz=1
       mesh%dx=mesh%box%width/(mesh%Nx+1)
       mesh%dy=1.0
       mesh%dz=1.0
    else
       print *,' STOP in new_mesh(): dimension=',mesh%dim,' not yet implemented!'
       stop
    end if
    mesh%dv=mesh%dx*mesh%dy*mesh%dz
    print *,"# mesh > dv=",mesh%dv
    ! -------------------------------------------------------------------------------
    !
    ! from now, we deal with the nodes belonging to mesh
    
    mesh%Ntot=mesh%Nx*mesh%Ny*mesh%Nz ! total number of nodes (active + unactive)
    !    allocate(mesh%n_neighbors(mesh%N))
    if(allocated(mesh%node)) deallocate(mesh%node)
    allocate(mesh%node(mesh%Ntot))

    mesh%multipole%lmax=2
    mesh%multipole%mmax=(mesh%multipole%lmax+1)**2

    if(allocated(mesh%multipole%rs)) deallocate(mesh%multipole%rs)
    allocate(mesh%multipole%rs(0:mesh%multipole%lmax))

    if(allocated(mesh%multipole%qlm)) deallocate(mesh%multipole%qlm)
    allocate(mesh%multipole%qlm(0:mesh%multipole%lmax))
    if(allocated(mesh%multipole%sph_harm_l)) deallocate(mesh%multipole%sph_harm_l)
    allocate(mesh%multipole%sph_harm_l(0:mesh%multipole%lmax))
    do i=0,mesh%multipole%lmax
       allocate(mesh%multipole%rs(i)%val(mesh%Ntot))
       allocate(mesh%multipole%qlm(i)%m(-i:i))
       allocate(mesh%multipole%sph_harm_l(i)%m(-i:i))
       do j=-i,i
          allocate(mesh%multipole%sph_harm_l(i)%m(j)%val(mesh%Ntot))
          allocate(mesh%multipole%qlm(i)%m(j)%val(1))
       end do
    end do
    do i=1,mesh%Ntot
       allocate(mesh%node(i)%list_neighbors(2*mesh%dim)) !
       mesh%node(i)%list_neighbors(:)=0
       mesh%node(i)%n_neighbors=0
       allocate(mesh%node(i)%list_bound(2*mesh%dim)) !
       mesh%node(i)%list_bound(:)=0
       mesh%node(i)%n_bound=0
    end do

    call set_nodes(mesh)
    ! 2*mesh%dim=2  @1D
    ! 2*mesh%dim=4  @2D
    ! 2*mesh%dim=6  @3D
    !allocate(mesh%list_neighbors(mesh%N,2*mesh%dim)) !
    !mesh%list_neighbors(:,:)=0
    !mesh%n_neighbors(:)=0
    ! n_bound -> number of inactive neighbors of a point
    ! default = 0
    ! max = 3 (corner)
    if(allocated(mesh%n_bound))     deallocate(mesh%n_bound)
    allocate(mesh%n_bound(mesh%Ntot))
    mesh%n_bound(:)=0
    ! list_bound -> idx of the inactive neighbors. It corresponds to
    !                       bound(:)
    if(allocated(mesh%list_bound))      deallocate(mesh%list_bound) !
    allocate(mesh%list_bound(mesh%Ntot,3)) !
    mesh%list_bound(:,:)=0
    ! number of element in the boundary surface
    mesh%nbound=8+4*(mesh%Nx+mesh%Ny+mesh%Nz)+&
         2*(mesh%Nx*mesh%Ny+mesh%Nx*mesh%Nz+&
         mesh%Ny*mesh%Nz)
    print *,'new_mesh > nbound=',mesh%nbound
    mesh%nbound=         2*(mesh%Nx*mesh%Ny+mesh%Nx*mesh%Nz+&
         mesh%Ny*mesh%Nz)
    print *,'new_mesh > nbound=',mesh%nbound
    if(allocated(mesh%bound))     deallocate(mesh%bound)
    allocate(mesh%bound(mesh%nbound))


    call compute_list_neighbors(mesh)
    call cart2sph(mesh)
  end subroutine new_mesh
  ! -------------------------------------------------------------------------------------------------------------------
  !
  !                          set_idx_list(mesh)
  !
  ! * There is Ntot nodes in the system ; some of them are
  !    "active" (we compute quantities at these nodes), others are
  !     or "unactive" (the values are fixed)
  ! * It is in this subroutine that we set up the state of a node  
  ! * Each node of the system is defined in the data type mesh by field mesh%node(idx)
  ! * Depending of the state of the node, its index is either in [1:meash%nactive] or
  !    in (mesh%nactive+1,mesh%Ntot]
  !    
  ! -------------------------------------------------------------------------------------------------------------------
  subroutine set_nodes(mesh)
    type(t_mesh)::mesh
    integer::i,j,k
    double precision::d,x,y,z
    logical::CompDomain
    character (len=1024) :: filename
    
    print *," Starting set_nodes()"
    if(allocated(mesh%ijk_to_idx))     deallocate(mesh%ijk_to_idx)
    allocate(mesh%ijk_to_idx(mesh%Nx,mesh%Ny,mesh%Nz))
    mesh%nactive=0
    mesh%nunactive=mesh%Ntot+1
    !    if((mesh%dim.eq.3).or.(mesh%dim.eq.2)) then   ! 3D
    do k=1,mesh%Nz
       z=k*mesh%dz
       do i=1,mesh%Nx
          x=i*mesh%dx
          do j=1,mesh%Ny
             y=j*mesh%dy
             CompDomain=.FALSE.
             select case (mesh%box%shape)
             case ("sphere")
                !                if(mesh%box%shape.eq."sphere") then
                d=((x-mesh%box%center(1))**2&
                     +(y-mesh%box%center(2))**2&
                     +(z-mesh%box%center(3))**2)**(0.5)
                if(d.le.mesh%box%radius) CompDomain=.TRUE.
             case ("cylinder")
                !                else                   if(mesh%box%shape.eq."cylinder") then
                d=((x-mesh%box%center(1))**2&
                     +(y-mesh%box%center(2))**2)**(0.5)
                if(d.le.mesh%box%radius) CompDomain=.TRUE.
             case("cube")
                !                   else                   if(mesh%box%shape.eq."cube") then
                if((abs(x-mesh%box%center(1)).le.mesh%box%radius).and.&
                     (abs(y-mesh%box%center(2)).le.mesh%box%radius).and.&
                     (abs(z-mesh%box%center(3)).le.mesh%box%radius)) &
                     CompDomain=.TRUE.
             case default
                !                   else
                print *,' STOP in set_idx_list(): undefined ',mesh%box%shape,'  shape!'
                stop
                   !                end if
             end select
             if(CompDomain) then
                   mesh%nactive=mesh%nactive+1
                   mesh%ijk_to_idx(i,j,k)%n=mesh%nactive
                   mesh%ijk_to_idx(i,j,k)%active=.TRUE.
                   mesh%node(mesh%nactive)%i=i
                   mesh%node(mesh%nactive)%j=j
                   mesh%node(mesh%nactive)%k=k
                   mesh%node(mesh%nactive)%q(1)=x
                   mesh%node(mesh%nactive)%q(2)=y
                   mesh%node(mesh%nactive)%q(3)=z
                   mesh%node(mesh%nactive)%usefull_unactive=.FALSE.
                   mesh%node(mesh%nactive)%active=.TRUE.
                else
                   mesh%nunactive=mesh%nunactive-1
                   mesh%ijk_to_idx(i,j,k)%n=mesh%nunactive
                   mesh%ijk_to_idx(i,j,k)%active=.FALSE.
                   mesh%node(mesh%nunactive)%i=i
                   mesh%node(mesh%nunactive)%j=j
                   mesh%node(mesh%nunactive)%k=k
                   mesh%node(mesh%nunactive)%q(1)=x
                   mesh%node(mesh%nunactive)%q(2)=y
                   mesh%node(mesh%nunactive)%q(3)=z
                   mesh%node(mesh%nunactive)%usefull_unactive=.FALSE.
                   mesh%node(mesh%nunactive)%active=.FALSE.
                end if
                mesh%ijk_to_idx(i,j,k)%q(1)=x
                mesh%ijk_to_idx(i,j,k)%q(2)=y
                mesh%ijk_to_idx(i,j,k)%q(3)=z
             end do
          end do
       end do

       !    else
       !       print *,' STOP in set_nodes(): dimension=',mesh%dim,' not yet implemented!'
       !       stop
       !    end if
       
       
       write(filename,'(a)') 'domain.xyz'
       open(unit=1,file=filename,form='formatted',status='unknown')
       write(1,*) mesh%Ntot
       write(1,*)
       do k=1,mesh%Nz
          do i=1,mesh%Nx
          do j=1,mesh%Ny
             if(mesh%ijk_to_idx(i,j,k)%active) then
                write(1,*) 'Cu ',&
                     mesh%ijk_to_idx(i,j,k)%q(1),&
                     mesh%ijk_to_idx(i,j,k)%q(2),&
                     mesh%ijk_to_idx(i,j,k)%q(3)
             else
                write(1,*) 'H ',&
                     mesh%ijk_to_idx(i,j,k)%q(1),&
                     mesh%ijk_to_idx(i,j,k)%q(2),&
                     mesh%ijk_to_idx(i,j,k)%q(3)
             end if
          end do
       end do
    end do
    close(1)
    print *,"# set_nodes > total number of nodes (active + unactive)" ,mesh%Ntot
    print *,"# set_nodes > ",mesh%nactive," actives nodes"
    print *,"# set_nodes > ",mesh%Ntot-mesh%nunactive+1," unactives nodes"
    

    ! do i=1,mesh%N
    !    print *,mesh%node(i)%active,mesh%node(i)%usefull_unactive,mesh%node(i)%n_bound
    !    if(mesh%node(i)%usefull_unactive) then
    !       mesh%n_usefull_unactive=mesh%n_usefull_unactive+1
    !    end if
    ! end do

    ! stop


!      mesh%N=mesh%nactive
!      stop

    print *,"End of set_nodes()"
  end subroutine set_nodes
  ! -----------------------------------------------
  !
  !          free_mesh(m)
  !
  ! -----------------------------------------------
  subroutine free_mesh(m)
    implicit none
    type(t_mesh) :: m
    !deallocate(m%n_neighbors)
   ! deallocate(m%list_neighbors) !
    deallocate(m%n_bound)
    deallocate(m%list_bound) !
  end subroutine free_mesh
  ! -----------------------------------------------
  !
  !        compute_list_neighbors(m)
  !
  ! -----------------------------------------------
  subroutine compute_list_neighbors(m)
    implicit none
    type(t_mesh) :: m
    integer::i,j,k,nn,idx

    if(m%dim.eq.3) then   ! 3D
       do nn=1,m%nactive
          i=m%node(nn)%i
          j=m%node(nn)%j
          k=m%node(nn)%k
          if(i.gt.1) then
             if(m%ijk_to_idx(i-1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i-1,j,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i-1,j,k)%n
                m%node(m%ijk_to_idx(i-1,j,k)%n)%usefull_unactive=.TRUE.

             end if
          end if
          if(i.lt.m%Nx) then
             if(m%ijk_to_idx(i+1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i+1,j,k)%n
                m%node(m%ijk_to_idx(i+1,j,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(j.gt.1) then
             if(m%ijk_to_idx(i,j-1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j-1,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j-1,k)%n
                m%node(m%ijk_to_idx(i,j-1,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(j.lt.m%Ny) then
             if(m%ijk_to_idx(i,j+1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j+1,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j+1,k)%n
                m%node(m%ijk_to_idx(i,j+1,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(k.gt.1) then
             if(m%ijk_to_idx(i,j,k-1)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j,k-1)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j,k-1)%n
                m%node(m%ijk_to_idx(i,j,k-1)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(k.lt.m%Nz) then
             if(m%ijk_to_idx(i,j,k+1)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j,k+1)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j,k+1)%n
                m%node(m%ijk_to_idx(i,j,k+1)%n)%usefull_unactive=.TRUE.
             end if
          end if
       end do
    else    if(m%dim.eq.2) then       ! 2D
       k=1
       do nn=1,m%nactive
          i=m%node(nn)%i
          j=m%node(nn)%j
          if(i.gt.1) then
             if(m%ijk_to_idx(i-1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i-1,j,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i-1,j,k)%n
                m%node(m%ijk_to_idx(i-1,j,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(i.lt.m%Nx) then
             if(m%ijk_to_idx(i+1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i+1,j,k)%n
                m%node(m%ijk_to_idx(i+1,j,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(j.gt.1) then
             if(m%ijk_to_idx(i,j-1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j-1,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j-1,k)%n
                m%node(m%ijk_to_idx(i-1,j,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
          if(j.lt.m%Ny) then
             if(m%ijk_to_idx(i,j+1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j+1,k)%n
             else
                m%node(nn)%n_bound=m%node(nn)%n_bound+1
                m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i,j+1,k)%n
                m%node(m%ijk_to_idx(i,j+1,k)%n)%usefull_unactive=.TRUE.
             end if
          end if
       end do
       else if(m%dim.eq.1) then      ! 1D
          k=1
          j=1
          do nn=1,m%nactive
             i=m%node(nn)%i
             if(i.gt.1) then
                if(m%ijk_to_idx(i-1,j,k)%active) then
                   m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                   m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i-1,j,k)%n
                else
                   m%node(nn)%n_bound=m%node(nn)%n_bound+1
                   m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i-1,j,k)%n
                   m%node(m%ijk_to_idx(i-1,j,k)%n)%usefull_unactive=.TRUE.
                end if
             end if
             if(i.lt.m%Nx) then
                if(m%ijk_to_idx(i+1,j,k)%active) then
                   m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                   m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
                else
                   m%node(nn)%n_bound=m%node(nn)%n_bound+1
                   m%node(nn)%list_bound(m%node(nn)%n_bound)=m%ijk_to_idx(i+1,j,k)%n
                   m%node(m%ijk_to_idx(i+1,j,k)%n)%usefull_unactive=.TRUE.
                end if
             end if
          end do
       else
       print *,' STOP in compute_list_neighbors(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
    
    m%n_usefull_unactive=0
    do i=1,m%Ntot
!       print *,m%node(i)%active,m%node(i)%usefull_unactive,m%node(i)%n_bound
       if(m%node(i)%usefull_unactive) then
          m%n_usefull_unactive=m%n_usefull_unactive+1
       end if
    end do
    print *,"# n usefull unactive nodes = ",m%n_usefull_unactive
!    stop
  end subroutine compute_list_neighbors
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! update_bound(idx,i,j,k,di,dj,dk,m)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine update_bound(idx,i,j,k,di,dj,dk,m)
  !        implicit none
  !        integer :: idx,i,j,k,nn,di,dj,dk
  !        type(t_mesh)::m
  !        nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
  !        m%n_bound(nn)=m%n_bound(nn)+1
  !        m%list_bound(nn,m%n_bound(nn))=idx
  !        m%bound(idx)%q(1)=m%dx*(i+di)
  !        m%bound(idx)%q(2)=m%dy*(j+dj)
  !        m%bound(idx)%q(3)=m%dz*(k+dk)
  !        m%bound(idx)%d=sqrt((m%bound(idx)%q(1)-m%box%center(1))**2+&
  !             (m%bound(idx)%q(2)-m%box%center(2))**2+&
  !             (m%bound(idx)%q(3)-m%box%center(3))**2)
  !        idx=idx+1
  !      end subroutine update_bound



  ! -----------------------------------------------
  !
  !        cart2sph(mesh)
  !
  ! -----------------------------------------------
  subroutine cart2sph(mesh)
    implicit none
    type(t_mesh)::mesh
    double precision::x,y,z,r,theta,phi
    integer::i,l,m
    do i=1,mesh%Ntot
       x=mesh%node(i)%q(1)-mesh%box%center(1)
       y=mesh%node(i)%q(2)-mesh%box%center(2)
       z=mesh%node(i)%q(3)-mesh%box%center(3)
       r=sqrt(x*x+y*y+z*z)
       mesh%node(i)%r=r
       theta=calc_theta(x,y,z)
       phi=calc_phi(x,y)
       mesh%node(i)%theta=theta
       mesh%node(i)%phi=phi
       do l=0,mesh%multipole%lmax
          if(l.eq.0) then
             mesh%multipole%rs(l)%val(i)=1.0
          else
             mesh%multipole%rs(l)%val(i)=r**l
          end if
          do m=-l,l
             mesh%multipole%sph_harm_l(l)%m(m)%val(i)=func_sph_harm(l,m,theta,phi)
          end do
       end do
    end do
  end subroutine cart2sph

  ! -----------------------------------------------
  !
  !        calc_theta(x,y,z)
  !
  ! -----------------------------------------------
  function calc_theta(x,y,z)
    implicit none
    double precision::x,y,z,calc_theta
    calc_theta=atan(sqrt(x*x+y*y)/z);
    if(z.lt.0)calc_theta=pi+calc_theta
        if((x.eq.0.0).and.(y.eq.0.0).and.(z.eq.0.0)) calc_theta=0.0
    return
  end function calc_theta
  
  ! -----------------------------------------------
  !
  !        calc_phi(x,y)
  !
  ! -----------------------------------------------
  function calc_phi(x,y)
    implicit none
    double precision::x,y,calc_phi
    calc_phi=atan(y/x);
    if((x.lt.0.0).and.(y.gt.0.0)) calc_phi=pi+calc_phi
    if((x.lt.0.0).and.(y.lt.0.0)) calc_phi=pi+calc_phi
    if((x.gt.0.0).and.(y.lt.0.0)) calc_phi=2*pi+calc_phi
    if((x.eq.0.0).and.(y.eq.0.0)) calc_phi=0.0
    return
  end function calc_phi
  ! --------------------------------------------------------------------------------------
  !
  !  func_sph_harm(l,m,theta,phi) --> Ylm(theta,phi) COMPLEX
  !
  ! --------------------------------------------------------------------------------------
  function func_sph_harm(l,morig,theta,phi)
    implicit none
    integer::l,m,i,morig
    double precision::theta,phi,fac1,fac2,fac3
    double complex::func_sph_harm
    double precision::costheta

    m=abs(morig)
    !    double precision,external::plgndr
    costheta=cos(theta)
    fac1=1.0
    do i=2,l-m
       fac1=fac1*i
    end do
    fac2=fac1
    do i=l-m+1,l+m
       fac2=fac2*i
    end do
    fac3=sqrt((2*l+1)*fac1/(4*pi*fac2))*plgndr(l,m,costheta)
    
    if(morig.lt.0) then
       fac3=fac3*(-1)**m
       func_sph_harm=cmplx(fac3*cos(m*phi),-fac3*sin(m*phi))
    else
       func_sph_harm=cmplx(fac3*cos(m*phi),fac3*sin(m*phi))
    end if

    
    
  end function func_sph_harm
  ! --------------------------------------------------------------------------------------
  !
  !      Legendre polynomials, associated (spherical harmonics) [6.8]
  !      From http://sciold.ui.ac.ir/~sjalali/nrf/
  !      Compute the associated Legendre polynomial Plm(x)
  !      l and m are integers satisfying 0 <= m <= l, while x lies in
  !      the range -1 <= x <= 1
  ! --------------------------------------------------------------------------------------
  FUNCTION plgndr(l,m,x)  
    implicit none
    integer::l,m
    double precision::x,plgndr
    double precision::pmm,somx2,fact,pmmp1,pll
    integer::i,ll

    if((m.lt.0).or.(m.gt.l).or.(abs(x).gt.1.0)) then
       print *,"### ERROR in plgndr"
       print *,"### l=",l
       print *,"### m=",m
       print *,"### abs(x)=",abs(x)
       call exit()
    end if

    pmm=1.0
    if(m.gt.0) then
       somx2=sqrt((1.0-x)*(1.0+x))
       fact=1.0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.0
       end do
    end if
    if(l.eq.m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.(m+1)) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          end do
          plgndr=pll
       end if
    end if
    return
  end FUNCTION plgndr

end module mesh_mod
