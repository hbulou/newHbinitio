module mesh_mod
  use global
  use param_mod
  implicit none
contains
  ! -----------------------------------------------
  !
  !       new_mesh(m,param)
  !
  ! -----------------------------------------------
  subroutine new_mesh(mesh,param)
    implicit none
    type(t_mesh)::mesh
    type(t_param)::param
    integer::i
    mesh%dim=param%dim
    mesh%box%width=param%box%width
    mesh%box%shape=param%box%shape
    mesh%box%radius=param%box%radius
    mesh%box%center(1)=param%box%center(1)*mesh%box%width
    mesh%box%center(2)=param%box%center(2)*mesh%box%width
    mesh%box%center(3)=param%box%center(3)*mesh%box%width
    mesh%Nx=param%Nx
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

    ! -------------------------------------------------------------------------------
    !
    ! from now, we deal with the nodes belonging to mesh
    
    mesh%N=mesh%Nx*mesh%Ny*mesh%Nz ! total number of nodes (active + unactive)
!    allocate(mesh%n_neighbors(mesh%N))
    allocate(mesh%node(mesh%N))
    do i=1,mesh%N
       allocate(mesh%node(i)%list_neighbors(2*mesh%dim)) !
       mesh%node(i)%list_neighbors(:)=0
       mesh%node(i)%n_neighbors=0
       allocate(mesh%node(i)%list_bound(2*mesh%dim)) !
       mesh%node(i)%list_bound(:)=0
       mesh%node(i)%n_bound=0
    end do

    call set_idx_list(mesh)
    ! 2*mesh%dim=2  @1D
    ! 2*mesh%dim=4  @2D
    ! 2*mesh%dim=6  @3D
    !allocate(mesh%list_neighbors(mesh%N,2*mesh%dim)) !
    !mesh%list_neighbors(:,:)=0
    !mesh%n_neighbors(:)=0
    ! n_bound -> number of inactive neighbors of a point
    ! default = 0
    ! max = 3 (corner)
    allocate(mesh%n_bound(mesh%N))
    mesh%n_bound(:)=0
    ! list_bound -> idx of the inactive neighbors. It corresponds to
    !                       bound(:)
    allocate(mesh%list_bound(mesh%N,3)) !
    mesh%list_bound(:,:)=0
    ! number of element in the boundary surface
    mesh%nbound=8+4*(mesh%Nx+mesh%Ny+mesh%Nz)+&
         2*(mesh%Nx*mesh%Ny+mesh%Nx*mesh%Nz+&
         mesh%Ny*mesh%Nz)
    print *,'new_mesh > nbound=',mesh%nbound
    mesh%nbound=         2*(mesh%Nx*mesh%Ny+mesh%Nx*mesh%Nz+&
         mesh%Ny*mesh%Nz)
    print *,'new_mesh > nbound=',mesh%nbound
    allocate(mesh%bound(mesh%nbound))


    call compute_list_neighbors(mesh)
  end subroutine new_mesh
  ! -----------------------------------------------
  !
  ! set_idx_list(mesh)
  !
  ! -----------------------------------------------
    subroutine set_idx_list(mesh)
      type(t_mesh)::mesh
      integer::i,j,k
      double precision::d,x,y,z
      logical::CompDomain
      allocate(mesh%ijk_to_idx(mesh%Nx,mesh%Ny,mesh%Nz))
      mesh%nactive=0
      mesh%nunactive=mesh%N+1
      if((mesh%dim.eq.3).or.(mesh%dim.eq.2)) then   ! 3D
         do k=1,mesh%Nz
            z=(k-1)*mesh%dz
            do i=1,mesh%Nx
               x=(i-1)*mesh%dx
               do j=1,mesh%Ny
                  y=(j-1)*mesh%dy
                  CompDomain=.FALSE.
                  if(mesh%box%shape.eq."sphere") then
                     d=((x-mesh%box%center(1))**2&
                          +(y-mesh%box%center(2))**2&
                          +(z-mesh%box%center(3))**2)**(0.5)
                     if(d.le.mesh%box%radius) CompDomain=.TRUE.
                  else                   if(mesh%box%shape.eq."cylinder") then
                     d=((x-mesh%box%center(1))**2&
                          +(y-mesh%box%center(2))**2)**(0.5)
                     if(d.le.mesh%box%radius) CompDomain=.TRUE.
                  else                   if(mesh%box%shape.eq."cube") then
                     if((abs(x-mesh%box%center(1)).le.mesh%box%radius).and.&
                          (abs(y-mesh%box%center(2)).le.mesh%box%radius).and.&
                          (abs(z-mesh%box%center(3)).le.mesh%box%radius)) &
                          CompDomain=.TRUE.
                  else
                     print *,' STOP in set_idx_list(): undefined ',mesh%box%shape,'  shape!'
                     stop
                  end if

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
                  end if
                  mesh%ijk_to_idx(i,j,k)%q(1)=x
                  mesh%ijk_to_idx(i,j,k)%q(2)=y
                  mesh%ijk_to_idx(i,j,k)%q(3)=z
               end do
            end do
         end do
      else
         print *,' STOP in set_idx_list(): dimension=',mesh%dim,' not yet implemented!'
         stop
      end if

      open(unit=1,file='domain.xyz',form='formatted',status='unknown')
      write(1,*) mesh%N
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
      mesh%N=mesh%nactive
!      stop
    end subroutine set_idx_list
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
             end if
          end if
          if(i.lt.m%Nx) then
             if(m%ijk_to_idx(i+1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
             end if
          end if
          if(j.gt.1) then
             if(m%ijk_to_idx(i,j-1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j-1,k)%n
             end if
          end if
          if(j.lt.m%Ny) then
             if(m%ijk_to_idx(i,j+1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j+1,k)%n
             end if
          end if
          if(k.gt.1) then
             if(m%ijk_to_idx(i,j,k-1)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j,k-1)%n
             end if
          end if
          if(k.lt.m%Nz) then
             if(m%ijk_to_idx(i,j,k+1)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j,k+1)%n
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
             end if
          end if
          if(i.lt.m%Nx) then
             if(m%ijk_to_idx(i+1,j,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
             end if
          end if
          if(j.gt.1) then
             if(m%ijk_to_idx(i,j-1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j-1,k)%n
             end if
          end if
          if(j.lt.m%Ny) then
             if(m%ijk_to_idx(i,j+1,k)%active) then
                m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i,j+1,k)%n
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
                end if
             end if
             if(i.lt.m%Nx) then
                if(m%ijk_to_idx(i+1,j,k)%active) then
                   m%node(nn)%n_neighbors=m%node(nn)%n_neighbors+1
                   m%node(nn)%list_neighbors(m%node(nn)%n_neighbors)=m%ijk_to_idx(i+1,j,k)%n
                end if
             end if
          end do
       else
       print *,' STOP in compute_list_neighbors(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
  end subroutine compute_list_neighbors
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! update_bound(idx,i,j,k,di,dj,dk,m)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine update_bound(idx,i,j,k,di,dj,dk,m)
         implicit none
         integer :: idx,i,j,k,nn,di,dj,dk
         type(t_mesh)::m
         nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
         m%n_bound(nn)=m%n_bound(nn)+1
         m%list_bound(nn,m%n_bound(nn))=idx
         m%bound(idx)%q(1)=m%dx*(i+di)
         m%bound(idx)%q(2)=m%dy*(j+dj)
         m%bound(idx)%q(3)=m%dz*(k+dk)
         m%bound(idx)%d=sqrt((m%bound(idx)%q(1)-m%box%center(1))**2+&
              (m%bound(idx)%q(2)-m%box%center(2))**2+&
              (m%bound(idx)%q(3)-m%box%center(3))**2)
         idx=idx+1
       end subroutine update_bound

end module mesh_mod
