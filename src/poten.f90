module poten
  use global
  !  use IO
  implicit none
contains
    ! --------------------------------------------------------------------------------------
  !
  !             init_pot()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_pot(molecule)
    implicit none
    type(t_molecule)::molecule


    if(allocated(molecule%pot%ext))     deallocate(molecule%pot%ext)
    allocate(molecule%pot%ext(molecule%mesh%Ntot))
    if(allocated(molecule%pot%hartree))     deallocate(molecule%pot%hartree)
    allocate(molecule%pot%hartree(molecule%wf%nwfc,molecule%mesh%Ntot))
    if(allocated(molecule%pot%Ehartree))     deallocate(molecule%pot%Ehartree)
    allocate(molecule%pot%Ehartree(molecule%wf%nwfc))
    if(allocated(molecule%pot%Vx))     deallocate(molecule%pot%Vx)
    allocate(molecule%pot%Vx(molecule%mesh%Ntot))
    molecule%pot%hartree=0.0
    molecule%pot%Vx=0.0
    if(allocated(molecule%pot%perturb))     deallocate(molecule%pot%perturb)
    allocate(molecule%pot%perturb(molecule%mesh%Ntot))
    molecule%pot%perturb=0.0
    if(allocated(molecule%pot%tot))     deallocate(molecule%pot%tot)
    allocate(molecule%pot%tot(molecule%mesh%Ntot))
    call Vext2(molecule%mesh,molecule%pot%ext)
    call Vperturb(molecule%mesh,molecule%pot)
!    molecule%pot%tot=molecule%pot%ext+&
!         molecule%pot%hartree !+molecule%pot%perturb
  end subroutine init_pot
  ! --------------------------------------------------------------------------------------
  !
  !              Vext()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Vext(m,pot_ext)
    implicit none
    type(t_mesh) :: m
    double precision :: pot_ext(:)
    double precision :: pts(3),rsqr
    
    !    character (len=1024) :: filename
    integer :: i,j,k,nn
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           3D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(m%dim.eq.3) then
       do k=1,m%Nz
          pts(3)=k*m%dz
          do i=1,m%Nx
             pts(1)=i*m%dx
             do j=1,m%Ny
                pts(2)=j*m%dy
                rsqr=(pts(1)-m%box%center(1))**2&
                     +(pts(2)-m%box%center(2))**2&
                     +(pts(3)-m%box%center(3))**2
                nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx
!                pot_ext(nn)=0.5*1.0*rsqr
                pot_ext(nn)=-1.0/sqrt(rsqr)
             end do
          end do
       end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           2D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if(m%dim.eq.2) then
       open(unit=1,file="pot_ext.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          do j=1,m%Ny
             pts(2)=j*m%dy
             rsqr=(pts(1)-m%box%center(1))**2&
                  +(pts(2)-m%box%center(2))**2
             nn=j+(i-1)*m%Ny
             pot_ext(nn)=0.5*1.0*rsqr
             write(1,*) a02ang*pts(1),a02ang*pts(2),Ha2eV*pot_ext(nn)
          end do
       end do
       close(1)
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !           1D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else    if(m%dim.eq.1) then
       open(unit=1,file="pot_ext.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          rsqr=(pts(1)-m%box%center(1))**2
          pot_ext(i)=.5*1.0*rsqr
          write(1,*) a02ang*pts(1),Ha2eV*pot_ext(i)
       end do
       close(1)
    else
       print *,' STOP in Vext(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
       !stop
  end subroutine Vext
  ! --------------------------------------------------------------------------------------
  !
  !              Vext2()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Vext2(m,pot_ext)
    implicit none
    type(t_mesh) :: m
    double precision :: pot_ext(:)
    double precision :: pts(3),rsqr
    integer :: i,j,k,nn
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           3D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(m%dim.eq.3) then
       do nn=1,m%nactive
          pts(1)=m%node(nn)%q(1)
          pts(2)=m%node(nn)%q(2)
          pts(3)=m%node(nn)%q(3)
          rsqr=(pts(1)-m%box%center(1))**2&
               +(pts(2)-m%box%center(2))**2&
               +(pts(3)-m%box%center(3))**2
          pot_ext(nn)=-1.0/sqrt(rsqr)
       end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !           2D
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else if(m%dim.eq.2) then
       do nn=1,m%nactive
          pts(1)=m%node(nn)%q(1)
          pts(2)=m%node(nn)%q(2)
          rsqr=(pts(1)-m%box%center(1))**2+(pts(2)-m%box%center(2))**2
          pot_ext(nn)=0.0005*1.0*rsqr
       end do
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !
       !           1D
       !
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    else    if(m%dim.eq.1) then
       !
       ! Harmonic potential
       !
       do nn=1,m%nactive
          pts(1)=m%node(nn)%q(1)
          rsqr=(pts(1)-m%box%center(1))**2
          pot_ext(nn)=.5*1.0*rsqr*0
!          pot_ext(nn)=.5*0.01*rsqr
       end do
    else
       print *,' STOP in Vext(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
       !stop
  end subroutine Vext2
  ! --------------------------------------------------------------------------------------
  !
  !              Vperturb()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Vperturb(m,pot)
    implicit none
    type(t_mesh) :: m
    type(t_potential)::pot
    double precision :: pts(3),rsqr

    double precision,parameter::pi=4.0*atan(1.0)
    double precision :: invsig
    double precision :: facperturb
    integer :: i,j,nn

    facperturb=m%perturb%Intensity/sqrt(2*pi*m%perturb%sigma**2)
    invsig=0.5/m%perturb%sigma**2
    select case (m%dim)
       ! 3D case
    case (3)
       open(unit=1,file="pot_perturb.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          do j=1,m%Ny
             pts(2)=j*m%dy
             rsqr=(pts(1)-m%perturb%location(1))**2&
                  +(pts(2)-m%perturb%location(2))**2
             nn=j+(i-1)*m%Ny
             pot%perturb(nn)=facperturb*exp(-invsig*rsqr)
             write(1,*) pts(1),pts(2),pot%perturb(nn)
          end do
       end do
       close(1)
       ! 2D case
    case (2)
       open(unit=1,file="pot_perturb.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          do j=1,m%Ny
             pts(2)=j*m%dy
             rsqr=(pts(1)-m%perturb%location(1))**2&
                  +(pts(2)-m%perturb%location(2))**2
             nn=j+(i-1)*m%Ny
             pot%perturb(nn)=facperturb*exp(-invsig*rsqr)
             write(1,*) pts(1),pts(2),pot%perturb(nn)
          end do
       end do
       close(1)
       ! 1D case
    case (1)
       open(unit=1,file="pot_perturb.dat",form='formatted',status='unknown')
       do i=1,m%Nx
          pts(1)=i*m%dx
          rsqr=(pts(1)-m%perturb%location(1))**2
          pot%perturb(i)=facperturb*exp(-invsig*rsqr)
          write(1,*) pts(1),pot%perturb(i)
       end do
       close(1)
       ! other cases
    case default
       print *,' STOP in Vperturb(): dimension=',m%dim,' not yet implemented!'
    end select
  end subroutine Vperturb
  ! --------------------------------------------------------------------------------------
  !
  !              save_wavefunction(mesh,molecule)
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_potential(param,molecule)
    implicit none
    type(t_param)::param
    type(t_molecule)::molecule

    integer::nn
    character (len=1024) :: filename
    write(filename,'(a,a)') param%prefix(:len_trim(param%prefix)),'/pot_ext.agr'    
    open(unit=1,file=filename,form='formatted',status='unknown')
    write(1,*) "@with g0"
    write(1,*) '@    xaxis  label "x (ang.)"'
    write(1,*) '@    yaxis  label "Energy (eV)"'
    write(1,*) "@target G0.S0"
    write(1,*) "@type xy"
    do nn=1,molecule%mesh%nactive
       write(1,*) molecule%mesh%node(nn)%q(1),molecule%pot%ext(nn)
    end do
    write(1,*) "&"
    close(1)
  end subroutine save_potential

end module poten
