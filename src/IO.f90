module IO
  use global
  use tools
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              make_filename()
  !
  ! --------------------------------------------------------------------------------------
    subroutine make_filename(idx,pref,filename)
      implicit none
      integer::idx
      character (len=1024) :: filename,pref

      if(idx.le.9) then
         write(filename,'(A,A,I1,A)') trim(pref),'0000000000',idx,'.png'
      else        if(idx.le.99) then
         write(filename,'(A,A,I2,A)') trim(pref),'000000000',idx,'.png'
      else        if(idx.le.999) then
         write(filename,'(A,A,I3,A)') trim(pref),'00000000',idx,'.png'
      else        if(idx.le.9999) then
         write(filename,'(A,A,I4,A)') trim(pref),'000000',idx,'.png'
      else        if(idx.le.99999) then
         write(filename,'(A,A,I5,A)') trim(pref),'00000',idx,'.png'
      else        if(idx.le.999999) then
         write(filename,'(A,A,I6,A)') trim(pref),'0000',idx,'.png'
      else        if(idx.le.9999999) then
         write(filename,'(A,A,I7,A)') trim(pref),'000',idx,'.png'
      else        if(idx.le.99999999) then
         write(filename,'(A,A,I8,A)') trim(pref),'00',idx,'.png'
      else        if(idx.le.999999999) then
         write(filename,'(A,A,I9,A)') trim(pref),'0',idx,'.png'
      else
         print *,"Too large number of steps !"
         call exit()
      end if
    end subroutine make_filename
  ! --------------------------------------------------------------------------------------
  !
  !              SAVE_AGR()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_agr(idxmov)
    implicit none
    integer::idxmov
    character (len=1024) :: filename,pref

    open(unit=1,file="plot_wfc.bfile",form='formatted',status='unknown')
    write(1,*) 'ARRANGE (2,1,.1,.6,.6,ON,ON,ON)'
    write(1,*) 'FOCUS G0'
    write(1,*) 'VIEW 0.1,0.1,0.9,0.5'
    write(1,*) 'READ BLOCK "tdse.dat"'
    write(1,*) 'BLOCK xy "1:4"'
    write(1,*) 's0 LINEWIDTH 2.0'
    write(1,*) 'READ BLOCK "wfc.dat"'
    write(1,*) 'BLOCK xy "1:2"'
    write(1,*) 's1 LINEWIDTH 2.0'
    write(1,*) 'BLOCK xy "1:3"'
    write(1,*) 's2 LINEWIDTH 2.0'
    write(1,*) 'BLOCK xy "1:4"'
    write(1,*) 's3 LINEWIDTH 2.0'
    write(1,*) 's3 line color 8'
    write(1,*) 's4 line color 10'
    write(1,*) 's4 LINEWIDTH 4.0'
    write(1,*) 'READ BLOCK "pot_ext.dat"'
    write(1,*) 'BLOCK xy "1:2"'
    write(1,*) 'BLOCK xy "1:3"'
    write(1,*) 'WORLD YMIN -.4'
    write(1,*) 'WORLD YMAX .4'
    write(1,*) 'FOCUS G1'
    write(1,*) 'VIEW 0.1,0.6,0.9,0.9'
    write(1,*) 'READ BLOCK "fft.dat"'
    write(1,*) 'BLOCK xy "1:4"'


    
    write(pref,'(A)') 'output'
    call make_filename(idxmov,pref,filename)
    
    write(1,'(A,A,A)') 'PRINT TO "',trim(filename),'"'

    !write(1,*) trim(filename)
    !print *,filename
    !call exit()
    write(1,*) 'HARDCOPY DEVICE "PNG"'
    write(1,*) 'PRINT'
    close(1)
    !       call system("xmgrace -batch plot_wfc.bfile -nosafe -hardcopy")
    call execute_command_line("xmgrace -batch plot_wfc.bfile -nosafe -hardcopy",WAIT=.TRUE.)
    idxmov=idxmov+1
  end subroutine save_agr
  ! --------------------------------------------------------------------------------------
  !
  !              SAVE_CUBE()
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_cube_3D(data,filename,m)
    implicit none
    double precision :: data(:),val
!    integer :: idxmin,idxmax
    type(t_mesh) :: m
    character (len=1024) :: filename
    
    character(len=*),parameter :: FMT1='(I5,3F12.6)'
    integer :: i,j,k,nn,ifield

    print *,"# save_cube_3D > saving ",trim(filename)
    
    open(unit=1,file=filename,form='formatted',status='unknown')
    write(1,*) ' Cubefile created from Hbinitio.f90 calculation'
    write(1,*) ' H. Bulou, November 2018'
    write(1,FMT1) 1,-.5*m%box%center(1),-.5*m%box%center(2),-.5*m%box%center(3)
    write(1,FMT1) m%Nx,m%dx,0.0,0.0
    write(1,FMT1) m%Ny,0.0,m%dy,0.0
    write(1,FMT1) m%Nz,0.0,0.0,m%dz
    write(1,'(I5,4F12.6)') 1,1.0,-.5*m%box%center(1),-.5*m%box%center(2),-.5*m%box%center(3)
    do k=1,m%Nz
       ifield=0
       do i=1,m%Nx
          do j=1,m%Ny
             nn=m%ijk_to_idx(i,j,k)%n
             val=data(nn)
!             if(.not.(m%node(nn)%active)) val=0.0
             write(1,'(E13.5)',advance='no') val
             ifield=ifield+1
             if (mod(ifield,6).eq.0) then
                ifield=0
                write(1,*)
             end if
          end do
       end do
       write(1,*)
    end do
    close(1)
!    print *,"# save_cube_3D > ok"
  end subroutine save_cube_3D
  ! --------------------------------------------------------------------------------------
  !
  !              save_wavefunction(param,mesh,V,molecule)
  !
  ! --------------------------------------------------------------------------------------
  subroutine save_wavefunction(param,mesh,V,molecule)
    implicit none
    type(t_param)::param
    type(t_mesh)::mesh
    double precision::V(:,:)
    type(t_molecule)::molecule
    
    integer :: i,j,k,nn
    double precision::val
    character (len=1024) :: filecube
    
    do i=1,param%nvec_to_cvg
       call norm(mesh,V(:,i))
       call dcopy(molecule%mesh%nactive,V(:,i),1,molecule%wf%wfc(:,i),1)

       select case(mesh%dim)
          ! case 3D
       case(3)
          write(filecube,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/evec',i,'.cube'
          call save_cube_3D(V(:,i),filecube,mesh)
          ! case 2D
       case(2)
          write(filecube,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/evec',i,'.dat'
          open(unit=1,file=filecube,form='formatted',status='unknown')
          do j=1,mesh%Nx
             do k=1,mesh%Ny
                nn=mesh%ijk_to_idx(j,k,1)%n
                val=0.0
                if(mesh%node(nn)%active) val=V(nn,i)

                write(1,*) j*mesh%dx,k*mesh%dy,val
             end do
          end do
          close(1)
          ! case 1D
       case(1)
          write(filecube,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/evec',i,'.dat'
!          write(filecube,'(a,i0,a)') 'evec',i,'.dat'
          open(unit=1,file=filecube,form='formatted',status='unknown')
          do j=1,mesh%nactive
             write(1,*) j*mesh%dx,V(j,i)
          end do
          close(1)
       case default
          print *,' STOP in main(): dimension=',mesh%dim,' not yet implemented!'
          stop
       end select
    end do

  end subroutine save_wavefunction
  ! --------------------------------------------------------------------------------------
  !
  !                             save_config()
  !
  ! subroutine to save the configuration of the calculation in order to restart it
  ! later if necessary
  ! --------------------------------------------------------------------------------------
  subroutine save_config(V,m,nvecmin,param)
    implicit none
    type(t_mesh)::m
    type(t_param)::param
    double precision :: V(:,:)
    integer::nvecmin,i,j
    character (len=1024)::filename

    write(filename,'(a,a)') param%prefix(:len_trim(param%prefix)),'/evectors.dat'
    print *,"# save_config > saving ",trim(filename)
    open(unit=1,file=filename,form='formatted',status='unknown')
    write(1,*) "# nmesh=",m%nactive
    write(1,*) "# nvec=",nvecmin
    do i=1,m%nactive
       write(1,*) (V(i,j),j=1,nvecmin)
    end do
    close(1)
  end subroutine save_config
  ! --------------------------------------------------------------------------------------
  !
  !              read_config()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_config(V,m,nvecmin)
    implicit none
    type(t_mesh)::m
    double precision :: V(:,:)
    integer::nvecmin,i,j
    logical :: file_exists
    integer::npoints, nvec
    character(len=128)::junk
    
    INQUIRE(FILE="evectors.dat", EXIST=file_exists)
    if(file_exists) then
       open(unit=1,file="evectors.dat",form='formatted',status='unknown')
       read(1,*)  junk,junk,npoints
       read(1,*)    junk,junk,nvec
       print *,"npoints=",npoints," nvec=",nvec
       do i=1,m%nactive
          read(1,*) (V(i,j),j=1,nvecmin)
       end do
       close(1)
    else
       print *,"### ERROR ### evectors.dat doesn't exist"
       stop
    end if
  end subroutine read_config

end module IO
