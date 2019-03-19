program ppHibinitio
  implicit none
  character(len=1024)::filename
  character(len=1024)::line
  integer :: i,j,k,Nx,Ny,Nz,ifield,scan(6)
  double precision ::dx,dy,dz
  double precision,allocatable:: data(:,:,:)
  character(len=32) :: arg

  do i = 1, iargc() 
     call getarg(i, arg) 
     if (arg.eq.'-f')  then
        call getarg(i+1, filename) 
        write (6,*) filename
     end if
     if (arg.eq.'-scan')  then
        do j=1,6
           call getarg(i+j, arg)
           read(arg,*) scan(j)
        end do
        print *,scan        
     end if
  end do

!  read *,filename
!  call getarg(1, filename) 
 ! print *,filename
  open(unit=1,file=filename,form='formatted')
  do i=1,3
     read(1,*) ; 
  end do
  read(1,'(I5,F12.6,24X)') Nx,dx ;  print *,Nx,dx
  read(1,'(I5,12x,F12.6,12X)') Ny,dy ;  print *,Ny,dy
  read(1,'(I5,24x,F12.6)') Nz,dz ;  print *,Nz,dz
  read(1,*)
  allocate(data(Nx,Ny,Nz))
  
  do k=1,Nz
     ifield=0
     do i=1,Nx
        do j=1,Ny
           !nn=j+(i-1)*m%Ny+(k-1)*m%Ny*m%Nx               
           read(1,'(E13.5)',advance='no',end=999) data(i,j,k)
           !print *,i,j,k,data(i,j,k)
           ifield=ifield+1
           if (mod(ifield,6).eq.0) then
              ifield=0
              read(1,*)
           end if
        end do
     end do
     read(1,*)
  end do
  
  close(1)
999 continue

  open(unit=1,file='plot.dat',form='formatted',status='unknown')
    
  do k=scan(5),scan(6)
     do j=scan(3),scan(4)
        do i=scan(1),scan(2)
!           write(1,*) i*dx,j*dy,k*dz,data(i,j,k)
           write(1,*) i*dx,data(i,j,k)
        end do
     end do
  end do
  close(1)
  deallocate(data)  
end program ppHibinitio
