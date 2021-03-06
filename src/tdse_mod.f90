module tdse_mod
    use global
    use tools
    use IO
    use FFT_mod
    implicit none
contains
    subroutine tdse(molecule,param,wfc_ini)
      implicit none
      type(t_molecule)::molecule
      type(t_param)::param
      double complex,allocatable::wfc_ini(:)
      
      double complex,allocatable::wfc(:,:)
      double complex,allocatable::in_wfc(:)
      double complex,allocatable::fft_wfc(:)
      
      double precision::x,sigsqr,sig,k0
      double precision::dxsqr,dt,dxsqrdt,c,d,normesqr
      double precision::a,b,aa,bb
      double precision::g,h
      character (len=1024) :: filename
      integer::i,oldt,newt,itmp,time_loop,idxmov

      double precision::lambda,x1,epsilonx1,x0,dk
      double precision,allocatable::mu(:)
      double complex,allocatable::e(:)
      double complex,allocatable::f(:)
      double complex,allocatable::Omega(:)

      
      dk=1.0/(molecule%mesh%nactive*molecule%mesh%dx)

      
!      double precision::a0
!      a0=0.529
      
      oldt=1
      newt=2
      allocate(e(molecule%mesh%nactive))
      allocate(f(molecule%mesh%nactive))
      allocate(wfc(molecule%mesh%nactive,2))
      allocate(fft_wfc(molecule%mesh%nactive))
      allocate(in_wfc(molecule%mesh%nactive))
      allocate(Omega(molecule%mesh%nactive))
      allocate(mu(molecule%mesh%nactive))

      wfc(:,oldt)=wfc_ini
      
      !sig=0.05*molecule%mesh%box%radius
      !sig=1.6*ang2a0
      !print *, "sigma=",sig*a02ang,' ang.'

      !lambda=1.6*ang2a0
      !k0=0*2*pi/(.05*molecule%mesh%box%width)
      !sigsqr=sig*sig
      !x0=molecule%mesh%box%width/4
      !do i=1,molecule%mesh%nactive
      !   x=i*molecule%mesh%dx-x0
      !   wfc(i,oldt)=cmplx(exp(-x*x/sigsqr)*cos(2*pi*x/lambda),exp(-x*x/sigsqr)*sin(2*pi*x/lambda))
      !end do
      call norm(molecule%mesh,dreal(wfc(:,oldt)))
      write(filename,'(a)') 'tdse4.dat'
      open(unit=1,file=filename,form='formatted',status='unknown')
      do i=1,molecule%mesh%nactive
         write(1,*) i*molecule%mesh%dx,dreal(wfc(i,oldt)),dimag(wfc(i,oldt)),abs(wfc(i,oldt))
      end do
      close(1)



      
!      epsilonx1=0.5
!      x1=60.0
!      do i=1,molecule%mesh%nactive
!         x=i*molecule%mesh%dx-2*x0
!         if(abs(x).lt.10) then
!            molecule%pot%tot(i)=1.0
!         else
!            molecule%pot%tot(i)=0
!         end if


      !      call exit()
      !      end do
      open(unit=1,file="pot_ext.dat",form='formatted',status='unknown')
      do i=1,molecule%mesh%nactive
         mu(i)=0.0
         if(i*molecule%mesh%dx.ge.25) mu(i)=1.0*(i*molecule%mesh%dx-25)
         !         mu(i)=1.0/(molecule%mesh%nactive*molecule%mesh%dx-(i-1)*molecule%mesh%dx)
         write(1,*) i*molecule%mesh%dx,molecule%pot%tot(i),mu(i)
      end do
      close(1)


      
      dt=1.0e-4
      dxsqr=molecule%mesh%dx**2
      dxsqrdt=dxsqr/dt
      print *,"dx=",molecule%mesh%dx
      print *,"dt=",dt
      print *,"dx^2/dt=",dxsqrdt
      ! !!!!!!!!!!!!!!!!!!!!!
      !
      ! e calculation
      !
      ! !!!!!!!!!!!!!!!!!!!!!
      aa=2+2*dxsqr*molecule%pot%tot(i)
      bb=-4*dxsqrdt-2*dxsqr*mu(i)
      i=1 ;     e(i)=cmplx(aa,bb)
      do i=2,molecule%mesh%nactive-1
         a=dreal(e(i-1))
         b=dimag(e(i-1))
         normesqr=a**2+b**2
         aa=2+2*dxsqr*molecule%pot%tot(i)-a/normesqr
         bb=b/normesqr-4*dxsqrdt-dxsqr*mu(i)
         e(i)=cmplx(aa,bb)
      end do
!      open(unit=1,file="e.dat",form='formatted',status='unknown')
!      do i=1,molecule%mesh%nactive
!         write(1,*) i*molecule%mesh%dx,dreal(e(i)),dimag(e(i)),abs(e(i))
!      end do
!      close(1)

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! Time loop
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      idxmov=1
      do time_loop=1,molecule%param%tdse%nstep
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! Omega calculation
         !
         ! !!!!!!!!!!!!!!!!!!!!!!!!!


         i=1
         Omega(i)=cmplx(&
              (2+2*dxsqr*molecule%pot%tot(i))*dreal(wfc(i,oldt))+(2*dxsqr*mu(i)-4*dxsqrdt)*dimag(wfc(i,oldt))-dreal(wfc(i+1,oldt)),&
              (2+2*dxsqr*molecule%pot%tot(i))*dimag(wfc(i,oldt))+(4*dxsqrdt-2*dxsqr*mu(i))*dreal(wfc(i,oldt))-dimag(wfc(i+1,oldt)))
         do i=2,molecule%mesh%nactive-1
            Omega(i)=cmplx(&
                 (2+2*dxsqr*molecule%pot%tot(i))*dreal(wfc(i,oldt))+(2*dxsqr*mu(i)-4*dxsqrdt)*dimag(wfc(i,oldt))&
                 -dreal(wfc(i+1,oldt))&
                 -dreal(wfc(i-1,oldt)),&
                 (2+2*dxsqr*molecule%pot%tot(i))*dimag(wfc(i,oldt))+(4*dxsqrdt-2*dxsqr*mu(i))*dreal(wfc(i,oldt))&
                 -dimag(wfc(i+1,oldt))&
                 -dimag(wfc(i-1,oldt)))
         end do
         
         ! open(unit=1,file="Omega.dat",form='formatted',status='unknown')
         ! do i=1,molecule%mesh%nactive-1
         !    write(1,*) i*molecule%mesh%dx,dreal(Omega(i)),dimag(Omega(i)),abs(Omega(i))
         ! end do
         ! close(1)
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! f calculation
         !
         ! !!!!!!!!!!!!!!!!!!!!!!!!!
         i=1 ;     f(i)=Omega(i)
         do i=2,molecule%mesh%nactive-1
            a=dreal(e(i-1))
            b=dimag(e(i-1))
            normesqr=a**2+b**2
            c=dreal(f(i-1))
            d=dimag(f(i-1))
            f(i)=Omega(i)+cmplx(c*a/normesqr+b*d/normesqr,d*a/normesqr-c*b/normesqr)
         end do
         
         ! open(unit=1,file="f.dat",form='formatted',status='unknown')
         ! do i=1,molecule%mesh%nactive-1
         !    write(1,*) i*molecule%mesh%dx,dreal(f(i)),dimag(f(i)),abs(f(i))
         ! end do
         ! close(1)
         
         
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         ! new wfc calculation
         !
         ! !!!!!!!!!!!!!!!!!!!!!!!!!
         
         i=molecule%mesh%nactive-1
         a=dreal(e(i))
         b=dimag(e(i))
         normesqr=a**2+b**2
         c=dreal(f(i))
         d=dimag(f(i))
         wfc(i,newt)=cmplx(-(c*a+b*d)/normesqr,-(a*d-c*b)/normesqr)
         do i=molecule%mesh%nactive-2,1,-1
            a=dreal(e(i))
            b=dimag(e(i))
            normesqr=a**2+b**2
            c=dreal(f(i))
            d=dimag(f(i))
            g=dreal(wfc(i+1,newt))
            h=dimag(wfc(i+1,newt))
            wfc(i,newt)=cmplx((a*(g-c)+b*(h-d))/normesqr,(a*(h-d)-b*(g-c))/normesqr)
         end do


         
         if(mod(time_loop,molecule%param%tdse%freq_save).eq.0) then

            !write(filename,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/wfc',time_loop,'.dat'
           ! print *,"saving ",trim(filename)
           ! open(unit=1,file=filename,form='formatted',status='unknown')
           ! do i=1,molecule%mesh%nactive-1
           !    write(1,*) i*molecule%mesh%dx,dreal(wfc(i,newt)),dimag(wfc(i,newt)),abs(wfc(i,newt))
           ! end do
           ! close(1)

            in_wfc=wfc(:,newt)
            call FFT(in_wfc,fft_wfc,molecule%mesh%nactive)
!            write(filename,'(a,a,i0,a)') param%prefix(:len_trim(param%prefix)),'/fft_wfc',time_loop,'.dat'
!            open(unit=1,file=filename,form='formatted',status='unknown')
            open(unit=1,file="fft.dat",form='formatted',status='unknown')
            do i=molecule%mesh%nactive/2+1,molecule%mesh%nactive
               write(1,*) (i-molecule%mesh%nactive-1)*dk,dreal(fft_wfc(i)),dimag(fft_wfc(i)),abs(fft_wfc(i))
            end do
            do i=1,molecule%mesh%nactive/2
               write(1,*) i*dk,dreal(fft_wfc(i)),dimag(fft_wfc(i)),abs(fft_wfc(i))
            end do
            close(1)
            

            
            open(unit=1,file="wfc.dat",form='formatted',status='unknown')
            do i=1,molecule%mesh%nactive-1
               write(1,*) i*molecule%mesh%dx,dreal(wfc(i,newt)),dimag(wfc(i,newt)),abs(wfc(i,newt))
            end do
            close(1)
            call  save_agr(idxmov)
         end if
         
         itmp=oldt ; oldt=newt ; newt=itmp
         
       end do
       
       
    end subroutine tdse
  end module tdse_mod
