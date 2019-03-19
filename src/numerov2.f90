module numerov_mod_dev
  use global
  use poten
  use mesh_mod
  use conjugate_gradient_mod
  implicit none
contains
    ! --------------------------------------------------------------------------------------
  !
  !             Numerov()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov(molecule,cvg,param)
    implicit none
    type(t_molecule)::molecule
    type(t_cvg)::cvg
    type(t_param)::param


    integer :: i,j,n_nodes_wanted,idxwfc,nwfc
    double precision,allocatable::r(:)
    integer::iloop
    character (len=1024) :: filename
!    double precision::EHartree,EX
    double precision,parameter::pi=3.1415926


    ! ---------------------------------------------------
    !
    ! building the radial mesh
    !
    ! ---------------------------------------------------
    allocate(r(molecule%mesh%nactive))
    do i=1,molecule%mesh%nactive
       r(i)=i*molecule%mesh%dx
    end do
    print *,"# numerov > dx=",molecule%mesh%dx," first point @ r=",r(1)
    ! V=rR=r^(l+1)*summation
    iloop=1
    cvg%cvg=.FALSE.
    cvg%total_nrj%dnrj=2*cvg%ETA
    molecule%pot%hartree=0.0
    molecule%pot%Vx=0.0

    ! ---------------------------------------------------------------------------
    !
    !                Main loop of numerov()
    !
    ! ---------------------------------------------------------------------------
    do while ((iloop.le.param%loopmax).and.(.not.cvg%cvg))
       print *,'Numerov> ----------------------------------------------------------------------------'
       print *,'Numerov> iloop=',iloop
       ! ================================
       !  updating the potential
       ! ================================
       do i=1,molecule%mesh%nactive
          molecule%pot%tot(i)=-param%Z/r(i)+0.5*param%lorb*(param%lorb+1)/r(i)**2
       end do
       if(param%hartree) then
          do i=1,molecule%mesh%nactive
             molecule%pot%tot(i)=molecule%pot%tot(i)+molecule%pot%hartree(i)
          end do
       end if
       if(param%exchange) then
          do i=1,molecule%mesh%nactive
             molecule%pot%tot(i)=molecule%pot%tot(i)+molecule%pot%Vx(i)
          end do
       end if
       
       open(unit=1,file='pot.dat',form='formatted',status='unknown')
       do i=2,molecule%mesh%nactive
          write(1,*) r(i),molecule%pot%tot(i),-param%Z/r(i),0.5*param%lorb*(param%lorb+1)/r(i)**2,&
               molecule%pot%hartree(i),molecule%pot%Vx(i)
       end do
       close(1)
       print *,'Numerov> Potential extrema ',minval(molecule%pot%tot),maxval(molecule%pot%tot)


       ! ==============================================
       ! computing the radial part of the wave function for different nodes (n)
       !    --> calling numerov_step()
       !==============================================
       nwfc=cvg%nvec_to_cvg
       print *,'Numerov> ',nwfc,' wavefunction to converge'
       do i=1,nwfc
          n_nodes_wanted=i       
          idxwfc=i
          call numerov_step(n_nodes_wanted,molecule,r,cvg,&
               param,idxwfc,molecule%pot,molecule%mesh)
       end do
       
       ! =========================================
       ! from u(r) -> the electron density rho(r) we need to compute
       !        the Exchange potential
       ! =========================================
       molecule%rho=0.0
       do j=1,param%noccstate
          do i=1,molecule%mesh%nactive
             molecule%rho(i)=-(molecule%wf%wfc(i,j)/r(i))**2
          end do
       end do
       !==========================================
       ! computing the exchange potential and the exchange nrj
       !==========================================
       if(param%exchange) then
          do i=1,molecule%mesh%nactive
             molecule%pot%Vx(i)=-(3.0*molecule%wf%wfc(i,1)**2/(2*(pi*r(i))**2))**(1.0/3.0)
          end do
          molecule%pot%EX=  0.0
          do i=1,molecule%mesh%nactive-1
             molecule%pot%EX=molecule%pot%EX+&
                  0.5*molecule%mesh%dx*(&
                  molecule%pot%Vx(i)*molecule%wf%wfc(i,1)**2+&
                  molecule%pot%Vx(i+1)*molecule%wf%wfc(i+1,1)**2)
          end do
       end if
       
       !         call serie(molecule,mesh,r)
       !        open(unit=1,file=trim(param%prefix)//'/rho.dat',form='formatted',status='unknown')
       !       do i=1,mesh%nactive
       !         write(1,*) r(i),molecule%rho(i)/r(i)
       !     end do
       !    close(1)
       
       !==========================================
       ! computing the Hartree potential -> calling hartree_cg()
       !==========================================
       if(param%hartree) then
          do i=1,molecule%mesh%nactive
             molecule%pot%hartree(i)=molecule%pot%hartree(i)*r(i)
          end do

          call Hartree_cg(molecule,molecule%mesh,r,molecule%pot)
          do i=1,molecule%mesh%nactive
             molecule%pot%hartree(i)=molecule%pot%hartree(i)/r(i)
          end do
          ! Hartree energy
          molecule%pot%EHartree=  0.0
          do i=1,molecule%mesh%nactive-1
             molecule%pot%EHartree=molecule%pot%EHartree+&
                  0.5*molecule%mesh%dx*(&
                  molecule%pot%hartree(i)*molecule%wf%wfc(i,1)**2+&
                  molecule%pot%hartree(i+1)*molecule%wf%wfc(i+1,1)**2)
          end do
       end if
       ! ==============================
       ! check the convergence
       ! ==============================
       if(iloop.eq.1) then
          cvg%total_nrj%last=2*molecule%wf%eps(1)
          if(param%hartree) cvg%total_nrj%last=cvg%total_nrj%last-molecule%pot%EHartree
          if(param%exchange) cvg%total_nrj%last=cvg%total_nrj%last+0.5*molecule%pot%EX
          
          cvg%total_nrj%previous=cvg%total_nrj%last
          open(unit=1,file=param%filenrj,form='formatted',status='unknown',access='append')
          write(1,'(A4,A16,A16,A16,A16,A16)') "#loop","ev(1)","EHartree","EX","total_nrj","dnrj";
          close(1)
       else
          cvg%total_nrj%previous=cvg%total_nrj%last
          
          cvg%total_nrj%last=2*molecule%wf%eps(1)
          if(param%hartree) cvg%total_nrj%last=cvg%total_nrj%last-molecule%pot%EHartree
          if(param%exchange) cvg%total_nrj%last=cvg%total_nrj%last+0.5*molecule%pot%EX
          
          
          cvg%total_nrj%dnrj=cvg%total_nrj%last-cvg%total_nrj%previous
       end if
       
       if(abs(cvg%total_nrj%dnrj).lt.cvg%ETA) cvg%cvg=.TRUE.
       open(unit=1,file=param%filenrj,form='formatted',status='unknown',access='append')
       write(1,'(I4,F16.8,F16.8,F16.8,F16.8,E16.8)') iloop,molecule%wf%eps(1),&
            molecule%pot%EHartree,molecule%pot%EX,&
            cvg%total_nrj%last,cvg%total_nrj%dnrj
       close(1)
       
       
       do i=1,nwfc
          print *,'Numerov - ev > ',molecule%wf%eps(i),' Ha=',27.211*molecule%wf%eps(i),' eV'
       end do
       
       
       write(filename,'(a,a,i0,a)') trim(param%prefix),'/wfc',iloop,'.dat'
       open(unit=1,file=filename,form='formatted',status='unknown')
       do i=1,molecule%mesh%nactive
          write(1,*) r(i),(molecule%wf%wfc(i,j),j=1,nwfc)
       end do
       close(1)
       
       !       write(filename,'(a,a,i0,a)') trim(param%prefix),'/potential',iloop,'.dat'
       write(filename,'(a,a)') trim(param%prefix),'/potential.dat'
       open(unit=1,file=filename,form='formatted',status='unknown')
       do i=1,molecule%mesh%nactive
          write(1,*) r(i),molecule%pot%tot(i),molecule%pot%hartree(i),&
               molecule%pot%VX(i),molecule%rho(i)
       end do
       close(1)
       
       iloop=iloop+1
    end do
    
  end subroutine numerov
  ! --------------------------------------------------------------------------------------
  !
  !             numerov_step()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov_step(n_nodes_wanted,molecule,r,cvg,param,idxwfc,pot,mesh)
    integer :: n_nodes_wanted
    type(t_mesh)::mesh
    type(t_potential)::pot
    type(t_param)::param
    type(t_cvg)::cvg
    double precision::r(:)
    type(t_molecule)::molecule
    integer :: iloop
    integer::n_nodes,idxwfc
    double precision::eps
    double precision,allocatable::Q(:),Vin(:),Vout(:),sqrd
    integer :: impt,i
    double precision::emin,emax,dVin,dVout,Iout,Iin,deps,epsold,facsign,eta
    logical,parameter :: outward=.TRUE.,inward=.FALSE.
    ! first we search an eigenenergy close to the soution
    ! by considering the number of nodes of the wavefunction
    allocate(Q(molecule%mesh%nactive))
    allocate(Vin(molecule%mesh%nactive))
    allocate(Vout(molecule%mesh%nactive))
    facsign=(0.5*mod(n_nodes_wanted,2)-1)
    sqrd=mesh%dx**2
    
    emax=min(0.0,maxval(pot%tot))
    emin=1e10
    print *,"====================================================================="
    print *,"=                           Numerov_step                            ="
    print *,"====================================================================="
    print *,' Numerov_step> Potential @ infinite=',emax,maxval(pot%tot)    
    print *," Numerov_step> Searching a wfc with ",n_nodes_wanted," nodes" 

    eps=pot%tot(molecule%mesh%nactive/2)
    n_nodes=-1
    
    ! two steps
    !  1 - first, starting from the highest possible energy, we decrease it until
    !       we search an eigenvalue for which the wfc has n_node_wanted.
    !           Then eps corresponds to the a lower bound for the energy
    iloop=1
    do while       ((iloop.le.param%loopmax).and.(.not.(n_nodes.eq.n_nodes_wanted)))
       call compute_Q(Q,molecule%mesh%nactive,eps,r,pot)
       Vout(1)=0.001
       call numerov_integrate(outward,Q,Vout,molecule%mesh%nactive,sqrd)
       n_nodes=count_nodes(Vout,molecule%mesh%nactive)

       write(*,'(A,I0,A,F8.4,A,I0)',advance='no') " Numerov_step> phase 1 > eps(",iloop,")=",eps," n_node(s)=",n_nodes
       if(n_nodes.gt.n_nodes_wanted) then
          emax=eps
          facsign=eps-1.0
       else if(n_nodes.le.n_nodes_wanted) then
          emin=eps
          facsign=eps+1.0
       end if
       if(emin.lt.emax) then
          eps=0.5*(emin+emax)
       else
          eps=facsign
       end if
       write(*,'(A,F10.4,A,F10.4,A,F8.4)') ' [emin,emax]=[',emin,',',emax,']  new eps=',eps
       iloop=iloop+1
       if(iloop.eq.param%loopmax) then
          print *,"ERROR in numerov_step"
          stop
       end if
       
    end do
    write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(1) eps(',iloop,')=',eps,&
         ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
    !  2 - then starting from emin, we increase it until we reach an wfc with n_node_wanted+1
    !           Then eps corresponds to the a higher bound for the energy
    iloop=1
    do while           ((iloop.le.param%loopmax).and.(.not.(n_nodes.eq.(n_nodes_wanted+1))))
       eps=0.5*(emin+emax)
       print *,'phase 2 > eps =',eps
       if(eps.ge.0.0) stop
       call compute_Q(Q,molecule%mesh%nactive,eps,r,pot)
       Vout(1)=0.001
       call numerov_integrate(outward,Q,Vout,molecule%mesh%nactive,sqrd)

       open(unit=1,file="Vout.dat",form='formatted',status='unknown')
       do i=1,mesh%nactive
          write(1,*) r(i),Vout(i)
       end do
       close(1)




       n_nodes=count_nodes(Vout,molecule%mesh%nactive)
!       write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(2) eps(',iloop,')=',eps,&
!            ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
       if(n_nodes.gt.(n_nodes_wanted+1)) then
          emax=eps
          facsign=eps-1.0
       else if(n_nodes.lt.(n_nodes_wanted+1)) then
          emin=eps
          facsign=eps+1.0
       end if
       if(emin.eq.emax) then
          print *,'!!!! ERROR in numerov_step(): emin=emax=',emin
          stop
          end if
!       if(emin.lt.emax) then
!          eps=0.5*(emin+emax)
!       else
!          eps=facsign
!       end if

       iloop=iloop+1
       if(iloop.eq.param%loopmax) then
          print *,"ERROR in numerov_step"
          stop
       end if

    end do


!    write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(3) eps(',iloop,')=',eps,&
!         ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'
!    stop

    ! Energy levels of Hydrogen atom
    ! n=1 -13.59
    ! n=2   -3.40
    ! n=3  -1.51
    ! n=4  -0.85
    ! n=5  -0.54

    eta=2*cvg%ETA
    iloop=1

    do while((iloop.le.param%loopmax).and.(abs(eta).gt.cvg%ETA))
       eps=0.5*(emin+emax)
       !print *,"new eps=",eps
       call compute_Q(Q,molecule%mesh%nactive,eps,r,pot)
       impt=matching_point(Q,molecule%mesh%nactive)
       Vout(1)=0.001
       call numerov_integrate(outward,Q,Vout,molecule%mesh%nactive,sqrd)
       Vin(molecule%mesh%nactive)=0.001
       call numerov_integrate(inward,Q,Vin,molecule%mesh%nactive,sqrd)
       n_nodes=count_nodes(Vout,molecule%mesh%nactive)

        Iout=0.0
        do i=1,impt-1
           Iout=Iout+0.5*mesh%dx*(Vout(i)**2+Vout(i+1)**2)
        end do
        Iin=0.0
        do i=impt,molecule%mesh%nactive-1
           Iin=Iin+0.5*mesh%dx*(Vin(i)**2+Vin(i+1)**2)
        end do
        dVin=0.5*(Vin(impt+1)-Vin(impt-1))/mesh%dx
        dVout=0.5*(Vout(impt+1)-Vout(impt-1))/mesh%dx
        deps=dVin/Vin(impt)-dVout/Vout(impt)
        deps=deps/(Iout/Vout(impt)**2+Iin/Vin(impt)**2)


       if(iloop.gt.1) eta=eps-epsold
       epsold=eps
       !write(*,'(A,I4,A,F8.4,A,I4,A,F8.4,A,F8.4,A)') '(2) eps(',iloop,')=',eps,&
       !     ' number of node(s)=',n_nodes,' [emin,emax]=[',emin,',',emax,']'

       !      write(*,'(A,I4,A,F8.4,A,F8.4,A,F8.4,A,E12.6,A,E12.6,A,I4,A,E12.6,A,E12.6)') 'epsmax(',iloop,')=',eps,&
       !          ' [emin,emax]=[',emin,',',emax,'] eta=',eta,' deps=',deps,' impt=',impt,&
       !         ' dVout=',dVout/Vout(impt),' dVin=',dVin/Vin(impt)


       if(n_nodes.le.(n_nodes_wanted)) then
          emin=eps
       else
          emax=eps
       end if



       iloop=iloop+1
       if(iloop.eq.param%loopmax) then
          print *,"ERROR in numerov_step"
          stop
       end if
    end do
    


    molecule%wf%eps(idxwfc)=eps
    do i=1,impt
       molecule%wf%wfc(i,idxwfc)=Vout(i)/Vout(impt)
    end do
    do i=impt+1,mesh%nactive
       molecule%wf%wfc(i,idxwfc)=Vin(i)/Vin(impt)
    end do
    call norm(mesh,molecule%wf%wfc(:,idxwfc))
    print *,'Numerov_step> wfc(',idxwfc,')=',sqrt(simpson(mesh,molecule%wf%wfc(:,idxwfc)))
    
    deallocate(Q)
    deallocate(Vin)
    deallocate(Vout)
      

    end subroutine numerov_step
      ! --------------------------------------------------------------------------------------
    !
    !             compute_Q()
    !
    ! --------------------------------------------------------------------------------------
    subroutine compute_Q(Q,N,eps,r,pot)
      double precision::Q(:),r(:),eps
      type(t_potential)::pot
    integer::N,i
!    Q(1)=10000.0
    do i=1,N
       Q(i)=2*(eps-pot%tot(i))
    end do
    open(unit=1,file="Q.dat",form='formatted',status='unknown')
    do i=1,N
       write(1,*) r(i),Q(i)
    end do
    close(1)
    if(maxval(Q).lt.0.0) then
       print *,'!!!! ERROR in compute_Q() !!!!'
       print *,'Q<=',maxval(Q),'<0'
       stop
    end if
    
  end subroutine compute_Q
  ! --------------------------------------------------------------------------------------
  !
  !             Numerov_integrate()
  !
  ! --------------------------------------------------------------------------------------
  subroutine numerov_integrate(outward,Q,y,N,sqrd)
    logical :: outward
    double precision::Q(:),y(:)
    double precision::t(3),sqrd
    integer :: N,i

    if(outward) then
!       print *,'outward'
       t(1)=1.0+sqrd*Q(2)/12.0
       t(2)=2*(1.0-5.0*sqrd*Q(1)/12.0)
       y(2)=t(2)*y(1)/t(1)
       do i=2,N-1
          t(1)=1.0+sqrd*Q(i+1)/12.0
          t(2)=2*(1.0-5.0*sqrd*Q(i)/12.0)
          t(3)=1.0+sqrd*Q(i-1)/12.0
          y(i+1)=(t(2)*y(i)-t(3)*y(i-1))/t(1)
       end do
    else
!       print *,'inward'
       t(1)=1.0+sqrd*Q(N-1)/12.0
       t(2)=2*(1.0-5.0*sqrd*Q(N)/12.0)
       y(N-1)=t(2)*y(N)/t(1)
       do i=N-1,2,-1
          t(1)=1.0+sqrd*Q(i-1)/12.0
          t(2)=2*(1.0-5.0*sqrd*Q(i)/12.0)
          t(3)=1.0+sqrd*Q(i+1)/12.0
          y(i-1)=(t(2)*y(i)-t(3)*y(i+1))/t(1)
       end do
    end if
  end subroutine numerov_integrate

    ! --------------------------------------------------------------------------------------
    !
    !             matching_point()
    !
    ! --------------------------------------------------------------------------------------
    function matching_point(Q,N)
!      function matching_point(Q,N,eps,r,pot)
      implicit none
      integer::N
      double precision::Q(:)!,r(:),eps
      !type(t_potential)::pot
      integer::matching_point,i
      !      call compute_Q(Q,N,eps,r,pot)
      matching_point=-1
      do i=2,N-1
         if((Q(i)*Q(i+1)).le.0) then
            !     print *,i,pot%tot(i),Q(i),Q(i+1)
            if(matching_point.lt.0)          matching_point=i
         end if
      end do
    end function matching_point
  ! --------------------------------------------------------------------------------------
  !
  !             n_nodes()
  !
  ! --------------------------------------------------------------------------------------
  function count_nodes(V,N)
    double precision::V(:)
    integer::N,i
    integer::count_nodes
    count_nodes=1
    do i=1,N-1
       if((V(i)*V(i+1)).le.0) count_nodes=count_nodes+1
    end do
  end function count_nodes
  ! --------------------------------------------------------------------------------------
  !
  !             serie()
  ! It allows us to compute the value of R(@r=0) from u(0)/0
  ! --------------------------------------------------------------------------------------
  subroutine serie(molecule,mesh,r)
    implicit none
    type(t_mesh)::mesh
    type(t_molecule)::molecule
    integer, parameter :: nmax=5
    double precision, dimension(nmax) :: a
    double precision:: r(:)
    double precision, dimension(nmax,nmax) :: b
    double precision :: som
    integer :: i, j,info, lda, ldb, nrhs
    integer, dimension(nmax) :: ipiv
      
      

    do i=1,nmax
       do j=1,nmax
          b(i,j)=r(1+i)**j  
       end do
       a(i)=molecule%wf%wfc(1+i,1)
      end do
      nrhs = 1 ! number of right hand sides in b
      lda = nmax  ! leading dimension of a
      ldb = nmax  ! leading dimension of b

      call dgesv(nmax, nrhs, b, lda, ipiv, a, ldb, info)
!      print *,'a=',a

!      print *,'-----------------------------------------'
!      print *,'r ','u ','som'
!      print *,'-----------------------------------------'
      do j=1,nmax
         som=0
         do i=1,nmax
            som=som+a(i)*(r(1+j)**i)
         end do
!         print *,r(1+j),molecule%wf%wfc(1+j,1),som
      end do


!    open(unit=1,file='R.dat',form='formatted',status='unknown')
!    write(1,*) 0.0,a(1)
!    do i=2,mesh%nactive
!       write(1,*) r(i),molecule%wf%wfc(i,1)/r(i)
!    end do
!    close(1)



      ! Note: the solution is returned in b
      ! and a has been changed.

    end subroutine serie
    ! --------------------------------------------------------------------------------------
  !
  !             Hartree_cg() (conjugate gradient)
  ! 
  ! --------------------------------------------------------------------------------------
  subroutine Hartree_cg(molecule,mesh,r,pot)
    implicit none
    type(t_potential)::pot
    double precision::r(:)
    type(t_molecule)::molecule
    type(t_mesh)::mesh

    integer :: i
    double precision :: charge_inf
!    double precision,allocatable::source(:)
 !   double precision,allocatable::U(:)
    double precision,allocatable::b(:)

    
    
    charge_inf=0.5*mesh%dx*(molecule%wf%wfc(1,1)**2+molecule%wf%wfc(mesh%nactive,1)**2)
    do i=1,mesh%nactive-1
       charge_inf=charge_inf+0.5*mesh%dx*(molecule%wf%wfc(i,1)**2+molecule%wf%wfc(i+1,1)**2)
    end do
!    charge_inf =2*charge_inf 
    print *,"charge @ infinite=",charge_inf 



    allocate(b(mesh%nactive))
    do i=1,mesh%nactive
!       b(i)=-molecule%rho(i)
       b(i)=-molecule%wf%wfc(i,1)**2/r(i)
    end do


    do i=1,mesh%nactive
       pot%hartree(i)=pot%hartree(i)*r(i)
    end do
    b(mesh%nactive)=b(mesh%nactive)-charge_inf/mesh%dx**2
    
    call Conjugate_gradient_3D(-b,pot%hartree,mesh%nactive,mesh%dx,mesh)    

    open(unit=1,file='hartree.dat',form='formatted',status='unknown')
    write(1,*) r(1)-mesh%dx,0.0
    do i=1,mesh%nactive
!       pot%hartree(i)=pot%hartree(i)/r(i)
       write(1,*) r(i),pot%hartree(i)
    end do
    write(1,*) r(mesh%nactive)+mesh%dx,charge_inf  !/(r(mesh%nactive)+mesh%dx)
    close(1)
    
    deallocate(b)
  end subroutine Hartree_cg


  ! --------------------------------------------------------------------------------------
  !
  !             Hartree_sd() (steepest descent)
  ! 
  ! --------------------------------------------------------------------------------------
  ! subroutine Hartree2()
  !   implicit none
  !   integer ::n,i,j
  !   double precision::alpha,dr,rmax
  !   double precision,allocatable::source(:),b(:)
  !   double precision,allocatable::U(:),grad(:),y(:)
  !   double precision,allocatable::r(:)
  !   integer,parameter :: seed = 86456
  !   double precision, external :: ddot
  !   call srand(seed)
  !   n=1000
  !   rmax=50.0
  !   dr=rmax/(n-1)
  !   allocate(source(n))
  !   allocate(r(n))
  !   allocate(b(n))
  !   allocate(U(n))
  !   allocate(grad(n))
  !   allocate(y(n))
  !   do i=1,n
  !      r(i)=(i-1)*dr
  !   end do
  !   do i=1,n
  !      source(i)=-4.0*r(i)*exp(-2.0*r(i))
  !   end do

    
  !   U(1)=0
  !   do i=2,n-1
  !      U(i)=rand()
  !   end do
  !   U(n)=1.0

  !   do i=2,n-1
  !      b(i)=source(i)
  !   end do
  !   b(2)=source(2)-U(1)/dr**2
  !   b(n-1)=source(n-1)-U(n)/dr**2

  !   grad(1)=0.0
  !   grad(n)=0.0
  !   y(1)=0.0
  !   y(n)=0.0
  !   do j=1,1000000
  !      grad(2)=(U(3)-2*U(2))/dr**2-b(2)
  !      do i=3,n-2
  !         grad(i)=(U(i+1)-2*U(i)+U(i-1))/dr**2-b(i)
  !      end do
  !      grad(n-1)=(-2*U(n-1)+U(n-2))/dr**2-b(n-1)


  !      y(2)=(grad(3)-2*grad(2))/dr**2
  !      do i=3,n-2
  !         y(i)=(grad(i+1)-2*grad(i)+grad(i-1))/dr**2
  !      end do
  !      y(n-1)=(-2*grad(n-1)+grad(n-2))/dr**2
       
  !      alpha=ddot(n,grad,1,grad,1)/ddot(n,grad,1,y,1)
  !      print *,'alpha=',alpha ,ddot(n,grad,1,y,1),ddot(n,grad,1,grad,1)
       
  !      !    call dscal(mesh%nactive-2,-alpha,grad,1)
  !      call daxpy(n,-alpha,grad,1,U,1)

       
  !   end do
  !   U(1)=0.0
  !   U(n)=1.0


  !   open(unit=1,file='b.dat',form='formatted',status='unknown')
  !   do i=1,n
  !      write(1,*) r(i),source(i),-(r(i)+1)*exp(-2.0*r(i))+1,U(i)
  !   end do
    
  ! end subroutine Hartree2
  subroutine Hartree_sd(molecule,mesh,r)
    implicit none
    double precision::r(:)
    type(t_molecule)::molecule
    type(t_mesh)::mesh
    double precision,allocatable::source(:),U(:),grad(:),y(:),b(:)
    double precision::alpha,q
    integer::i,j
    double precision, external :: ddot
    integer,parameter :: seed = 86456
    character (len=1024) :: filesave

    allocate(source(molecule%mesh%nactive))
    allocate(U(molecule%mesh%nactive))
    allocate(b(molecule%mesh%nactive))
    allocate(grad(molecule%mesh%nactive))
    allocate(y(molecule%mesh%nactive))
    call srand(seed)

    q=0.0
    do i=1,mesh%nactive-1
       q=q+0.5*mesh%dx*(molecule%wf%wfc(i,1)**2+molecule%wf%wfc(i+1,1)**2)
    end do
    print *,"q=",q 

    do i=2,mesh%nactive-1
       U(i)=rand()
    end do
    U(1)=0.0
    U(mesh%nactive)=q


    do i=1,mesh%nactive
       source(i)=-molecule%wf%wfc(i+1,1)**2/r(i+1)
    end do


    do i=1,mesh%nactive
       b(i)=source(i)
    end do
    b(2)=b(2)-U(1)/mesh%dx**2
    b(mesh%nactive-1)=b(mesh%nactive-1)-U(mesh%nactive)/mesh%dx**2



    
    grad(1)=0.0
    grad(mesh%nactive)=0.0
    y(1)=0.0
    y(mesh%nactive)=0.0

    
    do j=1,100000000
       grad(2)=(U(3)-2*U(2))/mesh%dx**2-b(2)
       do i=3,mesh%nactive-2
          grad(i)=(U(i+1)-2*U(i)+U(i-1))/mesh%dx**2-b(i)
       end do
       grad(mesh%nactive-1)=(-2*U(mesh%nactive-1)+U(mesh%nactive-2))/mesh%dx**2-b(mesh%nactive-1)
       
       y(2)=(grad(3)-2*grad(2))/mesh%dx**2
       do i=3,mesh%nactive-2
          y(i)=(grad(i+1)-2*grad(i)+grad(i-1))/mesh%dx**2
       end do
       y(mesh%nactive-1)=(-2*grad(mesh%nactive-1)+grad(mesh%nactive-2))/mesh%dx**2
       
       alpha=ddot(mesh%nactive,grad,1,grad,1)/ddot(mesh%nactive,grad,1,y,1)
       print *,'alpha=',alpha,ddot(mesh%nactive,grad,1,y,1),ddot(mesh%nactive,grad,1,grad,1)
       
       !    call dscal(mesh%nactive-2,-alpha,grad,1)
       call daxpy(mesh%nactive,-alpha,grad,1,U,1)

       if(mod(j,10000).eq.0) then
          write(filesave,'(a,i0,a)') 'b',j,'.dat'
          open(unit=1,file=filesave,form='formatted',status='unknown')
          do i=1,mesh%nactive
             write(1,*) r(i),source(i-1),U(i)
          end do
       end if
       
    end do


!    open(unit=1,file='b.dat',form='formatted',status='unknown')
!    do i=1,mesh%nactive
!       write(1,*) r(i),source(i-1),U(i)
!    end do

!    write(1,*) r(mesh%nactive),source(i-1),x(i)
!    close(1)

    
    deallocate(y)
    deallocate(grad)
    deallocate(source)
    deallocate(U)
  end subroutine Hartree_sd
end module numerov_mod_dev
