module param_mod
  use time_tracking
  use global
  use mesh_mod
  use poten
  use davidson_mod
  use tools
  use experiment
  use FFT_mod
  use tdse_mod
  use numerov_mod_dev
  use pseudopotential
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              parse_line()
  !
  ! --------------------------------------------------------------------------------------
  subroutine parse_line(param,field,nfield,end_loop,nmol,molecule,time_spent,syst)
    implicit none
    type(t_system)::syst
    type(t_param)::param
    character (len=32)::field(32)
    logical::end_loop
    integer::nmol
    type(t_molecule),allocatable:: molecule(:)
    type(t_molecule),allocatable:: junk(:) 
    type(t_time) :: time_spent   
    integer::i,nfield,ifield,idxwfc,idxmol,j
    double complex, allocatable::SH(:,:)
    integer::l,m
    double precision,allocatable::junk_wfc(:)
    double complex,allocatable::tdse_wfc(:)
    double precision,allocatable::coeff(:)
    double precision::norm,k0,dk

    double precision::r0,sig,Intens,x,y,z,x0,y0,z0
    double precision::r1(4),I2(3)
    double precision,allocatable::r2(:,:)
    double precision::theta,phi
    double precision, external :: ddot
    character (len=1024) :: filename
    double complex,allocatable::in_wfc(:)
    double complex,allocatable::fft_wfc(:)
    type(t_pseudo),allocatable :: pp(:)

    select case (trim(field(1)))
    case("iprint_level") 
       read(field(2),*) syst%iprint_level
    case("box_radius") 
       read(field(2),*) param%box%radius
    case("box_width") 
       read(field(2),*) param%box%width
    case("box_shape") 
       read(field(2),*) param%box%shape
    case("box_center") 
       read(field(2),*) param%box%center(1)
       read(field(3),*) param%box%center(2)
       read(field(4),*) param%box%center(3)
    case("dimension") 
       read(field(2),*) param%dim
    case("ETA") 
       read(field(2),*) param%ETA
    case("exchange") 
       read(field(2),*) param%exchange
    case("extrap_add") 
       read(field(2),*) param%extrap_add
    case("extrapol") 
       read(field(2),*) param%extrapol
    case("hartree") 
       read(field(2),*) param%hartree
    case("init_wf") 
       read(field(2),*) param%init_wf
    case("perturb_intensity") 
       read(field(2),*) param%perturb%Intensity
    case("perturb_sigma") 
       read(field(2),*) param%perturb%sigma
    case("perturb_shape") 
       read(field(2),*) param%perturb%shape
    case("perturb_location") 
       read(field(2),*) param%perturb%location(1)
       read(field(3),*) param%perturb%location(2)
       read(field(4),*) param%perturb%location(3)
    case("loopmax") 
       read(field(2),*) param%loopmax
    case("noccstate") 
       read(field(2),*) param%noccstate
    case("nvecmax") 
       read(field(2),*) param%nvecmax
    case("nvecmin") 
       read(field(2),*) param%nvecmin
    case("nvec_to_cvg") 
       read(field(2),*) param%nvec_to_cvg
       deallocate(param%list_idx_to_cvg)
       allocate(param%list_idx_to_cvg(param%nvec_to_cvg))
       do i=1,param%nvec_to_cvg
          param%list_idx_to_cvg(i)=i
       end do
       deallocate(param%occupation)
       allocate(param%occupation(param%nvec_to_cvg))
       do i=1,param%nvec_to_cvg
          param%occupation(i)=1.0
       end do
    case("Nx") 
       read(field(2),*) param%Nx
    case("occupation") 
       do i=1,param%nvec_to_cvg
          read(field(1+i),*) param%occupation(i)
       end do
    case("prefix") 
       read(field(2),*) param%prefix
    case("restart") 
       read(field(2),*) param%restart
    case("scheme") 
       read(field(2),*) param%scheme
    case("Zato") 
       read(field(2),*) param%Z
    case("lorb") 
       read(field(2),*) param%lorb
       ! ---------------------------------------------------------------
       !
       !                    CREATE MOLECULE
       !
       ! --------------------------------------------------------------
    case("create")
       read(field(3),*) idxmol       
       molecule(idxmol)%mesh%box%center(1)=molecule(idxmol)%mesh%box%center(1)*&
            molecule(idxmol)%mesh%box%width
       molecule(idxmol)%mesh%box%center(2)=molecule(idxmol)%mesh%box%center(2)*&
            molecule(idxmol)%mesh%box%width
       molecule(idxmol)%mesh%box%center(3)=molecule(idxmol)%mesh%box%center(3)*&
            molecule(idxmol)%mesh%box%width
       call new_mesh(molecule(idxmol)%mesh)
       call init_pot(molecule(idxmol))

       print *,molecule(idxmol)%wf%nwfc
       allocate(molecule(idxmol)%wf%eps(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%epsprev(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%deps(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%wfc(0:molecule(idxmol)%mesh%nactive,molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%n(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%l(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%m(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%wf%occ(molecule(idxmol)%wf%nwfc))
       allocate(molecule(idxmol)%rho(molecule(idxmol)%mesh%Ntot))


       allocate(molecule(idxmol)%cvg%wfc(molecule(idxmol)%wf%nwfc))
       do i=1,molecule(idxmol)%cvg%nwfc
          molecule(idxmol)%cvg%wfc%cvg=.FALSE.
       end do
       allocate(molecule(idxmol)%cvg%list_idx(molecule(idxmol)%cvg%nvec_to_cvg))
       do i=1,molecule(idxmol)%cvg%nvec_to_cvg
          molecule(idxmol)%cvg%list_idx(i)=param%list_idx_to_cvg(i)
       end do
       
       allocate(molecule(idxmol)%numerov%Q(molecule(idxmol)%mesh%nactive))
       allocate(molecule(idxmol)%numerov%Vout(molecule(idxmol)%mesh%nactive))
       allocate(molecule(idxmol)%numerov%Vin(molecule(idxmol)%mesh%nactive))
       allocate(molecule(idxmol)%numerov%r(molecule(idxmol)%mesh%nactive))
       allocate(molecule(idxmol)%numerov%rho(molecule(idxmol)%mesh%nactive))
       allocate(molecule(idxmol)%numerov%rhoold(molecule(idxmol)%mesh%nactive))
       ! ---------------------------------------------------------------
       !
       !                    SETTING MOLECULE
       !
       ! --------------------------------------------------------------
    case("set")
       print *,syst%nmol," molecule(s) are defined"
       if(.not.(allocated(molecule))) then
          nmol=1
          allocate(molecule(nmol))
       end if
       idxmol=1       
       if(nfield.gt.1) then
          read(field(3),*) idxmol
          do i=4,nfield
             select case (trim(field(i)))
             case("dimension") 
                read(field(i+1),*) molecule(idxmol)%mesh%dim
                print *,"Dimension= ",molecule(idxmol)%mesh%dim
             case("width")
                read(field(i+1),*) molecule(idxmol)%mesh%box%width
                print *,"width= ",molecule(idxmol)%mesh%box%width
             case("shape")
                read(field(i+1),*) molecule(idxmol)%mesh%box%shape
             case("radius")
                read(field(i+1),*) molecule(idxmol)%mesh%box%radius
             case("box_center") 
                read(field(i+1),*) molecule(idxmol)%mesh%box%center(1)
                read(field(i+2),*) molecule(idxmol)%mesh%box%center(2)
                read(field(i+3),*) molecule(idxmol)%mesh%box%center(3)
             case("N")
                read(field(i+1),*) molecule(idxmol)%mesh%Nx
              case("norbital")
                read(field(i+1),*) molecule(idxmol)%wf%nwfc
              case("mixing")
                read(field(i+1),*) molecule(idxmol)%mixing
             case("nvec_to_cvg")
                read(field(i+1),*) molecule(idxmol)%cvg%nvec_to_cvg
             case("ETA")
                read(field(i+1),*) molecule(idxmol)%cvg%ETA
             case("nloopmax")
                read(field(i+1),*) molecule(idxmol)%param%loopmax
             case("occupation") 
                do j=1,molecule(idxmol)%wf%nwfc
                   print *,i,j,field(i+j)
                   read(field(i+j),*) molecule(idxmol)%wf%occ(j)
                end do
             case("tdse_nstep")
                read(field(i+1),*) molecule(idxmol)%param%tdse%nstep
             case("tdse_freq_save")
                read(field(i+1),*) molecule(idxmol)%param%tdse%freq_save
             case("numerov_Z")
                read(field(i+1),*) molecule(idxmol)%numerov%Z
             case("numerov_nmax")
                read(field(i+1),*) molecule(idxmol)%numerov%nmax
             case("hartree") 
                read(field(i+1),*) molecule(idxmol)%param%hartree
             case("exchange") 
                read(field(i+1),*) molecule(idxmol)%param%exchange
             end select
          end do
       end if


    
       !call init_param_molecule(molecule(idxmol),param)

       !call new_molecule(molecule(idxmol),molecule(idxmol)%param)
       !call new_molecule(molecule(idxmol),param)
       !       print *, "# Changing molecule  ",idxmol
       
    case("molecule") 
       select case (trim(field(2)))
       case("new")
          print *,"----------------------------------"
          print *,"     NEW  MOLECULE                "
          print *,"----------------------------------"
          
          nmol=nmol+1
          print *,allocated(molecule)
          if(.not.(allocated(molecule))) then
             nmol=1
             allocate(molecule(nmol))
             call init_param_molecule(molecule(nmol),param)
             call new_molecule(molecule(nmol),molecule(nmol)%param)
             call print_param(molecule(nmol)%param)
             call print_param(param)
             !call exit()
          else
             print *,"Not yet implemented"
             call exit()
             !    allocate(junk(nmol-1))
             !    call move_alloc(molecule,junk)
             !    allocate(molecule(nmol))
             !    call move_alloc(junk,molecule)
          end if
          print *, "# creating a molecule -> ",nmol," molecule(s)"
       end select
    case ("davidson")
       print *,"----------------------------------"
       print *,"          DAVIDSON            "
       print *,"----------------------------------"
       idxmol=1
       if(nfield.lt.1) then
          do i=2,nfield
             if(field(i).eq."molecule") read(field(i+1),*) idxmol
          end do
       end if
       print *,nfield
       !call exit()
       call davidson(molecule(idxmol)%param,molecule(idxmol)%mesh,&
            molecule(idxmol)%cvg,molecule(idxmol),molecule(idxmol)%pot,time_spent)    



       allocate(fft_wfc(molecule(idxmol)%mesh%nactive))
       allocate(in_wfc(molecule(idxmol)%mesh%nactive))
       do j=1,5
          do i=1,molecule(idxmol)%mesh%nactive
             in_wfc(i)=cmplx(molecule(idxmol)%wf%wfc(i,j),0.0)
          end do
          call FFT(in_wfc,fft_wfc,molecule(idxmol)%mesh%nactive)
          
          dk=1.0/(molecule(idxmol)%mesh%nactive*molecule(idxmol)%mesh%dx)
          write(filename,'(a,i0,a)') 'fft',j,'.dat'
          open(unit=1,file=filename,form='formatted',status='unknown')
          do i=molecule(idxmol)%mesh%nactive/2+1,molecule(idxmol)%mesh%nactive
             write(1,*) (i-molecule(idxmol)%mesh%nactive-1)*dk,dreal(fft_wfc(i)),dimag(fft_wfc(i)),abs(fft_wfc(i))
          end do
          do i=1,molecule(idxmol)%mesh%nactive/2
             write(1,*) (i-1)*dk,dreal(fft_wfc(i)),dimag(fft_wfc(i)),abs(fft_wfc(i))
          end do
          close(1)
       end do
       deallocate(fft_wfc)
       deallocate(in_wfc)
!       call exit()
    case("experiment")
       print *,"---------------------------------------------------------------"
       print *,"                        Experimental "
       print *,"---------------------------------------------------------------"
       molecule(nmol)%rho=0
       print *,molecule(nmol)%mesh%box%center,molecule(nmol)%mesh%box%width
       print *,molecule(nmol)%mesh%dx
       sig=0.05*molecule(nmol)%mesh%box%width
       do i=1,molecule(nmol)%mesh%Ntot
          x=0.5*(&
               (molecule(nmol)%mesh%node(i)%q(1)-molecule(nmol)%mesh%box%center(1))**2+&
               (molecule(nmol)%mesh%node(i)%q(2)-molecule(nmol)%mesh%box%center(2))**2+&
               (molecule(nmol)%mesh%node(i)%q(3)-molecule(nmol)%mesh%box%center(3))**2&
          )/sig**2
          molecule(nmol)%rho(i)=exp(-x)
          print *,i,x,y,z,molecule(nmol)%rho(i)
       enddo
       norm=sum(molecule(nmol)%rho)*molecule(nmol)%mesh%dv
       molecule(nmol)%rho=molecule(nmol)%rho/norm
       norm=sum(molecule(nmol)%rho)*molecule(nmol)%mesh%dv
       print *,sum(molecule(nmol)%rho)*molecule(nmol)%mesh%dv
       write(filename,'(a)') 'rho.cube'
       call save_cube_3D(molecule(nmol)%rho,filename,molecule(nmol)%mesh)

       open(unit=1,file="scan.dat",form='formatted',status='unknown')
       y=molecule(nmol)%mesh%box%center(2)
       z=molecule(nmol)%mesh%box%center(3)
       do i=-10,10
          x=molecule(nmol)%mesh%box%center(1)+(i)*molecule(nmol)%mesh%dx
          write(1,*) x,interpolate2(x,y,z,molecule(nmol)%mesh,molecule(nmol)%rho)
       end do
       close(1)

       allocate(r2(molecule(nmol)%mesh%nactive,3))
       do i=1,molecule(nmol)%mesh%nactive
          r2(i,1)=(molecule(nmol)%mesh%node(i)%q(1)-x0)*molecule(nmol)%rho(i)
          r2(i,2)=(molecule(nmol)%mesh%node(i)%q(2)-y0)*molecule(nmol)%rho(i)
          r2(i,3)=(molecule(nmol)%mesh%node(i)%q(3)-z0)*molecule(nmol)%rho(i)
       end do
       I2(1)=sum(r2(:,1))*molecule(nmol)%mesh%dv
       I2(2)=sum(r2(:,2))*molecule(nmol)%mesh%dv
       I2(3)=sum(r2(:,3))*molecule(nmol)%mesh%dv       
!       print *,I2(1),I2(2),I2(3)

       phi=0.0
       open(unit=1,file="phi.dat",form='formatted',status='unknown')
       x0=molecule(nmol)%mesh%box%center(1)
       y0=molecule(nmol)%mesh%box%center(2)
       z0=molecule(nmol)%mesh%box%center(3)       
       do i=0,360
          theta=1.0*i
          x=molecule(nmol)%mesh%box%center(1)+0.5*molecule(nmol)%mesh%box%radius*cos(theta*pi/180.0)
          y=molecule(nmol)%mesh%box%center(2)+0.5*molecule(nmol)%mesh%box%radius*sin(theta*pi/180.0)
          z=molecule(nmol)%mesh%box%center(3)       
          r1(1)=x-x0
          r1(2)=y-y0
          r1(3)=z-z0
          r1(4)=sqrt(r1(1)*r1(1)+r1(2)*r1(2)+r1(3)*r1(3))
          write(1,*) x,y,z,r1(4),norm/r1(4),(r1(1)*I2(1)+r1(2)*I2(2)+r1(3)*I2(3))/r1(4)**3,&
               norm/r1(4)+(r1(1)*I2(1)+r1(2)*I2(2)+r1(3)*I2(3))/r1(4)**3

!          print *,theta*pi/180,sph_harm(0,0,theta*pi/180,phi)
       end do
       close(1)

       deallocate(r2)




       write(filename,'(a)') 'Y00.cube'
       call save_cube_3D(dreal(molecule(nmol)%mesh%multipole%sph_harm_l(0)%m(0)%val),filename,molecule(nmol)%mesh)
       write(filename,'(a)') 'Y10.cube'
       call save_cube_3D(dreal(molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(0)%val),filename,molecule(nmol)%mesh)
       write(filename,'(a)') 'Y1+.cube'
       call save_cube_3D(&
            dimag(molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(-1)%val&
            +molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(1)%val)/sqrt(2.0),&
            filename,molecule(nmol)%mesh)
       write(filename,'(a)') 'Y1-.cube'
       call save_cube_3D(&
            dreal(molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(-1)%val&
            -molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(1)%val)/sqrt(2.0),&
            filename,molecule(nmol)%mesh)
!       do i=1,molecule(nmol)%mesh%Ntot
!          print *,molecule(nmol)%mesh%multipole%sph_harm_l(2)%m(1)%val(i)&
!               +molecule(nmol)%mesh%multipole%sph_harm_l(2)%m(3)%val(i),&
!               molecule(nmol)%mesh%multipole%sph_harm_l(2)%m(2)%val(i)
!       end do
!        print *,"HHHHHHH"
!        do l=0,molecule(nmol)%mesh%multipole%lmax
!           do m=-l,l
!              molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)=0.0
!              do i=1,molecule(nmol)%mesh%nactive
!                 molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)=&
!                      molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)+&
!                      molecule(nmol)%mesh%multipole%sph_harm_l(l)%m(m)%val(i)*&
!                      molecule(nmol)%mesh%multipole%rs(l)%val(i)*molecule(nmol)%rho(i)
!                 !          print *,&
!                 !               molecule(nmol)%mesh%multipole%sph_harm_l(l)%m(m)%val(i),&
!                 !               molecule(nmol)%mesh%multipole%rs(l)%val(i),&
!                 !               molecule(nmol)%rho(i)
!              end do
!              molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)=&
!                   molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)*&
!                   molecule(nmol)%mesh%dv
!           end do
!        end do
!        do i=molecule(nmol)%mesh%nunactive,molecule(nmol)%mesh%Ntot
!           molecule(nmol)%pot%hartree(i)=0.0
!           do l=0,molecule(nmol)%mesh%multipole%lmax
!              do m=-l,l
!                 molecule(nmol)%pot%hartree(i)=&
!                      molecule(nmol)%pot%hartree(i)+&
!                      molecule(nmol)%mesh%multipole%qlm(l)%m(m)%val(1)*&
!                      molecule(nmol)%mesh%multipole%sph_harm_l(l)%m(m)%val(i)/&
!                      ((2*l+1)*molecule(nmol)%mesh%node(i)%r**(l+1))
!              end do
!           end do
!           molecule(nmol)%pot%hartree(i)=4*pi*molecule(nmol)%pot%hartree(i)
! !          print *,i,molecule(nmol)%pot%hartree(i)
!        end do

!        write(filename,'(a)') 'hartree.cube'
!        call save_cube_3D(molecule(nmol)%pot%hartree,filename,molecule(nmol)%mesh)



!        print *,molecule(nmol)%mesh%dv*&
!             sum(molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(1)%val*molecule(nmol)%mesh%multipole%sph_harm_l(1)%m(1)%val)
       
!        call integrate_Yl1m1_Yl2m2(molecule(nmol)%mesh,0,0,0,0)
!        call integrate_Yl1m1_Yl2m2(molecule(nmol)%mesh,1,0,0,0)
!        call integrate_Yl1m1_Yl2m2(molecule(nmol)%mesh,1,0,1,0)
!        call integrate_Yl1m1_Yl2m2(molecule(nmol)%mesh,1,1,0,0)
       
       call exit()
    case ("numerov")
       print *,"----------------------------------"
       print *,"          NUMEROV             "
       print *,"----------------------------------"
       print *,"Dimension= ",molecule(nmol)%mesh%dim
       print *,"Ntot= ",molecule(nmol)%mesh%Ntot
       print *,"dx= ",molecule(nmol)%mesh%dx
       print *,"Z= ",molecule(nmol)%numerov%Z
       print *,"nmax= ",molecule(nmol)%numerov%nmax
       call       numerov_new(molecule(nmol),syst)
       call exit()
    case ("pseudopotential")
       print *,"----------------------------------"
       print *,"          PSEUDOPOTENTIAL             "
       print *,"----------------------------------"
       allocate(pp(1))
       print *,"nfield= ",nfield
       if(nfield.gt.1) then
          do i=2,nfield
             if(field(i).eq."file") then
                read(field(i+1),*) pp(1)%file
       
             end if
          end do
       end if
       print *,'pp file= ',trim(pp(1)%file)
       call read_pp(pp(1))
       !call       numerov_new(molecule(nmol),syst)
       !deallocate(pp)
       call exit()
    case("tdse")
       print *,"---------------------------------------------------------------"
       print *,"          Time-Dependent Schrodinger Equation"
       print *,"---------------------------------------------------------------"
       allocate(tdse_wfc(molecule(nmol)%mesh%nactive))
       allocate(junk_wfc(molecule(nmol)%mesh%nactive))
       tdse_wfc=0.0
       junk_wfc=0.0
       !
       ! - 1 - building the wave packets
       !
       r0=10.
       sig=1.0
       Intens=1.0

       call wave_packet(r0,sig,Intens,molecule(nmol),junk_wfc)


        open(unit=10,file="tdse.dat",form='formatted')
        do i=1,molecule(nmol)%mesh%nactive
           write(10,*) i*molecule(nmol)%mesh%dx,molecule(nmol)%pot%ext(i),&
                junk_wfc(i),&
                (molecule(nmol)%wf%wfc(i,j),j=1,param%nvec_to_cvg)
        end do
        close(10)
       ! !
       ! ! - 2 - projecting the wave packet onto the eigenstates of the potential
       ! !
       ! allocate(coeff(param%nvec_to_cvg))
       ! do  idxwfc=1,param%nvec_to_cvg
       !    coeff(idxwfc)=-molecule(nmol)%mesh%dx*ddot(molecule(nmol)%mesh%nactive,&
       !         junk_wfc,1,&
       !         molecule(nmol)%wf%wfc(:,idxwfc),1)
       !    call daxpy(molecule(nmol)%mesh%nactive,&
       !         coeff(idxwfc),molecule(nmol)%wf%wfc(:,idxwfc),1,&
       !         junk_wfc,1)  ! tdse_wfc+coeff*wfc(idxwfc) -> tdse_wfc
       !    norm=sqrt(molecule(nmol)%mesh%dx*ddot(molecule(nmol)%mesh%nactive,&
       !         junk_wfc,1,&
       !         junk_wfc,1))
          
       !    print *,"(coeff,norm)=",coeff(idxwfc),norm
       ! end do
       ! !
       ! ! saving the residual
       ! !
       ! open(unit=10,file="tdse1.dat",form='formatted')
       ! do i=1,molecule(nmol)%mesh%nactive
       !    write(10,*) i*molecule(nmol)%mesh%dx,junk_wfc(i)
       ! end do
       ! close(10)
       ! !
       ! ! - 3 - building the new wave packet from the projection coefficients
       ! !
       ! junk_wfc=0.0
       ! do  idxwfc=1,param%nvec_to_cvg
       !    call daxpy(molecule(nmol)%mesh%nactive,&
       !         -coeff(idxwfc),molecule(nmol)%wf%wfc(:,idxwfc),1,&
       !         junk_wfc,1)  ! tdse_wfc+coeff*wfc(idxwfc) -> tdse_wfc
       ! end do
       
       ! open(unit=10,file="tdse2.dat",form='formatted')
       ! do i=1,molecule(nmol)%mesh%nactive
       !    write(10,*) i*molecule(nmol)%mesh%dx,junk_wfc(i)
       ! end do
       ! close(10)
       !
       ! - 4 - propagation of the wave packet
       !

       k0=2*PI/1.0
       do i=1,molecule(nmol)%mesh%nactive
          tdse_wfc(i)=cmplx(junk_wfc(i)*cos(k0*i*molecule(nmol)%mesh%dx),&
               junk_wfc(i)*sin(k0*i*molecule(nmol)%mesh%dx))
       end do
!       open(unit=10,file="tdse3.dat",form='formatted')
!       do i=1,molecule(nmol)%mesh%nactive
!          write(10,*) i*molecule(nmol)%mesh%dx,dreal(tdse_wfc(i)),dimag(tdse_wfc(i))
!       end do
!       close(10)
       print *,'Starting TDSE scheme'
       call tdse(molecule(nmol),param,tdse_wfc)
       
       
       
!       deallocate(coeff)
       deallocate(tdse_wfc)
    case("end")
       end_loop=.TRUE.
       
       !       call exit()
       ! ---------------------------------------------------------------
       !
       !                    COMMAND PART
       !
       ! --------------------------------------------------------------
    case("cmd")
       print *,'<----- ',field(2)
       select case (trim(field(2)))
          !
          ! to create a new molecule
          !
          ! ----------------------------------------------------------------------
          !
          ! to compute the wavefunction with davidson
          !
          ! ----------------------------------------------------------------------
          ! ----------------------------------------------------------------------
          !
          ! to operation
          !
          ! ----------------------------------------------------------------------
       case ("operation")
          print *,"nfield=",nfield
          allocate(junk_wfc(molecule(nmol)%mesh%nactive))
          junk_wfc=0.0
          allocate(coeff(1))
          do ifield=3,nfield-1,2
             coeff(1)=str2real(field(ifield))
             idxwfc=str2int(field(ifield+1))
             print *,"(nfield,coeff,idxwfc)=",nfield,coeff(1),idxwfc
             call daxpy(molecule(nmol)%mesh%nactive,&
                  coeff(1),molecule(nmol)%wf%wfc(:,idxwfc),1,&
                  junk_wfc,1)  ! junk_wfc+coeff*wfc(idxwfc) -> junk_wfc
          end do
          deallocate(coeff)
          open(unit=10,file="operation.dat",form='formatted')
          do i=1,molecule(nmol)%mesh%nactive
             write(10,*) i*molecule(nmol)%mesh%dx,junk_wfc(i),&
                  ( molecule(nmol)%wf%wfc(i,str2int(field(ifield+1))), ifield=3,nfield-1,2)
          end do
          close(10)
          deallocate(junk_wfc)
          ! ----------------------------------------------------------------------
          !
          !              TDSE
          !
          ! ----------------------------------------------------------------------
          ! ----------------------------------------------------------------------
          !
          ! to get info
          !
          ! ----------------------------------------------------------------------
       case ("info")
          call print_param(param)
          print *,nmol," molecule(s)"
       end select  ! field(2)
       !
       !
       !
    end select
  end subroutine parse_line
  ! --------------------------------------------------------------------------------------
  !
  !              wave_packet()
  !
  ! --------------------------------------------------------------------------------------
  subroutine wave_packet(r0,sig,Intens,molecule,wp)
    implicit none
    double precision::r0,sig,Intens
    type(t_molecule):: molecule
    double precision::wp(:)
    integer :: i
    do i=1,molecule%mesh%nactive
       wp(i)=gauss(i*molecule%mesh%dx,r0,sig,Intens)
    end do
  end subroutine wave_packet
  ! --------------------------------------------------------------------------------------
  !
  !              read_param()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_param(syst)
    implicit none
    integer::nmol
    type(t_molecule),allocatable:: molecule(:)
    type(t_system) :: syst


    type(t_param)::param
    character (len=1024)::line,redline
    character (len=1024)::line2
    character (len=32)::field(32)
    integer::nfield,debidx,endidx
    logical::exist_file
    integer::lline,eqidx,i,di
    logical::eol,end_loop

    syst%nmol=0
    call init_param(param)

    print *,'Reading ',trim(syst%inputfile)
    open(unit=2,file=syst%inputfile,form='formatted')
    end_loop=.FALSE.
    do while((.not.(is_iostat_end(param%ieof))).and.(.not.(end_loop)))
       read(2,'(A)') line
       call line_parser(line,nfield,field)
       print *,nfield,' --> ',(trim(field(i)),i=1,nfield)
       call parse_line(param,field,nfield,end_loop,syst%nmol,syst%molecule,syst%time_spent,syst)
       print *,"end_loop=",end_loop
    end do
    close(2)    

    call exit()
    read(*,*)  
  end subroutine read_param
  ! --------------------------------------------------------------------------------------
  !
  !              init_param()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_param(param)
    implicit none
    type(t_param)::param
    double precision,parameter::pi=4.0*atan(1.0)
    param%ieof=0
    param%loopmax=1000
    param%prefix='./'
    param%scheme='numerov'
    param%restart=.FALSE.
    param%init_wf=.TRUE.
    param%extrapol=.FALSE.
    param%extrap_add=10
    param%nvecmin=20
    param%nvecmax=41
    param%Nx=30
    param%noccstate=1
    param%nvec_to_cvg=1
    allocate(param%list_idx_to_cvg(param%nvec_to_cvg))
    allocate(param%occupation(param%nvec_to_cvg))
    param%list_idx_to_cvg(1)=1
    param%ETA=1.0e-3
    param%dim=1
    param%box%width=pi/sqrt(2.0)
    param%box%radius=.5*param%box%width
    param%box%shape='cube'
    param%box%center(1)=0.5
    param%box%center(2)=0.5
    param%box%center(3)=0.5
    param%perturb%Intensity=1.0
    param%perturb%sigma=1.0
    param%perturb%Intensity=1.0
    param%perturb%location(1)=0.0
    param%perturb%location(2)=0.0
    param%perturb%location(3)=0.0
    param%perturb%shape='gaussian'
    param%hartree=.FALSE.
    param%exchange=.FALSE.
    param%Z=1.0
    param%lorb=0
    param%tdse%nstep=1000
    param%tdse%freq_save=1000
  end subroutine init_param


  subroutine init_param_molecule(molecule,param)
    implicit none
    type(t_param)::param
    type(t_molecule):: molecule
    double precision,parameter::pi=4.0*atan(1.0)
    molecule%param%ieof=    param%ieof
    molecule%param%loopmax=param%loopmax
    molecule%param%prefix=param%prefix
    molecule%param%scheme=param%scheme
    molecule%param%restart=param%restart
    molecule%param%init_wf=param%init_wf
    molecule%param%extrapol=    param%extrapol
    molecule%param%extrap_add=    param%extrap_add
    molecule%param%nvecmin=    param%nvecmin
    molecule%param%nvecmax=    param%nvecmax
    molecule%param%Nx=    param%Nx
    molecule%param%noccstate=    param%noccstate
    molecule%param%nvec_to_cvg=    param%nvec_to_cvg
    if(allocated(molecule%param%list_idx_to_cvg)) deallocate(molecule%param%list_idx_to_cvg)
    allocate(molecule%param%list_idx_to_cvg(molecule%param%nvec_to_cvg))
    if(    allocated(molecule%param%occupation)) deallocate(molecule%param%occupation)
    allocate(molecule%param%occupation(molecule%param%nvec_to_cvg))
    molecule%param%list_idx_to_cvg(1)=param%list_idx_to_cvg(1)
    molecule%param%ETA=    param%ETA
    molecule%param%dim=param%dim
    molecule%param%box%width=    param%box%width
    molecule%param%box%radius=molecule%param%box%width
    molecule%param%box%shape=    param%box%shape
    molecule%param%box%center(1)=    param%box%center(1)
    molecule%param%box%center(2)=    param%box%center(2)
    molecule%param%box%center(3)=    param%box%center(3)
    molecule%param%perturb%Intensity=    param%perturb%Intensity
    molecule%param%perturb%sigma =   param%perturb%sigma
    molecule%param%perturb%Intensity=    param%perturb%Intensity
    molecule%param%perturb%location(1)=    param%perturb%location(1)
    molecule%param%perturb%location(2)=    param%perturb%location(2)
    molecule%param%perturb%location(3)=    param%perturb%location(3)
    molecule%param%perturb%shape=    param%perturb%shape
    molecule%param%hartree=    param%hartree
    molecule%param%exchange=    param%exchange
    molecule%param%Z=    param%Z
    molecule%param%lorb=    param%lorb
    molecule%param%tdse%nstep=param%tdse%nstep
    molecule%param%tdse%freq_save=param%tdse%freq_save
  end subroutine init_param_molecule
  ! --------------------------------------------------------------------------------------
  !
  !              print_param()
  !
  ! --------------------------------------------------------------------------------------

  subroutine print_param(param)
    implicit none
    integer::i
    type(t_param)::param
    print *,'----------------------------------------------------------------------------------'
    print *,'                     Parametrization (begin)'  
    print *,'----------------------------------------------------------------------------------'
    print *,'#filenrj=',trim(param%filenrj)
    print *,'#restart=',param%restart
    print *,'#scheme=',trim(param%scheme)
    print *,'#init_wf=',param%init_wf
    print *,'#extrapol=',param%extrapol
    print *,'#extrap_add=',param%extrap_add
    print *,'#loopmax=',param%loopmax
    print *,'#nvecmin=',param%nvecmin
    print *,'#nvecmax=',param%nvecmax
    print *,'#ETA=',param%ETA
    print *,'#nvec_to_cvg=',param%nvec_to_cvg
    print *,"#occupation= ",(param%occupation(i),i=1,param%nvec_to_cvg)
    print *,'#Zato=',param%Z
    print *,'#lorb=',param%lorb
    print *,'#hartree=',param%hartree
    print *,'#exchange=',param%exchange
    print *,'#box_width=',param%box%width
    print *,'#box_radius=',param%box%radius
    print *,'#box_shape=',param%box%shape
    print *,'#box_center=[',param%box%center(1),',',param%box%center(2),',',param%box%center(3),']'
    print *,'#Nx=',param%nx
    print *,'#noccstate=',param%noccstate
    print *,'#dh=',param%box%width/(param%Nx+1)
    print *,'#Dimension of the mesh=',param%dim
    print *,'#Perturbation shape=',param%perturb%shape
    print *,'#Magnitude of the perturbation=',param%perturb%Intensity
    print *,'#Spread of the perturbation=',param%perturb%sigma
    print *,'#Perturbation location=[',param%perturb%location(1),',',param%perturb%location(2),',',param%perturb%location(3),']'
    print *,'----------------------------------------------------------------------------------'
    print *,'                     Parametrization (end)'
    print *,'----------------------------------------------------------------------------------'

  end subroutine print_param
  ! --------------------------------------------------------------------------------------
  !
  !              line_parser()
  !
  ! --------------------------------------------------------------------------------------

  subroutine line_parser(line,nfield,field)
    implicit none
    character (len=1024)::line
    character (len=32)::field(32)
    integer::nfield,debidx,endidx,i
    
    nfield=0
    do while (.not.(scan(line,' ').eq.1))
       debidx=1
       endidx=index(line,' ')
       if(.not.(endidx.eq.1)) nfield=nfield+1
       field(nfield)=line(debidx:endidx-1)
       !write(field(nfield),'(A)') line(debidx:endidx-1)
       line=adjustl(line(endidx:))
    end do
  end subroutine line_parser
  ! --------------------------------------------------------------------------------------
  !
  !             new_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine new_molecule(molecule,param)
    implicit none
    type(t_molecule)::molecule
    type(t_param)::param

    logical::exist_file
    integer::i

    print *," Starting new_molecule()"

    if(param%dim.lt.3) param%box%center(3)=0.0
    if(param%dim.lt.2) param%box%center(2)=0.0

    print *,'#prefix=',param%prefix
    call system("mkdir -p "//param%prefix)
    write(param%filenameeigen,'(a,a)') param%prefix(:len_trim(param%prefix)),'/eigenvalues.dat'
    inquire (file=param%filenameeigen,exist=exist_file)
    if(exist_file) then
       call system("rm "//param%filenameeigen)
    end if
    open(unit=1,file=param%filenameeigen,form='formatted',status='unknown');   close(1)
    param%filenameeigen =param%filenameeigen(:len_trim(param%filenameeigen ))
    print *,'#filenameeigen=',trim(param%filenameeigen)

    write(param%filenrj,'(a,a)') trim(param%prefix),'/energy.dat'
    inquire (file=param%filenrj,exist=exist_file)
    if(exist_file) then
       call system("rm "//param%filenrj)
    end if
    open(unit=1,file=param%filenrj,form='formatted',status='unknown') ;close(1)


    call new_mesh(molecule%mesh)
    call init_pot(molecule)
    call save_potential(param,molecule)
    molecule%wf%nwfc=param%nvecmin   ! number of wfc min to cvg
    if(allocated(molecule%wf%eps))     deallocate(molecule%wf%eps)
    allocate(molecule%wf%eps(molecule%wf%nwfc))
    if(allocated(molecule%wf%epsprev))     deallocate(molecule%wf%epsprev)
    allocate(molecule%wf%epsprev(molecule%wf%nwfc))
    if(allocated(molecule%wf%deps))     deallocate(molecule%wf%deps)
    allocate(molecule%wf%deps(molecule%wf%nwfc))
    if(allocated(molecule%wf%wfc))     deallocate(molecule%wf%wfc)
    allocate(molecule%wf%wfc(molecule%mesh%nactive,molecule%wf%nwfc))
    if(allocated(molecule%rho))     deallocate(molecule%rho)
    allocate(molecule%rho(molecule%mesh%Ntot))

    molecule%cvg%nwfc=param%nvecmin
    if(allocated(molecule%cvg%wfc))     deallocate(molecule%cvg%wfc)
    allocate(molecule%cvg%wfc(molecule%cvg%nwfc))
    do i=1,molecule%cvg%nwfc
       molecule%cvg%wfc%cvg=.FALSE.
    end do
    molecule%cvg%nvec_to_cvg=param%nvec_to_cvg
    molecule%cvg%ETA=param%ETA
    if(allocated(molecule%cvg%list_idx))     deallocate(molecule%cvg%list_idx)
    allocate(molecule%cvg%list_idx(param%nvec_to_cvg))
    do i=1,param%nvec_to_cvg
       molecule%cvg%list_idx(i)=param%list_idx_to_cvg(i)
    end do



    print *,"End of new_molecule()"
  end subroutine new_molecule



end module param_mod
