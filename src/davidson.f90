module davidson_mod
  use global
  use time_tracking
  use mesh_mod
  use poten
  use param_mod
  use IO
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              INIT_BASIS_SET()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_basis_set(V,nvec,m,molecule,tag)
    implicit none
    integer :: nvec
    double precision,allocatable :: V(:,:)
    type(t_mesh)::m
    type(t_molecule)::molecule
    character(len=*)::tag

    double precision, external :: ddot
    double precision ::normloc
    integer :: i,j
    double precision,allocatable :: Vdump(:,:)
    type(t_GramSchmidt) :: GS
    integer,parameter :: seed = 86456
    
    allocate(Vdump(m%nactive,nvec))

    select case (tag)
    case ('rand')
       !    if(tag.eq.'rand') then
       print *,'Init_basis_set > Random initialization'
       call srand(seed)
       do i=1,nvec
          do j=1,m%nactive
             Vdump(j,i)=rand()
          end do
       end do
    case default
       !    else
       print *,'Init_basis_set >  wfc initialization'
       do i=1,nvec
          call dcopy(m%nactive,molecule%wf%wfc(:,i),1,Vdump(:,i),1) ! g->d
       end do
       !    end if
    end select

    do i=1,nvec
       normloc=ddot(m%nactive,Vdump(:,i),1,Vdump(:,i),1)
       normloc=1.0/sqrt(normloc)
       call dscal(m%nactive,normloc,Vdump(:,i),1)
    end do
    GS%nindep=1
    call GramSchmidt(V,Vdump,nvec,m,GS)

    deallocate(Vdump)
  end subroutine init_basis_set
  
  ! --------------------------------------------------------------------------------------
  !
  !             davidson()
  !
  ! --------------------------------------------------------------------------------------
  subroutine davidson(param,mesh,cvg,molecule,pot,time_spent)
    type(t_param)::param
    type(t_mesh)::mesh
    type(t_cvg)::cvg
    type(t_potential)::pot
    type(t_molecule)::molecule
    type(t_time)::time_spent
    
    integer :: nvec
    double precision,allocatable :: V(:,:) ! wavefunctions
!    double precision::normloc
 !   double precision, external :: ddot
    integer :: iloop !,lorb

    nvec=param%nvecmin
    allocate(V(mesh%nactive,nvec))
    call davidson_basis_init(param,mesh,nvec,V,molecule)
      
      
    iloop=1
    cvg%ncvg=0
    molecule%wf%epsprev(:)=0.0
    molecule%wf%deps(:)=0.0
    do while((iloop.le.param%loopmax).and.(cvg%ncvg.lt.cvg%nvec_to_cvg))
       write(*,'(A)') 'Main > #######################################'     
       write(*,'(A,I4,A)') 'Main > ############ scf loop=',iloop,' ############'
       write(*,'(A)') 'Main > #######################################'     
       call davidson_step(nvec,V,mesh,param%nvecmin,iloop,cvg,pot,time_spent,molecule,param)

       call compute_density(molecule,param,V,mesh)
!       call read_param(param)
!       cvg%ETA=param%ETA
       
    end do
    call save_config(V,mesh,param%nvecmin,param)
    
    
    !--------------------------------------------------------------------------
    !
    ! save the wavefunctions
    !
    !--------------------------------------------------------------------------
    call save_wavefunction(param,mesh,V,molecule)
    deallocate(V)
    print *,"# Davidson > end of davidson subroutine"
  end subroutine davidson
  ! --------------------------------------------------------------------------------------
  !
  !              DAVIDSON_BASIS_INIT()
  !
  ! --------------------------------------------------------------------------------------
    subroutine davidson_basis_init(param,mesh,nvec,V,molecule)
      implicit none
      type(t_param)::param
      type(t_mesh)::mesh
      integer::nvec
      double precision ,allocatable:: V(:,:) ! wavefunctions
      type(t_molecule)::molecule
          
      select case (param%init_wf)
      case (.TRUE.)
         print *,"Davidson > Initialisation of the wf"
         select case (param%restart)
         case (.TRUE.)
            print *,'Davidson> restart an old calculation'
            call read_config(V,mesh,nvec)
         case (.FALSE.)
            print *,"Davidson > new calculation"
            call init_basis_set(V,nvec,mesh,molecule,'rand')
         end select
      case (.FALSE.)
         call init_basis_set(V,nvec,mesh,molecule,'wfc')
         !      call check_ortho(V,nvec,mesh%nactive)
         !       stop
      end select
    end subroutine davidson_basis_init
  ! --------------------------------------------------------------------------------------
  !
  !              DAVIDSON_STEP()
  !
  ! --------------------------------------------------------------------------------------
  subroutine davidson_step(nvec,V,m,nvecmin,iloop,cvg,pot,time_spent,molecule,param)
    implicit none
    type(t_molecule)::molecule
    type(t_param)::param
    integer :: nvec,nvecmin,iloop
    double precision,allocatable :: V(:,:)   !,pot_ext(:)
    type(t_mesh) :: m
    type(t_cvg)::cvg
    type(t_time)::time_spent
    type(t_potential)::pot
    
    double precision,allocatable :: S(:)          ! eigenvalues
    double precision,allocatable :: T(:,:)        ! reduced matrix T
    double precision,allocatable :: VRitz(:,:)    ! Ritz's vectors
    double precision,allocatable :: residual(:,:) ! residual
    double precision,allocatable :: delta(:,:)    ! delta vectors
    double precision,allocatable :: Vnew(:,:)     ! Vnew
    integer :: i
    integer :: ndelta
    integer :: newnvec
    type(t_GramSchmidt)::GS
    ! T (reduced matrix) computing
    allocate(T(nvec,nvec))
    call cpu_time(time_spent%start_loc)
    call compute_T2(T,V,nvec,m,pot%tot)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> compute_T')
    
    ! Diagonatilzation of T
    allocate(S(nvec))

    call cpu_time(time_spent%start_loc)
    call diagonalization(S,T,nvec)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Diagonalization')


    do i=1,cvg%nwfc
       cvg%wfc(i)%dnrj=S(i)-cvg%wfc(i)%nrj
       cvg%wfc(i)%nrjold=cvg%wfc(i)%nrj
       cvg%wfc(i)%nrj=S(i)
    end do
    molecule%wf%eps(1:nvecmin)=S(1:nvecmin)
    molecule%wf%deps(:)=molecule%wf%eps(1:nvecmin)-molecule%wf%epsprev(:)
    molecule%wf%epsprev(:)=molecule%wf%eps(1:nvecmin)
    ! computation of the Ritz's vectors
    allocate(VRitz(m%nactive,nvec))

    call cpu_time(time_spent%start_loc)
    call Ritz(VRitz,V,T,nvec)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Diagonalization')

    
    ! computation of residual
    allocate(residual(m%nactive,nvec))
    

    call cpu_time(time_spent%start_loc)
    call compute_residual2(residual,VRitz,molecule%wf%eps,nvec,cvg,m,pot%tot)

    cvg%ncvg=0
    do i=1,cvg%nwfc
       if(cvg%wfc(i)%cvg) cvg%ncvg=cvg%ncvg+1
    end do
    do i=1,cvg%nwfc
       write(*,'(A,I3,A,F12.6,A,E12.2,A)',advance='no') 'Main > Eigenvalue(',i,'): ',cvg%wfc(i)%nrj,'(',cvg%wfc(i)%dnrj,')'
       write(*,'(A,I3,A,E12.4,A,E10.2)',advance='no') ' | r(',i,')= ',cvg%wfc(i)%resi,'/',cvg%ETA
       write(*,'(L2,I2)') cvg%wfc(i)%cvg,cvg%ncvg
!       if(cvg%wfc(i)%cvg) then
!          write(*,*) ' -> cvg'
!       else
!          write(*,*)
!       end if
    end do
    !call cpu_time(inter)
    open(unit=1,file=param%filenameeigen,form='formatted',status='unknown',access='append')
    write(1,*) iloop,molecule%wf%eps(1:nvecmin)
    close(1)


    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Residual')

    ! computation of delta
    allocate(delta(m%nactive,nvec))
    delta(:,:)=0.0
    !call cpu_time(inter)

    call cpu_time(time_spent%start_loc)
    call compute_delta(delta,residual,molecule%wf%eps,nvec,cvg,m,ndelta,pot%tot)
    call cpu_time(time_spent%end_loc)
    call time_tracking_write(iloop,time_spent,'Davidson -> Delta')

    
    deallocate(V)
    allocate(V(m%nactive,nvec+ndelta))
    V(:,1:nvec)=VRitz(:,:)
    print *,'ndelta=',ndelta
    print *,'nvec=',nvec
    V(:,nvec+1:nvec+ndelta)=delta(:,:ndelta)
    allocate(Vnew(m%nactive,nvec+ndelta))
    newnvec=nvec+ndelta
    GS%nindep=nvec
    !call cpu_time(inter)
    call GramSchmidt(Vnew,V,newnvec,m,GS)    
    ! call cpu_time(inter2);     call dbg(iloop,inter,inter2,'GS')
    
    print *,'Main > ',GS%ndep,newnvec
    
    deallocate(V)
    if(newnvec.le.param%nvecmax) then
       nvec=newnvec
    else
       print *,'Main > restart from nvecmin'
       nvec=nvecmin
    end if
    allocate(V(m%nactive,nvec))
    V(:,:)=Vnew(:,1:nvec)
    
    !call check_ortho(V,nvec,m%nactive)
    print *,'Main > New size of the basis ',nvec
    iloop=iloop+1
    deallocate(S)
    deallocate(T)
    deallocate(VRitz)
    deallocate(residual)
    deallocate(delta)
    deallocate(Vnew)
  end subroutine davidson_step
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_T()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_T(T,V,nvec,m,pot_tot)
    implicit none
    double precision,allocatable :: V(:,:),T(:,:),pot_tot(:)
    integer :: nvec
    type(t_mesh)::m
    
    integer :: i,j,k,l
    double precision :: deltasqr,acc
    double precision, parameter::alpha=0.0
!    double precision::beta
    deltasqr=m%dx**2
    !$OMP PARALLEL private(acc) 
    !$OMP DO 
    do j=1,nvec
       do i=1,nvec ! Tij
          T(i,j)=0.0
          do k=1,m%nactive
             acc=(m%dim/deltasqr+pot_tot(k))*V(k,j) ! the potential will be added here
             do l=1,m%node(k)%n_neighbors
                acc=acc-0.5*V(m%node(k)%list_neighbors(l),j)/deltasqr
!                print *, omp_get_thread_num(),i,j,k,l
             end do
             T(i,j)=T(i,j)+V(k,i)*acc
!             print *,' -> ',omp_get_thread_num(),i,j,T(i,j)
          end do
!          print *,omp_get_thread_num(),i,j,T(i,j)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
!    stop

  end subroutine compute_T
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_T()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_T2(T,V,nvec,m,pot_tot)
    implicit none
    double precision,allocatable :: V(:,:),T(:,:),pot_tot(:)
    integer :: nvec
    type(t_mesh)::m
    
    integer :: i,j,l,nn
    double precision :: deltasqr,acc
    double precision, parameter::alpha=0.0
!    double precision::beta
    deltasqr=m%dx**2
    !$OMP PARALLEL private(acc) 
    !$OMP DO 
    do j=1,nvec
       do i=1,nvec ! Tij
          T(i,j)=0.0
          do nn=1,m%nactive
             acc=(m%dim/deltasqr+pot_tot(nn))*V(nn,j) ! the potential will be added here
             do l=1,m%node(nn)%n_neighbors
                acc=acc-0.5*V(m%node(nn)%list_neighbors(l),j)/deltasqr
!                print *, omp_get_thread_num(),i,j,k,l
             end do
             T(i,j)=T(i,j)+V(nn,i)*acc
!             print *,' -> ',omp_get_thread_num(),i,j,T(i,j)
          end do
!          print *,omp_get_thread_num(),i,j,T(i,j)
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
!    stop

  end subroutine compute_T2
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_DELTA()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_delta(delta,r,lambda,nvec,cvg,m,ndelta,pot_tot)
    ! INPUT: the residual |r>, the Ritz's vectors |VRitz>, the eigenvalues lambda
    ! OUTPUT : the correction |delta> to improve the  Ritz's vectors so that to
    !          minimize the residual
    implicit none
    type(t_mesh)::m
    double precision :: lambda(:),r(:,:),delta(:,:),pot_tot(:)
    integer :: nvec,ndelta
    type(t_cvg)::cvg
    
    double precision, external :: ddot
    double precision, parameter::alpha=1.0,beta=0.0
    double precision, allocatable :: normloc
!    double precision, allocatable :: Dinv(:,:)
    integer :: i,j
    double precision :: deltasqr

    deltasqr=m%dx**2
    print *,'Delta > ---------------------'
    print *,'Delta > --- compute_delta ---'
    print *,'Delta > ---------------------'
    delta(:,:)=0.0
   ! allocate(Dinv(m%nactive,m%nactive))
  !  Dinv(:,:)=0.0
    ndelta=0
    do i=1,nvec
       ! delta is computed only if the eigenvectors isn't converged
       if(.not.(cvg%wfc(i)%cvg)) then
          ndelta=ndelta+1
          do j=1,m%nactive
             delta(j,ndelta)=r(j,ndelta)/((m%dim/deltasqr+pot_tot(j))-lambda(i))
!             Dinv(j,j)=1.0/((3.0/deltasqr+pot_tot(j))-lambda(i))
          end do
          ! see Victor Eijkhout in "Introduction to scientific and technical computing" edited by Willmore et al
          ! Chap 15 Libraries for Linear Algebra
          ! to get a comprehensive way to use dgemv
 !         call dgemv('N',m%nactive,m%nactive,alpha,Dinv,m%nactive,r(:,ndelta),1,beta,delta(:,ndelta),1)
          !norm=sqrt(ddot(m%nactive,delta(:,i),1,delta(:,i),1))
          !write(*, '(A10,I4,A2,E12.6)',advance='no') ' Delta > delta(',i,')=',norm
          !delta(:,ndelta)=delta(:,ndelta)+VRitz(:,ndelta)
          normloc=1.0/sqrt(ddot(m%nactive,delta(:,ndelta),1,delta(:,ndelta),1))
          call dscal(m%nactive,normloc,delta(:,ndelta),1)
       end if
    end do
    !deallocate(Dinv)
    print *,'Delta > ',ndelta,' new vector(s)'
  end subroutine compute_delta
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_RESIDUAL()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_residual(r,VRitz,S,nvec,cvg,m,pot_tot)
    implicit none
    type(t_mesh)::m
    integer :: nvec
    double precision :: r(:,:),VRitz(:,:),S(:),pot_tot(:)
    type(t_cvg) :: cvg

    integer :: i,j,k
    double precision, external :: ddot        
    double precision :: deltasqr

    
    print *,'Residual > ------------------------'
    print *,'Residual > --- compute residual ---'
    print *,'Residual > ------------------------'
    deltasqr=m%dx**2
    r(:,:)=0.0
    cvg%ncvg=0
    do j=1,nvec
       do i=1,m%nactive
          r(i,j)=(m%dim/deltasqr+pot_tot(i))*VRitz(i,j)
          do k=1,m%node(i)%n_neighbors
             r(i,j)=r(i,j)-0.5*VRitz(m%node(i)%list_neighbors(k),j)/deltasqr
          end do
          r(i,j)=r(i,j)-S(j)*VRitz(i,j)
       end do
       if(j.le.cvg%nwfc) then
          cvg%wfc(j)%resi=m%dv*ddot(m%nactive,r(:,j),1,r(:,j),1)
          if (cvg%wfc(j)%resi.le.cvg%ETA) then
             cvg%wfc(j)%cvg=.TRUE.
             if(j.le.5) cvg%ncvg=cvg%ncvg+1
          end if
       end if
    end do
  end subroutine compute_residual
  ! --------------------------------------------------------------------------------------
  !
  !              COMPUTE_RESIDUAL()
  !
  ! --------------------------------------------------------------------------------------
  subroutine compute_residual2(r,VRitz,S,nvec,cvg,m,pot_tot)
    implicit none
    type(t_mesh)::m
    integer :: nvec
    double precision :: r(:,:),VRitz(:,:),S(:),pot_tot(:)
    type(t_cvg) :: cvg

    integer :: nn,j,k
    double precision, external :: ddot        
    double precision :: deltasqr

    
    print *,'Residual > ------------------------'
    print *,'Residual > --- compute residual ---'
    print *,'Residual > ------------------------'
    deltasqr=m%dx**2
    r(:,:)=0.0
    cvg%ncvg=0
    do j=1,nvec
       do nn=1,m%nactive
          r(nn,j)=(m%dim/deltasqr+pot_tot(nn))*VRitz(nn,j)
          do k=1,m%node(nn)%n_neighbors
             r(nn,j)=r(nn,j)-0.5*VRitz(m%node(nn)%list_neighbors(k),j)/deltasqr
          end do
          r(nn,j)=r(nn,j)-S(j)*VRitz(nn,j)
       end do
       if(j.le.cvg%nwfc) then
          cvg%wfc(j)%resi=m%dv*ddot(m%nactive,r(:,j),1,r(:,j),1)
          if (cvg%wfc(j)%resi.le.cvg%ETA) then
             cvg%wfc(j)%cvg=.TRUE.
             if(j.le.5) cvg%ncvg=cvg%ncvg+1
          end if
       end if
    end do
  end subroutine compute_residual2
  ! --------------------------------------------------------------------------------------
  !
  !              RITZ()
  !
  ! --------------------------------------------------------------------------------------
  subroutine Ritz(Vout,Vin,y,nvec)
    implicit none
    double precision :: Vin(:,:),Vout(:,:),y(:,:)
    integer :: nvec

    integer :: i,j
    print *,'Ritz > --------------'
    print *,'Ritz > --- Ritz() ---'
    print *,'Ritz > --------------'
    do i=1,nvec
       Vout(:,i)=y(1,i)*Vin(:,1)
       do j=2,nvec
          Vout(:,i)=Vout(:,i)+y(j,i)*Vin(:,j)
       end do
    end do
  end subroutine Ritz

end module davidson_mod
