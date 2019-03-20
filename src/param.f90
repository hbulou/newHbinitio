module param_mod
  use global
  implicit none
contains
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
    param%box%radius=param%box%width
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
  end subroutine init_param
  ! --------------------------------------------------------------------------------------
  !
  !              parse_line()
  !
  ! --------------------------------------------------------------------------------------
  subroutine parse_line(param,field,end_loop)
    implicit none
    type(t_param)::param
    character (len=32)::field(32)
    logical::end_loop
    
    integer::i

    select case (field(1))
    case("box_radius >") 
       read(field(2),*) param%box%radius
    case("box_width >") 
       read(field(2),*) param%box%width
    case("box_shape >") 
       read(field(2),*) param%box%shape
    case("box_center >") 
       read(field(2),*) param%box%center(1)
       read(field(3),*) param%box%center(2)
       read(field(4),*) param%box%center(3)
    case("dimension >") 
       read(field(2),*) param%dim
    case("ETA >") 
       read(field(2),*) param%ETA
    case("exchange >") 
       read(field(2),*) param%exchange
    case("extrap_add >") 
       read(field(2),*) param%extrap_add
    case("extrapol >") 
       read(field(2),*) param%extrapol
    case("hartree >") 
       read(field(2),*) param%hartree
    case("init_wf >") 
       read(field(2),*) param%init_wf
    case("perturb_intensity >") 
       read(field(2),*) param%perturb%Intensity
    case("perturb_sigma >") 
       read(field(2),*) param%perturb%sigma
    case("perturb_shape >") 
       read(field(2),*) param%perturb%shape
    case("perturb_location >") 
       read(field(2),*) param%perturb%location(1)
       read(field(3),*) param%perturb%location(2)
       read(field(4),*) param%perturb%location(3)
    case("loopmax >") 
       read(field(2),*) param%loopmax
    case("noccstate >") 
       read(field(2),*) param%noccstate
    case("nvecmax >") 
       read(field(2),*) param%nvecmax
    case("nvecmin >") 
       read(field(2),*) param%nvecmin
    case("nvec_to_cvg >") 
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
    case("Nx >") 
       read(field(2),*) param%Nx
    case("occupation >") 
       do i=1,param%nvec_to_cvg
          read(field(1+i),*) param%occupation(i)
       end do
    case("prefix >") 
       read(field(2),*) param%prefix
    case("restart >") 
       read(field(2),*) param%restart
    case("scheme >") 
       read(field(2),*) param%scheme
    case("Zato >") 
       read(field(2),*) param%Z
    case("lorb >") 
       read(field(2),*) param%lorb
    case("cmd >") 
       if(field(2).eq."end") then
          end_loop=.TRUE.
       end if
    end select


  end subroutine parse_line
  ! --------------------------------------------------------------------------------------
  !
  !              read_param()
  !
  ! --------------------------------------------------------------------------------------
  subroutine read_param(param)
    implicit none
    type(t_param)::param

    character (len=1024)::line,redline
    character (len=1024)::line2
    character (len=32)::field(32)
    integer::nfield,debidx,endidx
    logical::exist_file
    integer::lline,eqidx,i,di
    logical::eol,end_loop

    call init_param(param)


    print *,'Reading ',param%inputfile
    open(unit=1,file=param%inputfile,form='formatted')
    end_loop=.FALSE.
    do while((.not.(is_iostat_end(param%ieof))).and.(.not.(end_loop)))
       read(1,'(A)') line
       call line_parser(line,nfield,field)
       print *,nfield,' --> ',(trim(field(i)),i=1,nfield)
       call parse_line(param,field,end_loop)
    end do
    close(1)    


    if(param%dim.lt.3) param%box%center(3)=0.0
    if(param%dim.lt.2) param%box%center(2)=0.0
    

    print *,'#prefix=',param%prefix
    call system("mkdir "//param%prefix)
    write(param%filenameeigen,'(a,a)') param%prefix(:len_trim(param%prefix)),'/eigenvalues.dat'
    inquire (file=param%filenameeigen,exist=exist_file)
    if(exist_file) then
       call system("rm "//param%filenameeigen)
    end if
    open(unit=1,file=param%filenameeigen,form='formatted',status='unknown');   close(1)
!    write(filename,'(a,a)') param%prefix(:len_trim(param%prefix)),'/evectors.dat'
 !   print *,filename
  !  open(unit=1,file=filename,form='formatted',status='unknown')
   ! close(1)
    !stop
    param%filenameeigen =param%filenameeigen(:len_trim(param%filenameeigen ))
    print *,'#filenameeigen=',trim(param%filenameeigen)

    write(param%filenrj,'(a,a)') trim(param%prefix),'/energy.dat'
    inquire (file=param%filenrj,exist=exist_file)
    if(exist_file) then
       call system("rm "//param%filenrj)
    end if
    open(unit=1,file=param%filenrj,form='formatted',status='unknown') ;close(1)
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
    
  end subroutine read_param


  ! --------------------------------------------------------------------------------------
  !
  !              line_parser()
  !
  ! --------------------------------------------------------------------------------------

      subroutine line_parser(line,nfield,field)
        implicit none
        character (len=1024)::line,redline
        character (len=32)::field(32)
        integer::eqidx,lline,nfield,debidx,endidx,i
        eqidx=index(line,"=")
        nfield=1
        field(nfield)=line(1:eqidx-1)//' >'
        nfield=nfield+1
        redline=line(eqidx+1:)
        lline=len_trim(redline)
        if(lline.gt.0) then
           !nfield=1
           !print *,lline,trim(redline)
           debidx=1
           do i=1,len(trim(redline))
              !print *,redline(i:i)
              if(redline(i:i).eq.' ') then
                 endidx=i-1
                 field(nfield)=trim(redline(debidx:endidx))
                 !print *,'                 >>>>',debidx,endidx,redline(debidx:endidx)
                 debidx=endidx+1
                 nfield=nfield+1
              end if
              field(nfield)=trim(redline(debidx:))
           end do
        end if
      end subroutine line_parser
    

end module param_mod
