module molecule_mod
  use global
  use IO
  use poten
  use mesh_mod
  implicit none
  
contains
    ! --------------------------------------------------------------------------------------
  !
  !             new_molecule()
  !
  ! --------------------------------------------------------------------------------------
  subroutine new_molecule(molecule,param)
    implicit none
    type(t_molecule)::molecule
    type(t_param)::param
    integer::i
    
    print *,"> Entering new_molecule()"
    
    call new_mesh(molecule%mesh)
    call init_pot(molecule)
    call save_potential(param,molecule%mesh,molecule)

    molecule%wf%nwfc=param%nvecmin   ! number of wfc min to cvg
    allocate(molecule%wf%eps(molecule%wf%nwfc))
    allocate(molecule%wf%epsprev(molecule%wf%nwfc))
    allocate(molecule%wf%deps(molecule%wf%nwfc))
    allocate(molecule%wf%wfc(molecule%mesh%nactive,molecule%wf%nwfc))
    allocate(molecule%rho(molecule%mesh%Ntot))

    molecule%cvg%nwfc=param%nvecmin
    allocate(molecule%cvg%wfc(molecule%cvg%nwfc))
    do i=1,molecule%cvg%nwfc
       molecule%cvg%wfc%cvg=.FALSE.
    end do
    molecule%cvg%nvec_to_cvg=param%nvec_to_cvg
    molecule%cvg%ETA=param%ETA
    allocate(molecule%cvg%list_idx(param%nvec_to_cvg))
    do i=1,param%nvec_to_cvg
       molecule%cvg%list_idx(i)=param%list_idx_to_cvg(i)
    end do

    if(allocated(molecule%numerov%Q)) deallocate(molecule%numerov%Q)
    allocate(molecule%numerov%Q(molecule%mesh%nactive))
    if(allocated(molecule%numerov%r)) deallocate(molecule%numerov%r)
    allocate(molecule%numerov%r(molecule%mesh%nactive))
    print *,"#################### nactive= ",molecule%mesh%nactive

    
  end subroutine new_molecule

end module molecule_mod
