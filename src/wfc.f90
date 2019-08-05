module wfc
  use global
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !             init_wfc()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_wfc(molecule)
    implicit none
    type(t_molecule)::molecule
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

  end subroutine init_wfc
end module wfc
