module cvg
  use global
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !             init_cvg()
  !
  ! --------------------------------------------------------------------------------------
  subroutine init_cvg(molecule)
    implicit none
    type(t_molecule)::molecule
    integer::i
    molecule%cvg%nwfc=molecule%cvg%nvec_to_cvg
    if(allocated(molecule%cvg%wfc))     deallocate(molecule%cvg%wfc)
    allocate(molecule%cvg%wfc(molecule%cvg%nwfc))
    do i=1,molecule%cvg%nwfc
       molecule%cvg%wfc%cvg=.FALSE.
    end do
    if(allocated(molecule%cvg%list_idx))     deallocate(molecule%cvg%list_idx)
    allocate(molecule%cvg%list_idx(molecule%cvg%nvec_to_cvg))
    do i=1,molecule%cvg%nvec_to_cvg
       molecule%cvg%list_idx(i)=i
    end do

  end subroutine init_cvg
end module cvg
