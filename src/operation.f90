module operation_mod
  use global
  use mesh_mod
  use param_mod
  use IO
  implicit none
contains
  ! --------------------------------------------------------------------------------------
  !
  !              Operation()
  !
  ! --------------------------------------------------------------------------------------
  subroutine operation(param,mesh)
    implicit none
    type(t_param)::param
    type(t_mesh)::mesh
    integer::nvec
    double precision ,allocatable:: V(:,:) ! wavefunctions
    nvec=param%nvecmin
    allocate(V(mesh%nactive,nvec))
    call read_config(V,mesh,nvec)
    deallocate(V)
  end subroutine operation
end module operation_mod
