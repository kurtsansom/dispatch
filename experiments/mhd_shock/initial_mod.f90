!===============================================================================
!> A dummy `initial_mod` to override the default one.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: initial_t
    logical:: mhd=.true.
  contains
    procedure:: init
    procedure:: condition
  end type
CONTAINS

!===============================================================================
!> Stub -- not used
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
END SUBROUTINE init

!===============================================================================
!> Stub -- not used
!===============================================================================
SUBROUTINE condition (self, m, ff, idx)
  class(initial_t):: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: ff
  class(mesh_t), pointer, dimension(:):: m
END SUBROUTINE condition

END MODULE initial_mod
