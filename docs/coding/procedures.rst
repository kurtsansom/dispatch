Modules
--------------

Modules should follow this template:::

  !===============================================================================
  !> Module description ...
  !===============================================================================
  MODULE some_mod
    USE io_mod
    USE ...
    implicit none
    private
    type, public:: some_t
    contains
      procedure:: init
      procedure:: update
    end type
    type(some_t):: some
  CONTAINS
  
  !===============================================================================
  !> Procedure description ...
  !===============================================================================  
  SUBROUTINE init (self)
    class(some_t):: self
    ...
  END SUBROUTINE init 

  !===============================================================================
  !> Procedure description ...
  !===============================================================================  
  SUBROUTINE updatet (self)
    class(some_t):: self
    !..............................................................................
    integer, save:: itimer=0
    !------------------------------------------------------------------------------
    call trace%begin('some_t%update, itimer=itimer)
    ...
    call trace%end (itimer)
  END SUBROUTINE update

  END MODULE some_mod


.. toctree::
   :maxdepth: 3
