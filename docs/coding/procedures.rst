Procedures
--------------

Procedures should follow this template:::

  !===============================================================================
  !> Procedure description ...
  !===============================================================================  
  SUBROUTINE updatet (self, aux)
    class(some_t):: self
    integer:: aux
    !..............................................................................
    integer:: local_variables, ...
    real(KindScalarVar), dimension(:,:,:):: d, ...
    !------------------------------------------------------------------------------
    call trace%begin('some_t%update)
    ...
    !------------------------------------------------------------------------------
    ! Comment block
    !------------------------------------------------------------------------------
    ...
    var = expression                                    ! in-line comment
    ...
    call trace%end ()
  END SUBROUTINE update

or, if the procedures should be timed:::

  !===============================================================================
  !> Procedure description ...
  !===============================================================================  
  SUBROUTINE updatet (self, aux)
    class(some_t):: self
    integer:: aux
    !..............................................................................
    integer:: local_variables, ...
    real(KindScalarVar), dimension(:,:,:):: d, ...
    integer, save:: itimer
    !------------------------------------------------------------------------------
    call trace%begin('some_t%update, itimer=itimer)
    ...
    call trace%end (itimer)
  END SUBROUTINE update


.. toctree::
   :maxdepth: 3
