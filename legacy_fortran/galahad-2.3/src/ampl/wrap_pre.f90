! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
!! Implements an F90 wrapper around the solvers implemented in the
!!
!!            G  A  L  A  H  A  D
!!
!! optimization library, to make functions within Fortran 90 modules
!! visible to the linker when linking with a C object file.
!!
!! Grants access to
!!
!!   PRE
!!     Quadratic Program Preprocessor,
!!
!! D. Orban@ECE                                   Chicago, March 2003
!!
!! Thanks to the f2py project      http://cens.ioc.ee/projects/f2py2e
!!===================================================================

!! ===============================
!! Main wrapper around USE_PRE( )
!! ===============================

Subroutine Wrap_Use_Pre( setup_use_pre )

  Use GALAHAD_USEPRE_double      ! Main PRE driver module
  External setup_use_pre         ! Declared in galahad.c

  ! Pass pointer to subroutine USE_PRE to calling function
  Call setup_use_pre( USE_PRE )

  Return

End Subroutine Wrap_Use_Pre

!!=================================
!! End of wrapper
!! ================================
