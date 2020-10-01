! THIS VERSION: GALAHAD 2.2 - 21/04/2008 AT 16:00 GMT.

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*     SYMBOLS  M O D U L E    *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Nick Gould and Philippe Toint

!  History -
!   originally released pre GALAHAD Version 1.0. December 1st 2000
!   update released with GALAHAD Version 2.0. February 16th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_SYMBOLS

!  This module provides the list of all symbolic names that are common to
!  the GALAHAD modules. It is just intended as a dictionnary for use in other
!  modules.
!
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

      IMPLICIT NONE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!                       The public symbolic names

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!                            General
!-------------------------------------------------------------------------------

!     General integers

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_1                      =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_2                      =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_3                      =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_4                      =   4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_5                      =   5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_6                      =   6
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_7                      =   7
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_8                      =   8
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_9                      =   9
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_10                     =  10
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_11                     =  11
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_12                     =  12

!     Matrix storage schemes

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DIAGONAL               =  -3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DENSE                  =  -2 
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SPARSE_BY_ROWS         =  -1 
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_COORDINATE             =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ELEMENTAL              =   2

!     Constraint and variable status
!     We must have that ELIMINATED < ACTIVE < other active status < FREE

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_INACTIVE               =  -2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_STRUCTURAL             =  -1 
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ELIMINATED             =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ACTIVE                 =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FIXED                  =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_RANGE                  =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_UPPER                  =   4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_LOWER                  =   5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FREE                   =   6

!     Sign conventions

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_POSITIVE               =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NEGATIVE               =  -1

!     Special values

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ALL_ZEROS              =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ALL_ONES               =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_GENERAL                =   2

!     Print levels

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SILENT                 =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_TRACE                  =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ACTION                 =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DETAILS                =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DEBUG                  =   4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_CRAZY                  =   5

!     Checking levels

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_BASIC                  =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SEVERE                 =   2

!     CUTEr problem types

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_UNCONSTRAINED          =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_CONSTRAINED            =   1

!     Exit conditions
!     The idea is to reserve the codes from 0 to -20 for generic GALAHAD
!     return conditions, while the rest of the range is tool specific.

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SUCCESS                =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_MEMORY_FULL            =  -1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FILE_NOT_OPENED        =  -2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_COULD_NOT_WRITE        =  -3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_TOO_FEW_BITS_PER_BYTE  =  -4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_PROGRESS_IMPOSSIBLE    =  -5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_MAX_ITERATIONS_REACHED =  -6
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NOT_INITIALIZED        =  -7
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_WRONG_N                =  -8
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_WRONG_M                =  -9
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SORT_TOO_LONG          =  -10
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NOT_DIAGONAL           =  -11

!  New exit conditions (0 to -30; others will be package specific)

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ok                      = 0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_allocate          = - 1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_deallocate        = - 2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_restrictions      = - 3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_bad_bounds        = - 4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_primal_infeasible = - 5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_dual_infeasible   = - 6
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_unbounded         = - 7
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_no_center         = - 8
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_analysis          = - 9
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_factorization     = - 10
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_solve             = - 11
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_uls_analysis      = - 12
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_uls_factorization = - 13
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_uls_solve         = - 14
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_preconditioner    = - 15
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_ill_conditioned   = - 16
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_tiny_step         = - 17
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_max_iterations    = - 18
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_cpu_limit         = - 19
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_inertia           = - 20
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_file              = - 21
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_io                = - 22
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_upper_entry       = - 23
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_sort              = - 24
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_error_input_status      = - 25
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_warning_on_boundary     = - 30

!     Final status for files produced by a module

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_KEEP                   =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DELETE                 =   1

!     Miscellaneous
!     Note: AUTOMATIC must be different from NEWTON, GAUSS_NEWTON and all
!           the other values in this paragraph. NONE must not overlap any
!           of the matrix stotage schemes.

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_UNDEFINED              = -100
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NONE                   =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_USER_DEFINED           =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NEVER                  =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_INITIAL                =   4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_ALWAYS                 =   5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_AUTOMATIC              =  10


!-------------------------------------------------------------------------------
!                            PRESOLVE
!-------------------------------------------------------------------------------

!     PRESOLVE: termination strategies

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_REDUCED_SIZE           =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FULL_PRESOLVE          =   2

!     PRESOLVE: policy wrt dual variables and multipliers (inactive)

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FORCE_TO_ZERO          =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_LEAVE_AS_IS            =   1

!     PRESOLVE: final status of the bounds on the variables

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_TIGHTEST               =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_NON_DEGENERATE         =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_LOOSEST                =   2

!-------------------------------------------------------------------------------
!                            FILTRANE
!-------------------------------------------------------------------------------

!      Model type

       INTEGER, PUBLIC, PARAMETER :: GALAHAD_GAUSS_NEWTON          =  0
       INTEGER, PUBLIC, PARAMETER :: GALAHAD_NEWTON                =  1

!      Subproblem accuracy

       INTEGER, PUBLIC, PARAMETER :: GALAHAD_ADAPTIVE              =  0
       INTEGER, PUBLIC, PARAMETER :: GALAHAD_FULL                  =  1

!      Margin types (must be different from GALAHAD_FIXED)

       INTEGER, PUBLIC, PARAMETER :: GALAHAD_CURRENT               =  0
       INTEGER, PUBLIC, PARAMETER :: GALAHAD_SMALLEST              =  1

!      Automatic model criteria

       INTEGER, PUBLIC, PARAMETER :: GALAHAD_BEST_FIT              =  0
       INTEGER, PUBLIC, PARAMETER :: GALAHAD_BEST_REDUCTION        =  1
       
!-------------------------------------------------------------------------------
!                            LANCELOT
!-------------------------------------------------------------------------------

!     First derivative approximations

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_EXACT                  =   0
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_FORWARD                =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_CENTRAL                =   2

!     Second derivative approximations

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_BFGS                   =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DFP                    =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_PSB                    =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SR1                    =   4
      
!     Linear solver

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_CG                     =   1
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_DIAGONAL_CG            =   2
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_USERS_CG               =   3
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_EXPANDING_BAND_CG      =   4
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_MUNKSGAARD_CG          =   5
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_SCHNABEL_ESKOW_CG      =   6
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_GMPS_CG                =   7 
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_BAND_CG                =   8
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_LIN_MORE_CG            =   9
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_MULTIFRONTAL           =   11
      INTEGER, PUBLIC, PARAMETER :: GALAHAD_MODIFIED_MULTIFRONTAL  =   12

!-------------------------------------------------------------------------------
!                           PRECONDITIONERS
!-------------------------------------------------------------------------------

!     Note: the values in this series should be different from GALAHAD_NONE 
!           GALAHAD_and USER_DEFINED

      INTEGER, PUBLIC, PARAMETER :: GALAHAD_BANDED                 =   1

   END MODULE GALAHAD_SYMBOLS

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*   END SYMBOLS  M O D U L E  *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
