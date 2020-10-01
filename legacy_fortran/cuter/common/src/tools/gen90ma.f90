!     ( Last modified on 23 Dec 2000 at 22:01:38 )
Program GENMA
  !
  Use CUTEr_precis
  Use CUTEr_interfaces
  Use Generic_Driver
  !
  !  Generic package driver (example) for applying package GEN90 to problems
  !  from SIF files. This driver also demonstrates how to dynamically
  !  allocate arrays to be used with CUTEr.
  !
  !  D. Orban, August 2002, strongly inspired by Philippe's original driver.
  !
  Implicit None
  Integer :: N, M, NLIN, NEQ, NBNDS, EXITCODE
  Integer, Parameter :: INSPEC = 46, INPUT = 47, IOUT = 6
  Real( Kind = wp ) :: DUMMY
  Real( Kind = wp ), Dimension( : ), Allocatable :: X, BL, BU, V, CL, CU
  Real, Dimension( 2 ) :: CPU( 2 )
  Real, Dimension( 7 ) :: CALLS( 7 )
  Character( len = 10 ) ::  PNAME
  Character( len = 10 ), Dimension( : ), Allocatable :: VNAMES, GNAMES
  Logical :: EFIRST, LFIRST, NVFRST
  Logical, Dimension( : ), Allocatable :: EQUATN, LINEAR
  Logical ::  CONSTRAINED
  !
  !  Open the Spec file for the method (typically called METHOD.SPC)
  !
  Call GENSPC( INSPEC, 'GEN.SPC' )
  !
  !  Open the relevant problem file.
  !
  Open ( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
  Rewind INPUT
  !
  !  Get problem dimensions and determine which tools to use
  !
  CONSTRAINED = .False.
  Call CDIMEN( INPUT, N, M )
  If( M > 0 ) Then
     CONSTRAINED = .True.
  Else If( M < 0 ) Then
     Write( 6, '(A)' ) 'Error reading OUTSDIF.d'
     Stop
  Endif
  !
  !  Set up parameters
  !
  EFIRST = .True. ; LFIRST = .False. ; NVFRST = .False.
  !
  !  Set up SIF data from the problem file
  !
  Allocate( X( N ), BL( N ), BU( N ) )
  If( CONSTRAINED ) Then
     Allocate( V( M+1 ), CL( M+1 ), CU( M+1 ), EQUATN( M+1 ), LINEAR( M+1 ) )
     Call CSETUP( INPUT, IOUT, N, M, X, BL, BU, N, EQUATN, &
          LINEAR, V, CL, CU, M+1, EFIRST, LFIRST, NVFRST )
  Else
     Allocate( EQUATN( 1 ), LINEAR( 1 ) )
     Call USETUP( INPUT, IOUT, N, X, BL, BU, N )
  Endif
  !
  !  Obtain problem/variables/constraints names.
  !
  Allocate( VNAMES( N ) )
  If( CONSTRAINED ) Then
     Allocate( GNAMES( M ) )
    Call CNAMES( N, M, PNAME, VNAMES, GNAMES )
  Else
     Call UNAMES( N, PNAME, VNAMES )
  Endif
!
!  Obtain info on the problem
!
  NLIN  = 0 ; NEQ   = 0 ; NBNDS = 0
  If( CONSTRAINED ) Then
     Call GETINFO( N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
  Else
     EQUATN( 1 ) = .False.
     LINEAR( 1 ) = .False.
     Call GETINFO( N, 1, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
  Endif
  !
  !  Call the optimizer.
  !
  Call GEN( DUMMY )
  EXITCODE = 0
  !
  !  Close the problem file
  !
  Close( INPUT )
  !
  !  Write the standard statistics (of which some may be irrelevant)
  !
  !    CALLS( 1 ): number of calls to the objective function
  !    CALLS( 2 ): number of calls to the objective gradient
  !    CALLS( 3 ): number of calls to the objective Hessian
  !    CALLS( 4 ): number of Hessian times vector products
  !           --constrained problems only--
  !    CALLS( 5 ): number of calls to the constraint functions
  !    CALLS( 6 ): number of calls to the constraint gradients
  !    CALLS( 7 ): number of calls to the constraint Hessians
  !           -----------------------------
  !
  !    CPU( 1 ) : CPU time (in seconds) for USETUP or CSETUP
  !    CPU( 2 ) : CPU time ( in seconds) since the end of USETUP or CSETUP
  !
  !  Note that each constraint function is counted separately.
  !  Evaluating all the constraints thus results in PNC evaluations, where
  !  PNC is the number of constraints in the problem.  Note that PNC does not
  !  include repetitions for constraints having full ranges.
  
  !  (N, is the dimension of the problem, M is the number of constraints,
  !   DUMMY is the final value of the objective function)
  !
  If( CONSTRAINED ) Then
     Call CREPRT( CALLS, CPU )      
  Else
     Call UREPRT( CALLS, CPU )
  Endif
  Write ( IOUT, 2000 ) PNAME, N, M, NLIN, NEQ, M-NEQ, NBNDS,       &
       CALLS( 1 ), CALLS( 2 ), CALLS( 3 ), CALLS( 5 ), CALLS( 6 ), &
       CALLS( 7 )
  Write ( IOUT, 2001 ) EXITCODE, DUMMY, CPU( 1 ), CPU( 2 )
  !
  !  Free allocated memory
  !
  Deallocate( X, BU, BL, VNAMES, EQUATN, LINEAR )
  If( CONSTRAINED ) Deallocate( V, CL, CU, GNAMES )
  !
  !  Exit
  !
  Stop
  !
  !  Non-executable statements.
  !
  !  The following is the complete standard statistics output format: select
  !  the items that are relevant to the type of problems solved and adapt the
  !  name of the code. It is broken in two to comply with compilers
  !  which want to see no more than 19 continuation lines.
  !
2000 Format( /, 24('*'), ' CUTEr statistics ', 24('*') // &
          ,' Code used                :  GEN90',    / &
          ,' Variant                  :  name of a variant, if needed',/ &
          ,' Problem                  :  ', A10,    / &
          ,' # variables              =      ', I10 / &
          ,' # constraints            =      ', I10 / &
          ,' # linear constraints     =      ', I10 / &
          ,' # equality constraints   =      ', I10 / &
          ,' # inequality constraints =      ', I10 / &
          ,' # bounds                 =      ', I10 / &
          ,' # objective functions    =        ', F8.2 / &
          ,' # objective gradients    =        ', F8.2 / &
          ,' # objective Hessians     =        ', F8.2 / &
          ,' # constraints functions  =        ', F8.2 / &
          ,' # constraints gradients  =        ', F8.2 / &
          ,' # constraints Hessians   =        ', F8.2 )
2001 Format(                                          &
          ' Exit code                =      ', I10 / &
          ,' Final f                  = ', E15.7 / &
          ,' Set up time              =      ', 0P, F10.2, ' seconds'/ &
          ' Solve time               =      ', 0P, F10.2, ' seconds'// &
          66('*') / )
End Program GENMA

