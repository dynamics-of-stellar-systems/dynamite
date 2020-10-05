!     ( Last modified on 23 Dec 2000 at 22:01:38 )
Program COBMA
  !
  !  COBYLA test driver for problems derived from SIF files.
  !
  !  A. R. Conn and Ph. Toint (based upon Nick Gould's vf13ma.f)
  !  January 1995.
  !
  !  Fortran 90/95 version, D. Orban, December 2006
  !
  Use CUTEr_precis
  Use CUTEr_interfaces
  Use CUTEr_Data

  Implicit None

  Type( CUTEr_Data_type ) :: Data

  Integer :: M, MAXFUN, LW, LIW
  Integer :: IPRINT, I, MGEQ, NFIX
  Integer, Dimension(:), Allocatable :: IW

  Real( Kind = wp ), Dimension(:), Allocatable :: W
  Real( Kind = wp ) :: RHOBEG, RHOEND, F
  Real( Kind = sp ), Dimension(2) :: CPU
  Real( Kind = sp ), Dimension(7) :: CALLS
  Integer :: ierr
  Integer, Parameter :: INPUT = 55, INDR = 46, IOUT = 6
  !
  !  Open the relevant file.
  !
  Open( INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
  Rewind INPUT
  !
  !  Initialize problem data structure
  !
  Data%allocate_H = .False.       ! No need for Hessian of objective/Lagrangian
  Data%allocate_J = .False.       ! No need for Jacobian of constraints
  Call CUTEr_Data_Initialize( Data, INPUT )
  !
  !  Set up the data structures necessary to hold the problem functions.
  !
  Call CSETUP( INPUT, IOUT, Data%n, Data%m, Data%x, Data%x_l, Data%x_u, &
       Data%n, Data%equation, Data%linear, Data%y, Data%c_l, Data%c_u,  &
       Data%m, .True., .False., .False. )
  Close( INPUT )
  !
  !  Allocate temporary work arrays
  !
  LW  = Data%n * ( 3 * Data%n + 2 * Data%m + 11 ) + 4 * Data%m + 6
  Allocate( W( LW ), STAT = ierr )
  If( ierr /= 0 ) Then
     Write(6,*) 'Error allocating memory'
  End If
  LIW = Data%n + 1
  Allocate( IW( LIW ), STAT = ierr )
  If( ierr /= 0 ) Then
     Write(6,*) 'Error allocating memory'
  End If
  !
  !  Count the number of general equality constraints and ignore them by
  !  shifting the remaining constraints at the beginning of the constraint list
  !
  MGEQ = 0
  Do I = 1, Data%m
     If ( Data%equation( I ) ) MGEQ = MGEQ + 1
  End Do
  If ( MGEQ > 0 ) Then
     Write( 6, 3090 ) MGEQ
     Do I = Data%m, MGEQ + 1, -1
        Data%c_u( I - MGEQ ) = Data%c_u( I )
        Data%c_l( I - MGEQ ) = Data%c_l( I )
     End Do
  End If
  M = Data%m - MGEQ
  !
  !  If constraints have both lower and upper bounds, they have to be
  !  included twice!
  !
  Do I = 1, Data%m - MGEQ
     If ( Data%c_l( I ) > -INFTY .And. Data%c_u( I ) < INFTY ) M = M + 1
  End Do
  !
  !  Include any simple bounds.
  !
  NFIX = 0
  Do I = 1, Data%n
     If ( Data%x_l( I ) ==  Data%x_u( I ) ) Then
        NFIX = NFIX + 1
     Else
        If ( Data%x_l( I ) > -INFTY ) M = M + 1
        If ( Data%x_u( I ) <  INFTY ) M = M + 1
     End If
  End Do
  If ( NFIX > 0 ) Write( 6, 3020 ) NFIX
  !
  !  Open the Spec file for the method.
  !
  Open( INDR, FILE = 'COBYLA.SPC', FORM = 'FORMATTED', STATUS = 'OLD')
  Rewind INDR
  !
  !  Read input Spec data.
  !
  !  RHOBEG = the size of the simplex initially
  !  RHOEND = the size of the simplex at termination
  !  MAXFUN = the maximum number of function calls allowed.
  !  IPRINT   should be set to 0, 1, 2 or 3, it controls the amount of printing
  !
  !  Set up algorithmic input data.
  !
  Read ( INDR, 1000 ) RHOBEG, RHOEND, MAXFUN, IPRINT
  Close ( INDR )
  !
  !  Evaluate the objective function and constraints.
  !
  Call CALCFC( Data, F, MGEQ )
  !
  !  Perform the minimization.
  !
  Call COBYLA( Data%n, M, Data%x, RHOBEG, RHOEND, IPRINT, MAXFUN, W, IW)
  !
  !  Output report
  !
  Call CREPRT ( CALLS, CPU )
  !
  Call CNAMES( Data%n, Data%m, Data%pname, Data%vnames, Data%cnames )
  Call CALCFC( Data, F, MGEQ )
  Write( 6, 2110 ) ( I, Data%vnames( I ), Data%x( I ), &
       Data%x_l( I ), Data%x_u( I ), I = 1, Data%n )
  If ( Data%m > 0 ) Write( 6, 2120 ) ( I, Data%cnames( I ), Data%c( I ), &
       Data%c_l( I ), Data%c_u( I ), Data%linear( I ), I = 1, Data%m )
  Write( 6, 2000 ) Data%pname, Data%n, Data%m, CALLS(1), CALLS(5), F, &
       CPU(1), CPU(2)
  !
  !  Clean-up data structures
  !
  Call CUTEr_Data_Free( Data )
  Deallocate( IW, STAT = ierr )
  Deallocate( W, STAT = ierr )
  Stop
  !
  !  Non-executable statements
  !
2000 Format( /, 24('*'), ' CUTEr statistics ', 24('*') //, &
          ' Code used               :  COBYLA ',  /,    &
          ' Problem                 :  ', A10,    /,    &
          ' # variables             =      ', I10 /,    &
          ' # constraints           =      ', I10 /,    &
          ' # objective functions   =        ', F8.2 /, &
          ' # constraints functions =        ', F8.2 /, &
          ' Final f                 = ', E15.7 /,       &
          ' Set up time             =      ', 0P, F10.2, ' seconds' /, &
          ' Solve time              =      ', 0P, F10.2, ' seconds' //, &
          66('*') / )
1000 Format( D12.4, /, D12.4, /,I6, /, I6 )
2110 Format( /, ' The variables:', /, &
          '     I name          value    lower bound upper bound', &
          /, ( I6, 1X, A10, 1P, 3D12.4 ) )
2120 Format( /, ' The constraints:', /, &
          '     I name          value    lower bound upper bound', &
          ' linear? ', &
          /, ( I6, 1X, A10, 1P, 3D12.4, 5X, L1 ) )
3000 Format( /,'  ** Program CSETUP: array length ', A6, ' too small.', &
          /,'  -- Miminimization abandoned.', &
          /,'  -- Increase the parameter ', A6, ' by at least ', I8, &
            ' and restart.'  )
3020 Format( /,'  ** Warning from COBMA. **', &
          /,'     In the problem as stated , ', I6, &
            ' variables are fixed: they are changed to free.' )
3090 Format( /,'  ** Warning from COBMA. **', &
          /,'     The problem as stated includes ', I6, &
            ' equality constraints: they are ignored ' )
  !
  !  End of COBMA.
  !
End Program COBMA
Subroutine CALCFC( Data, F, MGEQ )
  !
  !  Evaluates the objective function value in a format compatible with COBYLA,
  !  but using the CUTEr tools.
  !
  !  A. R. Conn and Ph. Toint
  !  January 1995.
  !
  !  Fortran 90/95 version, D. Orban, December 2006
  !
  Use CUTEr_precis
  Use CUTEr_Data
  Implicit None

  Type( CUTEr_Data_type ), Intent( Inout ) :: Data
  Real( Kind = wp ), Intent( Inout ) :: F
  Integer, Intent( In ) :: MGEQ
  Integer :: I, MT
  !
  !  Evaluate the objective function and constraints.
  !
  Call CFN( Data%n, Data%m, Data%x, F, Data%m, Data%c )
  !
  !  If there are equality constraints, ignore them
  !  and shift all the inequality constraint values.
  !
  If ( MGEQ > 0 ) Then
     Do I = Data%m, MGEQ + 1, - 1
        Data%c( I - MGEQ ) = Data%c( I )
     End Do
  End If
  !
  !  If constraints have both lower and upper bounds, they have to
  !  be included twice! Reverse the signs of less-than-or-equal-to
  !  constraints.
  !
  MT = Data%m - MGEQ + 1
  Do I = 1, Data%m - MGEQ
     If ( Data%c_l( I ) > -INFTY .And. Data%c_u( I ) < INFTY ) Then
        Data%c( I )  = Data%c_u( I ) - Data%c( I )
        Data%c( MT ) = Data%c( I ) - Data%c_l( I )
        MT      = MT + 1
     Else If ( Data%c_l( I ) > -INFTY ) Then
        Data%c( I )  = Data%c( I ) - Data%c_l( I )
     Else If ( Data%c_u( I ) < INFTY ) Then
        Data%c( I )  = Data%c_u( I ) - Data%c( I )
     End If
  End Do
  !
  !  Include any simple bounds, including fixed variables.
  !
  Do I = 1, Data%n
     If ( Data%x_l( I ) /=  Data%x_u( I ) ) Then
        If ( Data%x_l( I ) > -INFTY ) Then
           Data%c( MT ) = Data%x( I ) - Data%x_l( I )
           MT      = MT + 1
        End If
        If ( Data%x_u( I ) < INFTY ) Then
           Data%c( MT ) = Data%x_u( I ) - Data%x( I )
           MT      = MT + 1
        End If
     End If
  End Do
  Return
  !
  !  End of CALCFC.
  !
End Subroutine CALCFC
