! CUTEr derived data type
Module CUTEr_Data

  Use CUTEr_precis

  Implicit None

  Private
  Public :: CUTEr_Data_Initialize, CUTEr_Data_Free

  Integer, Parameter :: error_device = 6

  !============================================================================

  Type, Public :: CUTEr_Data_type

     Integer :: n = 0          ! Number of variables
     Integer :: m = 0          ! Number of general constraints

     Character( Len = 10 ) :: pname    ! Problem name
     Character( Len = 10 ), Pointer, Dimension(:) :: vnames => NULL()
     Character( Len = 10 ), Pointer, Dimension(:) :: cnames => NULL()

     ! Vector of variables and lower and upper bounds
     Real( Kind = wp ), Pointer, Dimension(:) :: x   => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: x_l => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: x_u => NULL()

     ! Vector of dual variables associated to bound constraints
     Real( Kind = wp ), Pointer, Dimension(:) :: z => NULL()

     ! Vector of constraints and lower and upper bounds
     Real( Kind = wp ), Pointer, Dimension(:) :: c   => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: c_l => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: c_u => NULL()

     ! Lagrange multipliers associated with general constraints
     Real( Kind = wp ), Pointer, Dimension(:) :: y => NULL()

     ! Constraint indicators
     Logical, Pointer, Dimension(:) :: equation => NULL()
     Logical, Pointer, Dimension(:) :: linear   => NULL()

     ! Gradient of the Lagrangian
     Real( Kind = wp ), Pointer, Dimension(:) :: g

     ! Pointers for sparse Hessian storage (coordinate)
     Logical :: allocate_H = .False.
     Integer :: nnzh = 0
     Integer, Pointer, Dimension(:) :: H_row => NULL()
     Integer, Pointer, Dimension(:) :: H_col => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: H_val => NULL()

     ! Pointers for sparse Jacobian storage (coordinate)
     Logical :: allocate_J = .False.
     Integer :: nnzj = 0
     Integer, Pointer, Dimension(:) :: J_row => NULL()
     Integer, Pointer, Dimension(:) :: J_col => NULL()
     Real( Kind = wp ), Pointer, Dimension(:) :: J_val => NULL()

  End Type CUTEr_Data_type

Contains

  !============================================================================

  Subroutine CUTEr_Data_Initialize( problem, input )
    !
    ! Allocate main problem data structure
    !
    Implicit None
    !
    Type( CUTEr_Data_type ), Intent( Inout ) :: problem
    Integer, Intent( In ) :: input
    Integer :: ierr, nerr
    !
    ! Obtain problem dimensions
    ! 
    Call CDIMEN( input, problem%n, problem%m )
    !
    ! Trap invalid input
    !
    If( problem%n <= 0 .Or. problem%m < 0 ) Then
       Write( error_device, 10 ) problem%n, problem%m
       Call CUTEr_Data_Free( problem )
       Return
    Endif
    !
    ! Allocate arrays depending on n
    !
    ierr = 0
    nerr = 0
    Allocate( problem%vnames( problem%n ), STAT = ierr )
    If( ierr /= 0 ) Then
       nerr = nerr + 1
    End If
    Allocate( problem%x( problem%n ), STAT = ierr )
    If( ierr /= 0 ) Then
       nerr = nerr + 1
    End If
    Allocate( problem%x_l( problem%n ), STAT = ierr )
    If( ierr /= 0 ) Then
       nerr = nerr + 1
    End If
    Allocate( problem%x_u( problem%n ), STAT = ierr )
    If( ierr /= 0 ) Then
       nerr = nerr + 1
    End If
    Allocate( problem%z( problem%n ), STAT = ierr )
    If( ierr /= 0 ) Then
       nerr = nerr + 1
    End If
    !
    ! Allocate memory for Hessian if required
    !
    If( problem%allocate_H ) Then
       Call CDIMSH( problem%nnzh )
       Allocate( problem%H_row( problem%nnzh ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%H_col( problem%nnzh ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%H_val( problem%nnzh ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
    End If
    !
    ! Allocate arrays depending on m
    !
    If( problem%m > 0 ) Then
       Allocate( problem%cnames( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%c( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%c_l( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%c_u( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%y( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%equation( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       Allocate( problem%linear( problem%m ), STAT = ierr )
       If( ierr /= 0 ) Then
          nerr = nerr + 1
       End If
       !
       ! Allocate memory for Jacobian if required
       !
       If( problem%allocate_J ) Then
          Call CDIMSJ( problem%nnzj )
          Allocate( problem%J_row( problem%nnzj ), STAT = ierr )
          If( ierr /= 0 ) Then
             nerr = nerr + 1
          End If
          Allocate( problem%J_col( problem%nnzj ), STAT = ierr )
          If( ierr /= 0 ) Then
             nerr = nerr + 1
          End If
          Allocate( problem%J_val( problem%nnzj ), STAT = ierr )
          If( ierr /= 0 ) Then
             nerr = nerr + 1
          End If
       End If
    End If
    !
    ! Report errors if any
    !
    If( nerr > 0 ) Write( error_device, 20 ) nerr
    !
    ! Non-executable statements
    !
10  Format( 'CUTEr_Data_Initialize:: invalid problem dimensions:', &
         'n = ', I6, ', m = ', I6 )
20  Format( 'CUTEr_Data_Initialize:: ', I2, ' errors in memory allocation' )

  End Subroutine CUTEr_Data_Initialize

  !============================================================================

  Subroutine CUTEr_Data_Free( problem )

    ! Deallocates dynamically-allocated memory for problem storage

    Implicit None

    Type( CUTEr_Data_type ), Intent( Inout ) :: problem
    Integer :: ierr, nerr

    nerr = 0

    If( Associated( problem%vnames ) ) Then
       Deallocate( problem%vnames, STAT = ierr )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%cnames ) ) Then
       Deallocate( problem%cnames )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%x ) ) Then
       Deallocate( problem%x )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%x_l ) ) Then
       Deallocate( problem%x_l )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%x_u ) ) Then
       Deallocate( problem%x_u )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%z ) ) Then
       Deallocate( problem%z )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%c ) ) Then
       Deallocate( problem%c )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%c_l ) ) Then
       Deallocate( problem%c_l )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%c_u ) ) Then
       Deallocate( problem%c_u )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%y ) ) Then
       Deallocate( problem%y )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%equation ) ) Then
       Deallocate( problem%equation )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%linear ) ) Then
       Deallocate( problem%linear )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%g ) ) Then
       Deallocate( problem%g )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%H_row ) ) Then
       Deallocate( problem%H_row )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%H_col ) ) Then
       Deallocate( problem%H_col )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%H_val ) ) Then
       Deallocate( problem%H_val )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%J_row ) ) Then
       Deallocate( problem%J_row )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%J_col ) ) Then
       Deallocate( problem%J_col )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif
    If( Associated( problem%J_val ) ) Then
       Deallocate( problem%J_val )
       If( ierr /= 0 ) nerr = nerr + 1
    Endif

    Write( error_device, 100 ) nerr

100 Format( I9, ' errors encountered in freeing memory' )

  End Subroutine CUTEr_Data_Free

  !============================================================================

End Module CUTEr_Data
