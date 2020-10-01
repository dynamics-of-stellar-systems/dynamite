! THIS VERSION: GALAHAD 2.2 - 21/04/2008 AT 16:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ F D C    M O D U L E  -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.0.August 14th 2006

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_FDC_double

!     -----------------------------------------------------
!     |                                                   |
!     | Check if A x = c is consistent and if so which    |
!     | if any of the constraints are linearly dependeent |
!     |                                                   |
!     -----------------------------------------------------

!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SYMBOLS
      USE GALAHAD_SPACE_double
      USE GALAHAD_SMT_double
      USE GALAHAD_SILS_double
      USE GALAHAD_ROOTS_double
      USE GALAHAD_SPECFILE_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: FDC_initialize, FDC_read_specfile, FDC_find_dependent,          &
                FDC_terminate, SMT_type, SILS_factors

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: FDC_time_type
        REAL :: total, analyse, factorize
      END TYPE

      TYPE, PUBLIC :: FDC_control_type
        INTEGER :: error, out, print_level, indmin, valmin
        REAL ( KIND = wp ) :: pivot_tol, zero_pivot, max_infeas
        LOGICAL :: deallocate_error_fatal, space_critical, scale
        CHARACTER ( LEN = 30 ) :: prefix
        TYPE ( SILS_control ) :: CNTL
      END TYPE

      TYPE, PUBLIC :: FDC_inform_type
        INTEGER :: status, alloc_status, factorization_status
        INTEGER :: factorization_integer, factorization_real
        REAL ( KIND = wp ) :: non_negligible_pivot
        CHARACTER ( LEN = 80 ) :: bad_alloc
        TYPE ( FDC_time_type ) :: time
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      REAL ( KIND = wp ), PARAMETER :: half = 0.5_wp
      REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
      REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
      REAL ( KIND = wp ), PARAMETER :: epsmch = EPSILON( one )

   CONTAINS

!-*-*-*-*-*-   F D C _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE FDC_initialize( control, FACTORS )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for FDC. This routine should be called before
!  FDC_primal_dual
! 
!  --------------------------------------------------------------------
!
!  Arguments:
!
!  data     private internal data
!  control  a structure containing control information. Components are -
!
!  INTEGER control parameters:
!
!   error. Error and warning diagnostics occur on stream error 
!   
!   out. General output occurs on stream out
!   
!   print_level. The level of output required is specified by print_level
!   
!   indmin. An initial guess as to the integer workspace required by SILS
!
!   valmin. An initial guess as to the real workspace required by SILS
! 
!  REAL control parameters:
!
!   pivot_tol. The threshold pivot used by the matrix factorization when 
!    attempting to detect linearly dependent constraints.
!    See the documentation for SBLS for details
!
!   zero_pivot. Any pivots smaller than zero_pivot in absolute value will 
!    be regarded to be zero when attempting to detect linearly dependent 
!    constraints
!
!   max_infeas. The largest permitted residual
!
!  LOGICAL control parameters:
!
!
!  scale. If true, the rows of A will be scaled to have unit infinity norm.
!   Otherwise, no scaling will be applied

!   deallocate_error_fatal. If true, any array/pointer deallocation error
!     will terminate execution. Otherwise, computation will continue
!
!  CHARACTER control parameters:
!
!  prefix (len=30). All output lines will be prefixed by 
!    %prefix(2:LEN(TRIM(%prefix))-1)
!   where %prefix contains the required string enclosed in 
!   quotes, e.g. "string" or 'string'
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

      TYPE ( FDC_control_type ), INTENT( OUT ) :: control        
      TYPE ( SILS_factors ), OPTIONAL :: FACTORS

!  If required, initialize FACTORS

      IF ( PRESENT( FACTORS ) ) CALL SILS_initialize( FACTORS = FACTORS )

!  Set control parameters

!  Integer parameters

      control%error = 6
      control%out = 6
      control%print_level = 0
      control%indmin = 1000
      control%valmin = 1000

!  Real parameters

      control%pivot_tol = half
      control%zero_pivot = epsmch ** 0.75
      control%max_infeas = epsmch ** 0.33

!  Logical parameters

      control%scale = .FALSE.
      control%space_critical = .FALSE.
      control%deallocate_error_fatal  = .FALSE.

!  Character parameters

      control%prefix = '""                            '

!  Initalize SILS components

      CALL SILS_initialize( control = control%CNTL )
      control%CNTL%u = control%pivot_tol
      control%CNTL%tolerance = epsmch ** 2
      control%CNTL%pivoting = 1
      control%CNTL%lp = - 1
      control%CNTL%mp = - 1
      control%CNTL%wp = - 1
      control%CNTL%ordering = 3
!57V2 control%CNTL%ordering = 2
!57V3 control%CNTL%ordering = 5
!57V2 control%CNTL%scaling = 0
!57V2 control%CNTL%static_tolerance = zero
!57V2 control%CNTL%static_level = zero

      RETURN  

!  End of FDC_initialize

      END SUBROUTINE FDC_initialize

!-*-*-*-*-   F D C _ R E A D _ S P E C F I L E  S U B R O U T I N E   -*-*-*-*-

      SUBROUTINE FDC_read_specfile( control, device, alt_specname )

!  Reads the content of a specification file, and performs the assignment of 
!  values associated with given keywords to the corresponding control parameters

!  The defauly values as given by FDC_initialize could (roughly) 
!  have been set as:

! BEGIN FDC SPECIFICATIONS (DEFAULT)
!  error-printout-device                             6
!  printout-device                                   6
!  print-level                                       0
!  initial-integer-workspace                         1000
!  initial-real-workspace                            1000
!  maximum-permitted-infeasibility                   1.0D-5
!  pivot-tolerance-used-for-dependencies             0.5
!  zero-pivot-tolerance                              1.0D-12
!  scale-A                                           F
!  space-critical                                    F
!  deallocate-error-fatal                            F
! END FDC SPECIFICATIONS (DEFAULT)

!  Dummy arguments

      TYPE ( FDC_control_type ), INTENT( INOUT ) :: control        
      INTEGER, INTENT( IN ) :: device
      CHARACTER( LEN = 16 ), OPTIONAL :: alt_specname

!  Programming: Nick Gould and Ph. Toint, January 2002.

!  Local variables

      INTEGER, PARAMETER :: lspec = 11
      CHARACTER( LEN = 16 ), PARAMETER :: specname = 'FDC            '
      TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec

!  Define the keywords

      spec( : )%keyword = ''

!  Integer key-words

      spec(  1 )%keyword = 'error-printout-device'
      spec(  2 )%keyword = 'printout-device'
      spec(  3 )%keyword = 'print-level' 
      spec(  4 )%keyword = 'initial-integer-workspace'
      spec(  5 )%keyword = 'initial-real-workspace'

!  Real key-words

      spec(  6 )%keyword = 'maximum-permitted-infeasibility'
      spec(  7 )%keyword = 'pivot-tolerance-used-for-dependencies'
      spec(  8 )%keyword = 'zero-pivot-tolerance'

!  Logical key-words

      spec( 11 )%keyword = 'scale-A'
      spec(  9 )%keyword = 'space-critical'
      spec( 10 )%keyword = 'deallocate-error-fatal'

!  Read the specfile

      IF ( PRESENT( alt_specname ) ) THEN
        CALL SPECFILE_read( device, alt_specname, spec, lspec, control%error )
      ELSE
        CALL SPECFILE_read( device, specname, spec, lspec, control%error )
      END IF

!  Interpret the result

!  Set integer values

      CALL SPECFILE_assign_value( spec( 1 ), control%error,                    &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 2 ), control%out,                      &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 3 ), control%print_level,              &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 4 ), control%indmin,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 5 ), control%valmin,                   &
                                  control%error )

!  Set real values


      CALL SPECFILE_assign_value( spec( 6 ), control%max_infeas,               &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 7 ), control%pivot_tol,                &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 8 ), control%zero_pivot,               &
                                  control%error )

!  Set logical values


      CALL SPECFILE_assign_value( spec( 11 ), control%scale,                   &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 9 ), control%space_critical,           &
                                  control%error )
      CALL SPECFILE_assign_value( spec( 10 ),                                  &
                                  control%deallocate_error_fatal,              &
                                  control%error )

      RETURN

      END SUBROUTINE FDC_read_specfile

!-*-*-*-   F D C _ F I N D _ D E P E N D E N T  S U B R O U T I N E   -*-*-*-

      SUBROUTINE FDC_find_dependent( n, m, A_val, A_col, A_ptr, C, K,          &
                                     FACTORS, n_depen, C_depen, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Determine which, if any, of the equality constraints A x = c
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Arguments:
!
!   n is an INTEGER variable, which must be set by the user to the
!    number of columns of A. RESTRICTION: n >= 1
!                 
!   m is an INTEGER variable, which must be set by the user to the
!    number of rows of A. RESTRICTION: m >= 0
!        
!   A_val is a REAL array of length A_ptr( m + 1 ) - 1 that must be set by 
!    the user to the values of the components of A, stored by row, that is
!    entries for row i must directly preceed those in row i+1 (the order
!    within each row is unimportant).
!
!   A_col is an INTEGER array of length A_ptr( m + 1 ) - 1 that must be set by 
!    the user to the column indices of the components of A corresponding
!    to the values in A_val
!
!   A_ptr is an INTEGER array of length that must be set by the user to point 
!    to the position in A_val and A_col for the start of each row, as well
!    as one position past the last entry.
!         
!   C is a REAL array of length m, which is used to store the values of c
!
!   K is a structure of type SMT_type that need not be set by the user.
!
!   FACTORS is a structure of type SILS_factors that need not be set by the user.
!
!   n_depen is an INTEGER variable that gives the number of rows of A that
!    are deemed to be linearly dependent
!
!   C_depen is a INTEGER pointer array that will have been allocated to
!    be of length n_depen, and will contain the indices of the rows of A
!    that are deemed to be linearly dependent
!
!   control is a structure of type FDC_control_type that controls the 
!   execution of the subroutine and must be set by the user. Default values for
!   the elements may be set by a call to FDC_initialize. See FDC_initialize 
!   for details
!
!   inform is a structure of type FDC_inform_type that provides 
!    information on exit from FDC_find_dependent. The component status 
!    has possible values:
!  
!     0 Normal termination with a prediction of how many (and which) 
!       constraints are dependent
!
!    -1 An allocation error occured; the status is given in the component
!       alloc_status.
!
!    -2 A deallocation error occured; the status is given in the component
!       alloc_status.
!
!    -3 One of the restrictions n >= 0 and m >= 0 has been violated
!
!    -5 The constraints appear to be inconsistent
!
!    -9 The ordering failed. The return status from the factorization
!      package is given in inform%factorization_status
!    
!   -10 The factorization failed. The return status from the factorization
!      package is given in inform%factorization_status
!    
!  On exit from FDC_find_dependent, other components of inform give the 
!  following:
!
!     alloc_status = The status of the last attempted allocation/deallocation 
!     factorization_status = The return status from the factorization
!     factorization_integer = The total integer workspace required by the 
!              factorization.
!     factorization_real = The total real workspace required by the 
!              factorization.
!     nfacts = The total number of factorizations performed.
!     factorization_status = the return status from the matrix factorization
!              package.   
!     non_negligible_pivot = the smallest pivot which was not judged to be
!       zero when detecting linearly dependent constraints
!     bad_alloc = the name of the array for which an allocation/deallocation
!       error ocurred
!     time%total = the total time spent in the package.
!     time%analyse = the time spent analysing the required matrices prior to
!                  factorization.
!     time%factorize = the time spent factorizing the required matrices.
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( IN ) :: n, m
      INTEGER, INTENT( OUT ) :: n_depen
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      INTEGER, INTENT( IN ), DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_col
      REAL ( KIND = wp ), INTENT( IN ),                                        &
                          DIMENSION( A_ptr( m + 1 ) - 1 ) :: A_val
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( m ) :: C

!  pointers (which will be allocatable arrays in f2000!!)

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_depen
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( FDC_control_type ), INTENT( IN ) :: control        
      TYPE ( FDC_inform_type ), INTENT( INOUT ) :: inform

!  Local variables

      INTEGER :: A_ne, i, ii, l, nroots, out, pmax, pmin
      REAL :: time_start, time_now, time_record
      REAL ( KIND = wp ) :: root1, root2, dmax, dmin, dmax_allowed
      REAL ( KIND = wp ) :: big, res, res_max, rmax, rmin

      LOGICAL ::  twobytwo
      CHARACTER ( LEN = 80 ) :: array_name
      TYPE ( SILS_ainfo ) :: AINFO
      TYPE ( SILS_finfo ) :: FINFO
      TYPE ( SILS_sinfo ) :: SINFO
      TYPE ( SILS_control ) :: CNTL

!  Automatic arrays

      INTEGER, DIMENSION( n + m ) :: P
      REAL ( KIND = wp ), DIMENSION( m ) :: SCALE
      REAL ( KIND = wp ), DIMENSION( n + m ) :: SOL
      REAL ( KIND = wp ), DIMENSION( 2, n + m ) :: D

!  Initialize return information

      inform%alloc_status = 0
      inform%factorization_status = 0
      inform%factorization_integer = 0
      inform%factorization_real = 0
      inform%non_negligible_pivot = zero
      inform%bad_alloc = ''

!  Set initial timing breakdowns

      inform%time%total = 0.0 ; inform%time%analyse = 0.0
      inform%time%factorize = 0.0

      out = control%out
      CNTL = control%CNTL

      IF ( out > 0 .AND. control%print_level >= 1 ) WRITE( out,                &
         "( /, 7( ' -' ), ' test for rank defficiency', 7( ' - ' ) )" )

!  Check that the problem makes sense

      IF ( n < 0 .OR. m < 0 ) THEN
        inform%status = GALAHAD_error_restrictions ; RETURN ; END IF 

!  Check the case where n is zero

      IF ( n == 0 ) THEN
        IF ( MAXVAL( C ) <= control%max_infeas ) THEN
          DO i = 1, m - 1 
            C_depen( i ) = i
          END DO
          n_depen = m - 1
          inform%status = GALAHAD_ok
        ELSE
          inform%status = GALAHAD_error_primal_infeasible
        END IF
        RETURN
      END IF

!  Analyse the sparsity pattern

      CALL CPU_TIME( time_start ) 

!  Set the dimensions of K

      A_ne = A_ptr( m + 1 ) - 1
      K%n = n + m ; K%ne = A_ne + n

!  Allocate the arrays for the analysis phase

      array_name = 'fdc: data%K%row'
      CALL SPACE_resize_array( K%ne, K%row, inform%status,                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      array_name = 'fdc: data%K%col'
      CALL SPACE_resize_array( K%ne, K%col, inform%status,                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

      array_name = 'fdc: data%K%val'
      CALL SPACE_resize_array( K%ne, K%val, inform%status,                     &
             inform%alloc_status, array_name = array_name,                     &
             deallocate_error_fatal = control%deallocate_error_fatal,          &
             exact_size = control%space_critical,                              &
             bad_alloc = inform%bad_alloc, out = control%error )
      IF ( inform%status /= GALAHAD_ok ) RETURN

!  Put A into K

      DO i = 1, m
        ii = n + i
        IF ( control%scale ) THEN
          IF ( A_ptr( i + 1 ) - 1 >= A_ptr( i ) ) THEN
            SCALE( i ) = MAX( one,                                             &
              MAXVAL( ABS( A_val( A_ptr( i ) : A_ptr( i + 1 ) - 1 ) ) ) )
          ELSE
            SCALE( i ) = one
          END IF
        ELSE
          SCALE( i ) = one
        END IF
        DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
          K%row( l ) = ii ; K%col( l ) = A_col( l )
          K%val( l ) = A_val( l ) / SCALE( i )
        END DO
      END DO
!     dmax = MAXVAL( ABS( A_val( : A_ne ) ) )
!     dmax = MAX( dmax, one )
!     dmax = one
      dmax = ten ** ( - 2 )

!  Put the diagonal into K

      DO i = 1, n
        K%row( A_ne + i ) = i ; K%col( A_ne + i ) = i
      END DO
      K%val( A_ne + 1 : A_ne + n ) = dmax

!     write(78,*) K%n, K%ne
!     DO i = 1,  K%ne
!       write(78,"( 2I8, ES22.14 )" ) K%row( i ), K%col( i ), K%val( i )
!     END DO

!  Analyse the sparsity pattern of the matrix

      CALL SILS_analyse( K, FACTORS, CNTL, AINFO )

!  Record the storage requested

      inform%factorization_integer = AINFO%nirnec 
      inform%factorization_real = AINFO%nrlnec

!  Check for error returns

      inform%factorization_status = AINFO%flag
      IF ( AINFO%flag < 0 ) THEN
        IF ( control%error > 0 .AND. control%print_level >= 1 )                &
           WRITE( control%error,                                               &
           "( '   **  Error return ', I6, ' from SILS_analyse' )") AINFO%flag
        inform%status = GALAHAD_error_analysis ; RETURN
      ELSE IF ( AINFO%flag > 0 ) THEN 
        IF ( out > 0 .AND. control%print_level >= 1 )                          &
          WRITE( out,                                                          &
           "( /, ' ** Warning ', I8, ' from SILS_analyse' )" ) AINFO%flag
      END IF

      IF ( out > 0 .AND. control%print_level >= 2 ) WRITE( out,                &
          "( ' real/integer space required for factors ', 2I10 )" )            &
          AINFO%nrladu, AINFO%niradu

      CALL CPU_TIME( time_now ) ; inform%time%analyse = time_now - time_start
      IF ( out > 0 .AND. control%print_level >= 2 )                            &
        WRITE( out, "( ' ** analysis time = ', F10.2 ) " ) inform%time%analyse

!  Factorize the matrix

!     CNTL%lp = 6
!     CNTL%mp = 6
!     CNTL%wp = 6
!     CNTL%sp = 6
!     CNTL%ldiag = 3
      CNTL%pivoting = 1
      CNTL%liw = MAX( 2 * inform%factorization_integer, control%indmin )
      CNTL%la = MAX( 2 * inform%factorization_real, control%valmin )

      CALL CPU_TIME( time_record ) 
      CALL SILS_factorize( K, FACTORS, CNTL, FINFO )

      CALL CPU_TIME( time_now ) ; inform%time%factorize = time_now - time_record
      IF ( out > 0 .AND. control%print_level >= 2 )                            &
        WRITE( out, "( ' ** factorize time = ', F10.2 ) " ) inform%time%factorize

!  Record the storage required

      inform%factorization_integer = FINFO%nirbdu 
      inform%factorization_real = FINFO%nrlbdu 

!  Test that the factorization succeeded

      inform%factorization_status = FINFO%flag
      IF ( FINFO%flag < 0 ) THEN
        IF ( control%error > 0 .AND. control%print_level >= 1 )                &
           WRITE( control%error,                                               &
            "( '   **  Error return ', I6, ' from SILS_factorize' )") FINFO%flag
        inform%status = GALAHAD_error_factorization
        RETURN

!  Record warning conditions

      ELSE IF ( FINFO%flag > 0 ) THEN
        IF ( out > 0 .AND. control%print_level >= 2 ) THEN
          WRITE( out,                                                          &
           "( /, ' ** Warning ', I8, ' from SILS_factorize' )" ) FINFO%flag
          IF ( FINFO%flag == 4 ) WRITE( out, "( ' ** Matrix has ', I7,         &
         &   ' zero eigenvalues ' )" ) K%n - finfo%rank
        END IF 
      END IF 

!  Determine the block diagonal part of the factors and the pivot order

      CALL SILS_enquire( FACTORS, PIVOTS = P, D = D )

!  Compute the smallest and largest eigenvalues of the block diagonal factor

      n_depen = 0 ; twobytwo = .FALSE. 
      dmax = zero ; dmin = HUGE( one )
      dmax_allowed = zero
      big = one / MAX( control%zero_pivot, epsmch )

!  Loop over the diagonal blocks
  
      DO i = 1, finfo%rank
        IF ( twobytwo ) THEN
          twobytwo = .FALSE.
          CYCLE
        END IF
        IF ( i < finfo%rank ) THEN

!  A 2x2 block

          IF ( P( i ) < 0 ) THEN
            twobytwo = .TRUE.

            CALL ROOTS_quadratic( D( 1, i ) * D( 1, i + 1 ) - D( 2, i ) ** 2,  &
               - D( 1, i ) - D( 1, i + 1 ), one, epsmch, nroots, root1, root2 ) 
            rmax = MAX( ABS( root1 ), ABS( root2 ) )
            rmin = MIN( ABS( root1 ), ABS( root2 ) )
            dmax = MAX( rmax, dmax ) ; dmin = MIN( rmin, dmin )

            pmax = MAX( ABS( P( i ) ), ABS( P( i + 1 ) ) )
            pmin = MIN( ABS( P( i ) ), ABS( P( i + 1 ) ) )
            
            IF ( rmax >= big ) THEN
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE(  out, "( '2x2 block ', 2i7, ' eval = ', ES12.4 )" )     &
                 pmax - n, pmin - n, one / rmax
              END IF
            ELSE IF ( rmin == zero ) THEN
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE( out, "( '2x2 block ', 2i7, ' eval = infinity' )" )      &
                 pmax - n, pmin - n
              END IF
            END IF

            IF ( rmin >= big ) THEN
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE( out, "( '2x2 block ', 2i7, ' eval = ', ES12.4 )" )      &
                  pmin - n, pmax - n, one / rmin
              END IF
            ELSE IF ( rmax == zero ) THEN
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE( out, "( '2x2 block ', 2i7, ' eval = infinity' )" )      &
                  pmin - n, pmax - n 
              END IF
            END IF
            IF ( rmax < big .AND. rmin > zero )                                &
              dmax_allowed = MAX( rmax, dmax_allowed )

!  A 1x1 block

          ELSE
            rmax = ABS( D( 1, i ) ) ; rmin = rmax
            dmax = MAX( rmax, dmax ) ; dmin = MIN( rmin, dmin )
            IF ( rmax >= big ) THEN    
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE( out, "( '1x1 block ', i7, 7x, ' eval = ', ES12.4 )" )   &
                  P( i ) - n,  one / rmax
              END IF
            ELSE IF ( rmax == zero ) THEN
              n_depen = n_depen + 1
              IF ( out > 0 .AND. control%print_level >= 3 ) THEN
                WRITE( out, "( '1x1 block ', i7, 7x, ' eval = infinity' )" )   &
                  P( i ) - n
              END IF
            END IF
            IF ( rmax < big .AND. rmin > zero )                                &
              dmax_allowed = MAX( rmax, dmax_allowed )
          END IF
        ELSE

!  The final 1x1 block

          rmax = ABS( D( 1, i ) ) ; rmin = rmax
          dmax = MAX( rmax, dmax ) ; dmin = MIN( rmin, dmin )
          IF ( rmax >= big ) THEN    
            n_depen = n_depen + 1
            IF ( out > 0 .AND. control%print_level >= 3 ) THEN
              WRITE( out, "( '1x1 block ', i7, 7x, ' eval = ', ES12.4 )" )     &
                P( i ) - n,  one / rmax
            END IF
          ELSE IF ( rmax == zero ) THEN    
            n_depen = n_depen + 1
            IF ( out > 0 .AND. control%print_level >= 3 ) THEN
              WRITE( out, "( '1x1 block ', i7, 7x, ' eval = infinity ' )" )    &
                P( i ) - n
            END IF
          END IF
          IF ( rmax < big .AND. rmin > zero )                                  &
            dmax_allowed = MAX( rmax, dmax_allowed )
        END IF
      END DO

!  Any null blocks

      IF ( K%n > finfo%rank ) THEN
         n_depen = n_depen + K%n - finfo%rank
         dmax = HUGE( one )
      END IF

      DO i = finfo%rank + 1, K%n
        IF ( out > 0 .AND. control%print_level >= 3 )                          &
         WRITE( out, "( '1x1 block ', i7, 7x, ' eval = ', ES12.4 )" )          &
           P( i ) - n, zero
      END DO

      IF ( out > 0 .AND. control%print_level >= 1 ) THEN
        IF ( dmin == zero .OR. dmax == zero ) THEN
          WRITE( out, "( ' 1/ smallest,largest block eigenvalues =', 2ES12.4)")&
            dmin, dmax
        ELSE
          WRITE( out, "( ' smallest,largest block eigenvalues =',  2ES12.4 )") &
            one / dmax, one / dmin
          WRITE( out, "( ' smallest non-negligible eigenvalue =',  ES12.4 )")  &
              one / dmax_allowed
        END IF
        WRITE( out, "( I7, ' constraints appear to be dependent ' )" ) n_depen
      END IF

      IF ( dmax_allowed > zero )                                               &
        inform%non_negligible_pivot = one / dmax_allowed

!  Mark any dependent constraints for removal

      IF ( n_depen > 0 ) THEN

!  Allocate arrays to indicate which constraints have been freed

        array_name = 'fdc: data%C_depen'
        CALL SPACE_resize_array( n_depen, C_depen, inform%status,              &
               inform%alloc_status, array_name = array_name,                   &
               deallocate_error_fatal = control%deallocate_error_fatal,        &
               exact_size = control%space_critical,                            &
               bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= GALAHAD_ok ) RETURN

!  A second loop over the diagonal blocks
  
        n_depen = 0 ; twobytwo = .FALSE. 
        DO i = 1, finfo%rank
          IF ( twobytwo ) THEN
            twobytwo = .FALSE.
            CYCLE
          END IF
          IF ( i < finfo%rank ) THEN

!  A 2x2 block

            IF ( P( i ) < 0 ) THEN
              twobytwo = .TRUE.
  
              CALL ROOTS_quadratic( D( 1, i ) * D( 1, i + 1 ) - D( 2, i ) ** 2,& 
                - D( 1, i ) - D( 1, i + 1 ), one, epsmch, nroots, root1, root2 ) 

              IF ( ABS( root2 ) >= big .OR.                                    &
                   root1 == zero .OR. root2 == zero ) THEN
                IF ( ABS( root1 ) >= big .OR.                                  &
                     ( root1 == zero .AND. root2 == zero ) ) THEN
                  n_depen = n_depen + 1
                  C_depen( n_depen )                                           &
                    = MIN( ABS( P( i ) ), ABS( P( i + 1 ) ) ) - n
                  D( 1, i ) = zero ;  D( 2, i ) = zero ; D( 1, i + 1 ) = zero
                END IF
                n_depen = n_depen + 1
                C_depen( n_depen ) = MAX( ABS( P( i ) ), ABS( P( i + 1 ) ) ) - n
              END IF

            ELSE

!  A 1x1 block

              IF ( ABS( D( 1, i ) ) >= big .OR. D( 1, i ) == zero ) THEN    
                n_depen = n_depen + 1
                C_depen( n_depen ) = P( i ) - n
                D( 1, i ) = zero
              END IF
            END IF
          ELSE

!  The final 1x1 block

            IF ( ABS( D( 1, i ) ) >= big .OR. D( 1, i ) == zero ) THEN
              n_depen = n_depen + 1
              C_depen( n_depen ) = P( i ) - n
              D( 1, i ) = zero
            END IF
          END IF
        END DO

!  Any null blocks

        DO i = finfo%rank + 1, K%n
          n_depen = n_depen + 1
          C_depen( n_depen ) = P( i ) - n
        END DO

!  Reset "small" pivots to zero

        CALL SILS_alter_d( FACTORS, D, SINFO%flag )

!  Check to see if the constraints are consistent

        SOL( : n ) = zero
        SOL( n + 1 : K%n ) = C( : m ) / SCALE( : m )
        CALL SILS_solve( K, FACTORS, SOL, CNTL, SINFO )

        res_max = zero
        DO i = 1, m
          res = - C( i ) / SCALE( i )
          DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
            res = res + A_val( l ) * SOL( A_col( l ) ) / SCALE( i )
          END DO
          res_max = MAX( res_max, ABS( res ) )
        END DO

        IF ( res_max <= control%max_infeas ) THEN
          IF ( out > 0 .AND. control%print_level >= 1 ) WRITE( out,            &
            "( ' constraints are consistent: maximum infeasibility = ',        &
          &    ES12.4 )" ) res_max
        ELSE
          IF ( out > 0 .AND. control%print_level >= 1 ) WRITE( out,            &
            "( ' constraints are inconsistent: maximum infeasibility = ',      &
          &    ES12.4, /, 31X, ' is larger than control%max_infeas ' )" ) res_max
          inform%status = GALAHAD_error_primal_infeasible ; GO TO 800
        END IF

      END IF

!  Reset the pivot tolerance and zero pivoting threshold to their initial values

      inform%status = GALAHAD_ok

  800 CONTINUE
      CALL CPU_TIME( time_now ) ; inform%time%total = time_now - time_start
      IF ( out > 0 .AND. control%print_level >= 1 ) WRITE( out,                &
         "( 6( ' -' ), ' end of test for rank defficiency', 5( ' - ' ) )" )

      RETURN  
 
!  End of FDC_find_dependent

      END SUBROUTINE FDC_find_dependent

!-*-*-*-*-*-*-   F D C _ T E R M I N A T E   S U B R O U T I N E   -*-*-*-*-*

      SUBROUTINE FDC_terminate( K, FACTORS, C_depen, control, inform )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

!  Arguments:
!
!   K       see Subroutine FDC_find_dependent
!   FACTORS see Subroutine FDC_find_dependent
!   control see Subroutine FDC_initialize
!   inform  see Subroutine FDC_find_dependent

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_depen
      TYPE ( SMT_type ), INTENT( INOUT ) :: K
      TYPE ( SILS_factors ), INTENT( INOUT ) :: FACTORS
      TYPE ( FDC_control_type ), INTENT( IN ) :: control        
      TYPE ( FDC_inform_type ), INTENT( INOUT ) :: inform
 
!  Local variables

      CHARACTER ( LEN = 80 ) :: array_name

!  Deallocate all arrays allocated within SILS

      CALL SILS_finalize( FACTORS, control%CNTL, inform%alloc_status )
      IF ( inform%alloc_status /= GALAHAD_ok ) THEN
        inform%status = GALAHAD_error_deallocate
        inform%bad_alloc = ''
        IF ( control%deallocate_error_fatal ) RETURN
      END IF

!  Deallocate all remaining allocated arrays

      array_name = 'fdc: C_depen'
      CALL SPACE_dealloc_array( C_depen,                                       &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'fdc: K%row'
      CALL SPACE_dealloc_array( K%row,                                         &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'fdc: K%col'
      CALL SPACE_dealloc_array( K%col,                                         &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      array_name = 'fdc: K%val'
      CALL SPACE_dealloc_array( K%val,                                         &
         inform%status, inform%alloc_status, array_name = array_name,          &
         bad_alloc = inform%bad_alloc, out = control%error )
      IF ( control%deallocate_error_fatal .AND.                                &
           inform%status /= GALAHAD_ok ) RETURN

      RETURN

!  End of subroutine FDC_terminate

      END SUBROUTINE FDC_terminate 

!  End of module FDC

   END MODULE GALAHAD_FDC_double
