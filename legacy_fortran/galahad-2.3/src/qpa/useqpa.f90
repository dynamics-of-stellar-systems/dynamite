! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-  G A L A H A D   U S E Q P A    M O D U L E  -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authors: Nick Gould and Dominique Orban

!  History -
!   originally released pre GALAHAD Version 1.0. March 14th 2003
!   update released with GALAHAD Version 2.0. August 11th 2005

!  For full documentation, see
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_USEQPA_double

!    -----------------------------------------------
!    | CUTEr/AMPL interface to  QPA, a working-set |
!    ! algorithm for quadratic programming         |
!    -----------------------------------------------

      USE CUTEr_interface_double
!NOT95USE GALAHAD_CPU_time
      USE GALAHAD_SYMBOLS
      USE GALAHAD_QPT_double
      USE GALAHAD_QPA_double
      USE GALAHAD_SORT_double, only: SORT_reorder_by_rows
      USE GALAHAD_RAND_double
      USE GALAHAD_PRESOLVE_double
      USE GALAHAD_SPECFILE_double 
      USE GALAHAD_COPYRIGHT
      USE SCALING
      USE GALAHAD_SYMBOLS,                                                     &
          ACTIVE                => GALAHAD_ACTIVE,                             &
          TRACE                 => GALAHAD_TRACE,                              &
          DEBUG                 => GALAHAD_DEBUG,                              &
          GENERAL               => GALAHAD_GENERAL,                            &
          ALL_ZEROS             => GALAHAD_ALL_ZEROS

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: USE_QPA

    CONTAINS

!-*-*-*-*-*-*-*-*-*-   U S E _ Q P A  S U B R O U T I N E   -*-*-*-*-*-*-*-*-

     SUBROUTINE USE_QPA( input )

!  --------------------------------------------------------------------
!
!      Solve the l_1 quadratic program from CUTEr              
!                                                    
!        minimize     1/2 x(T) H x + g(T) x + f
!                       + rho_g min( A x - c_l , 0 )
!                       + rho_g max( A x - c_u , 0 )
!                       + rho_b min( x - x_l , 0 )
!                       + rho_b max( x - x_u , 0 )
!
!      using the GALAHAD package GALAHAD_QPA
!
!  --------------------------------------------------------------------

!  Dummy argument

     INTEGER, INTENT( IN ) :: input

!  Parameters

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
     REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
     REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
     REAL ( KIND = wp ), PARAMETER :: ten = 10.0_wp
     REAL ( KIND = wp ), PARAMETER :: stopr = ten ** ( - 10 )
     REAL ( KIND = wp ), PARAMETER :: infinity = ten ** 19
     REAL ( KIND = wp ), PARAMETER :: biginf = HUGE( one )

!    INTEGER, PARAMETER :: n_k = 100, k_k = 3, in = 28
!    REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( :, : ) :: k_val
!    CHARACTER ( len = 10 ) :: filename = 'k.val'

!  Scalars

     INTEGER :: n, m, ir, ic, ifail, nmax, mmax, mpn, la, lh, liw, iores
     INTEGER :: i, j, l, neh, nea, factorization_integer, factorization_real
     INTEGER :: status, mfixed, mdegen, iter, nfacts, nfixed, ndegen, mequal
     INTEGER :: alloc_stat, A_ne, H_ne, smt_stat
     INTEGER :: n_o, m_o, a_ne_o, h_ne_o
     INTEGER :: m_ref = 1000
     REAL :: time, timeo, times, timet, timep1, timep2, timep3, timep4
     REAL ( KIND = wp ) :: obj, qfval
     LOGICAL :: filexx, printo, printe, is_specfile
     CHARACTER ( LEN =  1 ) :: p_degen, d_degen
      
!  Specfile characteristics

     INTEGER, PARAMETER :: input_specfile = 34
     INTEGER, PARAMETER :: lspec = 21
     CHARACTER ( LEN = 16 ) :: specname = 'RUNQPA'
     TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec
     CHARACTER ( LEN = 16 ) :: runspec = 'RUNQPA.SPC'

!  The default values for QPA could have been set as:

! BEGIN RUNQPA SPECIFICATIONS (DEFAULT)
!  write-problem-data                        NO
!  problem-data-file-name                    QPA.data
!  problem-data-file-device                  26
!  write-initial-sif                         NO
!  initial-sif-file-name                     INITIAL.SIF
!  initial-sif-file-device                   51
!  initial-rho-g                             -1.0
!  initial-rho-b                             -1.0
!  scale-problem                             0
!  pre-solve-problem                         YES
!  write-presolved-sif                       NO
!  presolved-sif-file-name                   PRESOLVE.SIF
!  presolved-sif-file-device                 52
!  solve-problem                             YES
!  print-full-solution                       NO
!  write-solution                            NO
!  solution-file-name                        QPASOL.d
!  solution-file-device                      62
!  write-result-summary                      NO
!  result-summary-file-name                  QPARES.d
!  result-summary-file-device                47
! END RUNQPA SPECIFICATIONS

!  Default values for specfile-defined parameters

     INTEGER :: scale = 0
     INTEGER :: dfiledevice = 26
     INTEGER :: ifiledevice = 51
     INTEGER :: pfiledevice = 53
     INTEGER :: rfiledevice = 47
     INTEGER :: sfiledevice = 62
     LOGICAL :: write_problem_data   = .FALSE.
     LOGICAL :: write_initial_sif    = .FALSE.
     LOGICAL :: write_presolved_sif  = .FALSE.
     LOGICAL :: write_solution       = .FALSE.
     LOGICAL :: write_result_summary = .FALSE.
     CHARACTER ( LEN = 30 ) :: dfilename = 'QPA.data'
     CHARACTER ( LEN = 30 ) :: ifilename = 'INITIAL.SIF'
     CHARACTER ( LEN = 30 ) :: pfilename = 'PRESOLVE.SIF'
     CHARACTER ( LEN = 30 ) :: rfilename = 'QPARES.d'
     CHARACTER ( LEN = 30 ) :: sfilename = 'QPASOL.d'
     LOGICAL :: do_presolve = .TRUE.
     LOGICAL :: do_solve = .TRUE.
     LOGICAL :: fulsol = .FALSE. 

!  Output file characteristics

     INTEGER, PARAMETER :: out  = 6
     INTEGER :: errout = 6
     CHARACTER ( LEN =  5 ) :: state, solv = ' QPA '
     CHARACTER ( LEN = 10 ) :: pname

!  Arrays

     TYPE ( QPA_data_type ) :: QPA_data
     TYPE ( QPA_control_type ) :: QPA_control        
     TYPE ( QPA_inform_type ) :: QPA_inform
     TYPE ( QPT_problem_type ) :: prob
     TYPE ( SCALING_control_type ) :: control
     TYPE ( RAND_seed ) :: seed
     TYPE ( PRESOLVE_control_type ) :: PRE_control
     TYPE ( PRESOLVE_inform_type )  :: PRE_inform
     TYPE ( PRESOLVE_data_type ) :: PRE_data

!  Allocatable arrays

     CHARACTER ( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: VNAME, CNAME
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: SH, SA, X0, C
     LOGICAL, ALLOCATABLE, DIMENSION( : ) :: EQUATN, LINEAR
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW, C_stat, B_stat

     CALL CPU_TIME( time )

!  Determine the number of variables and constraints

     CALL CDIMEN( input, nmax, mmax )
     nmax = MAX( nmax, 1 ) ; mmax = MAX( mmax, 1 ) !  for fortran 77 CUTE tools

!  Allocate suitable arrays

     ALLOCATE( prob%X( nmax ), prob%X_l( nmax ), prob%X_u( nmax ),            &
               prob%G( nmax ), VNAME( nmax ), B_stat( nmax ),                 &
               STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'X', alloc_stat ; STOP
     END IF

     ALLOCATE( prob%C_l( mmax ), prob%C_u( mmax ), prob%Y( mmax ),            &
               CNAME( mmax ), EQUATN( mmax ), LINEAR( mmax ),                 &
               C_stat( mmax ), STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'C', alloc_stat ; STOP
     END IF

!  Set up the data structures necessary to hold the group partially
!  separable function.

     CALL CSETUP( input, out, n, m, prob%X, prob%X_l, prob%X_u, nmax, EQUATN, &
                  LINEAR, prob%Y, prob%C_l, prob%C_u, mmax, .FALSE.,          &
                  .FALSE., .FALSE. )
     DEALLOCATE( LINEAR )

!  Allocate derived types

     ALLOCATE( X0( n ), STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'X0', alloc_stat
       STOP
     END IF

     ALLOCATE( prob%C( m ), C( m ), STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'C', alloc_stat
       STOP
     END IF

!  Allocate and initialize dual variables.

     ALLOCATE( prob%Z( n ), STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'Z', alloc_stat
       STOP
     END IF
     prob%Z( : n ) = one

!  Determine the names of the problem, variables and constraints.

     CALL CNAMES( n, m, pname, VNAME, CNAME )
     WRITE( out, 2020 ) pname 
     WRITE( out, "( /, ' n   = ', i8, ' nmax = ', i8,                         &
                     & ' m   = ', i8, ' mmax = ', i8 )" ) n, nmax, m, mmax

!  Set up the initial estimate of the solution and
!  right-hand-side of the Kuhn-Tucker system.

!  Determine the constant terms for the problem functions.

     prob%X( : n ) = MIN( prob%X_u( : n ),                                    &
                          MAX( prob%X_l( : n ), prob%X( : n ) ) )

!  Set X0 to zero to determine the constant terms for the problem functions

     X0 = zero 

!  Evaluate the constant terms of the objective (obj) and constraint 
!  functions (C)

     CALL CFN( n, m, X0, obj, m, C( : m ) )
     DO i = 1, m 
       IF ( EQUATN( i ) ) THEN 
         prob%C_l( i ) = prob%C_l( i ) - C( i )
         prob%C_u( i ) = prob%C_l( i )
       ELSE
         prob%C_l( i ) = prob%C_l( i ) - C( i )
         prob%C_u( i ) = prob%C_u( i ) - C( i )
       END IF 
     END DO

!  Determine the number of nonzeros in the Jacobian

     CALL CDIMSJ( la )
     la = MAX( la, 1 )

!  Allocate arrays to hold the Jacobian

     ALLOCATE( prob%A%row( la ), prob%A%col( la ), prob%A%val( la ),          &
               STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, 2150 ) 'A', alloc_stat ; STOP
     END IF

!  Evaluate the linear terms of the constraint functions

     CALL CSGR( n, m, .FALSE., mmax, prob%Y, X0, nea, la, prob%A%val,         &
                prob%A%col, prob%A%row )

     DEALLOCATE( X0 )
      
!  Exclude zeros; set the linear term for the objective function

     A_ne = 0
     prob%G( : n )      = zero
     prob%gradient_kind = ALL_ZEROS
     DO i = 1, nea
       IF ( prob%A%val( i ) /= zero ) THEN
         IF ( prob%A%row( i ) > 0 ) THEN
           A_ne = A_ne + 1
           prob%A%row( A_ne ) = prob%A%row( i ) 
           prob%A%col( A_ne ) = prob%A%col( i )
           prob%A%val( A_ne ) = prob%A%val( i )
         ELSE
           prob%G( prob%A%col( i ) ) = prob%A%val( i )
           prob%gradient_kind        = GENERAL
         END IF  
       END IF
     END DO

!  Determine the number of nonzeros in the Hessian

     CALL CDIMSH( lh ) ; lh = MAX( lh, 1 )

!  Allocate arrays to hold the Hessian

     ALLOCATE( prob%H%row( lh ), prob%H%col( lh ), prob%H%val( lh ),          &
               STAT = alloc_stat )
     IF ( alloc_stat /= 0 ) THEN
       WRITE( out, "( ' nea = ', i8, ' la   = ', i8 )" ) nea, la
       WRITE( out, 2150 ) 'H', alloc_stat
       STOP
     END IF

!  Evaluate the Hessian of the Lagrangian function at the initial point.

     CALL CSH( n, m, prob%X, mmax, prob%Y, neh, lh, prob%H%val, prob%H%row,   &
               prob%H%col )
     WRITE( out, "( ' nea = ', i8, ' la   = ', i8,                            &
    &               ' neh  = ', i8, ' lh   = ', i8 )" ) nea, la, neh, lh

!  Remove Hessian out of range

     H_ne = 0
     DO l = 1, neh    
       i = prob%H%row( l ) ; j = prob%H%col( l )
       IF ( i < 1 .OR. i > n .OR. j < 1 .OR. j > n ) CYCLE
       H_ne = H_ne + 1 ; prob%H%val( H_ne ) = prob%H%val( l )
       IF ( i >= j ) THEN
         prob%H%row( H_ne ) = i
         prob%H%col( H_ne ) = j
       ELSE
         prob%H%row( H_ne ) = j
         prob%H%col( H_ne ) = i
       END IF
     END DO

     liw = MAX( m, n ) + 1
     ALLOCATE( prob%A%ptr( m + 1 ), prob%H%ptr( n + 1 ) )
     ALLOCATE( IW( liw ) )

!    ALLOCATE( k_val( n_k, n_k ) )
!    OPEN( in, FILE = filename, FORM = 'FORMATTED', STATUS = 'OLD' ) 
!    REWIND in
!    DO j = 1, n_k
!      DO i = 1, n_k
!         READ( in, "( ES24.16 )" ) k_val( i, j )
!      END DO
!    END DO
!    CLOSE( in )
!    DO l = 1, H_ne
!      i = MOD( prob%H%row( l ), n_k ) ; IF ( i == 0 ) i = n_k
!      j = MOD( prob%H%col( l ), n_k ) ; IF ( j == 0 ) j = n_k
!      IF ( prob%H%row( l ) <= k_k * n_k .AND.                                &
!           prob%H%col( l ) <= k_k * n_k ) THEN
!        IF ( ABS( prob%H%val( l ) - k_val( i, j ) ) > 0.0001 )               &
!          WRITE( 6, "( 2I6, 2ES22.14 )" )                                    &
!            prob%H%row( l ), prob%H%col( l ), prob%H%val( l ),  k_val( i, j )
!         prob%H%val( l ) = k_val( i, j )
!      ELSE
!        IF ( ABS( prob%H%val( l ) + k_val( i, j ) / k_k ) > 0.0001 )         &
!          WRITE( 6, "( 2I6, 2ES22.14 )" ) prob%H%row( l ), prob%H%col( l ),  &
!            prob%H%val( l ), - k_val( i, j ) / k_k
!         prob%H%val( l ) = - k_val( i, j ) / k_k
!      END IF
!    END DO
!    DEALLOCATE( k_val )

!  Transform A to row storage format

     IF ( A_ne /= 0 ) THEN
       CALL SORT_reorder_by_rows( m, n, A_ne, prob%A%row, prob%A%col, A_ne,   &
                                  prob%A%val, prob%A%ptr, m + 1, IW, liw,     &
                                  out, out, i )
     ELSE
       prob%A%ptr = 0
     END IF

!  Same for H

     IF ( H_ne /= 0 ) THEN
       CALL SORT_reorder_by_rows( n, n, H_ne, prob%H%row, prob%H%col, H_ne,   &
                                  prob%H%val, prob%H%ptr, n + 1, IW, liw,     &
                                  out, out, i )
     ELSE
       prob%H%ptr = 0
     END IF

!  Deallocate arrays holding matrix row indices

     DEALLOCATE( prob%A%row, prob%H%row )
     DEALLOCATE( IW )
     ALLOCATE( prob%A%row( 0 ), prob%H%row( 0 ) )

      prob%new_problem_structure = .TRUE.

!  Store the problem dimensions

     prob%n = n
     prob%m = m
     IF ( ALLOCATED( prob%H%type ) ) DEALLOCATE( prob%H%type )
     CALL SMT_put( prob%H%type, 'SPARSE_BY_ROWS', smt_stat )
     IF ( ALLOCATED( prob%A%type ) ) DEALLOCATE( prob%A%type )
     CALL SMT_put( prob%A%type, 'SPARSE_BY_ROWS', smt_stat )
     prob%f = obj ; prob%rho_g = 2 * m ; prob%rho_b = 2 * n

!  ------------------- problem set-up complete ----------------------

     CALL CPU_TIME( times )

!  ------------------ Open the specfile for runqpa ----------------

     INQUIRE( FILE = runspec, EXIST = is_specfile )
     IF ( is_specfile ) THEN
       OPEN( input_specfile, FILE = runspec, FORM = 'FORMATTED', STATUS = 'OLD' )

!   Define the keywords

       spec( 1 )%keyword = 'write-problem-data'
       spec( 2 )%keyword = 'problem-data-file-name'
       spec( 3 )%keyword = 'problem-data-file-device'
       spec( 4 )%keyword = 'write-initial-sif'
       spec( 5 )%keyword = 'initial-sif-file-name'
       spec( 6 )%keyword = 'initial-sif-file-device'
       spec( 7 )%keyword = 'initial-rho-g'
       spec( 8 )%keyword = 'initial-rho-b'
       spec( 9 )%keyword = 'scale-problem'
       spec( 10 )%keyword = 'pre-solve-problem'
       spec( 11 )%keyword = 'write-presolved-sif'
       spec( 12 )%keyword = 'presolved-sif-file-name'
       spec( 13 )%keyword = 'presolved-sif-file-device'
       spec( 14 )%keyword = 'solve-problem'
       spec( 15 )%keyword = 'print-full-solution'
       spec( 16 )%keyword = 'write-solution'
       spec( 17 )%keyword = 'solution-file-name'
       spec( 18 )%keyword = 'solution-file-device'
       spec( 19 )%keyword = 'write-result-summary'
       spec( 20 )%keyword = 'result-summary-file-name'
       spec( 21 )%keyword = 'result-summary-file-device'

!   Read the specfile

       CALL SPECFILE_read( input_specfile, specname, spec, lspec, errout )

!   Interpret the result

       CALL SPECFILE_assign_logical( spec( 1 ), write_problem_data, errout )
       CALL SPECFILE_assign_string ( spec( 2 ), dfilename, errout )
       CALL SPECFILE_assign_integer( spec( 3 ), dfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 4 ), write_initial_sif, errout )
       CALL SPECFILE_assign_string ( spec( 5 ), ifilename, errout )
       CALL SPECFILE_assign_integer( spec( 6 ), ifiledevice, errout )
       CALL SPECFILE_assign_real( spec( 7 ), prob%rho_g, errout )
       CALL SPECFILE_assign_real( spec( 8 ), prob%rho_b, errout )
       CALL SPECFILE_assign_integer( spec( 9 ), scale, errout )
       CALL SPECFILE_assign_logical( spec( 10 ), do_presolve, errout )
       CALL SPECFILE_assign_logical( spec( 11 ), write_presolved_sif, errout )
       CALL SPECFILE_assign_string ( spec( 12 ), pfilename, errout )
       CALL SPECFILE_assign_integer( spec( 13 ), pfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 14 ), do_solve, errout )
       CALL SPECFILE_assign_logical( spec( 15 ), fulsol, errout )
       CALL SPECFILE_assign_logical( spec( 16 ), write_solution, errout )
       CALL SPECFILE_assign_string ( spec( 17 ), sfilename, errout )
       CALL SPECFILE_assign_integer( spec( 18 ), sfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 19 ), write_result_summary, errout )
       CALL SPECFILE_assign_string ( spec( 20 ), rfilename, errout )
       CALL SPECFILE_assign_integer( spec( 21 ), rfiledevice, errout )
     END IF

     IF ( prob%rho_g <= zero ) prob%rho_g = 2 * m 
     IF ( prob%rho_b <= zero ) prob%rho_b = 2 * n

!  If required, print out the (raw) problem data

     IF ( write_problem_data ) THEN
       INQUIRE( FILE = dfilename, EXIST = filexx )
       IF ( filexx ) THEN
          OPEN( dfiledevice, FILE = dfilename, FORM = 'FORMATTED',            &
                STATUS = 'OLD', IOSTAT = iores )
       ELSE
          OPEN( dfiledevice, FILE = dfilename, FORM = 'FORMATTED',            &
                 STATUS = 'NEW', IOSTAT = iores )
       END IF
       IF ( iores /= 0 ) THEN 
         write( out, 2160 ) iores, dfilename
         STOP
       END IF

       WRITE( dfiledevice, "( 'n, m = ', 2I6, ' obj = ', ES12.4 )" )          &
         n, m, prob%f
       WRITE( dfiledevice, "( ' g ', /, ( 5ES12.4 ) )" ) prob%G( : n )
       WRITE( dfiledevice, "( ' x_l ', /, ( 5ES12.4 ) )" ) prob%X_l( : n )
       WRITE( dfiledevice, "( ' x_u ', /, ( 5ES12.4 ) )" ) prob%X_u( : n )
       WRITE( dfiledevice, "( ' c_l ', /, ( 5ES12.4 ) )" ) prob%C_l( : m )
       WRITE( dfiledevice, "( ' c_u ', /, ( 5ES12.4 ) )" ) prob%C_u( : m )
       WRITE( dfiledevice, "( ' A_ptr ', /, ( 10I6 ) )" ) prob%A%ptr( : m + 1 )
       WRITE( dfiledevice, "( ' A_col ', /, ( 10I6 ) )" ) prob%A%col( : A_ne )
       WRITE( dfiledevice, "( ' A_val ', /, ( 5ES12.4 ) )" )                  &
         prob%A%val( : A_ne )
       WRITE( dfiledevice, "( ' H_ptr ', /, ( 10I6 ) )" ) prob%H%ptr( : n + 1 )
       WRITE( dfiledevice, "( ' H_col ', /, ( 10I6 ) )" ) prob%H%col( : H_ne )
       WRITE( dfiledevice, "( ' H_val ', /, ( 5ES12.4 ) )" )                  &
         prob%H%val( : H_ne )

       CLOSE( dfiledevice )
     END IF

!  If required, append results to a file

     IF ( write_result_summary ) THEN
       INQUIRE( FILE = rfilename, EXIST = filexx )
       IF ( filexx ) THEN
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',            &
                STATUS = 'OLD', POSITION = 'APPEND', IOSTAT = iores )
       ELSE
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',            &
                STATUS = 'NEW', IOSTAT = iores )
       END IF
       IF ( iores /= 0 ) THEN 
         write( out, 2160 ) iores, rfilename
         STOP
       END IF
       WRITE( rfiledevice, 2180 ) pname
     END IF

!  Update control parameters if required.

     CALL SCALING_initialize( control )
     CALL QPA_initialize( QPA_data, QPA_control )
     IF ( is_specfile )                                                       &
       CALL QPA_read_specfile( QPA_control, input_specfile )
     
     control%print_level = QPA_control%print_level
     control%out         = QPA_control%out
     control%out_error   = QPA_control%error

     printo = out > 0 .AND. QPA_control%print_level > 0
     printe = out > 0 .AND. QPA_control%print_level >= 0

     IF ( printo ) CALL COPYRIGHT( out, '2002' )

!  Initalize random number seed

     CALL RAND_initialize( seed )

!  Set all default values, and override defaults if requested
 
     IF ( QPA_control%cold_start == 0 ) THEN
       IF ( m > 0 ) THEN
         mpn = MIN( m + n,                                                    &
                    COUNT( prob%X_l > - biginf .OR. prob%X_l < biginf )  +    &
                    COUNT( prob%C_l > - biginf .OR. prob%C_l < biginf ) )
       ELSE
         mpn = MIN( n, COUNT( prob%X_l > - biginf .OR. prob%X_l < biginf ) )
       END IF
       CALL RAND_random_integer( seed, mpn + 1, m_ref )
       m_ref = m_ref - 1
       C_stat = 0 ; B_stat = 0
       DO i = 1, m_ref
         DO
           CALL RAND_random_integer( seed, mpn, j )
           IF ( j > m ) THEN
             IF ( B_stat( j - m ) == 0 ) EXIT
             IF ( prob%X_l( j - m ) <= - biginf .AND.                         &
                  prob%X_l( j - m ) >= biginf ) EXIT
           ELSE
             IF ( C_stat( j ) == 0 ) EXIT
             IF ( prob%C_l( j ) <= - biginf .AND.                             &
                  prob%C_l( j ) >= biginf ) EXIT
           END IF
         END DO
         IF ( j > m ) THEN
           IF ( prob%X_l( j - m ) > - biginf ) THEN
             B_stat( j - m ) = - 1
           ELSE IF ( prob%X_l( j - m ) < biginf ) THEN
             B_stat( j - m ) = 1
           END IF
         ELSE
           IF ( prob%C_l( j ) > - biginf ) THEN
             C_stat( j ) = - 1
           ELSE IF ( prob%C_l( j ) < biginf ) THEN
             C_stat( j ) = 1
           END IF
         END IF
       END DO
     END IF

!  If required, scale the problem

     IF ( scale > 0 ) THEN
       ALLOCATE( SH( n ), SA( m ), STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
         IF ( printe ) WRITE( out, 2150 ) 'SH/SA', alloc_stat ; STOP
       END IF

!  Scale using K

       IF ( scale == 1 .OR. scale == 4 ) THEN
         IF ( printo ) WRITE( out, 2140 ) 'K'
         CALL SCALING_get_factors_from_K( n, m, prob%H%val, prob%H%col,       &
                                          prob%H%ptr, prob%A%val, prob%A%col, &
                                          prob%A%ptr, SH, SA, control, ifail )

!  Scale using A

       ELSE IF ( scale == 2 .OR. scale == 5 ) THEN
         IF ( printo ) WRITE( out, 2140 ) 'A'
         CALL SCALING_get_factors_from_A( n, m, prob%A%val, prob%A%col,       &
                                          prob%A%ptr, SH, SA, control, ifail )
       ELSE IF ( scale == 3 ) THEN
         SH = one ; SA = one
       END IF

!  Rescale A

       IF ( scale >= 3 ) THEN
         IF ( printo ) WRITE( out, 2170 )
         CALL SCALING_normalize_rows_of_A( n, m, prob%A%val, prob%A%col,      &
                                           prob%A%ptr, SH, SA )
       END IF

!  Apply the scaling factors

       CALL SCALING_apply_factors( n, m, prob%H%val, prob%H%col, prob%H%ptr,  &
                                   prob%A%val, prob%A%col, prob%A%ptr,        &
                                   prob%G, prob%X, prob%X_l, prob%X_u,        &
                                   prob%C_l, prob%C_u, prob%Y, prob%Z,        &
                                   infinity, SH, SA, .TRUE. )
     END IF

!  If the preprocessor is to be used, or the problem to be output,
!  allocate sufficient space

     IF ( write_initial_sif .OR. do_presolve ) THEN

       ALLOCATE( prob%X_status( n ), STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
         IF ( printe ) WRITE( out, 2150 ) 'X_status', alloc_stat
         STOP
       END IF
       prob%X_status = ACTIVE
       
       ALLOCATE( prob%C_status( m ), STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
         IF ( printe ) WRITE( out, 2150 ) 'X_status', alloc_stat ; STOP
       END IF
       prob%C_status =  ACTIVE
       
       ALLOCATE( prob%Z_l( n ), prob%Z_u( n ), STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
         IF ( printe ) WRITE( out, 2150 ) 'Z_lu', alloc_stat ; STOP
       END IF
       prob%Z_l( : n ) = - infinity
       prob%Z_u( : n ) =   infinity
       
       ALLOCATE( prob%Y_l( m ), prob%Y_u( m ), STAT = alloc_stat )
       IF ( alloc_stat /= 0 ) THEN
         IF ( printe ) WRITE( out, 2150 ) 'C_lu', alloc_stat ; STOP
       END IF
       prob%Y_l( : m ) = - infinity
       prob%Y_u( : m ) =   infinity

!  Writes the initial SIF file, if needed

       IF ( write_initial_sif ) THEN
         CALL QPT_write_to_sif( prob, pname, ifilename, ifiledevice,          &
                                .FALSE., .FALSE., infinity )
         IF ( .NOT. ( do_presolve .OR. do_solve ) ) STOP
       END IF

     END IF

!  Presolve

     IF ( do_presolve ) THEN

!      write( out, * ) ' >>> ===== starting presolve ===='
       
!      Allocate and initialize status variables.
       
       h_ne_o      = H_ne
       a_ne_o      = A_ne
       n_o         = prob%n
       m_o         = prob%m
       
!      prob%fgh_status = ACTIVE
       
       CALL CPU_TIME( timep1 )
       
!      set control variables
       
       CALL PRESOLVE_initialize( PRE_control, PRE_inform, PRE_data )
       IF ( is_specfile )                                                      &
         CALL PRESOLVE_read_specfile( input_specfile, PRE_control, PRE_inform )

       IF ( PRE_inform%status /= 0 ) STOP

       PRE_control%infinity   = QPA_control%infinity
       PRE_control%c_accuracy = ten * QPA_control%feas_tol
       PRE_control%z_accuracy = ten * QPA_control%feas_tol

!  Call the presolver

       CALL PRESOLVE_apply( prob, PRE_control, PRE_inform, PRE_data )
       IF ( PRE_inform%status < 0 ) THEN
         WRITE( out, "( '  ERROR return from PRESOLVE (exitc =', I6, ')' )" ) &
           PRE_inform%status
         STOP
       END IF
       
       CALL CPU_TIME( timep2 )
       
       A_ne = prob%A%ptr( prob%m + 1 ) - 1
       H_ne = prob%H%ptr( prob%n + 1 ) - 1
       IF ( printo ) WRITE( out, 2200 ) n_o, m_o, a_ne_o, h_ne_o,prob%n,      &
         prob%m, MAX( 0, A_ne ), MAX( 0, H_ne ), timep2 - timep1,             &
         PRE_inform%nbr_transforms

!  If required, write a SIF file containing the presolved problem

       IF ( write_presolved_sif ) THEN
         CALL QPT_write_to_sif( prob, pname, pfilename, pfiledevice,          &
                                .FALSE., .FALSE., QPA_control%infinity )
       END IF
     END IF

!  Solve the problem
!  =================

     IF ( do_solve ) THEN
       CALL CPU_TIME( timeo )
       CALL QPA_solve( prob, C_stat, B_stat, QPA_data, QPA_control, QPA_inform)
       CALL CPU_TIME( timet )

       status = QPA_inform%status ; iter = QPA_inform%iter
       nfacts = QPA_inform%nfacts
       factorization_integer = QPA_inform%factorization_integer 
       factorization_real = QPA_inform%factorization_real
       CALL QPA_terminate( QPA_data, QPA_control, QPA_inform )
       qfval = QPA_inform%obj
     ELSE
       timeo  = 0.0
       timet  = 0.0
       iter   = 0
       solv   = ' NONE'
       status = 0
       nfacts = 0
       factorization_integer = 0
       factorization_real    = 0.0
       qfval = prob%f
     END IF

!  Restore from presolve

     IF ( do_presolve ) THEN
       IF ( PRE_control%print_level >= DEBUG )                                &
         CALL QPT_write_problem( out, prob )
       
       CALL CPU_TIME( timep3 )
       CALL PRESOLVE_restore( prob, PRE_control, PRE_inform, PRE_data )
       IF ( PRE_inform%status /= 0 .AND. printo )                             &
         WRITE( out, " ( /, ' Warning: info%status following',                &
      &  ' PRESOLVE_restore is ', I5, / ) " ) PRE_inform%status
!      IF ( PRE_inform%status /= 0 ) STOP
       CALL PRESOLVE_terminate( PRE_control, PRE_inform, PRE_data )
       IF ( PRE_inform%status /= 0 .AND. printo )                             &
         WRITE( out, " ( /, ' Warning: info%status following',                &
      &    ' PRESOLVE_terminate is ', I5, / ) " ) PRE_inform%status
!      IF ( PRE_inform%status /= 0 ) STOP
       IF ( .NOT. do_solve ) STOP
       CALL CPU_TIME( timep4 )
       IF ( printo ) WRITE( out, 2210 )                                       &
         timep4 - timep3, timep2 - timep1 + timep4 - timep3
       qfval = prob%q
     ELSE
       PRE_control%print_level = TRACE
     END IF

!  If the problem was scaled, unscale it.

     IF ( scale > 0 ) THEN
       CALL SCALING_apply_factors( n, m, prob%H%val, prob%H%col, prob%H%ptr,  &
                                   prob%A%val, prob%A%col, prob%A%ptr,        &
                                   prob%G, prob%X, prob%X_l, prob%X_u,        &
                                   prob%C_l, prob%C_u, prob%Y, prob%Z,        &
                                   infinity, SH, SA, .FALSE., C = prob%C )
       DEALLOCATE( SH, SA )
     END IF

!  Print details of the solution obtained

     WRITE( out, 2010 ) status
     ndegen = 0  ; mdegen = 0 
     IF ( status == 0 .OR. status == - 8 ) THEN
       l = 4 ; IF ( fulsol ) l = n 

!  Print details of the primal and dual variables

       WRITE( out, 2090 ) 
       DO j = 1, 2 
         IF ( j == 1 ) THEN 
           ir = 1 ; ic = MIN( l, n ) 
         ELSE 
           IF ( ic < n - l ) WRITE( out, 2000 ) 
           ir = MAX( ic + 1, n - ic + 1 ) ; ic = n 
         END IF 
         DO i = ir, ic 
           state = ' FREE' 
           IF ( ABS( prob%X  ( i ) - prob%X_l( i ) ) < ten * stopr )          &
             state = 'LOWER'
           IF ( ABS( prob%X  ( i ) - prob%X_u( i ) ) < ten * stopr )          &
             state = 'UPPER'
           IF ( ABS( prob%X_l( i ) - prob%X_u( i ) ) <     1.0D-10 )          &
             state = 'FIXED'
           WRITE( out, 2050 ) i, VNAME( i ), state, prob%X( i ),              &
                              prob%X_l( i ), prob%X_u( i ), prob%Z( i )
!                             prob%X_l( i ), prob%X_u( i ), zero
         END DO 
       END DO 

!  Compute the number of fixed and degenerate variables.

       nfixed = 0 
       DO i = 1, n 
         IF ( ABS( prob%X( i ) - prob%X_l( i ) ) < stopr ) THEN
           nfixed = nfixed + 1 
!          IF ( ABS( prob%Z( i ) ) < ten * stopr ) ndegen = ndegen + 1 
         ELSE IF ( ABS( prob%X( i ) - prob%X_u( i ) ) < ten * stopr ) THEN
           nfixed = nfixed + 1 
!          IF ( ABS( prob%Z( i ) ) < ten * stopr ) ndegen = ndegen + 1 
         END IF 
       END DO 

!  Print details of the constraints.

       IF ( m > 0 ) THEN 

         WRITE( out, 2040 ) 
         l = 2  ; IF ( fulsol ) l = m 
         DO j = 1, 2 
           IF ( j == 1 ) THEN 
             ir = 1 ; ic = MIN( l, m ) 
           ELSE 
             IF ( ic < m - l ) WRITE( out, 2000 ) 
             ir = MAX( ic + 1, m - ic + 1 ) ; ic = m 
           END IF 
           DO i = ir, ic 
             state = ' FREE' 
             IF ( ABS( prob%C( I )   - prob%C_l( i ) ) < ten * stopr )        &
               state = 'LOWER' 
             IF ( ABS( prob%C( I )   - prob%C_u( i ) ) < ten * stopr )        &
               state = 'UPPER' 
             IF ( ABS( prob%C_l( i ) - prob%C_u( i ) ) <       stopr )        &
               state = 'EQUAL' 
             WRITE( out, 2130 ) i, CNAME( i ), STATE, prob%C( i ),            &
                                prob%C_l( i ), prob%C_u( i ), prob%Y( i ) 
           END DO 
         END DO 

!  Compute the number of equality, fixed inequality and degenerate constraints

         mequal = 0 ; mfixed = 0
         DO i = 1, m 
           IF ( ABS( prob%C( i ) - prob%C_l( i ) ) < ten * stopr .OR.         &
                ABS( prob%C( i ) - prob%C_u( i ) ) < ten * stopr ) THEN
             IF ( ABS( prob%C_l( i ) - prob%C_u( i ) ) < ten * stopr ) THEN 
                mequal = mequal + 1 
             ELSE 
                mfixed = mfixed + 1 
             END IF 
             IF ( ABS( prob%Y( i ) ) < stopr ) mdegen = mdegen + 1 
           END IF 
         END DO 
       END IF 
       WRITE( out, 2100 ) n, nfixed, ndegen 
       IF ( m > 0 ) THEN 
          WRITE( out, 2110 ) m, mequal, mdegen 
          IF ( m /= mequal ) WRITE( out, 2120 ) mfixed 
       END IF 
       WRITE( out, 2030 ) qfval, iter, nfacts,                                &
              factorization_integer, factorization_real 

!  If required, write the solution to a file

       IF ( write_solution ) THEN
         INQUIRE( FILE = sfilename, EXIST = filexx )
         IF ( filexx ) THEN
            OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',          &
                STATUS = 'OLD', IOSTAT = iores )
         ELSE
            OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',          &
                 STATUS = 'NEW', IOSTAT = iores )
         END IF
         IF ( iores /= 0 ) THEN 
           write( out, 2160 ) iores, sfilename
           STOP
         END IF

         WRITE( sfiledevice, 2250 ) pname, solv, qfval
         WRITE( sfiledevice, 2090 ) 

         DO i = 1, n 
           state = ' FREE' 
           IF ( ABS( prob%X( i )   - prob%X_l( i ) ) < ten * stopr )          &
             state = 'LOWER' 
           IF ( ABS( prob%X( i )   - prob%X_u( i ) ) < ten * stopr )          &
             state = 'UPPER' 
           IF ( ABS( prob%X_l( I ) - prob%X_u( I ) ) < stopr )                &
             state = 'FIXED' 
           WRITE( sfiledevice, 2050 ) i, VNAME( i ), STATE, prob%X( i ),      &
             prob%X_l( i ), prob%X_u( i ), prob%Z( i )
         END DO 
  
         IF ( m > 0 ) THEN 
           WRITE( sfiledevice, 2040 ) 
           DO i = 1, m 
             state = ' FREE' 
             IF ( ABS( prob%C( I ) - prob%C_l( i ) ) < ten * stopr )          &
               state = 'LOWER'
             IF ( ABS( prob%C( I ) - prob%C_u( i ) ) < ten * stopr )          &
               state = 'UPPER'
             IF ( ABS( prob%C_l( i ) - prob%C_u( i ) ) < stopr )              &
               state = 'EQUAL' 
             WRITE( sfiledevice, 2130 ) i, CNAME( i ), STATE, prob%C( i ),    &
               prob%C_l( i ), prob%C_u( i ), prob%Y( i )   
           END DO 
         END IF 
  
         WRITE( sfiledevice, 2030 ) qfval, iter, nfacts,                      &
                factorization_integer, factorization_real 
         CLOSE( sfiledevice ) 
       END IF 
     END IF 

     times = times - time ; timet = timet - timeo
     WRITE( out, 2060 ) times + timet 
     WRITE( out, 2070 ) pname 

!  Compare the variants used so far

     WRITE( out, 2080 ) solv, iter, nfacts, qfval, status, times,             &
                        timet, times + timet 

     IF ( write_result_summary ) THEN
       IF ( mdegen == 0 ) THEN ; p_degen = ' ' ; ELSE ; p_degen = 'P' ; END IF
       IF ( ndegen == 0 ) THEN ; d_degen = ' ' ; ELSE ; d_degen = 'D' ; END IF
       BACKSPACE( rfiledevice )
       WRITE( rfiledevice, 2190 ) pname, n, m, iter, nfacts, QPA_inform%obj,  &
         status, timet, p_degen, d_degen
     END IF

     DEALLOCATE( VNAME, CNAME, C )
     IF ( is_specfile ) CLOSE( input_specfile )

     RETURN

!  Non-executable statements

 2000 FORMAT( '      . .          .....  ..........',                          &
              '  ..........  ..........  .......... ' ) 
 2010 FORMAT( /,' Stopping with inform%status = ', I3 ) 
 2020 FORMAT( /, ' Problem: ', A10 )
 2030 FORMAT( /,' Final objective function value ', ES22.14, /,                &
          ' Total number of iterations = ',I6,' Number of factorizations = ',  &
          I6, //, I10, ' integer and ', I10, ' real words required',           &
          ' for the factorization' ) 
 2040 FORMAT( /,' Constraints : ', /, '                              ',        &
                '        <------ Bounds ------> ', /                           &
                '      # name       state    value   ',                        &
                '    Lower       Upper     Multiplier ' ) 
 2050 FORMAT( I7, 1X, A10, A6, 4ES12.4 ) 
 2060 FORMAT( /, ' Total time = ', 0P, F12.2 ) 
 2070 FORMAT( /, ' Problem: ', A10, //,                                        &
              '                                 objective',                    &
              '          < ------ time ----- > ', /,                           &
              ' Method  iterations   factors      value  ',                    &
              '   status setup   solve   total', /,                            &
              ' ------  ----------   -------    ---------',                    &
              '   ------ -----   -----   -----  ' ) 
 2080 FORMAT( A5, 2I10, 6X, ES12.4, I6, 0P, 3F8.2 ) 
 2090 FORMAT( /,' Solution : ', /, '                              ',          &
                '        <------ Bounds ------> ', /                          &
                '      # name       state    value   ',                       &
                '    Lower       Upper       Dual ' ) 
 2100 FORMAT( /, ' Of the ', I7, ' variables, ', 2X, I7,                      &
              ' are on bounds, &', I7, ' are dual degenerate' ) 
 2110 FORMAT( ' Of the ', I7, ' constraints, ', I7,' are equations, &', I7,   &
              ' are degenerate' ) 
 2120 FORMAT( ' Of the inequality constraints ', I6, ' are on bounds' ) 
 2130 FORMAT( I7, 1X, A10, A6, 4ES12.4 ) 
 2140 FORMAT( /, ' *** Problem will be scaled based on ', A1, ' *** ' )
 2150 FORMAT( ' Allocation error, variable ', A8, ' status = ', I6 )
 2160 FORMAT( ' IOSTAT = ', I6, ' when opening file ', A9, '. Stopping ' )
 2170 FORMAT( /, ' *** Further scaling applied to A *** ' )
 2180 FORMAT( A10 )
 2190 FORMAT( A10, 2I7, I8, I6, ES13.4, I6, 0P, F8.1, 1X, A1, A1 ) 
 2200 FORMAT( ' =%= old dimensions:  n = ', I7, ' m = ', I7,                   &
              ' a_ne = ', I9, ' h_ne = ', I9, /,                               &
              ' =%= new dimensions:  n = ', I7, ' m = ', I7,                   &
              ' a_ne = ', I9, ' h_ne = ', I9, /,                               &
              ' =%= preprocessing time =', F9.2,                               &
              '        number of transformations =', I10 )
 2210 FORMAT( ' === postprocessing time =', F9.2, /,                           &
              ' === processing time     =', F9.2 )
 2250 FORMAT( /, ' Problem:    ', A10, /, ' Solver :   ', A5,                  &
              /, ' Objective:', ES24.16 )

!  End of subroutine USE_QPA

     END SUBROUTINE USE_QPA

!  End of module USEQPA_double

   END MODULE GALAHAD_USEQPA_double
