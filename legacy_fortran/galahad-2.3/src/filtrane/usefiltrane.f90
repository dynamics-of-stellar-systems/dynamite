! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
! Updated 30/06/2003: Stop statement replaced by return

   MODULE GALAHAD_USEFILTRANE_double

!-------------------------------------------------------------------------------
!   U s e d   m o d u l e s   a n d   s y m b o l s
!-------------------------------------------------------------------------------

   USE GALAHAD_SYMBOLS,                                                        &
      OK                          => GALAHAD_SUCCESS,                          &
      MEMORY_FULL                 => GALAHAD_MEMORY_FULL,                      &
      SILENT                      => GALAHAD_SILENT,                           &
      TRACE                       => GALAHAD_TRACE,                            &
      DETAILS                     => GALAHAD_DETAILS,                          &
      COORDINATE                  => GALAHAD_COORDINATE,                       &
      USER_DEFINED                => GALAHAD_USER_DEFINED,                     &
      NONE                        => GALAHAD_NONE

   USE GALAHAD_NLPT_double      ! the NLP problem type

   USE GALAHAD_SPECFILE_double  ! the specfile tools

   USE GALAHAD_FILTRANE_double  ! the FILTRANE solver

!-------------------------------------------------------------------------------
!   A c c e s s 
!-------------------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE :: OK, MEMORY_FULL, SILENT, TRACE, DETAILS, COORDINATE,             &
              USER_DEFINED, NONE
              

   PUBLIC  :: USE_FILTRANE

 CONTAINS

!===============================================================================

   SUBROUTINE USE_FILTRANE( isif )

!  Reads a SIF problem and applies FILTRANE to it.

!  Argument
!  --------

   INTEGER, INTENT( IN ) :: isif

!         the number of the device on which the SIF problem file is opened.

!  Programming:  Ph. Toint and N. Gould, June 2003.

!===============================================================================


!-------------------------------------------------------------------------------
!   P r e c i s i o n
!-------------------------------------------------------------------------------

  INTEGER, PARAMETER :: sp = KIND( 1.0 )
  INTEGER, PARAMETER :: dp = KIND( 1.0D+0 )
  INTEGER, PARAMETER :: wp = dp

!-------------------------------------------------------------------------------
!   D e c l a r a t i o n s
!-------------------------------------------------------------------------------

  TYPE( NLPT_problem_type     ) :: problem
  TYPE( FILTRANE_control_type ) :: FILTRANE_control
  TYPE( FILTRANE_inform_type  ) :: FILTRANE_inform
  TYPE( FILTRANE_data_type    ) :: FILTRANE_data
  TYPE( SPECFILE_item_type    ), DIMENSION( 7 ) :: specs

  INTEGER, PARAMETER :: ispec = 55      ! SPECfile device number
  INTEGER, PARAMETER :: iout = 6        ! stdout
 
  REAL( KIND = wp ), PARAMETER :: INFINITY = (10.0_wp)**19

  INTEGER :: iostat
  LOGICAL :: filexx
  LOGICAL :: is_specfile

  INTEGER :: soldev = 57                ! solution file device number
  INTEGER :: sumdev = 58                ! summary file device number
  LOGICAL :: full_sol  = .FALSE.        
  LOGICAL :: write_sol = .FALSE.
  LOGICAL :: write_sum = .FALSE.
  INTEGER :: ierrout = 6                ! stderr
  CHARACTER ( LEN = 30 ) :: solfilename  = 'FILTRANE.sol'
  CHARACTER ( LEN = 30 ) :: sumfilename  = 'FILTRANE.sum'
  CHARACTER ( LEN = 16 ) :: specfilename = 'RUNFILT.SPC'
  CHARACTER ( LEN = 16 ) :: algo_name    = 'RUNFILT'

!-------------------------------------------------------------------------------
!  T h e   w o r k s
!-------------------------------------------------------------------------------

! Local variable

  INTEGER              :: nnzj, J_ne_plus_n, cuter_inform
  REAL, DIMENSION( 7 ) :: cuter_calls
  REAL, DIMENSION( 2 ) :: cuter_time

! Setup the current CUTEr problem

  CALL cuter_initialize( problem, isif, iout, cuter_inform )
  IF ( cuter_inform /= OK ) STOP

! Update J_ne to take the peculiarity of CUTEr into account.

  J_ne_plus_n = problem%J_ne + problem%n

! Initialize FILTRANE

  CALL FILTRANE_initialize( FILTRANE_control, FILTRANE_inform, FILTRANE_data )

! Open the specfile for RUNFILT and FILTRANE

  INQUIRE( FILE = specfilename, EXIST = is_specfile )
  IF ( is_specfile ) THEN
    OPEN( ispec, file = specfilename, form = 'FORMATTED', status = 'OLD',      &
          IOSTAT = iostat )
    IF ( iostat > 0 ) THEN
       WRITE( iout, 205 ) specfilename, ispec, iostat
       STOP
    END IF

! Define the keywords.

    specs( 1 )%keyword = 'print-full-solution'
    specs( 2 )%keyword = 'write-solution'
    specs( 3 )%keyword = 'solution-file-name'
    specs( 4 )%keyword = 'solution-file-device'
    specs( 5 )%keyword = 'write-result-summary'
    specs( 6 )%keyword = 'result-summary-file-name'
    specs( 7 )%keyword = 'result-summary-file-device'

! Read the specfile for RUNFILT.

    CALL SPECFILE_read( ispec, algo_name, specs, 7, ierrout )

! Interpret the result

    CALL SPECFILE_assign_logical( specs( 1 ), full_sol   , ierrout )
    CALL SPECFILE_assign_logical( specs( 2 ), write_sol  , ierrout )
    CALL SPECFILE_assign_string ( specs( 3 ), solfilename, ierrout )
    CALL SPECFILE_assign_integer( specs( 4 ), soldev     , ierrout )
    CALL SPECFILE_assign_logical( specs( 5 ), write_sum  , ierrout )
    CALL SPECFILE_assign_string ( specs( 6 ), sumfilename, ierrout )
    CALL SPECFILE_assign_integer( specs( 7 ), sumdev     , ierrout )

! Read the specfile for FILTRANE and close it.

    CALL FILTRANE_read_specfile( ispec, FILTRANE_control, FILTRANE_inform )
    CLOSE( ispec )
  END IF 

! Check the preconditioning and external product options, as these are
! not desired when using the CUTEr interface.

  IF ( FILTRANE_control%prec_used == USER_DEFINED ) THEN
     FILTRANE_control%prec_used = NONE
     WRITE( ierrout, 300 )
  END IF

  IF ( FILTRANE_control%external_J_products ) THEN
     FILTRANE_control%external_J_products = .FALSE.
     WRITE( ierrout, 301 )
  END IF

! Apply the solver in a reverse communication loop.

  DO 
  
     CALL FILTRANE_solve( problem, FILTRANE_control, FILTRANE_inform,          &
                          FILTRANE_data )

     SELECT CASE ( FILTRANE_inform%status )
  
     CASE ( 1, 2 )
  
        CALL CCFSG( problem%n, problem%m, problem%x, problem%m, problem%c,     &
                    nnzj, J_ne_plus_n, problem%J_val, problem%J_col,           &
                    problem%J_row, .TRUE. )
  
     CASE ( 3:5 )
  
        CALL CCFSG( problem%n, problem%m, problem%x, problem%m, problem%c,     &
                    nnzj, J_ne_plus_n, problem%J_val, problem%J_col,           &
                    problem%J_row, .FALSE. )
  
     CASE ( 6 )
  
        CALL CSGR(  problem%n, problem%m, .FALSE., problem%m, problem%y,       &
                    problem%x, nnzj, J_ne_plus_n, problem%J_val,               &
                    problem%J_col, problem%J_row )
  
     CASE ( 7 )

        WRITE( iout, 206 )
        EXIT

     CASE ( 8:11 )
  
        WRITE( iout, 206 )
        EXIT

     CASE ( 12:14 )

        WRITE( iout, 207 )
        EXIT

     CASE ( 15, 16 )
                    
        CALL CPROD( problem%n, problem%m, .NOT. FILTRANE_data%RC_newx,         &
                    problem%x, problem%m, problem%y, FILTRANE_data%RC_v,       &
                    FILTRANE_data%RC_Mv )

     CASE DEFAULT
  
        EXIT
  
     END SELECT
  
  END DO ! end of the reverse communication loop

! Get the CUTEr statistics.

  CALL CREPRT( cuter_calls, cuter_time )

! Terminate FILTRANE.

  FILTRANE_control%print_level = SILENT
  CALL FILTRANE_terminate( FILTRANE_control, FILTRANE_inform, FILTRANE_data )

! Output results.

  IF ( full_sol ) THEN
     CALL NLPT_write_problem( problem, iout, DETAILS )
  ELSE
     CALL NLPT_write_problem( problem, iout, TRACE )
  END IF

  WRITE( iout, 202  )FILTRANE_inform%nbr_iterations,                           &
                     FILTRANE_inform%nbr_cg_iterations
  WRITE( iout, 200 ) problem%pname, problem%n, problem%m,                      &
                     FILTRANE_inform%nbr_c_evaluations,                        &
                     FILTRANE_inform%nbr_J_evaluations,                        &
                     FILTRANE_inform%status, problem%f,                        &
                     cuter_time( 1 ), cuter_time( 2 )
  WRITE( iout, 100 ) 'FILTRANE', FILTRANE_inform%nbr_iterations,               &
                     FILTRANE_inform%nbr_cg_iterations, problem%f,             &
                     FILTRANE_inform%status, cuter_time( 1 ), cuter_time( 2 ), &
                     cuter_time( 1 ) + cuter_time( 2 )

! If required, write the solution to a file

  IF ( write_sol .AND. soldev > 0 ) THEN
     INQUIRE( FILE = solfilename, EXIST = filexx )
     IF ( filexx ) THEN
        OPEN( soldev, FILE = solfilename, FORM = 'FORMATTED',                  &
              STATUS = 'OLD', IOSTAT = iostat )
     ELSE
        OPEN( soldev, FILE = solfilename, FORM = 'FORMATTED',                  &
              STATUS = 'NEW', IOSTAT = iostat )
     END IF
     IF ( iostat /= 0 ) THEN 
        WRITE( iout, 205 ) solfilename, soldev, iostat
     ELSE
        CALL NLPT_write_problem( problem, soldev, DETAILS )
        CLOSE( soldev ) 
     END IF
  END IF 

! If required, write a result summary to a file

  IF ( write_sum .AND.sumdev > 0 ) THEN
     INQUIRE( FILE = sumfilename, EXIST = filexx )
     IF ( filexx ) THEN
        OPEN( sumdev, FILE = sumfilename, FORM = 'FORMATTED',                  &
              STATUS = 'OLD', IOSTAT = iostat )
     ELSE
        OPEN( sumdev, FILE = sumfilename, FORM = 'FORMATTED',                  &
              STATUS = 'NEW', IOSTAT = iostat )
     END IF
     IF ( iostat /= 0 ) THEN 
        WRITE( iout, 205 ) sumfilename, sumdev, iostat
     ELSE
        WRITE( sumdev, 208 ) problem%pname, problem%n, problem%m,              &
                             FILTRANE_inform%nbr_iterations,                   &
                             FILTRANE_inform%nbr_cg_iterations, problem%f,     &
                             FILTRANE_inform%status, cuter_time( 2 )
        CLOSE( sumdev ) 
     END IF
  END IF 

! Clean up the problem space

  CALL NLPT_cleanup( problem )

   RETURN

! Formats

100 FORMAT( /, '                         CG      objective',                   &
               '          < ------ time ----- > ', /,                          &
               ' Method  iterations  iterations    value  ',                   &
               '   status setup   solve   total', /,                           &
               ' ------  ----------  ----------  ---------',                   &
               '   ------ -----    ----   -----  ',/,                          &
                A8, 2I10, 3X, ES12.4, I6, 0P, 3F8.2 ) 
200 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //                       &
            ,' Code used               :  FILTRANE',    /                      &
            ,' Problem                 :  ', A10,    /                         &
            ,' # variables             =      ', I10 /                         &
            ,' # constraints           =      ', I10 /                         &
            ,' # constraints functions =      ', I10 /                         &
            ,' # constraints gradients =      ', I10 /                         &
             ' Exit code               =      ', I10 /                         &
            ,' Final f                 = ', 1pE15.7 /                          &
            ,' Set up time             =        ', 0P, F8.2, ' seconds' /      &
             ' Solve time              =        ', 0P, F8.2, ' seconds' //     &
            66('*') / )
202 FORMAT(1x,'Number of iterations = ',i6,' Number of CG iterations = ',i10)
205 FORMAT(1x,'ERROR: could not open file ',a,' as unit ',i2,' (IOSTAT = ', i6,&
           ')' )
206 FORMAT(1x,'ERROR: Jacobian products are requested to be internal.')
207 FORMAT(1x,'ERROR: preconditioner is requested to be internal.')
208 FORMAT(a10,2x,i10,1x,i10,1x,i10,1x,i10,1x,1pe11.3,1x,i3,1x,0p,f8.2)
300 FORMAT(1x,'WARNING: usefiltrane does not support USER_DEFINED ',           &
              'preconditioners.',/,'         Abandoning preconditioning.')
301 FORMAT(1x,'WARNING: usefiltrane does not support external Jacobian ',      &
              'products ',/, '         Using internal products.')

CONTAINS

!===============================================================================

      SUBROUTINE cuter_initialize( problem, isif, errout, inform_status )

!     Initializes the problem from its CUTEr description.

!     Arguments

      TYPE ( NLPT_problem_type ), INTENT( OUT ) :: problem
 
!            the problem;

      INTEGER, INTENT( IN ) :: isif

!            the device file for the OUTSDIF.d file

      INTEGER, INTENT( IN ) :: errout

!            the device number for error disagnostics;

      INTEGER, INTENT( OUT ) :: inform_status

!            the exit code.

!     Programming: Ph. Toint, November 2002.
!
!===============================================================================

! Local variables

  INTEGER :: j, iostat, nmax, mmax, J_size, n_free

!-------------------------------------------------------------------------------
! Initialize the exit status.
!-------------------------------------------------------------------------------

  inform_status = OK

!-------------------------------------------------------------------------------
! Set the infinity value.
!-------------------------------------------------------------------------------

  problem%infinity = INFINITY

! --------------------------------------------------------------------------
! Get the problem's dimensions
! --------------------------------------------------------------------------
  
  CALL CDIMEN( isif, nmax, mmax )

! --------------------------------------------------------------------------
! Allocate the problem's structure
! --------------------------------------------------------------------------
  
  ALLOCATE( problem%x( nmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 300 ) nmax
     RETURN
  END IF
  
  ALLOCATE( problem%x_l( nmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 301 )  nmax
     RETURN
  END IF
  
  ALLOCATE( problem%x_u( nmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 302 ) nmax
     RETURN
  END IF
  
  ALLOCATE( problem%x_status( nmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 303 ) nmax
     RETURN
  END IF

  ALLOCATE( problem%vnames( nmax ) )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 304 ) nmax
     RETURN
  END IF
  
  ALLOCATE( problem%equation( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 305 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%linear( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 306 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%c( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 307 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%c_l( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 308 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%c_u( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 309 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%y( mmax ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 310 ) mmax
     RETURN
  END IF
  
  ALLOCATE( problem%cnames( mmax ) )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 311 ) mmax
     RETURN
  END IF
  
! --------------------------------------------------------------------------
! CUTEr setup
! --------------------------------------------------------------------------
  
  CALL CSETUP( isif, errout, problem%n, problem%m, problem%x,                  &
               problem%x_l, problem%x_u, nmax, problem%equation,               &
               problem%linear, problem%y, problem%c_l, problem%c_u,            &
               mmax, .FALSE., .FALSE., .FALSE. )
  
  CALL CNAMES( problem%n, problem%m,                                           &
               problem%pname, problem%vnames, problem%cnames )
  
! --------------------------------------------------------------------------
! Allocate the Jacobian space.
! --------------------------------------------------------------------------
  
  CALL CDIMSJ( J_size )

  problem%J_ne = J_size - problem%n

  ALLOCATE( problem%J_val( J_size ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 312 )  J_size
     RETURN
  END IF
  
  ALLOCATE( problem%J_col( J_size ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 313 ) J_size
     RETURN
  END IF
  
  ALLOCATE( problem%J_row( J_size ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 314 ) J_size
     RETURN
  END IF
  problem%J_type = COORDINATE
  
!---------------------------------------------------------------------------
!    The gradient
!---------------------------------------------------------------------------
  
  ALLOCATE( problem%g( problem%n ), STAT = iostat )
  IF ( iostat /= 0 ) THEN
     inform_status = MEMORY_FULL
     WRITE( errout, 315 ) problem%n
     RETURN
   END IF

!-------------------------------------------------------------------------------
! Analyze the problem variables.
!-------------------------------------------------------------------------------

  n_free  = 0
  DO j = 1, problem%n
     IF ( problem%x_l( j ) <= - INFINITY .AND. problem%x_u( j ) >= INFINITY )  &
        n_free  = n_free + 1
  END DO

!-------------------------------------------------------------------------------
! If they are not all free, allocate the vector of dual variables
!-------------------------------------------------------------------------------

  IF ( n_free < problem%n  ) THEN

     ALLOCATE( problem%z( problem%n ), STAT = iostat )
     IF ( iostat /= 0 ) THEN
        inform_status = MEMORY_FULL
        WRITE( errout, 316 ) problem%n
        RETURN
     END IF
     problem%z = 0.0_wp

  END IF

!-------------------------------------------------------------------------------
! Deallocate some useless vectors, if no contraint is present
!-------------------------------------------------------------------------------

  IF ( problem%m <= 0 ) THEN
     DEALLOCATE( problem%y, problem%c, problem%c_l, problem%c_u,               &
                 problem%equation, problem%linear, problem%cnames )     
  END IF
 
  RETURN

!     Formats

300   FORMAT( 1x, 'ERROR: no memory for allocating x(',i6, ')')
301   FORMAT( 1x, 'ERROR: no memory for allocating x_l(',i6, ')')
302   FORMAT( 1x, 'ERROR: no memory for allocating x_u(',i6, ')')
303   FORMAT( 1x, 'ERROR: no memory for allocating x_status(',i6, ')')
304   FORMAT( 1x, 'ERROR: no memory for allocating vnames(',i6, ')')
305   FORMAT( 1x, 'ERROR: no memory for allocating equation(',i6, ')')
306   FORMAT( 1x, 'ERROR: no memory for allocating linear(',i6, ')')
307   FORMAT( 1x, 'ERROR: no memory for allocating c(',i6, ')')
308   FORMAT( 1x, 'ERROR: no memory for allocating c_l(',i6, ')')
309   FORMAT( 1x, 'ERROR: no memory for allocating c_u(',i6, ')')
310   FORMAT( 1x, 'ERROR: no memory for allocating y(',i6, ')')
311   FORMAT( 1x, 'ERROR: no memory for allocating cnames(',i6, ')')
312   FORMAT( 1x, 'ERROR: no memory for allocating J_val(',i10, ')')
313   FORMAT( 1x, 'ERROR: no memory for allocating J_col(',i10, ')')
314   FORMAT( 1x, 'ERROR: no memory for allocating J_row(',i10, ')')
315   FORMAT( 1x, 'ERROR: no memory for allocating g(',i6, ')')
316   FORMAT( 1x, 'ERROR: no memory for allocating z(',i6, ')')

   END SUBROUTINE cuter_initialize

!===============================================================================
!===============================================================================

   END SUBROUTINE USE_FILTRANE

   END MODULE GALAHAD_USEFILTRANE_double

