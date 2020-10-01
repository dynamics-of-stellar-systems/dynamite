! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

PROGRAM FILTRANE_MAIN

!-------------------------------------------------------------------------------
!   U s e d   m o d u l e s   a n d   s y m b o l s
!-------------------------------------------------------------------------------

   USE GALAHAD_NLPT_double      ! the NLP problem type

   USE GALAHAD_FILTRANE_double  ! the FILTRANE solver

   USE GALAHAD_SYMBOLS,                                                        &
      OK                          => GALAHAD_SUCCESS,                          &
      MEMORY_FULL                 => GALAHAD_MEMORY_FULL,                      &
      SILENT                      => GALAHAD_SILENT,                           &
      TRACE                       => GALAHAD_TRACE,                            &
!     ACTION                      => GALAHAD_ACTION,                           &
!     DETAILS                     => GALAHAD_DETAILS,                          &
!     DEBUG                       => GALAHAD_DEBUG,                            &
!     CRAZY                       => GALAHAD_CRAZY,                            &
      COORDINATE                  => GALAHAD_COORDINATE,                       &
      FIXED                       => GALAHAD_FIXED,                            &
      RANGE                       => GALAHAD_RANGE,                            &
      LOWER                       => GALAHAD_LOWER,                            &
      UPPER                       => GALAHAD_UPPER,                            &
      FREE                        => GALAHAD_FREE,                             &
      USER_DEFINED                => GALAHAD_USER_DEFINED,                     &
      NONE                        => GALAHAD_NONE

!-------------------------------------------------------------------------------
!   A c c e s s 
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!  PRIVATE :: OK, MEMORY_FULL_ SILENT, TRACE, ACTION,                          &
!             DETAILS, DEBUG, CRAZY, COORDINATE, FIXED, LOWER, UPPER,          &
!             RANGE, FREE, USER_DEFINED, NONE

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
  
  INTEGER, PARAMETER :: ispec = 55      ! SPECfile  device number
  INTEGER, PARAMETER :: isif  = 56      ! OUTSDIF.d device number
  INTEGER, PARAMETER :: iout = 6        ! stderr and stdout
  
  REAL( KIND = wp ), PARAMETER :: INFINITY = (10.0_wp)**19

  INTEGER :: iostat

!-------------------------------------------------------------------------------
!  T h e   w o r k s
!-------------------------------------------------------------------------------

! Local variable

  INTEGER              :: status, n_bounded, m_not_equal, nnzj, J_ne_plus_n,   &
                          cuter_inform
  REAL, DIMENSION( 7 ) :: cuter_calls
  REAL, DIMENSION( 2 ) :: cuter_time

! Open the SIF description file

  OPEN( isif, file = 'OUTSDIF.d', form = 'FORMATTED', status = 'OLD',          &
        IOSTAT = iostat )
  IF ( iostat > 0 ) THEN
     WRITE( iout, 201 )  isif
     STOP
  END IF
  REWIND( isif )

! Setup the current CUTEr problem

  CALL cuter_initialize( problem, isif, iout, cuter_inform )
  CLOSE( isif )
  IF ( cuter_inform /= OK ) STOP

! Update J_ne to take the peculiarity of CUTEr into account.

  J_ne_plus_n = problem%J_ne + problem%n

  CALL NLPT_write_stats( problem, iout )

! Initialize FILTRANE

  CALL FILTRANE_initialize( FILTRANE_control, FILTRANE_inform, FILTRANE_data )

! Read the FILTRANE spec file

  OPEN( ispec, file = 'FILTRANE.SPC', form = 'FORMATTED', status = 'OLD',      &
        IOSTAT = iostat )
  IF ( iostat > 0 ) THEN
     WRITE( iout, 205 )  ispec
     STOP
  END IF
  CALL FILTRANE_read_specfile( ispec, FILTRANE_control, FILTRANE_inform )
  CLOSE( ispec )

! Check the preconditioning and external product options, as these are
! not desired when using the CUTEr interface.

  IF ( FILTRANE_control%prec_used == USER_DEFINED ) THEN
     FILTRANE_control%prec_used = NONE
  END IF

  IF ( FILTRANE_control%external_J_products ) THEN
     FILTRANE_control%external_J_products = .FALSE.
  END IF

! Apply the solver in a reverse communication loop

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

! Output results

  CALL NLPT_write_problem( problem, iout, TRACE )
!  CALL NLPT_write_problem( problem, iout, DETAILS )

  WRITE( iout, 202  )FILTRANE_inform%nbr_iterations,                           &
                     FILTRANE_inform%nbr_cg_iterations
  WRITE( iout, 203 ) FILTRANE_inform%nbr_c_evaluations
  WRITE( iout, 204 ) FILTRANE_inform%nbr_J_evaluations
  WRITE( iout, 200 ) problem%pname, problem%n, problem%m,                      &
                     FILTRANE_inform%nbr_c_evaluations,                        &
                     FILTRANE_inform%nbr_J_evaluations,                        &
                     FILTRANE_inform%status, problem%f,                        &
                     cuter_time( 1 ), cuter_time( 2 )
  WRITE( iout, 100 ) problem%pname, 'FILTRANE', FILTRANE_inform%nbr_iterations,&
                     FILTRANE_inform%nbr_cg_iterations, problem%f,             &
                     FILTRANE_inform%status, cuter_time( 1 ), cuter_time( 2 ), &
                     cuter_time( 1 ) + cuter_time( 2 )

! Clean up the problem space

  CALL NLPT_cleanup( problem )

  STOP

! Formats

100 FORMAT( /, ' Problem: ', A10, //,                                          &
               '                         CG      objective',                   &
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
201 FORMAT(1x,'ERROR: could not open file OUTSDIF.d as unit ',i2)
202 FORMAT(1x,'Number of iterations = ',i6,' Number of CG iterations = ',i10)
203 FORMAT(1x,'Number of constraints evaluations = ',i10)
204 FORMAT(1x,'Number of Jacobian evaluations    = ',i10)
205 FORMAT(1x,'ERROR: could not open file FILTRANE.SPC as unit ',i2)
206 FORMAT(1x,'ERROR: Jacobian products are requested to be internal.')
207 FORMAT(1x,'ERROR: preconditioner is requested to be internal.')

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

  INTEGER :: i, j, iostat, nmax, mmax, J_size, n_free

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

!-------------------------------------------------------------------------------
! Deallocate some useless vectors, if no bound is present
!-------------------------------------------------------------------------------

  ELSE
     NULLIFY( problem%z )
  END IF

!-------------------------------------------------------------------------------
! Deallocate some useless vectors, if no contraint is present
!-------------------------------------------------------------------------------

  IF ( problem%m <= 0 ) THEN
     DEALLOCATE( problem%y, problem%c, problem%c_l, problem%c_u,               &
                 problem%equation, problem%linear, problem%cnames )     
  END IF
 
!-------------------------------------------------------------------------------
! Nullify the derivative pointers for the Hessian
!-------------------------------------------------------------------------------

  NULLIFY( problem%H_val, problem%H_row, problem%H_col, problem%H_ptr,         &
           problem%gL)

  RETURN

!     Formats

101   FORMAT( 1x, i4, 1x, a10, 3x,  ES12.4 )
102   FORMAT( 1x, i4, 1x, a10, 3x, 2ES12.4 )
103   FORMAT( 1x, i4, 1x, a10, 3x, 3ES12.4 )
202   FORMAT( 1x, i4, 1x, a10, 3x, 2ES12.4, ' linear' )
203   FORMAT( 1x, i4, 1x, a10, 3x, 3ES12.4, ' linear' )
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

END PROGRAM FILTRANE_main

