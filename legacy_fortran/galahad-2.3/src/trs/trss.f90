   PROGRAM GALAHAD_TRS_EXAMPLE  !  GALAHAD 2.2 - 05/06/2008 AT 13:30 GMT.
   USE GALAHAD_TRS_DOUBLE                          ! double precision version
   IMPLICIT NONE
   INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )       ! set precision
   REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp
   INTEGER, PARAMETER :: n = 3, h_ne = 4           ! problem dimensions
   INTEGER :: s
   REAL ( KIND = wp ), DIMENSION( n ) :: C, X
   TYPE ( SMT_type ) :: H
   TYPE ( TRS_data_type ) :: data
   TYPE ( TRS_control_type ) :: control        
   TYPE ( TRS_inform_type ) :: inform
   REAL ( KIND = wp ) :: f = one                ! constant term, f
   REAL ( KIND = wp ) :: radius = one           ! trust-region radius
   C = (/ 0.0_wp, 2.0_wp, 0.0_wp /)             ! linear term, c
!  C = (/ 0.0_wp, 2.0_wp, 0.0001_wp /)             ! linear term, c
!  C = (/ 5.0_wp, 0.0_wp, 4.0_wp /)             ! linear term, c
   CALL SMT_put( H%type, 'COORDINATE', s )      ! Specify co-ordinate
   ALLOCATE( H%val( h_ne ), H%row( h_ne ), H%col( h_ne ) )
   H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp /) ! Hessian, H
!  H%val = (/ 1.0_wp, 2.0_wp, 3.0_wp, 1.0_wp /) ! Hessian, H
   H%row = (/ 1, 2, 3, 3 /)                     ! NB lower triangle
   H%col = (/ 1, 2, 3, 1 /) ; H%ne = h_ne
   CALL TRS_initialize( data, control )         ! Initialize control parameters
!  control%print_level = 1
!  control%initial_multiplier = 4.0_wp
   CALL TRS_solve( n, radius, f, C, H, X, data, control, inform ) ! solve problem
   IF ( inform%status == 0 ) THEN !  Successful return
     WRITE( 6, "( 1X, I0, ' factorizations. Solution and Lagrange multiplier =',&
    &    2ES12.4 )" ) inform%factorizations, inform%obj, inform%multiplier
   ELSE  !  Error returns
     WRITE( 6, "( ' TRS_solve exit status = ', I0 ) " ) inform%status
   END IF
   CALL TRS_terminate( data, control, inform )  ! delete internal workspace
   END PROGRAM GALAHAD_TRS_EXAMPLE
