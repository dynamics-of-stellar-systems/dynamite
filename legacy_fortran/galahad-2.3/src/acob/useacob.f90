! THIS VERSION: GALAHAD 2.2 - 07/02/2008 AT 17:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D   U S E _ A C O B  -*-*-*-*-*-*-*-*-*-*-

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  Started: February 2nd 2008

   MODULE GALAHAD_USEACOB_double

!  This is the driver program for running ACOB for a variety of computing 
!  systems. It opens and closes all the files, allocate arrays, reads and 
!  checks data, and calls the appropriate minimizers

     USE GALAHAD_ACOB_double
     USE GALAHAD_SPECFILE_double 
     USE GALAHAD_COPYRIGHT
     USE GALAHAD_SPACE_double
     USE GALAHAD_CUTER_FUNCTIONS_double
     IMPLICIT NONE

     PRIVATE
     PUBLIC :: USE_ACOB

   CONTAINS

!-*-*-*-*-*-*-*-*-*-  U S E _ A C O B   S U B R O U T I N E  -*-*-*-*-*-*-*-

     SUBROUTINE USE_ACOB( input )

!  Dummy argument

     INTEGER, INTENT( IN ) :: input

!  Set precision

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!-------------------------------
!   D e r i v e d   T y p e s
!-------------------------------

     TYPE ( ACOB_control_type ) :: control
     TYPE ( ACOB_inform_type ) :: inform
     TYPE ( ACOB_data_type ) :: data
     TYPE ( NLPT_problem_type ) :: nlp
     TYPE ( NLPT_userdata_type ) :: userdata
     TYPE ( CUTER_FUNCTIONS_control_type ) :: cuter_control
     TYPE ( CUTER_FUNCTIONS_inform_type ) :: cuter_inform

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------

!  Problem input characteristics

     INTEGER :: iores, i
     LOGICAL :: filexx, is_specfile

!  Specfile characteristics

     INTEGER, PARAMETER :: input_specfile = 34
     INTEGER, PARAMETER :: lspec = 29
     CHARACTER ( LEN = 16 ) :: specname = 'RUNACOB'
     TYPE ( SPECFILE_item_type ), DIMENSION( lspec ) :: spec
     CHARACTER ( LEN = 16 ) :: runspec = 'RUNACOB.SPC'

!  Default values for specfile-defined parameters

     INTEGER :: dfiledevice = 26
     INTEGER :: rfiledevice = 47
     INTEGER :: sfiledevice = 62
     INTEGER :: wfiledevice = 59
     LOGICAL :: write_problem_data   = .FALSE.
     LOGICAL :: write_solution       = .FALSE.
!    LOGICAL :: write_result_summary = .FALSE.
     LOGICAL :: write_result_summary = .TRUE.
     CHARACTER ( LEN = 30 ) :: dfilename = 'ACOB.data'
     CHARACTER ( LEN = 30 ) :: rfilename = 'ACOBRES.d'
     CHARACTER ( LEN = 30 ) :: sfilename = 'ACOBSOL.d'
     CHARACTER ( LEN = 30 ) :: wfilename = 'ACOBSAVE.d'
     LOGICAL :: testal = .FALSE.
     LOGICAL :: dechk  = .FALSE.
     LOGICAL :: dechke = .FALSE.
     LOGICAL :: dechkg = .FALSE.
     LOGICAL :: not_fatal  = .FALSE.
     LOGICAL :: not_fatale = .FALSE.
     LOGICAL :: not_fatalg = .FALSE.
     LOGICAL :: getsca = .FALSE.
     INTEGER :: print_level_scaling = 0
     LOGICAL :: scale  = .FALSE.
     LOGICAL :: scaleg = .FALSE.
     LOGICAL :: scalev = .FALSE.
     LOGICAL :: get_max = .FALSE.
     LOGICAL :: warm_start = .FALSE.
     INTEGER :: istore = 0
!    LOGICAL :: one_norm = .TRUE.

!  Output file characteristics

     INTEGER :: out  = 6
     INTEGER :: errout = 6
     CHARACTER ( LEN =  6 ) :: solv = 'ACOB   '

!  ------------------ Open the specfile for runlpsqp ----------------

     INQUIRE( FILE = runspec, EXIST = is_specfile )
     IF ( is_specfile ) THEN
       OPEN( input_specfile, FILE = runspec, FORM = 'FORMATTED', STATUS = 'OLD' )

!   Define the keywords

       spec( 1 )%keyword  = 'write-problem-data'
       spec( 2 )%keyword  = 'problem-data-file-name'
       spec( 3 )%keyword  = 'problem-data-file-device'
       spec( 4 )%keyword  = ''
       spec( 5 )%keyword  = 'write-solution'
       spec( 6 )%keyword  = 'solution-file-name'
       spec( 7 )%keyword  = 'solution-file-device'
       spec( 8 )%keyword  = 'write-result-summary'
       spec( 9 )%keyword  = 'result-summary-file-name'
       spec( 10 )%keyword = 'result-summary-file-device'
       spec( 11 )%keyword = 'check-all-derivatives'
       spec( 12 )%keyword = 'check-derivatives'
       spec( 13 )%keyword = 'check-element-derivatives'
       spec( 14 )%keyword = 'check-group-derivatives'
       spec( 15 )%keyword = 'ignore-derivative-bugs'
       spec( 16 )%keyword = 'ignore-element-derivative-bugs'
       spec( 17 )%keyword = 'ignore-group-derivative-bugs'
       spec( 18 )%keyword = 'get-scaling-factors' 
       spec( 19 )%keyword = 'scaling-print-level' 
       spec( 20 )%keyword = 'use-scaling-factors' 
       spec( 21 )%keyword = 'use-constraint-scaling-factors' 
       spec( 22 )%keyword = 'use-variable-scaling-factors' 
       spec( 23 )%keyword = 'maximizer-sought' 
       spec( 24 )%keyword = 'restart-from-previous-point' 
       spec( 25 )%keyword = 'restart-data-file-name'
       spec( 26 )%keyword = 'restart-data-file-device'
       spec( 27 )%keyword = 'save-data-for-restart--every'
       spec( 28 )%keyword = ''
       spec( 29 )%keyword = ''

!   Read the specfile

       CALL SPECFILE_read( input_specfile, specname, spec, lspec, errout )

!   Interpret the result

       CALL SPECFILE_assign_logical( spec( 1 ), write_problem_data, errout )
       CALL SPECFILE_assign_string ( spec( 2 ), dfilename, errout )
       CALL SPECFILE_assign_integer( spec( 3 ), dfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 5 ), write_solution, errout )
       CALL SPECFILE_assign_string ( spec( 6 ), sfilename, errout )
       CALL SPECFILE_assign_integer( spec( 7 ), sfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 8 ), write_result_summary, errout )
       CALL SPECFILE_assign_string ( spec( 9 ), rfilename, errout )
       CALL SPECFILE_assign_integer( spec( 10 ), rfiledevice, errout )
       CALL SPECFILE_assign_logical( spec( 11 ), testal, errout )
       CALL SPECFILE_assign_logical( spec( 12 ), dechk, errout )
       CALL SPECFILE_assign_logical( spec( 13 ), dechke, errout )
       CALL SPECFILE_assign_logical( spec( 14 ), dechkg, errout )
       CALL SPECFILE_assign_logical( spec( 15 ), not_fatal, errout )
       CALL SPECFILE_assign_logical( spec( 16 ), not_fatale, errout )
       CALL SPECFILE_assign_logical( spec( 17 ), not_fatalg, errout )
       CALL SPECFILE_assign_logical( spec( 18 ), getsca, errout )
       CALL SPECFILE_assign_integer( spec( 19 ), print_level_scaling, errout )
       CALL SPECFILE_assign_logical( spec( 20 ), scale, errout )
       CALL SPECFILE_assign_logical( spec( 21 ), scaleg, errout )
       CALL SPECFILE_assign_logical( spec( 22 ), scalev, errout )
       CALL SPECFILE_assign_logical( spec( 23 ), get_max, errout )
       CALL SPECFILE_assign_logical( spec( 24 ), warm_start, errout )
       CALL SPECFILE_assign_string ( spec( 25 ), wfilename, errout )
       CALL SPECFILE_assign_integer( spec( 26 ), wfiledevice, errout )
       CALL SPECFILE_assign_integer( spec( 27 ), istore, errout )
     END IF
    
     IF ( dechk .OR. testal ) THEN ; dechke = .TRUE. ; dechkg = .TRUE. ; END IF
     IF ( not_fatal ) THEN ; not_fatale = .TRUE. ; not_fatalg = .TRUE. ; END IF
     IF ( scale ) THEN ; scaleg = .TRUE. ; scalev = .TRUE. ; END IF

!  If required, open a file for the results

     IF ( write_result_summary ) THEN
       INQUIRE( FILE = rfilename, EXIST = filexx )
       IF ( filexx ) THEN
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',             &
                STATUS = 'OLD', POSITION = 'APPEND', IOSTAT = iores )
       ELSE
          OPEN( rfiledevice, FILE = rfilename, FORM = 'FORMATTED',             &
                STATUS = 'NEW', IOSTAT = iores )
       END IF
       IF ( iores /= 0 ) THEN 
         write( errout, 2030 ) iores, rfilename
         STOP
       END IF
       READ( INPUT, "( /, I2, A8  )" ) iores, nlp%pname
       REWIND( input )
       WRITE( rfiledevice, "( A10 )" ) nlp%pname
     END IF

!  Initialize the problem data

     cuter_control%input = input ; cuter_control%error = errout
     CALL CUTER_initialize( nlp, cuter_control, cuter_inform, userdata )

!  Set copyright 

     IF ( out > 0 ) CALL COPYRIGHT( out, '2007' )

!  Set up data for next problem

     CALL ACOB_initialize( data, control )
     IF ( is_specfile ) CALL ACOB_read_specfile( control, input_specfile )

!  Solve the problem

     CALL ACOB_solve( nlp, control, inform, data,                              &
                      CUTER_eval_F, CUTER_eval_G, CUTER_eval_H, userdata )

!  If required, append results to a file

      IF ( write_result_summary ) THEN
        BACKSPACE( rfiledevice )
        WRITE( rfiledevice, "( A10, ES16.8, ES9.1, bn, I9, F12.2, I6 )" )      &
          nlp%pname, inform%obj, inform%norm_pg,                               &
          inform%iter, inform%time%total, inform%status
      END IF
      WRITE( errout, "( 'name        f               du-feas ',                &
     &                  '     its        time  stat' )" )
      WRITE( errout, "( A10, ES16.8, ES9.1, bn, I9, F12.2, I6 )" )             &
        nlp%pname, inform%obj, inform%norm_pg,                                 &
        inform%iter, inform%time%total, inform%status

!  If required, write the solution

     IF ( write_solution .AND.                                                 &
         ( inform%status == 0  .OR. inform%status == - 10 ) ) THEN

       INQUIRE( FILE = sfilename, EXIST = filexx )
       IF ( filexx ) THEN
          OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',             &
              STATUS = 'OLD', IOSTAT = iores )
       ELSE
          OPEN( sfiledevice, FILE = sfilename, FORM = 'FORMATTED',             &
               STATUS = 'NEW', IOSTAT = iores )
       END IF
       IF ( iores /= 0 ) THEN 
         write( out, 2030 ) iores, sfilename
         STOP
       END IF

       WRITE( sfiledevice, "( /, ' Problem:    ', A10, /, ' Solver :   ', A,   &
      &       /, ' Objective:', ES24.16 )" ) nlp%pname, solv, inform%obj

       WRITE( sfiledevice, 2000 )
       DO i = 1, nlp%n
         WRITE( sfiledevice, 2020 ) i, nlp%VNAMES( i ), nlp%X( i ),            &
           nlp%X_l( i ), nlp%X_u( i ), nlp%Z( i )
       END DO

    END IF

!  Close any opened files and deallocate arrays

     IF ( is_specfile ) CLOSE( input_specfile )
     CALL CUTER_terminate( nlp, cuter_inform, userdata )
     STOP

!  Non-executable statements

 2000 FORMAT( /,' Solution: ', /,'                        ',                   &
                '        <------ Bounds ------> ', /                           &
                '      # name          value   ',                              &
                '    Lower       Upper       Dual ' ) 
 2020 FORMAT( I7, 1X, A10, 4ES12.4 ) 
 2030 FORMAT( ' IOSTAT = ', I6, ' when opening file ', A9, '. Stopping ' )

!  End of subroutine USE_ACOB

     END SUBROUTINE USE_ACOB

!  End of module USEACOB_double

   END MODULE GALAHAD_USEACOB_double
