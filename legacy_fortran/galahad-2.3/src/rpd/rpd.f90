! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-*- G A L A H A D _ R P D   M O D U L E -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.0. January 22nd 2006

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_RPD_double

!     ----------------------------------------------
!     |                                            |
!     | Read data for the quadratic program        |
!     |                                            |
!     |    minimize     1/2 x(T) H x + g(T) x + f  |
!     |    subject to     c_l <= A x <= c_u        |
!     |                   x_l <=  x  <= x_u        |
!     |                                            |
!     | from a problem-data file                   |
!     |                                            |
!     ----------------------------------------------

      USE GALAHAD_SMT_double, only: SMT_put
      USE GALAHAD_QPT_double

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: RPD_read_problem_data, QPT_problem_type

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      INTEGER, PARAMETER :: input_line_length = 256

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

      TYPE, PUBLIC :: RPD_inform_type
        INTEGER :: status, alloc_status, io_status, line
        CHARACTER ( LEN = 1 ) :: bad_alloc
      END TYPE

   CONTAINS

!-*-*-   R P D _ R E A D _ P R O B L E M _ D A T A   S U B R O U T I N E   -*-*-

      SUBROUTINE RPD_read_problem_data( input, prob, inform )
      INTEGER, INTENT( IN ) :: input
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( RPD_inform_type ), INTENT( out ) :: inform

!  Read the problem-data file from unit input into the derived type prob

!  ****************************************************************************

!  The data should be input in a file on unit 5. The data is in
!  free format (blanks separate values),but must occur in the order 
!  given here with NO SEPARATING BLANK LINES. Each term in "quotes"
!  denotes a required value. Any strings beyond those required on a
!  given lines will be regarded as comments and ignored.

!  "problem name"
!  "number varibales"
!  "number general linear constraints"
!  "number of nonzeros in lower traingle of H"
!  "row" "column" "value" for each entry of H (if any), one triple on each line
!  "default value for entries in g"
!  "number of non-default entries in g"
!  "index" "value" for each non-default term in g (if any), one pair per line
!  "value of f"
!  "number of nonzeros in A"
!  "row" "column" "value" for each entry of A (if any), one triple on each line
!  "default value for entries in c_l"
!  "number of non-default entries in c_l"
!  "index" "value" for each non-default term in c_l (if any), one pair per line
!  "default value for entries in c_u"
!  "number of non-default entries in c_u"
!  "index" "value" for each non-default term in c_u (if any), one pair per line
!  "default value for entries in x_l"
!  "number of non-default entries in x_l"
!  "index" "value" for each non-default term in x_l (if any), one pair per line
!  "default value for entries in x_u"
!  "number of non-default entries in x_u"
!  "index" "value" for each non-default term in x_u (if any), one pair per line
!  "default value for starting value for variables x"
!  "number of non-default starting entries in x"
!  "index" "value" for each non-default term in x (if any), one pair per line
!  "default value for starting value for Lagrange multipliers y for constraints"
!  "number of non-default starting entries in y"
!  "index" "value" for each non-default term in y (if any), one pair per line
!  "default value for starting value for dual varibales z for simple bounds"
!  "number of non-default starting entries in z"
!  "index" "value" for each non-default term in z (if any), one pair per line

!  *****************************************************************************

!  Local variables

     INTEGER :: i, j, k, A_ne, H_ne, nnzx_0, nnzy_0, nnzz_0
     INTEGER :: nnzg, nnzc_l, nnzc_u, nnzx_l, nnzx_u, smt_stat
     REAL ( KIND = wp ) :: rv, default
     CHARACTER ( LEN = 10 ) :: pname
     CHARACTER ( LEN = input_line_length ) :: input_line, blank_line

     inform%line = 0
     inform%alloc_status = 0 ; inform%bad_alloc = ' '

     DO i = 1, input_line_length
       blank_line( i : i ) = ' '
     END DO

!  Determine the problem name

     pname = '          '
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) pname 
     ALLOCATE( prob%name( 10 ) )
     prob%name = TRANSFER( pname, prob%name )

!  Determine the number of variables and constraints

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%n
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%m
 
!  Allocate suitable arrays

     ALLOCATE( prob%X( prob%n ), prob%X_l( prob%n ), prob%X_u( prob%n ),       &
               prob%G( prob%n ), prob%Z( prob%n ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'X' ; RETURN
     END IF

     ALLOCATE( prob%C_l( prob%m ), prob%C_u( prob%m ), prob%Y( prob%m ),       &
               prob%C( prob%m ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'Y' ; RETURN
     END IF

!  Fill component H

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO

     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%H%ne
     ALLOCATE( prob%H%row( prob%H%ne ), prob%H%col( prob%H%ne ),               &
               prob%H%val( prob%H%ne ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'H' ; RETURN
     END IF

     H_ne = 0
     DO k = 1, prob%H%ne
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, j, rv
       IF ( rv == zero ) CYCLE
       H_ne = H_ne + 1 ; prob%H%val( H_ne ) = rv
       IF ( i >= j ) THEN
         prob%H%row( H_ne ) = i
         prob%H%col( H_ne ) = j
       ELSE
         prob%H%row( H_ne ) = j
         prob%H%col( H_ne ) = i
       END IF
     END DO
     prob%H%ne = H_ne
     IF ( ALLOCATED( prob%H%type ) ) DEALLOCATE( prob%H%type )
     CALL SMT_put( prob%H%type, 'COORDINATE', smt_stat )

!  Fill component g

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%G = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzg
     DO k = 1, nnzg
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%G( i ) = rv
     END DO

!  Fill component f

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%f

!  Fill component A

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%A%ne
     ALLOCATE( prob%A%row( prob%A%ne ), prob%A%col( prob%A%ne ),               &
               prob%A%val( prob%A%ne ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'A' ; RETURN
     END IF

     A_ne = 0
     DO k = 1, prob%A%ne
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, j, rv
       IF ( rv == zero ) CYCLE
       A_ne = A_ne + 1 ; prob%A%val( A_ne ) = rv
       prob%A%row( A_ne ) = i ; prob%A%col( A_ne ) = j
     END DO
     prob%A%ne = A_ne
     IF ( ALLOCATED( prob%A%type ) ) DEALLOCATE( prob%A%type )
     CALL SMT_put( prob%A%type, 'COORDINATE', smt_stat )

!  Fill component c_l

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%C_l = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzc_l
     DO k = 1, nnzc_l
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%C_l( i ) = rv
     END DO

!  Fill component c_u

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%C_u = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzc_u
     DO k = 1, nnzc_u
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%C_u( i ) = rv
     END DO

!  Fill component x_l

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%X_l = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzx_l
     DO k = 1, nnzx_l
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%X_l( i ) = rv
     END DO

!  Fill component x_u

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%X_u = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzx_u
     DO k = 1, nnzx_u
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%X_u( i ) = rv
     END DO

!  Fill component x

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%X = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzx_0
     DO k = 1, nnzx_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%X( i ) = rv
     END DO

!  Fill component y

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%Y = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzy_0
     DO k = 1, nnzy_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%Y( i ) = rv
     END DO

!  Fill component z

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%Z = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzz_0
     DO k = 1, nnzz_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%Z( i ) = rv
     END DO

     inform%status = 0
     RETURN  

!  Error returns

!  - end of file encountered

 930 CONTINUE
     inform%status = 3
     RETURN  

!  - other error encountered

 940 CONTINUE
     inform%status = 4
     RETURN  

!  End of RPD_read_problem_data

     END SUBROUTINE RPD_read_problem_data

!-*-*-*-*-*-   R P D _ I G N O R E _ S T R I N G   F U N C T I O N   -*-*-*-*-*-

     FUNCTION RPD_ignore_string( input_line )
     LOGICAL :: RPD_ignore_string

!  Ignore a string if it is (a) blank or (b) starts with "!", "%" or "#"

     CHARACTER ( LEN = input_line_length ), INTENT( IN ) :: input_line

!  Local variables

     INTEGER :: i, length_string

     length_string = LEN_TRIM( input_line )
     IF ( length_string <= 0 ) THEN
       RPD_ignore_string = .TRUE.
       RETURN
     END IF

     DO i = 1, length_string
       IF ( input_line( i : i ) == ' ' ) CYCLE
       IF ( input_line( i : i ) == '!' .OR. input_line( i : i ) == '#' .OR. &
            input_line( i : i ) == '%' .OR. input_line( i : i ) == '|' ) THEN
         RPD_ignore_string = .TRUE.
         RETURN
       END IF
       EXIT
     END DO
     RPD_ignore_string = .FALSE.

     RETURN

!  End of RPD_ignore_string

     END FUNCTION RPD_ignore_string

!  End of module RPD

   END MODULE GALAHAD_RPD_double
