C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      SUBROUTINE GEN( DUMMY )
C
CS  THIS IS THE SINGLE PRECISION VERSION OF THE GENERIC PACKAGE
CD  THIS IS THE DOUBLE PRECISION VERSION OF THE GENERIC PACKAGE
C
CS      REAL              DUMMY
CD      DOUBLE PRECISION  DUMMY
      WRITE(*,*)' ********************************'
      WRITE(*,*)' *                              *'
      WRITE(*,*)' *       HELLO FROM GEN!        *'
CS      WRITE(*,*)' *     (SINGLE PRECISION)       *'
CD      WRITE(*,*)' *     (DOUBLE PRECISION)       *'
      WRITE(*,*)' *                              *'
      WRITE(*,*)' ********************************'
      WRITE(*,*)' '
CS      DUMMY = 41.9999995555555E0
CD      DUMMY = 41.9999999999999D0
      WRITE( *, * ) ' OPTIMAL SOLUTION FOUND'
      WRITE( *, * ) ' THE ANSWER IS ', DUMMY
      RETURN
      END

      SUBROUTINE GENSPC( FUNIT, FNAME )

C     THIS IS A DUMMY ROUTINE TO READ A SPEC FILE
C     POSSIBLY, THIS ROUTINE CONTAINS PRECISION-DEPENDENT DIRECTIVES

      INTEGER     FUNIT, FERROR
      PARAMETER( FERROR = 6 )
      CHARACTER*7 FNAME

      OPEN( UNIT=FUNIT, FILE=FNAME, STATUS='UNKNOWN', ERR=100 )
      REWIND( FUNIT )

C     READ COMMANDS...

      CLOSE( FUNIT )
      RETURN

 100  WRITE( FERROR, '(A,A7)' ) 'Failure while reading ', FNAME
      RETURN

      END

      SUBROUTINE GETINFO(N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS)
C
C     Input/Output variables
C
      INTEGER N, M, NLIN, NEQ, NBNDS
CS      REAL             BL( N ), BU( N )
CS      REAL             INFTY
CS      PARAMETER      ( INFTY = 1.0E+20 )
CD      DOUBLE PRECISION BL( N ), BU( N )
CD      DOUBLE PRECISION INFTY
CD      PARAMETER      ( INFTY = 1.0D+20 )
      LOGICAL EQUATN( M ), LINEAR( M )
C
C     Local variables
C
      INTEGER I

      NLIN  = 0
      NEQ   = 0
      NBNDS = 0

      DO 200 I = 1, M
         IF( EQUATN( I ) ) NEQ  = NEQ  + 1
         IF( LINEAR( I ) ) NLIN = NLIN + 1
 200  CONTINUE

      DO 300 I = 1, N
         IF( BL( I ) .GT. -INFTY .OR. BU( I ) .LT. INFTY ) 
     .        NBNDS = NBNDS + 1
 300  CONTINUE

      END
