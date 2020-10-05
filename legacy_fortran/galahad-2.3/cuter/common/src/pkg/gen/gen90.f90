!     ( Last modified on 23 Dec 2000 at 22:01:38 )

Module Generic_Driver

  Use CUTEr_precis
  Implicit None
  Private
  Public :: GEN, GENSPC, GETINFO

Contains

  Subroutine GEN( DUMMY )

    Real( Kind = wp ) :: DUMMY

    Write( *, * ) ' ********************************'
    Write( *, * )' *                              *'
    Write( *, * )' *      HELLO FROM GEN90!       *'
!S      WRITE( *, * )' *     (SINGLE PRECISION)       *'
!D      WRITE( *, * )' *     (DOUBLE PRECISION)       *'
    Write( *, * )' *                              *'
    Write( *, * )' ********************************'
    Write( *, * )' '
!S      DUMMY = 41.9999995555555E0
!D      DUMMY = 41.9999999999999D0
    Write( *, * ) ' OPTIMAL SOLUTION FOUND'
    Write( *, * ) ' THE ANSWER IS ', DUMMY
    Return
  End Subroutine GEN

  Subroutine GENSPC( FUNIT, FNAME )

    ! This is a dummy routine to read a spec file
    ! possibly, this routine contains precision-dependent directives

    Integer :: FUNIT
    Integer, Parameter :: FERROR = 6
    Character( len = 7 ) :: FNAME

    Open( UNIT=FUNIT, FILE=FNAME, STATUS='UNKNOWN', ERR=100 )
    Rewind( FUNIT )

    !     READ COMMANDS...

    Close( FUNIT )
    Return

100 Write( FERROR, '(A,A7)' ) 'Failure while reading ', FNAME
    Return
  End Subroutine GENSPC

  Subroutine GETINFO(N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS)
    !
    ! Input/Output variables
    !
    Integer, Intent( IN  ) :: N, M
    Integer, Intent( OUT ) :: NLIN, NEQ, NBNDS
    Real( Kind = wp ), Dimension( N ), Intent( IN ) :: BL, BU
    Logical, Dimension( M ), Intent( IN ) :: EQUATN, LINEAR
!
!     Local variables
!
!S  Real( Kind = wp ), Parameter :: INFTY = 1.0E+20
!D  Real( Kind = wp ), Parameter :: INFTY = 1.0D+20
    Integer :: I

    NLIN  = 0 ; NEQ   = 0 ; NBNDS = 0

    Do I = 1, M
       If( EQUATN( I ) ) NEQ  = NEQ  + 1
       If( LINEAR( I ) ) NLIN = NLIN + 1
    End Do

    Do I = 1, N
       If( BL( I ) .Gt. -INFTY .Or. BU( I ) .Lt. INFTY ) NBNDS = NBNDS + 1
    End Do
  End Subroutine GETINFO

End Module Generic_Driver
