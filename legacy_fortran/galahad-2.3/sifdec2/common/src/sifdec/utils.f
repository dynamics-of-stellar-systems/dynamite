C     ( Last modified on 15 Mar 2001 at 22:28:00 )
      SUBROUTINE REORDA( NC    , NNZ   , IRN   , JCN   , A , IP, IW )
      INTEGER            NC    , NNZ
      INTEGER            IRN   ( NNZ  ), JCN   ( NNZ )
      INTEGER            IW    ( *    ), IP    ( *   )
      DOUBLE PRECISION   A     ( NNZ  )
C
C  SORT A SPARSE MATRIX FROM ARBITRARY ORDER TO COLUMN ORDER.
C  NICK GOULD. 7TH NOVEMBER, 1990.
C  FOR CGT PRODUCTIONS.
C
      INTEGER          I , J , K , L , IC
      INTEGER          NCP1  , ITEMP , JTEMP,  LOCAT
      DOUBLE PRECISION ANEXT , ATEMP
C
C  INITIALIZE THE WORKSPACE AS ZERO.
C
      NCP1       = NC + 1
      DO 10 J    = 1, NCP1
         IW( J ) = 0
   10 CONTINUE
C
C  PASS 1. COUNT THE NUMBER OF ELEMENTS IN EACH COLUMN.
C
      DO 20 K   = 1, NNZ
        J       = JCN( K )
        IW( J ) = IW( J ) + 1
C        IP( J ) = IP( J ) + 1
   20 CONTINUE
C
C  PUT THE POSITIONS WHERE EACH COLUMN BEGINS IN
C  A COMPRESSED COLLECTION INTO IP AND IW.
C
      IP( 1 )       = 1
      DO 30 J       = 2, NCP1
        IP( J )     = IW( J - 1 ) + IP( J - 1 )
        IW( J - 1 ) = IP( J - 1 )
   30 CONTINUE
C
C  PASS 2. REORDER THE ELEMENTS INTO COLUMN ORDER. 
C          FILL IN EACH COLUMN IN TURN.
C
      DO 70 IC = 1, NC
C
C  CONSIDER THE NEXT UNFILLED POSITION IN COLUMN IC.
C
        DO 60 K = IW( IC ), IP( IC + 1 ) - 1
C
C  THE ENTRY SHOULD BE PLACED IN COLUMN J.
C
          I       = IRN( K )
          J       = JCN( K )
          ANEXT   = A( K )
          DO 40 L = 1, NNZ
C
C  SEE IF THE ENTRY IS ALREADY IN PLACE.
C
             IF ( J .EQ. IC ) GO TO 50
             LOCAT = IW( J )
C          
C  AS A NEW ENTRY IS PLACED IN COLUMN J, INCREASE THE POINTER 
C  IW( J ) BY ONE.
C          
             IW( J  ) = LOCAT + 1
C
C  RECORD DETAILS OF THE ENTRY WHICH CURRENTLY OCCUPIES LOCATION LOCAT.
C
             ITEMP = IRN( LOCAT )
             JTEMP = JCN( LOCAT )
             ATEMP = A( LOCAT )
C
C  MOVE THE NEW ENTRY TO IT CORRECT PLACE. 
C
             IRN( LOCAT ) = I 
             JCN( LOCAT ) = J  
             A( LOCAT )   = ANEXT
C
C  MAKE THE DISPLACED ENTRY THE NEW ENTRY.
C
             I          = ITEMP
             J          = JTEMP
             ANEXT      = ATEMP
   40     CONTINUE
C
C  MOVE THE NEW ENTRY TO IT CORRECT PLACE. 
C
   50     CONTINUE
          JCN( K ) = J
          IRN( K ) = I
          A( K )   = ANEXT
   60   CONTINUE
   70 CONTINUE
      RETURN
C
C  END OF REORDA.
C
      END
