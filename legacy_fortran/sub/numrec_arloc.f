      SUBROUTINE HUNTPOS(xx,n,x,jbin)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Find the bin jbin in which to store a datapoint x, given an array
C xx(1:n). Based on the Numerical Recipes Routine hunt. At input jbin must
C be an initial guess.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      REAL*8 x,xx(n)
C
C Find the point jbin, such that x is between xx(jbin) and xx(jbin+1)
C
      Stop " functuon huntpos used !"
      CALL HUNT(xx,n,x,jbin)
C
C Find the correct bin
C
      IF (jbin.EQ.0) THEN
        jbin = 1
      ELSE IF (jbin.NE.n) THEN
        IF (ABS(x-xx(jbin)).GT.(0.5D0*ABS(xx(jbin+1)-xx(jbin)))) THEN
          jbin = jbin+1
        END IF
      END IF
C
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Numerical Recipes routine to locate a position in an array.
C Transformed to double precision.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE hunt(xx,n,x,jlo)
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END

