      FUNCTION ran1(init)
      IMPLICIT NONE
      INTEGER init,idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
C     "Minimal" random number generator of Park and Miller with Bays-Durham
C     shuffle and added safeguards. Returns a uniform random deviate between
C     0.0 and 1.0 (exclusive of the endpoint values). Call with init a
C     negative integer to initialize; thereafter, idum is not altered between
C     successive deviates in a sequence. RNMX should approximate the largest
C     floating value that is less than 1.
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy,idum
      DATA iv /NTAB*0/, iy /0/
      if (init.le.0.or.iy.eq.0) then
            init=max(-init,1)
            do j=NTAB+8,1,-1
                  k=init/IQ
                  init=IA*(init-k*IQ)-IR*k
                  if (init.lt.0) init=init+IM
                  if (j.le.NTAB) iv(j)=init
            enddo
            iy=iv(1)
            idum=init
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END

c       PROGRAM rantest
c       IMPLICIT NONE
c       INTEGER idnum, i
c       REAL ran1
c       idnum = -42
c       print *, ran1(idnum)
c       do i=1,10
c             print *, ran1(1)
c       enddo
c       END PROGRAM rantest
