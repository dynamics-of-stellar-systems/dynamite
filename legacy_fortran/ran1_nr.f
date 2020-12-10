      DOUBLE PRECISION FUNCTION ran1(init)
      IMPLICIT NONE
      INTEGER init,idum,IA,IM,IQ,IR,NTAB,NDIV,mode
C      DOUBLE PRECISION ran1,AM,EPS,RNMX
      DOUBLE PRECISION AM,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1d0/IM,IQ=127773,IR=2836,
C     *           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1e0-EPS)
     *           NTAB=32,NDIV=1+(IM-1)/NTAB)
C     Dynamite random number generator. First call:
C     init .le. 0 -> initialize random sequence
C     init .gt. 0 -> return next "random" number 0 < ran1 < 1
C     "Minimal" random number generator of Park and Miller with Bays-Durham
C     shuffle and added safeguards. Returns a uniform random deviate between
C     0.0 and 1.0 (exclusive of the endpoint values). Call with init a
C     negative integer to initialize; thereafter, idum is not altered between
C     successive deviates in a sequence. RNMX should approximate the largest
C     floating value that is less than 1.
C     Function has been adapted to double precision, otherwise taken from
C     Numerical Rexipes 2nd ed.
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy,idum
      DATA iv /NTAB*0/, iy /0/

      RNMX = 1d0-2.23d-16 ! 2.23d-16 is approx. epsilon(1d0)
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
