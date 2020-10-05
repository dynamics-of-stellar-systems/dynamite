      SUBROUTINE GETGAUHER (vm,sg,veltemp,Nvhist,Nvmax,dvhist,
     &                      hh,Nhermmax,gam,ingam)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate the first Nhermmax Gauss-Hermite moments from the histogram
C veltemp with spacing dvhist, using the values vm and sg for the weighting
C function. The results are returned in hh.
C
C If ingam=1, then gam is taken as input.
C If ingam=0, then gam is returned on output such that h0=1.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
C
      PARAMETER (pi = 3.14159265358979323846D0)
C
      DIMENSION veltemp(-Nvmax:Nvmax),hh(0:Nhermmax)
C
CCCCCCCCCCCCCCCCCCCC
C
C Initialize
C
 
      DO l=0,Nhermmax
        hh(l) = 0.0D0
      END DO
C
C Loop over the velocities
C
      DO i=-Nvhist,Nvhist
C
        vel = dvhist * DBLE(i)
        w   = (vel-vm)/sg
C
C Loop over the Gauss-Hermite moments
C
        DO l=0,Nhermmax
          hh(l) = hh(l) + (veltemp(i)*SDGAUSS(w)*H_POL(l,w))
        END DO
C
      END DO
C
C Normalize properly. The factor dvhist arises through the stepsize of the
C Euler integration. The value of gamma is determined by the constraint
C that h0=1.
C
      IF (ingam.EQ.0) THEN
        gam = hh(0) * 2.0D0 * SQRT(pi) * dvhist
      END IF
C
      DO l=0,Nhermmax
        hh(l) = hh(l) * 2.0D0 * SQRT(pi) * dvhist / gam
      END DO
C
      END


      SUBROUTINE GETGAUSSWEIGHTMOM
     &        (vm,sg,veltemp,Nvhist,Nvmax,dvhist,gw,Nmommax)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate the first Nmommax Gaussian weighted moments from the histogram
C veltemp with spacing dvhist, using the values vm and sg for the weighting
C function. The results are returned in hh.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
C
      PARAMETER (pi = 3.14159265358979323846D0)
C
      DIMENSION veltemp(-Nvmax:Nvmax),gw(0:Nmommax)
C
CCCCCCCCCCCCCCCCCCCC
C
C Initialize
C
      DO l=0,Nmommax
        gw(l) = 0.0D0
      END DO
C
C Loop over the velocities
C
      DO i=-Nvhist,Nvhist
C
        vel = dvhist * DBLE(i)
        w   = (vel-vm)/sg
C
C Loop over the moments
C
        DO l=0,Nmommax
          gw(l) = gw(l) + (veltemp(i)*SDGAUSS(w)*(vel**l))
        END DO
C
      END DO
C
C Normalize properly
C
      DO l=0,Nmommax
        gw(l) = gw(l) * dvhist
      END DO
C
      END


      SUBROUTINE GETTRUEMOM (veltemp,Nvhist,Nvmax,dvhist,
     &                       allmom,Nmommax)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate the first Nmommax moments from the histogram
C veltemp with spacing dvhist. The results are returned in allmom.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
C
      DIMENSION veltemp(-Nvmax:Nvmax),allmom(0:Nmommax)
C
CCCCCCCCCCCCCCCCCCCC
C
C Initialize
C
      DO l=0,Nmommax
        allmom(l) = 0.0D0
      END DO
C
C Loop over the velocities
C
      DO i=-Nvhist,Nvhist
C
        vel = dvhist * DBLE(i)
C
C Loop over the moments
C
        DO l=0,Nmommax
          allmom(l) = allmom(l) + (veltemp(i)*(vel**l))
        END DO
C
      END DO
C
C Normalize properly
C
      DO l=0,Nmommax
        allmom(l) = allmom(l) * dvhist
      END DO
C
      END


      SUBROUTINE GAUSSFIT (veltemp,Nvhist,Nvmax,dvhist,
     &                              gam,vm,sg,h12,NMC)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate the normalization gam, mean vm and dispersion sg of
C the best-fitting Gaussian to the histogram veltemp with spacing dvhist.
C
C A good initial guess is sought using a random number scheme, by drawing
C NMC (vm,sg) points.
C
C The program uses an amoeba scheme which should be very robust.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
C
      PARAMETER (eps = 1.0D-10    ,
     &           Nhistmax = 10000 )
C
      DIMENSION veltemp(-Nvmax:Nvmax),vvtemp(-Nhistmax:Nhistmax),
     &          allmom(0:2)
C
C Starting simpleces for routine AMOEBA, and the function values
C in the starting simplices.
C
      DIMENSION P(3,2),Y(3),help(2)
C
      COMMON /vpcur/ vvtemp,dvhistt,gamma,Nvhistt
      COMMON /truesigcur/ sg0
C
      EXTERNAL CHI2H_2D
C
CCCCCCCCCCCCCCCCCCCC
C
C Store the velocity profile in a common block
C
      Nvhistt = Nvhist
      dvhistt = dvhist
      DO i=-Nhistmax,Nhistmax
        IF (ABS(i).LE.Nvhist) THEN
          vvtemp(i) = veltemp(i)
        ELSE
          vvtemp(i) = 0.0D0
        END IF
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Get the lowest order true moments
C
      CALL GETTRUEMOM (veltemp,Nvhist,Nvmax,dvhist,allmom,2)
C
C Initialize the estimates for the mean and dispersion of the best
C Gaussian fit.
C
      vm0 = allmom(1)/allmom(0)
      sg0 = SQRT((allmom(2)/allmom(0))-(vm0*vm0))
C
CCCCCCCCCCCCCCCCCCCC
C
C Draw a certain number of random points in (vm,sg). Remember
C the one with the lowest CHI2H_2D.
C
      Ndraw   = NMC
      idum    = -100
      vm      = vm0
      sg      = sg0
      help(1) = vm
      help(2) = sg
      chi2min = CHI2H_2D(help)
C
51    DO i=1,Ndraw
        help(1) = vm0 + (2.0D0 * sg0 * ((2.0D0*randraw(idum))-1.0D0))
        help(2) = 2.0D0 * sg0 * randraw(idum)
        chi2cur = CHI2H_2D(help)
        IF (chi2cur.LT.chi2min) THEN
          vm = help(1)
          sg = help(2)
          chi2min = chi2cur
        END IF
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Initialize the starting simplex fits
C
      epsmal = 0.2D0
      pl1    = 1.0D0 + epsmal
      xmn1   = 1.0D0 - epsmal
C
      P(1,1) = vm + (epsmal*sg)
      P(1,2) = pl1 * sg
      P(2,1) = vm - (epsmal*sg)
      P(2,2) = xmn1 * sg
      P(3,1) = vm + (epsmal*sg)
      P(3,2) = xmn1 * sg
C
C Initialize
C
      DO i=1,3
        DO j=1,2
          help(j) = P(i,j)
        END DO
        Y(i) = CHI2H_2D(help)
      END DO
C
      CALL AMOEBA(P,Y,3,2,2,eps,CHI2H_2D,iter)
C
      vm  = (P(1,1)+P(2,1)+P(3,1)) / 3.0D0
      sg  = (ABS(P(1,2))+ABS(P(2,2))+ABS(P(3,2))) / 3.0D0
C
      help(1) = vm
      help(2) = sg
      h12 = CHI2H_2D(help) - 1.0D0
      gam = gamma
C
C Start from the beginning if the found minimum does not satisfy
C h1=h2=0 to high enough accuracy.
C
      IF ((h12.GT.100.0D0*eps).AND.(Ndraw.LE.3000)) THEN
        Ndraw   = Ndraw*5
        chi2min = h12 + 1.0D0
        WRITE (*,*) 'Redraw', Ndraw
        GOTO 51
      END IF
C
      END




      REAL*8 FUNCTION CHI2H_2D (y)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculates the chih^2 = (h1^2) + (h2^2) for a Gaussian
C with parameters Vgau = y(1), sig = |y(2)|,
C for the VP in the common block /vpcur/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
C
      PARAMETER (Nhistmax = 10000)
C
      DIMENSION y(2),harr(0:10)
C
      DIMENSION veltemp(-Nhistmax:Nhistmax)
      COMMON /vpcur/ veltemp,dvhist,gam,Nvhist
      COMMON /truesigcur/ sg0
C
CCCCCCCCCCCCCCCCCCCC
C
C Note: it may be necessary to avoid values of gamma and sigma too
C close to zero.
C
      vm  = y(1)
      sg  = MAX(0.1D0*sg0,ABS(y(2)))
C
      CALL GETGAUHER (vm,sg,veltemp,Nvhist,Nhistmax,dvhist,
     &                harr,2,gam,0)
C
      CHI2H_2D = 1.0D0 + SQRT(MAX(0.0D0,
     &         (harr(1)**2.0D0)+(harr(2)**2.0D0) ))
C
      END


      SUBROUTINE GAUSSFITfix (veltemp,Nvhist,Nvmax,dvhist,
     &                                  vm,sg,h12,NMC)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculate the normalization mean vm and dispersion sg of
C the best-fitting NORMALIZED Gaussian to the histogram veltemp with
C spacing dvhist.
C
C A good initial guess is sought using a random number scheme, by drawing
C NMC (vm,sg) points.
C
C The program uses an amoeba scheme which should be very robust.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
C
      PARAMETER (eps = 1.0D-10    ,
     &           Nhistmax = 10000 )
C
      DIMENSION veltemp(-Nvmax:Nvmax),vvtemp(-Nhistmax:Nhistmax),
     &          allmom(0:2)
C
C Starting simpleces for routine AMOEBA, and the function values
C in the starting simplices.
C
      DIMENSION P(3,2),Y(3),help(2)
C
      COMMON /vpcur_FIX/ vvtemp,dvhistt,Nvhistt
      COMMON /truesigcur_FIX/ sg0
C
      EXTERNAL CHI2H_2D_FIX
C
CCCCCCCCCCCCCCCCCCCC
C
C Store the velocity profile in a common block
C
      Nvhistt = Nvhist
      dvhistt = dvhist
      DO i=-Nhistmax,Nhistmax
        IF (ABS(i).LE.Nvhist) THEN
          vvtemp(i) = veltemp(i)
        ELSE
          vvtemp(i) = 0.0D0
        END IF
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Get the lowest order true moments
C
      CALL GETTRUEMOM (veltemp,Nvhist,Nvmax,dvhist,allmom,2)
C
C Initialize the estimates for the mean and dispersion of the best
C Gaussian fit.
C
      vm0 = allmom(1)/allmom(0)
      sg0 = SQRT((allmom(2)/allmom(0))-(vm0*vm0))
C
CCCCCCCCCCCCCCCCCCCC
C
C Draw a certain number of random points in (vm,sg). Remember
C the one with the lowest CHI2H_2D.
C
      Ndraw   = NMC
      idum    = -100
      vm      = vm0
      sg      = sg0
      help(1) = vm
      help(2) = sg
      chi2min = CHI2H_2D_FIX(help)
C
51    DO i=1,Ndraw
        help(1) = vm0 + (2.0D0 * sg0 * ((2.0D0*randraw(idum))-1.0D0))
        help(2) = 2.0D0 * sg0 * randraw(idum)
        chi2cur = CHI2H_2D_FIX(help)
        IF (chi2cur.LT.chi2min) THEN
          vm = help(1)
          sg = help(2)
          chi2min = chi2cur
        END IF
      END DO
C
CCCCCCCCCCCCCCCCCCCC
C
C Initialize the starting simplex fits
C
      epsmal = 0.2D0
      pl1    = 1.0D0 + epsmal
      xmn1   = 1.0D0 - epsmal
C
      P(1,1) = vm + (epsmal*sg)
      P(1,2) = pl1 * sg
      P(2,1) = vm - (epsmal*sg)
      P(2,2) = xmn1 * sg
      P(3,1) = vm + (epsmal*sg)
      P(3,2) = xmn1 * sg
C
C Initialize
C
      DO i=1,3
        DO j=1,2
          help(j) = P(i,j)
        END DO
        Y(i) = CHI2H_2D_FIX(help)
      END DO
C
      CALL AMOEBA(P,Y,3,2,2,eps,CHI2H_2D_FIX,iter)
C
      vm  = (P(1,1)+P(2,1)+P(3,1)) / 3.0D0
      sg  = (ABS(P(1,2))+ABS(P(2,2))+ABS(P(3,2))) / 3.0D0
C
      help(1) = vm
      help(2) = sg
      h12 = CHI2H_2D_FIX(help) - 1.0D0
C
C Start from the beginning if the found minimum does not satisfy
C h12 = 0 to high enough accuracy.
C
      IF ((h12.GT.100.0D0*eps).AND.(Ndraw.LE.3000)) THEN
        Ndraw   = Ndraw*5
        chi2min = h12 + 1.0D0
        WRITE (*,*) 'Redraw', Ndraw
        GOTO 51
      END IF
C
      END


      REAL*8 FUNCTION CHI2H_2D_FIX (y)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Calculates the chih^2 = (h1^2) + ((h0-SQRT(2)*h2)^2) for a Gaussian
C with parameters Vgau = y(1), sig = |y(2)|,
C for the VP in the common block /vpcur_FIX/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
C
      PARAMETER (Nhistmax = 10000)
C
      DIMENSION y(2),harr(0:10)
C
      DIMENSION veltemp(-Nhistmax:Nhistmax)
      COMMON /vpcur_FIX/ veltemp,dvhist,Nvhist
      COMMON /truesigcur_FIX/ sg0
C
CCCCCCCCCCCCCCCCCCCC
C
      vm  = y(1)
      sg  = MAX(0.1D0*sg0,ABS(y(2)))
C
      CALL GETGAUHER (vm,sg,veltemp,Nvhist,Nhistmax,dvhist,
     &                harr,2,1.0D0,1)
C
      CHI2H_2D_FIX = 1.0D0 + SQRT(MAX(0.0D0,
     &      (harr(1)**2.0D0) +
     &      ((1.0D0-harr(0)+(SQRT(2.0D0)*harr(2)))**2.0D0) ))
C
      END


      REAL*8 FUNCTION SDGAUSS (x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Returns the standard Gaussian as function of x
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
      PARAMETER (pi=3.14159265358979D0)
      SDGAUSS = (1.0D0/SQRT(2.0D0*pi)) * EXPP(-0.5D0*x*x)
      END


      REAL*8 FUNCTION EXPP(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Modified exponential function that will not underflow
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IF (x.GE.-500D0) THEN
        EXPP = EXP(x)
      ELSE
        EXPP = 0.0D0
      END IF
      END


      REAL*8 FUNCTION HE_POL (ll,x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Returns the value of the Hermite polynomial He_l(x) as defined
C in Appendix A of van der Marel & Franx.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
C
      LOGICAL firstc
      DATA firstc/.TRUE./
C
      DIMENSION pfac(0:20,0:20)
      SAVE firstc,pfac
C
CCCCCCCCCCCCCCCCCCCC
C
      IF (firstc) THEN
        DO l=0,20
          dl = DBLE(l)
          DO j=l,0,-1
            IF (2*((j+l)/2).EQ.j+l) THEN
              dj = DBLE(j)
              pfacln = (0.5D0*gammaln(dl+1.0D0)) -
     &             gammaln(dj+1.0D0) - gammaln(dl-dj+1.0D0) +
     &             gammaln(0.5D0*(dl-dj+1.0D0)) -
     &             gammaln(0.5D0) + (0.5D0*(dl-dj)*LOG(2.0D0))
              pfac(l,j) = ((-1.0D0)**((l-j)/2)) * EXPP(pfacln)
            ELSE
              pfac(l,j) = 0.0D0
            END IF
          END DO
        END DO
        firstc = .FALSE.
      END IF
C
CCCCCCCCCCCCCCCCCCCC
C
      HE_POL = 0.0D0
      DO j=ll,0,-1
        HE_POL = pfac(ll,j) + (x*HE_POL)
      END DO
C
      END


      REAL*8 FUNCTION H_POL (l,x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Returns the value of the Hermite polynomial H_l(x) as defined
C in Appendix A of van der Marel & Franx.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n) 
      H_POL = HE_POL(l,x*SQRT(2.0D0))
      END


      REAL*8 FUNCTION gammaln(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The logarithm of the gamma function
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      PARAMETER (pi=3.14159265358979D0)
      IF (x.GE.1.0D0) THEN
        gammaln = gammln(x)
      ELSE IF (x.GE.0.0D0) THEN
        z = 1.0D0 - x
        gammaln = (LOG((pi*z)/SIN(pi*z))) - gammln(2.0D0-x)
      ELSE
        PAUSE 'x < 0 in gammaln'
      END IF
      END


      REAL*8 FUNCTION gammln(x)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C From Numerical recipes. Modified to be double precision
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER j
      DIMENSION cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END


      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C From Numerical Recipes. Modified to be double precision
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER*4 iter,mp,ndim,np,NMAX,ITMAX
      REAL*8 ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=1000)
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER*4 i,ihi,ilo,inhi,j,m,n
      REAL*8 rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1     do 12 n=1,ndim
        sum=0.0D0
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
C
C Possibility to write intermediate results
C
C      WRITE (*,'(I5,3F20.12)') iter,(p(ilo,j),j=1,2),y(ilo)
C
      if (iter.ge.ITMAX) then
        write (*,*) 'ITMAX exceeded in amoeba'
        return
      endif
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0D0)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0D0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5D0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5D0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END


      REAL*8 FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C From Numerical Recipes. Modified to be double precision
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER*4 ihi,mp,ndim,np,NMAX
      REAL*8 fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
      REAL*8 fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.0D0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry)
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END


      FUNCTION randraw(idum)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Function to compute a sequence of pseudo-random numbers in the
C   interval (0,1). 'Minimal' random number generator of Park and Miller
C   with Bays-Durham shuffle and added safeguards. idum must be negative
C   and must not be alterated on successive calls. (See Numerical
C   Recipes for more details, 2nd. Edition, it corresponds to ran). It
C   is not useful for a number of calls > 10^8.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Declaration of local variables.
C
      INTEGER*4 idum, ia,im,iq,ir,ntab,ndiv
      REAL*8  randraw,am,epsran,rnmx
      PARAMETER (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,
     &           ntab=32,ndiv=1+(im-1)/ntab,epsran=1.2e-7,
     &           rnmx=1.0-epsran)
      INTEGER j,k,iv(ntab),iy
C
      SAVE iv,iy
      DATA iv /ntab*0/, iy /0/
C
      IF (idum.LE.0.or.iy.EQ.0) THEN
        idum=max(-idum, 1)
        DO j=ntab+8, 1, -1
          k=idum/iq
          idum=ia*(idum-k*iq)-ir*k
          IF (idum.LT.0) idum=idum+im
          IF (j.LE.ntab) iv(j)=idum
        END DO
        iy=iv(1)
      END IF
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      IF (idum.LT.0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      randraw=MIN(am*iy,rnmx)
C
      RETURN
      END
