!Josh Adams combined these routines to solve the incomplete Beta
!function for the calculation of the gNFW potential. The integral 
!needing solution is from Zhao'96, Eq A2 as B(a,b,x)=
!int(dt*t^(a-1)*(1-t)^(b-1),0,x). I have taken routines from numrec, 
!but removed some unnecessary normalization. We also want the complete 
!beta function, which can be simply related to the gamma function. 
!To avoid conflicts with the 
!original routines, I've prefixed zh_ to all the routines herein.
      FUNCTION zh_betai(a,b,x)
      use numeric_kinds
      Real(kind=dp):: zh_betai,a,b,x,zh_beta
!U    USES betacf,gammln
      Real(kind=dp):: bt,zh_betacf,zh_gammln
      if(x.lt.0.0_dp.or.x.gt.1.0_dp) print*, 'bad argument x in zh_betai'
      if(x.eq.0.0_dp)then
        bt=0.0_dp
      else if(x.eq.1.0_dp)then
        bt=zh_beta(a,b)
        return
      else
        bt=x**a*(1.0_dp-x)**b
      endif
      if(x.lt.(a+1.0_dp)/(a+b+2.0_dp).or.b.le.0.0_dp)then
!Can't use the inverse relation with b=0 b/c beta_comp will be infinite.
        zh_betai=bt*zh_betacf(a,b,x)/a
        return
      else
        zh_betai=zh_beta(a,b)-bt*zh_betacf(b,a,1.0_dp-x)/b
        return
      endif
      END
!---------------------------------------------------------------
      FUNCTION zh_beta(z,w)
      use numeric_kinds
      Real(kind=dp):: zh_beta,w,z
!U    USES gammln
      Real(kind=dp):: zh_gammln
      zh_beta=exp(zh_gammln(z)+zh_gammln(w)-zh_gammln(z+w))
      return
      END
!---------------------------------------------------------------
      FUNCTION zh_betacf(a,b,x)
      use numeric_kinds
      INTEGER(kind=i4b):: MAXIT
      Real(kind=dp):: zh_betacf,a,b,x,EPS,FPMIN
      PARAMETER (MAXIT=500,EPS=3.e-7_dp,FPMIN=1.e-30_dp)
      INTEGER(kind=i4b):: m,m2
      Real(kind=dp):: aa,c,d,del,h,qab,qam,qap
      qab=a+b
      qap=a+1.0_dp
      qam=a-1.0_dp
      c=1.0_dp
      d=1.0_dp-qab*x/qap
      if(abs(d).lt.FPMIN)d=FPMIN
      d=1.0_dp/d
      h=d
      do 11 m=1,MAXIT
        m2=2*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0_dp+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.0_dp+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.0_dp/d
        h=h*d*c
        aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.0_dp+aa*d
        if(abs(d).lt.FPMIN)d=FPMIN
        c=1.0_dp+aa/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp).lt.EPS)goto 1
11    continue
      print*, 'a or b too big, or MAXIT too small in zh_betacf', a, b, &
x, abs(del-1.0_dp)
1     zh_betacf=h
      return
      END
!---------------------------------------------------------------
      FUNCTION zh_gammln(xx)
      use numeric_kinds
      Real(kind=dp):: zh_gammln,xx
      INTEGER(kind=i4b):: j
      Real(kind=dp):: ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146_dp,-86.50532032941677_dp,&
      24.01409824083091_dp,-1.231739572450155_dp,.1208650973866179e-2_dp,&
      -.5395239384953e-5_dp,2.5066282746310005_dp/
      x=xx
      y=x
      tmp=x+5.5_dp
      tmp=(x+0.5_dp)*log(tmp)-tmp
      ser=1.000000000190015_dp
      do 11 j=1,6
        y=y+1.0_dp
        ser=ser+cof(j)/y
11    continue
      zh_gammln=tmp+log(stp*ser/x)
      return
      END
!---------------------------------------------------------------
