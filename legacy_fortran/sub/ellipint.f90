module ellipticalintegrals
! Code written by VdVen. From numerical recipes (?).
use numeric_kinds
implicit none
private

  private :: rf, rd ! Carslon circular and elliptic functions
  public :: ellf, elle ! Legendre ellptic integrals

contains

  !!! Carlson's function R_F(x,y,z) !!!
  SUBROUTINE rf(x,y,z,res)
    REAL (kind=dp), intent(in) :: x,y,z
    REAL (kind=dp), intent(out) :: res
    REAL (kind=dp), parameter :: ERRTOL=0.08,TINY=1.5e-38,BIG=3.E37,THIRD=1.0/3.0,&
         C1=1.0/24.0,C2=0.1,C3=3.0/44.0,C4=1.0/14.0  
    REAL (kind=dp) :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt 
    if(min(x,y,z) < 0.0.or.min(x+y,x+z,y+z) < TINY.or.max(x,y,z) > BIG)&
         stop "invalid arguments in rf"
    xt=x
    yt=y
    zt=z
    do 
    sqrtx=sqrt(xt)
    sqrty=sqrt(yt)
    sqrtz=sqrt(zt)
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    xt=0.25*(xt+alamb)
    yt=0.25*(yt+alamb)
    zt=0.25*(zt+alamb)
    ave=THIRD*(xt+yt+zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    if(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
    end do
    e2=delx*dely-delz**2
    e3=delx*dely*delz
    res=(1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
  END SUBROUTINE rf

  !!! Carlson's function R_D(x,y,z) !!!
  SUBROUTINE rd(x,y,z,res)
    REAL (kind=dp), intent(in) :: x,y,z
    REAL (kind=dp), intent(out) :: res
    REAL (kind=dp), parameter :: ERRTOL=0.05,TINY=1.5e-25,BIG=4.5E21,&
         C1=3.0/14.0,C2=1.0/6.0,C3=9.0/22.0,C4=3.0/26.0,C5=0.25*C3,C6=1.5*C4 
    REAL (kind=dp) :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,&
         sqrtx,sqrty,sqrtz,sum,xt,yt,zt 
    if(min(x,y) < 0.0.or.min(x+y,z) < TINY.or.max(x,y,z) > BIG)&
         stop "invalid arguments in rd"
    xt=x
    yt=y
    zt=z
    sum=0.0
    fac=1.0
    do
       sqrtx=sqrt(xt)
       sqrty=sqrt(yt)
       sqrtz=sqrt(zt)
       alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
       sum=sum+fac/(sqrtz*(zt+alamb))
       fac=0.25*fac
       xt=0.25*(xt+alamb)
       yt=0.25*(yt+alamb)
       zt=0.25*(zt+alamb)
       ave=0.2*(xt+yt+3.0*zt)
       delx=(ave-xt)/ave
       dely=(ave-yt)/ave
       delz=(ave-zt)/ave
       if(max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
    end do
    ea=delx*dely
    eb=delz*delz
    ec=ea-eb
    ed=ea-6.0*eb
    ee=ed+ec+ec
    res=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*&
         ec+delz*C4*ea)))/(ave*sqrt(ave))
  END SUBROUTINE rd


  SUBROUTINE ellf(pp,kk,res)
    REAL (kind=dp), intent(in) :: pp,kk
    REAL (kind=dp), intent(out) :: res
    REAL (kind=dp) :: s,dum
    s=sin(pp)
    call rf(cos(pp)**2,(1.0-s*kk)*(1.0+s*kk),1.0_dp,dum) 
    res=s*dum
  END SUBROUTINE ellf

  !!! elliptic integral of the second kind !!!
  !!! Def.: E(phi,k) = int(sqrt(1-k^2*sin^2theta),0..phi) !!!
  SUBROUTINE elle(pp,kk,res)
    REAL (kind=dp), intent(in) :: kk,pp
    REAL (kind=dp), intent(out) :: res
    REAL (kind=dp) :: s,cc,q,dum1,dum2
    s=sin(pp)
    cc=cos(pp)**2
    q=(1.0-s*kk)*(1.0+s*kk)
    call rf(cc,q,1.0_dp,dum1)
    call rd(cc,q,1.0_dp,dum2)
    res=s*(dum1-((s*kk)**2)*dum2/3.0)
  END SUBROUTINE elle
 
end module ellipticalintegrals