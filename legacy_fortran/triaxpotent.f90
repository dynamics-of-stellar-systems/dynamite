!#######################################################################
!
! Calculate the Triaxial Potential from an MGE fit.
!
! The potential is calculated seperately for the inner, middle and outer 
! regions.
!
! This Fortran 90 code is compatible with the F language subset
! http://www.fortran.com/imagine1
!
! USAGE:
!
!  Initial Setup : Call tp_setup ()
!  Potential     : Call tp_potent(x,y,z,potential)
!  Acceleration  : Call tp_accel (x,y,z,vz,vy,vz )
!
!
! ToDo:
!
! - None
!  
! Bugs:
! 
!
! HISTORY:
!
! V1.0 Written By Remco van den Bosch JULY/2004
! V1.1 RvdB 24/AUG/2004:
!      Fixed temporairy qintr**2 = qintr**4
!      Added check for A1+A2+A3=sqrt(1-q^2)/(p*q).
!      Simplified the constant for the point mass approx.
!      Added additional braces for readability.
! V1.2 Added ellf and elle into this file for F compatibility
! V1.3 Fixed math bug in the accin. 
!      Checked overlap regions
!      Added overlap region check on runtime
! V1.4 Removed the black hole and dark halo and moved them to dmpotent.f90
!
!######################################################################

module triaxpotent
  use numeric_kinds
  use ellipticalintegrals
  use initial_parameters
  implicit none
  private

  ! calculate potential \phi at (x,y,z) 
  public:: tp_potent
  
  ! calculate accel ax,ay,ax at (x,y,z)
  public:: tp_accel
  
  ! setup the constants for the potential
  public:: tp_setup
  
  ! The potentential estimators for the internal and middel part of the pot.
  private:: potin, potmid
  
  ! The acceleration estimators for the internal and middel part of the pot.
  private:: accmid,accin
  
  ! Integrator help functions
  private:: potfunc,axfunc,ayfunc,azfunc
  
  ! constants for the central approximation
  real (kind=dp),private,dimension(:),allocatable :: A1,A2,A3,F,V0
  
  ! deprojected information about the gaussians
  real (kind=dp),public ,dimension(:),allocatable :: sigintr_km, pintr
  real (kind=dp),public ,dimension(:),allocatable :: qintr, dens
  real (kind=dp),public ,dimension(:),allocatable :: triaxpar
  ! global helper variables for the integration
  real   (kind=dp ),private:: gx,gy,gz
  integer(kind=i4b),private:: gn  ! global gauss index number

  ! Boundary for the inner approximation (in fraction of gauss sigma)
  real (kind=dp),private,parameter :: inner_approx=0.0001
  ! Boundary for the inner approximation (in fraction of gauss sigma)
  real (kind=dp),private,parameter :: outer_approx=300.0

!######################################################################

contains

  !+++++++++++++++++++++++++++++++++++++
  subroutine tp_setup()
    use  ellipticalintegrals, only : elle,ellf
  !-------------------------------------
    integer(kind=i4b) :: i
    real   (kind=dp ) :: E,k,p,q,trmass,secth,cotph
    real   (kind=dp ) :: ix,iy,iz,ax,ay,az
    real   (kind=dp ),dimension(:),allocatable:: delp, &
                         nom1minq2, nomp2minq2, denom

    print*," Setting up triaxial potential routines"
    allocate ( delp      (ngauss_mge), nom1minq2(ngauss_mge), &
               nomp2minq2(ngauss_mge), denom    (ngauss_mge), & 
               sigintr_km(ngauss_mge), pintr    (ngauss_mge), &
               qintr     (ngauss_mge), triaxpar (ngauss_mge), &
               dens      (ngauss_mge), V0       (ngauss_mge))

    ! do the analytic triaxial deprojection under the MGE hypotesis
    ! (see e.g. Cappellari 2002)
    secth = 1.0_dp/cos(theta_view)
    cotph = 1.0_dp/tan(  phi_view)

    delp(:) = 1.0_dp - qobs(:)**2

    nom1minq2(:) = delp(:)*( 2.0_dp*cos(2.0_dp*psi_obs(:)) + &
         sin(2.0_dp*psi_obs(:))*(secth*cotph - cos(theta_view)*tan(phi_view)))

    nomp2minq2(:) = delp(:)*( 2.0_dp*cos(2.0_dp*psi_obs(:)) + &
         sin(2.0_dp*psi_obs(:))*(cos(theta_view)*cotph - secth*tan(phi_view)))

    denom(:)  = 2.0_dp*sin(theta_view)**2*( delp(:)*cos(psi_obs(:))*&
         (cos(psi_obs(:)) + secth*cotph*sin(psi_obs(:))) - 1.0_dp )

    ! These are temporary values of the squared intrinsic axial 
    ! ratios p^2 and q^2
    qintr(:) = (1.0_dp   - nom1minq2(:)/denom(:))
    pintr(:) = (qintr(:) + nomp2minq2(:)/denom(:))

    ! Quick check to see if we are not going to take the sqrt of 
    ! a negative number.
    if ( any(qintr(:) < 0.0_dp) .or. any (pintr(:) < 0.0_dp)) &
        stop "p^2 or q^2 is below 0."

    ! intrinsic axial ratios p and q 
    qintr(:) = sqrt(qintr(:))
    pintr(:) = sqrt(pintr(:))

    print*,"middle axis ratio P"
    print*,pintr(:)
    print*,"minor  axis ratio q"
    print*,qintr(:)

    ! triaxiality parameter T = (1-p^2)/(1-q^2)
    triaxpar(:)=(1.0_dp-pintr(:)**2)/(1.0_dp-qintr(:)**2)
    if ( any(triaxpar(:) < 0.0_dp) .or.  any(triaxpar(:) > 1.0_dp) ) &
            stop "No triaxial deprojection possible"
    
    print*,"triaxiality parameters:"
    print*,triaxpar(:)

    if ( any ( qintr(:) > pintr(:) )) stop "q>p"  
    if ( any ( pintr(:) > 1.0_dp   )) stop "p>1"

    ! intrinsic sigma (Cappellari 2002 eq 9.) 
    sigintr_km(:) = sigobs_km(:)*sqrt( qobs(:)/sqrt((pintr(:)*&
        cos(theta_view))**2+ (qintr(:)*sin(theta_view))**2*((pintr(:)*&
        cos(phi_view))**2 + sin(phi_view)**2) ) )

    print*,"Unitless Length of the projected major axis: U"
    print*,sigobs_km(:)/sigintr_km(:)

    ! density factor
    dens(:) = surf_km(:)*qobs(:)*sigobs_km(:)**2/&
             (sqrt(twopi_d)*pintr(:)*qintr(:)*sigintr_km(:)**3)

    ! Integration Constant
    V0(:) = 4.0_dp * Pi_d * grav_const_km * sigintr_km(:)**2 * pintr(:) &
         * qintr(:) * dens(:)
      
    ! Masses of the individual gausses.
    print*,'  * Masses of the individual gausses'
    print*, twopi_d*(surf_km(:)*qobs(:)*sigobs_km(:)**2)
    ! Total mass of the galaxy
    trmass = twopi_d*sum(surf_km(:)*qobs(:)*sigobs_km(:)**2)
    print*, "Total mass (Msun) projected:",trmass
    trmass = sqrt(twopi_d)**3*sum(dens(:)*pintr(:)*qintr(:)*sigintr_km(:)**3)
    print*, "Total mass (Msun) intrinsic:",trmass

    ! Setup constants for the computation of the potential for the inner region
    allocate(A1(ngauss_mge),A2(ngauss_mge),A3(ngauss_mge)&
            ,F(ngauss_mge))

    ! Compute the Constants for the inner approximation
    do i=1,ngauss_mge 
       p     = pintr(i)
       q     = qintr(i)
 
       ! Calculate the Elliptical integrals
       k     = sqrt ( (1.0_dp - p*p)/ (1.0_dp-q*q) )
       call ellf(acos(q),k,F(i))
       call elle(acos(q),k,E   )
       A1(i) = (F(i) - E) / (1.0_dp - p*p)

       A2(i) = ((1.0_dp-q*q)*E - (p*p-q*q) * F(i) - (q/p) *(1.0_dp-p*p) * &
               sqrt(1.0_dp-q*q))/((1.0_dp-p*p)*(p*p-q*q))

       A3(i) = ( (p/q) * sqrt(1.0_dp-q*q) - E ) / (p*p - q*q)

       ! According to Glenn a1+a2+a3 should be equal to sqrt(1-q**2)/(p*q)
       if ( abs ( (sqrt(1-q**2)/(p*q)) - A1(i) - A2(i) - A3(i) ) &
            > 1.0e-6) then 
          print*,"  * Failure to properly compute A1, A2 and A3 in tp_setup"
          print*,"  * gauss_n,A1,A2,A3: ",i,A1(i),A2(i),A3(i)
          print*,abs ( (sqrt(1-q**2)/(p*q)) - A1(i) - A2(i) - A3(i) )
          print*,p,q, (1.0_dp - p*p)
          print*, (F(i) - E),k,epsilon(1.0_dp)
          stop
       end if

    end do



    ! Test the accuracy of the approximation regimes
    print*,"  Testing the accuracy of the approximation regimes."

    do i=1,ngauss_mge 
       call potin  (i,inner_approx*sigintr_km(i),1.0_dp,1.0_dp,ax)
       call potmid (i,inner_approx*sigintr_km(i),1.0_dp,1.0_dp,ix)
       if (abs((ix-ax)/ix) > 1.0e-4) stop "Failed test 1"  
       call potin  (i,1.0_dp,1.0_dp,inner_approx*sigintr_km(i),ax)
       call potmid (i,1.0_dp,1.0_dp,inner_approx*sigintr_km(i),ix)
       if (abs((ix-ax)/ix) > 1.0e-4) stop "Failed test 2" 
       ax=sqrt(pi_d/2.0_dp) * V0(i)/outer_approx
       call potmid (i,outer_approx*sigintr_km(i),0.0_dp,0.0_dp,ix)
       if (abs((ix-ax)/ix) > 1.0e-4) stop "Failed test 3" 
       ax=sqrt(pi_d/2.0_dp) * V0(i)/outer_approx
       call potmid (i,1.0_dp,1.0_dp,outer_approx*sigintr_km(i),ix)
       if (abs((ix-ax)/ix) > 1.0e-4) stop "Failed test 4" 
       call accin  (i,inner_approx*sigintr_km(i)*0.95,0.2*inner_approx*sigintr_km(i),0.2*inner_approx*sigintr_km(i),ax,ay,az)
       call accmid (i,inner_approx*sigintr_km(i)*0.95,0.2*inner_approx*sigintr_km(i),0.2*inner_approx*sigintr_km(i),ix,iy,iz)
       if ( sqrt( (ix-ax)**2+(iy-ay)**2+(iz-az)**2)/sqrt(ix**2+iy**2+iz**2)>1.0e-3) stop "Failed test 5"
       call accin  (i,0.2*inner_approx*sigintr_km(i),0.2*inner_approx*sigintr_km(i),0.95*inner_approx*sigintr_km(i),ax,ay,az)
       call accmid (i,0.2*inner_approx*sigintr_km(i),0.2*inner_approx*sigintr_km(i),0.95*inner_approx*sigintr_km(i),ix,iy,iz)
       if ( sqrt( (ix-ax)**2+(iy-ay)**2+(iz-az)**2)/sqrt(ix**2+iy**2+iz**2)>1.0e-2) stop "Failed test 6"
    end do

    print*,"  Triaxial potential routines setup finished"

  end subroutine tp_setup


!+++++++++++++++++++++++++++++++++++
subroutine tp_potent(x,y,z,pot)
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(out):: pot
!------------------------------------
  real(kind=dp) :: d2,tpot
  integer(kind=i4b) :: i

  pot = 0.0_dp
  
  d2  = x*x + y*y + z*z
  
  do i=1,ngauss_mge
  
     if      (d2 < ( inner_approx * sigintr_km(i))**2) then
        ! central approximation
        call potin (i,x,y,z,tpot)
      else if (d2 < (outer_approx * sigintr_km(i))**2) then
        ! numerically compute the full integral
        call potmid(i,x,y,z,tpot)
      else
        ! point mass approximation:
        ! tpot= sqrt(2.0_dp*pi_d)**3 *sigintr_km(i) &
        !       /(4.0_dp*pi_d) * V0(i)/sqrt(d2)
        tpot= sqrt(pi_d/2.0_dp) * sigintr_km(i) * V0(i)/sqrt(d2)
      endif
     
     pot = pot + tpot
  end do


end subroutine tp_potent

!+++++++++++++++++++++++++++++++++++
subroutine tp_accel(x,y,z,vx,vy,vz)
  real(kind=dp), intent(in) ::  x, y, z
  real(kind=dp), intent(out):: vx,vy,vz
!------------------------------------
  real(kind=dp) :: d2,tx,ty,tz,t,t1,t2,t3
  integer(kind=i4b) :: i

  vx = 0.0_dp
  vy = 0.0_dp
  vz = 0.0_dp
  
  d2 = x*x + y*y + z*z
  
  do i=1,ngauss_mge
  
     if      (d2 < ( inner_approx * sigintr_km(i))**2) then
        ! central approximation
        call accin (i,x,y,z,tx,ty,tz)
     else if (d2 < (outer_approx * sigintr_km(i))**2) then
        ! numerically compute the full integral.
        call accmid(i,x,y,z,tx,ty,tz)
     else
        ! point mass approximation:
        ! t  = -sqrt(2.0*pi_d)**3 *sigintr_km(i) /(4.0*pi_d) &
        !      *  V0(i)/d2**(3.0_dp/2.0)
        t  = - sqrt(pi_d/2.0_dp) * sigintr_km(i) * V0(i)/d2**(3.0_dp/2.0)
        tx = x * t
        ty = y * t
        tz = z * t
     endif
     vx = vx + tx
     vy = vy + ty
     vz = vz + tz

  enddo

end subroutine tp_accel

!++++++++++++++++++++++++++++++++++++
subroutine potmid(n,x,y,z,res)
  integer(kind=i4b),intent(in ) :: n
  real   (kind=dp ),intent(in ) :: x,y,z
  real   (kind=dp ),intent(out) :: res
!------------------------------------
  real   (kind=dp )             :: tpot
  integer(kind=i4b),parameter   :: limit=20,leniw=limit*3,lenw=limit*46
  integer(kind=i4b)             :: ier,last
  integer(kind=i4b),dimension(leniw) ::iwork
  real   (kind=dp )             :: epsabs,epsrel,abserr
  real   (kind=dp ),dimension(limit*46) :: work

  !set global integral variables
  gx = x
  gy = y
  gz = z
  gn = n
  
  epsabs = 0.0e-0_dp
  epsrel = 1.0e-8_dp
  call dqxgs(potfunc,0.0_dp,1.0_dp,epsabs,epsrel,tpot,&
       abserr,ier,limit,leniw,lenw,last,iwork,work)
  if (ier /= 0 ) print*,"dqxgs error in subroutine potin:",ier

  res = V0(n) * tpot
end subroutine potmid

!++++++++++++++++++++++++++++++++++++
subroutine accmid(n,x,y,z,vx,vy,vz)
  integer(kind=i4b),intent(in) :: n
  real   (kind=dp ),intent(in) ::  x, y, z
  real   (kind=dp ),intent(out):: vx,vy,vz
!------------------------------------
  real   (kind=dp )            :: tx,ty,tz
  integer(kind=i4b),parameter  :: limit=20,leniw=limit*3,lenw=limit*46
  integer(kind=i4b)            :: ier,last
  integer(kind=i4b),dimension(leniw) ::iwork
  real   (kind=dp )            :: epsabs,epsrel,abserr
  real   (kind=dp ),dimension(limit*46)  :: work

  !set global integral variables
  gx = x
  gy = y
  gz = z
  gn = n
  
  epsabs = 0.0_dp
  epsrel = 1.0e-8_dp
  call dqxgs(axfunc,0.0_dp,1.0_dp,epsabs,epsrel,tx,&
       abserr,ier,limit,leniw,lenw,last,iwork,work)
  if (ier /= 0 ) print*,"dqxgs error in subroutine potin:",ier
  vx = V0(n) * tx

  epsabs = 0.0_dp
  epsrel = 1.0e-8_dp
  call dqxgs(ayfunc,0.0_dp,1.0_dp,epsabs,epsrel,ty,&
       abserr,ier,limit,leniw,lenw,last,iwork,work)
  if (ier /= 0 ) print*,"dqxgs error in subroutine potin:",ier
  vy = V0(n) * ty

  epsabs = 0.0_dp
  epsrel = 1.0e-8_dp
  call dqxgs(azfunc,0.0_dp,1.0_dp,epsabs,epsrel,tz,&
       abserr,ier,limit,leniw,lenw,last,iwork,work)
  if (ier /= 0 ) print*,"dqxgs error in subroutine potin:",ier
  vz = V0(n) * tz

end subroutine accmid


function potfunc(t) result(res)
  real (kind=dp), intent(in) :: t
  real (kind=dp)             :: res
  real (kind=dp)             :: a,d,e
  
  d=1.0_dp -  ((1.0_dp-pintr(gn)*pintr(gn))*t*t)
  e=1.0_dp -  ((1.0_dp-qintr(gn)*qintr(gn))*t*t)

  ! Integral part of formula 12 of Cappellari 2002.  
  a = exp( -t*t / ( 2.0_dp*sigintr_km(gn)**2) *(gx*gx +(gy*gy)/d +(gz*gz)/e))
  res= a / sqrt( d * e)
end function potfunc

function axfunc(t) result(res)
  real (kind=dp), intent(in) :: t
  real (kind=dp)             :: res
  real (kind=dp)             :: a,d,e
  
  d=1.0_dp -  ((1.0_dp-pintr(gn)*pintr(gn))*t*t)
  e=1.0_dp -  ((1.0_dp-qintr(gn)*qintr(gn))*t*t)
  
  ! Integral part of formula 12 of Cappellari 2002.
  a = exp( -t*t / ( 2.0_dp*sigintr_km(gn)**2) *(gx*gx +(gy*gy)/d +(gz*gz)/e))
  ! dF/dx
  res = - gx / sigintr_km(gn)**2 * t * t *  a / sqrt( d * e)
end function axfunc 

function ayfunc(t) result(res)
  real (kind=dp), intent(in) :: t
  real (kind=dp)             :: res
  real (kind=dp)             :: a,d,e  

  d=1.0_dp -  ((1.0_dp-pintr(gn)*pintr(gn))*t*t)
  e=1.0_dp -  ((1.0_dp-qintr(gn)*qintr(gn))*t*t)

  ! Integral part of formula 12 of Cappellari 2002.
  a = exp( -t*t / ( 2.0_dp*sigintr_km(gn)**2) *(gx*gx +(gy*gy)/d +(gz*gz)/e))
  ! dF/dy
  res = - gy / sigintr_km(gn)**2 * t * t / d *  a / sqrt( d * e)
end function ayfunc

function azfunc(t) result(res)
  real (kind=dp), intent(in) :: t
  real (kind=dp)             :: res
  real (kind=dp)             :: a,d,e
  
  ! integral part of formula 12 of Cappellari 2002.
  d=1.0_dp -  ((1.0_dp-pintr(gn)*pintr(gn))*t*t)
  e=1.0_dp -  ((1.0_dp-qintr(gn)*qintr(gn))*t*t)
  
  ! Integral part of formula 12 of Cappellari 2002.
  a  = exp( -t*t / ( 2.0_dp*sigintr_km(gn)**2) *(gx*gx +(gy*gy)/d +(gz*gz)/e))
  ! dF/dz
  res = - gz / sigintr_km(gn)**2 * t * t / e *  a / sqrt( d * e)
end function azfunc

!++++++++++++++++++++++++++++++++++++
subroutine potin(n,x,y,z,res)
  integer(kind=i4b),intent(in)::n
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(out):: res
!------------------------------------
  real(kind=dp) :: x2,y2,z2
  real(kind=dp) :: O1,O2
  real(kind=dp) :: A11,A22,A33,A12,A23,A31

  x2=x*x
  y2=y*y
  z2=z*z

  A12 = - ( A1(n) - A2(n) )       / ( 1.0_dp       - pintr(n)**2)
  A23 = - ( A2(n) - A3(n) )       / ( pintr(n)**2  - qintr(n)**2)
  A31 = - ( A3(n) - A1(n) )       / ( qintr(n)**2  - 1.0_dp     )

  A11 = (1.0_dp/3.0_dp) * (2.0_dp / ( 1.0_dp     ) - A12 - A31  )
  A22 = (1.0_dp/3.0_dp) * (2.0_dp / ( pintr(n)**2) - A23 - A12  )
  A33 = (1.0_dp/3.0_dp) * (2.0_dp / ( qintr(n)**2) - A31 - A23  )

  O1 = -1.0_dp/(2.0_dp*sigintr_km(n)**2) * ( A1(n)*x2 + A2(n)*y2 + A3(n)*z2 )

  O2 =  1.0_dp/(8.0_dp*sigintr_km(n)**4) * &
       (        A11*x2*x2 +        A22*y2*y2 +        A33*z2*z2 + &
         2.0_dp*A12*x2*y2 + 2.0_dp*A23*y2*z2 + 2.0_dp*A31*z2*x2 )

  res  =  V0(n) / sqrt(1.0_dp-qintr(n)**2) * ( F(n) + O1 + O2)

end subroutine potin

!++++++++++++++++++++++++++++++++++++
subroutine accin(n,x,y,z,vx,vy,vz)
  integer(kind=i4b),intent(in)::n
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(out):: vx,vy,vz
!------------------------------------
  real(kind=dp) :: x2,y2,z2
  real(kind=dp) :: A11,A22,A33,A12,A23,A31

  x2=x*x
  y2=y*y
  z2=z*z

  A12 = - ( A1(n) - A2(n) )       / ( 1.0_dp       - pintr(n)**2)
  A23 = - ( A2(n) - A3(n) )       / ( pintr(n)**2  - qintr(n)**2)
  A31 = - ( A3(n) - A1(n) )       / ( qintr(n)**2  - 1.0_dp     )

  A11 = (1.0_dp/3.0_dp) * (2.0_dp / ( 1.0_dp     ) - A12 - A31  )
  A22 = (1.0_dp/3.0_dp) * (2.0_dp / ( pintr(n)**2) - A23 - A12  )
  A33 = (1.0_dp/3.0_dp) * (2.0_dp / ( qintr(n)**2) - A31 - A23  )

  vx = -V0(n) / sqrt(1.0-qintr(n)**2) * x/sigintr_km(n)**2 * ( A1(n) - &
       1.0_dp/(2.0_dp*sigintr_km(n)**2) * (A11*x2+A12*y2+A31*z2))

  vy = -V0(n) / sqrt(1.0-qintr(n)**2) * y/sigintr_km(n)**2 * ( A2(n) - &
       1.0_dp/(2.0_dp*sigintr_km(n)**2) * (A12*x2+A22*y2+A23*z2))

  vz = -V0(n) / sqrt(1.0-qintr(n)**2) * z/sigintr_km(n)**2 * ( A3(n) - &
       1.0_dp/(2.0_dp*sigintr_km(n)**2) * (A31*x2+A23*y2+A33*z2))

end subroutine accin

end module triaxpotent
