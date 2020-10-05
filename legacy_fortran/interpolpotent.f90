module interpolpot
  use numeric_kinds

  implicit none
  private

  ! calculate potential phi at (x,y,z) 
  public:: ip_potent
  
  ! calculate accel ax,ay,ax at (x,y,z)
  public:: ip_accel
  
  ! setup the constants for the potential
  public:: ip_setup

  ! setup the constants for the potential
  public:: ip_stop

  integer(kind=i4b),private,parameter :: nRadius=640,nTheta=64,nPhi=64
  !integer(kind=i4b),private,parameter :: nRadius=768,nTheta=64,nPhi=64
  !integer(kind=i4b),private,parameter :: nRadius=8,nTheta=4,nPhi=4
  real (kind=dp),private,allocatable,dimension(:,:,:,:) :: grid

  real (kind=dp),private::thetastep,phistep,rlogstep,rlogminp	,rmin2,rmax2
  
  ! Counters to check if the interpolation range was suffiecient   
  integer(kind=i4b),private :: count_Rmin,count_Rmax

  ! BUGS: uses sig_intr_km from triaxpotent module  
  !       Big speed penalty for converting cartesion coordinates to spherical.
  !           (Optimal fix: integrate in spherical coordinates)  

contains

subroutine ip_setup()
  use dmpotent, only : dm_setup
  use initial_parameters, only:conversion_factor
  integer (kind=i4b) :: error

  error=1
  
  allocate(grid(3,nphi,ntheta,nradius))
  
  call dm_setup()
  call ip_read(error)
  if (error == 0) call ip_testaccuracy(error)
  if (error /= 0) then 
     call ip_setup_grid()
     call ip_save()
     call ip_testaccuracy(error)
     if (error /=0) stop "failed to setup interpolation grid"
  endif
 
  print*,"  * Potential interpolation setup" 
end subroutine ip_setup

subroutine ip_stop()
  use dmpotent, only : dm_stop
  print*,"  * Number of times the points where inside the interpolation:",count_Rmin
  print*,"  * Number of times the points where outside the interpolation:",count_Rmax
  deallocate(grid)
  print*,"  * Potential interpolation stopped" 
  call dm_stop()
end subroutine ip_stop

subroutine ip_setup_grid()
  use initial_parameters, only:sigobs_km,rlogmax,conversion_factor,rlogmin
  use dmpotent, only : dm_setup, dm_accel
  use triaxpotent, only : sigintr_km
  integer (kind=i4b) :: i,j,k
  real (kind=dp) :: x,y,z,theta,phi,r,vx,vy,vz,rlogpolmax,t1,t2,t3, &
       maxdiff,t4,meandiff,tinyval
  real(kind=sp) :: time1, time2
  
 
  print*," Setting up interpolating grid"
  call cpu_time(time1)

  ! choose minimum from the sigobs (because sigobs<sigintr)
  rmin2 = (minval(sigobs_km)/10.0_dp)**2
  ! choose maximum from sigintr    (because sigintr>sigobs)
  rmax2 = (maxval(sigintr_km)*6.0_dp)**2

  ! try to ensure the inner radius is small enough for smallest orbits
  rmin2 = minval( (/ ((10**rlogmin)*0.01)**2 , rmin2 /))

  ! let the grid go out far enough to encompass all orbits in the library
   rmax2 = maxval( (/ ((10.0_dp**rlogmax)*1.05)**2 , rmax2*2.0_dp /))
  rlogminp = log10(sqrt(rmin2))
  rlogpolmax = log10(sqrt(rmax2))

  thetastep = pio2_d/(ntheta-1)
  phistep   = pio2_d/(nphi-1)
  rlogstep  =(rlogpolmax-rlogminp)/(nradius-1)

  print*,"inner and outer radius of interpolation"
  print*,sqrt(rmin2)/conversion_factor,sqrt(rmax2)/conversion_factor

  grid(:,:,:,:)=log(tiny(grid(1,1,1,1))) ! ~ log(tiny)
  tinyval=tiny(grid(1,1,1,1))

  do i=1,nradius 
     r     = 10 ** (rlogminp + (i-1) * rlogstep)
     do j=1,ntheta
        theta = (j-1.0_dp)*thetastep
        if (j == 1     ) theta=(j-0.5_dp)*thetastep
        if (j == ntheta) theta=(j-1.1_dp)*thetastep
        z = r * cos (theta)
        do k=1,nphi
           phi   = (k-1.0_dp)*phistep
           if (k == 1   ) phi=(k-0.5_dp)*phistep
           if (k == nphi) phi=(k-1.1_dp)*phistep
           x = r * sin (theta) * cos (phi)
           y = r * sin (theta) * sin (phi)
           call dm_accel(x,y,z,vx,vy,vz)
           if (-vx > tinyval*x) grid(1,k,j,i)=log(-vx/x)
           if (-vy > tinyval*y) grid(2,k,j,i)=log(-vy/y)
           if (-vz > tinyval*z) grid(3,k,j,i)=log(-vz/z)
           !print*,grid(1,k,j,i),log(-vx/x),log(tiny(grid(1,1,1,1))),vx,x
        end do
     end do
  end do

  !initialize Counters
  count_Rmin=0
  count_Rmax=0	

  call cpu_time(time2)
  print*,"   Time spent initializing :",time2-time1,"seconds"
end subroutine ip_setup_grid

subroutine ip_testaccuracy(error)
  use dmpotent, only : dm_accel
  integer(kind=i4b),intent(out)::error
  integer (kind=i4b) :: i,j,k
  real (kind=dp) :: x,y,z,theta,phi,r,vx,vy,vz,rlogpolmax,t1,t2,t3, &
       maxdiff,t4,meandiff,len,diff

  !
  ! test the accuracy of the interpolation
  !
  error=0
  print*,"  * Testing accuracy of interpolation" 
 
  rlogpolmax = log10(sqrt(rmax2))

  maxdiff=0.0_dp
  meandiff=0.0_dp
  t4=0.0_dp
  call random_seed()
  do i=1,20000 

     call random_number(r)
     call random_number(theta)
     call random_number(phi)

     r  =  10 ** (rlogminp + (rlogpolmax-rlogminp) * r)
     theta=theta*pio2_d
     phi  =phi  *pio2_d

     x = r * sin (theta) * cos (phi)
     y = r * sin (theta) * sin (phi)
     z = r * cos (theta)
     
     call dm_accel(x,y,z,vx,vy,vz)
     call ip_accel(x,y,z,t1,t2,t3)
     
     len=sqrt(vx**2+vy**2+vz**2)
     diff=sqrt(abs(vx-t1)**2+ abs(vy-t2)**2+ abs(vz-t3)**2)/len
     maxdiff = max( maxdiff , diff)
     meandiff = meandiff +  diff
     if (diff > 1.0e-3) then
       write(unit=*,fmt="(1es15.3,4f15.3,4e15.3)"),diff,theta*180/pi_d,phi*180/pi_d,log10(r),rlogminp,&
              (vx-t1)/t1,(vy-t2)/t2,(vz-t3)/t3,len
     endif
     end do

  close (unit=30)
  meandiff = meandiff/(20000.0)
  print*,"  * Max relative acceleration error" ,maxdiff
  print*,"  * Mean relative acceleration error",meandiff
  if ( maxdiff >= 1.0e-2) error=1 
  if ( meandiff >= 1.0e-4) error=1 
  if (error == 1) print*,"  * Required accuracy not reached"
end subroutine ip_testaccuracy


!+++++++++++++++++++++++++++++++++++
subroutine ip_accel(x,y,z,vx,vy,vz)
  use dmpotent, only : dm_accel
  real(kind=dp), intent(in) ::  x, y, z
  real(kind=dp), intent(out):: vx,vy,vz
  !------------------------------------
  real(kind=dp) :: r2,rlog,theta,phi,tf,pf,rf
  real(kind=dp),dimension(3) :: acc
  integer(kind=i4b) :: t,p,r
  real(kind=dp),parameter :: I=1.0_dp
  real(kind=dp),external ::  pow

  r2=x*x+y*y+z*z
  if (r2 > rmin2 .and. r2 < rmax2 ) then
     theta=atan2(sqrt(x*x+y*y),abs(z))
     phi=atan2(abs(y),abs(x))
     rlog=0.5_dp * log10(r2) ! = log10(sqrt(r2))

     t=floor(theta         / thetastep) + 1
     p=floor(phi           / phistep  ) + 1
     r=floor((rlog-rlogminp)/ rlogstep ) + 1
     tf=theta         /thetastep - floor(theta         / thetastep)
     pf=phi           /phistep   - floor(phi           / phistep  )
     rf=(rlog-rlogminp)/rlogstep  - floor((rlog-rlogminp)/ rlogstep )

     !if (p>1) then
     !   if (t>1) then ! p>1 and t>1
           acc(:)= abs((pf-I) * (tf-I) * (rf-I)) * grid(:,p  ,t  ,r  ) + &
                   abs((pf-I) * (tf-I) * (rf  )) * grid(:,p  ,t  ,r+1) + &
                   abs((pf-I) * (tf  ) * (rf  )) * grid(:,p  ,t+1,r+1) + &
                   abs((pf-I) * (tf  ) * (rf-I)) * grid(:,p  ,t+1,r  ) + &
                   abs((pf  ) * (tf-I) * (rf-I)) * grid(:,p+1,t  ,r  ) + &
                   abs((pf  ) * (tf-I) * (rf  )) * grid(:,p+1,t  ,r+1) + &
                   abs((pf  ) * (tf  ) * (rf  )) * grid(:,p+1,t+1,r+1) + &
                   abs((pf  ) * (tf  ) * (rf-I)) * grid(:,p+1,t+1,r  ) 
     !   else !p>1 and t=1
     !       acc(:)=abs((pf-I) * (rf  )) * grid(:,p  ,t+1,r+1) + &
     !              abs((pf-I) * (rf-I)) * grid(:,p  ,t+1,r  ) + &
     !              abs((pf  ) * (rf  )) * grid(:,p+1,t+1,r+1) + &
     !              abs((pf  ) * (rf-I)) * grid(:,p+1,t+1,r  ) 
     !   endif
     !else !p=1
     !   if (t>1) then ! p=1 and t>1
     !      acc(:)= abs((tf-I) * (rf-I)) * grid(:,p+1,t  ,r  ) + &
     !              abs((tf-I) * (rf  )) * grid(:,p+1,t  ,r+1) + &
     !              abs((tf  ) * (rf  )) * grid(:,p+1,t+1,r+1) + &
     !              abs((tf  ) * (rf-I)) * grid(:,p+1,t+1,r  ) 
     !   else     !p=1 and t=1
     !      acc(:)= abs((rf  )) * grid(:,p+1,t+1,r+1) + &
     !              abs((rf-I)) * grid(:,p+1,t+1,r  ) 
     !   endif
     !
     !endif
     vx = -x * exp(acc(1)) 
     vy = -y * exp(acc(2)) 
     vz = -z * exp(acc(3))
  else 
	 if (r2 < rmin2) count_Rmin =count_Rmin + 1_i4b
	 if (r2 > rmax2) count_Rmax =count_Rmax + 1_i4b  
	 call dm_accel(x,y,z,vx,vy,vz)
  endif
end subroutine ip_accel

subroutine ip_potent(x,y,z,pot)
use dmpotent, only : dm_potent
!use triaxpotent, only : tp_potent

  real(kind=dp), intent(in) ::  x, y, z
  real(kind=dp), intent(out):: pot
  call dm_potent(x,y,z,pot)
end subroutine ip_potent

subroutine ip_save()
  print*,"  *  Writing interpolation table to disk."
  open  (unit=35, file="interpolgrid", status="replace",action="write",form="unformatted")
  write (unit=35) 3_i4b,nphi,ntheta,nradius
  write (unit=35) 6_i4b,thetastep,phistep,rlogstep,rlogminp,rmin2,rmax2
  write (unit=35) grid(:,:,:,:)
  close (unit=35)
end subroutine ip_save

subroutine ip_read(error)
  integer(kind=i4b),intent(out)::error
  integer(kind=i4b) :: t1,t2,t3,t4
  print*,"  *  Reading interpolation table from disk."
  open  (unit=34, file="interpolgrid", status="old",position="rewind", & 
       action="read",form="unformatted",iostat=error)
  if (error == 0 ) &
       read (unit=34,iostat=error) t1,t2,t3,t4
  if (error /= 0 .or. t1/=3 .or. t2/=nphi .or. t3/=ntheta &
       .or. t4/=nradius) error = 1
  if (error == 0 ) &
       read (unit=34,iostat=error) t1,thetastep,phistep,rlogstep,rlogminp,rmin2,rmax2
  if (error /= 0 .or. t1/=6) error=1
  if (error == 0 ) &
       read (unit=34,iostat=error) grid(:,:,:,:)
  close (unit=34)
  if (error /= 0 ) &
       print*,"  *  Interpolation file could not be read"

end subroutine ip_read

end module interpolpot
