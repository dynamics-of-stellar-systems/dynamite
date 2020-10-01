 module binmass
  use numeric_kinds
  use initial_parameters
  implicit none
  private

  ! Main program for calculating the aperture mases
  public :: binmass_Main

  ! helper functions
  private :: integrand,integrate,loopoverbins,aperture_boxed_readfile
  private :: read_binningfile,write_apermass

  real (kind=dp),private, dimension(9) :: binmass_intergrand_global

  ! global for storing the calculated aperture masses before writing
  integer(kind=i4b),private :: global_n_apermass=0
!  real (kind=dp),private, dimension(20000) :: global_apermass=0.0_dp
  real (kind=dp),private, dimension(400_i4b**2*6) :: global_apermass=0.0_dp
  
contains

function integrand(x) result(res)
  real (kind=dp), intent(in) :: x
  real (kind=dp)             :: res

  real (kind=dp)             :: y0,y1,angle,qi,sigma,surf_km,isotwist
  real (kind=dp)             :: psfwidth,psfweight

  real (kind=dp)             :: f0,f1,k,c,p,psigma,sb,qb,surcor
  real (kind=dp)             :: alpha,dens

   y0       =binmass_intergrand_global(1)
   y1       =binmass_intergrand_global(2)
   angle    =binmass_intergrand_global(3)
   qi       =binmass_intergrand_global(4)
   sigma    =binmass_intergrand_global(5)
   surf_km  =binmass_intergrand_global(6)
   isotwist =binmass_intergrand_global(7)
   psfwidth =binmass_intergrand_global(8)
   psfweight=binmass_intergrand_global(9)

   ! compute new gaussian size after convolving with the psf: 
   sb = sqrt( sigma**2.0_dp+psfwidth**2.0_dp)
   qb = sqrt((( sigma*sigma*qi*qi + psfwidth*psfwidth))/(sigma*sigma + psfwidth*psfwidth))
   surcor= surf_km * qi/qb *(sigma**2) /(sb**2)

   ! Angle is the angle form the PA to the X-axis ccw
   ! To rotate the grid to the PA rotate of -angle ccw
   ! Isotwist is the rotation from the PA to the major axis of the gaussian ccw
   ! so this can just be added.
   alpha = -1.0_dp*angle + isotwist
   
   ! Glenn formula
   !k= sqrt( 1.0_dp+qb*qb+(1.0_dp-qb*qb)*cos(2.0_dp*alpha))
   !f0= 1.0_dp/(2.0_dp*k*qb*sb) * ( k*k*y0 - (1.0_dp-qb*qb)*sin(2.0_dp*alpha)*x)
   !f1= 1.0_dp/(2.0_dp*k*qb*sb) * ( k*k*y1 - (1.0_dp-qb*qb)*sin(2.0_dp*alpha)*x)
   !res = -psfweight * surcor * qb * sb * sqrt(pi_d) & 
   !      / k * (derf(f0)-derf(f1)) * exp ( - (x*x)/(k*k*sb*sb))

   ! Cappellari formula from ic1459 paper (2002) Appendix B3, formula (B6 and B7)
   p  = sqrt (1+qb*qb + (1.0_dp-qb*qb) * cos(2.0_dp*alpha))
   f0 = 1.0_dp/(2*p*qb*sb) *  ((1.0_dp-qb*qb)*x*sin(2.0_dp*alpha) - p*p * y0)
   f1 = 1.0_dp/(2*p*qb*sb) *  ((1.0_dp-qb*qb)*x*sin(2.0_dp*alpha) - p*p * y1)
   res  = psfweight * surcor * qb * sb * sqrt(pi_d)/p &
        * (derf(f0)-derf(f1)) * exp(-(x/(p*sb))**2.0_dp)
   
end function integrand

subroutine integrate(x0,x1,res)
  real   (kind=dp ),intent(in) ::  x0,x1
  real   (kind=dp ),intent(out):: res
!------------------------------------
  integer(kind=i4b),parameter  :: limit=10000,leniw=limit*3,lenw=limit*46
  integer(kind=i4b)            :: ier,last
  integer(kind=i4b),dimension(leniw) ::iwork
  real   (kind=dp )            :: epsabs,epsrel,abserr
  real   (kind=dp ),dimension(limit*46)  :: work

  epsabs = 0.0_dp
  epsrel = 1.0e-5_dp
  call dqxgs(integrand,x0,x1,epsabs,epsrel,res,&
       abserr,ier,limit,leniw,lenw,last,iwork,work)
  if (ier /= 0 .and. ier /=2 ) print*,"dqxgs error in subroutine integrate:",ier

end subroutine integrate

subroutine loopoverbins(mnx,mny,xsize,ysize,xbins,ybins,angle,psfwidth,psfweight)
  use initial_parameters , only : totalmass,ngauss_mge,qobs,psi_view
real (kind=dp),intent(in) :: mnx,mny,xsize,ysize,angle
integer(kind=i4b),intent(in) :: xbins,ybins
real (kind=dp),intent(in) :: psfwidth(:),psfweight(:)
real (kind=dp)            :: bx(xbins+1),by(ybins+1),grid(xbins,ybins)
real (kind=dp)            :: res
integer(kind=i4b)         :: i,j,k,l

! determine the boundaries in the x direction
do i=1,xbins+1
   bx(i)=mnx+ (i-1.0_dp)*xsize / xbins
end do

! determine the boundaries in the y direction
do i=1,ybins+1
   by(i)=mny+ (i-1.0_dp)*ysize / ybins
end do

do i=1,xbins
   print*,i,xbins
   do j=1,ybins

      grid(i,j)=0.0_dp

      do k=1,ngauss_mge
         do l=1,size(psfwidth)
            binmass_intergrand_global(1)=by(j)
            binmass_intergrand_global(2)=by(j+1)
            binmass_intergrand_global(3)=angle
            binmass_intergrand_global(4)=qobs(k)
            binmass_intergrand_global(5)=sigobs_km(k)
            binmass_intergrand_global(6)=surf_km(k)
            binmass_intergrand_global(7)=psi_obs(k)-psi_view ! = isotwist
            binmass_intergrand_global(8)=psfwidth(l)
            binmass_intergrand_global(9)=psfweight(l)

            call integrate(bx(i),bx(i+1),res)
            grid(i,j)=grid(i,j)+res

         end do
      end do
   end do
end do

print*," Percentage of the light in this aperture"
print*,sum(grid(:,:))/totalmass*100.0

! normalising the grid
grid(:,:)=grid(:,:)/totalmass

!Do binning and store in the global variable
call read_binningfile(reshape(grid,(/xbins*ybins/)))

end subroutine loopoverbins

subroutine binmass_main()
  use initial_parameters    , only : conversion_factor
real (kind=dp) :: mnx,mny,xsize,ysize,angle
integer(kind=i4b) :: xbins,ybins,i,j,naper,npsfgauss
real (kind=dp),allocatable,dimension(:) :: psfwidth,psfweight

print*,"  * How many aperture file to integrate?"
read*,naper
print*,naper
do i=1,naper
   print*,"  * aperture :",i
   call aperture_boxed_readfile(mnx,mny,xsize,ysize,xbins,ybins,angle)

   print*,"  * How many Gaussians in the psf for this aperture?"
   read*,npsfgauss
   print*,npsfgauss
   allocate (psfwidth(npsfgauss),psfweight(npsfgauss))
   do j=1,npsfgauss
      read*,psfweight(j),psfwidth(j)
   end do 
   
   psfwidth(:)=psfwidth(:)* conversion_factor

   print*,psfwidth

   call loopoverbins(mnx,mny,xsize,ysize,xbins,ybins,angle,psfwidth,psfweight)
   deallocate(psfwidth,psfweight)
end do

call write_apermass()

end subroutine binmass_main

subroutine aperture_boxed_readfile(mnx,mny,xsize,ysize,xbins,ybins,angle)
  use initial_parameters , only : conversion_factor
  character(len = 80) :: file,string
  integer (kind=i4b) :: handle
  !----------------------------------------------------------------------
  real (kind=dp),intent(out) :: mnx,mny,xsize,ysize,angle
  integer(kind=i4b),intent(out) :: xbins,ybins
  print*, "  * Reading boxed aperture file."
  print*,"  * What's the filename of the aperture file ? :"
  read*, file

  open (unit=handle,file=file,action="read",status="old"&
       &,position="rewind")
  print*,"  * Checking type."
  read(unit=handle,fmt=*) string
  
  select case (string)
  case ("#counterrotation_polygon_aperturefile_version_1")
     stop "Aperture type not supported"
  case ("#counter_rotation_boxed_aperturefile_version_2")
     ! empty
  case default
     stop " Unkown aperture type"
  end select

  print*,"  *  Reading box info"
  print*,"  *  Order: begin(x,y)"
  read(unit=handle,fmt=*) mnx,mny
  print*,"      size(x,y) "
  read(unit=handle,fmt=*) xsize,ysize
  print*,"      rotation"
  read(unit=handle,fmt=*) angle
  angle=angle*(pi_d/180.0_dp)
  print*,"      bin(x,y)"
  read(unit=handle,fmt=*) xbins,ybins

  ! convert arcsec into km
  mnx=mnx*conversion_factor
  mny=mny*conversion_factor
  xsize=xsize*conversion_factor
  ysize=ysize*conversion_factor

  print*,"   Total bins " , xsize*ysize
  print*,"   begin      " , mnx,mny
  print*,"   size       " , xsize,ysize
  print*,"   rotation   " , angle/(pi_d/180.0_dp)
  print*,"   binx       " , xbins,ybins
  print*," "
  print*,"  * Finished reading aperture"
  close (unit=handle)

end subroutine aperture_boxed_readfile

subroutine read_binningfile(grid)
real    (kind=dp) ,dimension(:),intent(in) :: grid
integer (kind=i4b) :: i,bins,minst,maxst
character (len=256) :: string
integer (kind=i4b),dimension(size(grid)) :: bin_order

! start this array at 0. This is used for unused pixels
real    (kind=dp) ,dimension(0:size(grid)) :: binned

print*,"  * Aperture: ",i
print*,"  * Give the filename of the binning file."
read*,string
print*,"  * Opening: ",string
open (unit=30,file=string,action="read",status="old"&
     &,position="rewind")
read (unit=30,fmt=*) string
if (string /= "#Counterrotaton_binning_version_1") &
     stop " Wrong version of file"
read (unit=30,fmt=*) bins
print*,"  * pixels in this aperture:",bins

if (bins /= size(grid) )&
     stop " bin and aperture file do not agree on the amount of pixels"
read (unit=30,fmt=*) bin_order(:)
close (unit=30)

print*,"  * Bins in this aperture:",maxval(bin_order)

if (maxval(bin_order) > size(grid)) &
     stop " Cannot store the highest bins in the grid"

binned=0.0_dp

do i=1, size(grid)
   binned(bin_order(i))=binned(bin_order(i))+grid(i)
end do

minst=global_n_apermass+1
maxst=global_n_apermass+maxval(bin_order)

if (maxst > size(global_apermass)) &
     stop "  TOo few elements in static array global_apermass"

global_apermass(minst:maxst) = global_apermass(minst:maxst) &
                             + binned(1:maxval(bin_order))

global_n_apermass=maxst

print*,"  * Done binning"
end subroutine read_binningfile

subroutine write_apermass()
integer(kind=i4b) :: i
character (len=256) :: string
print*,"  * Give the file name of the aperture mass output file"
read*,string
print*,string

print*,"  * Writing file"
open (unit=28,file=string,action="write", status="new")
write (unit=28, fmt="(i5)") global_n_apermass
do i=1,global_n_apermass
   write (unit=28, fmt="(i5,es20.12)") i, global_apermass(i)
end do
close (unit=28)
end subroutine write_apermass

end module binmass

program triaxmassbin
  use binmass
  use initial_parameters
  use triaxpotent
  use numeric_kinds
  implicit none
 
  call iniparam()
  call tp_setup()
  
  call binmass_main()


end program Triaxmassbin
