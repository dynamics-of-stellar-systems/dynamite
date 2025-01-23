
! 04 okt 2005 RvdB
! - fixed bugs to allow for hermax =2
! - BUG add_regularization: Added second if to avoid evaluation
!   of (orbint(4,reg(1)) if it is undefined.
! 19 okt
! - Do not evenize orbitlibrary 2 per default
! 03 feb 2005
! - split donnls routine in two. Now there is version with swapfile and without
! 19 dec 2005 regularisation of box orbits is flipped in I3
!
! $Id: triaxnnls.f90,v 1.4 2011/10/25 08:48:45 bosch Exp $

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The modified or added lines are commented by (! BT)
! In this file, reading orbit libraries and making observables are modified.
! If  (Omega = 0) it will be the same as before, it will read box orbit and tube orbit, in addition
! flip the sign of velocity of tube orbits to have a counter-rotating orbit library as the third one.
! if (Omega != 0) flipping tube orbits will be canceled and just read retrograde orbits (in orblibbox).
! So for (Omega = 0) we have three orbit libraries and for (Omega != 0) we have only two-orbit libraries.

! by Behzad Tahmasebzadeh   February 2022
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! 2025/01/19 Thomas I. Maindl
! can read the orbit library and losvd histograms from either one file (likely
! (orblib.dat) or separate files (likely orblib_qgrid.dat,
! orblib_losvd_hist.dat).

module NNLSfit

use numeric_kinds
implicit none
private

character(len=256), private :: orbfile1_1, orbfile1_2, orbfile2_1, orbfile2_2, outroot

integer (kind=i4b),private :: hermax
integer (kind=i4b),private :: nconstr,massconstr,orbitsinfit,convecl
real    (kind=dp ),private :: velunit,rmsresid,vm,sg,intrmasserror,reguldelta

real(kind=dp),dimension(:),allocatable,private ::con,econ,enermass, lcut
real(kind=dp),dimension(:),allocatable,private :: apermass,orbweight
real(kind=dp),dimension(:,:),allocatable,private :: velmom,dvelmom,orbmat
real(kind=dp),dimension(:,:,:),allocatable,private :: intrinsicmass

integer(kind=i4b),dimension(:,:),allocatable,private :: orbint

integer(kind=i4b),dimension(:),allocatable,private :: orbinmat,orbtype

! boolean for the use of regularisation
logical , private :: use_reg

public :: readorblibs,readobservations,donnls,allpred,makekinem
public :: stopnnlsfit,add_regularization,FIND_AVERAGE_ORBTYPE
public :: ADD_REG_TO_ORBMAT,FIND_ORBINT,donnls_nosave
public :: donnls_galahad,reg_forcetoabel,donnls_abel

contains

  subroutine find_average_orbtype(orb,avg)
    integer (kind=i4b),dimension(:),intent(in) :: orb
    integer (kind=i4b),intent(out) :: avg
    integer (kind=i4b) :: i,o1,o3,o4,n

    o1=0
    o3=0
    o4=0

    n=size(orb)
    do i=1,n
       if (orb(i) == 1) o1=o1+1
       if (orb(i) == 3) o3=o3+1
    enddo

    avg=4
    if (o1 > o4 ) avg=1
    if (o3 > o1 .and. o3 > o4 ) avg=3

 end subroutine find_average_orbtype


  subroutine readorblibs()
    use initial_parameters , only : conversion_factor, Omega
    integer (kind=i4b) :: i,j,k,l,m,iv,l1,l2,l3

    integer (kind=i4b) :: t1,t2,t3,t4,norbtotal,orboffset
    integer (kind=i4b) ::smom1,slr1,sth1,sph1
    integer (kind=i4b) ::smom2,slr2,sth2,sph2
    integer (kind=i4b),dimension(2) :: norb,ndith,nvhist,nconstrl
    integer (kind=i4b),dimension(4) ::  orbnumbers
    integer(kind=i4b) :: ivmin,ivmax
    real(kind=dp),dimension(0:hermax) :: hh
    real(kind=dp),dimension(2) :: dvhist
    real(kind=dp),dimension(:),allocatable :: quad_lr &
                              ,quad_lth,quad_lph,veltmp
    integer(kind=i4b),dimension(:),allocatable :: orbtypes
    real(kind=dp),dimension(:,:,:,:),allocatable :: quad_light
    real(kind=dp),dimension(:,:),allocatable :: velhist
    real(kind=dp) :: tmp, veldiff_l, norm_veltmp, dvm, dsg, vm_ol,sg_ol, gam_ol, h12_ol

    integer (kind=i4b) :: naperture_cut
    real(kind=dp),dimension(:), allocatable ::  aperture_cut

    integer(kind=i4b) :: qgrid_unit1, qgrid_unit2, hist_unit1, hist_unit2, read_from

    print*, "  * What is the velocity scaling factor?"
    read*,velunit
    print*,velunit

    print*, "  * Reading orbit libraries"

        ! the orbit library can be in one or two files
        print *, "Give first orbitlibrary name, file 1:"
        read (unit=*, fmt="(a256)") orbfile1_1
        orbfile1_1 = trim(orbfile1_1)
        print *, orbfile1_1
        print *, "Give first orbitlibrary name, file 2:"
        read (unit=*, fmt="(a256)") orbfile1_2
        orbfile1_2 = trim(orbfile1_2)
        print *, orbfile1_2
        qgrid_unit1 = 27
        open (unit=qgrid_unit1, file=orbfile1_1, action="read", &
            status="old", form="unformatted", position="rewind")
        if (orbfile1_1 == orbfile1_2) then
            hist_unit1 = qgrid_unit1
        else
            hist_unit1 = 37
            open (unit=hist_unit1, file=orbfile1_2, action="read", &
                status="old", form="unformatted", position="rewind")
        end if
        read (unit=qgrid_unit1) norb(1), t1, t2, t3, ndith(1)
        print *, norb(1), t1, t2, t3, ndith(1)
        read (unit=qgrid_unit1) smom1, sph1, sth1, slr1
        print *, smom1, sph1, sth1, slr1
        smom1 = 16
        if (slr1*sth1*sph1 /= massconstr) &
            stop " The number of mass constraints is not the same"
        ! remember that N bins have N+1 boundaries
        allocate (quad_lr(slr1 + 1), quad_lth(sth1 + 1), quad_lph(sph1 + 1), &
                  quad_light(smom1, sph1, sth1, slr1))
        read (unit=qgrid_unit1) quad_lr(:)
        print *, "  * Intrinisic grid radii:"
        print *, quad_lr/conversion_factor
        read (unit=qgrid_unit1) quad_lth(:)
        read (unit=qgrid_unit1) quad_lph(:)
        read (unit=hist_unit1) nconstrl(1), nvhist(1), dvhist(1)

        print *, "Give second orbitlibrary name, file 1:"
        print *, "Remember that this library is not (velocity) mirrored"
        read (unit=*, fmt="(a256)") orbfile2_1
        orbfile2_1 = trim(orbfile2_1)
        print *, orbfile2_1
        print *, "Give second orbitlibrary name, file 2:"
        print *, "Remember that this library is not (velocity) mirrored"
        read (unit=*, fmt="(a256)") orbfile2_2
        orbfile2_2 = trim(orbfile2_2)
        print *, orbfile2_2
        qgrid_unit2 = 28
        open (unit=qgrid_unit2, file=orbfile2_1, action="read", &
            status="old", form="unformatted", position="rewind")
        if (orbfile2_1 == orbfile2_2) then
            hist_unit2 = qgrid_unit2
        else
            hist_unit2 = 38
            open (unit=hist_unit2, file=orbfile2_2, action="read", &
                status="old", form="unformatted", position="rewind")
        end if
        read (unit=qgrid_unit2) norb(2), t1, t2, t3, ndith(2)
        read (unit=qgrid_unit2) smom2, sph2, sth2, slr2
        !if ( smom1 /= smom2 .or. slr1 /= slr2 .or. sth1 /= sth2 .or. sph1 /= sph2 ) &
        !     stop "  Intrinsic grid is not the same size as the other library"

        read (unit=qgrid_unit2) quad_lr(:)
        read (unit=qgrid_unit2) quad_lth(:)
        read (unit=qgrid_unit2) quad_lph(:)
        read (unit=hist_unit2) nconstrl(2), nvhist(2), dvhist(2)

        print *, nconstrl(2), nvhist(2), dvhist(2)
        deallocate (quad_lr, quad_lth, quad_lph)

    if (nconstr /= nconstrl(1) .or. nconstr /= nconstrl(2)) &
         stop " number of apertures is not consistent"

    norbtotal=2*norb(1)+norb(2)
    ! in case of rotating frame box orbit library is our retrogtete orbits (we do not flip the tubr orbits)
    if (Omega /= 0.0_dp ) then
      norbtotal=norb(1)+norb(2)
    endif

    print*," Allocating Matrix :",size(con),norbtotal
    allocate(orbmat(size(con),norbtotal),orbinmat(norbtotal))
    allocate(orbint(4,norbtotal),orbtype(norbtotal))
    allocate(lcut(norbtotal))
    allocate(aperture_cut(nconstr))

    orbmat(:,:) = 0.0_dp

    orboffset=0

    do i=1,2 ! loop over orbit libraries
       print*," Reading orbitlibrary :",i
       allocate(velhist(-nvhist(i):nvhist(i),nconstr),&
                veltmp (-nvhist(i):nvhist(i)))
       allocate(orbtypes(ndith(i)**3))

            do j = 1, norb(i) ! loop over orbits
                print *, " Orbit :", j, orboffset
                if (i == 1) then
                    read_from = qgrid_unit1
                else
                    read_from = qgrid_unit2
                end if
                read (read_from) t1, orbnumbers(1:4)
                if (t1 /= j) stop " orbit number does not match"
                read (unit=read_from) orbtypes(:)
                read (unit=read_from) quad_light(:, :, :, :)

                if (i == 1) then
                    read_from = hist_unit1
                else
                    read_from = hist_unit2
                end if
                velhist(:, :) = 0.0_dp
                do k = 1, nconstr ! loop over apertures
                    read (unit=read_from) ivmin, ivmax
                    if (ivmin <= ivmax) &
                        read (unit=read_from) velhist(ivmin:ivmax, k)
                end do
                ! Loop over orbit with reverse angular momentum.
                ! Do not loop when i=2 (=second orbitlibray)

         if (Omega /= 0.0_dp ) then     ! (BT) in case of figure rotation we only have two orbit libraries)
            l1=1
            l2=1
            l3=1
         else
            l1=1
            l2=-1+(i-1)*2
            l3=-2
         end if

         do k=l1,l2,l3           ! (BT) do not flip tubes in case of figure rotation

                orboffset=orboffset+1
	         	lcut(orboffset) = 0.0
                ! Store which orbit is where:
                orbinmat(orboffset)=(j + (i-1) * norb(1)) * k
                ! store orbit integral identifier (E, I2, I3)
                orbint(:,orboffset)=orbnumbers(:) * (/ 1_i4b+(i-1_i4b)*(-2_i4b),k,1_i4b,1_i4b /)
                ! find average orbtype
                call find_average_orbtype(orbtypes,orbtype(orboffset))

                ! Store total mass equals 1 constraint
                orbmat(1,orboffset)=1.0_dp

                ! Store intrinsic masses in the matrix
                orbmat(2:massconstr+1,orboffset)=&
                     reshape(quad_light(1,:,:,:),(/massconstr/))

		! at the start of loop over apertures, set to cut of this orbits to all apertures to be 0
		naperture_cut = 0
		do l =1, nconstr
		  aperture_cut(l) = 0
		end do

                do l=1,nconstr ! loop over aperture
                   !orbmat(massconstr+l,orboffset)=sum(velhist(:,l))
                   if (sum(velhist(:,l)) < 0.0 ) stop "sum is less then zero"
                   if (sum(velhist(:,l)) > 0.0_dp ) then
                   ! Store aperture mass
                    orbmat(massconstr+1+l,orboffset)=sum(velhist(:,l))
                   ! Flip sign of velocity field of nessecary:
                   if (k==1) then

                      !if ( minval(orbtypes) > 3) then
                      !   ! Evenize orbit if all orbit elements are boxes.
                      !   !(The second orbit library should only contain boxes)
                      !!!!! Disabled. should mirror the intrinsic mass weights
                      !!!!! Too
                      !   do m=-nvhist(i),nvhist(i)
                      !      veltmp(m) = 0.5_dp*(velhist(m,l)+velhist(-m,l))
                      !   end do
                      !else
                         veltmp(:)=velhist(:,l)
                      !end if
                   else
                      do m=-nvhist(i),nvhist(i)
                         veltmp(m)=velhist(-m,l)
                      end do
                   end if

                   ! Set the first and last point in the histogram to zero before
                   !continuing. These cannot be used, because they correspond to all
                   ! velocities to infinity. We assume that their contribution the
                   ! GH moments are zero. This will in fact automatically be the case
                   ! if sigma_obs << velocity histogram size.
                   veltmp(-nvhist) = 0.0_dp
                   veltmp( nvhist) = 0.0_dp

                   ! Calculate the mean velocity and dispersion for this observation in
                   ! the program units
                   vm=velmom(l,1)/velunit
                   sg=velmom(l,2)/velunit

                   dvm=dvelmom(l,1)/velunit
                   dsg=dvelmom(l,2)/velunit
                   ! Calculate the Gauss-Hermite moments for the observed mean velocity
                   ! and velocity dispersion. in this calculation the value of gamma is
                   ! kept equal to unity. The routine (automatically) returns zero for
                   ! all gauss-Hermite moments if there is no weight in the histogram
                   ! at all.
                   call getgauher (vm,sg,veltmp,&
                        nvhist(i),nvhist(i),dvhist(i),hh,hermax,1.0_dp,1)
                   ! Loop over the Gauss-Hermite moments, and store in the matrix

                   ! Some attention is required here. To normalize the histogram w.r.t.
                   ! integration over velocity one must divide by (totweight*dvhist).
                   ! The matrix entry should be the Gauss-Hermite moment, multiplied by
                   ! the fraction of the time the orbit is seen at the current aperture.
                   ! This is totweight. Hence, to get the correct matrix entry, we
                   ! should still divide by dvhist.


		   !!!!!!!!!!!!!!!!!!!!!!!!!!
		   ! Lingcut
		   ! label the orbits being cut as lcut = 1
		   ! Fit gauss to the orbit histgram veltem
  		   !call gaussfit(veltmp,nvhist(i),nvhist(i),dvhist(i), gam_ol,vm_ol,sg_ol,h12_ol,5)
		   !vm_ol = meanvel(veltmp, nvhist(i), dvhist(i))
		   vm_ol = 0
		   norm_veltmp = 0

                   do m=-nvhist(i),nvhist(i)
                      vm_ol = vm_ol + veltmp(m) * m * dvhist(i)
		      norm_veltmp = norm_veltmp + veltmp(m)
                   end do
		   vm_ol = vm_ol / norm_veltmp

		   !print*, "vm, mean vm_ol", vm, vm_ol
		   !call gaussfit(veltmp,nvhist(i),nvhist(i),dvhist(i), gam_ol,vm_ol,sg_ol,h12_ol,5)
		   !print*, "gauss vm_ol", vm_ol
		   !veldiff_l = minval([abs(vm_ol - vm -dvm), abs(vm_ol - vm + dvm)])

		   veldiff_l = abs(vm_ol - vm)

		   ! cut 1
		   if ( abs(vm) / sg > 1.5 .and. veldiff_l / (sg) > 3.0 .and. vm*vm_ol < 0.0 )  then
		   !if (veldiff_l / (sg) > 4.0 .and. vm*vm_ol < 0.0 )  then
                      aperture_cut(l) = 1
                      naperture_cut = naperture_cut + 1
		   end if

		   ! cut 2  m7vc

		   !if (abs(vm)/sg > 2.5 .and. vm_ol*vm_ol > 0.7 * 0.7 * (sg*sg + vm *vm)  .and. vm*vm_ol < 0.0 )  then
		   !	hh(1) = 3.0
		   !	lcut(orboffset) = 1.0
		   !end if

		   ! set unreasonable value to hh(1) for the reversed orbits
		   !!!!!!!!!!!!!!!!!!!!!!!!!!

                   orbmat(massconstr+nconstr+l+1:massconstr+l+hermax*nconstr+1 &
                        :nconstr,orboffset) = hh(1:hermax)/dvhist(i)
                   end if
                end do

		!naperture_cut  ! recored the number of aperture that want the orbit to be cut
		!aperture_cut   ! size of aperture, is aperture_cut =1, then set hh(1) = 3.0

                if (naperture_cut > 1) then
                   lcut(orboffset) = 2.0
                   do l=1,nconstr ! loop over aperture
                      if (aperture_cut(l) > 0) then
                         orbmat(massconstr+nconstr+l+1,orboffset) = 3.0/dvhist(i)
				! reset hh1 = 3.0 for this apertures with orbit orboffset
                      endif
                   end do
                end if


           end do

       end do
       deallocate(velhist,veltmp,orbtypes)
    end do

        close (unit=qgrid_unit1)
        if (qgrid_unit1 /= hist_unit1) close (unit=hist_unit1)
        close (unit=qgrid_unit2)
        if (qgrid_unit2 /= hist_unit2) close (unit=hist_unit2)

    orbitsinfit=orboffset
    print*,"  * Orbits in fit :",orbitsinfit

    deallocate(quad_light)
    print*, " Done"
  end subroutine readorblibs

subroutine readobservations()
  use initial_parameters , only : iniparam_bar,nEner,nI2,nI3,orbit_dithering
  character (len=256) :: infile
  integer (kind=i4b) :: ios,Nmass1,nmass2,nmass3,i,j,k,tmp,ntemp,nherm
  integer (kind=i4b) :: cmin,cmax,iherm

  real(kind=dp) :: ftmp,minmass,masserror
  real(kind=dp),allocatable,dimension(:) :: kinerrscale
  real(kind=dp),dimension(:,:,:),allocatable :: tmass

  print*,"  * Read parameters.in"
  call iniparam_bar()

  print*," Give the regularization strenght"
  read (unit=*,fmt=*) reguldelta
  print*,reguldelta
  use_reg=.false.
  if (reguldelta > 0.0_dp) then
     print*,"  * Regularization Enabled"
     use_reg=.true.
     reguldelta=sqrt(reguldelta/(2.0_dp*nEner*nI2*nI3/(orbit_dithering**3)))
     print*,"  * Regularization Enabled. reguldelta = ",reguldelta
  endif
  if (.not. use_reg ) print*,"  * Regularisation Disabled"

  print*, "  * Give root file name for results"
  read (unit=*, fmt="(a256)") outroot
  print*, trim(outroot)
  print*, " "

  print*, "Give input file with intrinsic masses"
  print*, "This is mass mass_qgrid.dat"
  read (unit=*, fmt="(a256)") infile
  print*,trim(infile)

  open (unit=41,file=infile,status="old",action="read",position="rewind")

  read (unit=41,fmt=*) Nmass1,nmass2,nmass3  ! I3,I2,E
  print*," Number of intrinsic mass constraints:"
  print*,Nmass1,nmass2,nmass3,"=",nmass1*nmass2*nmass3
  allocate(intrinsicmass(Nmass1,nmass2,nmass3))
  allocate(tmass(Nmass1,nmass2,nmass3))
  read (unit=41,fmt=*) intrinsicmass(:,:,:)
  close(unit=41)

  massconstr=nmass1*nmass2*nmass3
  minmass=1.0_dp
  ! Check sanity of the masses
  do i=1,nmass1
     do j=1,nmass2
        do k=1,nmass3
           if (intrinsicmass(i,j,k) > 0.0_dp .and. intrinsicmass(i,j,k) < minmass ) &
                minmass = intrinsicmass(i,j,k)
           if (intrinsicmass(i,j,k) < 0.0_dp ) stop " grid mass is less then zero!"
        end do
     end do
  end do

! -----------------------------------------------
! Read mass_radmass.dat for regularization

  allocate(enermass(nEner/orbit_dithering))

  print*, "Give input file with intrinsic masses"
  print*, "This is mass mass_qgrid.dat"
  !read (unit=*, fmt="(a256)") infile
  infile='datfil/mass_radmass.dat'
  print*,infile
  open (unit=41,file=infile,status="old",action="read",position="rewind")

  read (unit=41,fmt=*) i
  print*," Number of energies:"
  print*,i,"=",nEner/orbit_dithering
  if (i/=nEner/orbit_dithering) stop 'number of energies dont match'
  read (unit=41,fmt=*) enermass(:)
  close(unit=41)

  print*,enermass(:)

!-----------------------------------------------------

  print*,"Give input file with projected mass"
  print*, "at the kinematical constraint points"
  read (unit=*, fmt="(a256)") infile
  print*,infile
  ! this is trueapweights.dat
  open (unit=43,file=infile,status="old",action="read",position="rewind",iostat=ios)

  ! read the file

  read (unit=43, fmt=*,iostat=ios) nconstr
  print*,"number of aperture constraints:",nconstr
  ! read the input data
  allocate(apermass(nconstr))
  do i=1,nconstr
     read (unit=43,fmt=*,iostat=ios) tmp, apermass(i)
     if (ios /= 0) stop  "error reading aperture masses"
  end do
  close (unit=43)


!------------------------------------------------------

  print*, "  Give the number of hermite moments to be fitted"
  read (unit=*, fmt=*) hermax
  allocate (kinerrscale(hermax))

  print*, "Give input file with observed kinematical data"
  print*, "at the constraint points"

  read (unit=*, fmt="(a256)") infile
  print*,infile
  open (unit=44,file=infile,status="old",action="read",position="rewind")
  read (unit=44, fmt=*) ntemp
  print*, "ntemp,nconstr",ntemp,nconstr,hermax
  if (ntemp /= nconstr) stop "input error in GET_KinDATA"

  allocate(velmom(nconstr,hermax),dvelmom(nconstr,hermax))
  ! read the input data

  velmom(:,:)=0.0_dp
  dvelmom(:,:)=1.0e31_dp

  print*,"  * reading kinematics",hermax
  do i=1,nconstr
     read (unit=44,fmt=*) tmp,&
          velmom(i,1),dvelmom(i,1), velmom(i,2),dvelmom(i,2),nherm ,&
          (velmom(i,j),dvelmom(i,j),j=3,minval( (/hermax, nherm/) ))
  end do


  print*,"Give the kinematic systematic error for each moment."
  read (unit=*,fmt=*) kinerrscale(:)
  print*,kinerrscale(:)
  print*,'stope1',minval(kinerrscale)
  if (minval(kinerrscale) < 0.0_dp) stop " scale is smaller then zero"

  do i=1,hermax
     dvelmom(:,i)=sqrt(dvelmom(:,i)**2 + kinerrscale(i)**2)
  end do

  ! this maximum error exist because this is what it sets zero error's to.
  print*,'stope2', maxval(dvelmom)
  if (maxval(dvelmom) > 1.0e32_dp) stop " maximum error on kinematics exceeded"
  print*,'stope3', minval(dvelmom)
  if (minval(dvelmom) .le. 0.0_dp) stop " min error on kin is less than zero"

  close (unit=44)

  deallocate(kinerrscale)
!---------------------------------------------------------

 convecl=nmass1*nmass2*nmass3 + 1 + nconstr*(hermax+1)

if (use_reg) convecl=convecl + &
     3*3*nEner*ni2*nI3/(orbit_dithering**3)

print*,"Total constraints allocated:",convecl

allocate(con(convecl),&
     econ(convecl))

con(:)=0.0_dp
econ(:)=0.0_dp



print*," Give the relative error on the intrinsic mass"
read (unit=*,fmt=*) intrmasserror
print*,intrmasserror

print*,"  * Adding total mass equals ",sum(intrinsicmass(:,:,:))," constraint"

cmin=1
cmax=cmin

econ(cmin) = minval((/ intrmasserror/10.0_dp, abs(1.0_dp- sum(intrinsicmass(:,:,:))) /))
 con(cmin) = sum(intrinsicmass(:,:,:))

print*,"  * Constraint on total mass is within: ",econ(cmin)

! Intrinsic mass bins

cmin=cmin+1
cmax=cmin+massconstr-1

con(cmin:cmax)=reshape(intrinsicmass(:,:,:),(/massconstr/))
econ(cmin:cmax)=abs(reshape(intrinsicmass(:,:,:),(/massconstr/))*intrmasserror)

! Check sanity of mass error.
do i=cmin,cmax
   !if ( econ(i) .le. 0.0_dp ) stop " intrinsic mass equal to or less then zero"
   if ( econ(i) .le. 0.0_dp ) econ(i)=1.0e-16
end do

print*," * Give the relative error on the projected mass in the apertures"
read (unit=*,fmt=*) masserror
print*,masserror

cmin=cmax+1
cmax=cmin+nconstr-1
 con(cmin:cmax)=    apermass(:)
econ(cmin:cmax)=abs(apermass(:)*masserror)


print*,"  * Storing Velocity and Sigma"
! V and sigma
do i=1,minval((/2_i4b,hermax/))
cmin=cmax+1
cmax=cmin+nconstr-1
 con(cmin:cmax)=velmom(:,i)*0.0_dp
econ(cmin:cmax)=sqrt(0.5_dp)*dvelmom(:,i)/velmom(:,2)*apermass(:)
end do

print*,"  * Storing Hermite moments"
! other Hermite moments
do i=3,hermax
cmin=cmax+1
cmax=cmin+nconstr-1
 con(cmin:cmax)=     velmom(:,i)*apermass(:)
econ(cmin:cmax)=abs(dvelmom(:,i)*apermass(:))
end do



end subroutine readobservations

subroutine add_regularization()

  real (kind=dp) :: rnorm
  integer (kind=i4b) :: i,mi3,iE,iI2,iI3
  integer (kind=i4b),dimension(3) :: reg,ints
  integer (kind=i4b),dimension(4,orbitsinfit) :: orbinttmp

  if (use_reg) then

  print*,"Adding regularization"

  if (reguldelta <= 0.0_dp) stop " zero regularization not implemented"

   ! add error to the whole range to avoid devide by zero.
  econ(nconstr*(hermax+1)+massconstr+1:convecl) = 1.0_dp

  !  store the orbint so we can change it temporarily
  orbinttmp=orbint

  ! Change the orbit integrals identifier of the box orbits so that they
  ! match up with the tube orbits.
  mI3= maxval ( orbint(3,:))

  ! flip box orbits so that there orbit numbers are aligned with the tube orbits
  where ( orbint(1,:) < 0 )
     orbint(3,:) = 1 + mI3*2 - orbint(3,:)
     orbint(1,:) = orbint(1,:) * (-1)
  end where

  ! move the flipped tube orbits (i2<0) so that there is no gap.
  ! This makes the it so that there are orbits with I2=0
  where ( orbint(2,:) < 0 )
     orbint(2,:) = orbint(2,:) + 1
  end where


  do i=1,orbitsinfit
     !print*,"orbit ",i,orbint(1:3,i)
     ! Energy reg
     reg(1) = find_orbint(orbint(1:3,i) + (/ -1_i4b, 0_i4b, 0_i4b /))
     reg(2) = i
     reg(3) = find_orbint(orbint(1:3,i) + (/  1_i4b, 0_i4b, 0_i4b /))
     ! regularize if 3 orbits are found
     if (minval(reg) > 0 ) call add_reg_to_orbmat(reg)

     ! I2 reg
     reg(1) = find_orbint(orbint(1:3,i) + (/  0_i4b,-1_i4b, 0_i4b /))
     reg(2) = i
     reg(3) = find_orbint(orbint(1:3,i) + (/  0_i4b, 1_i4b, 0_i4b /))
     ! regularize if 3 orbits are found
     if (minval(reg) > 0 )  call add_reg_to_orbmat(reg)

     ! I3 reg
     reg(1) = find_orbint(orbint(1:3,i) + (/  0_i4b, 0_i4b,-1_i4b /))
     reg(2) = i
     reg(3) = find_orbint(orbint(1:3,i) + (/  0_i4b, 0_i4b, 1_i4b /))
     ! do not regularize if the orbit is border long axis tube orbit that
     ! should be correlated with a box orbit
     if (minval(reg) > 0 .and. orbint(4,i) == 0  ) then
          if (orbint(4,reg(1)) ==0)  &
          call add_reg_to_orbmat(reg)
       endif
  end do

  !reverse the orbint change
  where ( orbint(2,:) <= 0 )
     orbint(2,:) = orbint(2,:) - 1
  end where

  ! change structure so that the minus-rotation orbits lines up with the boxes
  where ( orbint(3,:) < mi3-1 .or. orbint(3,:) > mI3+2 .or. &
          ( orbint(2,:) > 0 .and. orbint(3,:) <= mI3 ))
     orbint(1,:)=0
     orbint(2,:)=0
     orbint(3,:)=0
  end where

  orbint(2,:)=abs(orbint(2,:))

  do i=1,orbitsinfit
     ! I3 reg along the -Lz and box border
     reg(1) = find_orbint(orbint(1:3,i) + (/  0_i4b, 0_i4b,-1_i4b /))
     reg(2) = i
     reg(3) = find_orbint(orbint(1:3,i) + (/  0_i4b, 0_i4b, 1_i4b /))
     ! do not regularize if the orbit is border long axis tube orbit that
     ! should be correlated with a box orbit
     if (minval(reg) > 0 .and. orbint(4,i) == 0 ) then
        if (orbint(4,reg(1)) == 0) &
          call add_reg_to_orbmat(reg)
     endif
  end do

  ! restore orbint
  orbint=orbinttmp

  endif ! add regularisation
end subroutine add_regularization

subroutine add_reg_to_orbmat(orbs)
        ! Add orbs regularization constraint to orbmat.

        ! Regularisation is perturbed to keep the matrix from becoming singular
        ! This helps the minimizer converge, but does not change the solution

        integer (kind=i4b), intent(in),dimension(3) :: orbs
        integer (kind=i4b) :: i,row
        integer (kind=i4b),save :: nreg=0
        real    (kind=dp)  :: r

        call random_number(r)

        nreg=nreg+1
        row = nconstr*(hermax+1)+massconstr + 1 + nreg

        if (row > convecl) stop " array to small for regularization"

        con(row)  = 0.0_dp  + (r-0.5_dp)*reguldelta*1.0e-5
        econ(row) = reguldelta

        call random_number(r)
		r= 1.0_dp + (r-0.5_dp) * 1.0e-2
        orbmat(row,orbs(1)) = -1.0_dp / enermass(orbint(1,orbs(1))) * r
        call random_number(r)
		r= 1.0_dp + (r-0.5_dp) * 1.0e-2
        orbmat(row,orbs(2)) =  2.0_dp / enermass(orbint(1,orbs(2))) * r
        call random_number(r)
		r= 1.0_dp + (r-0.5_dp) * 1.0e-2
        orbmat(row,orbs(3)) = -1.0_dp / enermass(orbint(1,orbs(3))) * r

end subroutine add_reg_to_orbmat


function find_orbint(int) result (res)
  integer (kind=i4b) :: i
  integer (kind=i4b),dimension(3),intent(in) :: int
  integer (kind=i4b) :: res

  res = 0
  do i=1,orbitsinfit
     if (orbint(1,i) == int(1) .and. orbint(2,i) == int(2) .and. &
          orbint(3,i) == int(3) ) then
        res = i
        !   break
     endif
  enddo

end function find_orbint

subroutine donnls()
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Solve the NNLS problem for the orbital weights

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real (kind=dp), dimension(convecl) :: b, zz
  real (kind=dp), dimension(orbitsinfit) :: w
  !real (kind=dp), dimension(size(orbmat,1),size(orbmat,2)) :: orbtmp
  ! Working space and other arrays
  integer (kind=i4b) :: len, i, j, k
  real (kind=dp) :: rnorm
  integer (kind=i4b), dimension(orbitsinfit) :: index
  integer (kind=i4b) :: it, ior, mode
	character (len=260) :: scratchfile

  print*," Starting NNLS in function"

  ! first element is the total mass = 1 constraint.
  ! in qp this has zero error, but nnls cant cope with that.
  if (econ(1) <= 0.0_dp ) econ(1) = con(1)*1.0e-2_dp

  allocate(orbweight(orbitsinfit))

  !open (unit=36,file=trim(outroot)//"orbmat.tmp",status="replace",action="write",form="unformatted")
  open (unit=36,status="SCRATCH",action="readwrite",form="unformatted")
	!INQUIRE(unit=36,name=scratchfile)
	!print*,'scratchfile is ',scratchfile
  print*, "Storing nonzero ORBMAT elements on file..."
  len = 0
  do k=1,orbitsinfit
     do j=1,size(con)
        if (orbmat(j,k) /= 0.0_dp) then
           len = len + 1
           write (unit=36) j, k, orbmat(j,k)
        endif
     end do
  end do
  print*, "Nonzero elements (%):", len*100.0_sp/(orbitsinfit*size(con))
  !close (unit=36,STATUS='KEEP')

  b(:)=con(:)/econ(:)

   do it=1,size(con)
      if (econ(it) <= 0.0_dp) stop "Zero error"
           orbmat(it,:) = orbmat(it,:) / econ(it)
   end do


   print*, "doing NNLS fit ..."
   print*, " "

   ! Get the NNLS fit

   call nnls(orbmat(:,1:orbitsinfit),size(con),size(con),orbitsinfit,&
        b,orbweight,rmsresid,w,zz,index,mode)

   ! write some info to screen

   print*, "NNLS fit performed"
   write (unit=*, fmt="(es20.12,i6)") rmsresid, mode
   print*, " "

   print*, "Restoring nonzero ORBMAT from file to memory..."
   !open (unit=36,file=scratchfile,status="old",action="read",&
   !     form="unformatted",position="rewind")
	 rewind(unit=36)
   orbmat = 0.0_dp
   do i=1,len
      read (unit=36) j, k, orbmat(j,k)
   end do
   close (unit=36)

  !print*, "Shrinking orbmat.tmp"
  !open (unit=36,file=trim(outroot)//"orbmat.tmp",status="replace",action="write",form="unformatted")
  !write (unit=36) 0.0_dp
  !close (unit=36)


   print*,"sum orbweight :",sum(orbweight)
 end subroutine donnls

subroutine donnls_store_matrix_and_read_solution()
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Solve the NNLS problem for the orbital weights

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real (kind=dp), dimension(convecl) :: b, zz
  real (kind=dp), dimension(orbitsinfit) :: w
  !real (kind=dp), dimension(size(orbmat,1),size(orbmat,2)) :: orbtmp
  ! Working space and other arrays
  integer (kind=i4b) :: len, i, j, k
  real (kind=dp) :: rnorm
  integer (kind=i4b), dimension(orbitsinfit) :: index
  integer (kind=i4b) :: it, ior, mode

  print*," Starting NNLS in function"

  allocate(orbweight(orbitsinfit))

  open (unit=36,file="orbmat.gal",status="replace",action="write",form="unformatted")
  print*, "Storing Matrix in file..."
  write (unit=36) convecl,orbitsinfit
  write (unit=36) nconstr,massconstr,hermax
  write (unit=36) massconstr+nconstr,orbitsinfit
  write (unit=36) con(1:massconstr+nconstr),econ(1:massconstr+nconstr)
  write (unit=36) orbmat(1:massconstr+nconstr,:)
  write (unit=36) size(orbmat,1)-massconstr-nconstr,orbitsinfit
  write (unit=36) con(massconstr+nconstr+1:size(orbmat,1)),&
       econ(massconstr+nconstr+1:size(orbmat,1))
  write (unit=36) orbmat(massconstr+nconstr+1:size(orbmat,1),:)

  close (unit=36)

  open (unit=27,file='sol',action="read", &
         status="old",form="formatted",position="rewind")

  read (unit=27,fmt=*),orbweight(:)
  close (unit=27)

end subroutine donnls_store_matrix_and_read_solution

subroutine donnls_galahad(save,FitZeroMoment)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Solve the NNLS problem for the orbital weights

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  integer (kind=i4b),intent(in) :: save,FitZeroMoment
  !real (kind=dp), dimension(size(orbmat,1),size(orbmat,2)) :: orbtmp
	real(kind=dp),dimension(:,:),allocatable :: orbtmp
  ! Working space and other arrays
  integer (kind=i4b) :: len, i, j, k
  real (kind=dp),dimension(size(econ)) :: tmpcon,tmpecon
  real (kind=dp) :: rescale,ran,nc
  integer (kind=i4b), dimension(orbitsinfit) :: index
  integer (kind=i4b) :: it, ior, mode,a,b,c,d

   print*,"  * Solver function memory allocated."

  allocate(orbweight(orbitsinfit))

  if (save .eq. 1) then
    !open (unit=36,file=trim(outroot)//"_orbmat.tmp",status="replace",action="write",form="unformatted")
    open (unit=36,status="SCRATCH",action="readwrite",form="unformatted")
		print*, "Storing nonzero ORBMAT elements on file..."
    len = 0
    do k=1,orbitsinfit
       do j=1,size(con)
          if (orbmat(j,k) /= 0.0_dp) then
             len = len + 1
             write (unit=36) j, k, orbmat(j,k)
          endif
       end do
    end do
    print*, "Nonzero elements (%):", len*100.0_sp/(orbitsinfit*size(con))
    !close (unit=36)
	else
    allocate(orbtmp(size(orbmat,1),size(orbmat,2)))
    orbtmp=orbmat
  endif

  ! temporary store orbmat
  !orbtmp=orbmat
  tmpcon=con
  tmpecon=econ

  ! a-b is the observable part of the matrix
  a=1+1+massconstr+nconstr
  b=size(orbmat,1)
  ! c-d is the constraint part of the matrix
  c=1
  d=1+massconstr+nconstr

  ! if fitZeroMoment = 1 then we fit the zeroth velocity moment.
	! Otherwise it is a constraint.
  if (fitZeroMoment .eq. 1) then
		print*,' Using both the 3D and zeroth moments as constraints'
		a= a - nconstr
		d= d - nconstr
	endif

  ! FIXME: Scale such that the lowest mass in mass_qgrid*rescale gt 1e-5
  rescale=1.0e5

  orbmat(:,:)=orbmat(:,:)/rescale

  ! adding error to matrix
  con(a:b)=con(a:b)/econ(a:b)
   do it=a,b
		  if (econ(it) <= 0.0_dp)  print*,"Zero error in ",it,econ(it)
      if (econ(it) <= 0.0_dp)  stop "Zero error"
           orbmat(it,:) = orbmat(it,:) / econ(it)
   end do

  ! normalize the constraint part of the matrix
  do it=c,d
     if (abs(con(it)) .gt. 1.0d-40) then
		  !print*,it,con(it),econ(it)
		  !  call random_number(ran)
		  nc=1.0_dp !+ 0.02*ran -0.01
          orbmat(it,:) = orbmat(it,:) / con(it) * nc
          econ(it) = abs((abs(con(it))+econ(it))/abs(con(it)) * nc - nc)
          con(it)=nc
          !print*,it,con(it),econ(it)
     endif
  end do

!  Alternative normalisation the constraint part of the matrix
!  This normalization does not improve the solution convergence.
!do it=c,d
!   if (abs(con(it)) .gt. 0.0_dp) then
!		  !print*,it,con(it),econ(it)
!		  !  call random_number(ran)
!		  nc= 1.0_dp / sqrt (sum (orbmat(it,:)**2))
!        orbmat(it,:) = orbmat(it,:) * nc
!        econ(it) = econ(it) * nc
!        con(it)= con(it) * nc
!        !print*,it,con(it),econ(it),nc
!   endif
!end do


   call nnlsgal(orbmat(a:b,:), con(a:b),  &
                orbmat(c:d,:), con(c:d),econ(c:d),orbweight)
   orbweight=orbweight/rescale

   allocate(orbmat(size(con),size(orbtype)))

	 if (save .eq. 1) then
     print*, "Restoring nonzero ORBMAT from file to memory..."
     !open (unit=36,file=trim(outroot)//"_orbmat.tmp",status="old",action="read",&
     !   form="unformatted",position="rewind")
		 rewind(unit=36)
     orbmat = 0.0_dp
     do i=1,len
      read (unit=36) j, k, orbmat(j,k)
     end do
     close (unit=36)

     !print*, "Shrinking orbmat.tmp"
     !open (unit=36,file=trim(outroot)//"_orbmat.tmp",status="replace",action="write",form="unformatted")
     !write (unit=36) 0.0_dp
     !close (unit=36)
	 else
	    orbmat=orbtmp
      deallocate(orbtmp)
		endif
   con=tmpcon
   econ=tmpecon

   print*,sum(orbweight)
   print*,'  * Solving done'
end subroutine donnls_galahad


   subroutine nnlsgal(obsmat,cons,fixmat,conf,econf,orbweight)
   USE GALAHAD_QPB_DOUBLE
   USE GALAHAD_QPT_DOUBLE                           ! Double precision
   USE GALAHAD_PRESOLVE_DOUBLE                      ! Double precision
   USE GALAHAD_SYMBOLS
   USE GALAHAD_SMT_DOUBLE
   real(kind=dp),dimension(:),intent(in) ::conf,econf,cons
   real(kind=dp),dimension(:,:),intent(in) :: fixmat,obsmat
   real(kind=dp),dimension(:),intent(out) ::orbweight
   real(kind=dp),dimension(:,:),allocatable :: tmpmat
   REAL ( KIND = dp ), PARAMETER :: infinity = 10.0_dp ** 20
   REAL ( KIND = dp ), PARAMETER :: rescale=1.0_dp!e4
   TYPE ( QPT_problem_type ) :: p
   TYPE ( QPB_data_type ) :: data
   TYPE ( QPB_control_type ) :: control
   TYPE ( QPB_inform_type ) :: info
   TYPE ( PRESOLVE_control_type ) :: pcontrol
   TYPE ( PRESOLVE_inform_type )  :: inform
   TYPE ( PRESOLVE_data_type )    :: pdata

   INTEGER (kind=i4b) :: n, m, h_ne, a_ne ,i,j,t ,s
   INTEGER (kind=i4b),ALLOCATABLE, DIMENSION( : ) :: C_stat, B_stat
   integer :: stat
   real(kind=dp) :: mxv

   n = size(obsmat,2)
   if (size(obsmat,2) /= size(fixmat,2)) stop ' fixmat wrong size'
   if (size(cons) /= size(obsmat,1)) stop ' cons wrong size'
   m = size(fixmat,1)
   if (m /= size( conf)) stop '  conf wrong size'
   if (m /= size(econf)) stop ' econf wrong size'

   print*,"n,m,obs:",n,m,size(cons)
! start problem data
   ALLOCATE( p%G( n ), p%X_l( n ), p%X_u( n ) )
   ALLOCATE( p%C( m ), p%C_l( m ), p%C_u( m ) )
   ALLOCATE( p%X( n ), p%Y( m ), p%Z( n ) )
   ALLOCATE( B_stat( n ), C_stat( m ) )
   ! P%X0 needs to be allocated because QPB forgets to do it. (BUG in qpb)
   ALLOCATE( p%X0(n) )
   p%new_problem_structure = .TRUE.           ! new structure
   p%n = n ; p%m = m ; p%f = 5e5          ! dimensions & objective constant
	 print*,'call matmul'
   p%G = matmul(-transpose(obsmat),cons(:))
   print*,minval(p%G),maxval(p%G)
   p%C_l = conf(:)-econf(:)              ! constraint lower bound
   p%C_u = conf(:)+econf(:)              ! constraint upper bound
   p%X = 1.0_dp/real(n) ; p%Y = 0.0_dp ; p%Z = 0.0_dp ! start from zero
   p%X_l =  0.0_dp              ! variable lower bound
   p%X_u =  infinity ! 1.0_dp *rescale           ! variable upper bound
!   p%rho_g = 1.0_dp ; p%rho_b = 1.0_dp        ! initial penalty parameters


!  sparse co-ordinate storage format
   CALL smt_put( p%H%type, 'COORDINATE',s) ! Specify co-ordinate
   CALL smt_put( p%A%type, 'COORDINATE',s )
!  Make A
   a_ne=0
   mxv=maxval(abs(fixmat(:,:)))*0.0_dp
    do i=1,m
      do j=1,n
         if (abs(fixmat(i,j)) > mxv) a_ne=a_ne+1
      end do
   end do
   p%A%ne = a_ne
   print*,"sparsity A",a_ne*1.0/(n*1.0*m)*100.0
   ALLOCATE( p%A%val( a_ne ), p%A%row( a_ne ), p%A%col( a_ne ) )
   t=0
   do i=1,m
      do j=1,n
         if (abs(fixmat(i,j)) > mxv) then
            t=t+1
            p%A%val(t)=fixmat(i,j)
            p%A%row(t)=i
            p%A%col(t)=j
         endif
      end do
   end do

   !print*,'Converting A to row sparse'
   call  QPT_A_from_C_to_S(p,stat)
   IF ( stat /= 0 ) STOP ' Failed to Convert'


   print*,'Computing H'
   allocate(tmpmat(n,n))
   !FIXME: Replace this loop with the appropriate LAPACK/BLAS call.
   do i=1,n
      do j=1,i
	         tmpmat(i,j)=sum(obsmat(:,i)*obsmat(:,j))
      enddo
   enddo
   print*,'Made matrix'
   deallocate(orbmat)
   h_ne=0
   mxv=maxval(abs(tmpmat))*0.0_dp
   t=0
   do i=1,n
      do j=1,i
         if (abs(tmpmat(i,j)) > mxv .or. i==j ) H_ne=H_ne+1
      end do
   end do
   p%H%ne = h_ne
   print*,"sparsity H:",h_ne/((n*(n+1.0))/2.0)*100.0
   ALLOCATE( p%H%val( h_ne ), p%H%row( H_ne ), p%H%col( H_ne ) )
   t=0
   do i=1,n
      do j=1,i
         if (abs(tmpmat(i,j)) > mxv .or. i==j ) then
            t=t+1
            p%H%val(t)=tmpmat(i,j)
            p%H%row(t)=i
            p%H%col(t)=j
         endif
      end do
   end do
   deallocate(tmpmat)
   print*,'Converting H to row sparse'
   call  QPT_H_from_C_to_S(p,stat)
   IF ( stat /= 0 ) STOP ' Failed to Convert'

   print*,'Calling Solver'
   !  presolve
   CALL PRESOLVE_initialize( pcontrol, inform, pdata )
   IF ( inform%status /= 0 ) STOP
   pcontrol%print_level =  GALAHAD_TRACE              ! Ask for some output
   pcontrol%get_x=.true.  ! only need this one.
   pcontrol%c_accuracy = 1d-12 !1d-6    good
   !pcontrol%z_accuracy = 1d-12 !1d-6
   pcontrol%pivot_tol = 1d-20  !1d-10   good
   pcontrol%min_rel_improve = 1d-20 !1d-10
   pcontrol%max_growth_factor= 1d4 ! 1d8
   pcontrol%infinity = infinity
   !pcontrol%termination=GALAHAD_FULL_PRESOLVE
   !! apply presolving to reduce the problem
   !CALL PRESOLVE_apply( p, pcontrol, inform, pdata )
   IF ( inform%status /= 0 ) STOP

   CALL QPB_initialize( data, control)!,info )         ! Initialize control parameters
!   CALL QPB_initialize( data, control,info )         ! Initialize control parameters
   !   control%rho_g = 1.0_dp ; control%rho_b = 1.0_dp        ! initial penalty parameters
   control%infinity = infinity                  ! Set infinity
   control%print_level = 1

   !for GALAHAD 2.3.0000
   !control%restore_problem=0
   ! control%generate_sif_file=.true.
   ! control%sif_file_name=''
   ! control%print_level = 4
   control%maxit=n+m
   control%cg_maxit=(n+m+1)*4.0   ! *3.0 /2.0

   !galahad 2.4 for QP
   !control%scale=7
   !!control%presolve=.TRUE.
   !control%QPB_control%print_level=1
   !control%QPA_control%print_level=1
   !control%QPC_control%print_level=1
   !control%CQP_control%print_level=1
   !control%SCALE_control%print_level=1
   !control%PRESOLVE_control%print_level=1
   !control%QPB_control%maxit=n+m
   !control%QPB_control%cg_maxit=(n+m+1)*3.0 /2.0 !(N+M+1)*4
   !!control%CQP_control%restore_problem=0
   !!control%CQP_control%space_critical=.TRUE.
   !control%quadratic_programming_solver='qpb'

   CALL QPB_solve( p,data, control, info,C_stat,B_stat )
   !CALL QPC_solve( p,C_stat, B_stat,data, control, info )
   IF ( info%status == 0 ) THEN                 !  Successful return
	print*,sum(p%X)
   ELSE                                         !  Error returns
     WRITE( 6, "( ' QPB_solve exit status = ', I6 ) " ) info%status
   END IF

   ! restore the solved reduced problem to the original formulation
   !pcontrol%print_level =  1
   !pcontrol%get_q=.true.
   !pcontrol%get_z=.true.
   !pcontrol%get_c=.true.
   !pcontrol%get_y=.true.
   !pcontrol%get_x=.true.
   !pcontrol%get_c_bounds=.true.
   !CALL PRESOLVE_restore( p, pcontrol, inform, pdata )
   print*,'resolving done'
   IF ( inform%status /= 0 ) STOP

   control%print_level =0
   pcontrol%print_level =0
   CALL QPB_terminate( data, control, info )    !  delete internal workspace
   CALL PRESOLVE_terminate( pcontrol, inform, pdata )
   orbweight(:)=p%X(:)

   print*,'finished finding orbweight'
   ! Deallocate arrays
   DEALLOCATE( p%G, p%X_l, p%X_u )
   DEALLOCATE( p%C, p%C_l, p%C_u )
   DEALLOCATE( p%X, p%Y, p%Z)
   DEALLOCATE( p%H%ptr, p%A%ptr)
   DEALLOCATE( p%H%val, p%H%col)
   DEALLOCATE( p%A%val, p%A%col)
 END subroutine nnlsgal


subroutine donnls_nosave()
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Solve the NNLS problem for the orbital weights

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real (kind=dp), dimension(convecl) :: b, zz
  real (kind=dp), dimension(orbitsinfit) :: w
  ! real (kind=dp), dimension(size(orbmat,1),size(orbmat,2)) :: orbtmp
	real(kind=dp),dimension(:,:),allocatable :: orbtmp
  real (kind=dp),dimension(size(econ)) :: tmpcon,tmpecon
  ! Working space and other arrays
  integer (kind=i4b) :: len, i, j, k,c,d
  real (kind=dp) :: rnorm
  integer (kind=i4b), dimension(orbitsinfit) :: index
  integer (kind=i4b) :: it, ior, mode

  print*," Starting NNLS in function"

  allocate(orbweight(orbitsinfit))
	allocate(orbtmp(size(orbmat,1),size(orbmat,2)))
  orbtmp(:,:)=orbmat(:,:)
  tmpcon=con
  tmpecon=econ

  ! first element is the total mass = 1 constraint.
  ! in qp this has zero error, but nnls cant cope with that.
  if (econ(1) <= 0.0_dp ) econ(1) = con(1)*1.0e-2_dp

  b(:)=con(:)/econ(:)

   do it=1,size(con)
	 	  if (econ(it) <= 0.0_dp)  print*,"Zero error in ",it,econ(it)
      if (econ(it) <= 0.0_dp .and. it > 1) stop "Zero error"
      orbmat(it,:) = orbmat(it,:) / econ(it)
   end do

   ! --- save NNLS matrix to *_orbmat.out ---
   open (unit=36, file=trim(outroot)//"_orbmat.out", status="replace", action="write", form="formatted")
   print *, "Storing nonzero ORBMAT elements on file..."
   write (unit=36, fmt=*) orbitsinfit, size(con)
   len = 0
   do k = 1, orbitsinfit
      do j = 1, size(con)
         if (orbmat(j, k) /= 0.0_dp) then
            len = len + 1
         end if
         write (unit=36, fmt=*) orbmat(j, k)
      end do
   end do

   ! --- save NNLS rhs vector to *_orbmat.out file ---
   do j = 1, size(con)
      write (unit=36, fmt=*) b(j)
   end do

   print *, "Nonzero elements (%):", len*100.0_sp/(orbitsinfit*size(con))


   print*, "doing NNLS fit ..."
   print*, " "

   call nnls(orbmat(:,1:orbitsinfit),size(con),size(con),orbitsinfit,&
        b,orbweight,rmsresid,w,zz,index,mode)

   ! --- save NNLS solution to *_orbmat.out file ---
   do k = 1, orbitsinfit
      write (unit=36, fmt=*) orbweight(k)
   end do
   close (unit=36, STATUS='KEEP')


   ! write some info to screen
   print*, "NNLS fit performed"
   write (unit=*, fmt="(es20.12,i6)") rmsresid, mode
   print*, " "

   orbmat(:,:)=orbtmp(:,:)

   con=tmpcon
   econ=tmpecon
	 deallocate(orbtmp)

   print*,"sum orbweight :",sum(orbweight)
 end subroutine donnls_nosave


subroutine allpred()
  real(kind=dp),dimension(size(con)) :: conpr
  real(kind=dp) :: chimass,parchi,chisq,chitotmass,chireg
  integer (kind=i4b) :: b,e,s,i,j
  real (kind=dp) ,dimension(nconstr,hermax) :: modvel
  print*, " Making predictions"

 ! Calculate the predictions for each constraint by straight
 ! matrix*vector multiplication

 conpr(:) = matmul(orbmat(:,1:orbitsinfit),orbweight(:))

 open (unit=30,file=trim(outroot)//"_nnls.out", status="replace", &
      action="write")

 chimass = sum(((conpr(2:massconstr+1)-con(2:massconstr+1))/econ(2:massconstr+1))**2)
 write (unit=*,fmt="(A6,i6,2es20.6,f20.6)") "mass",&
      massconstr,chimass,sqrt(chimass),chimass/massconstr
 write (unit=30,fmt="(A6,i6,2es20.6,f20.6)") "mass",&
      massconstr,chimass,sqrt(chimass),chimass/massconstr

 do i=0,hermax
    b=massconstr+1+i*nconstr+1
    e=massconstr+1+(i+1)*nconstr
    parchi = sum(((conpr(b:e)-con(b:e))/econ(b:e))**2)
    write (unit=*,fmt="(2i6,2es20.6,f20.6)")  i,&
         nconstr,parchi,sqrt(parchi),parchi/nconstr
    write (unit=30,fmt="(2i6,2es20.6,f20.6)")  i,&
         nconstr,parchi,sqrt(parchi),parchi/nconstr

 end do

 ! Regularisation
 b=massconstr+(hermax+1)*nconstr+1+1
 e=size(con)
 chireg=0.0_dp
 if (use_reg) chireg = sum(((conpr(b:e)-con(b:e))/econ(b:e))**2)

 write (unit=*,fmt="(A6,i6,2es20.6,f20.6)") "reg",&
      e-b+1,chireg,sqrt(chireg),chireg/(e-b+1)
 write (unit=30,fmt="(A6,i6,2es20.6,f20.6)") "reg",&
      e-b+1,chireg,sqrt(chireg),chireg/(e-b+1)



 !! Total mass eq 1.0 constraint
 !chitotmass = (conpr(convecl)-con(convecl) / econ(convecl))**2
 !write (unit=*,fmt="(A7,i6,2es15.6,f15.6)") "totmass",&
 !     1,chitotmass,sqrt(chitotmass),chitotmass/nconstr
 !write (unit=30,fmt="(A7,i6,2es15.6,f15.6)") "totmass",&
 !     1,chitotmass,sqrt(chitotmass),chitotmass/nconstr

 ! Total chi2
 chisq = sum(((conpr(2:convecl)-con(2:convecl))/econ(2:convecl))**2)
 write (unit=*,fmt="(A8,i6,2es20.6,f20.6)")  "total",&
      size(con),chisq,sqrt(chisq),chisq/(size(con))
 write (unit=30,fmt="(A8,i6,2es20.6,f20.6)")  "total",&
      size(con),chisq,sqrt(chisq),chisq/(size(con))

! kinematical chi2
 b=massconstr+1+1*nconstr+1
 e=massconstr+1+(hermax+1)*nconstr
 chisq = sum(((conpr(b:e)-con(b:e))/econ(b:e))**2)
 write (unit=*,fmt="(A6,i6,2es20.6,f20.6)")  "kin",&
      e-b+1,chisq,sqrt(chisq),chisq/(e-b+1)
 write (unit=30,fmt="(A6,i6,2es20.6,f20.6)")  "kin",&
      e-b+1,chisq,sqrt(chisq),chisq/(e-b+1)



 close (unit=30)

!print*,"Total mass constraint:"
!print*,conpr(convecl),con(convecl) ,econ(convecl)

 open (unit=30,file=trim(outroot)//"_con.out", status="replace", &
      action="write")
 write (unit=30, fmt=*) size(con)
 do i=1,size(con)
    write (unit=30,fmt="(i6,3es15.6)") i,con(i),conpr(i),econ(i)
 end do
 close (unit=30)

! velmom(nconstr,hermax) modvel

!do i=1,minval((/2_i4b,hermax/))
!   modvel(:,i) = velmom(:,i) + velmom(:,2)/(apermass(:)*sqrt(0.5_dp))/ &
!        conpr(massconstr+nconstr*i+1:massconstr+nconstr*(i+1))
!end do
!
!do i=3,hermax
!
!   modvel(:,i) =  conpr(massconstr+nconstr*i+1:massconstr+nconstr*(i+1))/apermass(:)
!end do
!
!open (unit=30, file=trim(outroot)//"_kinem-hh-bugged.out", status="replace",&
!     action="write")
!write (unit=30, fmt=*) nconstr,hermax
!do i=1,nconstr
!   write (unit=30, fmt="(i8,30es13.5)") i,0.0_dp,0.0_dp,&
!        (modvel(i,j),velmom(i,j),dvelmom(i,j),j=1,hermax)
!end do
!close (unit=30)

open (unit=30,file=trim(outroot)//"_orb.out", status="replace", &
      action="write")
 write (unit=30, fmt=*) orbitsinfit
do i=1,orbitsinfit
   write (unit=30, fmt="(6i8,es13.5,es13.5,es13.5,es13.5,es13.5,es13.5,6i8)") i, orbint(1:4,i),orbtype(i),orbweight(i), lcut(i)
end do

close (unit=30)

deallocate (orbmat)
end subroutine allpred


subroutine makekinem(select)
  use initial_parameters , only : conversion_factor, Omega
  ! select defines the type of output
  integer (kind=i4b),intent(in) :: select
  integer (kind=i4b) :: i,j,k,l,m,orb_qualify,l1,l2,l3    ! behzad

  integer (kind=i4b) :: t1,t2,t3,t4,norbtotal,orboffset
  integer (kind=i4b) :: smom1,slr1,sth1,sph1
  integer (kind=i4b) :: smom2,slr2,sth2,sph2
  integer (kind=i4b),dimension(2) :: norb,ndith,nvhist,nconstrl
  integer(kind=i4b) :: ivmin,ivmax
  real(kind=dp),dimension(0:hermax) :: hh
  real(kind=dp),dimension(2) :: dvhist
  real(kind=dp),dimension(:),allocatable :: quad_lr,quad_lth,quad_lph,veltmp,veltmp1
  integer(kind=i4b),dimension(:),allocatable :: orbtypes
  real(kind=dp),dimension(:,:,:,:),allocatable :: quad_light,quad_store
  real(kind=dp),dimension(:,:  ),allocatable :: velhist,velstore
  real(kind=dp) :: vm,sg,h12,gam
  real (kind=dp) ,dimension(nconstr,hermax) :: modvel
  real (kind=dp) ,dimension(nconstr)        :: fitapermass
  character (len=260) :: outroot_kinem

        integer(kind=i4b) :: qgrid_unit1, qgrid_unit2, hist_unit1, hist_unit2, read_from

  print*, "  * Reading orbit libraries"

        qgrid_unit1 = 27
        open (unit=qgrid_unit1, file=orbfile1_1, action="read", &
            status="old", form="unformatted", position="rewind")
        if (orbfile1_1 == orbfile1_2) then
            hist_unit1 = qgrid_unit1
        else
            hist_unit1 = 37
            open (unit=hist_unit1, file=orbfile1_2, action="read", &
                status="old", form="unformatted", position="rewind")
        end if
        read (unit=qgrid_unit1) norb(1), t1, t2, t3, ndith(1)
        read (unit=qgrid_unit1) smom1, sph1, sth1, slr1
        smom1 = 16
        allocate (quad_lr(slr1 + 1), quad_lth(sth1 + 1), quad_lph(sph1 + 1), &
                  quad_light(smom1, sph1, sth1, slr1), quad_store(smom1, sph1, sth1, slr1))
        read (unit=qgrid_unit1) quad_lr(:)
        read (unit=qgrid_unit1) quad_lth(:)
        read (unit=qgrid_unit1) quad_lph(:)
        read (unit=hist_unit1) nconstrl(1), nvhist(1), dvhist(1)

        qgrid_unit2 = 28
        open (unit=qgrid_unit2, file=orbfile2_1, action="read", &
            status="old", form="unformatted", position="rewind")
        if (orbfile2_1 == orbfile2_2) then
            hist_unit2 = qgrid_unit2
        else
            hist_unit2 = 38
            open (unit=hist_unit2, file=orbfile2_2, action="read", &
                status="old", form="unformatted", position="rewind")
        end if
        read (unit=qgrid_unit2) norb(2), t1, t2, t3, ndith(2)
        read (unit=qgrid_unit2) smom2, sph2, sth2, slr2

        read (unit=qgrid_unit2) quad_lr(:)
        read (unit=qgrid_unit2) quad_lth(:)
        read (unit=qgrid_unit2) quad_lph(:)
        read (unit=hist_unit2) nconstrl(2), nvhist(2), dvhist(2)

        if (nvhist(1) /= nvhist(2) .or. dvhist(1) /= dvhist(2)) &
            stop " Velocity histograms of both libs are inconsistent"

        norbtotal = 2*norb(1) + norb(2)
    ! in case of rotating frame box orbit library is our retrogtete orbits (we do not flip the tubr orbits)
    if (Omega /= 0.0_dp ) then          !  (BT)
      norbtotal=norb(1)+norb(2)
    endif

  orboffset=0

  allocate(velhist  (-nvhist(1):nvhist(1),nconstr),&
       veltmp   (-nvhist(1):nvhist(1)        ),&
       veltmp1  (-nvhist(1):nvhist(1)        ),&
       velstore (-nvhist(1):nvhist(1),nconstr))

  velstore  (:,:)   = 0.0_dp
  velhist   (:,:)   = 0.0_dp
  quad_store(:,:,:,:) = 0.0_dp

        do i = 1, 2 ! loop over orbit libraries
            print *, " Reading orbitlibrary :", i

            allocate (orbtypes(ndith(i)**3))

            do j = 1, norb(i) ! loop over orbits
                if (i == 1) then
                    read_from = qgrid_unit1
                else
                    read_from = qgrid_unit2
                end if
                read (unit=read_from) t1, t2, t3, t4
                if (t1 /= j) stop " orbit number does not match"
                read (unit=read_from) orbtypes(:)
                read (unit=read_from) quad_light(:, :, :, :)
                velhist(:, :) = 0.0_dp

                if (i == 1) then
                    read_from = hist_unit1
                else
                    read_from = hist_unit2
                end if
                do k = 1, nconstr ! loop over apertures
                    read (unit=read_from) ivmin, ivmax
                    if (ivmin <= ivmax) &
                        read (unit=read_from) velhist(ivmin:ivmax, k)
                end do

        ! Loop over orbit with reverse angular momentum.
        ! Do not loop when i=2 (=second orbitlibray)
        !do k=1,-1+(i-1)*2,-2

         if (Omega /= 0.0_dp ) then              !  (BT)
            l1=1
            l2=1
            l3=1
         else
            l1=1
            l2=-1+(i-1)*2
            l3=-2
         end if

         do k=l1,l2,l3                              !  (BT)
           ! do not store twice if all the orbits are a box orbits
           !if (k==1 .or. (k==-1 .and. minval(orbtypes) < 4)) then
           orboffset=orboffset+1
           ! Check if we cruise the orbits the same way:
           if (orbinmat(orboffset) /=(j + (i-1) * norb(1)) * k) &
                stop "  ?!? Reading different orbit library?"

        ! In the normal case all orbits qualify.
         orb_qualify = 1
        select case (select)
           case (1)
               if (orbtype(orboffset) /= 1) orb_qualify = 0
           case (3)
      		   if (orbtype(orboffset) /= 3) orb_qualify = 0
		   case (4)
      		   if (orbtype(orboffset) /= 4) orb_qualify = 0
		   case (31)
      		   if (orbtype(orboffset) /= 3 .or. k == 1) orb_qualify = 0
		   case (11)
      		   if (orbtype(orboffset) /= 1 .or. k == 1) orb_qualify = 0
		   case (32)
      		   if (orbtype(orboffset) /= 3 .or. k == -1) orb_qualify = 0
		   case (12)
      		   if (orbtype(orboffset) /= 1 .or. k == -1) orb_qualify = 0
		 case default
               orb_qualify = 1
        end select

           ! Store data if this orbit has weight and is selected
           if (orbweight(orboffset) > 0.0_dp .and. orb_qualify == 1) then

              ! add the total light and flip sign of vx,vy,vz
              ! (light,x,y,z,vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx)
              quad_store(1 ,:,:,:)=   quad_store(1 ,:,:,:)+ orbweight(orboffset)*quad_light(1 ,:,:,:)
              quad_store(2 ,:,:,:)=   quad_store(2 ,:,:,:)+ orbweight(orboffset)*quad_light(2 ,:,:,:)*quad_light(1,:,:,:)
              quad_store(3 ,:,:,:)=   quad_store(3 ,:,:,:)+ orbweight(orboffset)*quad_light(3 ,:,:,:)*quad_light(1,:,:,:)
              quad_store(4 ,:,:,:)=   quad_store(4 ,:,:,:)+ orbweight(orboffset)*quad_light(4 ,:,:,:)*quad_light(1,:,:,:)
              quad_store(5 ,:,:,:)=   quad_store(5 ,:,:,:)+ orbweight(orboffset)*quad_light(5 ,:,:,:)*quad_light(1,:,:,:)*k
              quad_store(6 ,:,:,:)=   quad_store(6 ,:,:,:)+ orbweight(orboffset)*quad_light(6 ,:,:,:)*quad_light(1,:,:,:)*k
              quad_store(7 ,:,:,:)=   quad_store(7 ,:,:,:)+ orbweight(orboffset)*quad_light(7 ,:,:,:)*quad_light(1,:,:,:)*k
              quad_store(8 ,:,:,:)=   quad_store(8 ,:,:,:)+ orbweight(orboffset)*quad_light(8 ,:,:,:)*quad_light(1,:,:,:)
              quad_store(9 ,:,:,:)=   quad_store(9 ,:,:,:)+ orbweight(orboffset)*quad_light(9 ,:,:,:)*quad_light(1,:,:,:)
              quad_store(10,:,:,:)=   quad_store(10,:,:,:)+ orbweight(orboffset)*quad_light(10,:,:,:)*quad_light(1,:,:,:)
              quad_store(11,:,:,:)=   quad_store(11,:,:,:)+ orbweight(orboffset)*quad_light(11,:,:,:)*quad_light(1,:,:,:)
              quad_store(12,:,:,:)=   quad_store(12,:,:,:)+ orbweight(orboffset)*quad_light(12,:,:,:)*quad_light(1,:,:,:)
              quad_store(13,:,:,:)=   quad_store(13,:,:,:)+ orbweight(orboffset)*quad_light(13,:,:,:)*quad_light(1,:,:,:)
              quad_store(14,:,:,:)=   quad_store(14,:,:,:)+ orbweight(orboffset)*quad_light(14,:,:,:)
              quad_store(15,:,:,:)=   quad_store(15,:,:,:)+ orbweight(orboffset)*quad_light(15,:,:,:)
              quad_store(16,:,:,:)=   quad_store(16,:,:,:)+ orbweight(orboffset)*quad_light(16,:,:,:)

              do l=1,nconstr ! loop over aperture

                 ! Flip sign of velocity field of nessecary
                 if (k==1) then
                    veltmp(:)=velhist(:,l)
                 else
                    do m=-nvhist(i),nvhist(i)
                       veltmp(m)=velhist(-m,l)
                    end do
                 end if

                 velstore(:,l) = velstore(:,l) + veltmp(:)*orbweight(orboffset)

              end do
           end if
       end do
    end do
    deallocate(orbtypes)
 end do

        close (unit=qgrid_unit1)
        if (qgrid_unit1 /= hist_unit1) close (unit=hist_unit1)
        close (unit=qgrid_unit2)
        if (qgrid_unit2 /= hist_unit2) close (unit=hist_unit2)

print*, "  done reading orblibs"

! normalize quad_light
where (quad_store(1,:,:,:) /= 0.0_dp)
!  !  (light,x,y,z,vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx)
   quad_store( 2,:,:,:)=quad_store( 2,:,:,:)           /quad_store(1,:,:,:)/conversion_factor
   quad_store( 3,:,:,:)=quad_store( 3,:,:,:)           /quad_store(1,:,:,:)/conversion_factor
   quad_store( 4,:,:,:)=quad_store( 4,:,:,:)           /quad_store(1,:,:,:)/conversion_factor
   quad_store( 5,:,:,:)=quad_store( 5,:,:,:)*velunit   /quad_store(1,:,:,:)
   quad_store( 6,:,:,:)=quad_store( 6,:,:,:)*velunit   /quad_store(1,:,:,:)
   quad_store( 7,:,:,:)=quad_store( 7,:,:,:)*velunit   /quad_store(1,:,:,:)
   quad_store( 8,:,:,:)=quad_store( 8,:,:,:)*velunit**2/quad_store(1,:,:,:)
   quad_store( 9,:,:,:)=quad_store( 9,:,:,:)*velunit**2/quad_store(1,:,:,:)
   quad_store(10,:,:,:)=quad_store(10,:,:,:)*velunit**2/quad_store(1,:,:,:)
   quad_store(11,:,:,:)=quad_store(11,:,:,:)*velunit**2/quad_store(1,:,:,:)
   quad_store(12,:,:,:)=quad_store(12,:,:,:)*velunit**2/quad_store(1,:,:,:)
   quad_store(13,:,:,:)=quad_store(13,:,:,:)*velunit**2/quad_store(1,:,:,:)

end where

! write aperture mass

select case (select)
   case (1)
       outroot_kinem=trim(outroot)//'_t1'
   case (3)
	   outroot_kinem=trim(outroot)//'_t3'
   case (4)
	   outroot_kinem=trim(outroot)//'_t4'
   case (31)
	   outroot_kinem=trim(outroot)//'_t3p'
   case (11)
	   outroot_kinem=trim(outroot)//'_t1p'
   case (32)
	   outroot_kinem=trim(outroot)//'_t3n'
   case (12)
	   outroot_kinem=trim(outroot)//'_t1n'
 case default
       outroot_kinem=trim(outroot)
end select

print*,"writing "//trim(outroot_kinem)//"intrinsic_moments.out"

open (unit=30, file= trim(outroot_kinem)//"_intrinsic_moments.out", status="replace",&
     action="write")
write (unit=30, fmt=*) "smom1,sph1,sth1,slr1"
write (unit=30, fmt=*) smom1,sph1,sth1,slr1
write (unit=30, fmt=*) "phi boundaries"
write (unit=30, fmt="(30es13.5)") quad_lph
write (unit=30, fmt=*) "theta boundaries"
write (unit=30, fmt="(30es13.5)") quad_lth
write (unit=30, fmt=*) "radius boundaries in arcsec"
write (unit=30, fmt="(30es13.5)") quad_lr/conversion_factor
write (unit=30, fmt="(A120)") 'phi,theta,r,mgemass,fitmass,errmass,x,y,z '
write (unit=30, fmt="(A120)") '(in arcsec),vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx,orbit long, orbit short,boxes'

do i=1,sph1
   do j=1,sth1
      do k=1,slr1
         !  (light,x,y,z,vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx)
         write (unit=30, fmt="(3i8,30es13.5)") i,j,k,intrinsicmass(i,j,k),&
              quad_store(1,i,j,k),intrinsicmass(i,j,k)*intrmasserror,&
         (quad_store(l,i,j,k),l=2,smom1)
      end do
   end do
end do
close (unit=30)

do L=1,nconstr
   fitapermass(L)=sum(velstore(:,l))
   ! do not use the outer bins for the gh fit, since they contain
   ! the light outside the vel histogram
   veltmp1(-nvhist(1):nvhist(1))=velstore(-nvhist(1):nvhist(1),l)
   veltmp1(-nvhist(1))=0.0_dp
   veltmp1( nvhist(1))=0.0_dp
   !veltmp1(-nvhist(1)+1:nvhist(1)-1)=velstore(-nvhist(1)+1:nvhist(1)-1,l)
   ! normalize the histogram
   veltmp1(:)=veltmp1(:)/(sum(velstore(:,l))*dvhist(1))

   !call gaussfit(veltmp1(:),nvhist(1)-1,nvhist(1)-1,dvhist(1), gam,vm,sg,h12,5)
   call gaussfit(veltmp1(:),nvhist(1),nvhist(1),dvhist(1), gam,vm,sg,h12,5)
   ! Calculate the Gauss-Hermite moments for the observed mean velocity
   ! and velocity dispersion. in this calculation the value of gamma is
   ! kept equal to unity. The routine (automatically) returns zero for
   ! all gauss-Hermite moments if there is no weight in the histogram
   ! at all.
   call getgauher (vm,sg,veltmp1(:),nvhist(1),nvhist(1),&
        dvhist(1),hh,hermax,1.0_dp,1)

!   call getgauher (vm,sg,veltmp1(:),nvhist(1)-1,nvhist(1)-1,&
!        dvhist(1),hh,hermax,1.0_dp,1)
!   print*,hh

   if (hermax >= 1) modvel(l,1       )=vm * velunit
   if (hermax >= 2) modvel(l,2       )=sg * velunit
   if (hermax >= 3) modvel(l,3:hermax)=hh(3:hermax)

   ! Loop over the Gauss-Hermite moments, and store in the matrix

   ! Some attention is required here. To normalize the histogram w.r.t.
   ! integration over velocity one must divide by (totweight*dvhist).
   ! The matrix entry should be the Gauss-Hermite moment, multiplied by
   ! the fraction of the time the orbit is seen at the current aperture.
   ! This is totweight. Hence, to get the correct matrix entry, we
   ! should still divide by dvhist.

end do

deallocate(quad_lr,quad_lth,quad_lph)

print*,"writing kinem.out"

open (unit=30, file= trim(outroot_kinem)//"_kinem.out", status="replace",&
     action="write")
write (unit=30, fmt=*) nconstr,hermax
do i=1,nconstr
   write (unit=30, fmt="(i8,30es13.5)") i,apermass(i),fitapermass(i),&
        (modvel(i,j),velmom(i,j),dvelmom(i,j),j=1,hermax)
end do
close (unit=30)

open (unit=30, file= trim(outroot_kinem)//"_aphist.out", status="replace",&
     action="write")
write (unit=30, fmt="(2i8,es13.5)") nvhist(1),nconstr,dvhist(1)
write (unit=30, fmt="(6es13.5)") ((velstore(i,j),i=-nvhist(1),nvhist(1)),j=1,nconstr)
close (unit=30)

deallocate(velhist,veltmp,velstore)

deallocate(quad_light)

end subroutine makekinem

subroutine stopnnlsfit()
  deallocate(orbweight,orbint)
  !X Add all global allocatable arrays
end subroutine stopnnlsfit

subroutine reg_forcetoabel()
  real (kind=dp) :: rnorm,t1
  integer (kind=i4b) :: i,mi3,iE,iI2,iI3
  integer (kind=i4b),dimension(3) :: reg,ints
  integer (kind=i4b) :: nreg,row
  integer (kind=i4b),dimension(4,orbitsinfit) :: orbinttmp
  real    (kind=dp),dimension(7,orbitsinfit) :: abelorb

  if (use_reg) then

  print*,"Adding regularization"

  !stop ' FIXME: Read mass_radmass.dat file to get approx energy weights'

  if (reguldelta <= 0.0_dp) stop " zero regularization not implemented"

  ! add error to the whole range to avoid devide by zero.
  econ(nconstr*(hermax+1)+massconstr+1:convecl) = 1.0_dp

  !  store the orbint so we can change it temporarily
  orbinttmp=orbint

! real    (kind=dp),dimension(7,orbitsinfit) :: abelorb
open (unit=30,file="datfil/abel_orb.out", status="old", &
      action="read")
read (unit=30, fmt=*) i
read (unit=30, fmt=*) abelorb
close(unit=30)

!print*,abelorb(7,:)/sum(abelorb(7,:))
rnorm=sum(abelorb(7,:))

abelorb(6,:)=abelorb(7,:)/rnorm

rnorm=-1.0_dp!*reguldelta
  nreg=0

  do i=1,orbitsinfit
  if (abelorb(6,i) > rnorm) then
     nreg=nreg+1
     row = nconstr*(hermax+1)+massconstr + 1 + nreg

     if (row > convecl) stop " array to small for regularization"

     con(row)  = abelorb(6,i)
     econ(row) = maxval((/reguldelta,rnorm /))
     !print*,nreg,i,abelorb(6,i)
     orbmat(row,i) = 1.0_dp
  end if
  end do

  print*,nreg,sum(con(nconstr*(hermax+1)+massconstr + 2:row))

  ! restore orbint
  orbint=orbinttmp

  endif ! add regularisation
end subroutine reg_forcetoabel

subroutine donnls_abel()
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! Solve the NNLS problem for the orbital weights

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real (kind=dp), dimension(size(orbmat,1),size(orbmat,2)) :: orbtmp
  ! Working space and other arrays
  integer (kind=i4b) :: len, i, j, k
  real (kind=dp) :: tmpcon,tmpecon,t1
  integer (kind=i4b), dimension(orbitsinfit) :: index
  integer (kind=i4b) :: it, ior, mode,a,b,c,d
  real    (kind=dp),dimension(7,orbitsinfit) :: abelorb
  print*,"  * Solver function memory allocated."

  allocate(orbweight(orbitsinfit))

open (unit=30,file="datfil/abel_orb.out", status="old", &
      action="read")
read (unit=30, fmt=*) i
read (unit=30, fmt=*) abelorb
close(unit=30)

print*,sum(abelorb(7,:))
orbweight(:)=abelorb(7,:)!/sum(abelorb(7,:))


   print*,'  * Solving done'
end subroutine donnls_abel


end module NNLSfit

program triaxnnls
use nnlsfit
use omp_lib

integer  :: i

!call omp_set_num_threads(1)
call readobservations()
call readorblibs()
call add_regularization()
!call reg_forcetoabel() ! alternative regularisation

print*,"Which solver do you want?"
print*,"(0) Least squares with bounds and constraints FitZeroMoment"
print*,"(1) NNLS "
print*,"(2) NNLS and save the matrix to disk to save memory"
print*,"(3) Read solution from ./datfil/abel_orb.out"
print*,"(4) Least squares with bounds and constraints savetodisk,FitZeroMoment"
print*,"(5) Least squares with bounds and constraints ConstrainZeroMoment (Preferred)"
print*,"(6) Least squares with bounds and constraints savetodisk,ConstrainZeroMoment"
i=0
read (unit=*,fmt=*) i
print*,i
if (i < 0 .or. i > 6) stop " Input not valid."

print*,"Starting Solver"

if (i == 0) call donnls_galahad(0,1)
if (i == 1) call donnls_nosave()
if (i == 2) call donnls()
if (i == 3) call donnls_abel()
if (i == 4) call donnls_galahad(1,1)
if (i == 5) call donnls_galahad(0,0)
if (i == 6) call donnls_galahad(1,0)

call allpred()
call makekinem(0)
!call makekinem(1)
!call makekinem(3)
!call makekinem(4)
!call makekinem(11)
!call makekinem(12)
!call makekinem(31)
!call makekinem(32)

end program triaxnnls
