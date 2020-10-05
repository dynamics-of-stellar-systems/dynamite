module orbsolread
use numeric_kinds
private
public :: readorb


contains

subroutine readorb()
use output , only: out_file
use initial_parameters    , only : conversion_factor
integer(kind=i4b) :: i,j,k,hdl,ni1,ni2,ni3,ndith,slr,sth,sph,nconstr,i1,i2,iap
integer(kind=i4b) :: nvhist,smom
real(kind=dp) :: t1,t2,t3,t4,histwidth
integer(kind=i4b) :: ivmin,ivmax
integer(kind=i4b),dimension(:),allocatable :: orbtypes
real(kind=dp),dimension(:),allocatable :: quad_lr,quad_lth,quad_lph,veltmp1
real(kind=dp),dimension(:,:,:,:),allocatable :: quad_light
real(kind=dp),dimension(:,:),allocatable :: velhist
character (len=256) :: outroot

!real(kind=dp),dimension(:,:,:),allocatable :: quad_light

! Define the outout files name
i=LEN_TRIM(out_file)  
if (i-4 < 1) stop "output filename too short"
outroot=out_file(1:i-4)

hdl=27
open (unit=hdl,file=trim(outroot)//".dat",action="read", &
                status="old",form="unformatted",position="rewind")
print*,"opened"
read (unit=hdl) i,ni1,ni2,ni3 ,ndith
print*,i,ni1,ni2,ni3,ndith
allocate(orbtypes(ndith))
read (unit=hdl) smom,sph,sth,slr
print*,smom,slr,sth,sph
  ! remember that N bins have N+1 boundaries
allocate(quad_lr(slr+1),quad_lth(sth+1),quad_lph(sph+1),quad_light(smom,sph,sth,slr))
read (unit=hdl) quad_lr (:)
read (unit=hdl) quad_lth(:)
read (unit=hdl) quad_lph(:)
read (unit=hdl) nconstr, nvhist, histwidth
print*, nconstr, nvhist, histwidth
!nvhist=1
 allocate(velhist(-nvhist:nvhist,nconstr),veltmp1(-nvhist:nvhist))

i=1
print*,i
velhist(:,:) = 0.0_dp

read (unit=hdl) i1,ni1,ni2,ni3
read (unit=hdl) orbtypes(:)  
read (unit=hdl) quad_light(:,:,:,:)
!print*,quad_light(1,:,:,:)
do iap=1,nconstr
   read (unit=hdl) ivmin, ivmax
   !print*,ivmin,ivmax
   if (ivmin <= ivmax) &
        read (unit=hdl) velhist(ivmin:ivmax,iap)
end do
close(unit=hdl)

! write aperture mass

print*,"  * Writing aper_mass"
open (unit=31,file=trim(outroot)//"_apermass.out",action="write", status="replace")
write (unit=31, fmt="(i5)") nconstr
do i=1,nconstr
   write (unit=31, fmt="(i5,es20.12)") i,sum(velhist(:,i))
   !write (unit=31, fmt="(es20.12)") sum(velhist(:,i))
end do
close(unit=31)


!open (unit=28,file=trim(outroot)//".out",action="write", &
!                status="replace")
!write(unit=28,fmt=*) nvhist*2+1,nconstr, histwidth
!write(unit=28,fmt=*) velhist
!close(unit=28)


! converting units from km to arcsec
quad_lr=quad_lr/conversion_factor
quad_light(2:4,:,:,:)=quad_light(2:4,:,:,:)/conversion_factor

! write aperture mass

print*,"writing intrinsic_moments.out"

open (unit=30, file= trim(outroot)//"_intrinsic_moments.out", status="replace",&
     action="write")
write (unit=30, fmt=*) "smom,sph,sth,slr"
write (unit=30, fmt=*) smom,sph,sth,slr
write (unit=30, fmt=*) "phi boundaries"
write (unit=30, fmt="(30es13.5)") quad_lph 
write (unit=30, fmt=*) "theta boundaries"
write (unit=30, fmt="(30es13.5)") quad_lth 
write (unit=30, fmt=*) "radius boundaries in arcsec"
write (unit=30, fmt="(30es13.5)") quad_lr/conversion_factor
write (unit=30, fmt=*) "phi,theta,r, 0.0, mass,0.0,x,y,z (in arcsec),vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx"

do i=1,sph
   do j=1,sth
      do k=1,slr
         !  (light,x,y,z,vx,vy,vz,xv2,vy2,vz2,vxvy,vyvz,vzvx)
         write (unit=30, fmt="(3i8,30es13.5)") i,j,k,0.0_dp &
              , quad_light(1,i,j,k),0.0_dp &
              ,(quad_light(l,i,j,k),l=2,smom)
      end do
   end do
end do
close (unit=30)

open (unit=30, file= trim(outroot)//"_aphist.out", status="replace",&
     action="write")
write (unit=30, fmt="(2i8,es13.5)") nvhist,nconstr,histwidth
write (unit=30, fmt="(6es13.5)") ((velhist(i,j),i=-nvhist,nvhist),j=1,nconstr)
close (unit=30)

deallocate(velhist,veltmp1)

deallocate(quad_light)


end subroutine readorb

end module orbsolread



module Orbitcombine
  use numeric_kinds
  implicit none
  private

  integer (kind=i4b),private,dimension(:,:),allocatable :: orbinfo

! setup/run/stop the program.
public :: setup,run,stob

private :: read_orbsol
contains

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine setup()
    use integrator         , only : integrator_setup
    use projection         , only : projection_setup,projection_change_direction
    use quadrantgrid       , only : qgrid_setup
    use aperture_routines  , only : aperture_setup
    use histograms         , only : histogram_setup
    use psf                , only : psf_setup
    use output             , only : output_setup
  !----------------------------------------------------------------------
    character(len = 80) :: string
    print *,"  ** Start Orbit combine Setup"
    print *,"  * Give setup version info: [U for unspecified]"
    read *, string
    if (string == "#Orbitcombine_setupfile_version_2" .or. string=="U") then
       print*,"  * Setupfile is Version 1"
       call integrator_setup()
       call projection_setup()
       call projection_change_direction()
       call qgrid_setup()
       call psf_setup()
       call aperture_setup()
       call histogram_setup()
       call output_setup()
       call read_orbsol_relative_precision()
    else
       print*, "This version is not understood by this program"
       STOP "program terminated in high_level:setup"
    end if

    print *,"  ** Setup Finished"

  end subroutine setup


subroutine read_orbsol ()
  use integrator  , only : integrator_points,integrator_dithering
  character (len=30) :: infil
  integer (kind=i4b)  :: i,j,norb,nparticles 
  real    (kind=dp )  :: precision,t1,t2

  real (kind=dp)    ,dimension(:),allocatable :: orbwght

  print*, "  * Give name of file with the orbitsolution (*_orb.out)"
  read (unit=*, fmt="(a30)") infil
  print*,"   ",infil

  open (unit=40,file=infil, status="OLD", action="read",position="rewind")
  read (unit=40, fmt=*) norb
  allocate(orbinfo(6,norb),orbwght(norb))
  do i=1,norb
     read (unit=40, fmt="(6i8,es13.5)") j, orbinfo(1,i),orbinfo(2,i),orbinfo(3,i),orbinfo(4,i),orbinfo(5,i),orbwght(i)
     !orbint(1:4,i),orbtype(i),orbweight(i)
  end do
  close (unit=40)

 
  print*, "  * max amount of particles that could be generated (int64):", huge(1_i4b)
  print*, "  * Give the number of particles to generate:"
  read *, nparticles
  print*,"   ",nparticles

  if (nparticles <= 1_i4b) stop 

  orbinfo(6,:)=nint(orbwght(:)*nparticles/integrator_dithering**3)

 
  !t2=sum(orbwght(:))/precision*integrator_points
  !t1=t2-sum(dble(orbinfo(6,:)))
  t1=nparticles - sum (orbinfo(6,:))*(integrator_dithering**3) 

  print*,"  * Total orbweight (~1.00):",sum(orbwght(:))
  print*,"  * Total number of photons:", sum(dble(orbinfo(6,:)))*integrator_dithering**3
	print*,t1
  print*,"  * Relative Precision of the model :",(nparticles-t1)/t1
  !if (abs(t1/t2*integrator_dithering**3*8.0_dp) > 0.01) stop "precison less than 1 procent"
  print*,"  * Absolute number of photons difference:",t1
  print*,"  * symmetry number:",(integrator_dithering**3* 8.0)
  !rint*,(orbwght(:)/precision*integrator_points)*integrator_dithering**3

  !print*,epsilon(sum(dble(orbinfo(6,:)))*integrator_dithering**3*8.0_dp),exponent(sum(dble(orbinfo(6,:)))*integrator_dithering**3*8.0_dp)

  deallocate(orbwght)
end subroutine read_orbsol


subroutine read_orbsol_relative_precision ()
  use integrator  , only : integrator_points,integrator_dithering
  character (len=30) :: infil
  integer (kind=i4b)  :: i,j,norb
  real    (kind=dp )  :: precision,t1,t2,maxpre

  real (kind=dp)    ,dimension(:),allocatable :: orbwght

  print*, "  * Give name of file with the orbitsolution (*_orb.out)"
  read (unit=*, fmt="(a30)") infil
  print*,"   ",infil

  open (unit=40,file=infil, status="OLD", action="read",position="rewind")
  read (unit=40, fmt=*) norb
  allocate(orbinfo(6,norb),orbwght(norb))
  do i=1,norb
     read (unit=40, fmt="(6i8,es13.5)") j, orbinfo(1,i),orbinfo(2,i),orbinfo(3,i),orbinfo(4,i),orbinfo(5,i),orbwght(i)
     !orbint(1:4,i),orbtype(i),orbweight(i)
  end do
  close (unit=40)

  maxpre=maxval(orbwght(:))/real(huge(i))*integrator_points
  print*, "  * Maximum percision is:",maxpre
  read *, precision
  print*,"   ",precision
 
  if (precision <= maxpre) stop 


  t1=sum(orbwght(:))*real(huge(i))*integrator_points*integrator_dithering**3*8.0_dp
  t2=sum(orbwght(:))*real(huge(i))*integrator_points*integrator_dithering**3*8.0_dp+1.0_dp


  orbinfo(6,:)=nint(orbwght(:)/precision*integrator_points)


  print*,"Check to see if sumation works",t2,t2,t2-t1,t2-2.0_dp-t1
  t2=sum(orbwght(:))/precision*integrator_points
  t1=t2-sum(dble(orbinfo(6,:)))

  
  print*,"  * Total orbweight (~1.00):",sum(orbwght(:))
  print*,"  * Total number of photons:", sum(dble(orbinfo(6,:)))*integrator_dithering**3*8.0_dp
  print*,"  * Relative Precision of the model :",t1/t2*integrator_dithering**3*8.0_dp
  !if (abs(t1/t2*integrator_dithering**3*8.0_dp) > 0.01) stop "precison less than 1 procent"
  print*,"  * Absolute number of photons difference:",t1

  !rint*,(orbwght(:)/precision*integrator_points)*integrator_dithering**3

  !print*,epsilon(sum(dble(orbinfo(6,:)))*integrator_dithering**3*8.0_dp),exponent(sum(dble(orbinfo(6,:)))*integrator_dithering**3*8.0_dp)

  deallocate(orbwght)
 !  stop
end subroutine read_orbsol_relative_precision

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine run()
    use histograms  , only : histogram_reset,hist_thesame, &
         histogram_velbin,histogram_store
    use projection  , only : project,projection_symmetry
    use integrator  , only : integrator_integrate,integrator_points,vyini,&
         integrator_current,vyini,ini_integ,nEner,nI2,nI3,integrator_dithering
    use output      , only : output_write
    use quadrantgrid, only : qgrid_reset,qgrid_store
    use psf         , only : psf_n,psf_gaussian
    use aperture         , only : aperture_n,aperture_psf
    use aperture_routines, only : aperture_find
    use output , only: out_file
    !----------------------------------------------------------------------
    logical :: done,first,alldone
    real (kind=dp)    ,dimension(:,:),allocatable :: pos,vel
    real (kind=dp)    ,dimension(:,:),allocatable :: proj,vec_gauss
    real (kind=dp)    ,dimension(:  ),allocatable :: losvel
    integer (kind=i4b),dimension(:  ),allocatable :: velb,poly	
    integer (kind=i4b)                            :: ap,i,pf,no
    integer (kind=i4b)  :: type
    real    (kind=dp )  :: t1,t2
    character (len=256) :: outroot
    !print*,orbinfo(6,:)

    call histogram_reset()
    call qgrid_reset()      
	i=LEN_TRIM(out_file)  
	if (i-4 < 1) stop "output filename too short"
    outroot =      out_file(1:i-4)
!	print*," writing particles file: ",trim(outroot)//"_particles.out"
!	open (unit=33,file=trim(outroot)//"_particles.out",action="write", status="replace")
	
	  no=(nEner*nI2*nI3/integrator_dithering**3)
    print*,no
    call compute_orbits(1,2,no*2)
    print*,"  * Starting second set of orbits"
    vyini(:)=-vyini(:)
    call compute_orbits(2,2,no*2)

    print*,"  * Starting box orbits"

    call ini_integ()
    call compute_orbits(no*2+1,1,no*3)
    close(unit=33)

    call output_write()

    print*,"  ** Finished Orbit Calculations"

  end subroutine run

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine compute_orbits(orbstart,orbdelta,orbend)
    use histograms  , only : histogram_reset,hist_thesame, &
         histogram_velbin,histogram_store
    use projection  , only : project,projection_symmetry
    use integrator  , only : integrator_integrate,integrator_points,vyini,&
         integrator_current,nEner,integrator_set_current
    use output      , only : output_write
    use quadrantgrid, only : qgrid_reset,qgrid_store
    use psf         , only : psf_n,psf_gaussian
    use aperture         , only : aperture_n,aperture_psf
    use aperture_routines, only : aperture_find

    !----------------------------------------------------------------------
    integer (kind=i4b),intent(in) :: orbstart,orbdelta,orbend
    logical :: done,first,alldone
    real (kind=dp)    ,dimension(:,:),allocatable :: pos,vel
    real (kind=dp)    ,dimension(:,:),allocatable :: proj,vec_gauss
    real (kind=dp)    ,dimension(:  ),allocatable :: losvel
    integer (kind=i4b),dimension(:  ),allocatable :: velb,poly	
    integer (kind=i4b)                            :: ap,i,pf,orbn
    integer (kind=i4b)  :: type
    real    (kind=dp )  :: t1,t2

    orbn=orbstart
    integrator_current=0
    alldone=.false.
    print*,"  ** Starting Orbit Calculations"
    do  ! for each orbit
      
       call cpu_time(t1)
       
       integrator_points=orbinfo(6,orbn)
       orbn=orbn+orbdelta

       integrator_current=integrator_current+1
       if (integrator_points > 0) then 
          print*," This",(orbn-orbstart)/orbdelta," orbit has ",integrator_points,"photons"
          print*,orbn,integrator_current
       allocate (pos(integrator_points,3),vel(integrator_points,3))
       allocate (proj(integrator_points*projection_symmetry,2))
       allocate (vec_gauss(integrator_points*projection_symmetry,2))
       allocate (losvel(integrator_points*projection_symmetry))
       allocate (velb(integrator_points*projection_symmetry))
       allocate (poly(integrator_points*projection_symmetry))
       first=.true.
       call integrator_set_current(integrator_current-1)
       do ! for all dithers
          call integrator_integrate(pos,vel,type,done,first,alldone)
          first=.false.
          if (done .or. alldone) exit                        
		  do i=1,integrator_points 
 !           write (unit=33, fmt="(6es20.12)") pos(i,:),vel(i,:)    
          end do
          call qgrid_store(pos(:,:),vel(:,:),type)
          first = .true.
          do ! for all projections
             call project(type,pos,vel,proj,losvel,done,first)
             if (done) exit
             first = .false.

             if (hist_thesame) call histogram_velbin(1,losvel,velb)
             do i=1,psf_n
                if (.not. hist_thesame) call histogram_velbin(i,losvel,velb)
                   call psf_gaussian(i,proj,vec_gauss)
                do ap=1,aperture_n
                   if ( i == aperture_psf(ap) ) then
                      call aperture_find  (ap,vec_gauss,poly)
                      call histogram_store(ap,poly,velb,size(proj,1))
                   end if
                end do
             end do
          end do
       end do
       deallocate (pos,vel,proj,vec_gauss,losvel,velb,poly)
       call cpu_time(t2)
       print*,"  * Time spent one orbit:",t2-t1," seconds"
       end if
       if (alldone .or. orbn > orbend) exit
    end do

    print*,"  ** Finished Orbit Calculations"

  end subroutine compute_orbits
 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine stob()
    use integrator        , only : integrator_stop
    use projection        , only : projection_stop
    use histograms        , only : histogram_stop
    use psf               , only : psf_stop
    use aperture_routines , only : aperture_stop
    use output            , only : output_close
    use orbsolread        , only : readorb
  !----------------------------------------------------------------------

    call output_close    ()

    call readorb()

    call integrator_stop ()
    call projection_stop ()
    call aperture_stop   ()
    call psf_stop        ()
    call histogram_stop  ()

    deallocate(orbinfo)

    
  end subroutine stob

end module Orbitcombine

!######################################################################
!######################################################################
!######################################################################

program counterrotation
  use numeric_kinds
  use Orbitcombine, only : setup,run,stob
  implicit none
  real (kind=dp) :: t1,t2
  print*,"  ** Triaxial Orbit library by Remco C.E. van den Bosch <bosch@strw.leidenuniv.nl>"
  print*,"  * Jan 2007 "
  print*,"  * $Id: modelgen.f90,v 1.1 2010/03/17 12:56:57 bosch Exp $"  

  call setup()
  print*,"Starting orbit integration routines"
  call cpu_time(t1)
  call run()
  call cpu_time(t2)
  call stob()
  print*,"  ** Program Finished"
  print*,"  * Time spent calculating :",t2-t1," seconds"
end program counterrotation

!######################################################################
!######################################################################
!######################################################################
 
