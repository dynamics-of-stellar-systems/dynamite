
module Orbitcombine
  use numeric_kinds
  implicit none
  private

  integer (kind=i4b),private,dimension(:,:),allocatable :: orbinfo
  real (kind=dp)    ,dimension(:),allocatable :: orbwght      
  character (len=256) :: outroot      
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
       !call projection_change_direction()
       !call qgrid_setup()
       !call psf_setup()
       !call aperture_setup()
       !call histogram_setup()
       !call output_setup()  
	   print*,'bla'
       call read_orbsol()
    else
       print*, "This version is not understood by this program"
       STOP "program terminated in high_level:setup"
    end if

    print *,"  ** Setup Finished"

  end subroutine setup


subroutine read_orbsol ()
  use integrator  , only : integrator_points,integrator_dithering,rcirc     
  use initial_parameters, only : totalmass,conversion_factor
  character (len=30) :: infil
  integer (kind=i4b)  :: i,j,norb,nparticles 
  real    (kind=dp )  :: precision,t1,t2,rem,t3
  real (kind=dp)    ,dimension(:),allocatable :: massscale

  print*, "  * Give name of file with the orbitsolution (*_orb.out)"
  read (unit=*, fmt="(a256)") infil
  print*,"   ",infil

  open (unit=40,file=infil, status="OLD", action="read",position="rewind")
  read (unit=40, fmt=*) norb
  allocate(orbinfo(6,norb),orbwght(norb),massscale(norb))
  do i=1,norb
     !read (unit=40, fmt="(6i8,es13.5)") j, orbinfo(1,i),orbinfo(2,i),orbinfo(3,i),orbinfo(4,i),orbinfo(5,i),orbwght(i)
     read (unit=40, fmt=*) j, orbinfo(1,i),orbinfo(2,i),orbinfo(3,i),orbinfo(4,i),orbinfo(5,i),orbwght(i)
     !orbint(1:4,i),orbtype(i),orbweight(i)
  end do
  close (unit=40)   

   i=LEN_TRIM(infil)  
   if (i-4 < 1) stop "output filename too short"
   outroot =      infil(1:i-4)
 
  print*, "  * max amount of particles that could be generated (int64):", huge(1_i4b)
  print*, "  * Give the number of particles to generate:"
  read *, nparticles
  print*,"   ",nparticles

  if (nparticles <= 1_i4b) stop 
 
  print*,size(rcirc) ,norb/3 ,(norb/3*(integrator_dithering**3))
  !massscale(1:norb/3)= &
  !  rcirc(125/2:(norb/3*(integrator_dithering**3)):integrator_dithering**3) &
  !  /conversion_factor
  !massscale(Norb/3  +1:Norb/3*2)=massscale(1:norb/3)
  !massscale(Norb/3*2+1:Norb)=massscale(1:norb/3)
  !massscale(:)=(massscale(:)+0.01)/(massscale(:)+100)     
  !t1=sum(massscale(:))      

  print*,'Equal mass particles'  
  massscale(:)=sum(orbwght(:)) !1.0_dp
  !t1=1.0_dp   
  t1=sum(orbwght(:)) 
 
 print*,"  * Total orbweight (~1.00): ",sum(orbwght(:))
  ! Nint rounds to nearest whole number
  !orbinfo(6,:)=nint(orbwght(:)*nparticles/massscale(:)*t1)
  orbinfo(6,:)=nint(orbwght(:)*nparticles/massscale(:)*t1)


  ! Enforce exactly equal mass particles, and add the remainder mass to the
  ! future orbits.
  rem=0.0_dp
  do i=1,norb
	 t3= orbinfo(6,i) * massscale(i) / nparticles /t1
     rem = rem + orbwght(i) - t3
     !if (rem .ge. 4.0_dp/nparticles/t1*massscale(i)) then
	 !  rem = rem - 1.0_dp/nparticles/t1*massscale(i)
	 !  orbinfo(6,i) = orbinfo(6,i) + 1 
	 !print*,i,rem,1.0_dp/nparticles
	 !endif
     orbwght(i) = orbinfo(6,i) * massscale(i) / nparticles /t1
  end do  
  print*,'remaining mass remainder',rem

  ! ensure all the mass in the particles
  !where ( orbinfo(6,:) .lt. 1 .and. orbwght(:) .gt. 1.0e-20)    orbinfo(6,:) = 1

  t2=0.0_dp
  do i=1,norb
     if ( orbinfo(6,i) .gt. 0 )  t2=t2+orbwght(i)
  end do

  !print*,massscale/conversion_factor
  print*,conversion_factor
  t1=nparticles - sum (orbinfo(6,:)) 

  print*,"  * Total orbweight (~1.00): ",sum(orbwght(:))
  print*,"  * Total mass in realization: ",t2
  print*,"  * relative Mass percision of the model: ",(t2-sum(orbwght))/sum(orbwght)    
  print*,"  * absolute Mass percision of the model: ",(t2-sum(orbwght))*totalmass
  print*,"  * Total number of photons:", sum(dble(orbinfo(6,:)))
  print*,"  * Relative Precision of the model :",(t1)/nparticles
  print*,"  * Absolute number of photons difference:",t1
  print*,"  * symmetry number:",(integrator_dithering**3* 8.0)  

  open (unit=33,file=trim(outroot)//"_part_masses.out",action="write", status="replace")
  do i=1, norb
  write (unit=33, fmt=*) orbwght(i), massscale(i), orbinfo(6,i)
  end do
  !stop
end subroutine read_orbsol



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
    use initial_parameters, only : xmbh,grav_const_km, conversion_factor 	
    !use output , only: out_file
    !----------------------------------------------------------------------
    logical :: done,first,alldone
    real (kind=dp)    ,dimension(:,:),allocatable :: pos,vel
    real (kind=dp)    ,dimension(:,:),allocatable :: proj,vec_gauss
    real (kind=dp)    ,dimension(:  ),allocatable :: losvel
    integer (kind=i4b),dimension(:  ),allocatable :: velb,poly	
    integer (kind=i4b)                            :: ap,i,pf,no
    integer (kind=i4b)  :: type
    real    (kind=dp )  :: t1,t2
 
	print*," writing particles file: ",trim(outroot)//"_particles.out"
	open (unit=33,file=trim(outroot)//"_particles.out",action="write", status="replace")

    ! Define units 
    write (unit=33, fmt="(A40,1es20.12)"),"Gravitional constant in km^3/(s^2 Msun):",   grav_const_km
    ! Parsec unit is not needed for N-body calculations. We only use it for MGE->pot
    !write (unit=33, fmt="(A20,1es20.12)"),"Parsec in km",   parsec_km    
    write (unit=33, fmt="(A40,1es20.12)"),"To convert from km to arcsec use:",1.0_dp/conversion_factor 

    write (unit=33,fmt="(A40)"), "First particle is the black hole"
    write (unit=33,fmt="(A40)"), "Rows: Mass, x, y, z, vx, vy, vz"
    write (unit=33,fmt="(A40)"), "Units: Msun, km, km, km, km/s, km/s, km/s"
    ! Black hole is the first particle
    write (unit=33, fmt="(7es20.12)") xmbh,0.0,0.0,0.0,0.0,0.0,0.0
 
	  no=(nEner*nI2*nI3/integrator_dithering**3)
    print*,no
    call compute_orbits(1,2,no*2)
    print*,"  * Starting second set of orbits"
    vyini(:)=-vyini(:)
    call compute_orbits(2,2,no*2)

    print*,"  * Starting box orbits"
    call ini_integ()
    call compute_orbits(no*2+1,1,no*3)  

    write (unit=33,fmt="(A20)"), "The END_OF_FILE"
    close(unit=33)

    !call output_write()

    print*,"  ** Finished Orbit Calculations"

  end subroutine run

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine compute_orbits(orbstart,orbdelta,orbend)
    use histograms  , only : histogram_reset,hist_thesame, &
         histogram_velbin,histogram_store
    use projection  , only : project,projection_symmetry
    use integrator  , only : integrator_integrate,integrator_points,vyini,&
         integrator_current,nEner,integrator_set_current,integrator_dithering
    use output      , only : output_write
    use quadrantgrid, only : qgrid_reset,qgrid_store
    use psf         , only : psf_n,psf_gaussian
    use aperture         , only : aperture_n,aperture_psf
    use aperture_routines, only : aperture_find
    use initial_parameters, only : totalmass
    !----------------------------------------------------------------------
    integer (kind=i4b),intent(in) :: orbstart,orbdelta,orbend
    logical :: done,first,alldone
    real (kind=dp)    ,dimension(:,:),allocatable :: pos,vel
    real (kind=dp)    ,dimension(:,:),allocatable :: postot,veltot
    real (kind=dp)    ,dimension(:,:),allocatable :: proj,vec_gauss
    real (kind=dp)    ,dimension(:  ),allocatable :: losvel
    integer (kind=i4b),dimension(:  ),allocatable :: velb,poly	
    integer (kind=i4b)                            :: ap,i,ri,pf,orbn,nphot ,j
    integer (kind=i4b)  :: type, indtot,index_project
    real    (kind=dp )  :: t1,t2,rand,weight

    orbn=orbstart
    integrator_current=0
    alldone=.false.
    print*,"  ** Starting Orbit Calculations"
    allocate (postot(integrator_points*8*integrator_dithering**3,3))
	allocate (veltot(integrator_points*8*integrator_dithering**3,3))
    allocate (pos(integrator_points,3),vel(integrator_points,3))
    allocate (proj(integrator_points*projection_symmetry,2))
    allocate (vec_gauss(integrator_points*projection_symmetry,2))
    allocate (losvel(integrator_points*projection_symmetry))
    allocate (velb(integrator_points*projection_symmetry))
    allocate (poly(integrator_points*projection_symmetry))

    do  ! for each orbit
      
       call cpu_time(t1)    
       indtot=0
       nphot=orbinfo(6,orbn)
       weight=orbwght(orbn)*totalmass
       orbn=orbn+orbdelta
       integrator_current=integrator_current+1
       if (nphot > 0) then 
       	print*," This",(orbn-orbstart)/orbdelta," orbit has ",nphot,"photons"
	    !print*,orbn,integrator_current
	   
	   first=.true.
       call integrator_set_current(integrator_current-1)
	   index_project=1_i4b   
		postot=0.0_dp
		veltot=0.0_dp
		j=1
       do ! for all dithers
          call integrator_integrate(pos,vel,type,done,first,alldone)
          first=.false.
          if (done .or. alldone) exit    
          if (j .eq. integrator_dithering**3/2+1) then          
           do i=1,integrator_points 
          write (unit=33, fmt="(7es20.12)") 1.0,pos(i,:),vel(i,:) 
           enddo                                     
          endif
		! project overwrite the same vector everytime
!          call project_part(type,pos,vel,postot,veltot,index_project)
           j=j+1
        end do
       
	
!	   ! Fisher-Yates shuffle  
!       do i=1,nphot  
!		   call random_number(rand)
!		   ! random number between 1 and integrator_points*8*integrator_dithering**3- (i-1) 
!		   ri=int(rand*(integrator_points*8*integrator_dithering**3- (i-1))) + 1
!		   write (unit=33, fmt="(7es20.12)") weight/nphot,postot(ri,:),veltot(ri,:)
!		   ! Swap elements to make ensure they do not get drawn again.
!		   postot(size(postot,1)-(i-1),:)=postot(ri,:) 
!		   veltot(size(postot,1)-(i-1),:)=veltot(ri,:)
!		end do



       call cpu_time(t2)
       print*,"  * Time spent one orbit:",t2-t1," seconds"
       end if
       if (alldone .or. orbn > orbend) exit
    end do           
    deallocate (postot,veltot,pos,vel,proj,vec_gauss,losvel,velb,poly)
    print*,"  ** Finished Orbit Calculations"

  end subroutine compute_orbits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine project_part(type,pos,vel,postot,veltot,ind)
  use initial_parameters, only : theta_view, phi_view
    real (kind=dp),intent(in ),dimension(:        ,:)   :: pos
    real (kind=dp),intent(in ),dimension(size(pos,1),3) :: vel
    real (kind=dp),intent(out),dimension(:,:)           :: postot
    real (kind=dp),intent(out),dimension(:,:)           :: veltot
    integer (kind=i4b),intent(in)                       :: type  
	integer (kind=i4b),intent(inout)                    :: ind
  !----------------------------------------------------------------------
    real (kind=dp)               :: t1, t2, t3,theta,phi
    integer (kind=i4b)          :: i ,k ,l
    ! Signs of the (vx,vy,vz) for each Projection and type of Orbit
    real (kind=dp),dimension(3,8,5),parameter :: vsgn= reshape((/  &
    ! X tubes
    1 , 1 , 1    ,-1 , 1 , 1  , -1 , 1 ,-1  ,  1 , 1 , -1 , & 
    1 ,-1 , 1    ,-1 ,-1 , 1  , -1 ,-1 ,-1  ,  1 ,-1 , -1 , & 
    ! Y tubes
    1 , 1 , 1    , 1 , 1 ,-1  ,  1 ,-1 ,-1  ,  1 ,-1 , 1 , & 
   -1 , 1 , 1    ,-1 , 1 ,-1  , -1 ,-1 ,-1  , -1 ,-1 , 1 , & 
    ! Z tubes
    1 , 1 , 1    , 1 ,-1 , 1  , -1 ,-1 , 1  , -1 , 1 , 1 , & 
    1 , 1 ,-1    , 1 ,-1 ,-1  , -1 ,-1 ,-1  , -1 , 1 ,-1 , & 
    ! Boxed 
    1 , 1 , 1    ,-1 , 1 , 1  , -1 ,-1 , 1  ,  1 ,-1 , 1 , & 
    1 , 1 ,-1    ,-1 , 1 ,-1  , -1 ,-1 ,-1  ,  1 ,-1 ,-1 , & 
    ! Stochastic
    1 , 1 , 1    ,-1 , 1 , 1  , -1 ,-1 , 1  ,  1 ,-1 , 1 , & 
    1 , 1 ,-1    ,-1 , 1 ,-1  , -1 ,-1 ,-1  ,  1 ,-1 ,-1 /),(/3,8,5/))

    !Signs of the x,y,z for each projection  :psgn( [x,y,z], project )  
    real (kind=dp),dimension(3,8),parameter :: psgn= reshape((/  &   
    1 , 1 , 1   , -1 , 1 , 1  , -1 , -1 , 1 ,  1 , -1 , 1 , & 
    1 , 1 ,-1   , -1 , 1 ,-1  , -1 , -1 ,-1 ,  1 , -1 ,-1 /),(/3,8/))
    
    ! check orbit type
    if (type > 5 .or. type < 1 ) stop "project_n: Wrong orbit type" 
  
    ! Use sign matrix for the symmetries.
    ! Using the inverse (transpose) of the projection (eq. 4) of Thesis Ellen.
    
    do i=1,8
	  do l=1,3
	    postot(ind:ind+size(pos,1)-1,l)=pos(:,l)*psgn(l,i)
	    veltot(ind:ind+size(pos,1)-1,l)=vel(:,l)*vsgn(l,i,type)
      end do
      ind=ind+size(pos,1)
	end do
    ! y'
    !t1 = -cos(theta)*cos(phi) * psgn(1,n) 
    !t2 = -cos(theta)*sin(phi) * psgn(2,n) 
    !t3 =  sin(theta)          * psgn(3,n)
    !proj(:,2) = t1 * pos(:,1) + t2 * pos(:,2) + t3 * pos(:,3)

    ! v_LOS
    !t1 = sin(theta)*cos(phi)  * vsgn(1,n,type) 
    !t2 = sin(theta)*sin(phi)  * vsgn(2,n,type) 
    !t3 = cos(theta)           * vsgn(3,n,type)
    !losvel(:) = t1 * vel(:,1) + t2 * vel(:,2) + t3 * vel(:,3)



  end subroutine project_part
 

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine stob()
    use integrator        , only : integrator_stop
    use projection        , only : projection_stop
    use histograms        , only : histogram_stop
    use psf               , only : psf_stop
    use aperture_routines , only : aperture_stop
    use output            , only : output_close
  !----------------------------------------------------------------------

    call output_close    ()
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
  print*,"  * $Id: partgen.f90,v 1.1 2010/03/17 12:56:57 bosch Exp $"  

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
 
