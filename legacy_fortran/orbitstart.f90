module orbitstart
  use numeric_kinds
  use interpolpot
  implicit none
  private

  ! accuracy of integrator
  real (kind=dp),private :: integrator_accuracy=1d-4

  ! temporary pos and velocity array for the dense output of the integrator
  real (kind=dp),private,allocatable,dimension(:,:) :: vel_t,pos_t

  ! functions and variables copied from original orblib
  private :: derivs,soloutob,findtubeorbitwidth,calc_startpos,findtube

  public :: runorbitstart

  private :: findReq,make_startpoints,find_type,solouttype
  character (len=256), private  :: orbitstart_status, orbitstart_begin, orbitstart_beginbox
contains

  subroutine runorbitstart()

    use initial_parameters, only: nEner,nI2,nI3,rlogmin,rlogmax,&
         psi_view,orbit_dithering,conversion_factor
    use triaxpotent, only: pintr, triaxpar

    real    (kind=dp ), dimension(nEner):: rCirc,vcirc,tcirc,Epot
    real    (kind=dp ), dimension(nI2)  :: Theta,rib,rob,focplane
    real    (kind=dp ), dimension(nEner,nI2) :: boundout,boundin,boundmid
    real    (kind=dp ), dimension(nEner,nI2) :: noreggrid

    integer (kind=i4b), dimension(nEner):: irregular
    integer (kind=i4b), dimension(nEner,Ni2,Ni3) :: grid
    real    (kind=dp ), dimension(6):: startpos
    integer (kind=i4b) :: i,j,k
    real    (kind=dp ) :: pot,vx,vy,vz
    real    (kind=dp ) ::tube
    real    (kind=dp ) ::t1,t2,t3
    integer            ::i1,i2

 	read (unit=*, fmt="(a256)") orbitstart_status
	read (unit=*, fmt="(a256)") orbitstart_begin
	read (unit=*, fmt="(a256)") orbitstart_beginbox  
	
	print*,"Maximum triaxility:",maxval(triaxpar)

	! define the orbits logarithmically.
	do i=1,nEner
   		rcirc(i)=rLogMin + (rLogMax-rLogMin) *&
        	(i-1.0)/(Nener-1.0)
	end do
 
    rcirc(:) = 10.0_dp ** rcirc(:)
    print*,'Radiia of the energies in arcsec:'
    print*,rcirc(:)/conversion_factor

    ! calculate the circular velocities at the radii        
    do i=1,nEner
       ! estimate the tube is at 0.5_dp           
	   ! [Later in the program this will be computed correctly]
       t1=0.5_dp

       call ip_accel(rcirc(i)*t1,0.0_dp,0.0_dp,vx,vy,vz)
       vcirc(i) = sqrt(       rcirc(i)*t1*abs(vx)  )
       tcirc(i) = 2.0_dp*pi_d*rcirc(i)*t1/vcirc(i)
       call ip_potent(rcirc(i),0.0_dp,0.0_dp,Epot(i))
    end do

    print*,'Potentials of the energies'
    print*,epot(:)

    ! time for a approx. circular orbit.
    tcirc(:)= 2.0_dp*pi_d*rcirc(:)*t1/vcirc(:)

    ! The angular space in x,z start space is open.
    ! This is done to avoid duplicate orbits at x=0
    ! and to avoid the (un-integratable) centrophilic box orbits at z=0 
    do i=1,nI2
       Theta(i) = Pi_d/2.0_dp * (0.5+ i-1 ) / nI2 !( nI2+1.0_dp )
    end do

    print*,"Finding R_eq"
    do j=1,NEner
       do i=1,nI2
          call findReq(Rcirc(j),Epot(j),theta(i),0.0_dp,boundout(j,i))
       end do
    end do

    call find_innerboundary(boundin,irregular,boundout,tcirc,epot,theta)

    do i=1,nEner
       ! Use the position of the tube to define the 'circular' velocity
       t1=boundin(i,Ni2)/rcirc(i)
       call ip_accel(rcirc(i)*t1,0.0_dp,0.0_dp,vx,vy,vz)
       vcirc(i) = sqrt( rcirc(i)*t1*abs(vx)  )
       tcirc(i)= 2.0_dp*pi_d*rcirc(i)*t1/vcirc(i)
    end do

    
    ! check if galaxy model is triaxial.
    if ( maxval (triaxpar(:)) > 0.01_dp .or. minval(triaxpar(:)) > 0.99_dp) then 
       ! galaxy is triaxial
       print*,"galaxy is triaxial"
	
       call find_outerboundary(boundin,boundmid,boundout,irregular,tcirc,epot,theta)

       ! mark all the energys inside the last "irregular energy" as irregular
       j=0
       do i=1,Nener
          if (irregular(i) == 1) j=i
       end do
       print*,"last irregular energy is :",j

       ! round irregular energy up to nearest dithering boundary
       ! (to make sure each orbital dithering is similar)
       j= floor((j-1.0_dp+orbit_dithering)/orbit_dithering)*orbit_dithering
       print*,"moved for dithering to:",j

       if (j > 0 ) irregular(1:j) = 1
    else
       ! special settings for the axisymmetric case.
       print*,"galaxy is assumed to be axisymmetric"
       irregular(:)=0
       boundmid=boundout
    endif

    
    call find_unregorbits(boundout,boundmid,irregular,noreggrid)

    call make_startpoints(boundin,boundmid,boundout,irregular,noreggrid, &
         theta,epot,tcirc,rcirc,vcirc)

    call make_boxstartpoints(epot,tcirc,rcirc,vcirc) 

    ! write a status file, confirming that we finished succesfully
    open (unit=31,file=orbitstart_status,status="new",action="write",&
           position="rewind")
    write (unit=31, fmt=*) " finished "
    close (unit=31)

end subroutine  runorbitstart

subroutine  find_unregorbits(boundout,boundmid,irregular,noreggrid)
 use initial_parameters, only: nEner,nI2,nI3
  real    (kind=dp ),intent(in), dimension(nEner,nI2) :: boundmid,boundout
  real    (kind=dp ),intent(out), dimension(nEner,nI2) :: noreggrid
  integer (kind=i4b),dimension(Nener),intent(in) :: irregular
  integer (kind=i4b) :: i,j,k,noreg

! don't allow long axis tubes at the border to be regularized with 
! the box orbits 

do i=1,nEner
   noreg=0
   do j=nI2,1,-1
      if ( abs(boundmid(i,j)-boundout(i,j))/boundout(i,j) .gt. 1.0e-5  &
           .and. irregular(i) .eq. 0 )  &
           noreg=1
      noreggrid(i,j) = noreg
   end do
end do
end subroutine find_unregorbits
      
subroutine make_startpoints(boundin,boundmid,boundout,irregular &
                           ,noreggrid,theta,epot,tcirc,rcirc,vcirc)
 use initial_parameters, only: nEner,nI2,nI3
 real    (kind=dp ),intent(in),dimension(nEner):: tcirc,Epot,rcirc,vcirc
 real    (kind=dp ),intent(in), dimension(nI2)  :: Theta
 real    (kind=dp ),intent(in), dimension(nEner,nI2) :: boundin,boundmid
 real    (kind=dp ),intent(in), dimension(nEner,nI2) :: boundout,noreggrid
 integer (kind=i4b),dimension(Nener),intent(in) :: irregular
 real    (kind=dp ),dimension(nEner,nI2) :: bndin,bndmid
 real    (kind=dp ) :: t1
 real    (kind=dp ), dimension(6):: startpos
 integer (kind=i4b) :: i,j,k
 integer (kind=i4b) :: noreg,LastIrregE

 print*," Writing ", orbitstart_begin

 bndin(:,:)=boundin(:,:)
 bndmid(:,:)=boundmid(:,:)
 LastIrregE=0
 do i=1,NEner
	! If this energy is classified as irregular sample the whole X-Z plane.
    if (irregular(i) .eq. 1) then 
	  bndin (i,:)=0.0_dp
      bndmid(i,:)=boundout(i,:)
	  LastIrregE =i ! Keep track of the last irregular energy
    endif
 end do
 
  open (unit=31,file=orbitstart_begin,status="new",action="write" , &
      position="rewind")
 write (unit=31, fmt=*) nEner,nI2,nI3

 do i=1,NEner
    ! reset the do_not_regularize counter for each energy
    do j=1,nI2
       do k=1,nI3
          ! open boundaries:
          !t1= bndin(i,j) + (bndmid(i,j)-bndin(i,j)) * (k-0.5_dp)/nI3
		  ! nearly closed boundaries  (fully close boundaries do not work. RvdB 27032010):
          t1= bndin(i,j) + (bndmid(i,j)-bndin(i,j)) * (k-0.9_dp)/(nI3-0.8_dp)
          call calc_startpos(t1,theta(j),Epot(i), startpos)
          noreg=0
          
          ! do not regularize the orbits on the boundmid
          if (k .eq. nI3 .and. noreggrid(i,j) .eq. 1 ) noreg=1 
          ! if this is the last irregular energy do not regularize
          if (maxval(irregular(:)) .eq. i)  noreg=1

          write (unit=31, fmt="(3I5,9ES30.10,I4)") i,j,k,startpos,rcirc(i),tcirc(i),vcirc(i),noreg
       enddo
    enddo
 enddo
 close (unit=31)
end subroutine make_startpoints

subroutine make_boxstartpoints(epot,tcirc,rcirc,vcirc)
 use initial_parameters, only: nEner,nI2,nI3
 real    (kind=dp ),intent(in),dimension(nEner):: tcirc,Epot,rcirc,vcirc
 real    (kind=dp ) :: t1,theta,phi,rorbit,xp,yp,zp
 integer (kind=i4b) :: i,j,k

 print*," Writing ", orbitstart_beginbox
  open (unit=31,file=orbitstart_beginbox,status="new",action="write",&
      position="rewind")
 write (unit=31, fmt=*) nEner,nI2,nI3
 do i=1,NEner
    do j=1,nI2
       do k=1,nI3
          ! open boundaries:
          Theta = Pi_d/2.0_dp * ( j -0.5_dp) / nI2
          Phi   = Pi_d/2.0_dp * ( k -0.5_dp) / nI3
    
         call findReq(Rcirc(i),Epot(i),theta,phi,rorbit)
    
         xp=Rorbit*sin(theta)*cos(phi)   ! x
         yp=Rorbit*sin(theta)*sin(phi)   ! y 
         zp=Rorbit*cos(theta)            ! z
         write (unit=31, fmt="(3I5,9ES30.10,I4)") i,j,k,xp,yp,zp,0.0_dp & 
                    ,0.0_dp,0.0_dp,rcirc(i),tcirc(i),vcirc(i),0
       enddo
    enddo
 enddo
 close (unit=31)
end subroutine make_boxstartpoints

subroutine find_innerboundary(boundin,irregular,boundout,tcirc,epot,theta)
  use initial_parameters, only: nEner,nI2,nI3,rlogmin,rlogmax
  real    (kind=dp ),intent(in),dimension(nEner):: tcirc,Epot
  real    (kind=dp ),intent(in), dimension(nI2)  :: Theta
  real    (kind=dp ),intent(in), dimension(nEner,nI2) :: boundout
  real    (kind=dp ),intent(out), dimension(nEner,nI2) :: boundin
  integer (kind=i4b),dimension(Nener),intent(out) :: irregular
  real    (kind=dp ), dimension(6):: startpos
  integer (kind=i4b) :: i,j
  integer (kind=i4b),dimension(nI2) :: otype
  real    (kind=dp ) :: pot,vx,vy,vz,rbu,rbi,rg,acc
  integer (kind=i4b) :: orbtype,foundxtubes
  real    (kind=sp ) ::t1,t2,t3
  integer (kind=i2b) ::i1,i2
  
  print*,"Finding inner tubes"
  
  if (nI2 <= 3) stop "nI2 is smaller then 4"
	
  irregular(:)=0
  do j=1,nEner
	
	 ! Setup bounds for first guess.
     foundxtubes=0    
     i=nI2
     rbi=boundout(j,i)*0.11
     rbu=boundout(j,i)*0.89
     rg =boundout(j,i)*0.50
     if (j>1) rg=boundin(j-1,i)/boundout(j-1,i)*boundout(j,i)
     rg= max(min(rg,boundout(j,i)*0.88),boundout(j,i)*0.12)

     ! Find first 
     call findtube(rbi, rg,rbu,Epot(j),tcirc(j),theta(i),2,boundin(j,i))
     call find_type(boundin(j,i),theta(i),Epot(j),tcirc(j),otype(i))
     t1=boundin(j,i)/boundout(j,i)
     print*,j,i,t1,nEner,otype(i)

     ! find second
     i=nI2-1
     rg=boundin(j,i+1)
     rg= max(min(boundin(j,i+1) ,boundout(j,i)*0.99),boundout(j,i)*0.02)
     rbi=max(rg - boundout(j,i)*0.10,boundout(j,i)*0.01)
     rbu=min(rg + boundout(j,i)*0.10,boundout(j,i)     )
     call findtube(rbi,rg,rbu,Epot(j),tcirc(j),theta(i),2,boundin(j,i))
     call find_type(boundin(j,i),theta(i),Epot(j),tcirc(j),otype(i))
     t1=boundin(j,i)/boundout(j,Ni2)
     print*,j,i,t1,ni2,otype(i)
     
     do i=nI2-2,1,-1
        ! Use previous radius as guess.
        rg=boundin(j,i+1)
        
        ! set sane boundaries
        rg =min(max(rg ,boundout(j,i)*0.11),boundout(j,i)*0.99)
        rbi=max(rg - boundout(j,nI2)*0.08,boundout(j,i)*0.10)
        rbu=min(rg + boundout(j,nI2)*0.08,boundout(j,i)*1.00)

        call findtube(rbi,rg,rbu,Epot(j),tcirc(j),theta(i),2,boundin(j,i))
        if (abs(boundin(j,i)-rbi) < boundin(j,i)*1.0e-2  .or. &
            abs(boundin(j,i)-rbu) < boundin(j,i)*1.0e-2 ) then
           print*,'tube moved to edge of allowed boundary. irregular set'
           irregular(j)=1
           !exit
        endif
        call find_type(boundin(j,i),theta(i),Epot(j),tcirc(j),otype(i))
        ! If we detect boxorbits as a type orbit, define irregular energy
        if (otype(i) == 1) then
            foundxtubes=1
        endif
        t1=boundin(j,i)/boundout(j,nI2)
        t2=rbu/boundout(j,i)
        print*,j,i,t1,t2,otype(i)
      
     enddo
     
     if (foundxtubes == 0) irregular(j)=1
     if (irregular(j) == 1) print*,"irregular energy"
  end do
end subroutine find_innerboundary


subroutine find_outerboundary(boundin,boundmid,boundout,irregular,&
                              tcirc,epot,theta)
  use initial_parameters, only: nEner,nI2,nI3,rlogmin,rlogmax
  real    (kind=dp ),intent(in),  dimension(nEner):: tcirc,Epot
  real    (kind=dp ),intent(in),  dimension(nI2)  :: Theta
  real    (kind=dp ),intent(in),  dimension(nEner,nI2) :: boundout
  real    (kind=dp ),intent(out), dimension(nEner,nI2) :: boundmid
  real    (kind=dp ),intent(in ), dimension(nEner,nI2) :: boundin
  integer (kind=i4b),dimension(Nener),intent(inout)    :: irregular
  real    (kind=dp ), dimension(6):: startpos
  integer (kind=i4b) :: i,j,k,l,o1,o2,o3,notubes
  integer (kind=i4b),dimension(nI2) :: otype
  real    (kind=dp ) :: pot,vx,vy,vz,rbu,rbi,rg,acc,bp,r,rel_rbi
  integer (kind=i4b) :: nodice
  real    (kind=dp ) ::t1,t2,t3
  integer (kind=i2b) ::i1,i2

  print*,"Finding outer tubes"
  if (nI2 <= 3) stop "nI2 is smaller then 4"

  ! if we cannot find the outer thin tubes use equipotential.
  boundmid(:,:)=boundout(:,:)

  do j=1,NEner
     notubes=0
     ! find box orbits at last theta first and work backwards
     i=nI2
                    
     ! rel_rbi will be the relative radius at which tubes convert box orbits. 
     ! But first we underestimate it by the rel. radius of the thin short orbit tubes
     rel_rbi=(boundin(j,i))/boundout(j,i)

     ! Find the radius at which (short) axis tubes convert to box orbits
     do k=1,nI3*3
	            
	      ! we keep increasing r until we hit boxes
	      r=boundin(j,i) + (boundout(j,i)-(boundin(j,i))) *k / (nI3*3+1)
	      call find_type(r,theta(i),Epot(j),tcirc(j),o1)          
	      
	      if (o1 .eq. 3) rel_rbi= r/boundout(j,i) ! Have not found boxes (yet), so we set
	                                              ! so we move up rel_rbi

          ! if we find boxes we are done.	 
		  if ((o1 .eq. 1 .or. o1 .eq. 4) .and. k .ge. 2 ) exit
     end do	                        

     ! if we couldn't find any box orbits then there is no point in looking for thin tubes,
     ! as further in would not add more orbit resolution 
      if ((o1 .eq. 3 .and. k .ge. nI3))  then  
       print*,'Short to box conversion at',rel_rbi,', which is to close enough to the outer boundary'
       ! to not search for intermediate boundary'
       
	   irregular(j)=0
	   notubes=1

	 else
      
      if (o1 .eq. 4 .or. o1 .eq. 5 ) then 
     	! found boxes. Now find the transitition to x tubes in the theta direction.    
 
        ! r is the radius at which boxes are found in the previous step.
        ! bp is the relative radius at which the boxes are found. 
     	bp=r/boundout(j,i)      

     	do i=nI2-1,1,-1          
	
	        ! See if the're tubes at this rel. radius.
        	call find_type(boundout(j,i)*bp,theta(i),Epot(j),tcirc(j),o1)
        	if ( o1 == 1 ) exit ! x tubes found; done
	    	print "(3I5 , 4F7.3)", -j,i,o1,boundin(j,i)/boundout(j,Ni2),rel_rbi*boundout(j,i)/boundout(j,Ni2), &
	               boundmid(j,i)/boundout(j,Ni2),boundout(j,i)/boundout(j,Ni2)
     	enddo
        i=minval( (/ i+1 , nI2-1/) )
      endif  
      endif
    
     if ( notubes == 0 .and. i > 1) then
	    ! Found boxes and long axis tubes. Now find the thin x-axis tubes.    
   
        ! Find the thin long axis tubes is not very robust. I sometimes wonder 
        ! why we bother. The current approach tries to be robust by creating the
        ! smallest possible search boundaries, based on the previous estimate and 
        ! the expectation that the radius does not decrease faster than 3 times
        ! the maxiumum possible avg. linear speed. RvdB 27032010  
	
	    boundmid(j,i)=boundout(j,i)
	
        if ( i > 2) then 

           do k=i-1,1,-1        
	
	          ! inner boundary is either the relative tube-to-box conversion radius or
	          ! the thin tube orbit radius.
              rbi=maxval( (/boundout(j,k)*rel_rbi,boundin(j,:) /))

              ! Outer boundary is the equipotential
              rbu=boundout(j,k)	- 1.0e-6 ! Make sure we do not go over the equipotential
                                         ! as then the V_y velocity becomes undefined
                                         ! Rvdb 27032010


              ! boundmid should always shrink, due to galaxy flattening, hence we can assume 
              ! that boundmid(j,k) < boundmid(j,k+1)
              rbu=minval((/boundmid(j,k+1),rbu/))

              ! guess that the new boundmid is the relative boundmid of the previous step.
              rg =maxval((/minval((/boundmid(j,k+1)/boundout(j,k+1)*boundout(j,k),rbu/)),rbi/))

              ! In the most extreme case boundmid could touch rel_rbi at the end,
              ! this we allow boundmid to go towards that in linear steps.
              ! (This is not perfect, but works well.)  
              r=((1-rel_rbi)/Ni3*3)                  
              rbi=maxval((/rg*(1-r),rbi/)) 
                         

              !print "(3F7.3)" ,rbi/boundout(j,k),rg/boundout(j,k),rbu/boundout(j,k)
              ! Sanity check on the boundaries
              if (rbi .ge. rbu) stop 'inner boundary bigger than outer boundary'
              if (rg .le. rbi .or. rg .ge. rbu) rg = 0.5*(rbi+rbu) 

              ! find the thin long axis tubes radius 
			  call findtube(rbi,rg,rbu,Epot(j),tcirc(j),theta(k),1,r)
              boundmid(j,k)=r

              call find_type(r,theta(k),Epot(j),tcirc(j),o1)

			  print "(3I5 , 4F7.3)" ,j,k,o1,boundin(j,k)/boundout(j,Ni2),rel_rbi*boundout(j,k)/boundout(j,Ni2), &
			  boundmid(j,k)/boundout(j,Ni2),boundout(j,k)/boundout(j,Ni2) 
			
           end do
        endif
     endif

  end do

end subroutine find_outerboundary

subroutine calc_startpos(r,theta,eeq,Y)
  real(kind=dp),intent(in) :: r,theta,eeq
  real(kind=dp),intent(out), dimension(6) :: Y
  real(kind=dp) :: E
  Y(1)=R*sin(theta)            ! x
  Y(2)=0.0_dp                  ! y 
  Y(3)=R*cos(theta)            ! z
  call ip_potent(Y(1),Y(2),Y(3),E) 
  Y(4)=0.0_dp                  ! vx
  Y(5)=2.0_dp *(E- Eeq)        ! vy^2
  if ( Y(5) >= 1.0e-300_dp) Y(5)=sqrt(Y(5)) 
  if ( Y(5) <  0.0_dp .or. ISNAN(Y(5))) then               
	print*,r,theta,eeq,e,2*(e-eeq),Y(1:6), 'Starting velocity is negative or NAN' ,sqrt(2.0_dp *(E*1.0e-8))                
	Y(5)= sqrt(2.0_dp *(E*1.0e-12))! set vy to a really small positive value
  endif
  Y(6)=0.0_dp                  ! vz
end subroutine calc_startpos

  subroutine findReq(Req,E,th,ph,R)  
    real    (kind=dp ),intent(in ):: Req,E,th,ph
    real    (kind=dp ),intent(out):: R
    real    (kind=dp )            :: pot,emn,emx,rmn,rmx
    real    (kind=dp )            :: sth,cth,sph,cph
    integer (kind=i4b)            :: i

    ! Find the position R for which pot(E,th) = E
    
    ! th is defined from the other side:
    sth = sin(th)
    cth = cos(th)

    sph = sin(ph)
    cph = cos(ph)
    
    ! bisection
    rmx = 1.1_dp * Req
    rmn = 0.01_dp * Req
    
    call ip_potent( Rmx*sth*cph,Rmx*sth*sph,Rmx*cth,Emx) 
    call ip_potent( Rmn*sth*cph,Rmn*sth*sph,Rmn*cth,Emn)

    i=0

    do
       R = 0.5_dp*(Rmn+Rmx)
      
       call ip_potent(R*sth*cph,R*sth*sph,R*cth,pot)
        if ( abs(( E - pot )/E) < 1.0e-7_dp ) exit
       
       if (pot > E ) then 
          Rmn = 0.5_dp*(Rmn+Rmx)
          Emn = pot
       else
          Rmx = 0.5_dp*(Rmn+Rmx)
          Emx = pot
       end  if

       i=i+1

       if ( i > 60000_i4b ) then 
	   print*,E,emn,emx 
	   print*,R,rmn,rmx ,1.5_dp * Req ,0.01_dp * Req 
	   stop "Can not find R_eq"
       endif

    end do

  end subroutine findReq

 subroutine  findtube(Rin,Rmid,Rout,Epot,tcirc,theta,plane,R)  
    real    (kind=dp ),intent(in ):: Rin,Rmid,Rout,Epot,theta,tcirc
    integer (kind=i4b),intent(in)  :: plane
    real    (kind=dp ),intent(out):: R
    real    (kind=dp ),parameter  :: gldn=0.61803399_dp
    real    (kind=dp )            :: Tube,r0,r1,r2,r3,t1,t2
    real    (kind=dp )            :: sth,cth
    real    (kind=sp )            :: t4,t5,t6,t7
    integer (kind=i4b)            :: i

    ! golden Ratio search
    R0 = Rin
    R3 = Rout 

    if (abs(Rout-Rmid) > abs (Rmid-Rin)) then
       R1=Rmid
       R2=Rmid+(1.0_dp-gldn)*(Rout-Rmid)
    else
       R2=Rmid
       R1=Rmid-(1.0_dp-gldn)*(Rmid-Rin)
    endif
    call findtubeorbitwidth(R1,theta,Epot,tcirc,plane,T1)
    call findtubeorbitwidth(R2,theta,Epot,tcirc,plane,T2)
    do   
       if (abs(R3-R0) < 1.0e-4 *(abs(R1)+abs(R2))) exit 
       if (T2 < T1 ) then 
          R0=R1
          R1=R2
          R2=gldn*R1+(1.0-gldn)*R3
          T1=T2
          call findtubeorbitwidth(R2,theta,Epot,tcirc,plane,T2)
       else
          R3=R2
          R2=R1
          R1=gldn*R2+(1.0-gldn)*R0
          T2=T1
          call findtubeorbitwidth(R1,theta,Epot,tcirc,plane,T1)
       end if
     end do
    if (T1 < T2) then
       R=R1
    else
       R=R2
    endif
  end subroutine findtube


subroutine findtubeorbitwidth(R,theta,Eeq,tcirc,plane,Tube)
    real    (kind=dp ),intent(in ) :: R,theta,Eeq,tcirc
    real    (kind=dp), intent(out) :: tube
    integer (kind=i4b),intent(in)  :: plane
    !----------------------------------------------------------------------
    integer (kind=i4b),parameter :: intsteps=400 
	integer (kind=i4b),parameter :: N=6,nrdens=3,LWORK=11*N+8*NRDENS+21
    integer (kind=i4b),parameter :: LIWORK = NRDENS+21
    integer (kind=i4b) :: IOUT,IDID,ITOL
    real    (kind=dp ) :: X,Xend,TOL,E,tvx,tvy,tvz
    real    (kind=dp ), dimension(N) :: Y,RTOL,ATOL
    real    (kind=dp ), dimension(LWORK):: WORK
    integer (kind=i4b), dimension(LIWORK) :: IWORK
    real    (kind=dp ), dimension(2) :: RPAR
    integer (kind=i4b), dimension(1) :: IPAR
    real    (kind=dp ), dimension(intsteps) :: Rad
    ! Setting up the array for storing the points.
    allocate (pos_t(intsteps,3))
    ! --- DIMENSION OF THE SYSTEM
    IDID=0
    ! --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
    IOUT=2
    ! --- INITIAL VALUES
    X=0.0
    call calc_startpos(R,theta,eeq, Y)

    ! --- ENDPOINT OF INTEGRATION
    ! integration time 800 * tcirc
    Xend = 500.0_dp * intsteps * tcirc
    ! Stepsize. Make sure enough steps are in the integration by adding
    ! room for one extra step. RvdB, DK 16/06/03
    RPAR(2) = 0.0 ! unused 
    RPAR(1) = R   ! used for the relative bisection accuracy
    IPAR(1) = plane   ! find crossings of coordinate ( x = Y(1),etc)
    ! --- REQUIRED TOLERANCE
    ! Use same relative tolerance TOL for (R,z,VR,Vz)
    ! but ensures a fixed absolute tolerance of
    ! TOL*rcirc on (R,z) and TOL*vcirc on (VR,Vz).
    ! The integrator keeps the local error on Y(I)
    ! below RTOL(I)*ABS(Y(I)) + ATOL(I).
    TOL     = integrator_accuracy
    ITOL    = 1 ! Tolerances are vectors
    RTOL(:) = TOL
    ATOL(:) = 1e-8_dp !TOL* (/R,R,R,sqrt(R*abs( Y(5))),sqrt(R*abs( Y(5))),&
         !sqrt(R*abs( Y(5))) /)
    
    ! --- VALUES FOR PARAMETERS
    ! Default values are used when IWORK or WORK are zero
    IWORK(:)=0
     WORK(:)=0.0_dp
    ! Maximum of allowed steps (just really high) (100000)
    IWORK(1)=100000000_i4b
    !number of dense components needed. (only the positions in this cases)
    IWORK(5)=nrdens
    IWORK(21)=1
    IWORK(22)=2
    IWORK(23)=3
    !stiffnes detection (negative --> do not try to detect)
    IWORK(4)=-1
    ! maximal stepsize
    !WORK(6) = XEND/(4000+1.0_dp)
    ! Print no messages:
    IWORK(3)=-1
    ! --- CALL OF THE SUBROUTINE DOPRI8 ( The dop853 integrator.)
    CALL DOP853(N,derivs,X,Y,XEND, RTOL,ATOL,ITOL, SOLOUTOB,IOUT, &
         &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    select case (IDID)
    case (-1)
       stop "integrator:  INPUT IS NOT CONSISTENT,"
    case (-2)
       stop "integrator:  LARGER NMAX IS NEEDED,"
    case (-3)
       stop "integrator:  STEP SIZE BECOMES TOO SMALL."
    case (-4)
       stop "integrator:  PROBLEM IS PROBABLY STIFF (INTERRUPTED)."
    case default
       ! Integrating went ok!
    end select

    if (plane .eq. 1) then
       rad(:)=pos_t(:,2)**2.0_dp + pos_t(:,3)**2.0_dp
    endif

    if (plane .eq. 2) then
       rad(:)=pos_t(:,1)**2.0_dp + pos_t(:,3)**2.0_dp
    endif

    if (plane .eq. 3) then
       rad(:)=pos_t(:,1)**2.0_dp + pos_t(:,2)**2.0_dp
    endif

    tube  = maxval(sqrt(rad(:))) - minval(sqrt(rad(:)))
	
    deallocate(pos_t)

  end subroutine findtubeorbitwidth


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine derivs (N,xin,yin,dydx,RPAR,IPAR)
    use interpolpot, only : ip_accel
    integer           ,intent(in   )              :: N
    real    (kind=dp ),intent(in   )              :: xin
    real    (kind=dp ),intent(in   ),dimension(6) :: yin
    real    (kind=dp ),intent(out  ),dimension(6) :: dydx
    real    (kind=dp ),intent(inout),dimension(2) :: RPAR
    integer (kind=i4b),intent(inout),dimension(1) :: IPAR
  !----------------------------------------------------------------------

    ! subroutine which returns the right-hand side derivatives.
    !   x      = t
    !   yin(1) = x
    !   yin(2) = y
    !   yin(3) = z
    !   yin(4) = dx/dt
    !   yin(5) = dy/dt
    !   yin(6) = dz/dt

    ! First calculate the true accelerations at the given position
    dydx(1:3) = yin(4:6)
    call ip_accel(yin(1),yin(2),yin(3),dydx(4),dydx(5),dydx(6))

  end subroutine derivs

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE SOLOUTOB (NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
    ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
    ! --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
    integer          ,intent(in)    :: NR,N,ND
    real   (kind=dp ),intent(in)    :: XOLD,X
    integer          ,intent(inout) :: IRTRN
    real   (kind=dp ),intent(in)   ,dimension(N)    :: Y
    real   (kind=dp ),intent(in)   ,dimension(8*ND) :: CON
    integer(kind=i4b),intent(in)   ,dimension(ND)   :: ICOMP
    real   (kind=dp ),intent(inout),dimension(2)    :: RPAR
    integer(kind=i4b),intent(inout),dimension(1)    :: IPAR
    !--------------------------------------------------------------------
    real   (kind=dp ),save          :: XOUT,yold
    integer(kind=i4b),save          :: count=0
    real   (kind=dp )               :: contd8,R
    real   (kind=dp) :: Xmx,Xmn,Ymn,Ymx,Yt,Xt
    integer (kind=i4b) :: biscount,tp
    tp=ipar(1) ! find crossing of this plane. x=1, y=2 and x=3
    R=RPAR(1)
    IF (NR == 1) THEN
       count=0
    ELSE
       if (count >= size(pos_t,1)) IRTRN=-1
       IF (y(tp) * yold < 0.0_dp .and. count < size(pos_t,1)) THEN
       
          if (y(tp)> 0.0_dp) then
             Xmx = X
             Xmn = Xold
          else          
             Xmx = Xold
             Xmn = X
          endif
          Ymx=CONTD8(tp,Xmx,CON,ICOMP,ND)
          Ymn=CONTD8(tp,Xmn,CON,ICOMP,ND)
          biscount=0
          do
             Xt=(Xmx+Xmn)/2.0_dp
             Yt=CONTD8(tp,Xt,CON,ICOMP,ND)
             biscount=biscount+1
             if ( abs( Yt - 0.0_dp ) < R*1.0e-4 .or. biscount > 40) exit

             if (Yt < 0.0_dp ) then 
                Xmn = Xt
                Ymn = Yt
             else
                Xmx = Xt
                Ymx = Yt
             end  if
          end do
          if (biscount < 40) then
             ! only store value if the bisection converged.
             count = count+1
             pos_t(count,1) = CONTD8(1,Xt,CON,ICOMP,ND)
             pos_t(count,2) = CONTD8(2,Xt,CON,ICOMP,ND)
             pos_t(count,3) = CONTD8(3,Xt,CON,ICOMP,ND)
          end if
       END IF
    END IF
    yold=y(tp)
  END SUBROUTINE SOLOUTOB

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine find_type(R,theta,Eeq,tcirc,orbtype)
    real    (kind=dp ),intent(in ) :: R,theta,Eeq,tcirc
    integer (kind=i4b), intent(out) :: orbtype
  !----------------------------------------------------------------------

    integer (kind=i4b),parameter :: intsteps=5000 
	integer (kind=i4b),parameter :: N=6,NRDENS=N,LWORK=11*N+8*NRDENS+21
    integer (kind=i4b),parameter :: LIWORK = NRDENS+21
    integer (kind=i4b) :: IOUT,IDID,ITOL
    real    (kind=dp ) :: X,Xend,TOL,E,tvx,tvy,tvz
    real    (kind=dp ), dimension(N) :: Y,RTOL,ATOL
    real    (kind=dp ), dimension(LWORK)  :: WORK
    integer (kind=i4b), dimension(LIWORK) :: IWORK
    real    (kind=dp ), dimension(2) :: RPAR
    integer (kind=i4b), dimension(1) :: IPAR
    real (kind=dp) :: lxmax,lxmin,lymax,lymin,lzmax,lzmin,lxc,lyc,lzc
    real (kind=dp),dimension(intsteps) :: t
    real (kind=dp),parameter :: nul=0.0_dp

    ! setting up the array for storing the points.
    allocate (pos_t(intsteps,3),vel_t(intsteps,3))
    ! --- DIMENSION OF THE SYSTEM
    IDID=0
    ! --- OUTPUT ROUTINE (AND DENSE OUTPUT) IS USED DURING INTEGRATION
    IOUT=2
    ! --- INITIAL VALUES
    X=0.0
    call calc_startpos(R,theta,eeq, Y)

    ! --- ENDPOINT OF INTEGRATION
    XEND =  100.0_dp * tcirc
    ! Stepsize. Make sure enough steps are in the integration by adding
    ! room for a couple of  extra steps. RvdB 19/12/04
    RPAR(2) = XEND/(intsteps+4)
    RPAR(1) = 0.0 ! Unused
    ! --- REQUIRED TOLERANCE
    ! Use same relative tolerance TOL for (R,z,VR,Vz)
    ! but ensures a fixed absolute tolerance of
    ! TOL*rcirc on (R,z) and TOL*vcirc on (VR,Vz).
    ! The integrator keeps the local error on Y(I)
    ! below RTOL(I)*ABS(Y(I)) + ATOL(I).
    TOL     = integrator_accuracy
    ITOL    = 1 ! Tolerances are vectors
    RTOL(:) = TOL
    ATOL(:) = 1e-8_dp !TOL * (/R,R,R,  &
         !sqrt(R*abs(y(5))),sqrt(R*abs(y(5))),sqrt(R*abs(y(5))) /)
      
    ! Work arrays's
    ! --- VALUES FOR PARAMETERS
    ! Default values are used when IWORK or WORK are zero
    IWORK(:)=0
     WORK(:)=0.0_dp
    ! Maximum of allowed steps (just really high) (100000)
    IWORK(1)=100000000_i4b
    !number of dense components needed. (all in our case)
    IWORK(5)=NRDENS
    !stiffnes detection (negative --> do not try to detect)
    IWORK(4)=-1

    ! --- CALL OF THE SUBROUTINE DOPRI8 ( The dop853 integrator.)
    CALL DOP853(N,derivs,X,Y,XEND, RTOL,ATOL,ITOL, SOLOUTTYPE,IOUT, &
         &                  WORK,LWORK,IWORK,LIWORK,RPAR,IPAR,IDID)
    select case (IDID)
    case (-1)
       stop "integrator:  INPUT IS NOT CONSISTENT,"
    case (-2)
       stop "integrator:  LARGER NMAX IS NEEDED,"
    case (-3)
       stop "integrator:  STEP SIZE BECOMES TOO SMALL."
    case (-4)
       stop "integrator:  PROBLEM IS PROBABLY STIFF (INTERRUPTED)."
    case default
       ! Integrating went ok!
    end select

    !  Lx = y*Vz-z*Vy
    t= pos_t(:,2) * vel_t(:,3) - pos_t (:,3) * vel_t(:,2)
    lxc=maxval(t) * minval(t)

    !  Ly = z*Vx-x*Vz
    t= pos_t(:,3) * vel_t(:,1) - pos_t(:,1) * vel_t(:,3)
    lyc=maxval(t) * minval(t)

    !  Lz = x*Vy-y*Vx
    t= pos_t(:,1) * vel_t(:,2) - pos_t (:,2) * vel_t(:,1)
    lzc=maxval(t) * minval(t)

    ! assume orbit is chaotic, unless proven otherwise.
    orbtype=5
    
    ! X tube
    if (lxc > nul .and. lyc < nul .and. lzc < nul ) orbtype=1
    ! Y tube
    if (lxc < nul .and. lyc > nul .and. lzc < nul ) orbtype=2
    ! Z tube
    if (lxc < nul .and. lyc < nul .and. lzc > nul ) orbtype=3
    ! Box
    if (lxc < nul .and. lyc < nul .and. lzc < nul ) orbtype=4

    deallocate(pos_t,vel_t)

  end subroutine find_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE SOLOUTtype (NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
    ! --- PRINTS SOLUTION AT EQUIDISTANT OUTPUT-POINTS
    ! --- BY USING "CONTD8", THE CONTINUOUS COLLOCATION SOLUTION
    integer          ,intent(in)    :: NR,N,ND
    real   (kind=dp ),intent(in)    :: XOLD,X
    integer          ,intent(inout) :: IRTRN
    real   (kind=dp ),intent(in)   ,dimension(N)    :: Y
    real   (kind=dp ),intent(in)   ,dimension(8*ND) :: CON
    integer(kind=i4b),intent(in)   ,dimension(ND)   :: ICOMP
    real   (kind=dp ),intent(inout),dimension(2)    :: RPAR
    integer(kind=i4b),intent(inout),dimension(1)    :: IPAR
  !----------------------------------------------------------------------
    real   (kind=dp ),save          :: XOUT
    integer(kind=i4b),save          :: count=0
    real   (kind=dp )               :: contd8,step,rnd
    step=RPAR(2)
    IF (NR == 1) THEN
       XOUT=X+step
       count=0
    ELSE
       do
          ! Make sure count > integrator_points to make sure we do not
          ! go out of bounds on pos_t(:,:). RvdB, DK 16/06/03
          IF (X < XOUT .or. count >= size(pos_t,1)) exit
          count = count+1
          pos_t(count,1) = CONTD8(1,XOUT,CON,ICOMP,ND)
          pos_t(count,2) = CONTD8(2,XOUT,CON,ICOMP,ND)
          pos_t(count,3) = CONTD8(3,XOUT,CON,ICOMP,ND)
          vel_t(count,1) = CONTD8(4,XOUT,CON,ICOMP,ND)
          vel_t(count,2) = CONTD8(5,XOUT,CON,ICOMP,ND)
          vel_t(count,3) = CONTD8(6,XOUT,CON,ICOMP,ND)
          XOUT = XOUT + step
       end do
    END IF
  END SUBROUTINE SOLOUTtype




end module orbitstart

program orbitstart_run
use orbitstart
use initial_parameters
use interpolpot
use triaxpotent
call iniparam()
call ip_setup()
call runorbitstart()
call ip_stop()
end program orbitstart_run
