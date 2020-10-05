module intrinic_mass
  use numeric_kinds
  implicit none
  private
  
  private :: intrin_spher_grid_func, intrin_spher_grid

  public ::  intrin_spher,intrin_radii

  real (kind=dp),private, dimension(6) :: intrin_spher_grid_func_global 
  character (len=512), public :: orblib_filename
contains
  
  function intrin_spher_grid_func(dimen,coordinates)  result(answ)
    integer (kind=i4b), intent(in):: dimen
    real (kind=dp), intent (in), dimension(2) :: coordinates
    real (kind=dp) :: answ,c
    
    real (kind=dp) :: sth, cth  
    real (kind=dp) :: sph, cph
    real (kind=dp) :: P ,Q,sigma,r0,r1,rho0
    
    real (kind=dp) :: angular_part,res,term1,term2,rad1,rad2

    ! These come from the Global parameter
    P      = intrin_spher_grid_func_global(1)
    Q      = intrin_spher_grid_func_global(2)
    sigma  = intrin_spher_grid_func_global(3)
    r0     = intrin_spher_grid_func_global(4)
    r1     = intrin_spher_grid_func_global(5)
    rho0   = intrin_spher_grid_func_global(6)
    
    sth=sin(coordinates(1))
    cth=cos(coordinates(1))
    sph=sin(coordinates(2))
    cph=cos(coordinates(2))
    
    c=sqrt(2.)*sigma/sqrt(  (cph**2. + (sph/p)**2. )*sth*sth + (cth/q)**2.)
     
    term1 =  r0 * exp(-(r0/c)**2) - 0.5*c*sqrt(pi_d)*derf( r0/c)
    term2 =  r1 * exp(-(r1/c)**2) - 0.5*c*sqrt(pi_d)*derf( r1/c) 

    answ = 0.5 * rho0 *c * c * ( term1 - term2 ) * sth 

    ! the orbit library folds the 8 symmetries into one octant: 
    answ = answ *8.0_dp

  end function intrin_spher_grid_func
  
  
!!! calculate multi-dimensional integrals !!!
  subroutine intrin_spher_grid(low,up,res)
    integer (kind=i4b) ,parameter:: dim = 2
    real (kind=dp), intent(in), dimension(dim) :: low,up
    interface
       function integrand(dimen,coordinates)  result(answ)
         use numeric_kinds
         implicit none
         integer (kind=i4b), intent(in):: dimen
         real (kind=dp), intent (in), dimension(dimen) :: coordinates
         real (kind=dp) :: answ
       end function integrand
    end interface
    real (kind=dp), intent(out) :: res
    real (kind=dp), dimension(dim) :: lowL,upL
    logical :: diffcoord
    integer (kind=i4b) :: ii, intgsign
    real (kind=dp) :: dum
    ! definitions that are required for the integrator
    integer (kind=i4b) :: minpts,maxpts,ifail
    integer (kind=i4b) ::lenwrk
    real (kind=dp):: acc,eps
    real (kind=dp), dimension(:), allocatable:: wrkstr  
    !
    minpts=1
    maxpts=100000*dim
    lenwrk=(dim+2)*(1+maxpts/(2**dim+2*dim*dim+2*dim+1))
    ifail=1
    eps = 1.0e-8
    allocate(wrkstr(lenwrk))
    ! check integration limits
    lowL=low
    upL=up
    diffcoord=.true.
    intgsign=1
    do ii=1,dim
       if(lowL(ii) == upL(ii)) then
          print*,"lower and upper limit are equal"
          diffcoord=.false.
       end if
       if(lowL(ii)>upL(ii)) then
          print*,"changing order of integration"
          dum=lowL(ii)
          lowL(ii)=upL(ii)
          upL(ii)=dum
          intgsign=-intgsign
       end if
    end do
    if(diffcoord) then 
       call d01fcf(dim,lowL,upL,minpts,maxpts,intrin_spher_grid_func,eps,acc,&
            lenwrk,wrkstr,res,ifail)  
       res=intgsign*res
    else 
       res=0.0e0
    end if
    if (ifail /= 0 ) print*,ifail,res,acc
    if (ifail /= 0 .and. ifail /= 2) stop "integrator failed" 
  end subroutine intrin_spher_grid
  
  subroutine intrin_spher()
    use initial_parameters , only : totalmass,ngauss_mge,qobs
    use triaxpotent        , only : pintr, qintr,sigintr_km,dens
    real (kind=dp),dimension(:),allocatable :: quad_lr
    real (kind=dp),dimension(:),allocatable  :: quad_lth
    real (kind=dp),dimension(:),allocatable  :: quad_lph
    real (kind=dp),dimension(:,:,:),allocatable  :: quad_grid
    real (kind=dp),dimension(2) :: low,up
    integer (kind=i4b) :: i,j,k,l,quad_nr,quad_nth,quad_nph,quad_mom
    integer (kind=i4b) :: norbits,ni1,ni2,ni3,ndith
    character (len=512) :: mass_qgrid_filename
    real (kind=dp) :: res

    open (unit=27,file=orblib_filename,action="read", &
         status="old",form="unformatted",position="rewind")
    read (unit=27) norbits,ni1,ni2,ni3,ndith
    read (unit=27) quad_mom,quad_nph,quad_nth,quad_nr

    print*, quad_mom,quad_nr, quad_nth,quad_nph      
    ! remember that N bins have N+1 boundaries
    allocate (quad_grid(quad_nph,quad_nth,quad_nr))
    allocate (quad_lr (quad_nr +1))
    allocate (quad_lth(quad_nth+1))
    allocate (quad_lph(quad_nph+1))

    read (unit=27) quad_lr (:)
    read (unit=27) quad_lth(:)
    read (unit=27) quad_lph(:)
    close(unit=27)
    
    print*,quad_lth
    print*,quad_lph

    quad_grid(:,:,:)=0.0_dp

    do i=1,quad_nr
       do j=1,quad_nth
          do k=1,quad_nph
             do l=1,ngauss_mge
                intrin_spher_grid_func_global(1) = pintr(l)      ! P
                intrin_spher_grid_func_global(2) = qintr(l)      ! Q
                intrin_spher_grid_func_global(3) = sigintr_km(l) ! sigma
                intrin_spher_grid_func_global(4) = quad_lr(i)    ! r0
                intrin_spher_grid_func_global(5) = quad_lr(i+1)  ! r1
                intrin_spher_grid_func_global(6) = dens(l)       ! rho0
                low    = (/ quad_lth(j  ),quad_lph(k  ) /)
                up     = (/ quad_lth(j+1),quad_lph(k+1) /)

                call intrin_spher_grid(low,up,res)  

                quad_grid(k,j,i)=quad_grid(k,j,i)+res 
             end do
          end do
       end do
    end do

    print*,"procent of the Mass inside the projected grid"
    print*,sum(quad_grid(:,:,:))/totalmass*100.0_dp

    quad_grid(:,:,:) = quad_grid(:,:,:)/totalmass

    print*,"Give the filename of the output mass_qgrid file"
    read (unit=*, fmt="(a512)"), mass_qgrid_filename
    print*," Writing ", mass_qgrid_filename
    open (unit=28,file=mass_qgrid_filename,action="write", &
         status="new")
    write (unit=28,fmt=*) size(quad_lph)-1,size(quad_lth)-1,size(quad_lr)-1
    write (unit=28,fmt=*) quad_grid(:,:,:)
    close (unit=28)

  end subroutine intrin_spher

  subroutine intrin_radii()
    use initial_parameters , only : totalmass,ngauss_mge,qobs, & 
             rlogmin,rlogmax,nener,orbit_dithering 
    use triaxpotent        , only : pintr, qintr,sigintr_km, &
                              dens
    real (kind=dp),dimension(nener/orbit_dithering +1) :: quad_lr
    real (kind=dp),dimension(2)  :: quad_lth
    real (kind=dp),dimension(2)  :: quad_lph
    real (kind=dp),dimension(10000) :: radmass
    real (kind=dp),dimension(2) :: low,up
    integer (kind=i4b) :: i,j,k,l,quad_nr,quad_nth,quad_nph
    integer (kind=i4b) :: norbits,ni1,ni2,ni3,ndith,nr
    character (len=512) :: mass_radmass_filename
    real (kind=dp) :: res
    
    nr=nener/orbit_dithering 
    if (nr > size(radmass) ) stop 'increase arraysize radmass'
    print*,nr
    do i=1,nr
       quad_lr(i)=10.0_dp**(rLogMin + (rLogMax-rLogMin) * (i-2.0)/(nr-1.0))
    end do
    quad_lr(1)=0.0

    quad_lr(nr+1)=10.0_dp**rlogmax*100.0

    quad_lth(1)=0.00
    quad_lph(1)=0.00
    
    quad_lth(2)=pio2_d
    quad_lph(2)=pio2_d
    
    do i=1,nr
       do j=1,1
          do k=1,1
             do l=1,ngauss_mge
                intrin_spher_grid_func_global(1) = pintr(l)      ! P
                intrin_spher_grid_func_global(2) = qintr(l)      ! Q
                intrin_spher_grid_func_global(3) = sigintr_km(l) ! sigma
                intrin_spher_grid_func_global(4) = quad_lr(i)    ! r0
                intrin_spher_grid_func_global(5) = quad_lr(i+1)  ! r1
                intrin_spher_grid_func_global(6) = dens(l)       ! rho0
                low    = (/ quad_lth(j  ),quad_lph(k  ) /)
                up     = (/ quad_lth(j+1),quad_lph(k+1) /)

                call intrin_spher_grid(low,up,res)  
                print*,i,l,res/totalmass*100.0_dp
                radmass(i)=radmass(i)+res
             end do
          end do
       end do
    end do

    print*,"procent of the Mass inside the radial shells"
    print*,sum(radmass(:))/totalmass*100.0_dp
    radmass(:) = radmass(:)/totalmass
    print*,"radmass"
    print*,radmass(1:nener/orbit_dithering)

    print*,"Give the filename of the output mass file"            
    read (unit=*, fmt="(a512)"), mass_radmass_filename
    print*," Writing ", mass_radmass_filename
    open (unit=28,file=mass_radmass_filename,action="write", &    
         status="new")
    write (unit=28,fmt=*) nener/orbit_dithering 
    write (unit=28,fmt=*) radmass(1:nener/orbit_dithering )
    close (unit=28)

  end subroutine intrin_radii
  
end module intrinic_mass

program triaxmass
  use intrinic_mass
  use initial_parameters
  use triaxpotent
  use numeric_kinds
  implicit none

  call iniparam()

  print*,"Give the filename of the orbit library"
  read (unit=*, fmt="(a512)"), orblib_filename  

  call tp_setup()

  call intrin_radii()
  call intrin_spher()


 end program Triaxmass
