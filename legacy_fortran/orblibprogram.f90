!######################################################################
!
! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands
!
! HISTORY:
!
! V1.4: MC Verification of intrinsic velocity moments calculation.
!     Fixed significant bug in mergrid_store. Slightly revised 
!     implementation
!     to make it easier to check. Independent test of the results 
!      still needed
! V2.0: RvdB. Fork for the Triaxial orbit library code.
! V2.0.1 : Changed quadrant grid bins to have more bins.
! V2.0.2 : fix 10**rlogmax in bin setup of qgrid_setup
!          add internal moments
! V2.0.3 : change the zeroth moment grid
! V2.0.4 : add zero psf routine
!          fixed projection
! V3.0.0 : Make intrinisic grid size recipe to be able to numerically  
!          Fit the intrinsic mass.
!     RvdB, Leiden, oktober/2005
!
!######################################################################

! $Id: orblibprogram.f90,v 1.1 2010/03/17 12:56:57 bosch Exp $

! Written by Remco van den Bosch <bosch@strw.leidenuniv.nl>
! Sterrewacht Leiden, The Netherlands


!######################################################################
!######################################################################
!######################################################################

program counterrotation
  use numeric_kinds
  use high_level, only : setup,run,stob
  implicit none
  real (kind=dp) :: t1,t2
  print*,"  ** Triaxial Orbit library by Remco C.E. van den Bosch <bosch@strw.leidenuniv.nl>"
  print*,"  * Jan 2007 "
  print*,"  * $Id: orblibprogram.f90,v 1.1 2010/03/17 12:56:57 bosch Exp $"  

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
