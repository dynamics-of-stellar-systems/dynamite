! The sub program to be insert in orblib-f.f90
! It calculate and stores the probability (pij) for discrete data i contributed by orbit j
! First we only include los velocity, so the discrete data has (x,y,losv), not proper motion at the start



module pij_calculate
  private
  real (kind=dp),private,allocatable,dimension(:) :: pj
! store pij orbits by orbits, pj
contains

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine pj_setup
    use read_discrete, only: num_disc
    allocate (pj(num_disc))
  end subroutine pj_setup

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine pij_cal(proj,losvel)
    use read_discrete, only: num_disc, x_disc, y_disc, losv_disc, losve_disc, s_pos, s_v ! read in the number of data, and ,x,y,losv,losv_err, smoothing in position, smoothing in velocity
! refer to qgrid quadrant_light

  end subroutine pij_cal

end module pij_calculate
 
