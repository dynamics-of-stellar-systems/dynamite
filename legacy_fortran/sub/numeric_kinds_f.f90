!#######################################################################
! This Fortran 90 code is compatible with the F language subset
! a compiler can be found from http://www.fortran.com/imagine1
!#######################################################################

module numeric_kinds

    implicit none

    integer, parameter, public :: &
        i4b = selected_int_kind(9), &
        i2b = selected_int_kind(4), &
        i1b = selected_int_kind(2), &
        sp = kind(1.0), &
        dp = selected_real_kind(2*precision(1.0_sp)), &
        lgt = kind(.true.)
    real(kind=dp), parameter, public :: &
        pi = 3.141592653589793238462643383279502884197_sp, &
        pio2 = 1.57079632679489661923132169163975144209858_sp, &
        twopi = 6.283185307179586476925286766559005768394_sp
    real(kind=dp), parameter, public :: &
        pi_d = 3.141592653589793238462643383279502884197_dp, &
        pio2_d = 1.57079632679489661923132169163975144209858_dp, &
        twopi_d = 6.283185307179586476925286766559005768394_dp

end module numeric_kinds
