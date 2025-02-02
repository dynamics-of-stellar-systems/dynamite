program orbitstart_run
    use orbitstart
    use initial_parameters
    use interpolpot
    use triaxpotent

    implicit none

    integer :: r_seed
    double precision :: ran1, r_num

! read random number seed, any value <= 0 results in stochastic seed
    read (unit=*, fmt=*) r_seed
    if (r_seed <= 0) then
        call random_seed()
        call random_number(r_num)
        r_seed = 2147483647*r_num
        r_num = ran1(-r_seed)
        write (*, *) "Using ran1 with random seed ", r_seed
    else
        r_num = ran1(-r_seed)
        write (*, *) "Using ran1 with given seed ", r_seed
    end if

    call iniparam()
    call ip_setup()
    call runorbitstart()
    call ip_stop()
end program orbitstart_run
