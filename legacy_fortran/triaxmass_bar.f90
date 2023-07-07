program triaxmass
    use intrinic_mass
    use initial_parameters
    use triaxpotent
    use numeric_kinds
    implicit none

    call iniparam_bar()

    print *, "Give the filename of the orbit library"
    read (unit=*, fmt="(a512)"), orblib_filename

    call tp_setup_bar()

    call intrin_radii()
    call intrin_spher()

end program Triaxmass
