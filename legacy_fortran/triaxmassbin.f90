program triaxmassbin
    use binmass
    use initial_parameters
    use triaxpotent
    use numeric_kinds
    implicit none

    call iniparam()
    call tp_setup()

    call binmass_main()

end program Triaxmassbin
