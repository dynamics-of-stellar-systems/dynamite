program triaxmassbin
    use binmass
    use initial_parameters
    use triaxpotent
    use numeric_kinds
    implicit none

    call iniparam_bar()
    call tp_setup_bar()

    call binmass_main()

end program Triaxmassbin
