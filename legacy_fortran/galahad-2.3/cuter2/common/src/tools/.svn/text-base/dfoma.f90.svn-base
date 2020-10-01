!
!  Main CUTEr driver for DFO 2.0.0.
!  Original version by K. C. Dang, 2008
!  Fortran 90 version by D. Orban, 2009
!
Program DFOMA
  Use CUTEr_precis
  Use CUTEr_interfaces
  Implicit None
  !
  !  Variable declarations
  !
  Integer :: N, M, NCLIN, NCNLN, NLIN, NEQ, NBNDS
  Real(Kind = wp), Dimension(:), Allocatable :: X0, BL, BU, V, CL, CU
  Logical, Dimension(:), Allocatable :: EQUATN, LINEAR
  Logical :: EFIRST, LFIRST, NVFRST, IFINIV, CONSTRAINED
  Character(Len = 10) :: PNAME
  Character(Len = 10), Dimension(:), Allocatable :: VNAMES, GNAMES
  !     - Variables for I/O
  Integer, Parameter :: IOUT=6, INPUT=47, INSPEC=198, REPRTOUT=1812
  !     - Variables for algorithm parameter
  Integer :: NX, MAXIT, MAXNF,STPCRTR, IPRINT, SCALE
  !      LOGICAL IFINTV
  Real(Kind = wp) :: DELMIN, DELTA, CNSTOL, PP, STPTHR
  !     - Variables for CUTEr report
  Real :: CPU(2), CALLS(7)
  !     - Variables for working space
  Integer :: LDA
  Integer :: IT, NF, INFO
  Integer :: I
  Real(Kind = wp) :: F0
  Real(Kind = wp), Dimension(:), Allocatable :: X, FX, C, CONX, LB, UB, ALIN
  !
  !  Open data file.
  !
  Open(INPUT, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD')
  Rewind INPUT
  !
  !  Allocate working vectors.
  !
  CONSTRAINED = .False.
  Call CDIMEN(INPUT, N, M)
  If( M > 0 ) Then
     CONSTRAINED = .True.
  Else If( M < 0 ) Then
     Write( 6, '(A)' ) 'Error reading OUTSDIF.d'
     Stop
  Endif
  !
  !  Set up SIF data.
  !
  EFIRST = .False. ; LFIRST = .True. ; NVFRST = .False.
  Allocate(X0(N)) ; Allocate(BL(N)) ; Allocate(BU(N))
  If( CONSTRAINED ) Then
     Allocate(C(M))
     Allocate(V(M), CL(M), CU(M), EQUATN(M), LINEAR(M) )
     Call CSETUP( INPUT, IOUT, N, M, X0, BL, BU, N, EQUATN, &
          LINEAR, V, CL, CU, M, EFIRST, LFIRST, NVFRST )
  Else
     Call USETUP( INPUT, IOUT, N, X0, BL, BU, N )
  Endif
  !
  !  Obtain problem, variables and constraint names.
  !
  Allocate(VNAMES(N))
  If( CONSTRAINED ) Then
     Allocate(GNAMES(M))
    Call CNAMES(N, M, PNAME, VNAMES, GNAMES)
  Else
     Call UNAMES(N, PNAME, VNAMES)
  Endif
  !
  !  Obtain info on the problem
  !
  NLIN  = 0 ; NEQ   = 0 ; NBNDS = 0
  If( CONSTRAINED ) Then
     Call GETINFO( N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
  Else
     Call GETINFO( N, 0, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS )
  Endif
  !
  !  Treat all constraints as derivative free (ignore linear constraints)
  !
  NCLIN = 0
  NCNLN = 0
  LDA = 1
  Allocate(LB(N + NCLIN + NCNLN + M))
  Allocate(UB(N + NCLIN + NCNLN + M))
  Call DCOPY(N, BL, 1, LB(1:N), 1)
  Call DCOPY(N, BU, 1, UB(1:N), 1)
  If( CONSTRAINED ) Then
     Call DCOPY(M, CL, 1, LB(N+1:N+M), 1)
     Call DCOPY(M, CU, 1, UB(N+1:N+M), 1)
  Endif
  !
  !  Read algorithm specification
  !
  Open(INSPEC, FILE = 'DFO.SPC', FORM = 'FORMATTED', STATUS='OLD')
  Rewind INSPEC
  Read(INSPEC, 1000) NX, MAXIT, MAXNF, STPCRTR, DELMIN, STPTHR, CNSTOL, &
       DELTA, PP, SCALE, IPRINT
  Close(INSPEC)
  If( NX >= 2 ) Then
     Write(IOUT,*) 'NX >= 2 not currently supported; Check spec file.'
     Goto 999
  Endif
  !
  !  Allocate space for initial points supplied
  !
  Allocate(FX(NX))
  Allocate(X(N*NX))
  If( CONSTRAINED ) Allocate(CONX(M*NX))
  Allocate(ALIN(1))     ! Not used
  IFINIV = (NX >= 2)
  !
  !  Evaluate initial objective and constraint values
  !
  If( CONSTRAINED ) Then
     Call CFN(N, M, X0, F0, M, C)
  Else
     Call UFN(N, X0, F0)
  Endif
  If( NX == 1 ) Then
     Call DCOPY(N, X0, 1, X(1:N), 1)
     FX(1) = F0
     If( CONSTRAINED ) Call DCOPY(M, C, 1, CONX(1:M), 1)
  Else
     If( CONSTRAINED ) Then
        Do I = 1, NX
           Call CFN(N, M, X((I-1)*N + 1:I*N), FX(I), M, CONX((I-1)*M + 1:I*M))
        End Do
     Else
        Do I = 1, NX
           Call UFN(N, X((I-1)*N + 1:I*N), FX(I))
        End Do
     Endif
  Endif
  !
  !  Call main DFO routine
  !
  Call DFO(N, NX, X, N, FX, CONX, IFINIV, M, C ,NCLIN, NCNLN, LB, UB ,  &
       ALIN , LDA , VNAMES, PNAME, GNAMES, IT, NF, INFO, MAXIT,  MAXNF, &
       STPCRTR, DELMIN, STPTHR, CNSTOL, DELTA, PP, SCALE,  IOUT , IPRINT)
  !
  !  Write out statistics
  !
  Call CREPRT(CALLS, CPU)
  Write(IOUT, 2000) PNAME, N, M, NLIN, NEQ, M-NEQ, NBNDS, CALLS(1), CALLS(5)
  Write(IOUT, 2001) INFO, FX(1), CPU(1), CPU(2)
  !
  !  Write select statistics to file
  !
  Open(REPRTOUT, FILE = 'cuter.log', FORM = 'FORMATTED')
  Write(REPRTOUT, 2002) PNAME, N, M, IT, Int(CALLS(1)), Int(CALLS(5)), &
       INFO, FX(1), F0, CPU(1), CPU(2)
  Close(REPRTOUT)
999 Continue
  Close(INPUT)
  !
  !  Free allocated space
  !
  Deallocate(X0)
  Deallocate(BL)
  Deallocate(BU)
  Deallocate(VNAMES)
  If( CONSTRAINED ) Then
     Deallocate(EQUATN)
     Deallocate(LINEAR)
     Deallocate(V)
     Deallocate(CL)
     Deallocate(CU)
     Deallocate(GNAMES)
     Deallocate(C)
     Deallocate(CONX)
  Endif
  Deallocate(LB)
  Deallocate(UB)
  Deallocate(X)
  Deallocate(FX)
  Deallocate(ALIN)
  !
  !  Non-executable statements
  !
1000 Format(4(I10, /), 5(D10.3, /), 1(I10, /), I10)
2000 Format( /, 24('*'), ' CUTEr statistics ', 24('*') // &
          ,' Code used                :  DFO',    / &
          ,' Problem                  :  ', A10,    / &
          ,' # variables              =      ', I10 / &
          ,' # constraints            =      ', I10 / &
          ,' # linear constraints     =      ', I10 / &
          ,' # equality constraints   =      ', I10 / &
          ,' # inequality constraints =      ', I10 / &
          ,' # bounds                 =      ', I10 / &
          ,' # objective functions    =        ', F8.2 / &
          ,' # constraints functions  =        ', F8.2 )

2001 Format(                                          &
          ' Exit code                =      ', I10 / &
          ,' Final f                  = ', D15.7 / &
          ,' Set up time              =      ', 0P, F10.2, ' seconds'/ &
          ' Solve time               =      ', 0P, F10.2, ' seconds'// &
          66('*') / )

2002 Format(' Code used               : RNAME    : C : DFO',/ &
       ,' Problem                 : PNAME    : C : ', A15 ,/ &
       ,' # variables             : NVAR     : I : ', I15 ,/ &
       ,' # constraints           : NCON     : I : ', I15 ,/  &
       ,' # iterations            : NITER    : I : ', I15 ,/  &
       ,' # objective functions   : NFEVAL   : I : ', I15 ,/ &
       ,' # objective gradients   : NCEVAL   : I : ', I15 ,/  &
       ,' Exit code               : EXITCODE : I : ', I15 ,/ &
       ,' Final f                 : FVAL     : F : ', E15.7 ,/ &
       ,' Initial f               : FZERO    : F : ', E15.7 ,/ &
       ,' Set up time (in second) : PTIME    : F : ', 0P,F15.7 ,/ &
       ,' Solve time (in second)  : STIME    : F : ', 0P,F15.7 ,/)
End Program DFOMA

!==============================================================================

Subroutine FUN(N, M, X, F, C, IFERR)
  !
  !  Evaluate objective and constraint values at X
  !
  Use CUTEr_precis
  Implicit None

  Integer, Intent(In) :: N, M
  Real(Kind=wp), Dimension(N), Intent(In) :: X
  Real(Kind=wp), Dimension(M), Intent(Out) :: C
  Real(Kind=wp), Intent(Out) :: F
  Logical, Intent(Out) :: IFERR
  !Intrinsic :: ISNAN           ! Only in GFortran >= 4.3

  Integer :: i

  IFERR = .False.

  If( M > 0 ) Then
     Call CFN(N,M,X,F,M,C)
  Else
     Call UFN(N,X,F)
  Endif
  If( F /= F ) Then             ! If( ISNAN(F) ) Then
     IFERR = .True.
     Write(6,*) 'CUTEr : Function value is NaN!'
     Write(6,*) 'X = ', X
     Goto 3000
  Endif
  Do i = 1, M
      If( C(i) /= C(i) ) Then   ! If( ISNAN(C(i)) ) Then
        IFERR = .True.
        Write(6,*) 'CUTEr : Constraint value is NaN!'
        Write(6,*) 'X = ', X
        Goto 3000
     Endif
  End Do
3000 Continue
  Return
End Subroutine FUN

!==============================================================================

  Subroutine GETINFO(N, M, BL, BU, EQUATN, LINEAR, NLIN, NEQ, NBNDS)
    !
    ! Input/Output variables
    !
    Use CUTEr_precis
    Implicit None
    Integer, Intent( IN  ) :: N, M
    Integer, Intent( OUT ) :: NLIN, NEQ, NBNDS
    Real( Kind = wp ), Dimension( N ), Intent( IN ) :: BL, BU
    Logical, Dimension( M ), Intent( IN ) :: EQUATN, LINEAR
    !
    !     Local variables
    !
    Integer :: I

    NLIN  = 0 ; NEQ   = 0 ; NBNDS = 0

    Do I = 1, M
       If( EQUATN( I ) ) NEQ  = NEQ  + 1
       If( LINEAR( I ) ) NLIN = NLIN + 1
    End Do

    Do I = 1, N
       If( BL(I) > -INFTY .And. BU(I) < INFTY ) Then
          NBNDS = NBNDS + 2
       Else If( BL(I) > -INFTY .Or. BU(I) < INFTY ) Then
          NBNDS = NBNDS + 1
       Endif
    End Do
  End Subroutine GETINFO

!==============================================================================
