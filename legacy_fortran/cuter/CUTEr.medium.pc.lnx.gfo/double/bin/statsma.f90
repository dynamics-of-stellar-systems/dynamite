!     ( Last modified on Thu Aug 11 15:46:30 EDT 2005 )
Program STATSMA
  !
  Use CUTEr_precis
  Use CUTEr_interfaces
  !
  !  Statistics-collecting package for gathering information on
  !  variables and constraints from SIF problems.
  !
  !  D. Orban, August 2005
  !
  Implicit None

  ! Information pertaining to variables
  Type :: CUTEr_var_type
     Integer :: nvar
     Integer :: nfixed, nbelow, nabove, n2sided, nfree, nbnds
  End Type CUTEr_var_type

  ! Information pertaining to constraints
  Type :: CUTEr_con_type
     Integer :: ncon
     Integer :: nlle, nlge, nlrange, nleq, nlin
     Integer :: nnlle, nnlge, nnlrange, nnleq, nnlin
  End Type CUTEr_con_type

  Type :: CUTEr_db_type
     Integer :: nfr, nfx, n1s, n2s  ! # free, fixed, 1-sided, 2-sided bounds
     Integer :: m1sl, m2sl, meql ! # 1 and 2-sided ineq, equalities (linear)
     Integer :: m1sg, m2sg, meqg ! # 1 and 2-sided ineq, equalities (general)
     Character( Len = 86 ) :: var_str
     Character( Len = 76 ) :: lcon_str, gcon_str
     Character( Len = 99 ) :: classf_db_str
  End Type CUTEr_db_type

  ! Dynamic allocation flag
  Integer :: alloc_stat

  Integer, Parameter :: input = 47, iout = 6
  Integer :: i

  Real( Kind = wp ), Dimension( : ), Allocatable :: x, bl, bu, v, cl, cu
  Real, Dimension( 2 ) :: CPU( 2 )
  Real, Dimension( 7 ) :: CALLS( 7 )
  Character( len = 10 ) :: pname
  Logical :: efirst, lfirst, nvfrst
  Logical, Dimension( : ), Allocatable :: equatn, linear
  Logical :: constrained
  !
  Type( CUTEr_var_type ) :: vars
  Type( CUTEr_con_type ) :: cons
  Type( CUTEr_db_type  ) :: db
  !
  !  Open the relevant problem file.
  !
  Open ( input, FILE = 'OUTSDIF.d', FORM = 'FORMATTED', STATUS = 'OLD' )
  Rewind input
  !
  !  Get problem dimensions and determine which tools to use
  !
  Call CDIMEN( input, vars%nvar, cons%ncon )
  If( cons%ncon < 0 ) Then
     Close( input )
     Write( iout, '(A)' ) 'Error reading OUTSDIF.d'
     Write( iout, '(A21,1X,I6,1X,A3)' ) &
          'Number of constraints', cons%ncon, '< 0'
     Stop
  Endif
  constrained = (cons%ncon > 0)
  !
  !  Set up parameters
  !
  efirst = .False. ; lfirst = .True. ; nvfrst = .False.
  !
  !  Allocate arrays to hold problem data
  !
  Allocate( x( vars%nvar ),  STAT = alloc_stat )
  If( alloc_stat /= 0 ) Then
     Write( iout, 3000 ) 'X', vars%nvar
     Goto 900
  Endif
  Allocate( bl( vars%nvar ), STAT = alloc_stat )
  If( alloc_stat /= 0 ) Then
     Write( iout, 3000 ) 'BL', vars%nvar
     Goto 900
  Endif
  Allocate( bu( vars%nvar ), STAT = alloc_stat )
  If( alloc_stat /= 0 ) Then
     Write( iout, 3000 ) 'BU', vars%nvar
     Goto 900
  Endif
  If( constrained ) Then
     Allocate( v( cons%ncon ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'V', cons%ncon
        Goto 900
     End If
     Allocate( cl( cons%ncon ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'CL', cons%ncon
        Goto 900
     End If
     Allocate( cu( cons%ncon ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'CU', cons%ncon
        Goto 900
     End If
     Allocate( equatn( cons%ncon ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'EQUATN', cons%ncon
        Goto 900
     End If
     Allocate( linear( cons%ncon ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'LINEAR', cons%ncon
        Goto 900
     End If
     !
     !  If all ok, initialize problem data
     !
     Call CSETUP( input, iout, vars%nvar, cons%ncon,             &
                  x, bl, bu, vars%nvar, equatn, linear,          &
                  v, cl, cu, cons%ncon+1, efirst, lfirst, nvfrst )
  Else
     Allocate( equatn( 1 ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'EQUATN', 1
        Goto 900
     End If
     Allocate( linear( 1 ), STAT = alloc_stat )
     If( alloc_stat /= 0 ) Then
        Write( iout, 3000 ) 'LINEAR', 1
        Goto 900
     End If
     Call USETUP( input, iout, vars%nvar, x, bl, bu, vars%nvar )
  Endif
  !
  !  Obtain problem name.
  !
  Call PBNAME( vars%nvar, pname )
  !
  !  Initialize data on variables
  !
  vars%nfixed  = 0 ; vars%nbelow  = 0
  vars%nabove  = 0 ; vars%n2sided = 0
  vars%nbnds   = 0 ; vars%nfree   = 0 
  !
  !  Initialize data on constraints
  !
  cons%nlle    = 0 ; cons%nlge     = 0
  cons%nleq    = 0 ; cons%nlin     = 0
  cons%nnlle   = 0 ; cons%nnlge    = 0
  cons%nnleq   = 0 ; cons%nnlin    = 0
  cons%nlrange = 0 ; cons%nnlrange = 0
  !
  !  Obtain info on the variables
  !
  Do i = 1, vars%nvar
     If( (bl(i) > -INFTY .Or. bu(i) < INFTY) .And. bl(i) /= bu(i) ) &
          vars%nbnds = vars%nbnds + 1
     If( bl(i) > -INFTY .And. bu(i) >= INFTY ) Then
        vars%nbelow = vars%nbelow + 1
     Else If( bl(i) <= -INFTY .And. bu(i) < INFTY ) Then
        vars%nabove = vars%nabove + 1
     Else If( bl(i) > -INFTY .And. bu(i) < INFTY ) Then
        If( bl(i) == bu(i) ) Then
           vars%nfixed = vars%nfixed + 1
        Else
           vars%n2sided = vars%n2sided + 1
        End If
     Else
        vars%nfree = vars%nfree + 1
     End If
  End Do
  !
  !  Obtain info on the constraints
  !
  Do i = 1, cons%ncon
     If( linear(i) ) Then
        ! Process linear constraint
        If( equatn(i) ) Then
           cons%nleq = cons%nleq + 1
        Else
           If( cl(i) > -INFTY .And. cu(i) >= INFTY ) Then
              cons%nlge = cons%nlge + 1
           Else If( cl(i) <= -INFTY .And. cu(i) < INFTY ) Then
              cons%nlle = cons%nlle + 1
           Else If( cl(i) > -INFTY .And. cu(i) < INFTY ) Then
              cons%nlrange = cons%nlrange + 1
           Else
              Write( iout, * ) 'Oops!'
              Goto 900
           End If
        End If
     Else
        ! Process nonlinear constraint
        If( equatn(i) ) Then
           cons%nnleq = cons%nnleq + 1
        Else
            If( cl(i) > -INFTY .And. cu(i) >= INFTY ) Then
              cons%nnlge = cons%nnlge + 1
           Else If( cl(i) <= -INFTY .And. cu(i) < INFTY ) Then
              cons%nnlle = cons%nnlle + 1
           Else If( cl(i) > -INFTY .And. cu(i) < INFTY ) Then
              cons%nnlrange = cons%nnlrange + 1
           Else
              Write( iout, * ) 'Oops!'
              Goto 900
           End If
        End If
     End If
  End Do
  !
  !  Assemble classification strings
  !  We use the scheme
  !
  !  * CLASSIFICATION  probname XXXn-XX-n-m
  !  * CLASSIFICATION  probname VARIABLES           nfx  n1s  n2s
  !  * CLASSIFICATION  probname LINEAR CONSTRAINTS  m1sl m2sl meql
  !  * CLASSIFICATION  probname GENERAL CONSTRAINTS m1sg m2sg meqg
  db%nfr = vars%nfree
  db%nfx = vars%nfixed
  db%n1s = vars%nbelow + vars%nabove
  db%n2s = vars%n2sided
  Write( db%var_str, '(A18,A8,A21,I9,3(1X,I9))' )              &
       '* CLASSIFICATION  ', pname(1:8), ' VARIABLES           ', &
       db%nfr, db%n1s, db%n2s, db%nfx

  db%m1sl = cons%nlle + cons%nlge
  db%m2sl = cons%nlrange
  db%meql = cons%nleq
  Write( db%lcon_str, '(A18,A8,A21,I9,1X,I9,1X,I9)' )             &
       '* CLASSIFICATION  ', pname(1:8), ' LINEAR CONSTRAINTS  ', &
       db%m1sl, db%m2sl, db%meql

  db%m1sg = cons%nnlle + cons%nnlge
  db%m2sg = cons%nnlrange
  db%meqg = cons%nnleq
  Write( db%gcon_str, '(A18,A8,A21,I9,1X,I9,1X,I9)' )             &
       '* CLASSIFICATION  ', pname(1:8), ' GENERAL CONSTRAINTS ', &
       db%m1sg, db%m2sg, db%meqg
  Write( db%classf_db_str, '(I9,9(1X,I9))' )     &
       db%nfr,  db%n1s,  db%n2s, db%nfx, &
       db%m1sl, db%m2sl, db%meql,        &
       db%m1sg, db%m2sg, db%meqg
  
  !Write( iout, * ) 'Database classification strings:'
  !Write( iout, '(A86)' ) db%var_str
  !Write( iout, '(A76)' ) db%lcon_str
  !Write( iout, '(A76)' ) db%gcon_str
  Write( iout, '(A99)' ) db%classf_db_str
  Goto 900
  !
  !  Display collected data
  !
  Write ( IOUT, 2000 ) pname, vars%nvar, vars%nbnds, vars%nbelow,     &
       vars%nabove, vars%n2sided, vars%nfixed, vars%nfree,            &
       cons%ncon,                                                     &
       cons%nlin, cons%nleq, cons%nlle, cons%nlge, cons%nlrange,      &
       cons%nnlin, cons%nnleq, cons%nnlle, cons%nnlge, cons%nnlrange, &
       cons%nleq + cons%nnleq,                                        &
       cons%nlle+cons%nlge+2*cons%nlrange+cons%nnlle+cons%nnlge+2*cons%nnlrange
  !
  !  Terminate
  !
900 Continue
  !
  !  Close the problem file
  !
  Close( input )
  !
  !  Free allocated memory
  !
  If( Allocated( x ) ) Deallocate( x )
  If( Allocated( bl ) ) Deallocate( bl )
  If( Allocated( bu ) ) Deallocate( bu )
  If( Allocated( equatn ) ) Deallocate( equatn )
  If( Allocated( linear ) ) Deallocate( linear )
  If( constrained ) Then
     If( Allocated( v ) ) Deallocate( v )
     If( Allocated( cl ) ) Deallocate( cl )
     If( Allocated( cu ) ) Deallocate( cu )
  End If
  !
  !  Exit
  !
  Stop
  !
  !  Non-executable statements.
  !
  !  The following is the complete standard statistics output format: select
  !  the items that are relevant to the type of problems solved and adapt the
  !  name of the code. It is broken in two to comply with compilers
  !  which want to see no more than 19 continuation lines.
  !
2000 Format( /, 24('='), ' Problem statistics ', 24('=') //, &
          ' Code used                :       STATS',    /, &
          ' Problem                  :      ', A10,     /, &
          ' # variables              =      ', I10 /, &
          ' #   bounded              =      ', I10 /, &
          ' #     below only         =      ', I10 /, &
          ' #     above only         =      ', I10 /, &
          ' #     below and above    =      ', I10 /, &
          ' #   fixed                =      ', I10 /, &
          ' #   free                 =      ', I10 /, &
          ' # constraints            =      ', I10 /, &
          ' #   linear               =      ', I10 /, &
          ' #     equalities         =      ', I10 /, &
          ' #     <= inequalities    =      ', I10 /, &
          ' #     >= inequalities    =      ', I10 /, &
          ' #     range              =      ', I10 /, &
          ' #   nonlinear            =      ', I10 /, &
          ' #     equalities         =      ', I10 /, &
          ' #     <= inequalities    =      ', I10 /, &
          ' #     >= inequalities    =      ', I10 /, &
          ' #     range              =      ', I10 /, &
          ' # equality constraints   =      ', I10 /, &
          ' # inequality constraints =      ', I10 /, &
          '   (ranges count as two)',              /, &
          68('=') )
3000 Format( /, 'Error allocating array ', A6, ', dim = ', I6 )

End Program STATSMA

