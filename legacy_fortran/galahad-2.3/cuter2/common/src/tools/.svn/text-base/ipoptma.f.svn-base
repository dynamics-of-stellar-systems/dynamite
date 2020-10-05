C Copyright (C) 2002, Carnegie Mellon University, Dominique Orban and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
      PROGRAM           IPOPTMA
C
C     IPOPT CUTEr driver.
C     D. Orban,  adapted from Andreas Wachter's CUTE driver.
C
      IMPLICIT NONE
      INTEGER IOUT
      PARAMETER( IOUT = 6 )
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/



      INTEGER LRW, LIW
      PARAMETER( LRW=0, LIW=0 )
      INTEGER IW
      DOUBLE PRECISION RW


      INTEGER N, M, NLB, NUB
      DOUBLE PRECISION X( CUTE_NMAX )
      INTEGER ILB( CUTE_NMAX )
      INTEGER IUB( CUTE_NMAX )
      DOUBLE PRECISION BNDS_L( CUTE_NMAX )
      DOUBLE PRECISION BNDS_U( CUTE_NMAX )
      DOUBLE PRECISION V_L( CUTE_NMAX )
      DOUBLE PRECISION V_U( CUTE_NMAX )
      DOUBLE PRECISION LAM( CUTE_MMAX )
      DOUBLE PRECISION C( CUTE_MMAX )
      DOUBLE PRECISION IPOPT
C
C     Algorithmic Parameters (INITPARAMS)
C
      INTEGER NARGS
      DOUBLE PRECISION ARGS( 50 )
      CHARACTER*20 CARGS( 50 )

      INTEGER ITER
      INTEGER IERR

      EXTERNAL EVAL_F
      EXTERNAL EVAL_C
      EXTERNAL EVAL_G
      EXTERNAL EVAL_A
      EXTERNAL EVAL_H
      EXTERNAL EVAL_HESSLAG_V
      EXTERNAL EVAL_HESSOBJ_V
      EXTERNAL EVAL_HESSCON_V
C
C     The following arrays could be used instead of common blocks to
C     communicate between driver routines and evaluation subroutines
C
      DOUBLE PRECISION DAT(1)
      INTEGER IDAT(1)

      INTEGER fevals, cevals
      COMMON /EVALS/ fevals, cevals

      REAL CALLS( 7 ), CPU( 2 )
      CHARACTER*10 PNAME
      CHARACTER*10 VNAMES( CUTE_NMAX ), GNAMES( CUTE_MMAX )
      DOUBLE PRECISION FinalF, DummyG, cnrm
      INTEGER IDAMAX
C
      fevals = 0
      cevals = 0
C
C     Get problem dimensions and initialize
C
      CALL CUTE_INIT(N, M, CUTE_NMAX, X, NLB, ILB, BNDS_L,
     .     NUB, IUB, BNDS_U)
C
C     Get problem name.
C
      CALL CNAMES( CUTE_N, CUTE_M, PNAME, VNAMES, GNAMES )
C
C     Set algorithmic parameters (None :)
C
      NARGS = 0
C
C     Call IPOPT
C
      FinalF = IPOPT(N, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L,
     1     V_U, LAM, C, LRW, RW, LIW, IW, ITER, IERR, EVAL_F, EVAL_C,
     2     EVAL_G, EVAL_A, EVAL_H, EVAL_HESSLAG_V, EVAL_HESSOBJ_V,
     3     EVAL_HESSCON_V, DAT, IDAT, NARGS, ARGS, CARGS)
C
C     Display CUTEr statistics
C
      CALL CREPRT( CALLS, CPU )
      WRITE ( IOUT, 2000 ) PNAME, CUTE_N, CUTE_M, CALLS(1), CALLS(2),
     .     CALLS(3), CALLS(4), CALLS(5), CALLS(6), CALLS(7),
     .     IERR, FinalF, CPU(1), CPU(2)
c
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *     ,/,' Code used               :  IPOPT',    /
     *     ,' Problem                 :  ', A10,    /
     *     ,' # variables             =      ', I10 /
     *     ,' # constraints           =      ', I10 /
     *     ,' # objective functions   =      ', E15.7 /
     *     ,' # objective gradients   =      ', E15.7 /
     *     ,' # objective Hessians    =      ', E15.7 /
     *     ,' # Hessian-vector prdct  =      ', E15.7 /
     *     ,' # constraints functions =      ', E15.7 /
     *     ,' # constraints gradients =      ', E15.7 /
     *     ,' # constraints Hessians  =      ', E15.7 /
     *     ,' Exit code               =      ', I10 /
     *     ,' Final f                 = ', E15.7 /
     *     ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ,' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     ,/,66('*') / )

 9999 CONTINUE
      END



C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_A(TASK, N, X, NZ, A, ACON, AVAR, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_a.F,v 1.2 2004/03/11 01:27:50 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute Jacobian of constraints to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   TASK      I    INT     =0: Obtain NZ
CP                         <>0: Compute Jacobian
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where A is to be evaluated
CP   NZ       I/O   INT    TASK = 0: O: number of nonzero elements
CP                         otherwise: number of nonzero elements
CP                                     (size of A, AVAR, ACON)
CP   A         O    DP     (only TASK<>0) values in Jacobian
CP   ACON      O    INT    (only TASK<>0) row indices
CP   AVAR      O    INT    (only TASK<>0) column indices
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CCS    CDIMSJ
CCS    CCFSG
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer N
      double precision X(N)
      integer NZ
      double precision A(NZ)
      integer ACON(NZ)
      integer AVAR(NZ)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision c(CUTE_MMAX)
      double precision cjac(CUTE_NZMAX)
      integer indvar(CUTE_NZMAX), indfun(CUTE_NZMAX)
      integer i, nnzj
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( TASK.eq.0 ) then
C
C     Call CDIMSJ to obtain number of nonzero elements
C
         CALL CDIMSJ( NZ )
C
C     Substract contribution of (dense) objective function gradient
C
         NZ = NZ - CUTE_N
         NZ = NZ + CUTE_NIQ
      else
C
C     Call CCFSG to obtain Jacobian for constraints
C
         call CCFSG(CUTE_N, CUTE_M, X, CUTE_MMAX, c, nnzj,
     1        NZ, A, AVAR, ACON, .TRUE.)
C
C     Augment entries for slacks
C
         do i = 1, CUTE_NIQ
            A   (nnzj+i) = -1.d0
            AVAR(nnzj+i) = CUTE_N + i
            ACON(nnzj+i) = CUTE_IIQ(i)
         enddo

      endif

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_C(N, X, M, C, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_c.F,v 1.2 2004/03/11 01:27:51 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute values of constraints to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      07/01/99   BUG: problems if ineq not first
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where G is to be evaluated
CP   M         I    INT    number of constraints
CP   C         O    DP     values of constraints
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CCS    CCFG
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision X(N)
      integer M
      double precision C(M)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision f, dummy
      integer i, liq, leq
      logical ineq
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Call CCFG to obtain constraint values, but without slacks
C
      call CCFG(CUTE_N, CUTE_M, X, M, C, .FALSE., 1, 1, dummy, .FALSE.)
C
C     Add entries for slack variables and constant terms
C
      liq = 1
      leq = 1
      do i = 1, M
         ineq = .false.
         if( liq.le.CUTE_NIQ ) then
            if( CUTE_IIQ(liq).eq.i ) then
               C(i) = C(i) - X(CUTE_N+liq)
               liq = liq + 1
               ineq = .true.
            endif
         endif
         if( .not.ineq ) then
            C(i) = C(i) - CUTE_CEQ(leq)
            leq = leq + 1
         endif
      enddo

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_F(N, X, F, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_f.F,v 1.2 2004/03/11 01:27:51 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute objective function value to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where F is to be evaluated
CP   F         O    DP     objective function value
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
C     COFG
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision X(N)
      double precision F
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision dummy(1)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Call COFG to obtain value of objective function
C
      call COFG( CUTE_N, X, F, dummy, .false.)

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_G(N, X, G, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_g.F,v 1.2 2004/03/11 01:27:51 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute gradient of objective function to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where G is to be evaluated
CP   G         O    DP     gradient of objective function
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
C     COFG
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision X(N)
      double precision G(N)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision f
      integer i
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Call COFG to obtain gradient of objective function
C
      call COFG( CUTE_N, X, f, G, .true.)
C
C     Add entries for slack variables
C
      do i = CUTE_N + 1, N
         G(i) = 0.d0
      enddo

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C

      subroutine EVAL_H(TASK, N, X, M, LAM, NNZH, HESS, IRNH, ICNH,
     1     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_h.F,v 1.2 2004/03/11 01:27:52 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute Jacobian of constraints to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      03/23/00
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   TASK      I    INT     =0: Obtain NZ
CP                         <>0: Compute Jacobian
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where A is to be evaluated
CP   NZ       I/O   INT    TASK = 0: O: number of nonzero elements
CP                         otherwise: number of nonzero elements
CP                                     (size of A, AVAR, ACON)
CP   A         O    DP     (only TASK<>0) values in Jacobian
CP   ACON      O    INT    (only TASK<>0) row indices
CP   AVAR      O    INT    (only TASK<>0) column indices
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CCS    CDIMSH
CCS    CSH
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer N
      integer M
      integer NNZH
      double precision LAM(M)
      double precision X(N)
      double precision HESS(NNZH)
      integer IRNH(NNZH)
      integer ICNH(NNZH)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision v(CUTE_MMAX)
      double precision h(CUTE_NZMAX)
      integer indirnh(CUTE_NZMAX), indicnh(CUTE_NZMAX)
      integer nnzh2

      integer NNZH_STORE
      save    NNZH_STORE
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( TASK.eq.0 ) then
C
C     Get number of nonzeros in Hessian of the Lagrangian
C
         CALL CDIMSH( NNZH )
         if( NNZH.gt.CUTE_NZMAX ) then
            write(*,*) 'CUTE_NZMAX = ',CUTE_NZMAX,' too small in cute_h'
            write(*,*) 'Increase to at least ', NNZH
            stop
         endif
         NNZH_STORE = NNZH
      else
C
C     Call CSH to obtain Hessian for constraints
C
         call CSH(CUTE_N, CUTE_M, X, M, LAM, nnzh2, NNZH, HESS,
     1        IRNH, ICNH)
      endif

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine EVAL_HESSLAG_V(TASK, N, X, M, LAM, VIN, VOUT,
     1     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: eval_hesslag_v.F,v 1.2 2004/03/11 01:27:52 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute product of Hessian of Lagrangian (individual constraints
CT    weighted by LAM) with vector
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      11/05/00
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   TASK      I    INT     =0: reevaluate Hessian entries (because X or LAM
CP                              changed
CP                         <>0: no need to reevaluate
CP   N         I    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   X         I    DP     point where Hessians are to be evaluated
CP   M         I    INT    number of constraints (including slack equations)
CP   LAM       I    DP     weights for constraints Hessians
CP   VIN       I    DP     vector to be mutliplied with Hessian
CP   VOUT      O    DP     resulting product
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CCS    CPROD
CCS    DCOPY
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer N
      integer M
      double precision LAM(M)
      double precision X(N)
      double precision VIN(N)
      double precision VOUT(N)
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      logical goth
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( TASK.eq.0 ) then
         goth = .false.
      else
         goth = .true.
      endif

      call CPROD(CUTE_N, CUTE_M, goth, X, M, LAM, VIN, VOUT)
C
C     slacks only appear linear
C
      if( N.gt.CUTE_N ) then
         call DCOPY(N-CUTE_N, 0d0, 0, VOUT(CUTE_N+1), 1)
      endif

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine CUTE_INIT(N, M, LX, X, NLB, ILB, BNDS_L,
     1                     NUB, IUB, BNDS_U)
C
C*******************************************************************************
C
C    $Id: cute_init.F,v 1.2 2004/03/11 01:27:50 andreasw Exp $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Initialize interface to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      07/19/99  initialize slacks based on C
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CCD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         O    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   M         O    INT    number of equality constraints
CP   LX        I    INT    actual declared length of X, LIB, BNDS_L, IUB, BNDS_U
CP   X         O    DP     starting point
CP   NLB       O    INT    number of lower bounds
CP   ILB       O    INT    indices for lower bounds
CP                            ( BNDS_L(i) is lower bound for X(ILB(i)) )
CP   BNDS_L    O    INT    values of lower bounds
CP   NUB       O    INT    number of upper bounds
CP   IUB       O    INT    indices for upper bounds
CP                            ( BNDS_U(i) is lower bound for X(IUB(i)) )
CP   BNDS_U    O    INT    values of upper bounds
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CCS    CSETUP
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer M
      integer LX
      double precision X(LX)
      integer NLB
      integer ILB(LX)
      double precision BNDS_L(LX)
      integer NUB
      integer IUB(LX)
      double precision BNDS_U(LX)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision bl(CUTE_NMAX), bu(CUTE_NMAX)
      double precision v(CUTE_MMAX), cl(CUTE_MMAX), cu(CUTE_MMAX)
      double precision c(CUTE_MMAX), f, dummy
      logical equatn(CUTE_MMAX), linear(CUTE_MMAX)
      integer cnr_input, iout, i
      logical efirst, lfirst, nvfrst, ex
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
C-------------------------------------------------------------------------------

C
C     Call CSETUP to obtain problem size and starting point
C
      cnr_input = 60
      iout = 6
      efirst = .false.
      lfirst = .false.
      nvfrst = .false.

      open(cnr_input,file='OUTSDIF.d',status='old')

      call CSETUP(cnr_input, iout, CUTE_N, CUTE_M,
     1     X, bl, bu, CUTE_NMAX, equatn, linear, v, cl, cu,
     2     CUTE_MMAX, efirst, lfirst, nvfrst)
      close(cnr_input)
C
C     Added this:  Compute C in order to initialize slacks better
C
      call CCFG(CUTE_N, CUTE_M, X, M, C, .false., 1, 1, dummy, .false.)
C
      M = CUTE_M
      N = CUTE_N
C
C     Obtain bounds on variables
C
      NLB = 0
      do i = 1, CUTE_N
         if( bl(i).gt.CUTE_NOLB ) then
            NLB = NLB + 1
            ILB(NLB) = i
            BNDS_L(NLB) = bl(i)
         endif
      enddo

      NUB = 0
      do i = 1, CUTE_N
         if( bu(i).lt.CUTE_NOUB ) then
            NUB = NUB + 1
            IUB(NUB) = i
            BNDS_U(NUB) = bu(i)
         endif
      enddo
C
C     Find inequalities and augment X
C
      CUTE_NIQ = 0
      CUTE_NEQ = 0
      do i = 1, CUTE_M
         if( .not.equatn(i) ) then
            CUTE_NIQ = CUTE_NIQ + 1
            CUTE_IIQ(CUTE_NIQ) = i
            N = N + 1
            if( N.gt.LX ) then
               write(6,*) 'Error in cute_init: LX = ',
     1                    LX,' is too small. Abort.'
               stop
            endif
C
C           This is kind of arbitrary: initialize slack to zero...
C
C            X(N) = 0.d0
            X(N) = C(i)
            if( cl(i).gt.CUTE_NOLB ) then
               NLB = NLB + 1
               ILB(NLB) = N
               BNDS_L(NLB) = cl(i)
            endif
            if( cu(i).lt.CUTE_NOUB ) then
               NUB = NUB + 1
               IUB(NUB) = N
               BNDS_U(NUB) = cu(i)
            endif
         else
            CUTE_NEQ = CUTE_NEQ + 1
            CUTE_CEQ(CUTE_NEQ) = cl(i)
         endif
      enddo
C
C     For basis selection, write indices of slack variables into file
C     SLACKS.DAT, or delete this file, if no slack variables
C
      if( CUTE_NIQ.gt.0 ) then
         open(10,file='SLACKS.DAT',status='unknown')
         do i = CUTE_N+1, CUTE_N+CUTE_NIQ
            write(10,1000) i
 1000       format(i10)
         enddo
         close(10)
      else
         inquire(file='SLACKS.DAT',exist=ex)
         if( ex ) then
            open(10,file='SLACKS.DAT',status='old')
            close(10,status='delete')
         endif
      endif

 9999 continue
      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
      subroutine EVAL_HESSCON_V(TASK, N, X, M, LAM, VIN, VOUT,
     1     DAT, IDAT)
C
C    $Id: eval_hesscon_v.F,v 1.2 2004/03/11 01:27:52 andreasw Exp $
C
C     TASK = 0: reevaluate Hessian (because X or LAM changed)
C            1: do not need to reevaluate - but X and LAM are still
C               set to the correct values
      implicit none
      integer TASK, N, M
      double precision LAM(M), VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

      double precision lam2(CUTE_MMAX)
      double precision tmp(CUTE_NMAX)

      if( M.gt.CUTE_MMAX .or. N.gt.CUTE_NMAX) then
         write(*,*) 'N or M too large in eval_hesscon_v.'
         stop
      endif

      call EVAL_HESSLAG_V(0, N, X, M, LAM, VIN, VOUT, DAT, IDAT)

      call DCOPY(CUTE_M, 0d0, 0, lam2, 1)
      call EVAL_HESSLAG_V(0, N, X, M, lam2, VIN, tmp, DAT, IDAT)
      call DAXPY(CUTE_N, -1d0, tmp, 1, VOUT, 1)

      return
      end
C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
      subroutine EVAL_HESSOBJ_V(TASK, N, X, M, VIN, VOUT, DAT, IDAT)
C
C    $Id: eval_hessobj_v.F,v 1.2 2004/03/11 01:27:52 andreasw Exp $
C
C     TASK = 0: reevaluate Hessian (because X or LAM changed)
C            1: do not need to reevaluate - but X and LAM are still
C               set to the correct values
      implicit none
      integer TASK, N, M
      double precision  VIN(N), X(N), VOUT(N)
      double precision DAT(*)
      integer IDAT(*)

C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C
C     PARAMETER definitions and COMMON block for CUTEr interface
C

CB    CUTE_NMAX       maximal number of variables
CB    CUTE_MMAX       maximal number of constraints
CB    CUTE_NZMAX      maximal number of nonzero elements
CB    CUTE_NOLB       constant of CUTEr to indicate -inf as lower bound
CB    CUTE_NOUB       constant of CUTEr to indicate +inf as upper bound

      INTEGER CUTE_NMAX, CUTE_MMAX, CUTE_NZMAX
CTOY  PARAMETER( CUTE_NMAX = 1000,  CUTE_MMAX = 1000  )
CMED  PARAMETER( CUTE_NMAX = 10000, CUTE_MMAX = 10000 )
CBIG  PARAMETER( CUTE_NMAX = 50000, CUTE_MMAX = 50000 )
CCUS  PARAMETER( CUTE_NMAX = 200000, CUTE_MMAX = 200000 )
CTOY  PARAMETER( CUTE_NZMAX = 100000  )
CMED  PARAMETER( CUTE_NZMAX = 200000  )
CBIG  PARAMETER( CUTE_NZMAX = 1000000 )
CCUS  PARAMETER( CUTE_NZMAX = 10000000 )
      DOUBLE PRECISION CUTE_NOLB, CUTE_NOUB
      PARAMETER( CUTE_NOLB = -1.0D+20, CUTE_NOUB =  1.0D+20 )

CB    CUTE_N          number of variables in CUTEr problem
CB    CUTE_M          number of constraints in CUTEr problem
CB    CUTE_NIQ        number of inequality constraints
CB    CUTE_NEQ        number of equality constraints
CB    CUTE_IIQ        indices of inequality constraints
CB                       constraint number CUTE_IIQ(i) is inequality constraint
CB                       (always ordered increasingly)
CB    CUTE_CEQ        constants for equality constraints
CB                       (ordered as constraints without inequalities)
CB
CB    Structure of reformulation:
CB
CB    X (variables after reformulation) has a first entries all variables
CB    from original CUTEr problem, followed by the slack variables for the
CB    inequality constraints.

      integer CUTE_N, CUTE_M, CUTE_NIQ, CUTE_NEQ
      integer          CUTE_IIQ(CUTE_MMAX)
      double precision CUTE_CEQ(CUTE_MMAX)

      common /cute_stuff/ CUTE_CEQ, CUTE_IIQ, CUTE_N, CUTE_M,
     1                    CUTE_NIQ, CUTE_NEQ
      save   /cute_stuff/

      double precision lam(CUTE_MMAX)

      call DCOPY(CUTE_M, 0d0, 0, lam, 1)

      call EVAL_HESSLAG_V(0, N, X, M, lam, VIN, VOUT, DAT, IDAT)

      return
      end
