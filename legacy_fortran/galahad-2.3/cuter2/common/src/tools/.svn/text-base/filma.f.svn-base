C     ( Last modified on 23 Dec 2000 at 22:01:38 )
      PROGRAM FILMA
      
C     --------------------------------------------------------------
C     Solve NLP problems of the form
C     
C          minimize    f(x)
C          subject to  l_j <= c_j(x) <= u_j  ,    j = 1 , ... , m
C                      l_i <=   x_i  <= u_i  ,    i = 1 , ... , n
C
C     The problems solved are defined using CUTE.
C     --------------------------------------------------------------
C
C  CUTEr interface by Roger Fletcher and Sven Leyffer (U. Dundee)

C     ... internal constants
      INTEGER    NMAX,      MMAX,      MAXA,        KMX
CTOY  PARAMETER (NMAX=500, MMAX=500, MAXA=20000,  KMX=200)
CMED  PARAMETER (NMAX=1500, MMAX=1500, MAXA=40000,  KMX=500)
CBIG  PARAMETER (NMAX=5200, MMAX=5200, MAXA=500000,  KMX=100)
CCUS  PARAMETER (NMAX=2500, MMAX=2500, MAXA=100000,  KMX=200)
      INTEGER    MLP,     MXWK,         MXIWK,        MAXF
CTOY  PARAMETER (MLP=50, MXWK=1000000, MXIWK=500000, MAXF=50)
CMED  PARAMETER (MLP=50, MXWK=2000000, MXIWK=500000, MAXF=50)
CBIG  PARAMETER (MLP=200, MXWK=10000000, MXIWK=500000, MAXF=50)
CCUS  PARAMETER (MLP=100, MXWK=5000000, MXIWK=250000, MAXF=50)
C     ... space for user workspaces (real and integer) (here for CUTE)
      INTEGER    LUSER,      LIUSER
      PARAMETER (LUSER=NMAX, LIUSER=MAXA)
      
C     ... internal variables -- scalars
      INTEGER N, M, IPRINT, I, J, IDUMMY,
     .        M_NLN, IFAIL, KMAX, NOUT, MAX_ITER
CS    REAL                RHO, F, FMIN, 
CD    DOUBLE PRECISION    RHO, F, FMIN, 
     .        CPU_START, CPU_END, CPU_TIM, SECONDS

C     ... internal variables -- arrays
      INTEGER LA(0:MAXA+MMAX+2),ISTAT(14),IUSER(LIUSER),LWS(MXIWK)
CS    REAL                A(MAXA),BLO(NMAX+MMAX),BUP(NMAX+MMAX),X(NMAX),
CD    DOUBLE PRECISION    A(MAXA),BLO(NMAX+MMAX),BUP(NMAX+MMAX),X(NMAX),
     .        C(MMAX), WS(MXWK),RSTAT(7),USER(LUSER),LAM(NMAX+MMAX),
     .        S(NMAX+MMAX)
      LOGICAL EQUATN(MMAX), LINEAR(MMAX)
      CHARACTER*10 XNAMES(NMAX),GNAMES(MMAX)
      CHARACTER    CSTYPE(MMAX)

C     ... common statements
CS    REAL                             INFTY, EPS
CD    DOUBLE PRECISION                 INFTY, EPS
      COMMON /NLP_EPS_INF/ INFTY, EPS

C     ... common to indicate initial penalty parameter & updating or not
CS    REAL                           GIVEN_MU
CD    DOUBLE PRECISION               GIVEN_MU
      LOGICAL                      UPDATE_MU
      COMMON /PENALTY_C/ GIVEN_MU, UPDATE_MU

C     ... upper bound on filter
CS    REAL                      UBD, TT
CD    DOUBLE PRECISION          UBD, TT
      COMMON /UBDC/ UBD, TT

C     ... problem name (& length)
      INTEGER         CHAR_L
      CHARACTER*10            PNAME
      COMMON /CPNAME/ CHAR_L, PNAME

C     ... default options set here
      DATA IPRINT, MAX_ITER, NOUT, RHO, IDUMMY /1, 1000, 6, 1.E1, 0/

C     ========================  procedure body  ======================

C     ... read filter.par parameter file (or use defaults)
      CALL READPAR(IPRINT,IDUMMY,IDUMMY,MAX_ITER,IDUMMY,IDUMMY,IDUMMY,
     .             IDUMMY,NOUT,RHO,IDUMMY)

C     ... initialization
      CALL INITIALIZE_NLP(N,M,NMAX,MMAX,NOUT,BLO,BUP,X,LAM,EQUATN,
     .                    LINEAR,CSTYPE,XNAMES,GNAMES,USER,IUSER)
      IF ((N.GT.NMAX).OR.(M.GT.MMAX)) THEN
         WRITE(NOUT,*) 'Invalid choics for n, m'
         WRITE(NOUT,*) 'n,nmax,m,mmax=',n,nmax,m,mmax
         STOP
      ENDIF
      
C     ... variable/constraint scale factors
      CALL READSCALE (N,M,XNAMES,GNAMES,PNAME,CHAR_L,S,IFAIL)

C     ... initialize penalty parameter (for monotoring only)
      UPDATE_MU = .TRUE.
      GIVEN_MU  = 1.E0

C     ... set fmin, kmax
      FMIN  = - INFTY
      KMAX  = MIN ( N , KMX )
      IFAIL = 0

C     ... call the main SQP routine
      CPU_start = seconds()
      CALL FILTERSQP (N,M,KMAX,MAXA,MAXF,MLP,MXWK,MXIWK,IPRINT,NOUT,
     .                IFAIL,RHO,X,C,F,FMIN,BLO,BUP,S,A,LA,WS,LWS,
     .                LAM,CSTYPE,USER,IUSER,MAX_ITER,ISTAT,RSTAT)
      CPU_END  = SECONDS()
      CPU_TIM = CPU_END - CPU_START
      IF (IPRINT.GE.1) THEN
         WRITE(NOUT,*)' CPU time for this solve.............',CPU_tim
      ENDIF

C     ... count number of nonlinear c/s (for output)
      M_NLN = 0
      DO J=1,M
         IF (CSTYPE(J).EQ.'N') M_NLN = M_NLN + 1
      ENDDO
C     ... save result in a compact format
      OPEN (UNIT=12,FILE='00temp.txt')
      WRITE(12,9001) PNAME, IFAIL, N, M, M_NLN, ISTAT(1)
      WRITE(12,9001) PNAME, (ISTAT(J),J=2,7)
      WRITE(12,9001) PNAME, (ISTAT(J),J=8,14)
      WRITE(12,9002) PNAME, F,(RSTAT(J),J=5,7)
      WRITE(12,9002) PNAME, RHO,(RSTAT(J),J=1,4),CPU_TIM
      CLOSE(12)

C     ... write solution onto a files
      OPEN(UNIT=1, FILE=PNAME(1:CHAR_L)//'.solution')
      WRITE(1,*) 'Problemname.........',PNAME
      WRITE(1,*) 'No. of variables....',N
      WRITE(1,*) 'No. of constraints..',M
      WRITE(1,*) 'Solution:'
      WRITE(1,*) '========='
      WRITE(1,*) '  F* = ', F
      WRITE(1,*)
      WRITE(1,'(2A)')' Name    |    Lower bd   |       X*      |',
     .               '  Upper bd    | Multiplier   | Scales'
      WRITE(1,'(2A)')'---------+---------------+---------------+',
     .               '--------------+--------------+---------'
      WRITE(1,7002) (XNAMES(I),BLO(I),X(I),BUP(I),LAM(I),S(I),I=1,N)
      WRITE(1,*)
      WRITE(1,'(2A)')' Name    |    Lower BD   |      C(X*)    |',
     .               '  Upper BD    | Multiplier   | Scales'
      WRITE(1,'(2A)')'---------+---------------+---------------+',
     .               '--------------+--------------+---------'
      WRITE(1,7003)(GNAMES(I),BLO(N+I),C(I),BUP(N+I),LAM(N+I),
     .              S(N+I),CSTYPE(I),I=1,M)
      WRITE(1,*)
      WRITE(1,*) '------------------------------------------------',
     .     '-----------------------------'

 7000 FORMAT(2A)
 7001 FORMAT(4(I4,3X,A10,'|')) 
 7002 FORMAT(A,5G15.7)
 7003 FORMAT(A,5G15.7,A)
 9001 FORMAT(A,7(1X,I6))
 9002 FORMAT(A,6(G10.4,1X))

      STOP
      END

