C**********************************************************
C                                                         *
C     Copyright (c) 2001-2004 by Ziena Optimization, Inc. *
C                   2001-2004 by Northwestern University  *
C     All Rights Reserved                                 *
C                                                         *
C**********************************************************
C
C     ----------------------------------------------------------------
C                   KNITRO DRIVER FOR CUTEr INTERFACE
C     ----------------------------------------------------------------
C
      program          driver
C
C     This driver uses reverse communication, and interfaces with 
C     CUTEr by calling the following subroutines: 
C
C       csetup ---- get the problem dimension n, the number m of 
C                   constraints, the bounds cl, cu, bl, bu. 
C 
C       cnames ---- get the problem name, problem_name.
C
C       cfn    ---- compute the objective function value f and the 
C                   constraints c.
C
C       csgr   ---- compute the gradient g of the objective function, 
C                   and the constraint gradients. 
C
C       csh    ---- compute the Hessian of the Lagrangian function.
C 
      implicit none

C**********************************************************************
C
C     We begin by defining some parameters.
C
C     These parameters can be made smaller or larger to adjust for
C     application problem and host machine memory sizes. Please see
C     the recommended sizes below.
C
      integer  max_n, max_m, max_nnzj, max_nnzh

C**********************************************************************
C
C     Recommended maximum array sizes  
C
C        n: number of variables
C        m: number of constraints
C
C     *  very small scale problems   n/m <=     200
C     *  small scale problems        n/m <=   1 000
C     *  medium scale problems       n/m <=  10 000
C     *  large scale problems        n/m <=  50 000
C     *  very large scale problems   n/m <= 250 000
C
C**********************************************************************
C
C     Very small scale problems
C
C      parameter (max_n       = 200)
C      parameter (max_m       = 200)
C      parameter (max_nnzj    = 40000)
C      parameter (max_nnzh    = 40000)

C     Small scale problems
C
C      parameter (max_n       = 1000)
C      parameter (max_m       = 1000)
C      parameter (max_nnzj    = 100000)
C      parameter (max_nnzh    = 100000)

C     Medium scale problems
C
C      parameter (max_n       = 10000)
C      parameter (max_m       = 10000)
C      parameter (max_nnzj    = 200000)
C      parameter (max_nnzh    = 200000)

C     Large scale problems
C
C      parameter (max_n       = 100000)
C      parameter (max_m       = 100000)
C      parameter (max_nnzj    = 2500000)
C      parameter (max_nnzh    = 2500000)

C     Very large scale problems
C
      parameter (max_n       = 250000)
      parameter (max_m       = 250000)
      parameter (max_nnzj    = 5000000)
      parameter (max_nnzh    = 5000000)
C
C**********************************************************************
C
C     Define some parameters used for reverse communication.
C
      integer           KTR_RC_INITIAL, KTR_RC_EVALFC, KTR_RC_EVALGA,
     $                  KTR_RC_EVALH, KTR_RC_EVALX0, KTR_RC_NEWPOINT
      parameter        (KTR_RC_INITIAL      =  0,
     $                  KTR_RC_EVALFC       =  1,
     $                  KTR_RC_EVALGA       =  2,
     $                  KTR_RC_EVALH        =  3,
     $                  KTR_RC_EVALX0       =  4,
     $                  KTR_RC_NEWPOINT     =  6)
C
C     Declare some variables, whose detailed description is at 
C     the end of this file.
C
      logical           equatn(max_m), ctype(max_m)
      integer           n, m, nnzh, nnzj, hcol(max_nnzh),
     $                  hrow(max_nnzh), indvar(max_nnzj),
     $                  indfun(max_nnzj), status, ftype
      real              CALLS( 7 ), CPU( 2 )
      character*20      problem_name, cdummy(max_n)
      double precision  f, x(max_n), c(max_m), cl(max_m), cu(max_m),
     $                  hess(max_nnzh), lambda(max_m+max_n),
     $                  ddummy(max_m), bl(max_n), bu(max_n),
     $                  fgrad(max_n), cjac(max_nnzj), vector(max_n)
C
C     Other local variables.
C     
      integer           i, j, ind, nnzg
      integer           gradopt, hessopt, errorflag
C
C**********************************************************************
C
C     Specify whether using exact Hessians or gradients.
C
      hessopt  = 1
      gradopt  = 1
C
C     Set up the problem. Get n, m, bl, bu, cl, cu.
C
      open(15, file='OUTSDIF.d', form='formatted')
      call csetup(15, 6, n, m, x, bl, bu, max_n, equatn, ctype,
     $            ddummy, cl, cu, max_m, .TRUE., .FALSE., .FALSE.)
      ftype = 0;
      if (n .gt. max_n) then
         print*,'ERROR: n > max_n; max_n must be at least ', n
         stop
      endif
      if (m .gt. max_m) then
         print*,'ERROR: m > max_m; max_m must be at least ', m
         stop
      endif      
C
C     Get the problem name.
C
      call cnames(n, m, problem_name, cdummy, cdummy)
C
C     We need to call the CUTEr routines 'cdimsj' and 'cdimsh'
C     before calling 'KTRsolveF' for the first time to get the
C     values nnzj and nnzh which are required on input.
C     (If we are not using the exact Hessian version of KNITRO 
C     then we do not call 'cdimsh' and nnzh is defined below.)
C  
      CALL CDIMSJ( nnzj )
      if (nnzj .gt. max_nnzj) then
         print*,'ERROR: nnzj > max_nnzj;',
     $        ' max_nnzj must be at least ', nnzj
         stop
      endif
      if (hessopt .eq. 1) then
C  
C     Exact Hessian version.
C
         CALL CDIMSH( nnzh )
      else
         nnzh = 0
      endif
      if (nnzh .gt. max_nnzh) then
         print*,'ERROR: nnzh > max_nnzh;',
     $        ' max_nnzh must be at least ', nnzh
         stop
      endif

      status = KTR_RC_INITIAL 
C
C             << beginning of loop >>
C
 111  continue
C
C     Compute function values.
C
      if ( (status .eq. KTR_RC_EVALFC) .or.
     $     (status .eq. KTR_RC_EVALX0) ) then
         call cfn(n, m, x, f, max_m, c)
      endif
C
C     Compute gradients.
C
      if ( (status .eq. KTR_RC_EVALGA) .or.
     $     (status .eq. KTR_RC_EVALX0) ) then
         call csgr(n, m, .FALSE., 1, ddummy, x, nnzj, max_nnzj,
     $        cjac, indvar, indfun)
C        Extract dense gradient from cjac
         do i=1, n
            fgrad(i) = 0.0d0
         end do
         ind = 1
         nnzg = 0
         do i=1, nnzj
            if (indfun(i) .eq. 0) then
               j = indvar(i)
               fgrad(j) = cjac(i)
               nnzg = nnzg + 1
            else 
C              This is OK since ind <= i always
               indfun(ind) = indfun(i) - 1
               indvar(ind) = indvar(i) - 1
               cjac(ind) = cjac(i)
               ind = ind + 1
            end if
         end do
         nnzj = nnzj - nnzg
      endif
C
C     Compute the Hessian hess of the Lagrangian function.
C
      if (status .eq. KTR_RC_EVALH) then
         if (hessopt .eq. 1) then
            call csh(n, m, x, m, lambda, nnzh, max_nnzh, 
     $           hess, hrow, hcol)
C           Decrement Hessian indices
            do i=1, nnzh 
               hrow(i) = hrow(i) - 1
               hcol(i) = hcol(i) - 1
            end do
         elseif (hessopt .eq. 5) then
            print*,'ERROR: hessopt=5 not operational in CUTEr interface'
            stop            
         endif
      endif
C
C     We now have a new solution estimate.  If callback feature 
C     enabled, provide callback routine(s) below.
C     
      if (status .eq. KTR_RC_NEWPOINT) then
C        The user may insert newpoint routines here if desired. 
      endif
C
C     Call the KNITRO solver.
C
      call KTRsolveF(f, ftype, n, x, bl, bu, fgrad, m, c, cl, cu,
     $               ctype, nnzj, cjac, indvar, indfun, lambda,
     $               nnzh, hess, hrow, hcol, vector, status,
     $               gradopt, hessopt)
C
C     If the stopping test has not been satisfied, re-enter the loop.
C
      if (status .gt. 0) goto 111
C
C                << end of the loop>>
C
      CALL CREPRT( CALLS, CPU )
      WRITE ( *, 2000 ) problem_name, n, m, CALLS(1),
     *     CALLS(2), CALLS(5), CALLS(6), status, f, CPU(1), CPU(2)
C
 2000 FORMAT( /, 24('*'), ' CUTEr statistics ', 24('*') //
     *     ,' Code used               :  KNITRO',    /
     *     ,' Problem                 :  ', A10,    /
     *     ,' # variables             =      ', I10 /
     *     ,' # constraints           =      ', I10 /
     *     ,' # objective functions   =      ', E15.7 /
     *     ,' # objective gradients   =      ', E15.7 / 
     *     ,' # constraints functions =      ', E15.7 /
     *     ,' # constraints gradients =      ', E15.7 /
     *     ' Exit code               =      ', I10 /
     *     ,' Final f                 = ', E15.7 /
     *     ,' Set up time             =      ', 0P, F10.2, ' seconds' /
     *     ' Solve time              =      ', 0P, F10.2, ' seconds' //
     *     66('*') / )

      stop
      end
C
C==================== The end of CUTEr driver =======================
