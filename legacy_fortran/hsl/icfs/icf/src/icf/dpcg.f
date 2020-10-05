      subroutine dpcg(n,a,adiag,acol_ptr,arow_ind,
     +               l,ldiag,lcol_ptr,lrow_ind,
     +               b,rtol,itmax,x,iters,info,p,q,r,z)
      integer n, itmax, iters, info
      integer acol_ptr(n+1), arow_ind(*), lcol_ptr(n+1), lrow_ind(*)
      double precision rtol
      double precision a(*), adiag(n), l(*), ldiag(n), b(n)
      double precision x(n)
      double precision p(n), q(n), r(n), z(n)
c     *********
c     
c     Subroutine dpcg
c     
c     Given a sparse symmetric matrix A in compressed row storage,
c     this subroutine uses a preconditioned conjugate gradiene method
c     to solve the system A*x = b.
c     
c     The subroutine statement is
c
c       subroutine dpcg(n,a,adiag,acol_ptr,arow_ind,
c                       l,ldiag,lcol_ptr,lrow_ind,
c                       b,rtol,itmax,x,iters,info,p,q,r,z)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of A.
c         On exit n is unchanged.
c
c       a is a double precision array of dimension *.
c         On entry a must contain the strict lower triangular part
c            of A in compressed column storage.
c         On exit a is unchanged.
c
c       adiag is a double precision array of dimension n.
c         On entry adiag must contain the diagonal elements of A.
c         On exit adiag is unchanged.
c
c       acol_ptr is an integer array of dimension n + 1.
c         On entry acol_ptr must contain pointers to the columns of A.
c            The nonzeros in column j of A must be in positions
c            acol_ptr(j), ... , acol_ptr(j+1) - 1.
c         On exit acol_ptr is unchanged.
c
c       arow_ind is an integer array of dimension *.
c         On entry arow_ind must contain row indices for the strict 
c            lower triangular part of A in compressed column storage.
c         On exit arow_ind is unchanged.
c
c       l is a double precision array of dimension nnz+n*p.
c         On entry l need not be specified.
c         On exit l contains the strict lower triangular part
c            of L in compressed column storage.
c
c       ldiag is a double precision array of dimension n.
c         On entry ldiag need not be specified.
c         On exit ldiag contains the diagonal elements of L.
c
c       lcol_ptr is an integer array of dimension n + 1.
c         On entry lcol_ptr need not be specified.
c         On exit lcol_ptr contains pointers to the columns of L.
c            The nonzeros in column j of L are in the
c            lcol_ptr(j), ... , lcol_ptr(j+1) - 1 positions of l.
c
c       lrow_ind is an integer array of dimension nnz+n*p.
c         On entry lrow_ind need not be specified.
c         On exit lrow_ind contains row indices for the strict lower
c            triangular part of L in compressed column storage.
c    
c       b is a double precision array of dimension n.
c         On entry b must contain the vector b.
c         On exit b is unchanged.
c
c       rtol is a double precision variable.
c         On entry rtol specifies the relative convergence test.
c         On exit rtol is unchanged
c     
c       itmax is an integer variable.
c         On entry itmax specifies the limit on the number of
c            conjugate gradient iterations.
c         On exit itmax is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x need not be specified.
c         On exit x contains the final conjugate gradient iterate.
c
c       iters is an integer variable.
c         On entry iters need not be specified.
c         On exit iters is set to the number of conjugate 
c            gradient iterations.
c     
c       info is an integer variable.
c         On entry info need not be specified.
c         On exit info is set as follows:
c
c             info = 1  The residual has the relative error
c                       specified by rtol.
c
c             info = 2  Negative curvature is encountered, and
c                       p is a direction of negative curvature.
c
c             info = 3  Failure to converge after itmax iterations.
c
c       p is a double precision work array of dimension n.
c     
c       q is a double precision work array of dimension n.
c     
c       r is a double precision work array of dimension n.
c     
c       z is a double precision work array of dimension n.
c     
c     Subprograms called
c
c       MINPACK-2  ......  dstrsol, dssyax
c
c       Level 1 BLAS  ...  daxpy, dcopy, ddot, dnrm2
c
c     MINPACK-2 Project. May 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      double precision zero, one
      parameter(zero=0.0d0,one=1.0d0)

      integer i
      double precision alpha, beta, bnorm, ptq, rho, rnorm, rtz

      external dstrsol, dssyax
      external daxpy, dcopy 
      double precision ddot, dnrm2

c     Initialize the iterate x and the residual r.

      iters = 0
      do i = 1, n
         x(i) = zero
      end do
      call dcopy(n,b,1,r,1)
      bnorm = dnrm2(n,b,1)
      rnorm = bnorm
      info = 1

c     Exit if b = 0.

      if (bnorm .eq. zero) return

c     Initialize z by solving L*L'*z = r.
         
      call dcopy(n,r,1,z,1)
      call dstrsol(n,l,ldiag,lcol_ptr,lrow_ind,z,'N')
      call dstrsol(n,l,ldiag,lcol_ptr,lrow_ind,z,'T')

c     Initialize p.

      call dcopy(n,z,1,p,1)

c     Initialize rho.

      rho = ddot(n,r,1,z,1) 

      do iters = 1, itmax

c        Compute q = A*p.

         call dssyax(n,a,adiag,acol_ptr,arow_ind,p,q)
         
c        Check for negative curvature.
         
         ptq = ddot(n,p,1,q,1)
         if (ptq .le. zero) then
            info = 2
            return
         end if

c        Update x.
         
         alpha = rho/ptq 
         call daxpy(n,alpha,p,1,x,1)

c        Update the residual.

         call daxpy(n,-alpha,q,1,r,1)
         rnorm = dnrm2(n,r,1)

c        Exit if the residual convergence test is satisfied.

         if (rnorm .le. rtol*bnorm) return
         
c        Solve L*L'*z = r.

         call dcopy(n,r,1,z,1)
         call dstrsol(n,l,ldiag,lcol_ptr,lrow_ind,z,'N')
         call dstrsol(n,l,ldiag,lcol_ptr,lrow_ind,z,'T')

c        Compute new rho.
         
         rtz = ddot(n,r,1,z,1) 

         beta = rtz/rho

c        Compute p = z + beta*p

         call dscal(n,beta,p,1)
         call daxpy(n,one,z,1,p,1)

c        Update rho.

         rho = rtz

      end do

      info = 3
      return

      end
