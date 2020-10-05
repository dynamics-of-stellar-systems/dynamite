      subroutine dssyax(n,a,adiag,jptr,indr,x,y)
      integer n
      integer jptr(n+1), indr(*)
      double precision a(*), adiag(n), x(n), y(n)
c     **********
c
c     Subroutine dssyax
c
c     This subroutine computes the matrix-vector product y = A*x, 
c     where A is a symmetric matrix with the strict lower triangular 
c     part in compressed column storage.
c
c     The subroutine statement is
c
c       subroutine dssyax(n,a,adiag,jptr,indr,x,y)
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
c       jptr is an integer array of dimension n + 1.
c         On entry jptr must contain pointers to the columns of A.
c            The nonzeros in column j of A must be in positions
c            jptr(j), ... , jptr(j+1) - 1.
c         On exit jptr is unchanged.
c
c       indr is an integer array of dimension *.
c         On entry indr must contain row indices for the strict 
c            lower triangular part of A in compressed column storage.
c         On exit indr is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x must contain the vector x.
c         On exit x is unchanged.
c
c       y is a double precision array of dimension n.
c         On entry y need not be specified.
c         On exit y contains the product A*x.
c
c     MINPACK-2 Project. May 1998.
c     Argonne National Laboratory.
c
c     **********
      double precision zero
      parameter (zero=0.0d0)

      integer i, j
      double precision sum

      do i = 1, n
         y(i) = adiag(i)*x(i)
      end do

      do j = 1, n
         sum = zero
         do i = jptr(j), jptr(j+1) - 1
            sum = sum + a(i)*x(indr(i))
            y(indr(i)) = y(indr(i)) + a(i)*x(j)
         end do
         y(j) = y(j) + sum
      end do

      return

      end
