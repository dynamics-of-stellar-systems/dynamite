      subroutine dstrsol(n,l,ldiag,jptr,indr,r,task)
      character*60 task
      integer n
      integer jptr(n+1), indr(*)
      double precision l(*), ldiag(n), r(n)
c     **********
c
c     Subroutine dstrsol
c
c     This subroutine solves the triangular systems L*x = r or L'*x = r.
c
c     The subroutine statement is
c
c       subroutine dstrsol(n,l,ldiag,jptr,indr,r,task)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of L.
c         On exit n is unchanged.
c
c       l is a double precision array of dimension *.
c         On entry l must contain the nonzeros in the strict lower
c            triangular part of L in compressed column storage.
c         On exit l is unchanged.
c
c       ldiag is a double precision array of dimension n.
c         On entry ldiag must contain the diagonal elements of L.
c         On exit ldiag is unchanged.
c
c       jptr is an integer array of dimension n + 1.
c         On entry jptr must contain pointers to the columns of A.
c            The nonzeros in column j of A must be in positions
c            jptr(j), ... , jptr(j+1) - 1.
c         On exit jptr is unchanged.
c
c       indr is an integer array of dimension *.
c         On entry indr must contain row indices for the strict 
c            lower triangular part of L in compressed column storage.
c         On exit indr is unchanged.
c
c       r is a double precision array of dimension n.
c         On entry r must contain the vector r.
c         On exit r contains the solution vector x.
c
c       task is a character variable of length 60.
c         On entry 
c            task(1:1) = 'N' if we need to solve L*x = r
c            task(1:1) = 'T' if we need to solve L'*x = r
c         On exit task is unchanged.
c
c     MINPACK-2 Project. May 1998.
c     Argonne National Laboratory.
c
c     **********
      double precision zero
      parameter(zero=0.0d0)

      integer j,k
      double precision temp

c     Solve L*x =r and store the result in r.

      if (task(1:1) .eq. 'N') then

         do j = 1, n
            temp = r(j)/ldiag(j)
            do k = jptr(j), jptr(j+1) - 1
               r(indr(k)) = r(indr(k)) - l(k)*temp
            end do
            r(j) = temp
         end do

         return

      end if

c     Solve L'*x =r and store the result in r.

      if (task(1:1) .eq. 'T') then

         r(n) = r(n)/ldiag(n)
         do j = n - 1, 1, -1
            temp = zero
            do k = jptr(j), jptr(j+1) - 1
               temp = temp + l(k)*r(indr(k))
            end do
            r(j) = (r(j) - temp)/ldiag(j)
         end do

         return

      end if

      end 
