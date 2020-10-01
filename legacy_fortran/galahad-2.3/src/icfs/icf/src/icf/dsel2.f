      subroutine dsel2(n,x,keys,k)
      integer n, k
      integer keys(n)
      double precision x(n)
c     **********
c
c     Subroutine dsel2
c
c     Given an array x of length n, this subroutine permutes
c     the elements of the array keys so that 
c    
c       abs(x(keys(i))) <= abs(x(keys(k))),  1 <= i <= k,
c       abs(x(keys(k))) <= abs(x(keys(i))),  k <= i <= n.
c
c     In other words, the smallest k elements of x in absolute value are
c     x(keys(i)), i = 1,...,k, and x(keys(k)) is the kth smallest element.
c     
c     The subroutine statement is
c
c       subroutine dsel2(n,x,keys,k)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of keys.
c         On exit n is unchanged.
c
c       x is a double precision array of length n.
c         On entry x is the array to be sorted.
c         On exit x is unchanged.
c
c       keys is an integer array of length n.
c         On entry keys is the array of indices for x.
c         On exit keys is permuted so that the smallest k elements
c            of x in absolute value are x(keys(i)), i = 1,...,k, and
c            x(keys(k)) is the kth smallest element.
c
c       k is an integer.
c         On entry k specifes the kth largest element.
c         On exit k is unchanged.
c
c     MINPACK-2 Project. March 1998.
c     Argonne National Laboratory.
c     William D. Kastak, Chih-Jen Lin, and Jorge J. More'.
c
c     **********
      integer i, l, lc, lp, m, p, p1, p2, p3, u
      integer swap

      if (n .le. 1 .or. k .le. 0 .or. k .gt. n) return

      u = n
      l = 1
      lc = n
      lp = 2*n

c     Start of iteration loop.

      do while (l .lt. u)

c        Choose the partition as the median of the elements in
c        positions l+s*(u-l) for s = 0, 0.25, 0.5, 0.75, 1.
c        Move the partition element into position l.

         p1 = (u+3*l)/4
         p2 = (u+l)/2
         p3 = (3*u+l)/4

c        Order the elements in positions l and p1.

         if (dabs(x(keys(l))) .gt. dabs(x(keys(p1)))) then
            swap = keys(l)
            keys(l) = keys(p1)
            keys(p1) = swap
            end if

c        Order the elements in positions p2 and p3.

         if (dabs(x(keys(p2))) .gt. dabs(x(keys(p3)))) then
            swap = keys(p2)
            keys(p2) = keys(p3)
            keys(p3) = swap
            end if

c        Swap the larger of the elements in positions p1
c        and p3, with the element in position u, and reorder
c        the first two pairs of elements as necessary.

         if (dabs(x(keys(p3))) .gt. dabs(x(keys(p1)))) then
            swap = keys(p3)
            keys(p3) = keys(u)
            keys(u) = swap
            if (dabs(x(keys(p2))) .gt. dabs(x(keys(p3)))) then
               swap = keys(p2)
               keys(p2) = keys(p3)
               keys(p3) = swap
               end if
         else
            swap = keys(p1)
            keys(p1) = keys(u)
            keys(u) = swap
            if (dabs(x(keys(l))) .gt. dabs(x(keys(p1)))) then
               swap = keys(l)
               keys(l) = keys(p1)
               keys(p1) = swap
               end if
            end if

c        If we define a(i) = abs(x(keys(i)) for i = 1,..., n, we have 
c        permuted keys so that 
c
c          a(l) <= a(p1), a(p2) <= a(p3), max(a(p1),a(p3)) <= a(u).
c
c        Find the third largest element of the four remaining
c        elements (the median), and place in position l.

         if (dabs(x(keys(p1))) .gt. dabs(x(keys(p3)))) then
            if (dabs(x(keys(l))) .le. dabs(x(keys(p3)))) then
               swap = keys(l)
               keys(l) = keys(p3)
               keys(p3) = swap
               end if
         else
            if (dabs(x(keys(p2))) .le. dabs(x(keys(p1)))) then
               swap = keys(l)
               keys(l) = keys(p1)
               keys(p1) = swap
            else
               swap = keys(l)
               keys(l) = keys(p2)
               keys(p2) = swap
               end if
            end if

c        Partition the array about the element in position l.

         m = l
         do i = l+1, u
            if (dabs(x(keys(i))) .lt. dabs(x(keys(l))))then
               m = m + 1
               swap = keys(m)
               keys(m) = keys(i)
               keys(i) = swap
               end if
         end do

c        Move the partition element into position m.

         swap = keys(l)
         keys(l) = keys(m)
         keys(m) = swap

c        Adjust the values of l and u.

         if (k .ge. m) l = m + 1
         if (k .le. m) u = m - 1

c        Check for multiple medians if the length of the subarray
c        has not decreased by 1/3 after two consecutive iterations.

         if (3*(u-l) .gt. 2*lp .and. k .gt. m) then

c           Partition the remaining elements into those elements
c           equal to x(m), and those greater than x(m). Adjust
c           the values of l and u.

            p = m
            do i = m+1, u
               if (dabs(x(keys(i))) .eq. dabs(x(keys(m)))) then
                  p = p + 1
                  swap = keys(p)
                  keys(p) = keys(i)
                  keys(i) = swap
                  end if
            end do
            l = p + 1
            if (k .le. p) u = p - 1
            end if

c        Update the length indicators for the subarray.

         lp = lc
         lc = u-l

      end do

      return

      end

