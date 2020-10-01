      subroutine ihsort(n,keys)
      integer n
      integer keys(n)
c     **********
c
c     Subroutine ihsort
c
c     Given an integer array keys of length n, this subroutine uses
c     a heap sort to sort the keys in increasing order.
c
c     This subroutine is a minor modification of code written by
c     Mark Jones and Paul Plassmann.
c
c     The subroutine statement is
c
c       subroutine ihsort(n,keys)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of keys.
c         On exit n is unchanged.
c
c       keys is an integer array of length n.
c         On entry keys is the array to be sorted.
c         On exit keys is permuted to increasing order.
c
c     MINPACK-2 Project. March 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      integer k, m, lheap, rheap, mid
      integer x
     
      if (n. le. 1) return
     
c     Build the heap.

      mid = n/2
      do k = mid, 1, -1
         x = keys(k)
         lheap = k
         rheap = n
         m = lheap*2
         do while (m .le. rheap)
            if (m .lt. rheap) then
               if (keys(m) .lt. keys(m+1)) m = m + 1
            endif
            if (x .ge. keys(m)) then
               m = rheap + 1
            else
               keys(lheap) = keys(m)
               lheap = m
               m = 2*lheap
            end if
         end do
      keys(lheap) = x
      end do
     
c     Sort the heap.

      do k = n, 2, -1
         x = keys(k)
         keys(k) = keys(1)
         lheap = 1
         rheap = k-1
         m = 2
         do while (m .le. rheap)
            if (m .lt. rheap) then
               if (keys(m) .lt. keys(m+1)) m = m+1
            endif
            if (x .ge. keys(m)) then
               m = rheap + 1
            else
               keys(lheap) = keys(m)
               lheap = m
               m = 2*lheap
            end if
         end do
      keys(lheap) = x
      end do
     
      return

      end
