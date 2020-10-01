      subroutine insort(n,keys)
      integer n
      integer keys(n)
c     **********
c
c     Subroutine insort
c
c     Given an integer array keys of length n, this subroutine uses
c     an insertion sort to sort the keys in increasing order.
c
c     The subroutine statement is
c
c       subroutine insort(n,keys)
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
      integer i, j, ind

      do j = 2, n
         ind = keys(j)
         i = j - 1
         do while (i .gt. 0 .and. keys(i) .gt. ind)
            keys(i+1) = keys(i)
            i = i - 1
         end do
         keys(i+1) = ind
      end do

      return

      end 
