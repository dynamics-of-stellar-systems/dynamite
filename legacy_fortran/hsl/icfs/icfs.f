      subroutine dicfs(n,nnz,a,adiag,acol_ptr,arow_ind,
     +                 l,ldiag,lcol_ptr,lrow_ind,
     +                 p,alpha,iwa,wa1,wa2)
      integer n, nnz, p
      integer acol_ptr(n+1), arow_ind(nnz)
      integer lcol_ptr(n+1), lrow_ind(nnz+n*p)
      integer iwa(3*n)
      double precision alpha
      double precision wa1(n), wa2(n)
      double precision a(nnz), adiag(n), l(nnz+n*p), ldiag(n)
c     *********
c
c     Subroutine dicfs
c
c     Given a symmetric matrix A in compressed column storage, this
c     subroutine computes an incomplete Cholesky factor of A + alpha*D,
c     where alpha is a shift and D is the diagonal matrix with entries
c     set to the l2 norms of the columns of A.
c
c     The subroutine statement is
c
c       subroutine dicfs(n,nnz,a,adiag,acol_ptr,arow_ind,
c                        l,ldiag,lcol_ptr,lrow_ind,
c                        p,alpha,iwa,wa1,wa2)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of A.
c         On exit n is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz is the number of nonzeros in the strict lower
c            triangular part of A.
c         On exit nnz is unchanged.
c
c       a is a double precision array of dimension nnz.
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
c       arow_ind is an integer array of dimension nnz.
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
c       p is an integer variable.
c         On entry p specifes the amount of memory available for the
c            incomplete Cholesky factorization.
c         On exit p is unchanged.
c
c       alpha is a double precision variable.
c         On entry alpha is the initial guess of the shift.
c         On exit alpha is final shift
c
c       iwa is an integer work array of dimesnion 3*n.
c
c       wa1 is a double precision work array of dimension n.
c
c       wa2 is a double precision work array of dimension n.
c
c     Subprograms called
c
c       MINPACK-2  ......  dicf
c
c     MINPACK-2 Project. October 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      integer nbmax
      parameter(nbmax=3)
      double precision alpham, nbfactor
      parameter(alpham=1.0d-3,nbfactor=512)
      double precision zero, one, two
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0)

      integer i, info, j, k, nb
      double precision alphas

      external dicf
      
c     Compute the l2 norms of the columns of A.

      do i = 1, n
         wa1(i) = adiag(i)**2
      end do
      do j = 1, n
         do i = acol_ptr(j), acol_ptr(j+1)-1
            k = arow_ind(i)
            wa1(j) = wa1(j) + a(i)**2
            wa1(k) = wa1(k) + a(i)**2
         end do
      end do
      do j = 1, n
         wa1(j) = sqrt(wa1(j))
      end do

c     Compute the scaling matrix D.

      do i = 1, n
         if (wa1(i) .gt. zero) then
            wa2(i) = one/sqrt(wa1(i))
         else
	    wa2(i) = one
         endif
      end do

c     Determine a lower bound for the step. 

      if (alpha .le. zero) then
         alphas = alpham
      else
         alphas = alpha
      end if

c     Compute the initial shift.

      alpha = zero
      do i = 1, n
         if (adiag(i) .eq. zero) then
            alpha = alphas
         else
            alpha = max(alpha,-adiag(i)*(wa2(i)**2))
         end if
      end do
      if (alpha .gt. zero) alpha = max(alpha,alphas)

c     Search for an acceptable shift. During the search we decrease 
c     the lower bound alphas until we determine a lower bound that 
c     is not acceptable. We then increase the shift.
c     The lower bound is decreased by nbfactor at most nbmax times.

      nb = 1
      do while (1 .eq. 1)

c        Copy the sparsity structure of A into L.

         do i = 1, n+1
            lcol_ptr(i) = acol_ptr(i)
         end do
         do i = 1, nnz
            lrow_ind(i) = arow_ind(i)
         end do

c        Scale A and store in the lower triangular matrix L.

         do j = 1, n
            ldiag(j) = adiag(j)*(wa2(j)**2) + alpha
         end do
         do j = 1, n
            do i = acol_ptr(j), acol_ptr(j+1)-1
               l(i) = a(i)*wa2(j)*wa2(arow_ind(i))
            end do
         end do

c        Attempt an incomplete factorization.

         call dicf(n,nnz,l,ldiag,lcol_ptr,lrow_ind,p,info,
     +        iwa(1),iwa(n+1),iwa(2*n+1),wa1)

c        If the factorization exists, then test for termination.
c        Otherwise increment the shift.

         if (info .ge. 0) then

c           If the shift is at the lower bound, reduce the shift.
c           Otherwise undo the scaling of L and exit.

            if (alpha .eq. alphas .and. nb .lt. nbmax) then
               alphas = alphas/nbfactor
               alpha = alphas
               nb = nb + 1
            else
               do i = 1, n
                  ldiag(i) = ldiag(i)/wa2(i)
               end do
               do j = 1, lcol_ptr(n+1)-1
                  l(j) = l(j)/wa2(lrow_ind(j))
               end do
               return
            end if
         else
            alpha = max(two*alpha,alphas)
         end if
      end do

      return

      end
 


 
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
      subroutine dicf(n,nnz,a,diag,col_ptr,row_ind,p,info,
     +                indr,indf,list,w)
      integer n, nnz, p, info
      integer col_ptr(n+1), row_ind(*)
      integer indr(n), indf(n), list(n)
      double precision a(*), diag(n)
      double precision w(n)
c     *********
c     
c     Subroutine dicf
c     
c     Given a sparse symmetric matrix A in compressed row storage,
c     this subroutine computes an incomplete Cholesky factorization.
c
c     Implementation of dicf is based on the Jones-Plassmann code.
c     Arrays indf and list define the data structure.
c     At the beginning of the computation of the j-th column,
c
c       For k < j, indf(k) is the index of a for the first 
c       nonzero l(i,k) in the k-th column with i >= j.
c
c       For k < j, list(i) is a pointer to a linked list of column 
c       indices k with i = row_ind(indf(k)).
c
c     For the computation of the j-th column, the array indr records
c     the row indices. Hence, if nlj is the number of nonzeros in the 
c     j-th column, then indr(1),...,indr(nlj) are the row indices. 
c     Also, for i > j, indf(i) marks the row indices in the j-th  
c     column so that indf(i) = 1 if l(i,j) is not zero.
c
c     The subroutine statement is
c
c       subroutine dicf(n,nnz,a,diag,col_ptr,row_ind,p,info,
c                       indr,indf,list,w)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of A.
c         On exit n is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz is the number of nonzeros in the strict lower
c            triangular part of A.
c         On exit nnz is unchanged.
c
c       a is a double precision array of dimension nnz+n*p.
c         On entry the first nnz entries of a must contain the strict 
c            lower triangular part of A in compressed column storage.
c         On exit a contains the strict lower triangular part
c            of L in compressed column storage.
c
c       diag is a double precision array of dimension n.
c         On entry diag must contain the diagonal elements of A.
c         On exit diag contains the diagonal elements of L.
c
c       col_ptr is an integer array of dimension n + 1.
c         On entry col_ptr must contain pointers to the columns of A.
c            The nonzeros in column j of A must be in positions
c            col_ptr(j), ... , col_ptr(j+1) - 1.
c         On exit col_ptr contains pointers to the columns of L.
c            The nonzeros in column j of L are in the
c            col_ptr(j), ... , col_ptr(j+1) - 1 positions of l.
c
c       row_ind is an integer array of dimension nnz+n*p.
c         On entry row_ind must contain row indices for the strict 
c            lower triangular part of A in compressed column storage.
c         On exit row_ind contains row indices for the strict lower
c            triangular part of L in compressed column storage. 
c
c       p is an integer variable.
c         On entry p specifes the amount of memory available for the
c            incomplete Cholesky factorization.
c         On exit p is unchanged.
c
c       info is an integer variable.
c         On entry info need not be specified.
c         On exit info = 0 if the factorization succeeds, and
c            info < 0 if the -info pivot is not positive.
c
c       indr is an integer work array of dimension n.
c     
c       indf is an integer work array of dimension n.
c
c       list is an integer work array of dimension n.
c
c       w is a double precision work array of dimension n.
c
c     Subprograms called
c
c       MINPACK-2  ......  dsel2, ihsort, insort
c
c       Level 1 BLAS  ...  daxpy, dcopy, ddot, dnrm2
c
c     MINPACK-2 Project. May 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      integer insortf
      parameter (insortf=20)
      double precision zero
      parameter (zero=0.0d0)

      integer i, ip, j, k, kth, nlj, newk, np, mlj
      integer isj, iej, isk, iek, newisj, newiej
      double precision lval

      external dsel2, ihsort, insort

      info = 0
      do j = 1, n
         indf(j) = 0
         list(j) = 0
      end do

c     Make room for L by moving A to the last n*p positions in a.

      np = n*p
      do  j = 1, n + 1
         col_ptr(j) = col_ptr(j) + np
      end do
      do j = nnz, 1, -1
         row_ind(np+j) = row_ind(j)
         a(np+j) = a(j)
      end do

c     Compute the incomplete Cholesky factorization.

      isj = col_ptr(1)
      col_ptr(1) = 1
      do j = 1, n

c        Load column j into the array w. The first and last elements 
c        of the j-th column of A are a(isj) and a(iej).

         nlj = 0
         iej = col_ptr(j+1) - 1
         do ip = isj, iej
            i = row_ind(ip)
            w(i) = a(ip)
            nlj = nlj + 1
            indr(nlj) = i
            indf(i) = 1
         end do

c        Exit if the current pivot is not positive.

         if (diag(j) .le. zero) then
            info = -j
            return
         end if
         diag(j) = sqrt(diag(j))

c        Update column j using the previous columns.

         k = list(j)
         do while (k .ne. 0) 
            isk = indf(k)
            iek = col_ptr(k+1) - 1

c           Set lval to l(j,k).

            lval = a(isk)

c           Update indf and list.

            newk = list(k)
            isk = isk + 1
            if (isk .lt. iek) then
               indf(k) = isk
               list(k) = list(row_ind(isk))
               list(row_ind(isk)) = k
            endif
            k = newk

c           Compute the update a(i,i) <- a(i,j) - l(i,k)*l(j,k).
c           In this loop we pick up l(i,k) for k < j and i > j.

            do ip = isk, iek
               i = row_ind(ip)
               if (indf(i) .ne. 0) then
                  w(i) = w(i) - lval*a(ip)
               else
                  indf(i) = 1
                  nlj = nlj + 1
                  indr(nlj) = i
                  w(i) = - lval*a(ip)
               end if
            end do
         end do

c        Compute the j-th column of L.

         do k = 1, nlj
            w(indr(k)) = w(indr(k))/diag(j)
         end do

c        Set mlj to the number of nonzeros to be retained.

         mlj = min(iej-isj+1+p,nlj)
         kth = nlj - mlj + 1

        if (nlj .ge. 1) then

c           Determine the kth smallest elements in the current 
c           column, and hence, the largest mlj elements.
          
            call dsel2(nlj,w,indr,kth)
          
c           Sort the row indices of the selected elements. Insertion
c           sort is used for small arrays, and heap sort for larger
c           arrays. The sorting of the row indices is required so that
c           we can retrieve l(i,k) with i > k from indf(k).
          
            if (mlj .le. insortf) then
               call insort(mlj,indr(kth))
            else
               call ihsort(mlj,indr(kth))
            end if
         end if

c        Store the largest elements in L. The first and last elements 
c        of the j-th column of L are a(newisj) and a(newiej).

         newisj = col_ptr(j)
         newiej = newisj + mlj -1
         do k = newisj, newiej
            a(k) = w(indr(k-newisj+kth))
            row_ind(k) = indr(k-newisj+kth)
         end do

c        Update the diagonal elements.

         do k = kth, nlj
            diag(indr(k)) = diag(indr(k)) - w(indr(k))**2
         end do

c        Update indf and list for the j-th column.

         if (newisj .lt. newiej) then
            indf(j) = newisj
            list(j) = list(row_ind(newisj))
            list(row_ind(newisj)) = j
         endif

c        Clear out elements j+1,...,n of the array indf.

         do k = 1, nlj
            indf(indr(k)) = 0
         end do

c        Update isj and col_ptr.

         isj = col_ptr(j+1)
         col_ptr(j+1) = newiej + 1

      end do

      return

      end

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

