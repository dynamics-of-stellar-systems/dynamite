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

