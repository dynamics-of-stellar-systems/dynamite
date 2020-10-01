      subroutine srtdat2(n,nnz,a,adiag,row_ind,col_ind,col_ptr,iwa)
      integer n, nnz
      integer row_ind(*), col_ind(*), col_ptr(n+1)
      integer iwa(n)
      double precision a(*), adiag(n)
c     **********
c
c     subroutine srtdat2
c
c     Given the non-zero elements of an n by n matrix A in coordinate
c     format, this subroutine constructs the compressed column storage
c     with the diagonal elements stored separately.
c
c     On entry it is assumed that the elements are specified by
c
c           row_ind(k),col_ind(k), k = 1,...,nnz.
c
c     On exit row_ind is permuted so that the corresponding column
c     indices of col_ind are in non-decreasing order. In addition, the 
c     array col_ptr is set so that the row indices for the off-diagonal 
c     elements in column j are
c
c           row_ind(k), k = col_ptr(j),...,col_ptr(j+1)-1.
c
c     The subroutine statement is
c
c       srtdat2(n,nnz,a,adiag,row_ind,col_ind,col_ptr,iwa)
c
c     where
c
c       n is an integer variable.
c         On entry n is the order of the matrix.
c         On exit n is unchanged.
c
c       nnz is an integer variable.
c         On entry nnz is the number of non-zeros entries
c            in the input matrix. 
c         On exit nnz is unchanged.
c
c       a is a double precision array of dimension nnz.
c         On entry a contains the A in coordinate format.
c         On exit a contains A in compressed column storage.
c
c       adiag is a double precision array of dimension n.
c         On entry adiag need not be specified.
c         On exit adiag contains the diagonal elements of A.
c
c       row_ind is an integer array of length nnz. 
c         On entry row_ind must contain the row indices of the non-zero 
c            elements of A in coordinate format.
c         On exit row_ind is permuted so that the corresponding
c            column indices of col_ind are in non-decreasing order.
c
c       col_ind is an integer array of length nnz. 
c         On entry col_ind must contain the column indices of the 
c            non-zero elements of A. 
c         On exit col_ind is permuted so that these indices are in 
c            non-decreasing order.
c
c       col_ptr is an integer array of length n + 1. 
c         On entry col_ptr need not be specified.
c         On exit col_ptr specifies the locations of the row indices in 
c            the output row_ind. The row indices for the off-diagonal
c            elements in column j are
c
c              row_ind(k), k = col_ptr(j),...,col_ptr(j+1)-1.
c
c            Note that col_ptr(1) is set to 1 and that col_ptr(n+1)-1
c            is the number of non-zero off-diagonal elements.
c
c       iwa is an integer work array of length n.
c
c     MINPACK-2 Project. May 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      integer i, j, k, l, m
      double precision temp

      do j = 1, n
         iwa(j) = 0
      end do

c     Separate the diagonal elements and set m to the number
c     of non-zero off-diagonal elements. 

      m = 0
      do k = 1, nnz
         i = row_ind(k) 
         j = col_ind(k)
         if (i .eq. j) then
            adiag(i) = a(k)
         else
            m = m + 1
            row_ind(m) = i 
            col_ind(m) = j
            a(m) = a(k)

c           Store in array iwa the counts of non-zeroes in the columns.

            iwa(j) = iwa(j) + 1
         end if
      end do

c     Set pointers to the start of the columns in row_ind.

      col_ptr(1) = 1
      do  j = 1, n
         col_ptr(j+1) = col_ptr(j) + iwa(j) 
         iwa(j) = col_ptr(j)
      end do

c     Begin in-place sort.

      k = 1
      do while (k .le. m)
         j = col_ind(k)
        
         if (k .ge. col_ptr(j)) then
        
c        Current element is in position. Now examine the
c        next element or the first un-sorted element in
c        the j-th group.
        
            k = max(k+1,iwa(j))
         else
        
c        Current element is not in position. Place element
c        in position and make the displaced element the
c        current element.
        
            l = iwa(j)
            iwa(j) = iwa(j) + 1
            i = row_ind(k)
            temp = a(k)
            row_ind(k) = row_ind(l)
            col_ind(k) = col_ind(l)
            a(k) = a(l)
            row_ind(l) = i
            col_ind(l) = j
            a(l) = temp
         end if
      end do

      return

      end

