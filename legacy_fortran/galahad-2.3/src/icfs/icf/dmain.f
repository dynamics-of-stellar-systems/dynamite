      program dmain

      character*10 fread, fwrite
      parameter (fread='icf.in',fwrite='icf.out')
      integer nnzmax, nmax, maxnp
      parameter (nnzmax=600000,nmax=30000,maxnp=25)
      integer nread, nwrite, mmunit
      parameter (nread=1,nwrite=2,mmunit=8)
c     **********
c
c     Sample driver for the incomplete Cholesky decomposition.
c
c     MINPACK-2 Project. July 1998.
c     Argonne National Laboratory.
c     Chih-Jen Lin and Jorge J. More'.
c
c     **********
      character *60 prob_name

      integer i, j, np, nnza, nnzl
      integer n, nnz, info, iters, maxiter
      integer p(maxnp)
      double precision ftime, stime
      double precision fval, rerr

      integer arow_ind(nnzmax), acol_ind(nnzmax), acol_ptr(nmax+1)
      integer lcol_ptr(nmax), lrow_ind(nnzmax)
      integer iw(3*nmax)
      double precision alpha, rtol 
      double precision a(nnzmax), adiag(nmax), l(nnzmax), ldiag(nmax)
      double precision w(4*nmax) 
      double precision x(nmax), b(nmax)

c     Variables for the Matrix Market format.

      integer nrows, ncols
      character rep*10, field*7, symm*19
      integer ival
      complex cval

      double precision ddot, dnrm2
      external mminfo, mmread
      external dicfs, dnrm2, dpcg, dssyax, dtimer, srtdat2

      open (nread,file=fread,status='old')
      open (nwrite,file=fwrite)

c     Read values of p for the ICF preconditioner.

      read (nread,*) np
      read (nread,*) (p(i),i=1,np)

      do while (1 .eq. 1)

         read (nread,*) prob_name
         if (prob_name .eq. 'STOP') stop
         write (*,*) ' '
         write (*,*) prob_name
         open (mmunit,file=prob_name)

c        Read the Matrix Market information.

         call mminfo(mmunit,rep,field,symm,nrows,ncols,nnz)

c        Check the Matrix Market information.

         if (rep .ne. 'coordinate' .or. field .ne. 'real' .or.
     +      symm .ne. 'symmetric' .or. nrows .gt. nmax) then
            print *,'Matrix Market input error'
            stop
         end if

c        Read the matrix elements.

         call mmread(mmunit,rep,field,symm,nrows,ncols,nnz,
     +        nnzmax,arow_ind,acol_ind,ival,a,cval)
         n = nrows

c        Change from coordinate storage to compresses column storage.

         call srtdat2(n,nnz,a,adiag,arow_ind,acol_ind,acol_ptr,iw)

c        Generate the right hand side.
c        We set b = A*e, where e is the vector of ones.
c        Note that the paper of Lin and More', "Incomplete Cholesky 
c        factorizations with limited memory", uses b = e.

         do i = 1, n
            b(i) = adiag(i)
         end do
         do i = 1, n
            do k = acol_ptr(i), acol_ptr(i+1)-1
               b(arow_ind(k)) = b(arow_ind(k)) + a(k)
               b(i) = b(i) + a(k)
            end do
         end do

c        Modify the matrix for the optimization problems.

         if (prob_name(8:11) .eq. 'dgl2' .or. 
     +       prob_name(8:14) .eq. 'nlmsurf' .or. 
     +       prob_name(8:13) .eq. 'jimack') then
            do i = 1, n
               adiag(i) = adiag(i)*(1.0d0 + 1.d-5)
            end do
         end if

         do j = 1, np

c           Check for sufficient memory.

            nnz = acol_ptr(n+1) - 1 
            if (nnz + p(j)*n .gt. nnzmax) then
               write (*,*) 'Parameter nnzmax is too small'
               stop
            end if

c           Compute the preconditioner L.
            
            call dtimer(stime)

            if (p(j) .lt. 0) then
               do i = 1, n
                  ldiag(i) = sqrt(adiag(i))
                  lcol_ptr(i) = 1
               end do
               lcol_ptr(n+1) = 1
            else
               alpha = 0.0d0
               call dicfs(n,nnz,a,adiag,acol_ptr,arow_ind,
     +                    l,ldiag,lcol_ptr,lrow_ind,
     +                    p(j),alpha,iw,w,w(n+1))
            end if

c           Solve Ax = b with the preconditioned 
c           conjugate gradient method.

            rtol = 1.0d-6
            maxiter = max(n,100)
            call dpcg(n,a,adiag,acol_ptr,arow_ind,l,ldiag,
     +           lcol_ptr,lrow_ind,b,rtol,maxiter,x,
     +           iters,info,w,w(n+1),w(2*n+1),w(3*n+1))

            if (info .eq. 2) then
               write (*,*) 'Indefinite matrix'
               call dssyax(n,a,adiag,acol_ptr,arow_ind,w,w(n+1))
               write (*,*) 'Negative curvature',
     +            ddot(n,w,1,w(n+1),1)/dnrm2(n,w,1)**2
            end if

            call dtimer(ftime)

c           Compute relative error and function value.

            call dssyax(n,a,adiag,acol_ptr,arow_ind,x,w)

            fval = 0.5*ddot(n,x,1,w,1) - ddot(n,b,1,x,1)
            do i = 1, n
               w(i) = w(i) - b(i)
            end do
            rerr = dnrm2(n,w,1)/dnrm2(n,b,1)

            nnzl = lcol_ptr(n+1) - 1
            nnza = acol_ptr(n+1) - 1

            write (nwrite,1000) prob_name(8:30), n, p(j), nnz, nnzl, 
     +         dble(nnzl+n)/dble(nnz+n),
     +         dble(nnza+n*p(j)-nnzl)/dble(nnza+n*p(j)),
     +         iters,alpha,info,rerr,fval,ftime-stime
         end do

      end do

      close(nread)
      close(nwrite)
      close(mmunit)

      stop

 1000 format (1x,a23/,
     +        1x,'Order of the matrix                     ',3x,i10/,
     +        1x,'Memory parameter p                      ',3x,i10/,
     +        1x,'Number of nonzeroes in A                ',3x,i10/,
     +        1x,'Number of nonzeroes in L                ',3x,i10/,
     +        1x,'Memory usage in L nnz(L)/nnz(A)         ',3x,f10.2/,
     +        1x,'Wasted memory in L                      ',3x,f10.2/,
     +        1x,'Number of conjugate gradient iterations ',3x,i10/,
     +        1x,'Final shift                             ',3x,d10.3/,
     +        1x,'Info value for pcg                      ',3x,i10/,
     +        1x,'Relative error in the residual          ',3x,d10.3/,
     +        1x,'Function value                '          ,3x,d20.8/,
     +        1x,'Time for icf and pcg                    ',3x,d10.3/)

      end

