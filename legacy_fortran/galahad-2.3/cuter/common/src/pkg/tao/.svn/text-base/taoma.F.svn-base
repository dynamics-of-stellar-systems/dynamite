c Driver program for solving a .SIF example problem with TAO 
c (after sifdecode has been used on the .SIF file)

      program cutedriver
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/tao_solver.h"

c *** for statistical information
      integer n,nfev,ngrad,nhes
      common /stats/ n,nfev,ngrad,nhes

      double precision h_array(1800000)
      integer rowindex(1800000), colindex(1800000)
      common /hessian/ h_array, rowindex, colindex
      
      TAO_SOLVER       tao
      TAO_APPLICATION  myapp
      Vec              x,g,xl,xu
      Mat              H
      integer          infile, outfile
      PetscTruth       flg
      integer          nfree, nbind
      double precision f,gnorm
      double precision x_array(1), xl_array(1), xu_array(1),g_array(1)
      integer          x_index, xl_index, xu_index, g_index

      integer          i,info,nnzh
      double precision gatol,grtol
      external FormFunctionGradient, FormHessian, getstats

c *** absolute and relative gradient tolerances
      gatol = 1.0d-5  
      grtol = 1.0d-10

      call PetscInitialize(PETSC_NULL_CHARACTER,info)
      call TaoInitialize(PETSC_NULL_CHARACTER,info)

      infile = 1
      outfile = 6 
c *** outfile = 6  -> standard output      
      open (infile,file='OUTSDIF.d',status='old')


c *** Get the dimension size and create the vectors
      call udimen(infile,n)
      call VecCreateSeq(PETSC_COMM_SELF,n,x,info)
      call VecDuplicate(x,g,info)
      call VecDuplicate(x,xl,info)
      call VecDuplicate(x,xu,info)

c *** Generate the initial point and project into [xl,xu]

c *** When getting a petsc vector array using fortran, the VecGetArray
c     function returns a double precision address x_array and an offset 
c     where the data begins (x_index).  So to get the nth data element 
c     of x, use x_array(x_index+n) (using 1-based arrays)
      call VecGetArray(x,x_array,x_index,info)
      call VecGetArray(xl,xl_array,xl_index,info)
      call VecGetArray(xu,xu_array,xu_index,info)
      call VecGetArray(g,g_array,g_index,info)

      call usetup(infile, outfile, n, x_array(x_index+1),                 &
     &            xl_array(xl_index+1), xu_array(xu_index+1), n)
      
      do i=1,n
         if (xl_array(xl_index+i) .gt. xu_array(xu_index+i)) then
            print *,'Bad vector bounds'
            stop
         
         elseif (x_array(x_index+i) .lt. xl_array(xl_index+i))  then
            x_array(x_index+i) = xl_array(xl_index+i)
         elseif (x_array(x_index+i) .gt. xu_array(xu_index+i)) then
            x_array(x_index+i) = xu_array(x_index+i)
         endif
      enddo


c *** Count the free and binding variables, gradient projection norm
      call uofg(n,x_array(x_index+1), f, g_array(g_index+1),.true.)

      call getstats(n,x_array(x_index+1),xl_array(xl_index+1),              &
     &             xu_array(xu_index+1),g_array(g_index+1),nfree,nbind,     &
     &             gnorm)


      call VecRestoreArray(g,g_array,g_index,info)
      call VecRestoreArray(x,x_array,x_index,info)
      call VecRestoreArray(xl,xl_array,xl_index,info)
      call VecRestoreArray(xu,xu_array,xu_index,info)


c *** Get the number of nonzeros and create the hessian
      call udimsh(nnzh)
      call MatCreateSeqAIJ(PETSC_COMM_SELF,n,n,2*nnzh/n,                &
     &                     PETSC_NULL_INTEGER,H,info)
      call MatSetOption(H,MAT_SYMMETRIC,info)


      write (outfile,1000) n, nfree, n-nfree, nbind

      call TaoCreate(PETSC_COMM_SELF,'tao_blmvm',tao,info)
      call TaoPetscApplicationCreate(PETSC_COMM_SELF,myapp,info)

      call TaoSetTolerances(tao,0.0, 0.0, 0.0, 0.0, info)
      call TaoSetGradientTolerances(tao,gatol, grtol, 0.0, info)
      call TaoSetMaximumFunctionEvaluations(tao,10000,info)
      call TaoSetMaximumIterates(tao,10000,info)

      
      call TaoSetPetscFunctionGradient(myapp,x,g,FormFunctionGradient,  &
     &                                 PETSC_NULL_OBJECT,info)

      call TaoSetPetscHessian(myapp,H,H,FormHessian,nnzh,               &
     & info)

      call TaoSetPetscVariableBounds(myapp,xl,xu,info)
      call TaoSetPetscInitialVector(myapp,x,info)


      nfev = 0
      ngrad = 0
      nhes = 0
      call TaoSetApplication(tao,myapp,info)
      call TaoSetFromOptions(tao,info)
      call TaoSolve(tao,info)

c *** Count the free and binding variables again

c *** When getting a petsc vector array using fortran, the VecGetArray
c     function returns a double precision address x_array and an offset 
c     where the data begins (x_index).  So to get the nth data element 
c     of x, use x_array(x_index+n) (using 1-based arrays)

      call VecGetArray(x,x_array,x_index,info)
      call VecGetArray(xl,xl_array,xl_index,info)
      call VecGetArray(xu,xu_array,xu_index,info)
      call VecGetArray(g,g_array,g_index,info)

      call uofg(n,x_array(x_index+1), f, g_array(g_index+1),.true.)

      call getstats(n,x_array(x_index+1),xl_array(xl_index+1),              &
     &             xu_array(xu_index+1),g_array(g_index+1),nfree,nbind,     &
     &             gnorm)

      call VecRestoreArray(x,x_array,x_index,info)
      call VecRestoreArray(xl,xl_array,xl_index,info)
      call VecRestoreArray(xu,xu_array,xu_index,info)
      call VecRestoreArray(g,g_array,g_index,info)


      write (outfile,2000) nfree, n-nfree, nbind, nfev, f, gnorm


      call TaoView(tao)


      if (gnorm .le. gatol .or. gnorm .le. grtol*abs(f)) then
         write (outfile,4000)
      else
         write (outfile,5000) gnorm, gatol
      endif
         

      call VecDestroy(x,info)
      call VecDestroy(g,info)
      call VecDestroy(xl,info)
      call VecDestroy(xu,info)
      call MatDestroy(H,info)

      call TaoDestroy(tao,info)
      call TaoApplicationDestroy(myapp,info)

      close(outfile)
      close(infile)

      call TaoFinalize(info)
      call PetscFinalize(info)

 1000 format (
     +        ' Number of variables                         ',i12,/,
     +        ' Number of free variables at x(start)        ',i12,/,
     +        ' Number of active variables at x(start)      ',i12,/,
     +        ' Number of binding variables at x(start)     ',i12,/)
 2000 format (
     +        ' Number of free variables at x(final)        ',i12,/,
     +        ' Number of active variables at x(final)      ',i12,/,
     +        ' Number of binding variables at x(final)     ',i12,/,
     +        ' Number of function evaluations              ',i12,/,
     +        ' Function value at final iterate          '   ,d15.8,/,
     +        ' Projected gradient at final iterate      '   ,d15.3,/)


 4000 format (' Exit message     CONVERGENCE: GTOL TEST SATISFIED',//)
 5000 format (' Exit message     GRADIENT NORM ',f9.6,' > GTOL ',f9.6//)


      end


c ***************************************************** 
c ***  FormHessian
c ***       Calculates the hessian matrix at the
c ***       given vector.
c ***  
c ***       This function is called from withing the tao solver
c *** 
c *** IN:  tao - pointer to tao solver structure (used in C)
c ***      X   - the pointer to the PETSc vector 
c *** 
c *** OUT: H - PETSc matrix holding hessian matrix
c ***      Hpre - not used (for preconditioning)
c ***      flag - not used
c ***      nnzh - max number of nonzeros in the hessian
c ***      info - 0 if everything is ok
c *****************************************************

      subroutine FormHessian(tao, X, H, Hpre, flag, nnzh, info)
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/tao_solver.h"


      TAO_SOLVER tao
      Vec        X
      Mat        H,Hpre
      MatStructure flag
      integer    nnzh, info
      

      double precision h_array(1800000)
      integer colindex(1800000),rowindex(1800000)
      common /hessian/ h_array, colindex, rowindex

      integer n,nfev,ngrad,nhes
      common /stats/ n,nfev,ngrad,nhes

      external ush

      double precision x_array(1)
      integer x_index
      integer i

      
      info = 0

      nhes = nhes + 1

c *** When getting a petsc vector array using fortran, the VecGetArray
c     function returns a double precision address x_array and an offset 
c     where the data begins (x_index).  So to get the nth data element 
c     of x, use x_array(x_index+n) (using 1-based arrays)

      call VecGetArray(X,x_array,x_index,info)

c *** Get sparse hessian at x
      call ush(n,x_array(x_index+1),nnzh,nnzh,h_array,rowindex,colindex)

      call VecRestoreArray(X,x_array,x_index,info)

c *** Conver the sparse hessian to the petsc matrix
      call MatZeroEntries(H)
      do i=1,nnzh
         call MatSetValue(H,rowindex(i)-1,colindex(i)-1,h_array(i),         &
     &                    INSERT_VALUES, info)

         if (rowindex(i) .ne. colindex(i)) then
            call MatSetValue(H,colindex(i)-1,rowindex(i)-1,h_array(i),        &
     &                    INSERT_VALUES, info)
         endif
      enddo

      call MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY,info)
      call MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY,info)

      return
      end

c ***************************************************** 
c ***  FormFunctionGradient
c ***       Calculates the function and gradient at the
c ***       given vector.
c ***  
c ***       This function is called from withing the tao solver
c *** 
c *** IN:  tao - pointer to tao solver structure (used in C)
c ***      X   - the pointer to the PETSc vector 
c ***      dummy - not used
c *** 
c *** OUT: f - value of the function at X
c ***      G - PETSc vector holding value of gradient at X
c ***      info - 0 if everything is ok
c *****************************************************
      
      subroutine FormFunctionGradient(tao,X,f,G,dummy,info)
      implicit none

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscmat.h"
#include "include/finclude/tao_solver.h"


      TAO_SOLVER tao
      Vec X,G
      double precision f
      integer dummy, info
      

      integer n,nfev,ngrad,nhes
      common /stats/ n,nfev,ngrad,nhes


      external uofg

      double precision x_array(1),g_array(1)
      integer          x_index,g_index


      info = 0
      nfev = nfev + 1
      ngrad = ngrad + 1

c *** When getting a petsc vector array using fortran, the VecGetArray
c     function returns a double precision address x_array and an offset 
c     where the data begins (x_index).  So to get the nth data element 
c     of x, use x_array(x_index+n) (using 1-based arrays)

      call VecGetArray(X,x_array,x_index,info)
      call VecGetArray(G,g_array,g_index,info)

      call uofg(n,x_array(x_index+1), f, g_array(g_index+1),.true.)

      call VecRestoreArray(X,x_array,x_index,info)
      call VecRestoreArray(G,g_array,g_index,info)
      return
      end
      
      
c ***************************************************** 
c ***  getstats :
c ***       calculates the number of free and binding 
c ***       variables at the given vector
c *** IN:  n - number of elements in vector
c ***      x - array of vector elements
c ***      xl - array of lower bounds
c ***      xu - array of upper bounds
c ***      g - gradient array
c *** OUT: nfree - number of free variables at x
c ***      nbind - number of binding variables at x
c ***      gnorm - norm of the gradient vector at x
c *****************************************************
      subroutine getstats(n,x,xl,xu,g,nfree,nbind,gnorm)
      implicit none

      integer          n,nfree,nbind
      double precision x(n), xl(n), xu(n), g(n)
      double precision gnorm

      integer i

      nfree=0
      nbind =0

c *** count the free and binding variables      
      do i=1,n
         if (xl(i) .lt. x(i)  .and.  x(i) .lt. xu(i)) then
            nfree = nfree + 1
         elseif ((x(i).eq.xl(i) .and. g(i) .ge. 0)        .or.              &
     &           (x(i) .eq. xu(i) .and.   g(i) .le. 0)     .or.             &
     &           xl(i) .eq. xu(i)) then
            nbind = nbind + 1
         endif
      enddo

            
c *** calculate the gradient projection norm      
      gnorm = 0.0d0
      do i = 1, n
         if (xl(i) .ne. xu(i)) then
            if (x(i) .eq. xl(i)) then
               gnorm = gnorm + min(g(i),0.0d0)**2
            else if (x(i) .eq. xu(i)) then
               gnorm = gnorm + max(g(i),0.0d0)**2
            else
               gnorm = gnorm + g(i)**2
            end if
         end if
      end do
      gnorm = sqrt(gnorm)
 
      return 
      end

      
      
