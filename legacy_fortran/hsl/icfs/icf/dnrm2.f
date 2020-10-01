      double precision function dnrm2(n,x,incx)
      integer n, incx
      double precision x(n)
c     **********
c
c     Function dnrm2
c
c     Given a vector x of length n, this function calculates the
c     Euclidean norm of x with stride incx.
c
c     The function statement is
c
c       double precision function dnrm2(n,x,incx)
c
c     where
c
c       n is an integer variable.
c         On entry n is the length of x.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x specifies the vector x.
c         On exit x is unchanged.
c
c       incx is an integer variable.
c         On entry incx specifies the stride.
c         On exit incx is unchanged.
c
c     MINPACK-2 Project. October 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick and Jorge J. More'.
c
c     **********
      double precision zero
      parameter (zero=0.0d0)

      integer i
      double precision scale

      scale = zero
      do 10 i = 1, n, incx
         scale = max(scale,abs(x(i)))
   10 continue

      dnrm2 = zero
      if (scale .eq. zero) return

      do 20 i = 1, n, incx
         dnrm2 = dnrm2 + (x(i)/scale)**2
   20 continue
      dnrm2 = scale*sqrt(dnrm2)

      end
