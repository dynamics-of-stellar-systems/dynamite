! *******************************************************************
! COPYRIGHT (c) 2001 Council for the Central Laboratory
!               of the Research Councils
! All rights reserved.
!
! None of the comments in this Copyright notice between the lines
! of asterisks shall be removed or altered in any way.
!
! This Package is intended for compilation without modification,
! so most of the embedded comments have been removed.
!
! ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
! SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
!
! Please note that for an ACADEMIC Licence:
!
! 1. The Packages may only be used for academic research or teaching
!    purposes by the Licensee, and must not be copied by the Licensee for
!    use by any other persons. Use of the Packages in any commercial
!    application shall be subject to prior written agreement between
!    Hyprotech UK Limited and the Licensee on suitable terms and
!    conditions, which will include financial conditions.
! 2. All information on the Package is provided to the Licensee on the
!    understanding that the details thereof are confidential.
! 3. All publications issued by the Licensee that include results obtained
!    with the help of one or more of the Packages shall acknowledge the
!    use of the Packages. The Licensee will notify the Numerical Analysis
!    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
! 4. The Packages may be modified by or on behalf of the Licensee
!    for such use in research applications but at no time shall such
!    Packages or modifications thereof become the property of the
!    Licensee. The Licensee shall make available free of charge to the
!    copyright holder for any purpose all information relating to
!    any modification.
! 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
!    direct or consequential loss or damage whatsoever arising out of
!    the use of Packages by the Licensee.
! *******************************************************************
!
! Original date 18 October 2001

! 12th July 2004 Version 1.0.0. Version numbering added.
! 28 February 2006. Version 1.1.0. To avoid memory leaks, ma48_initialize and
!        ma48_finalize leave pointer component null instead having zero size.
! 27 April 2006. Version 1.2.0. Component drop of type ma48_finfo changed to
!        integer. Calculation of INFO(4) in SUBROUTINE NONZER corrected for
!        the case of a triangular matrix.
! 28 June 2006. Version 1.3.0. Value of la in ma48_analyse increased to ensure
!        that the arrays irn, icn, val are big enough for the original matrix.
! 5 October 2006. Version 1.4.0. Argument factors of ma48_solve given
!        intent(in) instead of intent(inout).
! 1 December 2006. Version 1.4.1. Deallocates in ma48_finalize separated out.
! 14 December 2006. Version 1.5.1. Error return added for pivot failure during
!                 fast factorize.
! 14 December 2006. Version 2.0.0. Sizes of cntl, icntl, info, rinfo increased
!        for compatibility with new version of ma48.
!        Immediate return from MA48B/BD if LA is too small.
!        control%switch_mode add to control whether a switch to slow
!        mode made when fast mode given unsuitable pivot sequence.
!        Pointer components of the derived type made allocatable (ZD11)

module hsl_ma48_double
   use hsl_zd11_double
   implicit none
   integer, parameter, private :: wp = kind(0.0d0)
   private :: nonzer
   type ma48_factors
     private
      integer, allocatable :: keep(:)
      integer, allocatable :: irn(:)
      integer, allocatable :: jcn(:)
      real(wp), allocatable :: val(:)
      integer :: m
      integer :: n
      integer :: lareq
      integer :: partial
      integer :: ndrop
      integer :: first
   end type ma48_factors
   type ma48_control
      real(wp) :: multiplier
      real(wp) :: reduce
      real(wp) :: u
      real(wp) :: switch
      real(wp) :: drop
      real(wp) :: tolerance
      real(wp) :: cgce
      integer :: lp
      integer :: wp
      integer :: mp
      integer :: ldiag
      integer :: btf
      logical :: struct
      integer :: maxit
      integer :: factor_blocking
      integer :: solve_blas
      integer :: la
      integer :: maxla
      integer :: pivoting
      logical :: diagonal_pivoting
      integer :: fill_in
      logical :: switch_mode
   end type ma48_control
   type ma48_ainfo
      real(wp) :: ops
      integer :: flag
      integer :: more
      integer :: len_analyse
      integer :: len_factorize
      integer :: ncmpa
      integer :: rank
      integer :: drop
      integer :: struc_rank
      integer :: oor
      integer :: dup
      integer :: stat
      integer :: lblock
      integer :: sblock
      integer :: tblock
   end type ma48_ainfo
   type ma48_finfo
      real(wp) :: ops
      integer :: flag
      integer :: more
      integer :: size_factor
      integer :: len_factorize
      integer :: drop
      integer :: rank
      integer :: stat
   end type ma48_finfo
   type ma48_sinfo
      integer :: flag
      integer :: more
      integer :: stat
   end type ma48_sinfo
contains
   subroutine ma48_initialize(factors,control)
      type(ma48_factors), intent(out), optional :: factors
      type(ma48_control), intent(out), optional :: control
      if (present(factors)) then
        factors%n = 0
        factors%first = 0
      end if
      if (present(control)) then
          control%switch = 0.5d0
          control%u      = 0.01d0
          control%drop   = 0.0d0
          control%tolerance = 0.0d0
          control%cgce      = 0.5d0
          control%lp = 6
          control%wp = 6
          control%mp = 6
          control%ldiag = 2
          control%pivoting = 3
          control%diagonal_pivoting = .false.
          control%fill_in = 3
          control%maxit = 10
          control%struct = .false.
          control%factor_blocking = 32
          control%solve_blas  = 2
          control%btf = 1
          control%la  = 0
          control%maxla = huge(0)
          control%multiplier = 2.0d0
          control%reduce     = 2.0d0
          control%switch_mode = .false.
      end if
    end subroutine ma48_initialize
   subroutine ma48_analyse(matrix,factors,control,ainfo,finfo,perm,lastcol)
      type(zd11_type), Intent(in) :: matrix
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      type(ma48_ainfo), intent(out) :: ainfo
      type(ma48_finfo), intent(out), optional :: finfo
      integer, intent(in), optional :: perm(matrix%m+matrix%n)
      integer, intent(in), optional :: lastcol(matrix%n)
      integer, allocatable :: iwork(:)
      integer :: i,job,k,la,lkeep,m,n,ne,stat,icntl(20),info(20)
      real(wp):: rinfo(10),cntl(10)
      external ma48ad
      cntl(1) = control%switch
      cntl(2) = control%u
      cntl(3) = control%drop
      cntl(4) = control%tolerance
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(4) = control%pivoting
      icntl(5) = control%factor_blocking
      icntl(6) = control%btf
      if (control%struct) then
        icntl(7) = 1
      else
        icntl(7) = 0
      end if
      m = matrix%m
      n = matrix%n
      ne = matrix%ne
      if(matrix%m .lt. 1) then
         ainfo%flag = -1
         ainfo%more = matrix%m
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%M has the value', matrix%m
       return
      end if
      if(matrix%n .lt. 1) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         ainfo%flag = -2
         ainfo%more = matrix%n
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%N has the value', matrix%n
       return
      end if
      if(matrix%ne .lt. 0) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         ainfo%flag = -3
         ainfo%more = matrix%ne
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'MATRIX%NE has the value', matrix%ne
       return
      end if
      ainfo%flag = 0
      ainfo%more = 0
      ainfo%stat = 0
      if(matrix%ne .eq. 0) then
        ainfo%ops = 0.0d0
        ainfo%rank = 0
        ainfo%drop = 0
        factors%ndrop = 0
        ainfo%oor = 0
        ainfo%dup = 0
        ainfo%len_analyse = 0
        ainfo%len_factorize = 0
        ainfo%struc_rank = 0
        ainfo%rank = 0
        ainfo%ncmpa = 0
        factors%first = 1
        factors%lareq = 0
        factors%partial = 0
        factors%m = m
        factors%n = n
        if (control%struct) then
          ainfo%flag = -5
           if (control%ldiag>0 .and. control%lp>=0 ) &
                 write (control%lp,'(/a,i3/a,i5)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
                 'Matrix is structurally singular with rank ',ainfo%struc_rank
        else
          ainfo%flag = 4
          if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_ANALYSE: ', &
              'ainfo%flag is equal to ',ainfo%flag
        endif
        return
      endif
      lkeep = m+5*n+4*n/icntl(6)+7
      if(allocated(factors%keep)) then
         if(size(factors%keep)/=lkeep) then
            deallocate(factors%keep,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%keep(lkeep),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%keep(lkeep),stat=stat)
         if (stat/=0) go to 100
      end if
      la = max(control%la,ne,control%fill_in*ne)
      if(allocated(factors%val)) then
         if(la /= size(factors%val)) then
            deallocate(factors%irn,factors%jcn,factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate(factors%irn(la),factors%jcn(la),factors%val(la),stat=stat)
            if (stat/=0) go to 100
         end if
      else
         allocate(factors%irn(la),factors%jcn(la),factors%val(la),stat=stat)
         if (stat/=0) go to 100
      end if
      if (present(perm)) then
         allocate (iwork(max(m,n)),stat=stat)
         if (stat/=0) go to 100
         iwork = 0
         do i = 1,m
           k = perm(i)
           if (k.gt.m .or. k.lt.0 .or. iwork(k) .ne. 0) then
             ainfo%flag = -6
             ainfo%more = i
             if (control%ldiag>0 .and. control%lp>=0 ) &
               write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'Invalid row permutation'
               deallocate (iwork,stat=stat)
             return
           end if
           iwork(k) = i
         end do
         iwork = 0
         do i = 1,n
           k = perm(m+i)
           if (k.gt.n .or. k.lt.0 .or. iwork(k) .ne. 0) then
             ainfo%flag = -6
             ainfo%more = i
             if (control%ldiag>0 .and. control%lp>=0 ) &
               write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'Invalid column permutation'
               deallocate (iwork,stat=stat)
             return
           end if
           iwork(k) = i
         end do
         deallocate (iwork,stat=stat)
         if (stat/=0) go to 100
         job = 2
      else
         job = 1
      end if
      if (control%diagonal_pivoting) job = 3
      allocate (iwork(6*m+3*n),stat=stat)
      if (stat/=0) go to 100
      do
         icntl(8) = 0
         if (present(lastcol)) then
           icntl(8) = 1
           factors%partial = count(lastcol(1:n) == 0)
           iwork(1:n) = lastcol(1:n)
         end if
         if (present(perm)) factors%keep(1:m+n) = perm(1:m+n)
         factors%irn(1:ne) = matrix%row(1:ne)
         factors%jcn(1:ne) = matrix%col(1:ne)
         factors%val(1:ne) = matrix%val(1:ne)
         info(4) = 0
         call ma48ad(m,n,ne,job,la,factors%val,factors%irn,  &
                     factors%jcn, &
                     factors%keep,cntl,icntl,iwork,info,rinfo)
         if (info(1).eq.-4) then
           ainfo%flag = -5
           ainfo%struc_rank = info(10)
           if (control%ldiag>0 .and. control%lp>=0 ) &
                 write (control%lp,'(/a,i3/a,i5)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
                 'Matrix is structurally singular with rank ',info(10)
         endif
         if (info(1)==-3) then
            la = control%multiplier*la
            la = max(la,info(3))
            if (la>control%maxla) then
               ainfo%flag = -7
               if (control%ldiag>0 .and. control%lp>=0 ) &
                 write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
                 'Main array needs to be bigger than', control%maxla
               return
            end if
            ainfo%more = la
            deallocate (factors%val,factors%irn,factors%jcn,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%val(la),factors%irn(la),factors%jcn(la),stat=stat)
            if (stat/=0) go to 100
         else
            exit
         end if
      end do
      if (info(1).gt.0) then
        if (info(1).eq.1 .and. info(11).gt.0) ainfo%flag = ainfo%flag + 2
        if (info(1).eq.1 .and. info(12).gt.0) ainfo%flag = ainfo%flag + 1
        if (info(1).eq.2) ainfo%flag = ainfo%flag + 4
        if (info(1).eq.4) ainfo%flag = ainfo%flag + 8
      endif
      if (ainfo%more.gt.0) ainfo%flag = ainfo%flag + 16
      factors%m = m
      factors%n = n
      factors%lareq = info(4)
      factors%first = 1
      ainfo%ops    = rinfo(1)
      ainfo%rank   = info(5)
      ainfo%drop   = info(6)
      factors%ndrop   = info(6)
      ainfo%oor    = info(12)
      ainfo%dup    = info(11)
      ainfo%len_analyse = info(3)
      ainfo%len_factorize = info(4)
      ainfo%lblock = info(7)
      ainfo%sblock = info(8)
      ainfo%tblock = info(9)
      ainfo%struc_rank = info(10)
      ainfo%ncmpa  = info(2)
      deallocate (iwork, stat=stat)
      if (stat/=0) go to 100
      if (present(finfo)) call MA48_factorize(matrix,factors,control,finfo)
      if (ainfo%flag.gt.0) then
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_ANALYSE: ', &
              'ainfo%flag is equal to ',ainfo%flag
      end if
      return
  100  ainfo%flag = -4
       ainfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_ANALYSE with ainfo%flag = ',ainfo%flag, &
         'allocate or deallocate failed with STAT=',stat
   end subroutine ma48_analyse
   subroutine ma48_get_perm(factors,perm)
      type(ma48_factors), intent(in), optional :: factors
      integer, intent(out) :: perm(:)
      integer m,n
      m = factors%m
      n = factors%n
      perm(1:m+n) = factors%keep(1:m+n)
    end subroutine ma48_get_perm
   subroutine ma48_factorize(matrix,factors,control,finfo,fast,partial)
      type(zd11_type), intent(in) :: matrix
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      type(ma48_finfo), intent(out) :: finfo
      integer, optional, intent(in) :: fast,partial
      integer :: job,la,m,n,ne
      integer stat
      integer icntl(20),info(20)
      real(wp) cntl(10),rinfo(10)
      real(wp), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      integer, allocatable :: itemp(:)
      external ma48bd
      if (factors%first .le. 0) then
         finfo%flag = -10
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'No prior call to MA48_ANALYSE'
         return
      endif
      m = matrix%m
      n = matrix%n
      ne = matrix%ne
      allocate (iwork(2*m+2*n),work(m),stat=stat)
      if (stat/=0) go to 100
      if(factors%m/=matrix%m) then
       finfo%flag = -1
       finfo%more = matrix%m
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12,a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'MATRIX%M has the value', matrix%m,' instead of',factors%m
       return
      end if
      if(factors%n/=matrix%n) then
        finfo%flag = -2
        finfo%more = matrix%n
          if (control%ldiag>0 .and. control%lp>=0 ) &
            write (control%lp,'(/a,i3/a,i12,a,i12)') &
           'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
           'MATRIX%N has the value', matrix%n,' instead of',factors%n
        return
      end if
      if(matrix%ne .lt. 0) then
         finfo%flag = -3
         finfo%more = matrix%ne
         if (control%ldiag>0 .and. control%lp>=0 ) &
            write (control%lp,'(/a,i3/a,i12)') &
           'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
           'MATRIX%NE has the value', matrix%ne
       return
      end if
      finfo%flag = 0
      finfo%more = 0
      finfo%stat = 0
      if(matrix%ne .eq. 0) then
        finfo%ops = 0.0d0
        finfo%rank = 0
        finfo%drop = 0
        finfo%len_factorize = 0
        finfo%size_factor = 0
        factors%first = 2
        factors%lareq = 0
        factors%partial = 0
        finfo%flag = 4
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
          'Warning from MA48_FACTORIZE: ', &
          'finfo%flag is equal to ',finfo%flag
        return
      endif
      la = size(factors%val)
      if (la<factors%lareq) then
            la = factors%lareq
            deallocate (factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%val(la),stat=stat)
            if (stat/=0) go to 100
            allocate (itemp(la),stat=stat)
            if (stat/=0) go to 100
            itemp(1:ne) = factors%irn(1:ne)
            deallocate (factors%irn,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%irn(la),stat=stat)
            if (stat/=0) go to 100
            factors%irn(1:ne) = itemp(1:ne)
            deallocate (itemp,stat=stat)
            if (stat/=0) go to 100
      endif
      if (la>control%reduce*factors%lareq) then
        la = factors%lareq
        deallocate (factors%val,stat=stat)
        if (stat/=0) go to 100
        allocate (factors%val(la),stat=stat)
        if (stat/=0) go to 100
        allocate (itemp(la),stat=stat)
        if (stat/=0) go to 100
        itemp(1:ne) = factors%irn(1:ne)
        deallocate (factors%irn,stat=stat)
        if (stat/=0) go to 100
        allocate (factors%irn(la),stat=stat)
        if (stat/=0) go to 100
        factors%irn(1:ne) = itemp(1:ne)
        deallocate (itemp,stat=stat)
        if (stat/=0) go to 100
      endif
      cntl(1) = control%switch
      cntl(2) = control%u
      cntl(3) = control%drop
      cntl(4) = control%tolerance
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(4) = control%pivoting
      icntl(5) = control%factor_blocking
      icntl(6) = control%btf
      icntl(10) = 1
      icntl(11) = 0
      if (control%switch_mode) icntl(11) = 1
      job = 1
      if (present(fast) .and. factors%ndrop==0 .and. factors%first.gt.1) &
          job = 2
      icntl(8) = 0
      if (present(partial) .and. factors%ndrop==0 .and. factors%partial.gt.0) &
      then
        job = 3
        icntl(8) = 1
      endif
      do
         factors%val(1:ne) = matrix%val(1:ne)
         call ma48bd(factors%m,factors%n,matrix%ne,job,la,factors%val, &
                  factors%irn,factors%jcn, &
                  factors%keep,cntl,icntl,work,iwork, &
                  info,rinfo)
         if (info(1)>=0) exit
         if (info(1)==-3) then
            la = control%multiplier*la
            if (la>control%maxla) then
               finfo%flag = -7
               if (control%ldiag>0 .and. control%lp>=0 ) &
                 write (control%lp,'(/a,i3/a,i12)') &
                 'Error return from MA48_FACTORIZE with finfo%flag = ',  &
                 finfo%flag, &
                 'Main array needs to be bigger than', control%maxla
               return
            end if
            finfo%more = la
            deallocate (factors%val,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%val(la),stat=stat)
            if (stat/=0) go to 100
            allocate (itemp(la),stat=stat)
            if (stat/=0) go to 100
            itemp(1:ne) = factors%irn(1:ne)
            deallocate (factors%irn,stat=stat)
            if (stat/=0) go to 100
            allocate (factors%irn(la),stat=stat)
            if (stat/=0) go to 100
            factors%irn(1:ne) = itemp(1:ne)
            deallocate (itemp,stat=stat)
            if (stat/=0) go to 100
          else if (info(1)==-7) then
            finfo%flag = -11
            if (control%ldiag>0 .and. control%lp>=0 ) &
               write (control%lp,'(/a,i3/a,i12)') &
               'Error return from MA48_FACTORIZE with finfo%flag = ',&
               finfo%flag,' Matrix entries unsuitable for fast factorization'
            return
          end if
      end do
      if (info(1).eq.2)    finfo%flag = finfo%flag + 4
      if (finfo%more.gt.0) finfo%flag = finfo%flag + 16
      deallocate (work,iwork,stat=stat)
      if (stat/=0) go to 100
      if (info(1)>=0) then
        finfo%len_factorize = info(4)
        factors%lareq = info(4)
        finfo%rank   = info(5)
        finfo%drop   = info(6)
        finfo%ops    = rinfo(1)
        call nonzer(m,n,factors%keep,info)
        finfo%size_factor = info(2)+info(3)+info(4)
        factors%first = 2
      end if
      if (finfo%flag.gt.0) then
        if (control%ldiag>1 .and. control%wp>=0) &
          write (control%wp,'(/a/a,i5)') &
              'Warning from MA48_FACTORIZE: ', &
              'finfo%flag is equal to ',finfo%flag
      end if
      return
  100  finfo%flag = -4
       finfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_FACTORIZE with finfo%flag = ',finfo%flag, &
         'allocate or deallocate failed with STAT=',stat
   end subroutine ma48_factorize
   subroutine ma48_solve(matrix,factors,rhs,x,control,sinfo,trans, &
                         resid,error)
      type(zd11_type), intent(in) :: matrix
      type(ma48_factors), intent(in) :: factors
      real(wp), intent(in) :: rhs(:)
      real(wp), intent(out) :: x(:)
      type(ma48_control), intent(in) :: control
      type(ma48_sinfo), intent(out) :: sinfo
      integer, optional, intent(in) :: trans
      real(wp), optional, intent(out) :: resid(2)
      real(wp), optional, intent(out) :: error
      integer icntl(20),info(20),job,m,n,stat
      real(wp) cntl(10),err(3)
      logical trans48
      integer, allocatable :: iwork(:)
      real(wp), allocatable :: work(:)
      real(wp), allocatable :: rhswork(:)
      external ma48cd
      if (factors%first .le. 1) then
         sinfo%flag = -10
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
         'No prior call to MA48_FACTORIZE'
         return
      endif
      m = matrix%m
      n = matrix%n
      if(factors%m/=matrix%m) then
         sinfo%flag = -1
         sinfo%more = matrix%m
         if (control%ldiag>0 .and. control%lp>=0 ) &
           write (control%lp,'(/a,i3/a,i12,a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%M has the value', matrix%m,' instead of',factors%m
       return
      end if
      if(factors%n/=matrix%n) then
         sinfo%flag = -2
         sinfo%more = matrix%n
         if (control%ldiag>0 .and. control%lp>=0 ) &
           write (control%lp,'(/a,i3/a,i12,a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%N has the value', matrix%n,' instead of',factors%n
       return
      end if
      if(matrix%ne .lt. 0) then
         sinfo%flag = -3
         sinfo%more = matrix%ne
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
          'MATRIX%NE has the value', matrix%ne
       return
      end if
      sinfo%flag = 0
      sinfo%more = 0
      sinfo%stat = 0
      if(matrix%ne .eq. 0) then
        x = 0.0d0
        if (present(resid)) resid = 0.0d0
        if (present(error)) error = 0.0d0
        return
      endif
      trans48 = present(trans)
      allocate (iwork(max(m,n)),work(max(3*m+n,3*n+m)),stat=stat)
      if (stat/=0) go to 100
      cntl(5)  = control%cgce
      icntl(1) = -1
      icntl(2) = -1
      if (control%ldiag.gt.2) icntl(2) = control%mp
      icntl(3) = control%ldiag
      icntl(5) = 0
      if (control%solve_blas.gt.1) icntl(5) = 2
      icntl(9) = control%maxit
      stat = 0
      job = 1
      if (control%maxit .gt. 0) then
        if (present(resid)) then
          if (present(error)) then
            job = 4
          else
            if (control%maxit == 1) then
              job = 2
            else
              job = 3
            end if
          end if
        end if
      end if
      if (trans48) then
        allocate (rhswork(n),stat=stat)
        if (stat/=0) go to 100
        rhswork = rhs(1:n)
      else
        allocate (rhswork(m),stat=stat)
        if (stat/=0) go to 100
        rhswork = rhs(1:m)
      endif
      call ma48cd(factors%m,factors%n,trans48,job,size(factors%val), &
                  factors%val, &
                  factors%irn,factors%keep,cntl,icntl,rhswork,   &
                  x,err,work,iwork,info)
      sinfo%flag = info(1)
      if(sinfo%flag .lt. 0) then
         if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3)') &
          'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag
       return
      end if
      if (job>1) resid(1) = err(1)
      if (job>1) resid(2) = err(2)
      if (job>3) error    = err(3)
      deallocate (iwork,work,rhswork,stat=stat)
      if (stat==0) return
  100  sinfo%flag = -4
       sinfo%stat = stat
       if (control%ldiag>0 .and. control%lp>=0 ) &
         write (control%lp,'(/a,i3/a,i12)') &
         'Error return from MA48_SOLVE with sinfo%flag = ',sinfo%flag, &
         'allocate or deallocate failed with STAT=',stat
   end subroutine ma48_solve
   subroutine ma48_finalize(factors,control,info)
      type(ma48_factors), intent(inout) :: factors
      type(ma48_control), intent(in) :: control
      integer, intent(out) :: info
      info = 0
      if (allocated(factors%keep)) deallocate(factors%keep,stat=info)
      if (info==0 .and. allocated(factors%val))  deallocate &
                       (factors%irn,factors%jcn,factors%val,stat=info)
      if (info==0) return
      if (control%ldiag>0 .and. control%lp>=0 ) write (control%lp,'(/2a,i5)') &
         'Error return from MA48_finalize: ',&
         'allocate or deallocate failed with STAT=',info
    end subroutine ma48_finalize
      SUBROUTINE NONZER(M,N,KEEP,INFO)
      INTEGER M,N
      INTEGER KEEP(*),INFO(4)
      INTEGER IPTRL,IPTRU,IPTRD,IPTRO,NBLOCK,MBLOCK,KBLOCK
      INTEGER NB,KB,JB,J1,J2,NC,NR,LC
      IPTRL = M+N
      IPTRU = IPTRL + N
      IPTRD = IPTRU + N
      IPTRO = IPTRD + N + 1
      NBLOCK = IPTRO + N - 1
      MBLOCK = NBLOCK + 1
      KBLOCK = MBLOCK + 1
      NB = KEEP(KBLOCK+3)
      INFO(1) = KEEP(IPTRO+N+1) - 1
      INFO(2) = KEEP(IPTRO+N+1) - KEEP(IPTRO+1)
      INFO(3) = 0
      J2 = 0
      LC = 0
      DO 100 JB = 1,NB
        NC = KEEP(NBLOCK+3*JB)
        J1 = J2 + 1
        J2 = J1 + NC - 1
        IF (KEEP(MBLOCK+3*JB).LT.0) THEN
          INFO(3) = INFO(3) + KEEP(IPTRD+J2+1) - KEEP(IPTRD+J1)
        ELSE
          KB = JB
          NR = NC
          LC = J2
        ENDIF
  100 CONTINUE
      IF (LC.EQ.0) THEN
        INFO(4) = 0
      ELSE
        NC = NR
        IF (NB.EQ.1) THEN
          NR = M
          INFO(4) = 0
        ELSE
          INFO(4) = KEEP(KBLOCK+3*KB) - 1
        ENDIF
        INFO(4) = INFO(4) + KEEP(IPTRL+LC) +  &
                  MAX(((NC-KEEP(MBLOCK+3*KB))+(NR-KEEP(MBLOCK+3*KB))), &
                      ((NC-KEEP(MBLOCK+3*KB))*(NR-KEEP(MBLOCK+3*KB))))
      ENDIF
      end subroutine nonzer
  subroutine ma48_special_rows_and_cols(factors,rank,rows,cols,info)
      type(ma48_factors), intent(in) :: factors
      integer,intent(out) :: rank,info
      integer,intent(out),dimension(factors%m) :: rows
      integer,intent(out),dimension(factors%n) :: cols
      integer :: la
      integer, allocatable :: iwork(:)
      external ma51ad
      allocate (iwork(max(factors%m,factors%n)),stat=info)
      if (info/=0) then
        info = -1
        return
      end if
      la = size(factors%val)
      call ma51ad(factors%m,factors%n,la,factors%irn,factors%keep,rank, &
                  rows,cols,iwork)
      deallocate (iwork, stat=info)
      if (info/=0) info = -2
   end subroutine ma48_special_rows_and_cols
   subroutine ma48_determinant(factors,sgndet,logdet,info)
      type(ma48_factors), intent(in) :: factors
      integer,intent(out) :: sgndet,info
      real(wp),intent(out) :: logdet
      integer :: la
      integer, allocatable :: iwork(:)
      external ma51cd
      allocate (iwork(factors%n),stat=info)
      if (info/=0) then
        info = -1
        return
      end if
      la = size(factors%val)
      call ma51cd(factors%m,factors%n,la,factors%val,factors%irn, &
                  factors%keep,sgndet,logdet,iwork)
      deallocate (iwork, stat=info)
      if (info/=0) info = -2
  end subroutine ma48_determinant
end module hsl_ma48_double
