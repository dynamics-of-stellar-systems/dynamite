! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*- G A L A H A D _ U L S    M O D U L E  -*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   development started, April 3rd 2006
!   originally released with GALAHAD Version 2.0. April 25th 2006

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_ULS_double

!     ----------------------------------------
!     |                                      |
!     |  Provide a MA48-style interface for  |
!     |  MA33 to allow the solution of       |
!     |                                      |
!     |      Unsymmetric Linear Systems      |
!     |                                      |
!     ----------------------------------------

     USE GALAHAD_SMT_double

     IMPLICIT NONE
 
     PRIVATE
     PUBLIC :: ULS_initialize, ULS_analyse, ULS_solve, ULS_finalize,           &
               ULS_special_rows_and_cols, SMT_type
   
!--------------------
!   P r e c i s i o n
!--------------------

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  

!  Set other parameters

     REAL ( KIND = wp ), PRIVATE, PARAMETER :: zero = 0.0_wp
     REAL ( KIND = wp ), PRIVATE, PARAMETER :: point01 = 0.01_wp
     REAL ( KIND = wp ), PRIVATE, PARAMETER :: half = 0.5_wp
     REAL ( KIND = wp ), PRIVATE, PARAMETER :: one = 1.0_wp
     REAL ( KIND = wp ), PRIVATE, PARAMETER :: two = 2.0_wp
     REAL ( KIND = wp ), PRIVATE, PARAMETER :: eps = EPSILON( one )

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

     TYPE, PUBLIC :: ULS_factors
       PRIVATE
       LOGICAL :: got_factors = .FALSE.
! scalars and arrays required by MA33
       INTEGER :: n, licn, lirn, orig_m, orig_n
       REAL ( KIND = wp ) :: u
       INTEGER, DIMENSION( 2 ) :: IDISP
       INTEGER, DIMENSION( 10 ) :: ICNTL, ICNTL_MA23
       INTEGER, DIMENSION( 10 ) :: INFO, INFO_MA23
       REAL ( KIND = wp ), DIMENSION( 5 ) :: CNTL
       REAL ( KIND = wp ), DIMENSION( 5 ) :: RINFO
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: ICN          ! length licn
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IFIRST       ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IP           ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IPC          ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IPTR         ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IQ           ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IRN          ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LASTC        ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LASTR        ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LENC         ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LENR         ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LENRL        ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: NEXTC        ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: NEXTR        ! length n
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: LENOFF       ! length 1/n 
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: A ! length licn
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: W ! length n
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R ! length n
     END TYPE ULS_factors

     TYPE, PUBLIC :: ULS_control
       REAL ( KIND = wp ) :: multiplier ! Factor by which arrays sizes are to be
                         ! increased if they are too small                    NEW
       REAL ( KIND = wp ) :: reduce ! if previously allocated internal workspace 
                         !  arrays are greater than reduce times the currently
                         !  required sizes, they are reset to current requirments
                         !                                                    NEW
       REAL ( KIND = wp ) :: u     ! Pivot threshold
       REAL ( KIND = wp ) :: switch ! Density for switch to full code         NEW
       REAL ( KIND = wp ) :: drop   ! Drop tolerance
       REAL ( KIND = wp ) :: tolerance ! anything < this is considered zero
       REAL ( KIND = wp ) :: cgce  ! Ratio for required reduction using IR
       INTEGER :: lp     ! Unit for error messages
       INTEGER :: wp     ! Unit for warning messages                          NEW
       INTEGER :: mp     ! Unit for monitor output
       INTEGER :: ldiag  ! Controls level of diagnostic output                NEW
       INTEGER :: btf    ! Minimum block size for BTF ... >=N to avoid        NEW
       LOGICAL :: struct ! Control to abort if structurally singular       ABORT1
       INTEGER :: maxit ! Maximum number of iterations
       INTEGER :: factor_blocking ! Level 3 blocking in factorize             NEW
       INTEGER :: solve_blas ! Switch for using Level 1 or 2 BLAS in solve.   NEW
       INTEGER :: la     ! Initial size for real array for the factors.       NEW
       INTEGER :: maxla  ! Max. size for real array for the factors.          NEW
       INTEGER :: pivoting  ! Controls pivoting:
                            !  Number of columns searched.  Zero for Markowitz
                            !                                               NSRCH
       LOGICAL :: diagonal_pivoting  ! Set to 0 for diagonal pivoting         NEW
       INTEGER :: fill_in ! Initially fill_in * ne space allocated for factors
     END TYPE ULS_control

     TYPE, PUBLIC :: ULS_ainfo
       REAL ( KIND = wp ) :: ops   ! Number of operations in elimination      NO
       INTEGER :: flag   ! Flags success or failure case                     
       INTEGER :: more    ! More information on failure
       INTEGER :: len_analyse ! Size for analysis            
       INTEGER :: len_factorize  ! Size for factorize
       INTEGER :: ncmpa   ! Number of compresses
       INTEGER :: rank    ! Estimated rank
       INTEGER :: drop    ! Number of entries dropped
       INTEGER :: struc_rank ! Structural rank of matrix
       INTEGER :: oor     ! Number of indices out-of-range                    NO
       INTEGER :: dup     ! Number of duplicates                              NO
       INTEGER :: stat    ! STAT value after allocate failure
       INTEGER :: lblock  ! Size largest non-triangular block                 NO
       INTEGER :: sblock  ! Sum of orders of non-triangular blocks            NO
       INTEGER :: tblock  ! Total entries in all non-tringular blocks         NO
     END TYPE ULS_ainfo

     TYPE, PUBLIC :: ULS_finfo
       REAL ( KIND = wp ) :: ops   ! Number of operations in elimination      NO
       INTEGER :: flag    ! Flags success or failure case
       INTEGER :: more    ! More information on failure
       INTEGER :: size_factor  ! Number of words to hold factors
       INTEGER :: len_factorize  ! Size for subsequent factorization
       INTEGER :: drop    ! Number of entries dropped
       INTEGER :: rank    ! Estimated rank
       INTEGER :: stat    ! STAT value after allocate failure
     END TYPE ULS_finfo

     TYPE, PUBLIC :: ULS_sinfo
       INTEGER :: flag    ! Flags success or failure case
       INTEGER :: more    ! More information on failure
       INTEGER :: stat    ! STAT value after allocate failure
     END TYPE ULS_sinfo

!--------------------------------
!   I n t e r f a c e  B l o c k
!--------------------------------

     INTERFACE 
     
       SUBROUTINE MA33ID( ICNTL, CNTL )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  
       INTEGER, INTENT( OUT ), DIMENSION( 10 ) :: ICNTL
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 5 ) :: CNTL
       END SUBROUTINE MA33ID

       SUBROUTINE MA33AD( n, ICN, A, licn, LENR, LENRL, IDISP, IP, IQ, IRN,    &
                          lirn, LENC, IFIRST, LASTR, NEXTR, LASTC, NEXTC,      &
                          IPTR, IPC, u, iflag, ICNTL, CNTL, INFO, RINFO )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  
       INTEGER, INTENT( IN ) :: n, licn, lirn
       INTEGER, INTENT( OUT ) :: iflag
       REAL ( KIND = wp ), INTENT( IN ) :: u
       INTEGER, INTENT( INOUT ), DIMENSION( licn ) ::  ICN
       INTEGER, INTENT( INOUT ), DIMENSION( n ) ::  LENR
       INTEGER, INTENT( OUT ), DIMENSION( n ) ::  LENRL
       INTEGER, INTENT( INOUT ), DIMENSION( 2 ) ::  IDISP
       INTEGER, INTENT( INOUT ), DIMENSION( n ) ::  IP, IQ
       INTEGER, INTENT( OUT ), DIMENSION( lirn ) :: IRN
       INTEGER, INTENT( OUT ), DIMENSION( n ) ::  IPC, IPTR, LENC, IFIRST
       INTEGER, INTENT( OUT ), DIMENSION( n ) ::  LASTR, NEXTR, LASTC, NEXTC
       INTEGER, INTENT( OUT ), DIMENSION( 10 ) :: INFO
       INTEGER, INTENT( IN ), DIMENSION( 10 ) :: ICNTL
       REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( licn ) :: A
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( 5 ) :: RINFO
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( 5 ) :: CNTL
       END SUBROUTINE MA33AD

       SUBROUTINE MA33CD( n, ICN, A, licn, LENR, LENRL, LENOFF, IDISP,         &
                          IP, IQ, X, W, mtype, RINFO )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  
       INTEGER, INTENT( IN ) :: n, licn, mtype
       INTEGER, INTENT( IN ), DIMENSION( licn ) :: ICN
       INTEGER, INTENT( IN ), DIMENSION( n ) :: LENR, LENRL, LENOFF
       INTEGER, INTENT( IN ), DIMENSION( 2 ) :: IDISP
       INTEGER, INTENT( INOUT ), DIMENSION( n ) :: IP, IQ
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( licn ) :: A
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( n ) :: W
       REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: X
       REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( 5 ) :: RINFO
       END SUBROUTINE MA33CD

       SUBROUTINE MC20AD( nc, maxa, A, INUM, JPTR, JNUM, jdisp )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )  
       INTEGER, INTENT( IN ) :: nc, maxa, jdisp
       INTEGER, INTENT( INOUT ), DIMENSION( maxa ) :: INUM, JNUM
       INTEGER, INTENT( OUT ), DIMENSION( nc ) :: JPTR
       REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( maxa ) :: A
       END SUBROUTINE MC20AD

     END INTERFACE

   CONTAINS

!-*-*-*-*-*-   U L S _ I N I T I A L I Z E   S U B R O U T I N E   -*-*-*-*-*

     SUBROUTINE ULS_initialize( FACTORS, CONTROL )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
!  Default control data for ULS. This routine should be called before
!  first call to ULS_analyse
!
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

     TYPE( ULS_factors ), INTENT( out ), OPTIONAL :: FACTORS
     TYPE( ULS_control ), INTENT( out ), OPTIONAL :: CONTROL

     IF ( PRESENT( FACTORS ) ) THEN
       FACTORS%got_factors = .FALSE.
     END IF

     IF ( PRESENT( CONTROL ) ) THEN
       CONTROL%switch = half
       CONTROL%u = point01
       CONTROL%drop = zero
       CONTROL%tolerance = zero
       CONTROL%cgce = half
       CONTROL%lp = 6
       CONTROL%wp = 6
       CONTROL%mp = 6
       CONTROL%ldiag = 2
!      CONTROL%pivoting = 3
       CONTROL%pivoting = 32768
       CONTROL%diagonal_pivoting = .FALSE.
       CONTROL%fill_in = 3
       CONTROL%maxit = 10
       CONTROL%struct = .FALSE.
       CONTROL%factor_blocking = 32
       CONTROL%solve_blas = 2
       CONTROL%btf = 1
       CONTROL%la = 0
       CONTROL%maxla = HUGE( 0 )
       CONTROL%multiplier = two
       CONTROL%reduce = two
     END IF

     RETURN  

!  End of SILS_initialize

     END SUBROUTINE ULS_initialize

!-*-*-*-*-*-*-*-   U L S _ A N A L Y S E  S U B R O U T I N E   -*-*-*-*-*-*-

     SUBROUTINE ULS_analyse( MATRIX, FACTORS, CONTROL, AINFO, FINFO )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Analyse the sparsity pattern to obtain a good potential ordering
!  for the subsequent factorization
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

     TYPE ( SMT_type ), INTENT( IN ) :: MATRIX
     TYPE ( ULS_factors ), INTENT( INOUT ) :: FACTORS
     TYPE ( ULS_control ), INTENT( IN ) :: CONTROL
     TYPE ( ULS_ainfo ), INTENT( OUT ) :: AINFO
     TYPE ( ULS_finfo ), INTENT( OUT ) :: FINFO

!  Local variables

     INTEGER :: i, iflag, stat, lp, mp, nsrch, length, newpos
     INTEGER :: move, i1, ii, iend, j, j1, j2, jj, jay, knum, newj1
     REAL( KIND = wp ):: themax, u
     LOGICAL :: lblock

     lp = max( control%lp, control%wp )
     mp = control%mp

!  assign default informational parameters

     AINFO%ops = zero ; FINFO%ops = zero
     AINFO%flag = 0 ; AINFO%more = 0 ; AINFO%stat = 0
     FINFO%flag = 0 ; FINFO%more = 0 ; FINFO%stat = 0
     AINFO%len_analyse = 0 ; AINFO%len_factorize = 0
     FINFO%size_factor = 0 ; FINFO%len_factorize = 0
     AINFO%drop = 0 ; AINFO%rank = 0 ; AINFO%struc_rank = 0 ; AINFO%ncmpa = 0
     FINFO%drop = 0 ; FINFO%rank = 0
     AINFO%oor = 0 ; AINFO%dup = 0
     AINFO%lblock = 0 ; AINFO%sblock = 0 ; AINFO%tblock = 0

!  Check input parameters

     FACTORS%orig_m = MATRIX%m ; FACTORS%orig_n = MATRIX%n
     FACTORS%n = MAX( MATRIX%m, MATRIX%n )
     IF ( MATRIX%n < 1 ) THEN
       AINFO%flag = - 1
       IF ( lp /= 0 ) WRITE( lp, "( A, I0, A )" )                              &
         ' Error return from ULS_analyse because MATRIX%n = ', MATRIX%n, ' < 1'
       RETURN
     ELSE IF ( MATRIX%m < 1 ) THEN
       AINFO%flag = - 2
       IF ( lp /= 0 ) WRITE( lp, "( A, I0, A )" )                              &
         ' Error return from ULS_analyse because MATRIX%m = ', MATRIX%m, ' < 1'
       RETURN
     ELSE IF ( MATRIX%ne < 0 ) THEN
       AINFO%flag = - 3
       IF ( lp /= 0 ) WRITE( lp, "( A, I0, A )" )                              &
         ' Error return from ULS_analyse because MATRIX%ne = ', MATRIX%ne, ' < 1'
       RETURN
     ELSE IF ( MATRIX%ne == 0 ) THEN
       IF ( CONTROL%struct ) THEN
         AINFO%flag = - 5
          IF ( CONTROL%ldiag > 0 .AND. CONTROL%lp >= 0  )                      &
            WRITE( CONTROL%lp, '( /, A, I3, /, A, I5 )' )                      &
              ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag, &
              ' Matrix is structurally singular with rank ', AINFO%struc_rank
       ELSE
         AINFO%flag = 4
         IF ( CONTROL%ldiag > 1 .AND. CONTROL%wp >= 0 )                        &
           WRITE( CONTROL%wp, '( /, A, /, A, I5 )' )                           &
             ' Warning from ULS_ANALYSE: AINFO%flag is equal to ', AINFO%flag
       END IF
       RETURN
     END IF

!  Check that all the indices are within range

     DO i = 1, MATRIX%ne
       IF ( MATRIX%ROW( i ) <= 0 .OR. MATRIX%ROW( i ) > MATRIX%m .OR.          &
            MATRIX%COL( i ) <= 0 .OR. MATRIX%COL( i ) > MATRIX%n ) THEN
         AINFO%oor = AINFO%oor + 1
         IF ( control%wp /= 0 .AND. AINFO%oor <= 10 ) WRITE( control%wp,       &
            "( 1X, I0, A, ES12.4, /, 20X, A, I0, A, I0 )" )                    &
              i, 'th element with value ', MATRIX%VAL( i ),                    &
              ' is out of range with indices ', MATRIX%ROW( i ), ' , ',        &
             MATRIX%COL( i )
       END IF
     END DO

     IF ( AINFO%oor > 0 .AND. mp > 0 ) WRITE( mp, "( A, I0, A )" )             &
       ' Warning from ULS_analyse because', AINFO%oor,                         &
       ' indices found out of range'

!  Allocate workspace arrays prior to the analyse/factorise phases

     IF ( ALLOCATED( FACTORS%IFIRST ) ) THEN
       IF ( SIZE( FACTORS%IFIRST ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%IFIRST, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%IFIRST( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%IFIRST( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%IP ) ) THEN
       IF ( SIZE( FACTORS%IP ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%IP, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%IP( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%IP( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%IPC ) ) THEN
       IF ( SIZE( FACTORS%IPC ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%IPC, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%IPC( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%IPC( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%IPTR ) ) THEN
       IF ( SIZE( FACTORS%IPTR ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%IPTR, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%IPTR( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%IPTR( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%IQ ) ) THEN
       IF ( SIZE( FACTORS%IQ ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%IQ, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%IQ( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%IQ( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LASTC ) ) THEN
       IF ( SIZE( FACTORS%LASTC ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LASTC, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LASTC( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LASTC( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LASTR ) ) THEN
       IF ( SIZE( FACTORS%LASTR ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LASTR, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LASTR( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LASTR( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LENC ) ) THEN
       IF ( SIZE( FACTORS%LENC ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LENC, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LENC( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LENC( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LENR ) ) THEN
       IF ( SIZE( FACTORS%LENR ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LENR, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LENR( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LENR( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LENRL ) ) THEN
       IF ( SIZE( FACTORS%LENRL ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LENRL, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LENRL( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LENRL( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%NEXTC ) ) THEN
       IF ( SIZE( FACTORS%NEXTC ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%NEXTC, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%NEXTC( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%NEXTC( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%NEXTR ) ) THEN
       IF ( SIZE( FACTORS%NEXTR ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%NEXTR, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%NEXTR( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%NEXTR( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%LENOFF ) ) THEN
       IF ( SIZE( FACTORS%LENOFF ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%LENOFF, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%LENOFF( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%LENOFF( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     IF ( ALLOCATED( FACTORS%W ) ) THEN
       IF ( SIZE( FACTORS%W ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%W, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%W( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%W( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

!  assign control parameters for MA33

      FACTORS%ICNTL( 1 ) = lp
      IF ( control%struct ) THEN
        FACTORS%ICNTL( 2 ) = 1
      ELSE
        FACTORS%ICNTL( 2 ) = 0
      END IF
      FACTORS%ICNTL( 3 ) = 0
      FACTORS%ICNTL( 4 ) = 0
      FACTORS%ICNTL( 5 ) = control%pivoting
      FACTORS%ICNTL( 6 ) = 0

      FACTORS%CNTL( 1 ) = two
      FACTORS%CNTL( 2 ) = control%DROP

      lblock = control%btf < FACTORS%n .AND. MATRIX%m == MATRIX%n
      nsrch = control%pivoting
      IF ( nsrch <= 0 ) nsrch = FACTORS%n

!  Set initial estimates for the space required for the factorization

     FACTORS%licn = MAX( CONTROL%la, MAX( 1, CONTROL%fill_in ) * MATRIX%ne )
     FACTORS%lirn = MAX( CONTROL%la, MAX( 1, CONTROL%fill_in ) * MATRIX%ne )

!  Prepare for the factorization

     DO

!  Allocate space for the factors

       IF ( ALLOCATED( FACTORS%IRN ) ) THEN
         IF ( SIZE( FACTORS%IRN ) /= FACTORS%lirn ) THEN
           DEALLOCATE( FACTORS%IRN, stat = stat )
           IF ( stat /= 0 ) GO TO 100
           ALLOCATE( FACTORS%IRN( FACTORS%lirn ), stat = stat )
           IF ( stat /= 0 ) GO TO 100
         END IF
       ELSE
         ALLOCATE( FACTORS%IRN( FACTORS%lirn ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF

       IF ( ALLOCATED( FACTORS%ICN ) ) THEN
         IF ( SIZE( FACTORS%ICN ) /= FACTORS%licn ) THEN
           DEALLOCATE( FACTORS%ICN, stat = stat )
           IF ( stat /= 0 ) GO TO 100
           ALLOCATE( FACTORS%ICN( FACTORS%licn ), stat = stat )
           IF ( stat /= 0 ) GO TO 100
         END IF
       ELSE
         ALLOCATE( FACTORS%ICN( FACTORS%licn ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF

       IF ( ALLOCATED( FACTORS%A ) ) THEN
         IF ( SIZE( FACTORS%A ) /= FACTORS%licn ) THEN
           DEALLOCATE( FACTORS%A, stat = stat )
           IF ( stat /= 0 ) GO TO 100
           ALLOCATE( FACTORS%A( FACTORS%licn ), stat = stat )
           IF ( stat /= 0 ) GO TO 100
         END IF
       ELSE
         ALLOCATE( FACTORS%A( FACTORS%licn ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF

!  Copy the matrix into the space reserved for the factors

       FACTORS%IRN( 1 : MATRIX%ne ) = MATRIX%row( 1 : MATRIX%ne )
       FACTORS%ICN( 1 : MATRIX%ne ) = MATRIX%col( 1 : MATRIX%ne )
       FACTORS%A( 1 : MATRIX%ne ) = MATRIX%val( 1 : MATRIX%ne )

!  Switch from co-ordinate to row order

       CALL MC20AD( FACTORS%n, MATRIX%ne, FACTORS%A, FACTORS%ICN,              &
                    FACTORS%IPC, FACTORS%IRN, 0 )

!  Use LENR and IP as temporary workspace 

       FACTORS%LENR = 0
       FACTORS%IP = 0

!  Check for, and sum, duplicate entries

       move = 0
       themax = zero
       j1 = FACTORS%IPC( 1 )
       DO i = 1, FACTORS%n
         IF ( i < FACTORS%n ) THEN
           iend = FACTORS%IPC( i + 1 )
         ELSE
           iend = MATRIX%ne + 1
         END IF
         length = iend - j1
         IF ( length == 0 ) CYCLE
         j2 = iend - 1
         newj1 = j1 - move
         DO jj = j1, j2
           j = FACTORS%ICN( jj )
           themax = MAX( themax, ABS( FACTORS%A( jj ) ) )
           IF ( FACTORS%IP( j ) /= i ) THEN
             FACTORS%IP( j ) = i
             FACTORS%IQ( j ) = jj - move - newj1
             IF ( move /= 0 ) THEN
               newpos = jj - move
               FACTORS%A( newpos ) = FACTORS%A( jj )
               FACTORS%ICN( newpos ) = FACTORS%ICN( jj )
             END IF
             CYCLE
           END IF
           move = move + 1
           length = length - 1
           jay = FACTORS%IQ( j ) + newj1
           IF ( mp /= 0 ) WRITE( mp, "( A, I0, A, I0, A, ES12.4 )" )            &
            ' Duplicate element in position ', i, ', ', j, ' with value ',      &
             FACTORS%A( jj )
           FACTORS%A( jay ) = FACTORS%A( jay ) + FACTORS%A( jj )
           AINFO%dup = AINFO%dup + 1
           themax = MAX( themax, ABS( FACTORS%A( jay ) ) )
         END DO
         FACTORS%LENR( i ) = length
         j1 = iend
       END DO
       knum = MATRIX%ne - move

!  If required, permute to block-triangular form

       IF ( lblock ) THEN
         FACTORS%ICNTL_MA23( 1 ) = lp
         FACTORS%ICNTL_MA23( 2 ) = FACTORS%ICNTL( 2 )
         FACTORS%ICNTL_MA23( 3 : 10 ) = 0
         CALL MC23_threadsafe( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,&
                               FACTORS%LENR, FACTORS%IDISP, FACTORS%IP,        &
                               FACTORS%IQ, FACTORS%LENOFF, FACTORS%IFIRST,     &
                               FACTORS%LASTR, FACTORS%NEXTR,  FACTORS%LASTC,   &
                               FACTORS%IPC, FACTORS%LENC,                      &
                               FACTORS%ICNTL_MA23, FACTORS%INFO_MA23 )
 
         IF ( FACTORS%INFO_MA23( 1 ) < 0 ) THEN
           AINFO%flag = - 7
           IF ( FACTORS%INFO_MA23( 1 ) == - 1 ) AINFO%flag = - 5
           IF ( lp /= 0 ) WRITE( lp, "( A )" )                                 &
            ' Error return from ULS_analyse following call to MC23_threadsafe'
           RETURN
         END IF
       ELSE
         DO i = 1, knum
           ii = knum - i + 1
           newpos = FACTORS%licn - i + 1
           FACTORS%ICN( newpos ) = FACTORS%ICN( ii )
           FACTORS%A( newpos ) = FACTORS%A( ii )
         END DO
         FACTORS%IDISP( 1 ) = 1
         FACTORS%IDISP( 2 ) = FACTORS%licn - knum + 1
         DO i = 1, FACTORS%n
           FACTORS%IP( i ) = i
           FACTORS%IQ( i ) = i
         END DO
         FACTORS%LENOFF( 1 ) = - 1
       END IF

!  Factorize the matrix

       u = CONTROL%u
       IF ( nsrch > FACTORS%n ) THEN
         CALL MA33AD( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,         &
                      FACTORS%LENR, FACTORS%LENRL, FACTORS%IDISP, FACTORS%IP,  &
                      FACTORS%IQ, FACTORS%IRN, FACTORS%lirn, FACTORS%LENC,     &
                      FACTORS%IFIRST, FACTORS%LASTR, FACTORS%NEXTR,            &
                      FACTORS%LASTC, FACTORS%NEXTC, FACTORS%IPTR, FACTORS%IPC, &
                      u, iflag, FACTORS%ICNTL, FACTORS%CNTL,                   &
                      FACTORS%INFO, FACTORS%RINFO )
       ELSE
         CALL MA33AD( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,         &
                      FACTORS%LENR, FACTORS%LENRL, FACTORS%IDISP, FACTORS%IP,  &
                      FACTORS%IQ, FACTORS%IRN, FACTORS%lirn, FACTORS%LENC,     &
                      FACTORS%IFIRST, FACTORS%LASTR, FACTORS%NEXTR,            &
                      FACTORS%IPC, FACTORS%IPC, FACTORS%LASTC, FACTORS%IPC,    &
                      u, iflag, FACTORS%ICNTL, FACTORS%CNTL,                   &
                      FACTORS%INFO, FACTORS%RINFO )
       END IF

!  Record return information

       AINFO%flag = iflag
       AINFO%len_analyse = MAX( FACTORS%INFO( 4 ), FACTORS%INFO( 5 ), MATRIX%ne )
       AINFO%len_factorize = AINFO%len_analyse
       AINFO%ncmpa = FACTORS%INFO( 1 ) + FACTORS%INFO( 2 )
       AINFO%rank = FACTORS%INFO ( 3 )
       AINFO%drop = FACTORS%INFO( 6 )
       AINFO%struc_rank = AINFO%rank

       FINFO%flag = AINFO%flag
       FINFO%size_factor = AINFO%len_analyse
       FINFO%len_factorize = AINFO%len_analyse
       FINFO%drop = AINFO%drop
       FINFO%rank = AINFO%rank

!  Check for a succesful exit

       IF ( iflag >= 0 ) EXIT

!  Check for un-required singularity

       IF ( iflag == - 1 ) THEN
         AINFO%flag = - 5
         AINFO%struc_rank = FACTORS%INFO( 3 )
         IF ( CONTROL%ldiag>0 .AND. CONTROL%lp >= 0  )                         &
           WRITE ( CONTROL%lp, '( /, A, I3, /, A, I5 )' )                      &
            ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag,   &
            ' Matrix is structurally singular with rank ', FACTORS%INFO( 3 )
         IF ( CONTROL%struct ) THEN
           RETURN
         ELSE
           EXIT
         END IF
       END IF

       IF ( iflag == - 2 ) THEN
         AINFO%flag = - 5
         AINFO%rank = FACTORS%INFO( 3 )
         IF ( CONTROL%ldiag>0 .AND. CONTROL%lp >= 0  )                         &
           WRITE ( CONTROL%lp, '( /a, i3/a, i5 )' )                            &
            ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag,   &
            ' Matrix is singular with estimated rank ', FACTORS%INFO( 3 )
         EXIT
       END IF

!  If necessary, increase the space required for the factors

       IF ( iflag == - 3 .OR. iflag == - 6 ) THEN
         i = FACTORS%lirn * CONTROL%multiplier
         FACTORS%lirn = MAX( i, FACTORS%INFO( 4 ) )
         IF ( FACTORS%lirn > CONTROL%maxla ) THEN
           AINFO%flag = - 7
           IF ( CONTROL%ldiag>0 .and. CONTROL%lp>=0  )                         &
             WRITE ( CONTROL%lp, '( /a, i3/a, i12 )' )                         &
             ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag,  &
             ' Array FACTORS%IRN needs to be bigger than',  CONTROL%maxla
           RETURN
         END IF
       END IF

       IF ( iflag == - 4 .OR. iflag == - 5 .OR. iflag == - 6 ) THEN
         i = FACTORS%licn * CONTROL%multiplier
         FACTORS%licn = MAX( i, FACTORS%INFO( 5 ) )
         IF ( FACTORS%licn > CONTROL%maxla ) THEN
           AINFO%flag = - 7
           IF ( CONTROL%ldiag>0 .and. CONTROL%lp>=0  )                         &
             WRITE ( CONTROL%lp, '( /a, i3/a, i12 )' )                         &
             ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag,  &
             ' Array FACTORS%IRN needs to be bigger than',  CONTROL%maxla
           RETURN
         END IF
       END IF

!  The factorization is complete

     END DO

!  Reorder the off-diagonal blocks according to pivot permutation.

     i1 = FACTORS%IDISP( 1 ) - 1
     IF ( i1 /= 0 )                                                            &
       CALL MC22_threadsafe( FACTORS%n, FACTORS%ICN, FACTORS%A, i1,            &
                             FACTORS%LENOFF, FACTORS%IP, FACTORS%IQ,           &
                             FACTORS%IPC,FACTORS%LENC, FACTORS%IRN )
!  Record exit status

     AINFO%flag = 0
     IF ( AINFO%oor > 0 ) AINFO%flag = AINFO%flag + 1
     IF ( AINFO%dup > 0 ) AINFO%flag = AINFO%flag + 2
     IF ( FACTORS%n /= FINFO%rank ) AINFO%flag = AINFO%flag + 4
     IF ( AINFO%more > 0 ) AINFO%flag = AINFO%flag + 16

!  Assign additional workspace

     IF ( ALLOCATED( FACTORS%R ) ) THEN
       IF ( SIZE( FACTORS%R ) /= FACTORS%n ) THEN
         DEALLOCATE( FACTORS%R, stat = stat )
         IF ( stat /= 0 ) GO TO 100
         ALLOCATE( FACTORS%R( FACTORS%n ), stat = stat )
         IF ( stat /= 0 ) GO TO 100
       END IF
     ELSE
       ALLOCATE( FACTORS%R( FACTORS%n ), stat = stat )
       IF ( stat /= 0 ) GO TO 100
     END IF

     FACTORS%got_factors = .TRUE.
     RETURN

!  Error returns

 100 CONTINUE
     AINFO%flag = - 4
     AINFO%stat = stat
     IF ( CONTROL%ldiag > 0 .AND. CONTROL%lp >= 0 )                           &
       WRITE( CONTROL%lp, '( /, A, I3, /, A, I12 )' )                         & 
        ' Error return from ULS_ANALYSE with AINFO%flag = ', AINFO%flag,      & 
        ' ALLOCATE or DEALLOCATE failed with STAT=', stat

     FACTORS%got_factors = .FALSE.
     RETURN

!  End of ULS_analyse

     CONTAINS

!-*-*-*-*-  U L S _ M C 2 2 _ T H R E A D S A F E   S U B R O U T I N E  -*-*-*-

       SUBROUTINE MC22_threadsafe( n, ICN, A, nz, LENROW, IP, IQ, IW1,        &
                                   IW2, IW11 )

!  A threadsafe version of Iain Duff's HSL package MC22 to permute the
!  rows and columns of a sparse, unsymmetric matrix.

!  *******************************************************************
!  COPYRIGHT (c) 2006 Hyprotech UK and CCLRC 
!  All rights reserved.
! 
!  None of the comments in this Copyright notice between the lines
!  of asterisks shall be removed or altered in any way.
! 
!  This Package is intended for compilation without modification,
!  so most of the embedded comments have been removed.
! 
!  ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
!  Licence, see http://hsl.rl.ac.uk/archive/cou.html
! 
!  Please note that for an HSL ARCHIVE Licence:
! 
!  1. The Package must not be copied for use by any other person.
!     Supply of any part of the library by the Licensee to a third party
!     shall be subject to prior written agreement between AEA
!     Hyprotech UK Limited and the Licensee on suitable terms and
!     conditions, which will include financial conditions.
!  2. All information on the Package is provided to the Licensee on the
!     understanding that the details thereof are confidential.
!  3. All publications issued by the Licensee that include results obtained
!     with the help of one or more of the Packages shall acknowledge the
!     use of the Packages. The Licensee will notify the Numerical Analysis
!     Group at Rutherford Appleton Laboratory of any such publication.
!  4. The Packages may be modified by or on behalf of the Licensee
!     for such use in research applications but at no time shall such
!     Packages or modifications thereof become the property of the
!     Licensee. The Licensee shall make available free of charge to the
!     copyright holder for any purpose all information relating to
!     any modification.
!  5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
!     direct or consequential loss or damage whatsoever arising out of
!     the use of Packages by the Licensee.
!  *******************************************************************

!  Dummy arguments

       INTEGER, INTENT( IN )  :: n, nz
       INTEGER, INTENT( INOUT ), DIMENSION( nz ) :: ICN
       INTEGER, INTENT( INOUT ), DIMENSION( n ) :: LENROW
       INTEGER, INTENT( IN ), DIMENSION( n ) :: IP, IQ
       INTEGER, INTENT( OUT ), DIMENSION( n ) :: IW1, IW2
       INTEGER, INTENT( OUT ), DIMENSION( nz ) :: IW11
       REAL( KIND = wp ), INTENT( INOUT ), DIMENSION( nz ):: A

!  Local variables

       INTEGER :: i, ichain, iold, ipos, j2, jj, jnum, jval, length, newpos
       REAL( KIND = wp ) :: aval

       IF ( nz <= 0 .OR. n <= 0 ) RETURN

!  Set start of row i in IW1( i ) and LENROW( i ) in IW2( i )

       IW1( 1 ) = 1
       IW2( 1 ) = LENROW( 1 )
       DO i = 2, n
         IW1( i ) = IW1( i - 1 ) + LENROW( i - 1 )
         IW2( i ) = LENROW( i )
       END DO

!  Permute LENROW according to IP.  Set off-sets for new position of row iold 
!  in IW1( iold ) and put old row indices in IW11 in positions corresponding to 
!  the new position of this row in A/ICN.

       jj = 1
       DO i = 1, n
         iold = ABS( IP( i ) )
         length = IW2( iold )
         LENROW( i ) = length
         IF ( length == 0 ) CYCLE
         IW1( iold ) = IW1( iold ) - jj
         j2 = jj + length - 1
         IW11( jj : j2 ) = iold
         jj = j2 + 1
      END DO

!  Set inverse permutation to IQ in IW2

       DO i = 1, n
         IW2( ABS( IQ( i ) ) ) = i
       END DO

!  Permute A and ICN in place, changing to new column numbers.

!  Each pass through the main loop places a closed chain of column indices
!  in their new (and final) positions. This is recorded by setting the IW11 
!  entry to zero so that any which are subsequently encountered during this 
!  major scan can be bypassed

       DO i = 1, nz
         iold = IW11( i )
         IF ( iold == 0 ) CYCLE
         ipos = i
         jval = ICN( i )

!  If row iold is in same position after permutation leave it alone

         IF ( IW1( iold ) /= 0 ) THEN
           aval = A( i )

!  Each pass through this loop places one (permuted) column index in its 
!  final position  .. viz. ipos.

           DO ichain = 1, nz

!  newpos is the original position in A/ICN of the element to be placed in 
!  position ipos.  It is also the position of the next element in the chain

             newpos = ipos + IW1( iold )

!  Is the chain complete?

             IF ( newpos == i ) EXIT
             A( ipos ) = A( newpos )
             jnum = ICN( newpos )
             ICN( ipos ) = IW2( jnum )
             ipos = newpos
             iold = IW11( ipos )
             IW11( ipos ) = 0
           END DO
           A( ipos ) = aval
         END IF
         ICN( ipos ) = IW2( jval )
       END DO

       RETURN

!  End of subroutine MC22_threadsafe

       END SUBROUTINE MC22_threadsafe

!-*-*-*-*-  U L S _ M C 2 3 _ T H R E A D S A F E   S U B R O U T I N E  -*-*-*-

       SUBROUTINE MC23_threadsafe( n, ICN, A, licn, LENR, IDISP, IP, IQ,       &
                                   LENOFF, IW1, IW2, IW3, IW4, IW11, IW12,     &
                                   ICNTL, INFO )

!  A threadsafe version of Iain Duff's HSL package MC23 to permute a given
!  sparse, unsymmetric matrix.to block-triangular form

!  *******************************************************************
!  COPYRIGHT (c) 2006 Hyprotech UK and CCLRC 
!  All rights reserved.
! 
!  None of the comments in this Copyright notice between the lines
!  of asterisks shall be removed or altered in any way.
! 
!  This Package is intended for compilation without modification,
!  so most of the embedded comments have been removed.
! 
!  ALL USE IS SUBJECT TO LICENCE. For full details of an HSL ARCHIVE
!  Licence, see http://hsl.rl.ac.uk/archive/cou.html
! 
!  Please note that for an HSL ARCHIVE Licence:
! 
!  1. The Package must not be copied for use by any other person.
!     Supply of any part of the library by the Licensee to a third party
!     shall be subject to prior written agreement between AEA
!     Hyprotech UK Limited and the Licensee on suitable terms and
!     conditions, which will include financial conditions.
!  2. All information on the Package is provided to the Licensee on the
!     understanding that the details thereof are confidential.
!  3. All publications issued by the Licensee that include results obtained
!     with the help of one or more of the Packages shall acknowledge the
!     use of the Packages. The Licensee will notify the Numerical Analysis
!     Group at Rutherford Appleton Laboratory of any such publication.
!  4. The Packages may be modified by or on behalf of the Licensee
!     for such use in research applications but at no time shall such
!     Packages or modifications thereof become the property of the
!     Licensee. The Licensee shall make available free of charge to the
!     copyright holder for any purpose all information relating to
!     any modification.
!  5. Neither CCLRC nor Hyprotech UK Limited shall be liable for any
!     direct or consequential loss or damage whatsoever arising out of
!     the use of Packages by the Licensee.
!  *******************************************************************

!  Dummy arguments

       INTEGER, INTENT( IN ) :: licn, n
       INTEGER, INTENT( OUT ), DIMENSION( 2 ) :: IDISP
       INTEGER, INTENT( INOUT ), DIMENSION( 10 ) :: INFO
       INTEGER, INTENT( INOUT ), DIMENSION( 10 ) :: ICNTL
       INTEGER, INTENT( OUT ), DIMENSION( n ) :: IP, IQ, LENOFF, IW11, IW12
       INTEGER, INTENT( OUT ), DIMENSION( n ) :: IW1, IW2, IW3, IW4
       INTEGER, INTENT( INOUT ), DIMENSION( n ) :: LENR
       INTEGER, INTENT( INOUT ), DIMENSION( licn ) :: ICN
       REAL( KIND = wp ), INTENT( INOUT ), DIMENSION( licn ):: A

!  Local variables

       INTEGER :: i, i1, i2, ibeg, iblock, iend, ii, ilend, inew, iold, irowb, &
                  irowe, j, jj, jnew, jnpos, jold, k, leni, nz, large, lp,     & 
                  num, numnz
       LOGICAL :: abort
       EXTERNAL :: MC13ED, MC21BD

       INFO( 1 ) = 0
       INFO( 2 ) = 0
       num = 0
       large = 0
       lp = ICNTL( 1 )
       abort = ICNTL( 2 ) == 1

!  Set up pointers IW1 to the beginning of the rows and set LENOFF to LENR

       LENOFF = LENR
       IW11( 1 ) = 1
       DO i = 2, n
         IW11( i ) = IW11( i - 1 ) + LENR( i - 1 )
       END DO

!  IDISP( 1 ) points to the first position in A/ICN after the off-diagonal 
!  blocks and untreated rows

       IDISP( 1 ) = IW11( N ) + LENR( n )

!  Find the row permutation IP to make the diagonal zero-free

       CALL MC21BD( n, ICN, licn, IW11, LENR, IP, numnz, IW1, IW2, IW3, IW4 )

!  Possible error return for structurally singular matrices

       IF ( numnz /= n .AND. abort ) THEN
         IF ( lp > 0 ) WRITE( lp, "( /, A, /, 10X, A, I0 )" )                  &
           ' Error return from MC23_threadsafe because',                       &
           ' matrix is structurally singular,  rank = ',  numnz
         IDISP( 1 ) = - 1
         INFO( 1 ) = - 1
         INFO( 2 ) = numnz
         GO TO 900
       END IF

!  IW12 and LENR are permutations of IW11 and LENR/LENOFF suitable for entry
!  to MC13ED since matrix with these row pointer and length arrays has maximum 
!  number of non-zeros on the diagonal

       DO ii = 1, n
         i = IP( ii )
         IW12( ii ) = IW11( i )
         LENR( ii ) = LENOFF( i )
       END DO

!  Find symmetric permutation IQ to block lower triangular form - there
!  are num blocks

       CALL MC13ED( n, ICN, licn, IW12, LENR, IQ, IW4, num, IW1, IW2, IW3 )

!  Move the whole matrix to the end of the storage if it irreducible

       IF ( num == 1 ) THEN
         DO i = 1, n
           LENR( i ) = LENOFF( i )
           IP( i ) = i
           IQ( i ) = i
         END DO
         LENOFF( 1 ) = - 1

!  IDISP( 1 ) is the first position after the last element in the off-diagonal 
!  blocks and untreated rows

         nz = IDISP( 1 ) - 1
         IDISP( 1 ) = 1

!  IDISP( 2 ) is the position in A/ICN of the first element in the diagonal block

         IDISP( 2 ) = licn - nz + 1
         large = n
         IF ( nz == licn ) GO TO 900
         DO k = 1, nz
           j = nz - k + 1
           jj = licn - k + 1
           A( jj ) = A( j )
           ICN( jj ) = ICN( j )
         END DO
         GO TO 900
       END IF

!  Reorder the data structure: composite row permutation IP( i ) = IP( IQ( i ) )

       DO ii = 1, n
         i = IQ( ii )
         IW1( ii ) = IP( i )
       END DO

       DO i = 1, n
         IP( i ) = IW1( i )
       END DO

!  Run through blocks in reverse order separating diagonal blocks which are moved
!  to the end of the storage. Elements in off-diagonal blocks are left in place 
!  unless a compress is necessary

!  ibeg indicates the smallest value of j for which ICN( j ) has been set to zero
!  when element in position j was moved to the diagonal block part of storage

       ibeg = licn + 1

!  iend is the position of the first element of rows which are in diagonal blocks

       iend = licn + 1

!  large is the dimension of the largest block encountered so far

       large = 0
       DO k = 1, num
         iblock = num - k + 1

!  i1 and i2 are the first and last row (in permuted form) of block iblock

         i1 = IW4( iblock )
         i2 = n
         IF ( k /= 1 ) i2 = IW4( iblock + 1 ) - 1
         large = MAX( large, i2 - i1 + 1 )

!  Go through the rows of block iblock in reverse order

         DO ii = i1, i2
           inew = i2 - ii + i1

!  Deal with row inew in permuted form - row iold in original matrix.

           iold = IP( inew )

!  Check if there is space to move up diagonal block portion of row

           IF ( iend - IDISP( 1 ) < LENOFF( iold ) ) THEN

!  In-place compress moves separated off-diagonal elements and untreated rows 
!  to the front of storage

             jnpos = ibeg
             ilend = IDISP( 1 ) - 1
             IF ( ilend < ibeg ) GO TO 800
             DO j = ibeg, ilend
               IF ( ICN( j ) /= 0 ) THEN
                 ICN( jnpos ) = ICN( j )
                 A( jnpos ) = A( j )
                 jnpos = jnpos + 1
               END IF
             END DO
             IDISP( 1 ) = jnpos
             IF ( iend - jnpos < LENOFF( iold ) ) GO TO 800
             ibeg = licn + 1

!  Reset pointers to the beginning of the rows

             DO i = 2, n
               IW11( I ) = IW11( I - 1 ) + LENOFF( I-1 )
             END DO

!  Row iold is now split into diagonal and off-diagonal parts

           END IF
           irowb = IW11( iold )
           leni = 0
           irowe = irowb + LENOFF( iold ) - 1

!  Backward scan of whole of row iold in the original matrix.

           IF ( irowe >= irowb ) THEN
             DO jj = irowb, irowe
               j = irowe - jj + irowb
               jold = ICN( j )

!  IW2 holds the inverse permutation to IQ as set by MC13ED

               jnew = IW2( jold )

!  If jnew < i1 then element is in off-diagonal block and so is left in place

               IF ( jnew >= i1 ) THEN

!  Element is in diagonal block and is moved to the end of the storage

                 iend = iend - 1
                 A( iend ) = A( j )
                 ICN( iend ) = jnew
                 ibeg = MIN( ibeg, j )
                 ICN( j ) = 0
                 leni = leni + 1
               END IF
             END DO

             LENOFF( iold ) = LENOFF( iold ) - leni
           END IF
           LENR( INEW ) = LENI
         END DO

         IP( I2 ) = -IP( I2 )
       END DO

!  Reset IP( n ) to a positive value

       IP( n ) = - IP( n )

!  IDISP( 2 ) is the position of the first element in the diagonal blocks
 
       IDISP( 2 ) = iend
 
!  This compress moves all off-diagonal elements to the front of the array

       IF ( ibeg > licn ) GO TO 900
       jnpos = ibeg
       ilend = IDISP( 1 ) - 1
       DO j = ibeg, ilend
         IF ( ICN( j ) /= 0 ) THEN
           ICN( jnpos ) = ICN( j )
           A( jnpos ) = A( j )
           jnpos = jnpos + 1
         END IF
       END DO

!  IDISP( 1 ) is first position after last element of off-diagonal blocks.

       IDISP( 1 ) = jnpos
       GO TO 900

!  Error return

  800  CONTINUE
       IF ( lp > 0 ) WRITE( lp, "( /, A, /, 10X, A, I0 )"  )                   &
         ' Error return from MC23_threadsafe because',                         &
         ' licn not big enough increase by ', n

       IDISP( 1 ) = - 2
       INFO( 1 ) = - 2
       INFO( 2 ) = licn

!  This is the only return point

  900  CONTINUE
       INFO( 3 ) = numnz
       INFO( 4 ) = num
       INFO( 5 ) = large

       RETURN

!  End of subroutine MC23_threadsafe

       END SUBROUTINE MC23_threadsafe

     END SUBROUTINE ULS_analyse

!-*-*-*-*-*-*-*-*-*-   U L S _ S O L V E  S U B R O U T I N E   -*-*-*-*-*-*-*

     SUBROUTINE ULS_solve( MATRIX, FACTORS, RHS, X, CONTROL, SINFO, trans  )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Solve the linear system using the factors obtained in the factorization
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

     TYPE ( SMT_type ), INTENT( IN ) :: MATRIX
     TYPE ( ULS_factors ), INTENT( INOUT ) :: FACTORS
     REAL ( KIND = wp ), INTENT( IN ), DIMENSION( MATRIX%m ) :: RHS
     REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( MATRIX%n ) :: X
     TYPE ( ULS_control ), INTENT( IN ) :: CONTROL
     TYPE ( ULS_sinfo ), INTENT( OUT ) :: sinfo
     INTEGER, OPTIONAL, INTENT( IN ) :: trans

!  Local variables

     INTEGER :: i, j, l, iter, mtype
     REAL ( KIND = wp ) :: resid, resid_init, resid_old

     SINFO%more = 0
     SINFO%stat = 0

! Check whether this has been preceded by a call to ULS_factorize

     IF ( .NOT. FACTORS%got_factors ) THEN
       SINFO%flag = - 10
       WRITE( CONTROL%lp, '( /, A, I3, /, A, I12 )' )                          &
         ' Error return from ULS_SOLVE with sinfo%flag = ', SINFO%flag,        &
         ' No prior call to ULS_FACTORIZE'
       RETURN
     END IF

     mtype = 1
     IF ( PRESENT( trans ) ) mtype = 0
     DO iter = 0, CONTROL%maxit

!  Compute the residuals ...

       IF ( iter > 0 ) THEN
         FACTORS%R( : MATRIX%m ) = RHS ; FACTORS%R( MATRIX%m + 1 : ) = zero

!   ... b - A x

         IF ( mtype == 1 ) THEN
           DO l = 1, MATRIX%ne
             i = MATRIX%row( l )
             FACTORS%R( i )                                                    &
               = FACTORS%R( i ) - MATRIX%val( l ) * X(  MATRIX%col( l ) )
           END DO

!   ... b - A^T x

         ELSE
           DO l = 1, MATRIX%ne
             j = MATRIX%col( l )
             FACTORS%R( j )                                                    &
               = FACTORS%R( j ) - MATRIX%val( l ) * X(  MATRIX%row( l ) )
           END DO
         END IF

         resid_old = resid
         resid = MAXVAL( ABS( FACTORS%R ) )
!        WRITE( 6, "( ' residual ', ES12.4 )" ) resid

!  Check for convergence

         IF ( resid <= eps * resid_init ) EXIT

!  Check for stagnation

         IF ( resid > CONTROL%cgce * resid_old ) THEN
           SINFO%flag = - 8
           IF ( CONTROL%lp >= 0 )                                              &
             WRITE( CONTROL%lp, '( /, A, I3, /, A, I12 )' )                    &
             ' Error return from ULS_SOLVE with sinfo%flag = ', SINFO%flag,    &
             ' Iterative refinement diverging'
           RETURN
         END IF

!  Compute the solution to the corrector equation A x = r (or A^T x = r)
!  overwriting r with x

         CALL MA33CD( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,         &
                      FACTORS%LENR, FACTORS%LENRL, FACTORS%LENOFF,             &
                      FACTORS%IDISP, FACTORS%IP, FACTORS%IQ, FACTORS%R,        &
                      FACTORS%W, mtype, FACTORS%RINFO )

!  Update the estimate of the solution

         IF ( MATRIX%m <= MATRIX%n ) THEN
           X = X + FACTORS%R
         ELSE
           X = X + FACTORS%R( : MATRIX%n )
         END IF

!  Special case for first iteration

       ELSE
         resid_init = MAXVAL( ABS( RHS ) )
         resid = resid_init
!        WRITE( 6, "( ' residual ', ES12.4 )" ) resid_init
         IF ( MATRIX%m <= MATRIX%n ) THEN
           X( : MATRIX%m ) = RHS ; X( MATRIX%m + 1 : ) = zero
           CALL MA33CD( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,      &
                        FACTORS%LENR, FACTORS%LENRL, FACTORS%LENOFF,          &
                        FACTORS%IDISP, FACTORS%IP, FACTORS%IQ, X,             &
                        FACTORS%W, mtype, FACTORS%RINFO )
         ELSE
           FACTORS%R = RHS
           CALL MA33CD( FACTORS%n, FACTORS%ICN, FACTORS%A, FACTORS%licn,      &
                        FACTORS%LENR, FACTORS%LENRL, FACTORS%LENOFF,          &
                        FACTORS%IDISP, FACTORS%IP, FACTORS%IQ, FACTORS%R,     &
                        FACTORS%W, mtype, FACTORS%RINFO )
           X = FACTORS%R( : MATRIX%n )
         END IF  
       END IF
     END DO

     SINFO%flag = 0

     RETURN

!  End of ULS_solve

     END SUBROUTINE ULS_solve

!-*-*-*-*-*-*-*-   U L S _ F I N A L I Z E  S U B R O U T I N E   -*-*-*-*-*-

     SUBROUTINE ULS_finalize( FACTORS, CONTROL, info )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  DEALLOCATE all currently ALLOCATEd arrays
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

     TYPE ( ULS_factors ), INTENT( INOUT ) :: FACTORS
     TYPE ( ULS_control ), INTENT( IN ) :: CONTROL
     INTEGER, INTENT( OUT ) :: info

!  Local variables

     INTEGER :: alloc_stat

     info = 0
     FACTORS%got_factors = .FALSE.

     IF ( ALLOCATED( FACTORS%ICN ) ) THEN
       DEALLOCATE( FACTORS%ICN, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IFIRST ) ) THEN
       DEALLOCATE( FACTORS%IFIRST, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IP ) )     THEN
       DEALLOCATE( FACTORS%IP, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IPC ) )    THEN
       DEALLOCATE( FACTORS%IPC, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IPTR ) )   THEN
       DEALLOCATE( FACTORS%IPTR, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IQ ) )     THEN
       DEALLOCATE( FACTORS%IQ, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%IRN ) )    THEN
       DEALLOCATE( FACTORS%IRN, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LASTC ) )  THEN
       DEALLOCATE( FACTORS%LASTC, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LASTR ) )  THEN
       DEALLOCATE( FACTORS%LASTR, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LENC ) )   THEN
       DEALLOCATE( FACTORS%LENC, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LENR ) )   THEN
       DEALLOCATE( FACTORS%LENR, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LENRL ) )  THEN
       DEALLOCATE( FACTORS%LENRL, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%NEXTC ) )  THEN
       DEALLOCATE( FACTORS%NEXTC, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%NEXTR ) )  THEN
       DEALLOCATE( FACTORS%NEXTR, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%A ) )      THEN
       DEALLOCATE( FACTORS%A, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%LENOFF ) ) THEN
       DEALLOCATE( FACTORS%LENOFF, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%W ) )      THEN
       DEALLOCATE( FACTORS%W, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF
     IF ( ALLOCATED( FACTORS%R ) )      THEN
       DEALLOCATE( FACTORS%R, stat = alloc_stat )
       IF ( alloc_stat /= 0 ) info = alloc_stat
     END IF

     IF ( info /= 0 .AND.CONTROL%ldiag > 0 .and. CONTROL%lp >= 0  )            &
       WRITE( CONTROL%lp, '( /, 2A, I0 )' ) ' Error return from ULS_finalize:',&
        ' DEALLOCATE failed with STAT=', info
   
     RETURN

!  End of ULS_finalize

     END SUBROUTINE ULS_finalize

!-*-  U L S _ S P E C I A L _ R O W S _ A N D _ C O L S  S U B R O U T I N E  -*-

     SUBROUTINE ULS_special_rows_and_cols( FACTORS, rank, ROWS, COLS, info )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!  Identify rows and columns of A taken into account when solving systems
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

     TYPE ( ULS_FACTORS ), INTENT( IN ) :: FACTORS
     INTEGER,INTENT( OUT ) :: rank, info
     INTEGER,INTENT( OUT ), DIMENSION( FACTORS%orig_m ) :: ROWS
     INTEGER,INTENT( OUT ), DIMENSION( FACTORS%orig_n ) :: COLS

!  Local variables

     INTEGER :: i, j

     rank = FACTORS%INFO ( 3 )
     rank = 0
     DO i = 1, FACTORS%n
       j = FACTORS%IQ( i )
       IF ( j > 0 ) THEN
         rank = rank + 1
         ROWS( rank ) = ABS( FACTORS%IP( i ) )
         COLS( rank ) = j
       END IF
     END DO
     info = 0

     RETURN

!  End of ULS_special_rows_and_cols

     END SUBROUTINE ULS_special_rows_and_cols

   END MODULE GALAHAD_ULS_double

!  IKEEP -> LENR
!  IKEEP( : 2 ) -> IP
!  IKEEP( : 3 ) -> IQ
!  IKEEP( : 4 ) -> LENRL
!  IKEEP( : 5 ) -> LENOFF
!  IW -> IPC
!  IW( : 2 ) -> LENC
!  IW( : 3 ) -> IFIRST
!  IW( : 4 ) -> LASTR
!  IW( : 5 ) -> NEXTR
!  IW( : 6 ) -> LASTC
!  IW( : 7 ) -> NEXTC
!  IW( : 8 ) -> IPTR

