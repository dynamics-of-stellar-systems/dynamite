! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-*-  G A L A H A D _ Q P D  M O D U L E  -*-*-*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.0. August 10th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

!     -----------------------------------------------------
!     | Provides a generic derived type to hold and share |
!     | private data between the GALAHAD QP packages      |
!     |      NOT INTENDED FOR PUBLIC CONSUMPTION          |
!     -----------------------------------------------------

   MODULE GALAHAD_QPD_double

     USE GALAHAD_RAND_double, ONLY : RAND_seed
     USE GALAHAD_SILS_double, ONLY : SILS_factors, SILS_control,               &
                                     SILS_ainfo, SILS_finfo, SMT_type
     USE GALAHAD_SBLS_double, ONLY: SBLS_data_type
     USE GALAHAD_SCU_double, ONLY : SCU_matrix_type, SCU_info_type, SCU_data_type
     USE GALAHAD_QPP_double, QPD_dims_type => QPP_dims_type
   
     IMPLICIT NONE

     PRIVATE
     PUBLIC :: QPD_HX, QPD_AX

!--------------------
!   P r e c i s i o n
!--------------------

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  ==============================
!  The QPD_data_type derived type
!  ==============================

     TYPE, PUBLIC :: QPD_data_type

! -----------------
!  Scalar componets
! -----------------

!  Common scalar components

       INTEGER :: start_print, stop_print
       LOGICAL :: new_problem_structure

!  QPA scalar components

       INTEGER :: prec_hist
       LOGICAL :: auto_prec, auto_fact

!  QPB/LSQP scalar components

       INTEGER :: trans, remove_more_deps
       LOGICAL :: tried_to_remove_deps, save_structure

! -----------------------
!  Allocatable components
! -----------------------

!  Common allocatable components

       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Abycol_row
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Abycol_ptr
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Abycol_val
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RHS
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H_s
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: A_s

!  QPA & QPB allocatable components

       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: S
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRAD

!  QPA & LSQP allocatable components

       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IBREAK
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: BREAKP

!  QPB & LSQP allocatable components

       INTEGER, ALLOCATABLE, DIMENSION( : ) :: IW
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: K_colptr
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Index_C_freed
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Index_C_more_freed
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES_x
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: SOL_y
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES_y
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: BEST_y
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: SOL
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: BEST
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: HX
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRAD_L
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_trial
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_X_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_X_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: BARRIER_X
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DY_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_C_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DY_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_C_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: BARRIER_C
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: SCALE_C
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DELTA
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DZ_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DZ_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X0
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIAG_X
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIAG_C
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_freed
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_more_freed
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y_last
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z_last

!  LPB allocatable components

       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_last

!  QPA allocatable components

       INTEGER, ALLOCATABLE, DIMENSION( : ) :: SC
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: REF
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_up_or_low
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_up_or_low 
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: PERM
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: S_row
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: S_col
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: S_colptr
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: A_norms
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERT
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: VECTOR
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: B
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: S_perm
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DX
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: R_pcg
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_pcg
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: P_pcg
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: S_val
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: RES_print
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: DIAG

!  QPB allocatable components

       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Index_X_fixed
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: Index_C_fixed
       INTEGER, ALLOCATABLE, DIMENSION( : ) :: H_band_ptr
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_fixed
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_fixed
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRAD_X_phi
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: GRAD_C_phi

!  LSQP and WPC allocatable components

       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: COEF0
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: COEF1
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: COEF2
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: COEF3
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: COEF4
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DELTA_cor
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DZ_cor_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DZ_cor_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DY_cor_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DY_cor_u

!  WPC allocatable components

       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_X_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_X_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_Y_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_Y_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_Z_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_Z_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_C_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: PERTURB_C_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_Y_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_Y_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_Z_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DIST_Z_u
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MU_X_l
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MU_X_U
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MU_C_L
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MU_C_U
       REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: MU

! -----------------------
!  Derived type componets
! -----------------------

!  Common derived type components

       TYPE ( QPD_dims_type ) :: dims
       TYPE ( SMT_type ) :: K, H_sub, A_sub, C_sub
       TYPE ( SILS_factors ) :: FACTORS
       TYPE ( SILS_control ) :: CNTL
       TYPE ( QPP_control_type ) :: QPP_control
       TYPE ( QPP_inform_type ) :: QPP_inform
       TYPE ( QPP_map_type ) :: QPP_map

!  LSQP derived type components

       TYPE ( QPP_map_type ) :: QPP_map_fixed, QPP_map_freed
       TYPE ( QPP_map_type ) :: QPP_map_more_freed
       TYPE ( QPD_dims_type ) :: dims_save_fixed, dims_save_freed
       TYPE ( QPD_dims_type ) :: dims_save_more_freed

!  QPA derived type components

       TYPE ( RAND_seed ) :: seed
       TYPE ( SILS_control ) :: CNTLA
       TYPE ( SILS_ainfo ) :: AINFO
       TYPE ( SILS_finfo ) :: FINFO
       TYPE ( SCU_matrix_type ) :: SCU_mat
       TYPE ( SCU_info_type ) :: SCU_info
       TYPE ( SCU_data_type ) :: SCU_data

!  SBLS derived type components

       TYPE ( SBLS_data_type ) :: SBLS_data

     END TYPE

! -------------------------------------------
!  Subroutines shared between the QP packages
! -------------------------------------------

   CONTAINS

!-*-*-*-*-*-*-*-*-*-*-   Q P D _ H X  S U B R O U T I N E  -*-*-*-*-*-*-*-*-*-

      SUBROUTINE QPD_HX( dims, n, R, H_ne, H_val, H_col, H_ptr, X, op,          &
                         semibw, H_band_ptr )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ............................................
!      .                                          .
!      .  Perform the operation r := r op H * x   .
!         where op is + or -                      .
!      ............................................

!  Arguments:
!  =========
!
!   dims    see module GALAHAD_QPP
!   H_*     sparse storage by rows or band
!   X       the vector x
!   R       the result of adding H * x to r
!   semibw  if present, only those entries within a band of semi-bandwidth
!           semibw will be accessed
!   op      character string "+" or "-"

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      TYPE ( QPD_dims_type ), INTENT( IN ) :: dims
      INTEGER, INTENT( IN ) :: n, H_ne
      INTEGER, OPTIONAL, INTENT( IN ) :: semibw
      CHARACTER( LEN = 1 ), INTENT( IN ) :: op
      INTEGER, INTENT( IN ), DIMENSION( n + 1 ) :: H_ptr
      INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( n ) :: H_band_ptr
      INTEGER, INTENT( IN ), DIMENSION( H_ne ) ::  H_col
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( H_ne ) :: H_val
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: R

!  Local variables

      INTEGER :: hd_start, hd_end, hnd_start, hnd_end, i, j, l, type
      REAL ( KIND = wp ) :: xi, ri

!  For a banded portion of H

      IF ( PRESENT( semibw ) ) THEN

        IF ( op( 1 : 1 ) == '+' ) THEN
  
!  r <- r + H * x (commented out since it is not used at present)
  
!         DO type = 1, 6
    
!           SELECT CASE( type )
!           CASE ( 1 )
    
!             hd_start  = 1
!             hd_end    = dims%h_diag_end_free
!             hnd_start = hd_end + 1
!             hnd_end   = dims%x_free
    
!           CASE ( 2 )
    
!             hd_start  = dims%x_free + 1
!             hd_end    = dims%h_diag_end_nonneg
!             hnd_start = hd_end + 1
!             hnd_end   = dims%x_l_start - 1
    
!           CASE ( 3 )
    
!             hd_start  = dims%x_l_start
!             hd_end    = dims%h_diag_end_lower
!             hnd_start = hd_end + 1
!             hnd_end   = dims%x_u_start - 1
    
!           CASE ( 4 )
    
!             hd_start  = dims%x_u_start
!             hd_end    = dims%h_diag_end_range
!             hnd_start = hd_end + 1
!             hnd_end   = dims%x_l_end
    
!           CASE ( 5 )
    
!             hd_start  = dims%x_l_end + 1
!             hd_end    = dims%h_diag_end_upper
!             hnd_start = hd_end + 1
!             hnd_end   = dims%x_u_end
    
!           CASE ( 6 )
    
!             hd_start  = dims%x_u_end + 1
!             hd_end    = dims%h_diag_end_nonpos
!             hnd_start = hd_end + 1
!             hnd_end   = n
    
!           END SELECT
    
!  rows with a diagonal entry
    
!           hd_end = MIN( hd_end, n )
!           DO i = hd_start, hd_end
!             DO l = H_band_ptr( i ), H_ptr( i + 1 ) - 2
!               j = H_col( l )
!               R( j ) = R( j ) + H_val( l ) * X( i )
!               R( i ) = R( i ) + H_val( l ) * X( j )
!             END DO
!             R( i ) = R( i ) + H_val( H_ptr( i + 1 ) - 1 ) * X( i )
!           END DO
!           IF ( hd_end == n ) EXIT
    
!  rows without a diagonal entry
    
!           hnd_end = MIN( hnd_end, n )
!           DO i = hnd_start, hnd_end
!             DO l = H_band_ptr( i ), H_ptr( i + 1 ) - 1
!               j = H_col( l )
!               R( j ) = R( j ) + H_val( l ) * X( i )
!               R( i ) = R( i ) + H_val( l ) * X( j )
!             END DO
!           END DO
!           IF ( hnd_end == n ) EXIT
    
!         END DO
        ELSE
  
!  r <- r - H * x

          DO type = 1, 6
    
            SELECT CASE( type )
            CASE ( 1 )
    
              hd_start  = 1
              hd_end    = dims%h_diag_end_free
              hnd_start = hd_end + 1
              hnd_end   = dims%x_free
    
            CASE ( 2 )
    
              hd_start  = dims%x_free + 1
              hd_end    = dims%h_diag_end_nonneg
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_start - 1
    
            CASE ( 3 )
    
              hd_start  = dims%x_l_start
              hd_end    = dims%h_diag_end_lower
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_start - 1
    
            CASE ( 4 )
    
              hd_start  = dims%x_u_start
              hd_end    = dims%h_diag_end_range
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_end
    
            CASE ( 5 )
    
              hd_start  = dims%x_l_end + 1
              hd_end    = dims%h_diag_end_upper
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_end
    
            CASE ( 6 )
    
              hd_start  = dims%x_u_end + 1
              hd_end    = dims%h_diag_end_nonpos
              hnd_start = hd_end + 1
              hnd_end   = n
    
            END SELECT
    
!  rows with a diagonal entry
    
            hd_end = MIN( hd_end, n )
            DO i = hd_start, hd_end
              xi = X( i )
              ri = R( i )
              DO l = H_band_ptr( i ), H_ptr( i + 1 ) - 2
                j = H_col( l )
                R( j ) = R( j ) - H_val( l ) * xi
                ri = ri - H_val( l ) * X( j )
              END DO
              R( i ) = ri - H_val( H_ptr( i + 1 ) - 1 ) * xi
            END DO
            IF ( hd_end == n ) EXIT
    
!  rows without a diagonal entry
    
            hnd_end = MIN( hnd_end, n )
            DO i = hnd_start, hnd_end
              xi = X( i )
              ri = R( i )
              DO l = H_band_ptr( i ), H_ptr( i + 1 ) - 1
                j = H_col( l )
                R( j ) = R( j ) - H_val( l ) * xi
                ri = ri - H_val( l ) * X( j )
              END DO
              R( i ) = ri
            END DO
            IF ( hnd_end == n ) EXIT
    
          END DO
        END IF

!  For the whole of H

      ELSE
        IF ( op( 1 : 1 ) == '+' ) THEN
  
!  r <- r + H * x
  
          DO type = 1, 6
    
            SELECT CASE( type )
            CASE ( 1 )
    
              hd_start  = 1
              hd_end    = dims%h_diag_end_free
              hnd_start = hd_end + 1
              hnd_end   = dims%x_free
    
            CASE ( 2 )
    
              hd_start  = dims%x_free + 1
              hd_end    = dims%h_diag_end_nonneg
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_start - 1
    
            CASE ( 3 )
    
              hd_start  = dims%x_l_start
              hd_end    = dims%h_diag_end_lower
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_start - 1
    
            CASE ( 4 )
    
              hd_start  = dims%x_u_start
              hd_end    = dims%h_diag_end_range
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_end
    
            CASE ( 5 )
    
              hd_start  = dims%x_l_end + 1
              hd_end    = dims%h_diag_end_upper
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_end
    
            CASE ( 6 )
    
              hd_start  = dims%x_u_end + 1
              hd_end    = dims%h_diag_end_nonpos
              hnd_start = hd_end + 1
              hnd_end   = n
    
            END SELECT
    
!  rows with a diagonal entry
    
            hd_end = MIN( hd_end, n )
            DO i = hd_start, hd_end
              xi = X( i )
              ri = R( i )
              DO l = H_ptr( i ), H_ptr( i + 1 ) - 2
                j = H_col( l )
                R( j ) = R( j ) + H_val( l ) * xi
                ri = ri + H_val( l ) * X( j )
              END DO
              R( i ) = ri + H_val( H_ptr( i + 1 ) - 1 ) * xi
            END DO
            IF ( hd_end == n ) EXIT
    
!  rows without a diagonal entry
    
            hnd_end = MIN( hnd_end, n )
            DO i = hnd_start, hnd_end
              xi = X( i )
              ri = R( i )
              DO l = H_ptr( i ), H_ptr( i + 1 ) - 1
                j = H_col( l )
                R( j ) = R( j ) + H_val( l ) * xi
                ri = ri + H_val( l ) * X( j )
              END DO
              R( i ) = ri
            END DO
            IF ( hnd_end == n ) EXIT
    
          END DO
        ELSE
  
!  r <- r - H * x
  
          DO type = 1, 6
    
            SELECT CASE( type )
            CASE ( 1 )
    
              hd_start  = 1
              hd_end    = dims%h_diag_end_free
              hnd_start = hd_end + 1
              hnd_end   = dims%x_free
    
            CASE ( 2 )
    
              hd_start  = dims%x_free + 1
              hd_end    = dims%h_diag_end_nonneg
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_start - 1
    
            CASE ( 3 )
    
              hd_start  = dims%x_l_start
              hd_end    = dims%h_diag_end_lower
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_start - 1
    
            CASE ( 4 )
    
              hd_start  = dims%x_u_start
              hd_end    = dims%h_diag_end_range
              hnd_start = hd_end + 1
              hnd_end   = dims%x_l_end
    
            CASE ( 5 )
    
              hd_start  = dims%x_l_end + 1
              hd_end    = dims%h_diag_end_upper
              hnd_start = hd_end + 1
              hnd_end   = dims%x_u_end
    
            CASE ( 6 )
    
              hd_start  = dims%x_u_end + 1
              hd_end    = dims%h_diag_end_nonpos
              hnd_start = hd_end + 1
              hnd_end   = n
    
            END SELECT
    
!  rows with a diagonal entry
    
            hd_end = MIN( hd_end, n )
            DO i = hd_start, hd_end
              xi = X( i )
              ri = R( i )
              DO l = H_ptr( i ), H_ptr( i + 1 ) - 2
                j = H_col( l )
                R( j ) = R( j ) - H_val( l ) * xi
                ri = ri - H_val( l ) * X( j )
              END DO
              R( i ) = ri - H_val( H_ptr( i + 1 ) - 1 ) * xi
            END DO
            IF ( hd_end == n ) EXIT
    
!  rows without a diagonal entry
    
            hnd_end = MIN( hnd_end, n )
            DO i = hnd_start, hnd_end
              xi = X( i )
              ri = R( i )
              DO l = H_ptr( i ), H_ptr( i + 1 ) - 1
                j = H_col( l )
                R( j ) = R( j ) - H_val( l ) * xi
                ri = ri - H_val( l ) * X( j )
              END DO
              R( i ) = ri
            END DO
            IF ( hnd_end == n ) EXIT
    
          END DO
        END IF
      END IF
      RETURN

!  End of subroutine QPD_HX

      END SUBROUTINE QPD_HX

!-*-*-*-*-*-*-*-*-*-*-   Q P D _ A x  S U B R O U T I N E  -*-*-*-*-*-*-*-*-

      SUBROUTINE QPD_AX( dim_r, R, m, A_ne, A_val, A_col, A_ptr, dim_x, X, op )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Perform the operation r := r +/- A * x    .
!      .                     or r := r +/- A^T * x  .
!      .                                            .
!      ..............................................

!  Arguments:
!  =========

!   R      the result r of adding/subtracting A * x or A^T *x to/from r
!   X      the vector x
!   op     2 string character: possible values are
!          '+ '   r <- r + A * x
!          '+T'   r <- r + A^T * x
!          '- '   r <- r - A * x
!          '-T'   r <- r - A^T * x
 
! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

      INTEGER, INTENT( IN ) :: dim_x, dim_r, m, A_ne
      CHARACTER( LEN = 2 ), INTENT( IN ) :: op
      INTEGER, INTENT( IN ), DIMENSION( m + 1 ) :: A_ptr
      INTEGER, INTENT( IN ), DIMENSION( A_ne ) ::  A_col
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( dim_x ) :: X
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( A_ne ) :: A_val
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( dim_r ) :: R

!  Local variables

      INTEGER :: i, l
      REAL ( KIND = wp ) :: xi, ri

      IF ( op( 1 : 1 ) == '+' ) THEN

!  r <- r + A^T * x

        IF ( op( 2 : 2 ) == 'T' .OR. op( 2 : 2 ) == 't' ) THEN
          DO i = 1, m
            xi = X( i )
            DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
              R( A_col( l ) ) = R( A_col( l ) ) + A_val( l ) * xi
            END DO
          END DO

!  r <- r + A * x

        ELSE
          DO i = 1, m
            ri = R( i )
            DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
              ri = ri + A_val( l ) * X( A_col( l ) )
            END DO
            R( i ) = ri
          END DO
        END IF

      ELSE

!  r <- r - A^T * x

        IF ( op( 2 : 2 ) == 'T' .OR. op( 2 : 2 ) == 't' ) THEN
          DO i = 1, m
            xi = X( i )
            DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
              R( A_col( l ) ) = R( A_col( l ) ) - A_val( l ) * xi
            END DO
          END DO

!  r <- r - A * x

        ELSE
          DO i = 1, m
            ri = R( i )
            DO l = A_ptr( i ), A_ptr( i + 1 ) - 1
              ri = ri - A_val( l ) * X( A_col( l ) )
            END DO
            R( i ) = ri
          END DO
        END IF

      END IF
      RETURN

!  End of subroutine QPD_Ax

      END SUBROUTINE QPD_Ax

!  End of module GALAHAD_QPD_double

   END MODULE GALAHAD_QPD_double
