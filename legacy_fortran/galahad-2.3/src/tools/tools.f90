! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
! 30/06/2003: procedures _is_are and _s moved from QPA.
! 23/07/2003: invalid "END INTERFACE" arguments removed

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*                                   *-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*        TOOLS   M O D U L E        *-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*                                   *-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

!  Copyright reserved: Nick Gould and Philippe Toint, for GALAHAD productions
!  July 2001
!
!               +---------------------------------------------+
!               |                                             |
!               |         Provides various simple tools       |
!               |                                             |
!               +---------------------------------------------+
!

   MODULE GALAHAD_TOOLS_double

!-------------------------------------------------------------------------------
!   A c c e s s 
!-------------------------------------------------------------------------------

      IMPLICIT NONE

!     Make everything private by default

      PRIVATE

!     Make the tools public

      PUBLIC :: TOOLS_output_vector

!            Outputs a vector with a reasonable layout.

      INTERFACE TOOLS_output_vector
         MODULE PROCEDURE TOOLS_output_vector_integer, TOOLS_output_vector_real
      END INTERFACE

      PUBLIC :: TOOLS_output_matrix_C, TOOLS_output_matrix_S,                  &
                TOOLS_output_matrix_D, TOOLS_is_are, TOOLS_s

!            Outputs a matrix with a reasonable layout.

      INTERFACE TOOLS_output_matrix_C
         MODULE PROCEDURE TOOLS_output_matrix_real_C
      END INTERFACE

      INTERFACE TOOLS_output_matrix_S
         MODULE PROCEDURE TOOLS_output_matrix_real_S
      END INTERFACE

      INTERFACE TOOLS_output_matrix_D
         MODULE PROCEDURE TOOLS_output_matrix_real_D
      END INTERFACE


!-------------------------------------------------------------------------------
!   P r e c i s i o n
!-------------------------------------------------------------------------------

      INTEGER, PRIVATE, PARAMETER :: sp = KIND( 1.0E+0 )
      INTEGER, PRIVATE, PARAMETER :: dp = KIND( 1.0D+0 )
      INTEGER, PRIVATE, PARAMETER :: wp = dp

!-------------------------------------------------------------------------------
!   O t h e r s
!-------------------------------------------------------------------------------

      INTEGER, PRIVATE, PARAMETER :: OK            =   0

   CONTAINS

!==============================================================================
!==============================================================================


      SUBROUTINE TOOLS_output_vector_real( n, x, out )

!     Print x

!     Arguments:

      INTEGER, INTENT( IN ) :: n

!             The dimension of x.

      REAL ( KIND = wp ), DIMENSION( n ), INTENT( IN ) :: x

!             The vector to print.

      INTEGER, INTENT( IN ) :: out

!             The output device number

!     Programming: Ph. L. Toint, November 2002.

!==============================================================================

!     Local variable

      INTEGER :: j, i

      WRITE( out, 101 ) 
      j = 1
      DO i = 1, n / 5
         WRITE( out, 100 ) j, x( j:j+4 )
         j = j + 5
      END DO
      IF ( j <= n ) WRITE( out, 100 ) j, x( j:n )
      WRITE( out, 101 ) 

      RETURN

100   FORMAT( 1x, i4, 5( 1x, 1pE14.6 ) )
101   FORMAT( / )

      END SUBROUTINE TOOLS_output_vector_real

!==============================================================================
!==============================================================================

      SUBROUTINE TOOLS_output_vector_integer( n, ix, out )

!     Print ix

!     Arguments:

      INTEGER, INTENT( IN ) :: n

!             The dimension of ix.

      INTEGER, DIMENSION( n ), INTENT( IN ) :: ix

!             The vector to print.

      INTEGER, INTENT( IN ) :: out

!             The output device number

!     Programming: Ph. L. Toint, November 2002.

!==============================================================================

!     Local variables

      INTEGER :: j, i

      WRITE( out, 101 )
      j = 1
      DO i = 1, n / 10
         WRITE( out, 100 ) j, ix( j:j+9 )
         j = j + 10
      END DO
      IF ( j <= n ) WRITE( out, 100 ) j, ix( j:n )
      WRITE( out, 101 )

      RETURN

100   FORMAT( 1x, i4, 2x, 10( 1x, i5 ) )
101   FORMAT( / )

      END SUBROUTINE TOOLS_output_vector_integer

!==============================================================================
!==============================================================================

      SUBROUTINE TOOLS_output_matrix_real_C( nnz, A_val, A_row, A_col, out )

      INTEGER, INTENT( IN ) :: nnz

      REAL( KIND = wp ), DIMENSION( nnz ), INTENT( IN ) :: A_val

      INTEGER, DIMENSION( nnz ), INTENT( IN ) :: A_row

      INTEGER, DIMENSION( nnz ), INTENT( IN ) :: A_col

      INTEGER, INTENT( IN ) :: out

!     Programming: Ph. L. Toint, November 2002.

!==============================================================================

!     Local variables

      INTEGER :: k, kk

      WRITE( out, 102 )
      k = 0
      DO kk = 1, nnz / 3
         WRITE( out, 100 ) A_row( k + 1 ), A_col( k + 1 ), A_val( k + 1 ),     &
                           A_row( k + 2 ), A_col( k + 2 ), A_val( k + 2 ),     &
                           A_row( k + 3 ), A_col( k + 3 ), A_val( k + 3 )
         k = k + 3
      END DO
      IF ( k < nnz ) THEN
         SELECT CASE ( nnz - k )
         CASE ( 1 )
            WRITE ( out, 100 ) A_row( nnz ), A_col( nnz ), A_val( nnz )
         CASE ( 2 )
            WRITE ( out, 100 ) A_row( k + 1 ), A_col( k + 1 ), A_val( k + 1  ),&
                               A_row( nnz ), A_col( nnz ), A_val( nnz )
         END SELECT            
      END IF
      WRITE( out, 101 )

      RETURN

100   FORMAT( 2( 1x, i4), 2x, 1pE12.4, 2( 4x, 2( 1x, i4), 2x, 1pE12.4 ) )
101   FORMAT( / )
102   FORMAT( /,1x,'   i    j       value  ',2(5x,'   i    j       value  '),/)

      END SUBROUTINE TOOLS_output_matrix_real_C

!==============================================================================
!==============================================================================

      SUBROUTINE TOOLS_output_matrix_real_S( nnz, A_val, A_ptr, A_col, out )

      INTEGER, INTENT( IN ) :: nnz

      REAL( KIND = wp ), DIMENSION( nnz ), INTENT( IN ) :: A_val

      INTEGER, DIMENSION( nnz ), INTENT( IN ) :: A_ptr

      INTEGER, DIMENSION( nnz ), INTENT( IN ) :: A_col

      INTEGER, INTENT( IN ) :: out

!     Programming: Ph. L. Toint, November 2002.

!==============================================================================

!     Local variables

      INTEGER :: k, kk, i, i1, i2, i3

      WRITE( out, 102 )
      k = 0
      i = 1
      DO kk = 1, nnz / 3
         DO
            IF ( k + 1  == A_ptr( i + 1 ) ) EXIT
            i = i + 1
         END DO
         i1 = i
         DO
            IF ( k + 2  == A_ptr( i + 1 ) ) EXIT
            i = i + 1
         END DO
         i2 = i
         DO
            IF ( k + 3  == A_ptr( i + 1 ) ) EXIT
            i = i + 1
         END DO
         i3 = i
         WRITE( out, 100 ) i1, A_col( k + 1 ), A_val( k + 1 ),                 &
                           i2, A_col( k + 2 ), A_val( k + 2 ),                 &
                           i3, A_col( k + 3 ), A_val( k + 3 )
         k = k + 3
      END DO
      IF ( k < nnz ) THEN
         SELECT CASE ( nnz - k )
         CASE ( 1 )
            DO
               IF ( nnz  == A_ptr( i + 1 ) ) EXIT
               i = i + 1
            END DO
            WRITE ( out, 100 ) i, A_col( nnz ), A_val( nnz )
         CASE ( 2 )
            DO
               IF ( k + 1  == A_ptr( i + 1 ) ) EXIT
               i = i + 1
            END DO
            i1 = i
            DO
               IF ( nnz == A_ptr( i + 1 ) ) EXIT
               i = i + 1
            END DO
            i2 = i
            WRITE ( out, 100 ) i1, A_col( k + 1 ), A_val( k + 1  ),            &
                               i2, A_col( nnz ), A_val( nnz )
         END SELECT            
      END IF
      WRITE( out, 101 )

      RETURN

100   FORMAT( 2( 1x, i4), 2x, 1pE12.4, 2( 4x, 2( 1x, i4), 2x, 1pE12.4 ) )
101   FORMAT( / )
102   FORMAT(/,1x,'   i    j       value  ',2(5x,'   i    j       value  '),/ )

      END SUBROUTINE TOOLS_output_matrix_real_S

!==============================================================================
!==============================================================================

      SUBROUTINE TOOLS_output_matrix_real_D( nrow, ncol, A_val, sym, out )

      INTEGER, INTENT( IN ) :: nrow, ncol

      REAL( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: A_val

      LOGICAL, INTENT( IN ) :: sym

      INTEGER, INTENT( IN ) :: out

!     Programming: Ph. L. Toint, November 2002.

!==============================================================================

!     Local variables

      INTEGER :: k, kk, i, i1, i2, j, j1, j2, nval, lrow

      WRITE( out, 102 )
      nval = nrow * ncol
      k = 0
      i = 1
      j = 0
      IF ( sym ) THEN
         lrow = 1
      ELSE
         lrow = ncol
      END IF
      DO kk = 1, nval / 3
         j = j + 1
         IF ( j > lrow ) THEN
            i = i + 1
            IF ( sym ) lrow = i
            j = 1
         END IF
         i1 = i
         j1 = j
         j  = j + 1
         IF ( j > lrow ) THEN
            i = i + 1
            IF ( sym ) lrow = i
            j = 1
         END IF
         i2 = i
         j2 = j
         j  = j + 1
         IF ( j > lrow ) THEN
            i = i + 1
            IF ( sym ) lrow = i
            j = 1
         END IF
         WRITE( out, 100 ) i1, j1, A_val( k + 1 ),                             &
                           i2, j2, A_val( k + 2 ),                             &
                           i , j , A_val( k + 3 )
         k = k + 3
      END DO
      IF ( k < nval ) THEN
         SELECT CASE ( nval - k )
         CASE ( 1 )
            j = j + 1
            IF ( j > lrow ) THEN
               i = i + 1
               IF ( sym ) lrow = i
               j = 1
            END IF
            WRITE ( out, 100 ) i, j, A_val( nval )
         CASE ( 2 )
            j = j + 1
            IF ( j > lrow ) THEN
               i = i + 1
               IF ( sym ) lrow = i
               j = 1
            END IF
            i1 = i
            j1 = j
            j  = j + 1
            IF ( j > lrow ) THEN
               i = i + 1
               IF ( sym ) lrow = i
               j = 1
            END IF
            WRITE ( out, 100 ) i1, j1, A_val( k + 1  ),            &
                               i , j , A_val( nval )
         END SELECT            
      END IF
      WRITE( out, 101 )

      RETURN

100   FORMAT( 2( 1x, i4), 2x, 1pE12.4, 2( 4x, 2( 1x, i4), 2x, 1pE12.4 ) )
101   FORMAT( / )
102   FORMAT(/,1x,'   i    j       value  ',2(5x,'   i    j       value  '),/)

      END SUBROUTINE TOOLS_output_matrix_real_D

!==============================================================================
!==============================================================================

!  Function that returns "is" for 1 item, "are" for any other number of items.
!  This is used for i/o purposes
!  Nick Gould, 1999

      FUNCTION TOOLS_is_are( num )
      CHARACTER ( len = 3 ) :: TOOLS_is_are
      INTEGER :: num
      IF ( num /= 1 ) THEN
        TOOLS_is_are = 'are'
      ELSE
        TOOLS_is_are = 'is '
      END IF
      RETURN
      END FUNCTION TOOLS_is_are

!==============================================================================
!==============================================================================

!  Function that returns " " for 1 item, "s" for any other number of items
!  This is used for i/o purposes
!  Nick Gould, 1999

      FUNCTION TOOLS_s( num )
      CHARACTER ( len = 1 ) :: TOOLS_s
      INTEGER :: num
      IF ( num /= 1 ) THEN
        TOOLS_s = 's'
      ELSE
        TOOLS_s = ' '
      END IF
      RETURN
      END FUNCTION TOOLS_s

!==============================================================================
!==============================================================================

   END MODULE GALAHAD_TOOLS_double

!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*    END TOOLS  M O D U L E   *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*                             *-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*







