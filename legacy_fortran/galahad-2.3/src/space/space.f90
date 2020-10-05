! THIS VERSION: GALAHAD 2.2 - 22/04/2008 AT 14:30 GMT.

!-*-*-*-*-*-*-*-  G A L A H A D _ S P A C E   M O D U L E  *-*-*-*-*-*-*-*-*-*

!  This module contains simple routines for possibly changing the size
!  of allocatable or pointer arrays, and for deallocating them after use

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  January 27th 2005

   MODULE GALAHAD_SPACE_double

     USE GALAHAD_SYMBOLS

     IMPLICIT NONE     

     PRIVATE
     PUBLIC :: SPACE_resize_pointer, SPACE_resize_array,                       &
               SPACE_dealloc_pointer, SPACE_dealloc_array

!  Set precision

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  generic interfaces to cope with both pointer and allocatable 
!  real and integer arrays

     INTERFACE SPACE_resize_pointer
       MODULE PROCEDURE SPACE_resize_real_pointer,                             &
                        SPACE_resize_reallu_pointer,                           &
                        SPACE_resize_real2_pointer,                            &
                        SPACE_resize_integer_pointer,                          &
                        SPACE_resize_integerlu_pointer,                        &
                        SPACE_resize_integer2_pointer,                         &
                        SPACE_resize_logical_pointer,                          &
                        SPACE_resize_character_pointer
     END INTERFACE 

     INTERFACE SPACE_resize_array
       MODULE PROCEDURE SPACE_resize_real_array,                               &
                        SPACE_resize_reallu_array,                             &
                        SPACE_resize_real2_array,                              &
                        SPACE_resize_reallu2_array,                            &
                        SPACE_resize_integer_array,                            &
                        SPACE_resize_integerlu_array,                          &
                        SPACE_resize_integer2_array,                           &
                        SPACE_resize_logical_array,                            &
                        SPACE_resize_character_array
     END INTERFACE 

     INTERFACE SPACE_dealloc_pointer
       MODULE PROCEDURE SPACE_dealloc_real_pointer,                            &
                        SPACE_dealloc_real2_pointer,                           &
                        SPACE_dealloc_integer_pointer,                         &
                        SPACE_dealloc_integer2_pointer,                        &
                        SPACE_dealloc_logical_pointer,                         &
                        SPACE_dealloc_character_pointer
     END INTERFACE 

     INTERFACE SPACE_dealloc_array
       MODULE PROCEDURE SPACE_dealloc_real_array,                              &
                        SPACE_dealloc_real2_array,                             &
                        SPACE_dealloc_integer_array,                           &
                        SPACE_dealloc_integer2_array,                          &
                        SPACE_dealloc_logical_array,                           &
                        SPACE_dealloc_character_array
     END INTERFACE 

   CONTAINS

!  *-  S P A C E _ R E S I Z E _ R E A L _ P O I N T E R  S U B R O U T I N E   -*-

     SUBROUTINE SPACE_resize_real_pointer( len, point, status, alloc_status,   &
       deallocate_error_fatal, point_name, exact_size, bad_alloc, out )

!  Ensure that the real pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point ) /= len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point ) < len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point ) < len ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_real_pointer

     END SUBROUTINE SPACE_resize_real_pointer

! -  S P A C E _ R E S I Z E _ R E A L L U _ P O I N T E R  S U B R O U T I N E  

     SUBROUTINE SPACE_resize_reallu_pointer( l, u, point, status, alloc_status, &
       deallocate_error_fatal, point_name, exact_size, bad_alloc, out )

!  Ensure that the real pointer array "point" has bounds at least l and u

!  If exact_size is prsent and true, point is reallocated to have bounds l and u
!  Otherwise point is only reallocated if its bounds are insufficient to cover
!  l and u

!  Dummy arguments

     INTEGER, INTENT( IN ) :: l, u
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( LBOUND( point, 1 ) /= l .OR. UBOUND( point, 1 ) /= u ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( LBOUND( point, 1 ) > l .OR. UBOUND( point, 1 ) < u ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( LBOUND( point, 1 ) > l .OR. UBOUND( point, 1 ) < u ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( l : u ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_reallu_pointer

     END SUBROUTINE SPACE_resize_reallu_pointer

! -  S P A C E _ R E S I Z E _ R E A L 2 _ P O I N T E R  S U B R O U T I N E   -

     SUBROUTINE SPACE_resize_real2_pointer( len1, len2, point, status,         &
       alloc_status, deallocate_error_fatal, point_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the 2D real pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len1, len2
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), POINTER, DIMENSION( : , : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point, 1 ) /= len1 .OR.                                  &
                SIZE( point, 2 ) /= len2  ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point, 1 ) < len1 .OR.                                   &
                SIZE( point, 2 ) < len2  ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point, 1 ) < len1 .OR.                                     &
              SIZE( point, 2 ) < len2  ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len1, len2 ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_real2_pointer

     END SUBROUTINE SPACE_resize_real2_pointer

! - S P A C E _ R E S I Z E _ I N T E G E R _ P O I N T E R  S U B R O U T I N E 

     SUBROUTINE SPACE_resize_integer_pointer( len, point, status, alloc_status, &
       deallocate_error_fatal, point_name, exact_size, bad_alloc, out )
 
!  Ensure that the integer pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point ) /= len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point ) < len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point ) < len ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integer_pointer

     END SUBROUTINE SPACE_resize_integer_pointer

! - S P A C E _ R E S I Z E _ I N T E G E R L U _ P O I N T E R  SUBROUTINE -  

     SUBROUTINE SPACE_resize_integerlu_pointer( l, u, point, status,           &
      alloc_status, deallocate_error_fatal, point_name, exact_size,            &
      bad_alloc, out )

!  Ensure that the integer pointer array "point" has bounds at least l and u

!  If exact_size is prsent and true, point is reallocated to have bounds l and u
!  Otherwise point is only reallocated if its bounds are insufficient to cover
!  l and u

!  Dummy arguments

     INTEGER, INTENT( IN ) :: l, u
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( LBOUND( point, 1 ) /= l .OR. UBOUND( point, 1 ) /= u ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( LBOUND( point, 1 ) > l .OR. UBOUND( point, 1 ) < u ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( LBOUND( point, 1 ) > l .OR. UBOUND( point, 1 ) < u ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( l : u ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integerlu_pointer

     END SUBROUTINE SPACE_resize_integerlu_pointer

! -*-  S P A C E _ R E S I Z E _ I N T E G E R 2 _ P O I N T E R  SUBROUTINE -*-

     SUBROUTINE SPACE_resize_integer2_pointer( len1, len2, point, status,         &
       alloc_status, deallocate_error_fatal, point_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the 2D integer pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len1, len2
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, POINTER, DIMENSION( : , : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point, 1 ) /= len1 .OR.                                  &
                SIZE( point, 2 ) /= len2  ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point, 1 ) < len1 .OR.                                   &
                SIZE( point, 2 ) < len2  ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point, 1 ) < len1 .OR.                                     &
              SIZE( point, 2 ) < len2  ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len1, len2 ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integer2_pointer

     END SUBROUTINE SPACE_resize_integer2_pointer

! - S P A C E _ R E S I Z E _ L O G I C A L _ P O I N T E R  S U B R O U T I N E 

     SUBROUTINE SPACE_resize_logical_pointer( len, point, status, alloc_status, &
       deallocate_error_fatal, point_name, exact_size, bad_alloc, out )

!  Ensure that the logical pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     LOGICAL, POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point ) /= len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point ) < len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point ) < len ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_logical_pointer

     END SUBROUTINE SPACE_resize_logical_pointer

! S P A C E _ R E S I Z E _ C H A R A C T E R _ P O I N T E R S U B R O U T I N E

     SUBROUTINE SPACE_resize_character_pointer( len, point, status,            &
       alloc_status, deallocate_error_fatal, point_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the character pointer array "point" is of lenth at least len.

!  If exact_size is prsent and true, point is reallocated to be of size len. 
!  Otherwise point is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     CHARACTER( LEN = * ), POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( point ) /= len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( point ) < len ) THEN 
             CALL SPACE_dealloc_pointer( point, status, alloc_status,          &
                                         point_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( point ) < len ) THEN 
           CALL SPACE_dealloc_pointer( point, status, alloc_status,            &
                                       point_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate point to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( point( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )                 &
         bad_alloc = point_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( point_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_character_pointer

     END SUBROUTINE SPACE_resize_character_pointer

!  *-  S P A C E _ R E S I Z E _ R E A L _ A R R A Y  S U B R O U T I N E   -*-

     SUBROUTINE SPACE_resize_real_array( len, array, status, alloc_status,     &
       deallocate_error_fatal, array_name, exact_size, bad_alloc, out )

!  Ensure that the real allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array ) /= len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array ) < len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array ) < len ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,             &
                                     array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_real_array

     END SUBROUTINE SPACE_resize_real_array

! -  S P A C E _ R E S I Z E _ R E A L L U _ A R R A Y   S U B R O U T I N E  

     SUBROUTINE SPACE_resize_reallu_array( l, u, array, status, alloc_status,  &
       deallocate_error_fatal, array_name, exact_size, bad_alloc, out )

!  Ensure that the real allocatable array "array" has bounds at least l and u

!  If exact_size is prsent and true, array is reallocated to have bounds l and u
!  Otherwise array is only reallocated if its bounds are insufficient to cover
!  l and u

!  Dummy arguments

     INTEGER, INTENT( IN ) :: l, u
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( LBOUND( array, 1 ) /= l .OR. UBOUND( array, 1 ) /= u ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                         array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( LBOUND( array, 1 ) > l .OR. UBOUND( array, 1 ) < u ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                         array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( LBOUND( array, 1 ) > l .OR. UBOUND( array, 1 ) < u ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,              &
                                       array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( l : u ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_reallu_array

     END SUBROUTINE SPACE_resize_reallu_array

!  *-  S P A C E _ R E S I Z E _ R E A L 2 _ A R R A Y  S U B R O U T I N E   -*

     SUBROUTINE SPACE_resize_real2_array( len1, len2, array, status,           &
       alloc_status, deallocate_error_fatal, array_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the 2D real allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len1, len2
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array, 1 ) /= len1 .OR.                                  &
                SIZE( array, 2 ) /= len2  ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array, 1 ) < len1 .OR.                                   &
                SIZE( array, 2 ) < len2  ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array, 1 ) < len1 .OR.                                     &
              SIZE( array, 2 ) < len2  ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,              &
                                       array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len1, len2 ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_real2_array

     END SUBROUTINE SPACE_resize_real2_array

! -  S P A C E _ R E S I Z E _ R E A L L U 2 _ A R R A Y   S U B R O U T I N E  

     SUBROUTINE SPACE_resize_reallu2_array( l1l, l1u, l2, array,                &
       status, alloc_status, deallocate_error_fatal,                            &
       array_name, exact_size, bad_alloc, out )

!  Ensure that the real allocatable array "array" has bounds at least l1l 
!  and l1u for its first argument and l2 for its second

!  If exact_size is prsent and true, array is reallocated to have bounds l and u
!  Otherwise array is only reallocated if its bounds are insufficient to cover
!  l1l and l1u for ithe first argument and l2 for the second

!  Dummy arguments

     INTEGER, INTENT( IN ) :: l1l, l1u, l2
     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( LBOUND( array, 1 ) /= l1l .OR. UBOUND( array, 1 ) /= l1u .OR.  &
                SIZE( array, 2 ) /= l2 ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( LBOUND( array, 1 ) > l1l .OR. UBOUND( array, 1 ) < l1u .OR.    &
                SIZE( array, 2 ) < l2 ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( LBOUND( array, 1 ) > l1l .OR. UBOUND( array, 1 ) < l1u .OR.      &
              SIZE( array, 2 ) < l2 ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,              &
                                       array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( l1l : l1u , l2 ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_reallu2_array

     END SUBROUTINE SPACE_resize_reallu2_array

!  - S P A C E _ R E S I Z E _ I N T E G E R _ A R R A Y  S U B R O U T I N E -

     SUBROUTINE SPACE_resize_integer_array( len, array, status, alloc_status,  &
       deallocate_error_fatal, array_name, exact_size, bad_alloc, out )

!  Ensure that the integer allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array ) /= len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array ) < len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array ) < len ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,             &
                                     array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integer_array

     END SUBROUTINE SPACE_resize_integer_array

! -*-  S P A C E _ R E S I Z E _ I N T E G E R L U _ A R R A Y   SUBROUTINE  -*-

     SUBROUTINE SPACE_resize_integerlu_array( l, u, array, status,             &
       alloc_status, deallocate_error_fatal, array_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the integer allocatable array "array" has bounds at least l and u

!  If exact_size is prsent and true, array is reallocated to have bounds l and u
!  Otherwise array is only reallocated if its bounds are insufficient to cover
!  l and u

!  Dummy arguments

     INTEGER, INTENT( IN ) :: l, u
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( LBOUND( array, 1 ) /= l .OR. UBOUND( array, 1 ) /= u ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( LBOUND( array, 1 ) > l .OR. UBOUND( array, 1 ) < u ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( LBOUND( array, 1 ) > l .OR. UBOUND( array, 1 ) < u ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,              &
                                     array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( l : u ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integerlu_array

     END SUBROUTINE SPACE_resize_integerlu_array

!  -*-  S P A C E _ R E S I Z E _ I N T E G E R 2 _ A R R A Y  SUBROUTINE   -*-

     SUBROUTINE SPACE_resize_integer2_array( len1, len2, array, status,        &
       alloc_status, deallocate_error_fatal, array_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the 2D integer allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len1, len2
     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, ALLOCATABLE, DIMENSION( : , : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array, 1 ) /= len1 .OR.                                  &
                SIZE( array, 2 ) /= len2  ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array, 1 ) < len1 .OR.                                   &
                SIZE( array, 2 ) < len2  ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,            &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array, 1 ) < len1 .OR.                                     &
              SIZE( array, 2 ) < len2  ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,              &
                                       array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len1, len2 ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_integer2_array

     END SUBROUTINE SPACE_resize_integer2_array

!  -  S P A C E _ R E S I Z E _ L O G I C A L _ A R R A Y  S U B R O U T I N E  -

     SUBROUTINE SPACE_resize_logical_array( len, array, status, alloc_status,  &
       deallocate_error_fatal, array_name, exact_size, bad_alloc, out )

!  Ensure that the logical allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     LOGICAL, ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array ) /= len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array ) < len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array ) < len ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,             &
                                     array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_logical_array

     END SUBROUTINE SPACE_resize_logical_array

! - S P A C E _ R E S I Z E _ C H A R A C T E R _ A R R A Y  S U B R O U T I N E 

     SUBROUTINE SPACE_resize_character_array( len, array, status,              &
       alloc_status, deallocate_error_fatal, array_name, exact_size,           &
       bad_alloc, out )

!  Ensure that the character allocatable array "array" is of lenth at least len.

!  If exact_size is prsent and true, array is reallocated to be of size len. 
!  Otherwise array is only reallocated if its length is currently smaller 
!  than len

!  Dummy arguments

     INTEGER, INTENT( IN ) :: len
     INTEGER, INTENT( OUT ) :: status, alloc_status
     CHARACTER ( LEN = * ), ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     LOGICAL, OPTIONAL :: deallocate_error_fatal, exact_size
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

!  Local variable

     LOGICAL :: reallocate

!  Check to see if a reallocation (or initial allocation) is needed

     status = GALAHAD_ok ; alloc_status = 0 ; reallocate = .TRUE.
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array ) ) THEN
       IF ( PRESENT( exact_size ) ) THEN
         IF ( exact_size ) THEN
           IF ( SIZE( array ) /= len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         ELSE
           IF ( SIZE( array ) < len ) THEN 
             CALL SPACE_dealloc_array( array, status, alloc_status,           &
                                       array_name, bad_alloc, out )
           ELSE ; reallocate = .FALSE.
           END IF
         END IF
       ELSE
         IF ( SIZE( array ) < len ) THEN 
           CALL SPACE_dealloc_array( array, status, alloc_status,             &
                                     array_name, bad_alloc, out )
          ELSE ; reallocate = .FALSE.
          END IF
       END IF
     END IF

!  If a deallocation error occured, return if desired

     IF ( PRESENT( deallocate_error_fatal ) ) THEN
       IF ( deallocate_error_fatal .AND. alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     ELSE
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate ; RETURN
       END IF
     END IF

!  Reallocate array to be of length len, checking for error returns

     IF ( reallocate ) ALLOCATE( array( len ), STAT = alloc_status )
     IF ( alloc_status /= 0 ) THEN
       status = GALAHAD_error_allocate
       IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )                 &
         bad_alloc = array_name
       IF ( PRESENT( out ) ) THEN
         IF ( PRESENT( array_name ) ) THEN
           IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
         ELSE
           IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Allocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Allocation error status = ', I6 ) 

!  End of SPACE_resize_character_array

     END SUBROUTINE SPACE_resize_character_array

!-  S P A C E _ D E A L L O C _ R E A L _ P O I N T E R   S U B R O U T I N E  -*

     SUBROUTINE SPACE_dealloc_real_pointer( point, status, alloc_status,       &
                                            point_name, bad_alloc, out )

!  Deallocate the real pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )               &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_real_pointer

     END SUBROUTINE SPACE_dealloc_real_pointer

!-  S P A C E _ D E A L L O C _ R E A L 2 _ P O I N T E R  S U B R O U T I N E -*

     SUBROUTINE SPACE_dealloc_real2_pointer( point, status, alloc_status,      &
                                             point_name, bad_alloc, out )

!  Deallocate the rank-2 real pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), POINTER, DIMENSION( : , : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )               &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_real2_pointer

     END SUBROUTINE SPACE_dealloc_real2_pointer

!-  S P A C E _ D E A L L O C _ I N T E G E R _ P O I N T E R   SUBROUTINE  -*-

     SUBROUTINE SPACE_dealloc_integer_pointer( point, status, alloc_status,    &
                                               point_name, bad_alloc, out )

!  Deallocate the integer pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )              &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_integer_pointer

     END SUBROUTINE SPACE_dealloc_integer_pointer

! -*- S P A C E _ D E A L L O C _ I N T E G E R 2 _ P O I N T E R  SUBROUTINE -*-

     SUBROUTINE SPACE_dealloc_integer2_pointer( point, status, alloc_status,    &
                                                point_name, bad_alloc, out )

!  Deallocate the rank-2 integer pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, POINTER, DIMENSION( : , : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )               &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_integer2_pointer

     END SUBROUTINE SPACE_dealloc_integer2_pointer

!-*-  S P A C E _ D E A L L O C _ L O G I C A L _ P O I N T E R   SUBROUTINE -*-

     SUBROUTINE SPACE_dealloc_logical_pointer( point, status, alloc_status,    &
                                               point_name, bad_alloc, out )

!  Deallocate the logical pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     LOGICAL, POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )               &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_logical_pointer

     END SUBROUTINE SPACE_dealloc_logical_pointer

!-  S P A C E _ D E A L L O C _ C H A R A C T E R  _ P O I N T E R  SUBROUTINE -*

     SUBROUTINE SPACE_dealloc_character_pointer( point, status, alloc_status,  &
                                                 point_name, bad_alloc, out )

!  Deallocate the character pointer array "point"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     CHARACTER( LEN = * ), POINTER, DIMENSION( : ) :: point
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: point_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ASSOCIATED( point) ) THEN
       DEALLOCATE( point, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( point_name ) )               &
           bad_alloc = point_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( point_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( point_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_character_pointer

     END SUBROUTINE SPACE_dealloc_character_pointer

!-*-  S P A C E _ D E A L L O C _ R E A L _ A R R A Y   S U B R O U T I N E  -*-

     SUBROUTINE SPACE_dealloc_real_array( array, status, alloc_status,         &
                                          array_name, bad_alloc, out )

!  Deallocate the real allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )              &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_real_array

     END SUBROUTINE SPACE_dealloc_real_array

!-*-  S P A C E _ D E A L L O C _ R E A L 2 _ A R R A Y  S U B R O U T I N E -*-

     SUBROUTINE SPACE_dealloc_real2_array( array, status, alloc_status,        &
                                           array_name, bad_alloc, out )

!  Deallocate the rank-2 real allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : , : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )               &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_real2_array

     END SUBROUTINE SPACE_dealloc_real2_array

!- S P A C E _ D E A L L O C _ I N T E G E R _ A R R A Y   S U B R O U T I N E - 

     SUBROUTINE SPACE_dealloc_integer_array( array, status, alloc_status,      &
                                             array_name, bad_alloc, out )

!  Deallocate the integer allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )              &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_integer_array

     END SUBROUTINE SPACE_dealloc_integer_array

!- S P A C E _ D E A L L O C _ I N T E G E R 2 _ A R R A Y  S U B R O U T I N E -

     SUBROUTINE SPACE_dealloc_integer2_array( array, status, alloc_status,     &
                                              array_name, bad_alloc, out )

!  Deallocate the rank-2 integer allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     INTEGER, ALLOCATABLE, DIMENSION( : , : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )               &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_integer2_array

     END SUBROUTINE SPACE_dealloc_integer2_array

!-*-*-  S P A C E _ D E A L L O C _ L O G I C A L _ A R R A Y   SUBROUTINE -*-*-

     SUBROUTINE SPACE_dealloc_logical_array( array, status, alloc_status,      &
                                             array_name, bad_alloc, out )

!  Deallocate the logical allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     LOGICAL, ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )               &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_logical_array

     END SUBROUTINE SPACE_dealloc_logical_array

!-  S P A C E _ D E A L L O C _ C H A R A C T E R  _ A R R A Y  SUBROUTINE  -*

     SUBROUTINE SPACE_dealloc_character_array( array, status, alloc_status,    &
                                               array_name, bad_alloc, out )

!  Deallocate the character allocatable array "array"

!  Dummy arguments

     INTEGER, INTENT( OUT ) :: status, alloc_status
     CHARACTER( LEN = * ), ALLOCATABLE, DIMENSION( : ) :: array
     INTEGER, OPTIONAL :: out
     CHARACTER ( LEN = 80 ), OPTIONAL :: array_name
     CHARACTER ( LEN = 80 ), OPTIONAL :: bad_alloc

     status = GALAHAD_ok ; alloc_status = 0
     IF ( PRESENT( bad_alloc ) ) bad_alloc = ''
     IF ( ALLOCATED( array) ) THEN
       DEALLOCATE( array, STAT = alloc_status )
       IF ( alloc_status /= 0 ) THEN
         status = GALAHAD_error_deallocate
         IF ( PRESENT( bad_alloc ) .AND. PRESENT( array_name ) )              &
           bad_alloc = array_name
         IF ( PRESENT( out ) ) THEN
           IF ( PRESENT( array_name ) ) THEN
             IF ( out > 0 ) WRITE( out, 2900 ) TRIM( array_name ), alloc_status
           ELSE
             IF ( out > 0 ) WRITE( out, 2910 ) alloc_status
           END IF
         END IF
       END IF
     END IF
     RETURN

!  Non-executable statements

2900 FORMAT( ' ** Deallocation error for ', A, /, '     status = ', I6 ) 
2910 FORMAT( ' ** Deallocation error status = ', I6 ) 

!  End of subroutine SPACE_dealloc_character_array

     END SUBROUTINE SPACE_dealloc_character_array

!  End of module GALAHAD_SPACE_double

   END MODULE GALAHAD_SPACE_double
