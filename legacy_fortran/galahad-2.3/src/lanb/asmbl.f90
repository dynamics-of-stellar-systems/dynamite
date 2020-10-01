! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-  L A N C E L O T  -B-  ASMBL  M O D U L E  *-*-*-*-*-*-*-*

!  Nick Gould, for GALAHAD productions
!  Copyright reserved
!  January 25th 1995

   MODULE LANCELOT_ASMBL_double

     USE LANCELOT_EXTEND_double, ONLY: EXTEND_arrays
     IMPLICIT NONE

     PRIVATE
     PUBLIC :: ASMBL_save_type, ASMBL_assemble_hessian

!  Set precision

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  Set other parameters

     REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
     REAL ( KIND = wp ), PARAMETER :: one = 1.0_wp

!  ================================
!  The ASMBL_save_type derived type
!  ================================

     TYPE :: ASMBL_save_type
       LOGICAL :: ptr_status
       INTEGER, DIMENSION( 30 ) :: ICNTL
       INTEGER, DIMENSION( 20 ) :: INFO
       REAL ( KIND = wp ), DIMENSION( 5 ) :: CNTL
     END TYPE ASMBL_save_type

   CONTAINS

!-*-*-  L A N C E L O T  -B-  ASMBL_assemble_hessian  S U B R O U T I N E -*-

     SUBROUTINE ASMBL_assemble_hessian(                                        &
                      n , ng, nel   , ntotel, nvrels, nnza  , maxsel,          &
                      nvargp, nfree , IFREE , ISTADH, ICNA  ,                  &
                      ISTADA, INTVAR, IELVAR, IELING, ISTADG, ISTAEV,          &
                      ISTAGV, ISVGRP, A     , GUVALS, lnguvl, HUVALS,          &
                      lnhuvl, GVALS2, GVALS3, GSCALE, ESCALE, GXEQX ,          &
                      ITYPEE, INTREP, RANGE , iprint, error , out   ,          &
                      buffer, use_band, no_zeros, fixed_structure, nsemib,     &
                      status, alloc_status, bad_alloc,                         &
                      S, lirnh, ljcnh, lh, IRNH, JCNH, H,                      &
                      LINK_col, POS_in_H, llink, lpos,                         &
                      IW_asmbl, GRAD_el, W_el, W_in, H_el, H_in, skipg,        &
                      nnzh, maxsbw, DIAG, OFFDIA, KNDOFG )

!  Assemble the second derivative matrix of a groups partially separable
!  function in either co-ordinate or band format

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     INTEGER, INTENT( IN  ) :: n , ng, maxsel, nsemib, nvargp, nnza
     INTEGER, INTENT( IN  ) :: nvrels, ntotel, nfree , nel   , buffer
     INTEGER, INTENT( IN  ) :: lnguvl, lnhuvl, iprint, error , out
     INTEGER, INTENT( OUT ) :: status, alloc_status
     LOGICAL, INTENT( IN  ) :: use_band, no_zeros, fixed_structure, skipg
     CHARACTER ( LEN = 24 ) :: bad_alloc
     INTEGER, INTENT( IN  ), DIMENSION( n       ) :: IFREE
     INTEGER, INTENT( IN  ), DIMENSION( nnza    ) :: ICNA
     INTEGER, INTENT( IN  ), DIMENSION( ng  + 1 ) :: ISTADA, ISTADG, ISTAGV
     INTEGER, INTENT( IN  ), DIMENSION( nel + 1 ) :: INTVAR, ISTAEV, ISTADH
     INTEGER, INTENT( IN  ), DIMENSION( nvrels  ) :: IELVAR
     INTEGER, INTENT( IN  ), DIMENSION( ntotel  ) :: IELING
     INTEGER, INTENT( IN  ), DIMENSION( nvargp  ) :: ISVGRP
     INTEGER, INTENT( IN  ), DIMENSION( nel ) :: ITYPEE
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( nnza ) :: A
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lnguvl ) :: GUVALS
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( lnhuvl ) :: HUVALS
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS2
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GVALS3
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ng ) :: GSCALE
     REAL ( KIND = wp ), INTENT( IN  ), DIMENSION( ntotel ) :: ESCALE
     LOGICAL, INTENT( IN ), DIMENSION( ng  ) :: GXEQX
     LOGICAL, INTENT( IN ), DIMENSION( nel ) :: INTREP
     TYPE ( ASMBL_save_type ), INTENT( INOUT ) :: S

!---------------------------------------------------------------
!   D u m m y   A r g u m e n t s   f o r   W o r k s p a c e 
!--------------------------------------------------------------

     INTEGER, INTENT( INOUT ) :: lirnh, ljcnh, lh, llink, lpos
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: LINK_col
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: POS_in_H
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: IRNH
     INTEGER, ALLOCATABLE, DIMENSION( : ) :: JCNH
     REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: H

     INTEGER, INTENT( OUT ), DIMENSION( : ) :: IW_asmbl
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: GRAD_el
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_el
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: W_in
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: H_el
     REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( : ) :: H_in

!--------------------------------------------------
!   O p t i o n a l   D u m m y   A r g u m e n t s
!--------------------------------------------------

     INTEGER, INTENT( OUT ), OPTIONAL :: maxsbw, nnzh
     REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL,                              &
                                        DIMENSION( nfree ) :: DIAG
     REAL ( KIND = wp ), INTENT( OUT ), OPTIONAL,                              &
                                        DIMENSION( nsemib, nfree ) :: OFFDIA
     INTEGER, INTENT( IN ), OPTIONAL, DIMENSION( ng ) :: KNDOFG

!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------

     INTERFACE
       SUBROUTINE RANGE( ielemn, transp, W1, W2, nelvar, ninvar, ieltyp,       &
                         lw1, lw2 )
       INTEGER, INTENT( IN ) :: ielemn, nelvar, ninvar, ieltyp, lw1, lw2
       LOGICAL, INTENT( IN ) :: transp
       REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( IN  ), DIMENSION ( lw1 ) :: W1
       REAL ( KIND = KIND( 1.0D+0 ) ), INTENT( OUT ), DIMENSION ( lw2 ) :: W2
       END SUBROUTINE RANGE
     END INTERFACE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

     INTEGER :: i , ii, j , jj, k , ig, ip, l , newpt , nin
     INTEGER :: iell  , iel   , ihnext, nvarel, ig1   , listvs, listve
     INTEGER :: ielh  , inext , ijhess, irow  , jcol  , jcolst, istart
     INTEGER :: nlh   , ulh   , mlh
     REAL ( KIND = wp ) :: wki   , hesnew, gdash , g2dash, scalee
     CHARACTER ( LEN = 2 ), DIMENSION( 36, 36 ) :: MATRIX
!    CHARACTER ( LEN = 80 ) :: array

!  If a band storage scheme is to be used, initialize the entries within the
!  band as zero

     IF ( use_band ) THEN
       maxsbw = 0
       DIAG = zero ; OFFDIA = zero

!  If a co-ordinate scheme is to be used, allocate arrays to hold the link
!  list which points to the row numbers  which are used in the columns of
!  the assembled Hessian

!  LINK_col( . ) gives the link list. The list for column J starts
!                 in LINK_col( J ) and ends when LINK_col( K ) = - 1
!  POS_in_H( . ) gives the position in H of the current link

     ELSE
       IF ( .NOT. S%ptr_status ) THEN
         S%ptr_status = .TRUE.
         llink = MIN( lirnh, ljcnh, lh ) + nfree
         lpos = MIN( lirnh, ljcnh, lh ) + nfree
       
         ALLOCATE( LINK_col( llink ), STAT = alloc_status )
         IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'LINK_c' ; GO TO 980
         END IF
       
         ALLOCATE( POS_in_H( lpos ), STAT = alloc_status )
         IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'POS_in' ; GO TO 980
         END IF
       
         ALLOCATE( IRNH( lirnh ), STAT = alloc_status )
         IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'IRNH' ; GO TO 980
         END IF
       
         ALLOCATE( JCNH( ljcnh ), STAT = alloc_status )
         IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'JCNH' ; GO TO 980
         END IF
       
         ALLOCATE( H( lh ), STAT = alloc_status )
         IF ( alloc_status /= 0 ) THEN ; bad_alloc = 'H' ; GO TO 980
         END IF
       
       ELSE
       
         k = MIN( lirnh, ljcnh, lh ) + nfree
         IF ( llink < k ) THEN
           DEALLOCATE( LINK_col ) ; llink = k
           ALLOCATE( LINK_col( llink ), STAT = alloc_status )
           IF ( alloc_status /= 0 ) THEN
             bad_alloc = 'LINK_c'; GO TO 980 ; END IF
         END IF
       
         k = MIN( lirnh, ljcnh, lh ) + nfree
         IF ( lpos < k ) THEN 
           DEALLOCATE( POS_in_H ) ; lpos = k
           ALLOCATE( POS_in_H( lpos ), STAT = alloc_status )
           IF ( alloc_status /= 0 ) THEN 
             bad_alloc = 'POS_in'; GO TO 980 ; END IF
         END IF
       
       END IF
       LINK_col( : nfree ) = - 1 ; POS_in_H( : nfree ) = - 1
       newpt = nfree

!  Make an initial allocation of the space required to hold the Hessian

       nnzh = 0
     END IF

!  Renumber the free variables so that they are variables 1 to NFREE

     IW_asmbl( : n ) = 0
     DO i = 1, nfree
       IW_asmbl( IFREE( i ) ) = i
     END DO
     IF ( iprint >= 10 ) WRITE( out, 2060 ) nfree, ( IFREE( i ), i = 1, nfree )

!  ------------------------------------------------------
!  Form the rank-one second order term for the I-th group
!  ------------------------------------------------------

     DO ig = 1, ng
       IF ( skipg ) THEN
         IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
       IF ( GXEQX( ig ) ) CYCLE
       IF ( .NOT. fixed_structure .AND. GSCALE( ig ) == zero ) CYCLE
       IF ( iprint >= 100 ) WRITE( out, 2070 ) ig
       g2dash = GSCALE( ig ) * GVALS3( ig )
       IF ( iprint >= 100 ) WRITE( 6, * ) ' GVALS3( ig ) ', GVALS3( ig )
       IF ( no_zeros .AND. g2dash == zero ) CYCLE
       ig1 = ig + 1
       listvs = ISTAGV( ig )
       listve = ISTAGV( ig1 ) - 1

!  Form the gradient of the IG-th group

       GRAD_el( ISVGRP( listvs : listve ) ) = zero

!  Consider any nonlinear elements for the group

       DO iell = ISTADG( ig ), ISTADG( ig1 ) - 1
         iel = IELING( iell )
         k = INTVAR( iel )
         l = ISTAEV( iel )
         nvarel = ISTAEV( iel + 1 ) - l
         scalee = ESCALE( iell )
         IF ( INTREP( iel ) ) THEN

!  The IEL-th element has an internal representation

           nin = INTVAR( iel + 1 ) - k
           CALL RANGE ( iel, .TRUE., GUVALS( k : k + nin - 1 ),                &
                        H_el, nvarel, nin, ITYPEE( iel ), nin, nvarel )
           DO i = 1, nvarel
             j = IELVAR( l )
             GRAD_el( j ) = GRAD_el( j ) + scalee * H_el( i )
             l = l + 1
           END DO
         ELSE

!  The IEL-th element has no internal representation

           DO i = 1, nvarel
             j = IELVAR( l )
             GRAD_el( j ) = GRAD_el( j ) + scalee * GUVALS( k )
             k = k + 1
             l = l + 1
           END DO
         END IF
       END DO

!  Include the contribution from the linear element

       DO k = ISTADA( ig ), ISTADA( ig1 ) - 1
         j = ICNA( k )
         GRAD_el( j ) = GRAD_el( j ) + A( k )
       END DO

!  The gradient is complete. Form the J-TH column of the rank-one matrix

       DO l = listvs, listve
         jj = ISVGRP( l )
         j = IW_asmbl( jj )
         IF ( j == 0 ) CYCLE

!  Find the entry in row I of this column

         DO k = listvs, listve
           ii = ISVGRP( k )
           i = IW_asmbl( ii )
           IF ( i == 0 .OR. i > j ) CYCLE

!  Skip all elements which lie outside a band of width NSEMIB

           IF ( use_band ) maxsbw = MAX( maxsbw, j - i )
           IF ( j - i > nsemib ) CYCLE
           hesnew = GRAD_el( ii ) * GRAD_el( jj ) * g2dash
           IF ( iprint >= 100 ) WRITE( out, 2090 ) i, j, hesnew
           IF ( no_zeros .AND. hesnew == zero ) CYCLE

!  Obtain the appropriate storage location in H for the new entry

!  Case 1: band matrix storage scheme

           IF ( use_band ) THEN

!  The entry belongs on the diagonal

             IF ( i == j ) THEN
               DIAG( i ) = DIAG( i ) + hesnew

!  The entry belongs off the diagonal

             ELSE
               OFFDIA( j - i, i ) = OFFDIA( j - i, i ) + hesnew
!              H( nfree + j - i + nsemib * ( i - 1 ) ) =                      &
!                H( nfree + j - i + nsemib * ( i - 1 ) ) + hesnew
              END IF

!  Case 2: co-ordinate storage scheme

           ELSE
             istart = j
 150         CONTINUE
             inext = LINK_col( istart )
             IF ( inext == - 1 ) THEN

!  The ( I,J )-th location is empty. Place the new entry in this location
!  and add another link to the list

               nnzh = nnzh + 1
               IF ( nnzh > lh .OR. nnzh > lirnh .OR. nnzh > ljcnh ) THEN
                 nlh = lirnh ; ulh = nnzh - 1 ; mlh = nnzh
!                array = 'IRNH'
!                CALL EXTEND_array( array, lirnh, ulh,                         &
!                                   nlh, mlh, status, iprint, out )
                 CALL EXTEND_arrays( IRNH, lirnh, ulh, nlh, mlh,               &
                                     buffer, status, alloc_status )
                 IF ( status /= 0 ) THEN
                   bad_alloc = 'IRNH' ; GO TO 980 ; END IF
                 lirnh = nlh
                 nlh = ljcnh ; ulh = nnzh - 1 ; mlh = nnzh
!                array = 'JCNH'
!                CALL EXTEND_array( array, lirnh, ulh,                         &
!                                   nlh, mlh, status, iprint, out )
                 CALL EXTEND_arrays( JCNH, lirnh, ulh, nlh, mlh,               &
                                     buffer, status, alloc_status )
                 IF ( status /= 0 ) THEN
                   bad_alloc = 'JCNH' ; GO TO 980 ; END IF
                 ljcnh = nlh
                 nlh = lh ; ulh = nnzh - 1 ; mlh = nnzh
!                array = 'H'
!                CALL EXTEND_array( array, lh, ulh,                            &
!                                   nlh, mlh, status, iprint, out )
                 CALL EXTEND_arrays( H, lh, ulh, nlh, mlh, buffer,             &
                                     status, alloc_status )
                 IF ( status /= 0 ) THEN
                   bad_alloc = 'H' ; GO TO 980 ; END IF
                 lh = nlh
               END IF
               IRNH( nnzh ) = i ; JCNH( nnzh ) = j
               H( nnzh ) = hesnew
               IF ( newpt == llink ) THEN
                 nlh = llink
                 ulh = newpt; mlh = llink + 1
!                array = 'LINK_col'
!                CALL EXTEND_array( array, llink, ulh,                         &
!                                   nlh, mlh, status, iprint, out )
                 CALL EXTEND_arrays( LINK_col, llink, ulh, nlh, mlh,           &
                                     buffer, status, alloc_status )
                 IF ( status /= 0 ) THEN
                   bad_alloc = 'LINK_col' ; GO TO 980 ; END IF
                 llink = nlh
                 nlh = lpos
                 ulh = newpt; mlh = lpos + 1
!                array = 'POS_in_H'
!                CALL EXTEND_array( array, lpos, ulh,                          &
!                                   nlh, mlh, status, iprint, out )
                 CALL EXTEND_arrays( POS_in_H, lpos, ulh, nlh, mlh,            &
                                     buffer, status, alloc_status )
                 IF ( status /= 0 ) THEN
                   bad_alloc = 'POS_in_H' ; GO TO 980 ; END IF
                 lpos = nlh
               END IF
               newpt = newpt + 1
               LINK_col( istart ) = newpt
               POS_in_H( istart ) = nnzh
               LINK_col( newpt  ) = - 1
               POS_in_H( newpt ) = - 1
             ELSE

!  Continue searching the linked list for an entry in row I, column J

               IF ( IRNH( POS_in_H( istart ) ) == i ) THEN
                 ip = POS_in_H( istart )
                 H( ip ) = H( ip ) + hesnew
               ELSE
                 istart = inext
                 GO TO 150
               END IF
             END IF
           END IF
         END DO
       END DO
     END DO

!  Reset the workspace array to zero

     W_el( : maxsel ) = zero

! --------------------------------------------------------
! Add on the low rank first order terms for the I-th group
! --------------------------------------------------------

     DO ig = 1, ng
       IF ( skipg ) THEN
         IF ( KNDOFG( ig ) == 0 ) CYCLE ; END IF
       IF ( .NOT. fixed_structure .AND. GSCALE( ig ) == zero ) CYCLE
       IF ( iprint >= 100 ) WRITE( out, 2100 ) ig
       IF ( GXEQX( ig ) ) THEN
         gdash = GSCALE( ig )
       ELSE
         gdash = GSCALE( ig ) * GVALS2( ig )
         IF ( iprint >= 100 ) WRITE( 6, * ) ' GVALS2( ig )', GVALS2( ig )
         IF ( no_zeros .AND. gdash == zero ) CYCLE
       END IF
       ig1 = ig + 1

!  See if the group has any nonlinear elements

       DO iell = ISTADG( ig ), ISTADG( ig1 ) - 1
         iel = IELING( iell )
         listvs = ISTAEV( iel )
         listve = ISTAEV( iel + 1 ) - 1
         nvarel = listve - listvs + 1
         ielh = ISTADH( iel )
         ihnext = ielh
         scalee = ESCALE( iell )
         DO l = listvs, listve
           j = IW_asmbl( IELVAR( l ) )
           IF ( j /= 0 ) THEN

!  The IEL-th element has an internal representation. Compute the J-th column
!  of the element Hessian matrix

             IF ( INTREP( iel ) ) THEN

!  Compute the J-th column of the Hessian

               W_el( l - listvs + 1 ) = one

!  Find the internal variables

               nin = INTVAR( iel + 1 ) - INTVAR( iel )
               CALL RANGE ( iel, .FALSE., W_el, W_in, nvarel, nin,             &
                            ITYPEE( iel ), nvarel, nin )

!  Multiply the internal variables by the element Hessian

               H_in( : nin ) = zero

!  Only the upper triangle of the element Hessian is stored

                 jcolst = ielh - 1
                 DO jcol = 1, nin
                   ijhess = jcolst
                   jcolst = jcolst + jcol
                   wki = W_in( jcol ) * gdash
                   DO irow = 1, nin
                     IF ( irow <= jcol ) THEN
                       ijhess = ijhess + 1
                     ELSE
                       ijhess = ijhess + irow - 1
                     END IF
                     H_in( irow ) = H_in( irow ) + wki * HUVALS( ijhess )
                   END DO
                 END DO

!  Scatter the product back onto the elemental variables

                CALL RANGE ( iel, .TRUE., H_in, H_el, nvarel, nin,             &
                             ITYPEE( iel ), nin, nvarel )
                W_el( l - listvs + 1 ) = zero
              END IF

!  Find the entry in row I of this column

              DO k = listvs, l
                i = IW_asmbl( IELVAR( k ) )

!  Skip all elements which lie outside a band of width nsemib

                IF ( use_band .AND. i /= 0 )                                   &
                  maxsbw = MAX( maxsbw, ABS( j - i ) )
                IF ( ABS( i - j ) <= nsemib .AND. i /= 0 ) THEN

!  Only the upper triangle of the matrix is stored

                  IF ( i <= j ) THEN
                    ii = i
                    jj = j
                  ELSE
                    ii = j
                    jj = i
                  END IF

!  Obtain the appropriate storage location in H for the new entry

                  IF ( INTREP( iel ) ) THEN
                    hesnew = scalee * H_el( k - listvs + 1 )
                  ELSE
                    hesnew = scalee * HUVALS( ihnext ) * gdash
                  END IF
                  IF ( iprint >= 100 ) WRITE( 6, 2080 ) ii, jj, iel, hesnew

!  Case 1: band matrix storage scheme

                  IF ( use_band ) THEN

!  The entry belongs on the diagonal

                    IF ( ii == jj ) THEN
                      DIAG( ii ) = DIAG( ii ) + hesnew
                      IF ( k /= l ) DIAG( ii ) = DIAG( ii ) + hesnew

!  The entry belongs off the diagonal

                    ELSE
                      OFFDIA( jj - ii, ii ) = OFFDIA( jj - ii, ii ) + hesnew
                    END IF

!  Case 2: co-ordinate storage scheme

                  ELSE
                    IF ( .NOT. no_zeros .OR. hesnew /= zero ) THEN
                     istart = jj
  230                CONTINUE
                     inext = LINK_col( istart )
                     IF ( inext == - 1 ) THEN

!  The ( I,J )-th location is empty. Place the new entry in this location
!  and add another link to the list

                       nnzh = nnzh + 1
                       
                       IF ( nnzh > lh .OR. nnzh > lirnh .OR. nnzh > ljcnh ) THEN
                         nlh = lirnh ; ulh = nnzh - 1; mlh = nnzh
!                        array = 'IRNH'
!                        CALL EXTEND_array( array, lirnh, ulh, nlh, mlh, 
!                                           status, iprint, out, alloc_status )
                         CALL EXTEND_arrays( IRNH, lirnh, ulh, nlh, mlh,       &
                                             buffer, status, alloc_status )
                         IF ( status /= 0 ) THEN
                            bad_alloc = 'IRNH' ; GO TO 980 ; END IF
                         lirnh = nlh
                         nlh = ljcnh ; ulh = nnzh - 1; mlh = nnzh
!                        array = 'JCNH'
!                        CALL EXTEND_array( array, lirnh, ulh, nlh, mlh,       &
!                                           status, iprint, out )
                         CALL EXTEND_arrays( JCNH, lirnh, ulh, nlh, mlh,       &
                                             buffer, status, alloc_status )
                         IF ( status /= 0 ) THEN
                            bad_alloc = 'JCNH' ; GO TO 980 ; END IF
                         ljcnh = nlh
                         nlh = lh ; ulh = nnzh - 1 ; mlh = nnzh
!                        array = 'H'
!                        CALL EXTEND_array( array, lh, ulh, nlh, mlh, status,  &
!                                           iprint, out )
                         CALL EXTEND_arrays( H, lh, ulh, nlh, mlh,             &
                                             buffer, status, alloc_status )
                         IF ( status /= 0 ) THEN
                            bad_alloc = 'H' ; GO TO 980 ; END IF
                         lh = nlh
                       END IF   
                       IRNH( nnzh ) = ii ; JCNH( nnzh ) = jj
                       H( nnzh ) = hesnew
                       IF( k /= l .AND. ii == jj ) H( nnzh ) = hesnew + hesnew
                       IF ( newpt == llink ) THEN
                         nlh = llink
                         ulh = newpt; mlh = llink + 1
!                        array = 'LINK_col'
!                        CALL EXTEND_array( array, llink,  ulh, nlh, mlh,      &
!                                           status, iprint, out )
                         CALL EXTEND_arrays( LINK_col, llink,  ulh, nlh, mlh,  &
                                             buffer, status, alloc_status )
                         IF ( status /= 0 ) THEN
                           bad_alloc = 'LINK_col' ; GO TO 980 ; END IF
                         llink = nlh
                         nlh = lpos
                         ulh = newpt; mlh = lpos + 1
!                        array = 'POS_in_H'
!                        CALL EXTEND_array( array, lpos, ulh, nlh, mlh,        &
!                                           status, iprint, out )
                         CALL EXTEND_arrays( POS_in_H, lpos, ulh, nlh, mlh,    &
                                             buffer, status, alloc_status )
                         IF ( status /= 0 ) THEN
                           bad_alloc = 'POS_in_H' ; GO TO 980 ; END IF
                         lpos = nlh
                       END IF
                       newpt = newpt + 1
                       LINK_col( istart ) = newpt
                       POS_in_H( istart ) = nnzh
                       LINK_col( newpt  ) = - 1
                       POS_in_H( newpt ) = - 1
                     ELSE

! Continue searching the linked list for an entry in row I, column J

                       IF ( IRNH( POS_in_H( istart ) ) == ii ) THEN
                         ip = POS_in_H( istart )
                         H( ip ) = H( ip ) + hesnew
                         IF( k /= l .AND. ii == jj ) H( ip ) = H( ip ) + hesnew
                       ELSE
                         istart = inext
                         GO TO 230
                       END IF
                     END IF
                   END IF
                 END IF
               END IF
               ihnext = ihnext + 1
             END DO
           ELSE
             ihnext = ihnext + l - listvs + 1
           END IF
         END DO
       END DO
     END DO

!  ---------------------------------------
!  For debugging, print the nonzero values
!  ---------------------------------------

     IF ( iprint >= 10 ) THEN
       IF ( .NOT. use_band )                                                   &
         WRITE( out, 2000 ) ( IRNH( i ), JCNH( i ), H( i ), i = 1, nnzh )

!  For debugging, print the nonzero pattern of the matrix


       IF ( nfree <= 36 ) THEN
         MATRIX( : nfree, : nfree ) = '  '
         IF ( use_band ) THEN
           DO i = 1, nfree
             IF ( DIAG( i ) /= zero ) MATRIX( i, i ) = ' *'
             DO j = 1, MIN( nsemib, nfree - i )
               IF ( OFFDIA( j, i ) /= zero ) THEN
                  MATRIX( i + j, i ) = ' *'
                  MATRIX( i, i + j ) = ' *'
               END IF
             END DO
           END DO
         ELSE
           DO i = 1, nnzh
             IF ( IRNH( i ) > nfree ) THEN
               WRITE( out, * ) ' Entry out of bounds in ASMBL ',              &
                                ' row number = ', IRNH( i )
!              STOP
             END IF
             IF ( JCNH( i ) > nfree ) THEN
               WRITE( out, * ) ' Entry out of bounds in ASMBL ',              &
                                ' col number = ', JCNH( i )
!              STOP
             END IF
             MATRIX( IRNH( i ), JCNH( i ) ) = ' *'
             MATRIX( JCNH( i ), IRNH( i ) ) = ' *'
           END DO
         END IF
         WRITE( out, 2040 ) ( i, i = 1, nfree )
         DO i = 1, nfree
           WRITE( out, 2050 ) i, ( MATRIX( i, j ), j = 1, nfree )
         END DO
       END IF
     END IF

!  Successful return

     status = 0
     RETURN

!  Unsuccessful returns

 980 CONTINUE
     WRITE( error, 2990 ) alloc_status, bad_alloc
     RETURN

!  Non-executable statements

 2000  FORMAT( '    Row  Column    Value        Row  Column    Value ', /      &
               '    ---  ------    -----        ---  ------    ----- ', /      &
               ( 2I6, ES24.16, 2I6, ES24.16 ) )
 2040  FORMAT( /, 5X, 36I2 )
 2050  FORMAT( I3, 2X, 36A2 )
 2060  FORMAT( /, I6, ' free variables. They are ', 8I5, /, ( 14I5 ) )
 2070  FORMAT( ' Group ', I5, ' rank-one terms ' )
 2080  FORMAT( ' Row ', I6, ' Column ', I6, ' used from element ', I6,         &
               ' value = ', ES24.16 )
 2090  FORMAT( ' Row ', I6, ' column ', I6, ' used. Value = ', ES24.16 )
 2100  FORMAT( ' Group ', I5, ' second-order terms ' )
 2990  FORMAT( ' ** Message from -LANCELOT_ASMBL_assemble_hessian-', /,        &
               ' Allocation error (status = ', I6, ') for ', A24 )

!  End of subroutine ASMBL_assemble_hessian

     END SUBROUTINE ASMBL_assemble_hessian

!  End of module LANCELOT_ASMBL

   END MODULE LANCELOT_ASMBL_double

