! THIS VERSION: GALAHAD 2.2 - 22/02/2008 AT 12:00 GMT.

!-*-*-*-*-*  G A L A H A D _ C U T E R _ F U N C T I O N S  M O D U L E  *-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal authos: Daniel Robinson and Nick Gould
!
!  History -
!   originally released pre GALAHAD Version 2.2. February 22nd 2008
!
   MODULE GALAHAD_CUTER_FUNCTIONS_double

     USE GALAHAD_SMT_double
     USE GALAHAD_SPACE_double
     USE GALAHAD_NLPT_double, ONLY: NLPT_problem_type, NLPT_userdata_type,     &
                                    NLPT_cleanup
     USE CUTEr_interface_double

     IMPLICIT NONE

!---------------------
!   P r e c i s i o n
!---------------------

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

     PRIVATE
     PUBLIC :: CUTER_eval_F, CUTER_eval_FC, CUTER_eval_C,                      &
               CUTER_eval_G, CUTER_eval_GJ, CUTER_eval_J,                      &
               CUTER_eval_H, CUTER_eval_HPROD, CUTER_eval_JPROD,               &
               CUTER_eval_HL, CUTER_eval_HLPROD,                               &
               CUTER_initialize, CUTER_terminate,                              &
               NLPT_problem_type, NLPT_userdata_type

!------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n
!------------------------------------------------

     TYPE, PUBLIC :: CUTER_FUNCTIONS_control_type
        INTEGER :: input = 5
        INTEGER :: error = 6
        LOGICAL :: separate_linear_constraints = .FALSE.
      END TYPE

     TYPE, PUBLIC :: CUTER_FUNCTIONS_inform_type
        INTEGER :: status, alloc_status
        CHARACTER ( LEN = 80 ) :: bad_alloc
      END TYPE

!----------------------
!   P a r a m e t e r s
!----------------------

     REAL ( KIND = wp ), PARAMETER ::  zero  = 0.0_wp
     REAL ( KIND = wp ), PARAMETER ::  one   = 1.0_wp
     REAL ( KIND = wp ), PARAMETER ::  two   = 2.0_wp
     REAL ( KIND = wp ), PARAMETER ::  ten   = 10.0_wp
     REAL ( KIND = wp ), PARAMETER ::  small = ten ** ( -8 )
     REAL ( KIND = wp ), PARAMETER ::  huge  = ten ** ( 19 )
     
     INTEGER, PARAMETER :: loc_m = 1
     INTEGER, PARAMETER :: loc_n = 2
     INTEGER, PARAMETER :: loc_m_a = 3
     INTEGER, PARAMETER :: loc_nnzh = 4
     INTEGER, PARAMETER :: loc_irnh = 5
     INTEGER, PARAMETER :: loc_icnh = 6
     INTEGER, PARAMETER :: loc_h = 7
     INTEGER, PARAMETER :: loc_nnzj = 8
     INTEGER, PARAMETER :: loc_indfun = 9
     INTEGER, PARAMETER :: loc_indvar = 10
     INTEGER, PARAMETER :: loc_cjac = 11
         
!---------------------------------
!   I n t e r f a c e  b l o c k s
!---------------------------------

!    INTERFACE CUTER_eval_H
!      MODULE PROCEDURE CUTER_eval_H, CUTER_eval_HL
!    END INTERFACE

!    INTERFACE CUTER_eval_HPROD
!      MODULE PROCEDURE CUTER_eval_HPROD, CUTER_eval_HLPROD
!    END INTERFACE

   CONTAINS

!-*-*-*-*-*-*-*-   C U T E R _ e v a l _ F   S U B R O U T I N E  -*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_F( status, f, X, userdata )

!  Evaluate the objective function f(X)

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), INTENT( OUT ) :: f
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! local variables

   INTEGER :: m, n

!  Extract scalar addresses

   m = userdata%integer( loc_m )
   n = userdata%integer( loc_n )

   IF ( m > 0 ) THEN
     CALL CFN( n, m, X, f, m, userdata%real( : m ) )
   ELSE
     CALL UFN( n, X, f )
   END IF
 
   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_F

!-*-*-*-*-*-*-*-   C U T E R _ e v a l _ C   S U B R O U T I N E  -*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_C( status, C, X, userdata )

!  Evaluate the constraint functions C(X)

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: C
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! local variables

   INTEGER :: m, m_a, n
   REAL ( KIND = wp ) :: f

!  Extract scalar addresses

   m   = userdata%integer( loc_m )
   n   = userdata%integer( loc_n )
   m_a = userdata%integer( loc_m_a )

   CALL CFN( n, m, X, f, m, userdata%real( : m ) )
   
   C = userdata%real( m_a + 1 : m )
 
   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_C

!-*-*-*-*-*-*-*-   C U T E R _ e v a l _ F C   S U B R O U T I N E  -*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_FC( status, f, C, X, userdata )

!  Evaluate the objective function f(X) and constraint functions C(X)

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), INTENT( OUT ) :: f
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: C
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! local variables

   INTEGER :: m, m_a, n

!  Extract scalar addresses

   m   = userdata%integer( loc_m )
   m_a = userdata%integer( loc_m_a )
   n   = userdata%integer( loc_n )

   CALL CFN( n, m, X, f, m, userdata%real( : m ) )

   C = userdata%real( m_a + 1 : m )
   
   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_FC

!-*-*-*-*-*-*-*-  C U T E R _ e v a l _ G   S U B R O U T I N E  -*-*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_G( status, G, X, userdata )

!  Evaluate the gradient of the objective function G(X)

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: G
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! Local variables

   INTEGER :: m, n, nnzj, indfun, indvar, cjac, lcjac, l
   REAL ( KIND = wp ), DIMENSION( 1 ) :: Y_dummy = zero

!  Extract scalar and array addresses

   m = userdata%integer( loc_m ) 
   n = userdata%integer( loc_n ) 

   IF ( m > 0 ) THEN
     nnzj = userdata%integer( loc_nnzj ) 
     indfun = userdata%integer( loc_indfun ) 
     indvar = userdata%integer( loc_indvar ) 
     cjac = userdata%integer( loc_cjac ) 

     lcjac = nnzj
     CALL CSGR( n, m, .FALSE., 1, Y_dummy, X, nnzj, lcjac,                      &
                userdata%real( cjac + 1 : cjac + nnzj ),                        &
                userdata%integer( indvar + 1 : indvar + nnzj ),                 &
                userdata%integer( indfun + 1 : indfun + nnzj ) )

! Untangle A: separate the gradient terms from the constraint Jacobian

     G( : n ) = zero
     DO l = 1, nnzj
        IF ( userdata%integer( indfun + l ) == 0 ) THEN
           G( userdata%integer( indvar + l ) ) = userdata%real( cjac + l )
        END IF
     END DO
   ELSE
     CALL UGR( n, X, G )
   END IF

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_G

!-*-*-*-*-*-*-*-*-  C U T E R _ e v a l _ J   S U B R O U T I N E  -*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_J( status, Jval, X, userdata )

!  Evaluate the values of the constraint Jacobian Jval(X) for the nonzeros
!  corresponding to the sparse coordinate format set in CUTER_initialize

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: Jval
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   
! Local variables

   INTEGER :: m, m_a, n, nnzj, indfun, indvar, cjac, Jne, lcjac, l
   REAL ( KIND = wp ), DIMENSION( 1 ) :: Y_dummy = zero

!  Extract scalar and array addresses

   m      = userdata%integer( loc_m )
   m_a    = userdata%integer( loc_m_a )
   n      = userdata%integer( loc_n ) 
   nnzj   = userdata%integer( loc_nnzj ) 
   indfun = userdata%integer( loc_indfun ) 
   indvar = userdata%integer( loc_indvar ) 
   cjac   = userdata%integer( loc_cjac ) 

   lcjac = nnzj
   CALL CSGR( n, m, .FALSE., 1, Y_dummy, X, nnzj, lcjac,                       &
              userdata%real( cjac + 1 : cjac + nnzj ),                         &
              userdata%integer( indvar + 1 : indvar + nnzj ),                  &
              userdata%integer( indfun + 1 : indfun + nnzj ) )

! Untangle A: separate the constraint Jacobian from the objective gradient
!             and the linear constraints.

   Jne = 0
   DO l = 1, nnzj
     IF ( userdata%integer( indfun + l ) > m_a ) THEN
       Jne = Jne + 1
       Jval( Jne ) = userdata%real( cjac + l )
     END IF
   END DO

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_J

!-*-*-*-*-*-*-*-*-  C U T E R _ e v a l _ G J   S U B R O U T I N E  -*-*-*-*-*-

   SUBROUTINE CUTER_eval_GJ( status, G, Jval, X, userdata )

!  Evaluate the gradient of the objective function G(X) and the values
!  of the constraint Jacobian Jval(X) for the nonzeros corresponding to 
!  the sparse coordinate format set in CUTER_initialize

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: G, Jval
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   
! Local variables

   INTEGER :: m, m_a, n, nnzj, indfun, indvar, cjac, Jne, lcjac, l
   REAL ( KIND = wp ), DIMENSION( 1 ) :: Y_dummy = zero

!  Extract scalar and array addresses

   m      = userdata%integer( loc_m )
   m_a    = userdata%integer( loc_m_a )
   n      = userdata%integer( loc_n ) 
   nnzj   = userdata%integer( loc_nnzj ) 
   indfun = userdata%integer( loc_indfun ) 
   indvar = userdata%integer( loc_indvar ) 
   cjac   = userdata%integer( loc_cjac ) 

   lcjac = nnzj
   CALL CSGR( n, m, .FALSE., 1, Y_dummy, X, nnzj, lcjac,                       &
              userdata%real( cjac + 1 : cjac + nnzj ),                         &
              userdata%integer( indvar + 1 : indvar + nnzj ),                  &
              userdata%integer( indfun + 1 : indfun + nnzj ) )

! Untangle A: separate the gradient terms from the constraint Jacobian

   Jne = 0
   G( : n ) = zero
   DO l = 1, nnzj
     IF ( userdata%integer( indfun + l ) == 0 ) THEN
       G( userdata%integer( indvar + l ) ) = userdata%real( cjac + l )
     ELSEIF ( userdata%integer( indfun + l ) > m_a ) then
       Jne = Jne + 1
       Jval( Jne ) = userdata%real( cjac + l )
     END IF
   END DO

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_GJ

!-*-*-*-*-*-*-*-  C U T E R _ e v a l _ H    S U B R O U T I N E  -*-*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_H( status, Hval, X, userdata )
 
!  Evaluate the values of the Herssian of the objective function Hval(X)
!  for the nonzeros corresponding to the sparse coordinate format set in 
!  CUTER_initialize. 

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: Hval
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! Local variables

   INTEGER :: n, nnzh, irnh, icnh, lh

!  Extract scalar and array addresses

   n = userdata%integer( loc_n ) 
   nnzh = userdata%integer( loc_nnzh ) 
   irnh = userdata%integer( loc_irnh ) 
   icnh = userdata%integer( loc_icnh ) 

! Evaluate the Hessian

   lh = nnzh
   CALL USH( n, X, nnzh, lh, Hval,                                             &
             userdata%integer( irnh + 1 : irnh + nnzh ),                       &
             userdata%integer( icnh + 1 : icnh + nnzh ) )                

   status = 0 
   RETURN

   END SUBROUTINE CUTER_eval_H

!-*-*-*-*-*-*-*-  C U T E R _ e v a l _ H L   S U B R O U T I N E  -*-*-*-*-*-*-

   SUBROUTINE CUTER_eval_HL( status, Hval, X, Y, userdata )
 
!  Evaluate the values of the Herssian of the Lagrangian function Hval(X,Y) 
!  for the nonzeros corresponding to the sparse coordinate format set in 
!  CUTER_initialize. By convention, the Lagrangian function is f - sum c_i y_i

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: X, Y
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( OUT ) :: Hval
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata

! Local variables

   INTEGER :: m, m_a, n, nnzh, irnh, icnh, lh
   REAL ( KIND = wp ), DIMENSION( userdata%integer( loc_m ) ) :: Y_full

!  Extract scalar and array addresses

   m    = userdata%integer( loc_m ) 
   m_a  = userdata%integer( loc_m_a ) 
   n    = userdata%integer( loc_n ) 
   nnzh = userdata%integer( loc_nnzh ) 
   irnh = userdata%integer( loc_irnh ) 
   icnh = userdata%integer( loc_icnh ) 

! Evaluate the Hessian

   Y_full                = zero
   Y_full( m_a + 1 : m ) = Y

   lh = nnzh
   CALL CSH( n, m, X, m, - Y_full, nnzh, lh, Hval,                             &
             userdata%integer( irnh + 1 : irnh + nnzh ),                       &
             userdata%integer( icnh + 1 : icnh + nnzh ) )                

   status = 0 
   RETURN

   END SUBROUTINE CUTER_eval_HL

!-*-*-*-*-*  C U T E R _ e v a l _ J P R O D    S U B R O U T I N E  -*-*-*-*-

   SUBROUTINE CUTER_eval_JPROD( status, transpose, U, V, userdata, X )

!  Compute the Jacobian-vector product 
!    U = U + J(X) * V
!  (if transpose is .FALSE.) or
!    U = U + J(X)' * V 
!  (if transpose is .TRUE.). The Jacobian is as recorded from the last
!  point at which it was evaluated if X is absent.

   INTEGER, INTENT( OUT ) :: status
   LOGICAL, INTENT( IN ) :: transpose
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( INOUT ) :: U
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: V
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   REAL ( KIND = wp ), DIMENSION( : ), OPTIONAL, INTENT( IN ) :: X
   
! Local variables

   INTEGER :: m, m_a, n
   REAL( KIND = wp ), DIMENSION( userdata%integer( loc_m ) ) :: full_V

!  Extract scalar and array addresses

   m    = userdata%integer( loc_m )
   m_a  = userdata%integer( loc_m_a )
   n    = userdata%integer( loc_n ) 

   if ( transpose ) then
      full_V = zero
      full_V( m_a + 1 : m ) = V
   end if

   IF ( PRESENT( X ) ) THEN
     IF ( transpose ) THEN
       CALL CJPROD( n, m, .FALSE., transpose, X, full_V, m, userdata%real( : n ), n )
     ELSE
       CALL CJPROD( n, m, .FALSE., transpose, X, V, n, userdata%real( : m ), m )
     END IF
   ELSE
     IF ( transpose ) THEN
       CALL CJPROD( n, m, .TRUE., transpose, userdata%real( : n ), full_V, m,  &
                    userdata%real( : n ), n )
     ELSE
       CALL CJPROD( n, m, .TRUE., transpose, userdata%real( : n ), V, n,       &
                    userdata%real( : m ), m )
     END IF
   END IF
   IF ( transpose ) THEN
     U( : n ) = U( : n ) + userdata%real( : n )
   ELSE
     U( : m - m_a ) = U( : m - m_a ) + userdata%real( m_a + 1 : m )
   END IF

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_JPROD

!-*-*-*-*-  C U T E R _ e v a l _ H P R O D    S U B R O U T I N E  -*-*-*-*-*-

   SUBROUTINE CUTER_eval_HPROD( status, U, V, userdata, X )

!  Compute the product U = U + H(X) * V involving the Hessian of the objective 
!  H(X). If X is absent, the Hessian is as recorded at the last point at which 
!  it was evaluated.

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( INOUT ) :: U
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: V
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   REAL ( KIND = wp ), DIMENSION( : ), OPTIONAL, INTENT( IN ) :: X
   
! Local variables

   INTEGER :: n

!  Extract scalar and array addresses

   n = userdata%integer( loc_n ) 

   IF ( PRESENT( X ) ) THEN
     CALL UPROD( n, .FALSE., X, V, userdata%real( : n ) )
   ELSE
     CALL UPROD( n, .TRUE., userdata%real( : n ), V, userdata%real( : n ) )
   END IF

   U( : n ) = U( : n ) + userdata%real( : n )

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_HPROD

!-*-*-*-*-  C U T E R _ e v a l _ H L P R O D    S U B R O U T I N E  -*-*-*-*-*-

   SUBROUTINE CUTER_eval_HLPROD( status, U, V, userdata, X, Y )

!  Compute the product U = U + H(X,Y) * V involving the Hessian of the
!  Lagrangian H(X,Y). If X and/or Y are absent, the Hessian is as recorded 
!  at the last point at which it was evaluated. By convention, the
!  Lagrangian function is f - sum c_i y_i

   INTEGER, INTENT( OUT ) :: status
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( INOUT ) :: U
   REAL ( KIND = wp ), DIMENSION( : ), INTENT( IN ) :: V
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   REAL ( KIND = wp ), DIMENSION( : ), OPTIONAL, INTENT( IN ) :: X, Y
   
! Local variables

   INTEGER :: m, m_a, n
   REAL ( KIND = wp ), DIMENSION( 1 ) :: Y_dummy = zero
   REAL ( KIND = wp ), DIMENSION( userdata%integer( loc_m ) ) :: full_Y

!  Extract scalar and array addresses

   m   = userdata%integer( loc_m )
   m_a = userdata%integer( loc_m_a )
   n   = userdata%integer( loc_n )

   full_Y = zero
   full_Y( m_a + 1 : m  ) = Y

   IF ( PRESENT( X ) .AND. PRESENT( Y ) ) THEN
     CALL CPROD( n, m, .FALSE., X, m, - full_Y, V, userdata%real( : n ) )
   ELSE
     CALL CPROD( n, m, .TRUE., userdata%real( : n ), 1, Y_dummy, V,     &
                 userdata%real( : n ) )
   END IF
   U( : n ) = U( : n ) + userdata%real( : n )

   status = 0
   RETURN

   END SUBROUTINE CUTER_eval_HLPROD

!-*-*-  C U T E R _ i n i t i a l i z e   S U B R O U T I N E  -*-*-*-*

   SUBROUTINE CUTER_initialize( nlp, control, inform, userdata )
 
   TYPE ( NLPT_problem_type ), INTENT( OUT ) :: nlp
   TYPE ( NLPT_userdata_type ), INTENT( OUT ) :: userdata
   TYPE ( CUTER_FUNCTIONS_control_type ), INTENT( IN ) :: control
   TYPE ( CUTER_FUNCTIONS_inform_type ), INTENT( OUT ) :: inform

! local variables.
     
   INTEGER :: i, j, l, lcjac, lh, status, alloc_status, iend, rend
   INTEGER :: m, n, nnzj, nnzh
   INTEGER :: indfun, indvar, cjac, irnh, icnh, h
   REAL( KIND = wp ) :: f, f2, alpha, alpha_min
   REAL( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y, C_l, C_u, C, X
   REAL( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C2, lin_const
   CHARACTER( LEN = 10 ), ALLOCATABLE, DIMENSION( : ) :: full_CNAMES
   CHARACTER( LEN = 80 ) :: array_name
  
! get dimensions

   CALL CDIMEN( control%input, n, m )
   nlp%n   = n
   nlp%m   = m
   nlp%m_a = 0   ! dummy initialization.

!  Allocate sufficient space to hold the problem
     
   CALL SPACE_resize_array( nlp%n, nlp%X, inform%status, inform%alloc_status )
   IF ( inform%status /= 0 ) THEN 
      inform%bad_alloc = 'nlp%X' ; GO TO 910 ; END IF
      
   CALL SPACE_resize_array( nlp%n, nlp%X_l, inform%status, inform%alloc_status ) 
   IF ( inform%status /= 0 ) THEN 
       inform%bad_alloc = 'nlp%X_l' ; GO TO 910 ; END IF
            
   CALL SPACE_resize_array( nlp%n, nlp%X_u, inform%status, inform%alloc_status )
   IF ( inform%status /= 0 ) THEN 
       inform%bad_alloc = 'nlp%X_u' ; GO TO 910 ; END IF

   CALL SPACE_resize_array( nlp%n, nlp%Z, inform%status, inform%alloc_status )
   IF ( inform%status /= 0 ) THEN 
      inform%bad_alloc = 'nlp%Z' ; GO TO 910 ; END IF
      
   CALL SPACE_resize_array( nlp%n, nlp%G, inform%status, inform%alloc_status )
   IF ( inform%status /= 0 ) THEN 
      inform%bad_alloc = 'nlp%G' ; GO TO 910 ; END IF
      
   CALL SPACE_resize_array( nlp%n, nlp%VNAMES, inform%status, inform%alloc_status )
   IF ( inform%status /= 0 ) THEN 
       inform%bad_alloc = 'nlp%VNAMES' ; GO TO 910 ; END IF

!  -- The problem is constrained --

   IF ( m > 0 ) THEN

     CALL SPACE_resize_array( m, nlp%EQUATION, inform%status,inform%alloc_status)
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%EQUATION' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( m, nlp%LINEAR, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%LINEAR' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( m, Y, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Y' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( m, C_l, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Y' ; GO TO 910 ; END IF
     
     CALL SPACE_resize_array( m, C_u, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Y' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( m, full_CNAMES, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'full_CNAMES' ; GO TO 910 ; END IF  
     
     ! Get the problem with linear constraints possibly first.

     i = n ; j = m
     CALL CSETUP( control%input, control%error, n, m, nlp%X, nlp%X_l, nlp%X_u, i,&
                  nlp%EQUATION, nlp%LINEAR,                                      &
                  Y, C_l, C_u, j, .FALSE., control%separate_linear_constraints, .FALSE. )

     nlp%m_a = 0
     if ( control%separate_linear_constraints ) then
        do l = 1, m
           if ( nlp%LINEAR( l ) ) then
              nlp%m_a = nlp%m_a + 1
           else
              exit
           end if
        end do
     end IF

     nlp%m = m - nlp%m_a

     CALL SPACE_resize_array( nlp%m_a, nlp%Y_a, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Y_a' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m_a, nlp%A_l, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%A_l' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m_a, nlp%Ax, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Ax' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m_a, nlp%A_u, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%A_u' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m_a, nlp%ANAMES, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%ANAMES' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m, nlp%Y, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Y' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m, nlp%C, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%C' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m, nlp%C_l, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%C_l' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m, nlp%C_u, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%C_u' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( nlp%m, nlp%CNAMES, inform%status, inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%CNAMES' ; GO TO 910 ; END IF

!  Obtain the names of the problem, its variables and general constraints

     CALL CNAMES( n, m, nlp%pname, nlp%VNAMES, full_CNAMES )

!  Define the "corrected" separated vectors.

     nlp%Y_a = Y  ( 1 : nlp%m_a ) ;      nlp%Y   = Y  ( nlp%m_a + 1 : m )
     nlp%A_l = C_l( 1 : nlp%m_a ) ;      nlp%C_l = C_l( nlp%m_a + 1 : m )
     nlp%A_u = C_u( 1 : nlp%m_a ) ;      nlp%C_u = C_u( nlp%m_a + 1 : m )

     nlp%ANAMES = full_CNAMES( 1 : nlp%m_a )
     nlp%CNAMES = full_CNAMES( nlp%m_a + 1 : m )

!  Deallocate arrays no longer needed.

     array_name = 'cuter_functions : C_l'
     CALL SPACE_dealloc_array( C_l,                                  &
        inform%status, inform%alloc_status, array_name = array_name,           &
        bad_alloc = inform%bad_alloc, out = control%error )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'C_l' ; GO TO 920 ; END IF

     array_name = 'cuter_functions : C_u'
     CALL SPACE_dealloc_array( C_u,                                  &
        inform%status, inform%alloc_status, array_name = array_name,           &
        bad_alloc = inform%bad_alloc, out = control%error )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'C_u' ; GO TO 920 ; END IF

     array_name = 'cuter_functions : full_CNAMES'
     CALL SPACE_dealloc_array( full_CNAMES,                                  &
        inform%status, inform%alloc_status, array_name = array_name,           &
        bad_alloc = inform%bad_alloc, out = control%error )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'full_CNAMES' ; GO TO 920 ; END IF

! Set up sparsity structure for A, J, and H.  ( Assume co-ordinate storage )

! Determine number of non-zeros in the matrix of gradients of the
! objective function AND constraint functions.

     CALL CDIMSJ( nnzj )

!  Determine how many nonzeros are required to store the Hessian matrix of the
!  Lagrangian, when the matrix is stored as a sparse matrix in "co-ordinate"
!  format (only the lower triangular part is stored).

     CALL CDIMSH( nnzh )
  
!  Set starting addresses for workspace array partitions 

     irnh    = loc_cjac
     icnh    = irnh + nnzh
     indfun  = icnh + nnzh
     indvar  = indfun + nnzj
     iend    = indvar + nnzj
     
     h    = 0
     cjac = h + nnzh
     rend = MAX( cjac + nnzj, m, n )

! Allocate space to hold scalars/arrays needed for subsequent calls

     CALL SPACE_resize_array( iend, userdata%integer, inform%status,           &
                              inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'userdata%integer' ; GO TO 910 ; END IF

     CALL SPACE_resize_array( rend, userdata%real, inform%status,              &
                              inform%alloc_status )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'userdata%real' ; GO TO 910 ; END IF

! Record workspace partitions in userdata%integer.

     userdata%integer( loc_m )       = m
     userdata%integer( loc_n )       = n
     userdata%integer( loc_m_a )     = nlp%m_a
     userdata%integer( loc_nnzh )    = nnzh
     userdata%integer( loc_irnh )    = irnh
     userdata%integer( loc_icnh )    = icnh
     userdata%integer( loc_h )       = h
     userdata%integer( loc_nnzj )    = nnzj
     userdata%integer( loc_indfun )  = indfun
     userdata%integer( loc_indvar )  = indvar
     userdata%integer( loc_cjac )    = cjac

! Determine if there is a constant in the linear constraints.  Adjust the bounds if necessary.

     if ( nlp%m_a > 0 ) then

        CALL SPACE_resize_array( m, C, inform%status, inform%alloc_status )
        IF ( inform%status /= 0 ) then
           inform%bad_alloc = 'userdata%real' ; GO TO 910 ; END IF

        CALL SPACE_resize_array( m, C2, inform%status, inform%alloc_status )
        IF ( inform%status /= 0 ) then
           inform%bad_alloc = 'userdata%real' ; GO TO 910 ; END IF

        CALL SPACE_resize_array( nlp%m_a, lin_const, inform%status, inform%alloc_status )
        IF ( inform%status /= 0 ) then
           inform%bad_alloc = 'userdata%real' ; GO TO 910 ; END IF

        CALL SPACE_resize_array( n, X, inform%status, inform%alloc_status )
        IF ( inform%status /= 0 ) THEN 
            inform%bad_alloc = 'userdata%real' ; GO TO 910 ; END IF

        ! Make X feasible with respect to bounds.

        X = zero
        X = min( nlp%X_u, max( nlp%X_l, X) )

        alpha_min = two
        do i = 1, n
           if ( X(i) > zero ) then
              alpha = min( one + one/X(i), nlp%X_u(i)/X(i) )
           elseif ( X(i) < zero ) then
              alpha = min( one - one/X(i), nlp%X_l(i)/X(i) )
           else
              alpha = two
           end if
           alpha_min = min( alpha_min, alpha )
        end do
        alpha_min = max( alpha_min, 1.001_wp )
        
        call CFN( n, m, X, f, m, C )
        call CFN( n, m, alpha_min*X, f2, m, C2 )

        lin_const = alpha_min * C( : nlp%m_a ) - C2( : nlp%m_a )
        lin_const = lin_const / (alpha_min - one )
        
        do i = 1, nlp%m_a
           if ( nlp%A_l( i ) > - huge ) then
              nlp%A_l( i ) = nlp%A_l( i ) - lin_const( i )
           end if
           if ( nlp%A_u( i ) < huge ) then
              nlp%A_u( i ) = nlp%A_u( i ) - lin_const( i )
           end if
        end do
        
        array_name = 'cuter_functions : C'
        CALL SPACE_dealloc_array( C,                                           &
             inform%status, inform%alloc_status, array_name = array_name,      &
             bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) THEN 
           inform%bad_alloc = 'C' ; GO TO 920
        END IF

        array_name = 'cuter_functions : X'
        CALL SPACE_dealloc_array( X,                                           &
             inform%status, inform%alloc_status, array_name = array_name,      &
             bad_alloc = inform%bad_alloc, out = control%error )
        IF ( inform%status /= 0 ) THEN 
           inform%bad_alloc = 'X' ; GO TO 920
        END IF

     end IF

! Evaluate the Jacobian and Hessian and get sparsity pattern.

     lcjac = nnzj
     lh    = nnzh
     CALL CSGRSH( nlp%n, nlp%m, nlp%X, .FALSE., m, - Y, nnzj, lcjac,            &
                  userdata%real( cjac + 1 : cjac + nnzj ),                      &
                  userdata%integer( indvar + 1 : indvar + nnzj ),               &
                  userdata%integer( indfun + 1 : indfun + nnzj ),               &
                  nnzh, lh,                                                     &
                  userdata%real( h + 1 : h + nnzh ),                            &
                  userdata%integer( irnh + 1 : irnh + nnzh ),                   &
                  userdata%integer( icnh + 1 : icnh + nnzh ) )              

! get number of nonzeros in the linear constraints and Jacobian constraints only.

     nlp%J%ne = 0
     nlp%A%ne = 0
     DO l = 1, nnzj
        IF ( userdata%integer( indfun + l ) == 0 ) THEN
           ! Relax....objective gradient component.
        ELSEIF ( userdata%integer( indfun + l ) <= nlp%m_a ) THEN
           nlp%A%ne = nlp%A%ne + 1
        ELSE
           nlp%J%ne = nlp%J%ne + 1
        END IF
     END DO

!  Deallocate arrays no longer needed.

     array_name = 'cuter_functions : Y'
     CALL SPACE_dealloc_array( Y,                                  &
        inform%status, inform%alloc_status, array_name = array_name,           &
        bad_alloc = inform%bad_alloc, out = control%error )
     IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'Y' ; GO TO 920 ; END IF

!  Allocate arrays that are now of correct length.

     CALL SPACE_resize_array( nlp%A%ne, nlp%A%row, status, alloc_status )
     IF ( status /= 0 ) GO TO 990
     
     CALL SPACE_resize_array( nlp%A%ne, nlp%A%col, status, alloc_status )
     IF ( status /= 0 ) GO TO 990
     
     CALL SPACE_resize_array( nlp%A%ne, nlp%A%val, status, alloc_status )
     IF ( status /= 0 ) GO TO 990

     CALL SPACE_resize_array( nlp%J%ne, nlp%J%row, status, alloc_status )
     IF ( status /= 0 ) GO TO 990
     
     CALL SPACE_resize_array( nlp%J%ne, nlp%J%col, status, alloc_status )
     IF ( status /= 0 ) GO TO 990
     
     CALL SPACE_resize_array( nlp%J%ne, nlp%J%val, status, alloc_status )
     IF ( status /= 0 ) GO TO 990
   
!  Untangle J: separate the gradient terms from the linear constraints
!              and the general constraints in the Jacobian
     
     nlp%J%ne = 0
     nlp%A%ne = 0

     nlp%A%n = n
     nlp%A%m = nlp%m_a

     nlp%J%n = n
     nlp%J%m = m - nlp%m_a

     DO l = 1, nnzj
        IF ( userdata%integer( indfun + l ) == 0 ) THEN
           ! Relax....objective gradient component.
        ELSEIF ( userdata%integer( indfun + l ) <= nlp%m_a ) THEN
           nlp%A%ne = nlp%A%ne + 1
           nlp%A%row( nlp%A%ne ) = userdata%integer( indfun + l )
           nlp%A%col( nlp%A%ne ) = userdata%integer( indvar + l )
           nlp%A%val( nlp%A%ne ) = userdata%real( cjac + l )
        ELSE
           nlp%J%ne = nlp%J%ne + 1
           nlp%J%row( nlp%J%ne ) = userdata%integer( indfun + l ) - nlp%m_a
           nlp%J%col( nlp%J%ne ) = userdata%integer( indvar + l )
        END IF
     END DO

!  Define the storage type for J

     CALL SMT_put( nlp%A%type, 'COORDINATE', status )
     IF ( status /= 0 ) GO TO 992

     CALL SMT_put( nlp%J%type, 'COORDINATE', status )
     IF ( status /= 0 ) GO TO 992 

   ELSE

!  -- The problem is unconstrained --

!  Set up the correct data structures for subsequent computations

      i = n
      CALL USETUP( control%input, control%error, n, nlp%X, nlp%X_l, nlp%X_u, i )

!  Obtain the names of the problem and its variables

      CALL UNAMES( n, nlp%pname, nlp%VNAMES )

! Set up sparsity structure for H.  ( Assume co-ordinate storage )

!  Determine how many nonzeros are required to store the Hessian matrix
!  when the matrix is stored as a sparse matrix in "co-ordinate" format 
!  (only the lower triangular part is stored).

      CALL CDIMSH( nnzh )
  
!  Set starting addresses for workspace array partitions 

      irnh = loc_h
      icnh = irnh + nnzh
      iend = icnh + nnzh
      h = 0
      rend = MAX( h + nnzh, n )

! Allocate space to hold scalars/arrays needed for subsequent calls

      CALL SPACE_resize_array( iend, userdata%integer, inform%status,           &
                              inform%alloc_status )
      IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'userdata%integer' ; GO TO 910
      END IF

      CALL SPACE_resize_array( rend, userdata%real, inform%status,              &
                               inform%alloc_status )
      IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'userdata%real' ; GO TO 910
      END IF

! Record workspace partitions in userdata%integer

      userdata%integer( loc_m ) = m
      userdata%integer( loc_n ) = n
      userdata%integer( loc_nnzh ) = nnzh
      userdata%integer( loc_irnh ) = irnh
      userdata%integer( loc_icnh ) = icnh
      userdata%integer( loc_h ) = h

! Evaluate the Jacobian and Hessian and get sparsity pattern

      lh = nnzh
      CALL USH( nlp%n, nlp%X, nnzh, lh,                                        &
                userdata%real( h + 1 : h + nnzh ),                             &
                userdata%integer( irnh + 1 : irnh + nnzh ),                    &
                userdata%integer( icnh + 1 : icnh + nnzh ) )                


! DPR: Added to prevent errors.

      CALL SPACE_resize_array( 0, nlp%C, inform%status, inform%alloc_status )
      IF ( inform%status /= 0 ) then
         inform%bad_alloc = 'nlp%C' ; GO TO 910
      END IF

      CALL SPACE_resize_array( 0, nlp%Y, inform%status, inform%alloc_status )
      IF ( inform%status /= 0 ) then
         inform%bad_alloc = 'nlp%Y' ; GO TO 910
      END IF

      CALL SPACE_resize_array( 0, nlp%Ax, inform%status, inform%alloc_status )
      IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Ax' ; GO TO 910
      END IF
 
      CALL SPACE_resize_array( 0, nlp%Y_a, inform%status, inform%alloc_status )
      IF ( inform%status /= 0 ) THEN 
         inform%bad_alloc = 'nlp%Ax' ; GO TO 910
      END IF


   END IF

! Define reduced costs

   nlp%Z = zero

! Define the storage type for H

   CALL SMT_put( nlp%H%type, 'COORDINATE', status )
   IF ( status /= 0 ) GO TO 992

! Set the sparsity pattern for H

   nlp%H%n = n
   nlp%H%ne = nnzh

   CALL SPACE_resize_array( nlp%H%ne, nlp%H%row, status, alloc_status )
   IF ( status /= 0 ) GO TO 990
   
   CALL SPACE_resize_array( nlp%H%ne, nlp%H%col, status, alloc_status )
   IF ( status /= 0 ) GO TO 990
   
   CALL SPACE_resize_array( nlp%H%ne, nlp%H%val, status, alloc_status )
   IF ( status /= 0 ) GO TO 990  

!  Ensure that the lower triangle of H is stored

   DO l = 1, nlp%H%ne
     i = userdata%integer( irnh + l )
     j = userdata%integer( icnh + l )
     IF ( i > j ) THEN
       nlp%H%row( l ) = i
       nlp%H%col( l ) = j
     ELSE
       nlp%H%row( l ) = j
       nlp%H%col( l ) = i
     END IF
   END DO
   RETURN

! Error format statements

910 CONTINUE
   WRITE( control%error, "( ' ** ERROR - subroutine cuter_functions - ',       &
        &         'Allocation error (status = ', I6, ') for ', A24, / )" )     &
        inform%status, inform%bad_alloc
   RETURN

920 CONTINUE
   WRITE( control%error, "( ' ** ERROR - subroutine cuter_functions - ',       &
        &         'Deallocation error (status = ', I6, ') for ', A24, / )" )   &
        inform%status, inform%bad_alloc
   RETURN
   
990 CONTINUE
   inform%status = status
   WRITE( control%error, "( ' CUTER_initialize: deallocation error. ', I0,     &
  &   'status ', I0 ) " ) status, alloc_status
   RETURN
   
992 CONTINUE
   inform%status = status
   WRITE( control%error,"(' CUTER_initialize: error using subroutine SMT_put', &
  & ' - status= ', I0 ) " )  status, alloc_status
   
   RETURN

   END SUBROUTINE CUTER_initialize

!-*-*-*-*-*-   C U T E R  _ t e r m i n a t e   S U B R O U T I N E   -*-*-*-*-

   SUBROUTINE CUTER_terminate( nlp, inform, userdata )

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!      ..............................................
!      .                                            .
!      .  Deallocate internal arrays at the end     .
!      .  of the computation                        .
!      .                                            .
!      ..............................................

! =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

!  Dummy arguments

   TYPE ( NLPT_problem_type ), INTENT( INOUT ) :: nlp
   TYPE ( NLPT_userdata_type ), INTENT( INOUT ) :: userdata
   TYPE ( CUTER_FUNCTIONS_inform_type ), INTENT( OUT ) :: inform

   CALL NLPT_cleanup( nlp )
   CALL SPACE_dealloc_array( userdata%integer, inform%status,                 &
                             inform%alloc_status )
   CALL SPACE_dealloc_array( userdata%real, inform%status,                    &
                             inform%alloc_status )

   RETURN

!  End of subroutine CUTER_terminate

   END SUBROUTINE CUTER_terminate

END MODULE GALAHAD_CUTER_FUNCTIONS_double

