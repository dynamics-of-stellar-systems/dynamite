! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-  G A L A H A D _ S T R I N G   M O D U L E  *-*-*-*-*-*-*-*-*

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released GALAHAD Version 2.0. September 14th 2005

!  For full documentation, see 
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_STRING_double

!     ----------------------------
!    |  Set strings appropriate   |
!    !  for singular and pleural  |
!    |  forms of words along with |
!    |  other useful strings      |
!     ----------------------------

     IMPLICIT NONE     

     PRIVATE
     PUBLIC :: STRING_pleural, STRING_are, STRING_have, STRING_their,         &
               STRING_sign, STRING_choice

!  Set precision

     INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!  Set other parameters

     REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp

   CONTAINS

!-*-*-*-  G A L A H A D -  S T R I N G _ p l e u r a l   F U N C T I O N  -*-*-*-

     FUNCTION STRING_pleural( val )

!   Given an integer val, returns "s" if v /= 0, otherwise returns " "

     CHARACTER ( len = 1 ) :: STRING_pleural

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     INTEGER, INTENT( IN ) :: val

     IF ( val /= 1 ) THEN
       STRING_pleural = "s"
     ELSE
       STRING_pleural = " "
     END IF

     RETURN

!  End of function STRING_pleural

      END FUNCTION STRING_pleural

!-*-*-*-*-  G A L A H A D -  S T R I N G _ a r e   F U N C T I O N  -*-*-*-*-

     FUNCTION STRING_are( val )

!   Given an integer val, returns "s" if v /= 0, otherwise returns " "

     CHARACTER ( len = 3 ) :: STRING_are

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     INTEGER, INTENT( IN ) :: val

     IF ( val /= 1 ) THEN
       STRING_are = "are"
     ELSE
       STRING_are = "is "
     END IF

     RETURN

!  End of function STRING_are

      END FUNCTION STRING_are

!-*-*-*-*-  G A L A H A D -  S T R I N G _ h a v e   F U N C T I O N  -*-*-*-*-

     FUNCTION STRING_have( val )

!   Given an integer val, returns "s" if v /= 0, otherwise returns " "

     CHARACTER ( len = 4 ) :: STRING_have

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     INTEGER, INTENT( IN ) :: val

     IF ( val /= 1 ) THEN
       STRING_have = "have"
     ELSE
       STRING_have = "has "
     END IF

     RETURN

!  End of function STRING_have

      END FUNCTION STRING_have

!-*-*-*-  G A L A H A D -  S T R I N G _ t h e i r    F U N C T I O N  -*-*-*-

     FUNCTION STRING_their( val )

!   Given an integer val, returns "their" if v /= 1, otherwise returns "its"

     CHARACTER ( len = 5 ) :: STRING_their

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     INTEGER, INTENT( IN ) :: val

     IF ( val /= 1 ) THEN
       STRING_their = "their"
     ELSE
       STRING_their = "its  "
     END IF

     RETURN

!  End of function STRING_their

      END FUNCTION STRING_their

!-*-*-*-  G A L A H A D -  S T R I N G _ c h o i c e    F U N C T I O N  -*-*-*-

     FUNCTION STRING_choice( val, string1, string2 )

!   Given an integer val, returns string1 if v /= 1, otherwise returns string2

     CHARACTER ( len = 120 ) :: STRING_choice
 
!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     INTEGER, INTENT( IN ) :: val
     CHARACTER ( len = * ), INTENT( IN ) :: string1, string2

     IF ( val /= 1 ) THEN
       STRING_choice = string1
     ELSE
       STRING_choice = string2
     END IF

     RETURN

!  End of function STRING_choice

      END FUNCTION STRING_choice

!-*-*-*-*-  G A L A H A D -  S T R I N G _ s i g n   F U N C T I O N  -*-*-*-*-

     FUNCTION STRING_sign( val, show_plus )

!   Given a real number val, returns " " (or "+" if show_plus is true)
!   if val >= 0 and "-" if val < 0

     CHARACTER ( len = 1 ) :: STRING_sign

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------

     REAL ( KIND = wp ), INTENT( IN ) :: val
     LOGICAL, INTENT( IN ) :: show_plus

     IF ( val < zero ) THEN
       STRING_sign = "-"
     ELSE
       IF ( show_plus ) THEN
         STRING_sign = "+"
       ELSE
         STRING_sign = " "
       END IF
     END IF

     RETURN

!  End of function STRING_sign

      END FUNCTION STRING_sign

!  End of module GALAHAD_STRING

   END MODULE GALAHAD_STRING_double
