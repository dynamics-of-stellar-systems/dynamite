! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-  G A L A H A D   R U N L A N B _ S I F  *-*-*-*-*-*-*-*-*-*-

!  Nick Gould, Dominique Orban and Philippe Toint, for GALAHAD productions
!  Copyright reserved
!  March 14th 2003

   PROGRAM RUNLANB_SIF
   USE GALAHAD_USELANB_double

!  Main program for the SIF interface to LANCELOT B, an augmented Lagrangian
!  algorithm for nonlinear programming

!  Problem input characteristics

   INTEGER, PARAMETER :: input = 55
   CHARACTER ( LEN = 16 ) :: prbdat = 'OUTSDIF.d'

!  Open the data input file

   OPEN( input, FILE = prbdat, FORM = 'FORMATTED', STATUS = 'OLD'  )
   REWIND input

!  Call the SIF interface

   CALL USE_LANB( input )

!  Close the data input file 

   CLOSE( input  )
   STOP

!  End of RUNLANB_SIF

   END PROGRAM RUNLANB_SIF
