! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.

!-*-*-*-*-*-*-*-*-  G A L A H A D   R U N P R E _ S I F  *-*-*-*-*-*-*-*-*-*-

!  Nick Gould, Dominique Orban and Philippe Toint, for GALAHAD productions
!  Copyright reserved
!  March 14th 2003

   PROGRAM RUNPRE_SIF
   USE GALAHAD_USEPRE_double

!  Main program for the SIF/CUTEr interface to PRE, a preprocessing
!  algorithm for quadratic programs

!  Problem input characteristics

   INTEGER, PARAMETER :: input = 55
   CHARACTER ( LEN = 16 ) :: prbdat = 'OUTSDIF.d'

!  Open the data input file

   OPEN( input, FILE = prbdat, FORM = 'FORMATTED', STATUS = 'OLD'  )
   REWIND input

!  Call the CUTEr interface

   CALL USE_PRE( input )

!  Close the data input file 

   CLOSE( input  )
   STOP

!  End of RUNPRE_SIF

   END PROGRAM RUNPRE_SIF
