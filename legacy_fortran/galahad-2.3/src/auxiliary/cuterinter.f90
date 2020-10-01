! THIS VERSION: GALAHAD 2.1 - 22/03/2007 AT 09:00 GMT.
!  Updated 06/02/2004: missing interfaces to CCFSG and CCIFSG added, and

!-*-*-*-*-*-  CUTEr_interface M O D U L E  *-*-*-*-*-*-*-*

!  Nick Gould, for CUTEr productions
!  Copyright reserved
!  November 10th 2000

   MODULE CUTEr_interface_double

!  Interface blocks for CUTE fortran tools

     IMPLICIT NONE

     INTERFACE

!  Interface block for unconstrained tools

       SUBROUTINE UDIMEN( INPUT, N )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: INPUT
       INTEGER, INTENT( OUT ) :: N
       END SUBROUTINE UDIMEN
       
       SUBROUTINE USETUP( INPUT , IOUT  , N , X , BL, BU, NMAX )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: INPUT, IOUT, NMAX
       INTEGER, INTENT( OUT ) :: N
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( NMAX ) :: X, BL, BU
       END SUBROUTINE USETUP
       
       SUBROUTINE UNAMES( N, PNAME, VNAME )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N
       CHARACTER ( LEN = 10 ), INTENT( OUT ) :: PNAME
       CHARACTER ( LEN = 10 ), INTENT( OUT ), DIMENSION( N ) :: VNAME
       END SUBROUTINE UNAMES
       
       SUBROUTINE UVARTY( N, IVARTY )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N
       INTEGER, INTENT( OUT ) :: IVARTY( N )
       END SUBROUTINE UVARTY
       
       SUBROUTINE UFN( N, X, F )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N
       REAL ( KIND = wp ), INTENT( OUT ) :: F
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       END SUBROUTINE UFN
       
       SUBROUTINE UGR( N, X, G )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       END SUBROUTINE UGR
       
       SUBROUTINE UOFG( N, X, F, G, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N
       REAL ( KIND = wp ), INTENT( OUT ) :: F
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       LOGICAL, INTENT( IN ) :: GRAD
       END SUBROUTINE UOFG
       
       SUBROUTINE UDH( N, X, LH1, H )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, LH1
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH1, N ) :: H
       END
       
       SUBROUTINE UGRDH( N, X, G, LH1, H )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, LH1
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH1, N ) :: H
       END SUBROUTINE UGRDH
       
       SUBROUTINE UDIMSH( NNZH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: NNZH
       END SUBROUTINE UDIMSH
       
       SUBROUTINE USH( N, X, NNZH, LH, H, IRNH, ICNH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, LH
       INTEGER, INTENT( OUT ) :: NNZH
       INTEGER, INTENT( OUT ), DIMENSION( LH ) :: IRNH, ICNH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH ) :: H
       END SUBROUTINE USH
       
       SUBROUTINE UDIMSE( NE, NNZH, NZIRNH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: NE, NNZH, NZIRNH
       END SUBROUTINE UDIMSE
       
       SUBROUTINE UEH( N , X , NE    , IRNHI , LIRNHI, LE    ,                 &
                       IPRNHI, HI    , LHI   , IPRHI , BYROWS )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , LE    , LIRNHI, LHI 
       INTEGER, INTENT( IN ) :: NE
       LOGICAL, INTENT( IN ) :: BYROWS
       INTEGER, INTENT( OUT ), DIMENSION( LIRNHI ) :: IRNHI
       INTEGER, INTENT( OUT ), DIMENSION( LE ) :: IPRNHI, IPRHI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LHI ) :: HI
       END SUBROUTINE UEH
       
       SUBROUTINE UGRSH( N, X, G, NNZH, LH, H, IRNH, ICNH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, LH
       INTEGER, INTENT( OUT ) :: NNZH
       INTEGER, INTENT( OUT ), DIMENSION( LH ) :: IRNH, ICNH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH ) :: H
       END SUBROUTINE UGRSH
       
       SUBROUTINE UGREH( N , X , G, NE , IRNHI , LIRNHI, LE    ,               &
                         IPRNHI, HI    , LHI   , IPRHI , BYROWS )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, LE, LIRNHI, LHI 
       INTEGER, INTENT( OUT ) :: NE
       LOGICAL, INTENT( IN ) :: BYROWS
       INTEGER, INTENT( OUT ), DIMENSION( LIRNHI ) :: IRNHI
       INTEGER, INTENT( OUT ), DIMENSION( LE ) :: IPRNHI, IPRHI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LHI ) :: HI
       END SUBROUTINE UGREH
       
       SUBROUTINE UPROD( N, GOTH, X, P, RESULT )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N
       LOGICAL, INTENT( IN ) :: GOTH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X, P
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: RESULT
       END SUBROUTINE UPROD
       
       SUBROUTINE UBANDH( N, GOTH, X, NSEMIB, BANDH, LBANDH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, NSEMIB, LBANDH
       LOGICAL, INTENT( IN ) ::  GOTH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) ::  X
       REAL ( KIND = wp ), INTENT( OUT ),                                      &
                                DIMENSION( 0 : LBANDH, N ) ::  BANDH
       END SUBROUTINE UBANDH
   
!  Interface block for constrained tools

       SUBROUTINE CDIMEN( INPUT, N, M )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: INPUT
       INTEGER, INTENT( OUT ) :: N, M
       END SUBROUTINE CDIMEN
       
       SUBROUTINE CSETUP( INPUT , IOUT  , N , M , X , BL , BU   , NMAX,        &
                          EQUATN, LINEAR, V , CL, CU     , MMAX , EFIRST,      &
                          LFIRST, NVFRST )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  INPUT , IOUT  , NMAX  , MMAX
       INTEGER, INTENT( OUT ) ::  N, M
       LOGICAL, INTENT( IN ) ::  EFIRST, LFIRST, NVFRST
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( NMAX ) :: X, BL, BU
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( MMAX ) :: V, CL, CU
       LOGICAL, INTENT( OUT ), DIMENSION( MMAX ) :: EQUATN, LINEAR
       END SUBROUTINE CSETUP
       
       SUBROUTINE CNAMES( N, M, PNAME, VNAME, GNAME )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M
       CHARACTER ( LEN = 10 ), INTENT( OUT ) :: PNAME
       CHARACTER ( LEN = 10 ), INTENT( OUT ), DIMENSION( N ) :: VNAME
       CHARACTER ( LEN = 10 ), INTENT( OUT ), DIMENSION( M ) :: GNAME
       END SUBROUTINE CNAMES
       
       SUBROUTINE CVARTY( N, IVARTY )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N
       INTEGER, INTENT( OUT ) :: IVARTY( N )
       END SUBROUTINE CVARTY
       
       SUBROUTINE CFN( N , M , X , F , LC, C )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , M , LC
       REAL ( KIND = wp ), INTENT( OUT ) :: F
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LC ) :: C
       END SUBROUTINE CFN
       
       SUBROUTINE CGR( N , M , X     , GRLAGF, LV, V , G     , JTRANS,         &
                       LCJAC1, LCJAC2, CJAC  )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , M , LV    , LCJAC1, LCJAC2
       LOGICAL, INTENT( IN ) :: GRLAGF, JTRANS
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       REAL ( KIND = wp ), INTENT( OUT ),                                      &
                                DIMENSION( LCJAC1, LCJAC2 ) :: CJAC
       END SUBROUTINE CGR
       
       SUBROUTINE COFG( N, X, F, G, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N
       REAL ( KIND = wp ), INTENT( OUT ) :: F
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       LOGICAL, INTENT( IN ) :: GRAD
       END SUBROUTINE COFG
       
       SUBROUTINE CDIMSJ( NNZJ )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: NNZJ
       END SUBROUTINE CDIMSJ
       
       SUBROUTINE CSGR( N , M , GRLAGF, LV, V , X     , NNZJ  ,                &
                        LCJAC , CJAC  , INDVAR, INDFUN )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , M , LV, LCJAC
       INTEGER, INTENT( OUT ) :: NNZJ
       LOGICAL, INTENT( IN ) :: GRLAGF
       INTEGER, INTENT( OUT ), DIMENSION( LCJAC ) :: INDVAR, INDFUN
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LCJAC ) :: CJAC
       END SUBROUTINE CSGR
       
       SUBROUTINE CCFG( N, M, X, LC, C, JTRANS, LCJAC1, LCJAC2, CJAC, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LC, LCJAC1, LCJAC2
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LC ) :: C
       REAL ( KIND = wp ), INTENT( OUT ),                                      &
                                DIMENSION( LCJAC1, LCJAC2 ) :: CJAC
       LOGICAL, INTENT( IN ) :: JTRANS, GRAD
       END SUBROUTINE CCFG
       
       SUBROUTINE CSCFG( N, M, X, LC, C, NNZJ, LCJAC, CJAC,                    &
                         INDVAR, INDFUN, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LC, LCJAC
       INTEGER, INTENT( OUT ) :: NNZJ
       LOGICAL, INTENT( IN ) :: GRAD
       INTEGER, INTENT( OUT ), DIMENSION( LCJAC ) :: INDVAR, INDFUN
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LC ) :: C
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LCJAC ) :: CJAC
       END SUBROUTINE CSCFG
       
       SUBROUTINE CCFSG( N, M, X, LC, C, NNZJ, LCJAC, CJAC,                    &
                         INDVAR, INDFUN, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LC, LCJAC
       INTEGER, INTENT( OUT ) :: NNZJ
       LOGICAL, INTENT( IN ) :: GRAD
       INTEGER, INTENT( OUT ), DIMENSION( LCJAC ) :: INDVAR, INDFUN
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LC ) :: C
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LCJAC ) :: CJAC
       END SUBROUTINE CCFSG
       
       SUBROUTINE CCIFG( N, ICON, X, CI, GCI, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, ICON
       LOGICAL, INTENT( IN ) :: GRAD
       REAL ( KIND = wp ), INTENT( OUT ) :: CI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: GCI
       END SUBROUTINE CCIFG
       
       SUBROUTINE CCIFSG( N, ICON, X, CI, NNZGCI, LGCI, GCI, INDVAR, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, ICON, LGCI
       INTEGER, INTENT( OUT ) :: NNZGCI
       LOGICAL, INTENT( IN ) :: GRAD
       INTEGER, INTENT( OUT ), DIMENSION( LGCI ) :: INDVAR
       REAL ( KIND = wp ), INTENT( OUT ) :: CI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LGCI ) :: GCI
       END SUBROUTINE CCIFSG
       
       SUBROUTINE CSCIFG( N, ICON, X, CI, NNZGCI, LGCI, GCI, INDVAR, GRAD )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, ICON, LGCI
       INTEGER, INTENT( OUT ) :: NNZGCI
       LOGICAL, INTENT( IN ) :: GRAD
       INTEGER, INTENT( OUT ), DIMENSION( LGCI ) :: INDVAR
       REAL ( KIND = wp ), INTENT( OUT ) :: CI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LGCI ) :: GCI
       END SUBROUTINE CSCIFG
       
       SUBROUTINE CDH( N, M, X, LV, V, LH1, H )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, M, LV, LH1
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH1, N ) :: H
       END SUBROUTINE CDH
       
       SUBROUTINE CIDH( N, X, IPROB, LH1, H )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, IPROB, LH1
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH1, N ) :: H
       END SUBROUTINE CIDH
       
       SUBROUTINE CGRDH( N , M , X     , GRLAGF, LV , V, G     ,               &
                         JTRANS, LCJAC1, LCJAC2, CJAC  , LH1, H     )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , M , LV    , LH1   , LCJAC1, LCJAC2
       LOGICAL, INTENT( IN ) :: GRLAGF, JTRANS
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: G
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH1, N ) :: H
       REAL ( KIND = wp ), INTENT( OUT ),                                      &
                                DIMENSION( LCJAC1, LCJAC2 ) :: CJAC
       END SUBROUTINE CGRDH
       
       SUBROUTINE CDIMSH( NNZH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: NNZH
       END SUBROUTINE CDIMSH
       
       SUBROUTINE CSH( N, M, X, LV, V, NNZH, LH, H, IRNH, ICNH  )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LV, LH
       INTEGER, INTENT( OUT ) :: NNZH
       INTEGER, INTENT( OUT ), DIMENSION( LH ) :: IRNH, ICNH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH ) :: H
       END SUBROUTINE CSH
       
       SUBROUTINE CISH( N, X, IPROB, NNZH, LH, H, IRNH, ICNH  )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, IPROB, LH
       INTEGER, INTENT( OUT ) :: NNZH
       INTEGER, INTENT( OUT ), DIMENSION( LH ) :: IRNH, ICNH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH ) :: H
       END SUBROUTINE CISH
       
       SUBROUTINE CDIMSE( NE, NNZH, NZIRNH )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( OUT ) :: NE, NNZH, NZIRNH
       END SUBROUTINE CDIMSE
       
       SUBROUTINE CEH( N , M , X , LV, V , NE, IRNHI , LIRNHI, LE    ,         &
                       IPRNHI, HI    , LHI   , IPRHI , BYROWS )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LV, LE, LIRNHI, LHI 
       INTEGER, INTENT( OUT ) :: NE
       LOGICAL, INTENT( IN ) :: BYROWS
       INTEGER, INTENT( OUT ), DIMENSION( LIRNHI ) :: IRNHI
       INTEGER, INTENT( OUT ), DIMENSION( LE ) :: IPRNHI, IPRHI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LHI ) :: HI
       END SUBROUTINE CEH
       
       SUBROUTINE CSGRSH( N , M , X     , GRLAGF, LV, V , NNZJ  ,              &
                          LCJAC , CJAC  , INDVAR, INDFUN, NNZH  ,              &
                          LH, H , IRNH  , ICNH  )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) ::  N, M, LV, LCJAC , LH
       INTEGER, INTENT( OUT ) :: NNZJ, NNZH
       LOGICAL, INTENT( IN ) ::  GRLAGF
       INTEGER, INTENT( OUT ), DIMENSION( LCJAC ) :: INDVAR, INDFUN
       INTEGER, INTENT( OUT ), DIMENSION( LH ) :: IRNH, ICNH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LH ) :: H
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LCJAC ) :: CJAC
       END SUBROUTINE CSGRSH
       
       SUBROUTINE CSGREH( N , M , X     , GRLAGF, LV, V , NNZJ  , LCJAC ,      &
                          CJAC  , INDVAR, INDFUN, NE    , IRNHI , LIRNHI,      &
                          LE    , IPRNHI, HI    , LHI   , IPRHI , BYROWS )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N, M, LV, LCJAC, LE, LIRNHI, LHI 
       INTEGER, INTENT( OUT ) :: NE, NNZJ
       LOGICAL, INTENT( IN ) :: GRLAGF, BYROWS
       INTEGER, INTENT( IN ), DIMENSION( LCJAC ) :: INDVAR, INDFUN
       INTEGER, INTENT( OUT ), DIMENSION( LIRNHI ) :: IRNHI
       INTEGER, INTENT( OUT ), DIMENSION( LE ) :: IPRNHI, IPRHI
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LHI ) :: HI
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( LCJAC ) :: CJAC
       END SUBROUTINE CSGREH
       
       SUBROUTINE CPROD( N , M , GOTH  , X , LV, V , P , RESULT )
       INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
       INTEGER, INTENT( IN ) :: N , M , LV
       LOGICAL, INTENT( IN ) :: GOTH
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( N ) :: X, P
       REAL ( KIND = wp ), INTENT( IN ), DIMENSION( LV ) :: V
       REAL ( KIND = wp ), INTENT( OUT ), DIMENSION( N ) :: RESULT
       END SUBROUTINE CPROD
  
     END INTERFACE
         
!  End of module CUTEr_interface_double

   END MODULE CUTEr_interface_double
      

