
SUBROUTINE SIGAM_GP (YDGEOMETRY, YDCST, YDDYN, KLON, KLEV, PD, PT, PSP)

!**** *SIGAM_GP* - Solve hydrostatic operator in semi-implicit

!     Purpose.
!     --------
!           Operator gamma to compute p.

!**   Interface.
!     ----------
!        *CALL* *SIGAM_GP(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE       (output)
!        PT    : TEMPERATURE      (input)
!        PSP   : SURFACE PRESSURE (input)
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED
!        KFLEVG: NUMBER OF ELEMENTS IN A VERTICAL COLUMN

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!                 value of LVERTFE in the SI linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      R.Brozkova + NEC 03-Mar-2021: Optimization for vector (NEC)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KLON)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPHI(KLON,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KLON,0:KLEV)


REAL(KIND=JPRB) :: ZSPHIX(KLON, 0:KLEV)
INTEGER(KIND=JPIM) :: JLEV, JLON
CHARACTER(LEN=4):: CLOPER
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM_GP',0,ZHOOK_HANDLE)

ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER)
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------

CLOPER='IBOT'

IF (YDCVER%LVERTFE) THEN

  IF (YDCVER%LVFE_COMPATIBLE) CLOPER='INTG'

  DO JLEV=1,KLEV
    ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    DO JLON=1,KLON
      ZSPHI(JLON,JLEV)=-YDCST%RD*PT(JLON,JLEV)*SILNPR(JLEV)*ZDETAH
    ENDDO
  ENDDO

  ZSPHI(:,0)=0.0_JPRB
  ZSPHI(:,KLEV+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,YDCVER,CLOPER,'11',KLON,1,KLON,KLEV,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)

  DO JLEV=1,KLEV
    DO JLON=1,KLON
      PD(JLON,JLEV)=ZOUT(JLON,JLEV-1)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO

ELSE

  ZSPHIX(:, KLEV)=0.0_JPRB

  DO JLEV=KLEV,1,-1
    DO JLON=1,KLON,1
      ZSPHIX(JLON, JLEV-1)=ZSPHIX(JLON, JLEV)+YDCST%RD*PT(JLON,JLEV)*SILNPR(JLEV)
      PD(JLON,JLEV)=ZSPHIX(JLON, JLEV)+SIALPH(JLEV)*YDCST%RD*PT(JLON,JLEV)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO

ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SIGAM_GP',1,ZHOOK_HANDLE)

END SUBROUTINE SIGAM_GP
