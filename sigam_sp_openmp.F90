#if defined(_OPENACC)
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP,ZSPHI,ZOUT)
#else
SUBROUTINE SIGAM_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)
#endif

!**** *SIGAM_SP_OPENMP* - Solve hydrostatic operator in semi-implicit

!     Purpose.
!     --------
!           Operator gamma to compute p.

!**   Interface.
!     ----------
!        *CALL* *SIGAM_SP_OPENMP(...)

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KSPEC)
#if defined(_OPENACC)
REAL(KIND=JPRB),intent(inout) :: ZSPHI(0:KLEV+1,kspec)
REAL(KIND=JPRB),intent(inout) :: ZOUT(0:KLEV,kspec)
#else

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPHI(KSPEC,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KSPEC,0:KLEV)
#endif


REAL(KIND=JPRB) :: ZSPHIX(0:KLEV, KSPEC)
INTEGER(KIND=JPIM) :: JLEV, JSPEC
CHARACTER(LEN=4):: CLOPER
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2

real(kind=JPRB)::intermediaire

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM_SP_OPENMP',0,ZHOOK_HANDLE)

ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER)
!!!!!!ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------

CLOPER='IBOT'

IF (YDCVER%LVERTFE) THEN

  IF (YDCVER%LVFE_COMPATIBLE) CLOPER='INTG'
!!$acc data present(pt,psp,pd,klev,kspec)
!$acc data present(pt,psp,pd,klev)
!!$acc data present(pt,psp,pd)
!$acc data present(ydveta%vfe_rdetah,YDDYN%SILNPR,YDDYN%SIALPH,YDDYN%SIRPRG,ydveta,ydcst)
!!$acc data create(zsphi,zout,intermediaire) 
!!!!!!!$acc data create(zsphi,zout)
!$acc data present(zsphi,zout)
!!copy(kspec,klev)

IF (LHOOK) CALL DR_HOOK('SIGAM_transpose1',0,ZHOOK_HANDLE2)

#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH) default(none)
!$acc loop collapse(2) !gang
  DO JSPEC=1,KSPEC
  !!$acc loop vector
  DO JLEV=1,KLEV
    !v1 ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    !v2 ZDETAH=-YDVETA%VFE_RDETAH(JLEV)*YDCST%RD
    ZDETAH=-YDVETA%VFE_RDETAH(JLEV)*YDCST%RD*YDDYN%SILNPR(JLEV)
      !v1 ZSPHI(JSPEC,JLEV)=-YDCST%RD*PT(JLEV,JSPEC)*YDDYN%SILNPR(JLEV)*ZDETAH
      !v2 ZSPHI(JSPEC,JLEV)=PT(JLEV,JSPEC)*YDDYN%SILNPR(JLEV)*ZDETAH
      ZSPHI(JLEV,JSPEC)=PT(JLEV,JSPEC)*ZDETAH
    ENDDO
  ENDDO
!$acc end parallel
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH)
!$OMP DO SCHEDULE(STATIC)
  DO JLEV=1,KLEV
    !v1 ZDETAH=YDVETA%VFE_RDETAH(JLEV)
    !v2 ZDETAH=-YDVETA%VFE_RDETAH(JLEV)*YDCST%RD
    ZDETAH=-YDVETA%VFE_RDETAH(JLEV)*YDCST%RD*YDDYN%SILNPR(JLEV)
    !$acc loop vector
    DO JSPEC=1,KSPEC
      !v1 ZSPHI(JSPEC,JLEV)=-YDCST%RD*PT(JLEV,JSPEC)*YDDYN%SILNPR(JLEV)*ZDETAH
      !v2 ZSPHI(JSPEC,JLEV)=PT(JLEV,JSPEC)*YDDYN%SILNPR(JLEV)*ZDETAH
      ZSPHI(JSPEC,JLEV)=PT(JLEV,JSPEC)*ZDETAH
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif

IF (LHOOK) CALL DR_HOOK('SIGAM_transpose1',1,ZHOOK_HANDLE2)

IF (LHOOK) CALL DR_HOOK('SIGAM_cond_lim',0,ZHOOK_HANDLE2)
!!!!a verifier
#if defined(_OPENACC)
!$acc parallel private(jspec) default(none)
!$acc loop gang
do jspec=1,kspec
  ZSPHI(0,jspec)=0.0_JPRB
  ZSPHI(KLEV+1,jspec)=0.0_JPRB
enddo
!$acc end parallel
#else
  ZSPHI(:,0)=0.0_JPRB
  ZSPHI(:,KLEV+1)=0.0_JPRB
#endif
IF (LHOOK) CALL DR_HOOK('SIGAM_cond_lim',1,ZHOOK_HANDLE2)

  CALL VERDISINT(YDVFE,YDCVER,CLOPER,'11',KSPEC,1,KSPEC,KLEV,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)

IF (LHOOK) CALL DR_HOOK('SIGAM_transpose2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JLEV,JSPEC,intermediaire) default(none)
!$acc loop gang
  DO JSPEC=1,KSPEC
     intermediaire=PSP(JSPEC)*YDDYN%SIRPRG
    !$acc loop vector
    !DO JSPEC=1,KSPEC
     DO JLEV=1,KLEV
      PD(JLEV,JSPEC)=ZOUT(JLEV-1,JSPEC)+intermediaire
      !PD(JLEV,JSPEC)=ZOUT(JSPEC,JLEV-1)+PSP(JSPEC)*SIRPRG
    ENDDO
  ENDDO
!$acc END PARALLEL
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,intermediaire)
!$OMP DO SCHEDULE(STATIC)
  !DO JLEV=1,KLEV
  DO JSPEC=1,KSPEC
     intermediaire=PSP(JSPEC)*YDDYN%SIRPRG
    !DO JSPEC=1,KSPEC
     DO JLEV=1,KLEV
      PD(JLEV,JSPEC)=ZOUT(JSPEC,JLEV-1)+intermediaire
      !PD(JLEV,JSPEC)=ZOUT(JSPEC,JLEV-1)+PSP(JSPEC)*SIRPRG
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
IF (LHOOK) CALL DR_HOOK('SIGAM_transpose2',1,ZHOOK_HANDLE2)

!$acc end data
!$acc end data
!$acc end data
ELSE
  ZSPHIX(KLEV, :)=0.0_JPRB

!$OMP PARALLEL PRIVATE(JLEV,JSPEC)
!$OMP DO SCHEDULE(STATIC)
  DO JSPEC=1,KSPEC,1
    DO JLEV=KLEV,1,-1
      ZSPHIX(JLEV-1, JSPEC)=ZSPHIX(JLEV, JSPEC)+YDCST%RD*PT(JLEV,JSPEC)*YDDYN%SILNPR(JLEV)
      PD(JLEV,JSPEC)=ZSPHIX(JLEV, JSPEC)+YDDYN%SIALPH(JLEV)*YDCST%RD*PT(JLEV,JSPEC)+PSP(JSPEC)*YDDYN%SIRPRG
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
!!!!!!END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SIGAM_SP_OPENMP',1,ZHOOK_HANDLE)

END SUBROUTINE SIGAM_SP_OPENMP
