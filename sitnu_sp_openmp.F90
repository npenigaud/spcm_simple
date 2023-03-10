#if defined(_OPENACC)
SUBROUTINE SITNU_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP,ZSDIV,ZOUT)
#else
SUBROUTINE SITNU_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)
#endif

!**** *SITNU_SP_OPENMP*   - Continuity equation for semi-implicit.

!     Purpose.
!     --------
!           Evaluate operators Tau and Nu in semi-implicit.

!**   Interface.
!     ----------
!        *CALL* *SITNU_SP_OPENMP(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE
!        PT    : TEMPERATURE
!        PSP   : SURFACE PRESSURE

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   None.
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
!                 value of LVERTFE in the SI NH linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP + optimization
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      R.Brozkova + NEC: Mar 2021: Optimization for vector (NEC)
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
 
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(KSPEC,klev)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KSPEC)
REAL(KIND=JPRB),intent(inout) :: ZSDIV(kspec,0:KLEV+1)          !!!!!!!on garde deux branches, mêmes ordres car il y a eu deux transpositions pour ces tableaux intermédiaires
REAL(KIND=JPRB),intent(inout) :: ZOUT(kspec,0:KLEV)
#else
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KSPEC)
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KSPEC,0:KLEV+1)
REAL(KIND=JPRB) :: ZOUT(KSPEC,0:KLEV)
#endif


REAL(KIND=JPRB) :: ZSDIVX(0:KLEV, KSPEC)
INTEGER(KIND=JPIM) :: JLEV, JSPEC
REAL(KIND=JPRB) :: ZREC
REAL(KIND=JPRB) :: ZDETAH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2

real(kind=JPRB)::intermediaire

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SITNU_SP_OPENMP',0,ZHOOK_HANDLE)

ASSOCIATE( YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDCVER=>YDGEOMETRY%YRCVER) 
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SIDELP=>YDDYN%SIDELP, SILNPR=>YDDYN%SILNPR, SIRDEL=>YDDYN%SIRDEL, &
 & SIRPRN=>YDDYN%SIRPRN, SITLAF=>YDDYN%SITLAF, SITR=>YDDYN%SITR, YDDIMV=>YDGEOMETRY%YRDIMV)

!     ------------------------------------------------------------------

!*       1.    SUM DIVERGENCE AND COMPUTES TEMPERATURE.
!              ----------------------------------------

IF(YDGEOMETRY%YRCVER%LVERTFE) THEN
IF (LHOOK) CALL DR_HOOK('SITNU_transfert1',0,ZHOOK_HANDLE2)
!!$acc data present(pd,psp,pt,sidelp,sitlaf,sitr,ydveta,ydcst,sirprn,klev,kspec)
!!$acc data present(pd,psp,pt,sidelp,sitlaf,sitr,ydveta,ydcst,sirprn,klev)
!$acc data present(pd,psp,pt,YDDYN%SIDELP,YDDYN%SITLAF,YDDYN%SITR,YDGEOMETRY%YRVETA,ydcst,YDDYN%SIRPRN)
!!$acc data create(zsdiv,zout,intermediaire) 
!$acc data present(zsdiv,zout)
!!copy(kspec)
IF (LHOOK) CALL DR_HOOK('SITNU_transfert1',1,ZHOOK_HANDLE2)

IF (LHOOK) CALL DR_HOOK('SITNU_transpose1',0,ZHOOK_HANDLE2)

#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH) default(none)
!$acc loop gang
  DO JLEV=1,KLEV
    ZDETAH=YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*YDDYN%SIDELP(JLEV)
    !$acc loop vector
    DO JSPEC=1,KSPEC               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pour les tableaux non intermédiaires changement, donc pas de transposition openacc
      ZSDIV(JSPEC,jlev)=PD(JSPEC,jlev)*ZDETAH
    ENDDO
  ENDDO
!$acc END PARALLEL
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZDETAH)
!$OMP DO SCHEDULE(STATIC)
  DO JLEV=1,KLEV
    ZDETAH=YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
    DO JSPEC=1,KSPEC
      ZSDIV(JSPEC,JLEV)=PD(JLEV,JSPEC)*YDDYN%SIDELP(JLEV)*ZDETAH
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_transpose1',1,ZHOOK_HANDLE2)

  IF (KSPEC>=1) THEN

IF (LHOOK) CALL DR_HOOK('SITNU_cond_lim',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
    !$acc parallel private(jspec) default(none)
    !$acc loop gang
    do jspec=1,kspec
      ZSDIV(jspec,0)=0.0_JPRB
      ZSDIV(jspec,KLEV+1)=0.0_JPRB
    enddo
    !$acc end parallel
#else
    ZSDIV(1:KSPEC,0)=0.0_JPRB
    ZSDIV(1:KSPEC,KLEV+1)=0.0_JPRB
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_cond_lim',1,ZHOOK_HANDLE2)

    CALL VERDISINT(YDGEOMETRY%YRVFE,YDGEOMETRY%YRCVER,'ITOP','11',KSPEC,1,KSPEC,KLEV,ZSDIV,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

intermediaire=YDCST%RKAPPA*YDDYN%SITR
IF (LHOOK) CALL DR_HOOK('SITNU_transpose2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JLEV,JSPEC,ZREC) default(none)
!$acc loop gang
  DO JLEV=1,KLEV
    ZREC=(1.0_JPRB/YDDYN%SITLAF(JLEV))*intermediaire
    !$acc loop vector
    DO JSPEC=1,KSPEC
      PT(JSPEC,jlev)=ZOUT(JSPEC,jlev-1)*ZREC
    ENDDO
  ENDDO
!$acc END PARALLEL
#else
!$OMP PARALLEL PRIVATE(JLEV,JSPEC,ZREC)
!$OMP DO SCHEDULE(STATIC)
  !DO JLEV=1,KLEV
  DO JSPEC=1,KSPEC
    !DO JSPEC=1,KSPEC
    DO JLEV=1,KLEV
      ZREC=1.0_JPRB/YDDYN%SITLAF(JLEV)
      !PT(JLEV,JSPEC)=YDCST%RKAPPA*SITR*ZOUT(JSPEC,JLEV-1)*ZREC
      PT(JLEV,JSPEC)=intermediaire*ZOUT(JSPEC,JLEV-1)*ZREC
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_transpose2',1,ZHOOK_HANDLE2)


IF (LHOOK) CALL DR_HOOK('SITNU_calcul1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$acc parallel loop private(jspec) default(none)
  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZOUT(JSPEC,klev)*YDDYN%SIRPRN
  ENDDO
!$acc end parallel
#else
  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZOUT(JSPEC,KLEV)*YDDYN%SIRPRN
  ENDDO
#endif
IF (LHOOK) CALL DR_HOOK('SITNU_calcul1',1,ZHOOK_HANDLE2)

IF (LHOOK) CALL DR_HOOK('SITNU_transfert2',0,ZHOOK_HANDLE2)
!$acc end data
!$acc end data
IF (LHOOK) CALL DR_HOOK('SITNU_transfert2',1,ZHOOK_HANDLE2)
ELSE

  ZSDIVX(0, :)=0.0_JPRB

IF (LHOOK) CALL DR_HOOK('SITNU_calcul2',0,ZHOOK_HANDLE2)
!$OMP PARALLEL PRIVATE(JLEV,JSPEC)
!$OMP DO SCHEDULE(STATIC)
  DO JSPEC=1,KSPEC
    DO JLEV=1,KLEV
      ZSDIVX(JLEV, JSPEC)=ZSDIVX(JLEV-1, JSPEC)+PD(JLEV,JSPEC)*YDDYN%SIDELP(JLEV)
      PT(JLEV,JSPEC)=YDCST%RKAPPA*YDDYN%SITR*(YDDYN%SIRDEL(JLEV)*YDDYN%SILNPR(JLEV)*ZSDIVX(JLEV-1, JSPEC)&
       & +YDDYN%SIALPH(JLEV)*PD(JLEV,JSPEC))
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
IF (LHOOK) CALL DR_HOOK('SITNU_calcul2',0,ZHOOK_HANDLE2)

  DO JSPEC=1,KSPEC
    PSP(JSPEC)=ZSDIVX(KLEV, JSPEC)*YDDYN%SIRPRN
  ENDDO

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SITNU_SP_OPENMP',1,ZHOOK_HANDLE)

END SUBROUTINE SITNU_SP_OPENMP

