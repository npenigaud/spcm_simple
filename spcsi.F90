#if defined(_OPENACC)
SUBROUTINE SPCSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 & zsdiv,zhelp,zsp,zst,zsdivp,zspdivp,zsphi,zout,&
 ! --- INPUT OPTIONAL --------------------------------------------------------
 & PSPAUXG)
#else
SUBROUTINE SPCSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
 ! --- INPUT OPTIONAL --------------------------------------------------------
 & PSPAUXG)
#endif

!**** *SPCSI* - SPECTRAL SPACE SEMI-IMPLICIT COMPUTATIONS FOR HYD MODEL.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCSI(..)

!        Explicit arguments :
!        --------------------  KM      - Zonal wavenumber 
!                              KMLOC   - Zonal wavenumber (DM-local numbering)
!                              KSTA    - First column processed
!                              KEND    - Last column processed
!                              LDONEM  - T if only one m if processed
!                              PSPVORG - Vorticity columns
!                              PSPDIVG - Divergence columns
!                              PSPTG   - Temperature columns
!                              PSPSPG  - Surface Pressure
!                              PSPTNDSI_VORG - [D vor/Dt]_SI
!                              PSPTNDSI_DIVG - [D div/Dt]_SI
!                              PSPTNDSI_TG   - [D T/Dt]_SI

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-11-24 (before 1997 spcsi.F was part of spc.F)

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N.Wedi        08-Mar-2005 remove mass correction      
!      K.Yessad 09-Dec-2004: move mass correction in SPCMASCOR + cleanings.
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad (Aug 2009): remove LSITRIC option.
!      F. Voitus: add DDH diagnostics.
!      T.Wilhelmsson 09-09-25: Remove LFULLM requirement for LIMPF
!      K. Yessad (Feb 2012): tests on LL3D, LLDOSI in the caller, simplifications.
!      P. Marguinaud (Nov 2012): Fix unallocated array arguments
!      P. Marguinaud (Sep 2012) : Make PSPAUXG optional
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      O. Marsden (May 2016): Removed redundant geometry arguments
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : MYSETV
USE YOMDYN       , ONLY : TDYN
USE YOMLDDH      , ONLY : TLDDH
USE YOMRIP       , ONLY : TRIP
USE YOMCST       , ONLY : TCST

#if defined(_OPENACC)
use cublas
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(:,:)
 
#if defined(_OPENACC)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zsdiv(1:YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zhelp(1:YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zst(1:YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zsp(1:YDGEOMETRY%YRMP%NSPEC2V)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zsdivp(max(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF),YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zspdivp(max(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF),YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zsphi(YDGEOMETRY%YRMP%NSPEC2V,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) ::zout(YDGEOMETRY%YRMP%NSPEC2V,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PSPAUXG(:,:)
#else
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PSPAUXG(:,:)
!     ------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: ZSDIVP (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVP(:,:)

REAL(KIND=JPRB) :: ZSDIV  (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZHELP  (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZST    (YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND) 
REAL(KIND=JPRB) :: ZSP    (       KSTA:KEND)                   
#endif

REAL(KIND=JPRB) :: ZSDIVPL (YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)

!!!!#if defined(_OPENACC)
!!!!real(kind=jprb)::zsphi(kend-ksta+1,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
!!!!real(kind=jprb)::zout(kend-ksta+1,0:YDGEOMETRY%YRDIMV%NFLEVG)
!!!!#endif

INTEGER(KIND=JPIM) :: II, IN, IOFF, IS0, IS02, ISE, ISPCOL, JLEV, JN, JSP  

REAL(KIND=JPRB) :: ZBDT, ZBDT2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE,ZHOOK_HANDLE2

real(kind=jprb)::intermediaire
!!integer(kind=jpim)::compteurprint1,compteurprint2

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"
#include "abor1.intfb.h"
#include "sigam_sp_openmp.intfb.h"
#include "spcimpfsolve.intfb.h"
#include "sitnu_sp_openmp.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPCSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, RBTS2=>YDDYN%RBTS2, &
 & SIHEG=>YDDYN%SIHEG, SIHEG2=>YDDYN%SIHEG2, SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, &
 & SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & LRSIDDH=>YDLDDH%LRSIDDH, &
 & NPTRSV=>YDMP%NPTRSV, NPTRSVF=>YDMP%NPTRSVF, NSPEC2V=>YDMP%NSPEC2V, &
 & NSPEC2VF=>YDMP%NSPEC2VF, &
 & TDT=>YDRIP%TDT, &
 & SCGMAP=>YDSPGEOM%SCGMAP)
!     ------------------------------------------------------------------

!*       0.    TESTINGS.
!              ---------


IF (LIMPF .AND. .NOT.PRESENT(PSPAUXG)) THEN
  CALL ABOR1(' SPCSI: If LIMPF=T, argument PSPAUXG must be present!')
ENDIF

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------
!!faux cas test
IF (LRSIDDH) THEN
  ! DDH memory transfer
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=-PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=-PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG  (1:NFLEVG,KSTA:KEND)=-PSPTG  (1:NFLEVG,KSTA:KEND)
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

#if defined(_OPENACC)

#else
ALLOCATE(ZSDIVP(NFLEVG,MAX(NSPEC2V,NSPEC2VF)))
ALLOCATE(ZSPDIVP(NFLEVG,MAX(NSPEC2V,NSPEC2VF)))
#endif


!*        2.1  Preliminary initialisations.

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF
ISPCOL=KEND-KSTA+1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

!*        2.3  Computes right-hand side of Helmholtz equation.
#if defined(_OPENACC)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',0,ZHOOK_HANDLE2)
!!!$acc data present(YDLAP,nflevg,nsmax,sivp,rstret) 
!!$acc data copy(pspdivg,psptg,pspspg,pspauxg,zsdivpl,zspdivpl)
!!$acc data present(pspdivg,psptg,pspspg) copy(ispcol,zbdt,zbdt2)
!!copyin(ispcol,zbdt,zbdt2) 
!!create(zhelp,zsp,zst)
!!$acc data copyout(zsdiv,zsdivp,zspdivp,zhelp,zsp,zst)

!$acc data present(YDGEOMETRY,YDGEOMETRY%YRLAP,YDGEOMETRY%YRLAP%NVALUE,YDGEOMETRY%YRLAP%RLAPIN,YDGEOMETRY%YRLAP%RLAPDI,nflevg,nsmax,YDDYN,YDDYN%SIVP,rstret) !!copy(ispcol,zbdt,zbdt2)
!!$acc data present(YDGEOMETRY%YRLAP,nflevg,nsmax,YDDYN%SIVP,rstret) !!copy(ispcol,zbdt,zbdt2)
!!$acc data present(YDLAP,nflevg,nsmax,sivp,rstret)
!!$acc data copy(pspdivg,psptg,pspspg,pspauxg,zsdivpl,zspdivpl)
!$acc data present(pspdivg,psptg,pspspg)
!!$acc data copyout(zsdiv,zhelp,zsp,zst,zsdivp,zspdivp)
!!!!!!!$acc data create(zsdiv,zhelp,zsp,zst,zsdivp,zspdivp)
!$acc data present(zsdiv,zhelp,zsp,zst,zsdivp,zspdivp,zsphi,zout)

!!faux cas test
if (limpf) then 
!$acc enter data copyin(pspauxg)
end if
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts1',1,ZHOOK_HANDLE2)
#endif

!!IF( .NOT.LDONEM ) CALL GSTATS(1655,0) ! Main routines and loops in SIGAM chain are parallel
#if defined(_OPENACC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!point à revoir sur le passage de sous-tableaux
!!! CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,ISPCOL,ZSDIV,PSPTG(:,KSTA:KEND),PSPSPG(KSTA:KEND),ZSPHI,ZOUT)  passage taille
print *,"valeur de ksta : ",ksta
print *,"valeur de kend : ",kend
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,ISPCOL,ZSDIV,PSPTG(:,:),PSPSPG(KSTA:KEND),ZSPHI,ZOUT)
#else
CALL SIGAM_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,ISPCOL,ZSDIV,PSPTG(:,KSTA:KEND),PSPSPG(KSTA:KEND))
#endif

!!IF( .NOT.LDONEM ) CALL GSTATS(1655,1)

!!IF( .NOT.LDONEM ) CALL GSTATS(1656,0)

!!faux cas test
IF (LSIDG) THEN
  IF (KM > 0) THEN
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV,IN) default(none)
!$acc loop gang
    DO JLEV=1,NFLEVG
      !$acc loop vector
      DO JSP=KSTA,KEND
        IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
        ZSDIV(JSP,jlev)=YDGEOMETRY%YRLAP%RLAPIN(IN)*PSPDIVG(JSP,jlev)-ZBDT*ZSDIV(JSP,jlev)      
      ENDDO
    ENDDO
!$acc end parallel
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
    !$acc loop vector
      DO JLEV=1,NFLEVG
        IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=YDGEOMETRY%YRLAP%RLAPIN(IN)*PSPDIVG(JLEV,JSP)-ZBDT*ZSDIV(JLEV,JSP)      
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
#endif
  ELSE
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV,IN) default(none)
   !$acc loop gang
    DO JSP=KSTA,KEND
    !$acc loop vector
      DO JLEV=1,NFLEVG
        IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)  
      ENDDO
    ENDDO
!$acc end parallel
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
    !$acc loop vector
      DO JLEV=1,NFLEVG
        IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)  
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
#endif
  ENDIF
ELSE
!!parcouru cas test
  ! Case of No Stretching
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle1',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV,IN,intermediaire) default(none)
!$acc loop gang
  DO JLEV=1,NFLEVG
    !$acc loop vector
    DO JSP=KSTA,KEND
      IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
      intermediaire=ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)
      ZSDIV(JSP,jlev)=PSPDIVG(JSP,jlev)-intermediaire*ZSDIV(JSP,jlev)
    ENDDO
  ENDDO
!$acc end parallel
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV,IN)
  DO JSP=KSTA,KEND
  IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
  intermediaire=ZBDT*YDGEOMETRY%YRLAP%RLAPDI(IN)
    DO JLEV=1,NFLEVG

      ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-intermediaire*ZSDIV(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle1',1,ZHOOK_HANDLE2)
ENDIF

!        Add [F] * result to rhs of Helmholtz equation

!!faux cas test
IF (LIMPF) THEN
!!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV)
!$acc loop gang
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
#endif
  DO JSP=KSTA,KEND
  !$acc loop vector
    DO JLEV=1,NFLEVG
      ZSDIV(JLEV,JSP)=ZSDIV(JLEV,JSP) + PSPAUXG(JLEV,JSP)
    ENDDO
  ENDDO
#if defined(_OPENACC)
!$acc END PARALLEL 
#else
!$OMP END PARALLEL DO
#endif
ENDIF
!!IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

!*        2.4  Solve Helmholtz equation

!           Current space --> vertical eigenmodes space.

!!IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Call to SGEMMX Parallelised
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',0,ZHOOK_HANDLE2)
!!$acc update host(zsdiv,zsdivp)
!!CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP(:,KSTA:KEND),1,NFLEVG,&
 !!& NFLEVG,NFLEVG,ISPCOL) 
!!!$acc update device(zsdiv,zsdivp)


!!do compteurprint1=1,nflevg
!!  do compteurprint2=1,nflevg
!!     print *,simi(compteurprint1,compteurprint2)
!!  end do
!!  print *,"ligne"
!!end do


#if defined(_OPENACC)
!$acc host_data use_device(SIMI,ZSDIV,ZSDIVP)
!!!!!!CALL cublasDgemm ('N','N',NFLEVG,ISPCOL,NFLEVG,1.0_JPRB, &
!!!!!!  & SIMI,NFLEVG,ZSDIV,NFLEVG,0.0_JPRB,ZSDIVP(1,ksta),NFLEVG)
CALL cublasDgemm ('N','T',ispcol,nflevg,nflevg,1.0_JPRB, &
  & zsdiv,ispcol,simi,NFLEVG,0.0_JPRB,ZSDIVP(ksta,1),ispcol)
!$acc end host_data
!$acc wait
#else
CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP(:,KSTA:KEND),1,NFLEVG,&
 & NFLEVG,NFLEVG,ISPCOL)
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop1',1,ZHOOK_HANDLE2)

!!IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

!!faux cas test
IF (LSIDG) THEN

  !             Inversion of two tridiagonal systems (Helmholtz equation)
  !                --> (SIMI*DIVprim(t+dt)).

  !             Reorganisation of divergence

  IS0=YDGEOMETRY%YRLAP%NSE0L(KMLOC)
  IS02=0
  II=MIN(KM,1)+1
  
  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB
  
  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSDIVPL(:,JN,1:2)=ZSDIVP(:,ISE:ISE+1)
  ENDDO
  
  IF (KM > 0) THEN

    !               Inversion of a symmetric matrix.

    CALL MXTURS(NSMAX+1-KM,NFLEVG,NFLEVG,II,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)  
  ELSE

    !               Inversion of a non-symmetric matrix.

    CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,-2,.TRUE.,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)  
    CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,3,.FALSE.,&
     & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
     & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)  
  ENDIF

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSPDIVP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
  ENDDO

ELSE
!!parcouru cas test, dans le else

  !             Case with NO Stretching :
  !!!non traite pour le moment, on n'y passe pas
  IF (LIMPF) THEN

    !               Solve complex pentadiagonal system

    CALL SPCIMPFSOLVE(YDGEOMETRY,YDCST,YDRIP,YDDYN,.FALSE.,.FALSE.,LDONEM,ZSDIVP,ZSPDIVP)
  !!parcouru cas test
  ELSE
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle2',0,ZHOOK_HANDLE2)
    !                 Inversion of a diagonal system (Helmholtz equation)
    !                 --> (SIMI*DIVprim(t+dt)).

#if defined(_OPENACC)
    !$acc parallel private(JSP,JLEV,intermediaire) default(none)
    !$acc loop gang
    DO JLEV=1,NFLEVG	
      intermediaire=ZBDT2*YDDYN%SIVP(JLEV)
      !$acc loop vector
      DO JSP=KSTA,KEND
        ZSPDIVP(JSP,jlev)=ZSDIVP(JSP,jlev)/(1.0_JPRB-intermediaire*YDGEOMETRY%YRLAP%RLAPDI(YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
    !$acc end parallel
#else
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
         & /(1.0_JPRB-ZBDT2*YDDYN%SIVP(JLEV)*YDGEOMETRY%YRLAP%RLAPDI(YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)))  
      ENDDO
    ENDDO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle2',1,ZHOOK_HANDLE2)
  ENDIF
ENDIF

!           Vertical eigenmodes space --> current space.

!!IF( .NOT.LDONEM ) CALL GSTATS(1660,0) ! MXMAOP Calls SGEMMX in parallel region

IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop2',0,ZHOOK_HANDLE2)
#if defined(_OPENACC)
!$acc host_data use_device(SIMO,ZSPDIVP,PSPDIVG)
!!!!!!!CALL cublasDgemm ('N','N',NFLEVG,ISPCOL,NFLEVG,1.0_JPRB, &
!!!!!!!  & SIMO,NFLEVG,ZSPDIVP(1,KSTA),NFLEVG,0.0_JPRB,PSPDIVG(1,KSTA),NFLEVG)
CALL cublasDgemm ('N','T',ispcol,nflevg,NFLEVG,1.0_JPRB, &
  & ZSPDIVP(KSTA,1),ispcol,SIMO,NFLEVG,0.0_JPRB,PSPDIVG(KSTA,1),ispcol)
!$acc end host_data
!$acc wait
#else
CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP(:,KSTA:KEND),1,NFLEVG,PSPDIVG(:,KSTA:KEND),1,&
 & NFLEVG,NFLEVG,NFLEVG,ISPCOL)
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_mxmaop2',1,ZHOOK_HANDLE2)

!!IF( .NOT.LDONEM ) CALL GSTATS(1660,1)

!!non parcourue dans le cas test
IF (LSIDG) THEN

  !           ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GM**2 * DIVprim(t+dt)) .

  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB

  !           Reorganisation of ZSDIVP (Back to the USSR)

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZSDIVPL(:,JN,1:2)=PSPDIVG(:,ISE:ISE+1)
  ENDDO

  !        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).

  CALL MXPTMA(NSMAX+1-KM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
   & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
   & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
   & ZSDIVPL,ZSPDIVPL)  

  !           Reorganisation of ZSPDIVPL

  DO JN=KM,NSMAX
    ISE=KSTA+2*(JN-KM)
    ZHELP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
  ENDDO
!!!parcourue dans le cas test
ELSE
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle3',0,ZHOOK_HANDLE2)
  !       ZSPDIV=(DIVprim(t+dt)) --> ZSPDIVG=(GMBAR**2 * DIVprim(t+dt)) .

!!  IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
#if defined(_OPENACC)
intermediaire=RSTRET*RSTRET
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop gang
  DO JLEV=1,NFLEVG
    !$acc loop vector
    DO JSP=KSTA,KEND
      ZHELP(JSP,jlev)=PSPDIVG(JSP,jlev)*intermediaire
    ENDDO
  ENDDO
!$acc end parallel
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZHELP(JLEV,JSP)=PSPDIVG(JLEV,JSP)*RSTRET*RSTRET
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle3',1,ZHOOK_HANDLE2)
!!  IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

ENDIF

!       If LSIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or if not LSIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]

!!IF( .NOT.LDONEM ) CALL GSTATS(1657,0)  ! Main routines and loops in SITNU chain are parallel
#if defined(_OPENACC)
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,ISPCOL,ZHELP,ZST,ZSP,ZSPHI,ZOUT)
#else
CALL SITNU_SP_OPENMP(YDGEOMETRY,YDCST,YDDYN,NFLEVG,ISPCOL,ZHELP,ZST,ZSP)
#endif
!!IF( .NOT.LDONEM ) CALL GSTATS(1657,1)

!*       2.5  Increment Temperature and surface pressure

!!IF( .NOT.LDONEM ) CALL GSTATS(1656,0)
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle4',0,ZHOOK_HANDLE2)
!!!!!!!!!!!!!!!!!!!!possible de mettre les deux en parallele
#if defined(_OPENACC)
!$acc PARALLEL PRIVATE(JSP,JLEV) default(none)
!$acc loop gang
DO JLEV=1,NFLEVG
!$acc loop vector
  DO JSP=KSTA,KEND
    PSPTG(JSP,jlev)=PSPTG(JSP,jlev)-ZBDT*ZST(JSP,jlev)
  ENDDO
  
ENDDO
!$acc end parallel

!$acc PARALLEL PRIVATE(JSP) default(none)
DO JSP=KSTA,KEND
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
enddo
!$acc end parallel
#else
!$OMP PARALLEL DO PRIVATE(JSP,JLEV)
DO JSP=KSTA,KEND
!$acc loop vector
  DO JLEV=1,NFLEVG
    PSPTG(JLEV,JSP)=PSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
  ENDDO
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO
!$OMP END PARALLEL DO
#endif
IF (LHOOK) CALL DR_HOOK('SPCSI_boucle4',1,ZHOOK_HANDLE2)
!!IF( .NOT.LDONEM ) CALL GSTATS(1656,1)

#if defined(_OPENACC)
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts2',0,ZHOOK_HANDLE2)
if (limpf) then 
!$acc exit data copyout(pspauxg)
end if

!$acc end data !!!copyout(zsdiv,zhelp,zst,zsp)
!$acc end data !!!copy(pspdivg,psptg,pspspg,zspdivp,zsdivp)
!$acc end data !!present
IF (LHOOK) CALL DR_HOOK('SPCSI_transferts2',1,ZHOOK_HANDLE2)
#endif

#if defined(_OPENACC)
#else
DEALLOCATE(ZSDIVP)
DEALLOCATE(ZSPDIVP)
#endif

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

IF (LRSIDDH) THEN
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=&
   & PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND) + PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)&
   & + PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)&
   & + PSPTG(1:NFLEVG,KSTA:KEND)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCSI',1,ZHOOK_HANDLE)
END SUBROUTINE SPCSI
