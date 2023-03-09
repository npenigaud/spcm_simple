INTERFACE
#if defined(_OPENACC)
SUBROUTINE SPCSI(&
& YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
& PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
& PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
& zsdiv,zhelp,zsp,zst,zsdivp,zspdivp,zsphi,zout,&
& PSPAUXG)
#else
SUBROUTINE SPCSI(&
& YDGEOMETRY,YDCST,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
& PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
& PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,&
& PSPAUXG)
#endif
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
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN), TARGET    :: YDDYN
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
#endif
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PSPAUXG(:,:)
END SUBROUTINE SPCSI

END INTERFACE
