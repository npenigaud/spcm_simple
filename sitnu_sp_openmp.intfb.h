INTERFACE
SUBROUTINE SITNU_SP_OPENMP (YDGEOMETRY, YDCST, YDDYN, KLEV, KSPEC, PD, PT, PSP)
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPEC
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(KLEV,KSPEC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KSPEC)
END SUBROUTINE SITNU_SP_OPENMP

END INTERFACE
