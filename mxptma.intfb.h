INTERFACE
SUBROUTINE MXPTMA(KLX,KVX,KVXS,KIX,PA,PBI,PCI,PBS,PCS,PX,PY)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
INTEGER(KIND=JPIM),INTENT(IN)    :: KLX
INTEGER(KIND=JPIM),INTENT(IN)    :: KVXS
INTEGER(KIND=JPIM),INTENT(IN)    :: KIX
INTEGER(KIND=JPIM),INTENT(IN)    :: KVX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBI(KLX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCI(KLX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBS(KLX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCS(KLX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PX(KVXS,KLX,KIX)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PY(KVXS,KLX,KIX)
END SUBROUTINE MXPTMA

END INTERFACE
