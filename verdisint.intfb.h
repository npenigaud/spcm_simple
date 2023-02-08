INTERFACE
SUBROUTINE VERDISINT(YDVFE,YDCVER,CDOPER,CDBC,KPROMA,KSTART,KPROF,KFLEV,PIN,POUT,KCHUNK)
USE PARKIND1,ONLY : JPIM, JPRB, JPRD
USE YOMCVER ,ONLY : TCVER
USE YOMHOOK ,ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  ,ONLY : NULERR
USE YOMVERT ,ONLY : TVFE
TYPE(TVFE),TARGET ,INTENT(IN) :: YDVFE
TYPE(TCVER)       ,INTENT(IN) :: YDCVER
CHARACTER(LEN=*)  ,INTENT(IN) :: CDOPER, CDBC
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA, KSTART, KPROF, KFLEV
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KCHUNK
REAL(KIND=JPRB)   ,INTENT(IN) :: PIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT):: POUT(KPROMA,KFLEV+1)
END SUBROUTINE VERDISINT

END INTERFACE
