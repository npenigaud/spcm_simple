MODULE UTIL_TLEP_MOD

USE YEMLAP, ONLY : TLEP

INTERFACE SAVE
MODULE PROCEDURE SAVE_TLEP
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_TLEP
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_TLEP
END INTERFACE

INTERFACE WIPE
MODULE PROCEDURE WIPE_TLEP
END INTERFACE



CONTAINS

SUBROUTINE SAVE_TLEP (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TLEP), INTENT (IN), TARGET :: YD
LOGICAL :: LMVALUE, LNCPL2M, LNCPL2N, LNCPL4M, LNCPL4N, LNCPLM, LNCPLN, LNESM0, LNESM0G, LNESPZERO
LOGICAL :: LNPME, LNPNE, LRLEPDIM, LRLEPDIN, LRLEPINM, LRLEPINN
LNCPL2M = ALLOCATED (YD%NCPL2M)
WRITE (KLUN) LNCPL2M
IF (LNCPL2M) THEN
  WRITE (KLUN) LBOUND (YD%NCPL2M)
  WRITE (KLUN) UBOUND (YD%NCPL2M)
  WRITE (KLUN) YD%NCPL2M
ENDIF
LNCPL4M = ALLOCATED (YD%NCPL4M)
WRITE (KLUN) LNCPL4M
IF (LNCPL4M) THEN
  WRITE (KLUN) LBOUND (YD%NCPL4M)
  WRITE (KLUN) UBOUND (YD%NCPL4M)
  WRITE (KLUN) YD%NCPL4M
ENDIF
LNCPLM = ALLOCATED (YD%NCPLM)
WRITE (KLUN) LNCPLM
IF (LNCPLM) THEN
  WRITE (KLUN) LBOUND (YD%NCPLM)
  WRITE (KLUN) UBOUND (YD%NCPLM)
  WRITE (KLUN) YD%NCPLM
ENDIF
LNCPL2N = ALLOCATED (YD%NCPL2N)
WRITE (KLUN) LNCPL2N
IF (LNCPL2N) THEN
  WRITE (KLUN) LBOUND (YD%NCPL2N)
  WRITE (KLUN) UBOUND (YD%NCPL2N)
  WRITE (KLUN) YD%NCPL2N
ENDIF
LNCPL4N = ALLOCATED (YD%NCPL4N)
WRITE (KLUN) LNCPL4N
IF (LNCPL4N) THEN
  WRITE (KLUN) LBOUND (YD%NCPL4N)
  WRITE (KLUN) UBOUND (YD%NCPL4N)
  WRITE (KLUN) YD%NCPL4N
ENDIF
LNCPLN = ALLOCATED (YD%NCPLN)
WRITE (KLUN) LNCPLN
IF (LNCPLN) THEN
  WRITE (KLUN) LBOUND (YD%NCPLN)
  WRITE (KLUN) UBOUND (YD%NCPLN)
  WRITE (KLUN) YD%NCPLN
ENDIF
LRLEPDIN = ALLOCATED (YD%RLEPDIN)
WRITE (KLUN) LRLEPDIN
IF (LRLEPDIN) THEN
  WRITE (KLUN) LBOUND (YD%RLEPDIN)
  WRITE (KLUN) UBOUND (YD%RLEPDIN)
  WRITE (KLUN) YD%RLEPDIN
ENDIF
LRLEPINN = ALLOCATED (YD%RLEPINN)
WRITE (KLUN) LRLEPINN
IF (LRLEPINN) THEN
  WRITE (KLUN) LBOUND (YD%RLEPINN)
  WRITE (KLUN) UBOUND (YD%RLEPINN)
  WRITE (KLUN) YD%RLEPINN
ENDIF
LRLEPDIM = ALLOCATED (YD%RLEPDIM)
WRITE (KLUN) LRLEPDIM
IF (LRLEPDIM) THEN
  WRITE (KLUN) LBOUND (YD%RLEPDIM)
  WRITE (KLUN) UBOUND (YD%RLEPDIM)
  WRITE (KLUN) YD%RLEPDIM
ENDIF
LRLEPINM = ALLOCATED (YD%RLEPINM)
WRITE (KLUN) LRLEPINM
IF (LRLEPINM) THEN
  WRITE (KLUN) LBOUND (YD%RLEPINM)
  WRITE (KLUN) UBOUND (YD%RLEPINM)
  WRITE (KLUN) YD%RLEPINM
ENDIF
LNESM0 = ALLOCATED (YD%NESM0)
WRITE (KLUN) LNESM0
IF (LNESM0) THEN
  WRITE (KLUN) LBOUND (YD%NESM0)
  WRITE (KLUN) UBOUND (YD%NESM0)
  WRITE (KLUN) YD%NESM0
ENDIF
LNESPZERO = ALLOCATED (YD%NESPZERO)
WRITE (KLUN) LNESPZERO
IF (LNESPZERO) THEN
  WRITE (KLUN) LBOUND (YD%NESPZERO)
  WRITE (KLUN) UBOUND (YD%NESPZERO)
  WRITE (KLUN) YD%NESPZERO
ENDIF
LNESM0G = ALLOCATED (YD%NESM0G)
WRITE (KLUN) LNESM0G
IF (LNESM0G) THEN
  WRITE (KLUN) LBOUND (YD%NESM0G)
  WRITE (KLUN) UBOUND (YD%NESM0G)
  WRITE (KLUN) YD%NESM0G
ENDIF
LNPME = ALLOCATED (YD%NPME)
WRITE (KLUN) LNPME
IF (LNPME) THEN
  WRITE (KLUN) LBOUND (YD%NPME)
  WRITE (KLUN) UBOUND (YD%NPME)
  WRITE (KLUN) YD%NPME
ENDIF
LNPNE = ALLOCATED (YD%NPNE)
WRITE (KLUN) LNPNE
IF (LNPNE) THEN
  WRITE (KLUN) LBOUND (YD%NPNE)
  WRITE (KLUN) UBOUND (YD%NPNE)
  WRITE (KLUN) YD%NPNE
ENDIF
LMVALUE = ALLOCATED (YD%MVALUE)
WRITE (KLUN) LMVALUE
IF (LMVALUE) THEN
  WRITE (KLUN) LBOUND (YD%MVALUE)
  WRITE (KLUN) UBOUND (YD%MVALUE)
  WRITE (KLUN) YD%MVALUE
ENDIF
END SUBROUTINE

SUBROUTINE LOAD_TLEP (KLUN, YD)
USE PARKIND1, ONLY : JPRD

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TLEP), INTENT (OUT), TARGET :: YD
INTEGER :: IL1(1), IU1(1)
LOGICAL :: LMVALUE, LNCPL2M, LNCPL2N, LNCPL4M, LNCPL4N, LNCPLM, LNCPLN, LNESM0, LNESM0G, LNESPZERO
LOGICAL :: LNPME, LNPNE, LRLEPDIM, LRLEPDIN, LRLEPINM, LRLEPINN
REAL(KIND=JPRD), ALLOCATABLE :: ZTMP1 (:)
READ (KLUN) LNCPL2M
IF (LNCPL2M) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPL2M (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPL2M
ENDIF
READ (KLUN) LNCPL4M
IF (LNCPL4M) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPL4M (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPL4M
ENDIF
READ (KLUN) LNCPLM
IF (LNCPLM) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPLM (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPLM
ENDIF
READ (KLUN) LNCPL2N
IF (LNCPL2N) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPL2N (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPL2N
ENDIF
READ (KLUN) LNCPL4N
IF (LNCPL4N) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPL4N (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPL4N
ENDIF
READ (KLUN) LNCPLN
IF (LNCPLN) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NCPLN (IL1(1):IU1(1)))
  READ (KLUN) YD%NCPLN
ENDIF
READ (KLUN) LRLEPDIN
IF (LRLEPDIN) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%RLEPDIN (IL1(1):IU1(1)))
  ALLOCATE (ZTMP1(LBOUND(YD%RLEPDIN,1):UBOUND(YD%RLEPDIN,1)))
  READ (KLUN) ZTMP1
  YD%RLEPDIN = ZTMP1
  DEALLOCATE (ZTMP1)
ENDIF
READ (KLUN) LRLEPINN
IF (LRLEPINN) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%RLEPINN (IL1(1):IU1(1)))
  ALLOCATE (ZTMP1(LBOUND(YD%RLEPINN,1):UBOUND(YD%RLEPINN,1)))
  READ (KLUN) ZTMP1
  YD%RLEPINN = ZTMP1
  DEALLOCATE (ZTMP1)
ENDIF
READ (KLUN) LRLEPDIM
IF (LRLEPDIM) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%RLEPDIM (IL1(1):IU1(1)))
  ALLOCATE (ZTMP1(LBOUND(YD%RLEPDIM,1):UBOUND(YD%RLEPDIM,1)))
  READ (KLUN) ZTMP1
  YD%RLEPDIM = ZTMP1
  DEALLOCATE (ZTMP1)
ENDIF
READ (KLUN) LRLEPINM
IF (LRLEPINM) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%RLEPINM (IL1(1):IU1(1)))
  ALLOCATE (ZTMP1(LBOUND(YD%RLEPINM,1):UBOUND(YD%RLEPINM,1)))
  READ (KLUN) ZTMP1
  YD%RLEPINM = ZTMP1
  DEALLOCATE (ZTMP1)
ENDIF
READ (KLUN) LNESM0
IF (LNESM0) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NESM0 (IL1(1):IU1(1)))
  READ (KLUN) YD%NESM0
ENDIF
READ (KLUN) LNESPZERO
IF (LNESPZERO) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NESPZERO (IL1(1):IU1(1)))
  READ (KLUN) YD%NESPZERO
ENDIF
READ (KLUN) LNESM0G
IF (LNESM0G) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NESM0G (IL1(1):IU1(1)))
  READ (KLUN) YD%NESM0G
ENDIF
READ (KLUN) LNPME
IF (LNPME) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NPME (IL1(1):IU1(1)))
  READ (KLUN) YD%NPME
ENDIF
READ (KLUN) LNPNE
IF (LNPNE) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%NPNE (IL1(1):IU1(1)))
  READ (KLUN) YD%NPNE
ENDIF
READ (KLUN) LMVALUE
IF (LMVALUE) THEN
  READ (KLUN) IL1
  READ (KLUN) IU1
  ALLOCATE (YD%MVALUE (IL1(1):IU1(1)))
  READ (KLUN) YD%MVALUE
ENDIF
END SUBROUTINE


SUBROUTINE COPY_TLEP (YD, LDCREATED)

IMPLICIT NONE
TYPE (TLEP), INTENT (IN), TARGET :: YD
LOGICAL, OPTIONAL, INTENT (IN) :: LDCREATED
LOGICAL :: LLCREATED
LOGICAL :: LMVALUE, LNCPL2M, LNCPL2N, LNCPL4M, LNCPL4N, LNCPLM, LNCPLN, LNESM0, LNESM0G, LNESPZERO
LOGICAL :: LNPME, LNPNE, LRLEPDIM, LRLEPDIN, LRLEPINM, LRLEPINN

LLCREATED = .FALSE.
IF (PRESENT (LDCREATED)) THEN
  LLCREATED = LDCREATED
ENDIF
IF (.NOT. LLCREATED) THEN
  !$acc enter data create (YD)
  !$acc update device (YD)
ENDIF
LNCPL2M = ALLOCATED (YD%NCPL2M)
IF (LNCPL2M) THEN
  !$acc enter data create (YD%NCPL2M)
  !$acc update device (YD%NCPL2M)
  !$acc enter data attach (YD%NCPL2M)
ENDIF

LNCPL4M = ALLOCATED (YD%NCPL4M)
IF (LNCPL4M) THEN
  !$acc enter data create (YD%NCPL4M)
  !$acc update device (YD%NCPL4M)
  !$acc enter data attach (YD%NCPL4M)
ENDIF

LNCPLM = ALLOCATED (YD%NCPLM)
IF (LNCPLM) THEN
  !$acc enter data create (YD%NCPLM)
  !$acc update device (YD%NCPLM)
  !$acc enter data attach (YD%NCPLM)
ENDIF

LNCPL2N = ALLOCATED (YD%NCPL2N)
IF (LNCPL2N) THEN
  !$acc enter data create (YD%NCPL2N)
  !$acc update device (YD%NCPL2N)
  !$acc enter data attach (YD%NCPL2N)
ENDIF

LNCPL4N = ALLOCATED (YD%NCPL4N)
IF (LNCPL4N) THEN
  !$acc enter data create (YD%NCPL4N)
  !$acc update device (YD%NCPL4N)
  !$acc enter data attach (YD%NCPL4N)
ENDIF

LNCPLN = ALLOCATED (YD%NCPLN)
IF (LNCPLN) THEN
  !$acc enter data create (YD%NCPLN)
  !$acc update device (YD%NCPLN)
  !$acc enter data attach (YD%NCPLN)
ENDIF

LRLEPDIN = ALLOCATED (YD%RLEPDIN)
IF (LRLEPDIN) THEN
  !$acc enter data create (YD%RLEPDIN)
  !$acc update device (YD%RLEPDIN)
  !$acc enter data attach (YD%RLEPDIN)
ENDIF

LRLEPINN = ALLOCATED (YD%RLEPINN)
IF (LRLEPINN) THEN
  !$acc enter data create (YD%RLEPINN)
  !$acc update device (YD%RLEPINN)
  !$acc enter data attach (YD%RLEPINN)
ENDIF

LRLEPDIM = ALLOCATED (YD%RLEPDIM)
IF (LRLEPDIM) THEN
  !$acc enter data create (YD%RLEPDIM)
  !$acc update device (YD%RLEPDIM)
  !$acc enter data attach (YD%RLEPDIM)
ENDIF

LRLEPINM = ALLOCATED (YD%RLEPINM)
IF (LRLEPINM) THEN
  !$acc enter data create (YD%RLEPINM)
  !$acc update device (YD%RLEPINM)
  !$acc enter data attach (YD%RLEPINM)
ENDIF

LNESM0 = ALLOCATED (YD%NESM0)
IF (LNESM0) THEN
  !$acc enter data create (YD%NESM0)
  !$acc update device (YD%NESM0)
  !$acc enter data attach (YD%NESM0)
ENDIF

LNESPZERO = ALLOCATED (YD%NESPZERO)
IF (LNESPZERO) THEN
  !$acc enter data create (YD%NESPZERO)
  !$acc update device (YD%NESPZERO)
  !$acc enter data attach (YD%NESPZERO)
ENDIF

LNESM0G = ALLOCATED (YD%NESM0G)
IF (LNESM0G) THEN
  !$acc enter data create (YD%NESM0G)
  !$acc update device (YD%NESM0G)
  !$acc enter data attach (YD%NESM0G)
ENDIF

LNPME = ALLOCATED (YD%NPME)
IF (LNPME) THEN
  !$acc enter data create (YD%NPME)
  !$acc update device (YD%NPME)
  !$acc enter data attach (YD%NPME)
ENDIF

LNPNE = ALLOCATED (YD%NPNE)
IF (LNPNE) THEN
  !$acc enter data create (YD%NPNE)
  !$acc update device (YD%NPNE)
  !$acc enter data attach (YD%NPNE)
ENDIF

LMVALUE = ALLOCATED (YD%MVALUE)
IF (LMVALUE) THEN
  !$acc enter data create (YD%MVALUE)
  !$acc update device (YD%MVALUE)
  !$acc enter data attach (YD%MVALUE)
ENDIF

END SUBROUTINE

SUBROUTINE WIPE_TLEP (YD, LDDELETED)

IMPLICIT NONE
TYPE (TLEP), INTENT (IN), TARGET :: YD
LOGICAL, OPTIONAL, INTENT (IN) :: LDDELETED
LOGICAL :: LLDELETED
LOGICAL :: LMVALUE, LNCPL2M, LNCPL2N, LNCPL4M, LNCPL4N, LNCPLM, LNCPLN, LNESM0, LNESM0G, LNESPZERO
LOGICAL :: LNPME, LNPNE, LRLEPDIM, LRLEPDIN, LRLEPINM, LRLEPINN

LNCPL2M = ALLOCATED (YD%NCPL2M)
IF (LNCPL2M) THEN
  !$acc exit data detach (YD%NCPL2M)
  !$acc exit data delete (YD%NCPL2M)
ENDIF

LNCPL4M = ALLOCATED (YD%NCPL4M)
IF (LNCPL4M) THEN
  !$acc exit data detach (YD%NCPL4M)
  !$acc exit data delete (YD%NCPL4M)
ENDIF

LNCPLM = ALLOCATED (YD%NCPLM)
IF (LNCPLM) THEN
  !$acc exit data detach (YD%NCPLM)
  !$acc exit data delete (YD%NCPLM)
ENDIF

LNCPL2N = ALLOCATED (YD%NCPL2N)
IF (LNCPL2N) THEN
  !$acc exit data detach (YD%NCPL2N)
  !$acc exit data delete (YD%NCPL2N)
ENDIF

LNCPL4N = ALLOCATED (YD%NCPL4N)
IF (LNCPL4N) THEN
  !$acc exit data detach (YD%NCPL4N)
  !$acc exit data delete (YD%NCPL4N)
ENDIF

LNCPLN = ALLOCATED (YD%NCPLN)
IF (LNCPLN) THEN
  !$acc exit data detach (YD%NCPLN)
  !$acc exit data delete (YD%NCPLN)
ENDIF

LRLEPDIN = ALLOCATED (YD%RLEPDIN)
IF (LRLEPDIN) THEN
  !$acc exit data detach (YD%RLEPDIN)
  !$acc exit data delete (YD%RLEPDIN)
ENDIF

LRLEPINN = ALLOCATED (YD%RLEPINN)
IF (LRLEPINN) THEN
  !$acc exit data detach (YD%RLEPINN)
  !$acc exit data delete (YD%RLEPINN)
ENDIF

LRLEPDIM = ALLOCATED (YD%RLEPDIM)
IF (LRLEPDIM) THEN
  !$acc exit data detach (YD%RLEPDIM)
  !$acc exit data delete (YD%RLEPDIM)
ENDIF

LRLEPINM = ALLOCATED (YD%RLEPINM)
IF (LRLEPINM) THEN
  !$acc exit data detach (YD%RLEPINM)
  !$acc exit data delete (YD%RLEPINM)
ENDIF

LNESM0 = ALLOCATED (YD%NESM0)
IF (LNESM0) THEN
  !$acc exit data detach (YD%NESM0)
  !$acc exit data delete (YD%NESM0)
ENDIF

LNESPZERO = ALLOCATED (YD%NESPZERO)
IF (LNESPZERO) THEN
  !$acc exit data detach (YD%NESPZERO)
  !$acc exit data delete (YD%NESPZERO)
ENDIF

LNESM0G = ALLOCATED (YD%NESM0G)
IF (LNESM0G) THEN
  !$acc exit data detach (YD%NESM0G)
  !$acc exit data delete (YD%NESM0G)
ENDIF

LNPME = ALLOCATED (YD%NPME)
IF (LNPME) THEN
  !$acc exit data detach (YD%NPME)
  !$acc exit data delete (YD%NPME)
ENDIF

LNPNE = ALLOCATED (YD%NPNE)
IF (LNPNE) THEN
  !$acc exit data detach (YD%NPNE)
  !$acc exit data delete (YD%NPNE)
ENDIF

LMVALUE = ALLOCATED (YD%MVALUE)
IF (LMVALUE) THEN
  !$acc exit data detach (YD%MVALUE)
  !$acc exit data delete (YD%MVALUE)
ENDIF

LLDELETED = .FALSE.
IF (PRESENT (LDDELETED)) THEN
  LLDELETED = LDDELETED
ENDIF
IF (.NOT. LLDELETED) THEN
  !$acc exit data delete (YD)
ENDIF
END SUBROUTINE



END MODULE
