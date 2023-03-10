MODULE UTIL_TDIM_MOD

USE YOMDIM, ONLY : TDIM

INTERFACE SAVE
MODULE PROCEDURE SAVE_TDIM
END INTERFACE

INTERFACE LOAD
MODULE PROCEDURE LOAD_TDIM
END INTERFACE

INTERFACE COPY
MODULE PROCEDURE COPY_TDIM
END INTERFACE

INTERFACE WIPE
MODULE PROCEDURE WIPE_TDIM
END INTERFACE



CONTAINS

SUBROUTINE SAVE_TDIM (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TDIM), INTENT (IN), TARGET :: YD
LOGICAL :: LNDLUNL, LNDLUXL
WRITE (KLUN) YD%NDGLG
WRITE (KLUN) YD%NDGLL
WRITE (KLUN) YD%NDGNH
WRITE (KLUN) YD%NDGSUR
WRITE (KLUN) YD%NDGSAG
WRITE (KLUN) YD%NDGSAL
WRITE (KLUN) YD%NDGSAH
WRITE (KLUN) YD%NDGSAFPH
WRITE (KLUN) YD%NDGENG
WRITE (KLUN) YD%NDGENL
WRITE (KLUN) YD%NDGENH
WRITE (KLUN) YD%NDGENFPH
WRITE (KLUN) YD%NDGUNG
WRITE (KLUN) YD%NDGUXG
WRITE (KLUN) YD%NDGUNL
WRITE (KLUN) YD%NDGUXL
WRITE (KLUN) YD%NDLON
WRITE (KLUN) YD%NDSUR1
WRITE (KLUN) YD%NSTENCILWIDE
WRITE (KLUN) YD%NDLSUR
WRITE (KLUN) YD%NDLSM
WRITE (KLUN) YD%NDLUNG
WRITE (KLUN) YD%NDLUXG
LNDLUNL = ALLOCATED (YD%NDLUNL)
WRITE (KLUN) LNDLUNL
IF (LNDLUNL) THEN
  WRITE (KLUN) LBOUND (YD%NDLUNL)
  WRITE (KLUN) UBOUND (YD%NDLUNL)
  WRITE (KLUN) YD%NDLUNL
ENDIF
LNDLUXL = ALLOCATED (YD%NDLUXL)
WRITE (KLUN) LNDLUXL
IF (LNDLUXL) THEN
  WRITE (KLUN) LBOUND (YD%NDLUXL)
  WRITE (KLUN) UBOUND (YD%NDLUXL)
  WRITE (KLUN) YD%NDLUXL
ENDIF
WRITE (KLUN) YD%NPROMA
WRITE (KLUN) YD%NPROMM
WRITE (KLUN) YD%NPROMM9
WRITE (KLUN) YD%NPROMNH
WRITE (KLUN) YD%NGPBLKS
WRITE (KLUN) YD%LOPTPROMA
WRITE (KLUN) YD%NRESOL
WRITE (KLUN) YD%NSMAX
WRITE (KLUN) YD%NMSMAX
WRITE (KLUN) YD%NVARMAX
WRITE (KLUN) YD%NSEFRE
WRITE (KLUN) YD%NSPECG
WRITE (KLUN) YD%NSPEC2G
WRITE (KLUN) YD%NSPEC
WRITE (KLUN) YD%NSPEC2
WRITE (KLUN) YD%NSPEC2MX
WRITE (KLUN) YD%NCMAX
WRITE (KLUN) YD%NUMP
WRITE (KLUN) YD%NUMCP
END SUBROUTINE

SUBROUTINE LOAD_TDIM (KLUN, YD)

IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE (TDIM), INTENT (OUT), TARGET :: YD
INTEGER :: IL2(2), IU2(2)
LOGICAL :: LNDLUNL, LNDLUXL
READ (KLUN) YD%NDGLG
READ (KLUN) YD%NDGLL
READ (KLUN) YD%NDGNH
READ (KLUN) YD%NDGSUR
READ (KLUN) YD%NDGSAG
READ (KLUN) YD%NDGSAL
READ (KLUN) YD%NDGSAH
READ (KLUN) YD%NDGSAFPH
READ (KLUN) YD%NDGENG
READ (KLUN) YD%NDGENL
READ (KLUN) YD%NDGENH
READ (KLUN) YD%NDGENFPH
READ (KLUN) YD%NDGUNG
READ (KLUN) YD%NDGUXG
READ (KLUN) YD%NDGUNL
READ (KLUN) YD%NDGUXL
READ (KLUN) YD%NDLON
READ (KLUN) YD%NDSUR1
READ (KLUN) YD%NSTENCILWIDE
READ (KLUN) YD%NDLSUR
READ (KLUN) YD%NDLSM
READ (KLUN) YD%NDLUNG
READ (KLUN) YD%NDLUXG
READ (KLUN) LNDLUNL
IF (LNDLUNL) THEN
  READ (KLUN) IL2
  READ (KLUN) IU2
  ALLOCATE (YD%NDLUNL (IL2(1):IU2(1), IL2(2):IU2(2)))
  READ (KLUN) YD%NDLUNL
ENDIF
READ (KLUN) LNDLUXL
IF (LNDLUXL) THEN
  READ (KLUN) IL2
  READ (KLUN) IU2
  ALLOCATE (YD%NDLUXL (IL2(1):IU2(1), IL2(2):IU2(2)))
  READ (KLUN) YD%NDLUXL
ENDIF
READ (KLUN) YD%NPROMA
READ (KLUN) YD%NPROMM
READ (KLUN) YD%NPROMM9
READ (KLUN) YD%NPROMNH
READ (KLUN) YD%NGPBLKS
READ (KLUN) YD%LOPTPROMA
READ (KLUN) YD%NRESOL
READ (KLUN) YD%NSMAX
READ (KLUN) YD%NMSMAX
READ (KLUN) YD%NVARMAX
READ (KLUN) YD%NSEFRE
READ (KLUN) YD%NSPECG
READ (KLUN) YD%NSPEC2G
READ (KLUN) YD%NSPEC
READ (KLUN) YD%NSPEC2
READ (KLUN) YD%NSPEC2MX
READ (KLUN) YD%NCMAX
READ (KLUN) YD%NUMP
READ (KLUN) YD%NUMCP
END SUBROUTINE


SUBROUTINE COPY_TDIM (YD, LDCREATED)

IMPLICIT NONE
TYPE (TDIM), INTENT (IN), TARGET :: YD
LOGICAL, OPTIONAL, INTENT (IN) :: LDCREATED
LOGICAL :: LLCREATED
LOGICAL :: LNDLUNL, LNDLUXL

LLCREATED = .FALSE.
IF (PRESENT (LDCREATED)) THEN
  LLCREATED = LDCREATED
ENDIF
IF (.NOT. LLCREATED) THEN
  !$acc enter data create (YD)
  !$acc update device (YD)
ENDIF























LNDLUNL = ALLOCATED (YD%NDLUNL)
IF (LNDLUNL) THEN
  !$acc enter data create (YD%NDLUNL)
  !$acc update device (YD%NDLUNL)
  !$acc enter data attach (YD%NDLUNL)
ENDIF

LNDLUXL = ALLOCATED (YD%NDLUXL)
IF (LNDLUXL) THEN
  !$acc enter data create (YD%NDLUXL)
  !$acc update device (YD%NDLUXL)
  !$acc enter data attach (YD%NDLUXL)
ENDIF




















END SUBROUTINE

SUBROUTINE WIPE_TDIM (YD, LDDELETED)

IMPLICIT NONE
TYPE (TDIM), INTENT (IN), TARGET :: YD
LOGICAL, OPTIONAL, INTENT (IN) :: LDDELETED
LOGICAL :: LLDELETED
LOGICAL :: LNDLUNL, LNDLUXL
























LNDLUNL = ALLOCATED (YD%NDLUNL)
IF (LNDLUNL) THEN
  !$acc exit data detach (YD%NDLUNL)
  !$acc exit data delete (YD%NDLUNL)
ENDIF

LNDLUXL = ALLOCATED (YD%NDLUXL)
IF (LNDLUXL) THEN
  !$acc exit data detach (YD%NDLUXL)
  !$acc exit data delete (YD%NDLUXL)
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
