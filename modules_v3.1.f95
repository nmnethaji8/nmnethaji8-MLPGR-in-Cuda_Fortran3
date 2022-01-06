!!----------------------------DOMAINMOD----------------------------!!
MODULE DOMAINMOD
   IMPLICIT NONE

   TYPE :: DOMAINDATA
      INTEGER(KIND=4)::NP,NPS
      REAL(KIND=8),MANAGED,ALLOCATABLE::CX(:),CY(:),CZ(:)
      REAL(KIND=8),MANAGED,ALLOCATABLE::UX(:),UY(:),UZ(:)
      REAL(KIND=8),MANAGED,ALLOCATABLE::PR(:),ET(:),Z0(:)
   CONTAINS
      PROCEDURE :: INITCOOR
      PROCEDURE :: INITALL
      PROCEDURE :: ENDDOMAINDATA

   END TYPE DOMAINDATA

CONTAINS

   SUBROUTINE INITCOOR(THIS,INP)
      IMPLICIT NONE

      CLASS(DOMAINDATA),INTENT(INOUT)::THIS
      INTEGER(KIND=4),INTENT(IN)::INP

      THIS%NP=INP
      THIS%NPS=INP
      ALLOCATE(THIS%CX(INP), THIS%CY(INP), THIS%CZ(INP))
      ALLOCATE(THIS%ET(INP), THIS%Z0(INP))

   END SUBROUTINE INITCOOR

   SUBROUTINE INITALL(THIS,INP)
      IMPLICIT NONE

      CLASS(DOMAINDATA),INTENT(INOUT)::THIS
      INTEGER(KIND=4),INTENT(IN)::INP

      THIS%NP=INP
      THIS%NPS=INP
      ALLOCATE(THIS%CX(INP), THIS%CY(INP), THIS%CZ(INP))
      ALLOCATE(THIS%UX(INP), THIS%UY(INP), THIS%UZ(INP))
      ALLOCATE(THIS%PR(INP), THIS%ET(INP), THIS%Z0(INP))

   END SUBROUTINE INITALL

   SUBROUTINE ENDDOMAINDATA(THIS)
      IMPLICIT NONE

      CLASS(DOMAINDATA),INTENT(INOUT)::THIS

      IF(ALLOCATED(THIS%CX)) DEALLOCATE(THIS%CX)
      IF(ALLOCATED(THIS%CY)) DEALLOCATE(THIS%CY)
      IF(ALLOCATED(THIS%CZ)) DEALLOCATE(THIS%CZ)

      IF(ALLOCATED(THIS%UX)) DEALLOCATE(THIS%UX)
      IF(ALLOCATED(THIS%UY)) DEALLOCATE(THIS%UY)
      IF(ALLOCATED(THIS%UZ)) DEALLOCATE(THIS%UZ)

      IF(ALLOCATED(THIS%PR)) DEALLOCATE(THIS%PR)
      IF(ALLOCATED(THIS%ET)) DEALLOCATE(THIS%ET)
      IF(ALLOCATED(THIS%Z0)) DEALLOCATE(THIS%Z0)


   END SUBROUTINE ENDDOMAINDATA

END MODULE DOMAINMOD
!!--------------------------END DOMAINMOD--------------------------!!

!!---------------------------MLPGSTORAGE---------------------------!!

MODULE MLPGSTORAGE

   IMPLICIT NONE

   INTEGER,MANAGED,ALLOCATABLE:: LINKTAB(:),IVV(:)
   DOUBLEPRECISION,MANAGED,ALLOCATABLE::SKK(:)
   REAL(KIND=8),MANAGED,ALLOCATABLE,TARGET::PTCSR(:)
   REAL(KIND=8),MANAGED,ALLOCATABLE,TARGET::FBCSR(:),SKKCSR(:)
   INTEGER(KIND=4),MANAGED,ALLOCATABLE,TARGET::IVVCSR(:),JVVCSR(:)
   !REAL(KIND=8),MANAGED,ALLOCATABLE::WORK(:,:) !!(NODEID(0),NTHREADS)
   !LOGICAL(KIND=1),MANAGED,ALLOCATABLE::NWORK(:,:) !!(NODEID(0),NTHREADS)

END MODULE MLPGSTORAGE

!!-------------------------END MLPGSTORAGE-------------------------!!

!!----------------------------MLPGKINE-----------------------------!!

MODULE MLPGKINE

   IMPLICIT NONE
   REAL(KIND=8),MANAGED,ALLOCATABLE::UX(:,:),UY(:,:),UZ(:,:)
   REAL(KIND=8),MANAGED,ALLOCATABLE::PPX(:,:),PPY(:,:),PPZ(:,:)
   REAL(KIND=8),MANAGED,ALLOCATABLE::ROU(:),R(:),CC(:),R0(:)

END MODULE MLPGKINE

!!--------------------------END MLPGKINE---------------------------!!


!!-----------------------------PROBES------------------------------!!
MODULE PROBESMOD
   IMPLICIT NONE

   TYPE :: PROBETYP
      INTEGER(KIND=4)::NP, NVAR, FILE
      REAL(KIND=8),MANAGED,ALLOCATABLE::XYZ(:,:)
      REAL(KIND=8),MANAGED,ALLOCATABLE::DATA(:,:)
   CONTAINS
      PROCEDURE :: INITPROBE
   END TYPE PROBETYP
CONTAINS

   SUBROUTINE INITPROBE(P,NP,NVAR,FILE)
      IMPLICIT NONE
      CLASS(PROBETYP),INTENT(INOUT)::P
      INTEGER(KIND=4),INTENT(IN)::NP,NVAR,FILE

      P%NP = NP
      P%NVAR = NVAR
      P%FILE = FILE
      IF(ALLOCATED(P%XYZ)) DEALLOCATE(P%XYZ)
      IF(ALLOCATED(P%DATA)) DEALLOCATE(P%DATA)
      IF((P%NP>0).AND.(P%NVAR>0))THEN
         ALLOCATE(P%XYZ(NP,3), P%DATA(NP,NVAR))
      ENDIF
      P%XYZ=0D0
      P%DATA=0D0

   END SUBROUTINE INITPROBE
END MODULE PROBESMOD
!!---------------------------END PROBES----------------------------!!

!!--------------------------FLUXPLANEMOD---------------------------!!
MODULE FLUXPLANEMOD
   IMPLICIT NONE

   !! VERTICAL PLANE ACROSS WHICH
   !! MASS AND MOMENTUN FLUX IS MEASURED

   TYPE :: FLUXPLANE
      !! ENSURE NXY AND NZ ARE ODD AND >3 FOR SIMPSON'S INTEGRATION

      INTEGER(KIND=4)::NXY, NZ, NP
      REAL(KIND=8),MANAGED,ALLOCATABLE::XFS(:), YFS(:), ZFS(:)
      REAL(KIND=8),MANAGED,ALLOCATABLE::X(:), Y(:), Z(:), WEI(:)
      REAL(KIND=8),MANAGED,ALLOCATABLE::U(:), V(:), W(:), P(:)

      !! BOTTOM-LEFT AND TOP-RIGHT COORD
      REAL(KIND=8),MANAGED::BL(3),TR(3)

      !! OUTWARD UNIT NORMAL
      REAL(KIND=8),MANAGED::NORM(3)

   CONTAINS
      PROCEDURE :: SETFLUXPLANE
      PROCEDURE :: GENPLANEPOI
      PROCEDURE :: CALCMASFLUX

   END TYPE FLUXPLANE

CONTAINS

   SUBROUTINE SETFLUXPLANE(F, NXY, NZ, BL, TR)

      IMPLICIT NONE

      CLASS(FLUXPLANE),INTENT(INOUT)::F
      INTEGER(KIND=4),INTENT(IN)::NXY, NZ
      REAL(KIND=8),INTENT(IN)::BL(3), TR(3)

      INTEGER(KIND=4)::I, J
      REAL(KIND=8)::R1, R2, DX, DY, DR, TANGENT(3)

      F%NXY = NXY
      F%NZ = NZ

      ! DO NOT USE NXY AND NZ ANYWHERE
      ! USE F%NXY AND F%NZ, BECAUSE THE VALUES MAY BE MODIFIED

      !! FOR SIMPSON'S INTEGRATION
      ! ENSURE GT 3
      IF(F%NXY .LT. 3) F%NXY = 3
      IF(F%NZ .LT. 3) F%NZ = 3

      ! ENSURE ODD
      IF( MOD(F%NXY,2) .EQ. 0 ) F%NXY = F%NXY + 1
      IF( MOD(F%NZ,2) .EQ. 0 ) F%NZ = F%NZ + 1

      IF(F%NXY .NE. NXY)THEN
         WRITE(8, '(" [INF] FLUXPLANE OBJECT NXY MODIFIED FOR SIMPSON INTEG")')
         WRITE(8, '(" [---] FROM ",I10," TO ",I10)')NXY, F%NXY
      ENDIF

      IF(F%NZ .NE. NZ)THEN
         WRITE(8, '(" [INF] FLUXPLANE OBJECT NZ MODIFIED FOR SIMPSON INTEG")')
         WRITE(8, '(" [---] FROM ",I10," TO ",I10)')NZ, F%NZ
      ENDIF

      F%NP = F%NXY * F%NZ

      ALLOCATE(F%XFS(F%NXY), F%YFS(F%NXY), F%ZFS(F%NXY))
      ALLOCATE(F%X(F%NP), F%Y(F%NP), F%Z(F%NP), F%WEI(F%NP))
      ALLOCATE(F%U(F%NP), F%V(F%NP), F%W(F%NP), F%P(F%NP))

      F%BL = BL
      F%TR = TR

      !! CALCULATE OUTWARD UNIT NORMAL FOR VERTICAL PLANE
      DX = F%TR(1) - F%BL(1)
      DY = F%TR(2) - F%BL(2)
      DR = DSQRT(DX**2 + DY**2)

      F%NORM(1) = -DY/DR
      F%NORM(2) = DX/DR
      F%NORM(3) = 0D0

      TANGENT(1) = DX/DR
      TANGENT(2) = DY/DR
      TANGENT(3) = 0D0

      R1 = DR / (F%NXY - 1)

      DO I = 0, (F%NXY-1)
         F%XFS(I+1) = F%BL(1) + (R1*I)*TANGENT(1)
         F%YFS(I+1) = F%BL(2) + (R1*I)*TANGENT(2)
      ENDDO

   END SUBROUTINE SETFLUXPLANE

   SUBROUTINE GENPLANEPOI(F)
      IMPLICIT NONE

      CLASS(FLUXPLANE),INTENT(INOUT)::F
      INTEGER(KIND=4)::I, J, I2
      REAL(KIND=8)::X, Y, Z0, Z, DXY, DZ, W1, W2

      IF( (MOD(F%NXY,2) .EQ. 0) .OR. (MOD(F%NZ,2) .EQ. 0) ) THEN
         WRITE(8, '(" [ERR] ENSURE F%NXY AND F%NZ ARE BOTH ODD FOR SIMPSON INTEG")')
         WRITE(8, '(" [---] F%NXY, F%NZ ")') F%NXY, F%NZ
         STOP
      ENDIF

      W1 = F%TR(1) - F%BL(1)
      W2 = F%TR(2) - F%BL(2)
      DXY = DSQRT(W1**2 + W2**2) / (F%NXY - 1)
      Z0 = F%BL(3)

      DO I = 1, F%NXY
         X = F%XFS(I)
         Y = F%YFS(I)
         DZ = (F%ZFS(I) - Z0) / (F%NZ - 1)

         ! XY WEIGHT SIMPSON 3 POINT
         IF( (I.EQ.1) .OR. (I.EQ.F%NXY) )THEN
            W1 = 1D0
         ELSEIF( MOD(I,2).EQ.0 )THEN
            W1 = 4D0
         ELSE
            W1 = 2D0
         ENDIF
         W1 = W1*DXY/3D0

         DO J = 1, F%NZ
            I2 = (I-1)*F%NZ + J
            F%X(I2) = X
            F%Y(I2) = Y
            F%Z(I2) = Z0 + (J-1)*DZ

            ! Z WEIGHT SIMPSON 3 POINT
            IF( (J.EQ.1) .OR. (J.EQ.F%NZ) )THEN
               W2 = 1D0
            ELSEIF( MOD(J,2).EQ.0 )THEN
               W2 = 4D0
            ELSE
               W2 = 2D0
            ENDIF
            W2 = W2*DZ/3D0

            F%WEI(I2) = ABS(W1*W2)
         ENDDO
      ENDDO

   END SUBROUTINE GENPLANEPOI

   SUBROUTINE CALCMASFLUX(F, RHOW, MASFLUX)
      IMPLICIT NONE

      CLASS(FLUXPLANE),INTENT(INOUT)::F
      REAL(KIND=8),INTENT(IN)::RHOW
      REAL(KIND=8),INTENT(OUT)::MASFLUX

      INTEGER(KIND=4)::I
      REAL(KIND=8)::VN

      MASFLUX = 0D0
      DO I = 1, F%NP
         VN = F%U(I)*F%NORM(1) + F%V(I)*F%NORM(2) + F%W(I)*F%NORM(3)
         MASFLUX = MASFLUX + (F%WEI(I) * VN * RHOW)
      ENDDO

   END SUBROUTINE CALCMASFLUX

END MODULE FLUXPLANEMOD
!!------------------------END FLUXPLANEMOD-------------------------!!

!!---------------------------VOLPLANEMOD---------------------------!!
MODULE VOLPLANEMOD
   IMPLICIT NONE

   !! Horizontal plane XY
   !! define the XY points for which the FS is obtained
   !! The volume under the FS is calculated then using SIMPSON integ

   TYPE:: VOLPLANE
      INTEGER(KIND=4)::NX, NY, NP
      REAL(KIND=8),MANAGED::BL(2), TR(2) !XY ONLY
      REAL(KIND=8),MANAGED,ALLOCATABLE::XFS(:), YFS(:), ZFS(:), WEI(:)

   CONTAINS
      PROCEDURE :: SETVOLPLANE
   END TYPE VOLPLANE

CONTAINS

   SUBROUTINE SETVOLPLANE(F, NX, NY, BL, TR)
      IMPLICIT NONE

      CLASS(VOLPLANE),INTENT(INOUT)::F
      INTEGER(KIND=4),INTENT(IN)::NX, NY
      REAL(KIND=8),INTENT(IN)::BL(2), TR(2)

      INTEGER(KIND=4)::I, J, I2
      REAL(KIND=8)::DX, DY, W1, W2

      F%NX = NX
      F%NY = NY
      F%BL = BL
      F%TR = TR

      !! FOR SIMPSON'S INTEGRATION
      ! ENSURE GT 3
      IF(F%NX .LT. 3) F%NX = 3
      IF(F%NY .LT. 3) F%NY = 3

      ! ENSURE ODD
      IF( MOD(F%NX,2) .EQ. 0 ) F%NX = F%NX + 1
      IF( MOD(F%NY,2) .EQ. 0 ) F%NY = F%NY + 1

      IF(F%NX .NE. NX)THEN
         WRITE(8, '(" [INF] VOLPLANE OBJECT NX MODIFIED FOR SIMPSON INTEG")')
         WRITE(8, '(" [---] FROM ",I10," TO ",I10)')NX, F%NX
      ENDIF

      IF(F%NY .NE. NY)THEN
         WRITE(8, '(" [INF] VOLPLANE OBJECT NY MODIFIED FOR SIMPSON INTEG")')
         WRITE(8, '(" [---] FROM ",I10," TO ",I10)')NY, F%NY
      ENDIF

      F%NP = F%NX * F%NY
      ALLOCATE(F%XFS(F%NP), F%YFS(F%NP), F%ZFS(F%NP), F%WEI(F%NP))

      DX = (TR(1) - BL(1)) / (F%NX - 1)
      DY = (TR(2) - BL(2)) / (F%NY - 1)

      DO I = 1, F%NX

         ! X WEIGHT FOR SIMPSON 3 POINT
         IF( (I.EQ.1) .OR. (I.EQ.F%NX) )THEN
            W1 = 1D0
         ELSEIF( MOD(I,2).EQ.0 )THEN
            W1 = 4D0
         ELSE
            W1 = 2D0
         ENDIF
         W1 = W1*DX/3D0

         DO J = 1, F%NY
            ! Y WEIGHT FOR SIMPSON 3 POINT
            IF( (J.EQ.1) .OR. (J.EQ.F%NY) )THEN
               W2 = 1D0
            ELSEIF( MOD(J,2).EQ.0 )THEN
               W2 = 4D0
            ELSE
               W2 = 2D0
            ENDIF
            W2 = W2*DY/3D0

            I2 = (I-1)*F%NY + J
            F%XFS(I2) = BL(1) + DX*(I-1)
            F%YFS(I2) = BL(2) + DY*(J-1)
            F%WEI(I2) = ABS(W1*W2)
         ENDDO
      ENDDO

   END SUBROUTINE SETVOLPLANE
END MODULE VOLPLANEMOD
!!-------------------------END VOLPLANEMOD-------------------------!!

!!---------------------------NEIGHNODES----------------------------!!
MODULE NEIGHNODES
   IMPLICIT NONE
   TYPE MLPGCONN
      INTEGER(KIND=4),MANAGED:: I(0:300)
      !REAL(KIND=8),MANAGED,ALLOCATABLE :: I(:)
      !INTEGER(KIND=4),MANAGED,ALLOCATABLE :: I(:)
   END TYPE MLPGCONN

   TYPE(MLPGCONN),MANAGED,ALLOCATABLE:: NLINK(:)

   TYPE MLPGCONN2
      INTEGER(KIND=4),ALLOCATABLE :: I(:)
   END TYPE MLPGCONN2

END MODULE NEIGHNODES
!!-------------------------END NEIGHNODES--------------------------!!
