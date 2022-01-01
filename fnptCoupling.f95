MODULE FNPTCPLMOD
    IMPLICIT NONE

    TYPE :: FNPTCPLTYP

    INTEGER(KIND=4)::FILE, NX, NY, NN
    REAL(KIND=8)::DDL, X0, RLXLEN
    CHARACTER(LEN=256)::FILENAME
    INTEGER(KIND=4)::MLNP
    INTEGER(KIND=4),MANAGED,ALLOCATABLE::MLPOI(:)
    REAL(KIND=8),MANAGED,ALLOCATABLE::MLVAL(:,:)
    REAL(KIND=8),MANAGED,ALLOCATABLE::FNCO(:,:), FNVA(:,:), FNFS(:,:)

    CONTAINS

    PROCEDURE :: SETCOUPLING
    PROCEDURE :: READFNPT

    END TYPE FNPTCPLTYP
    
    CONTAINS

    !!------------------------SETCOUPLING------------------------!!
    SUBROUTINE SETCOUPLING(F,NODN,NODEID,COORX)
        IMPLICIT NONE

        CLASS(FNPTCPLTYP),INTENT(INOUT)::F
        INTEGER(KIND=4),INTENT(IN)::NODN, NODEID(-7:NODN)
        REAL(KIND=8),INTENT(IN)::COORX(NODN) !COORX(:,1)

        INTEGER(KIND=4)::IX,II,IY
        REAL(KIND=8)::TMPR1

        ALLOCATE( F%FNCO(F%NN,2), F%FNVA(F%NN,6), F%FNFS(F%NX,2) )

        F%MLNP=0
        !! Wavemaker nodes
        DO IX=NODEID(-2)+1,NODEID(-3)
            F%MLNP = F%MLNP + 1
        ENDDO

        !! Wall nodes
        DO IX=NODEID(-3)+1,NODEID(-1)
            IF(NODEID(IX).EQ.2) CYCLE !DONT INCLUDE BOTTOM NODES
            IF(COORX(IX) .LT. F%RLXLEN)THEN
                F%MLNP = F%MLNP + 1
            ENDIF
        ENDDO

        !! Fluid nodes
        DO IX=1,NODEID(-2)
            IF(COORX(IX) .LT. F%RLXLEN)THEN
                F%MLNP = F%MLNP + 1
            ENDIF
        ENDDO

        ALLOCATE( F%MLPOI( F%MLNP ), F%MLVAL( F%MLNP, 7 ) )
        F%MLVAL=0D0

        ! FN_VA(:,1:6) and MP_VA(:,1:7) description
        ! FN_VA(:,1)    =   MP_VA(:,1)    =   pressure
        ! FN_VA(:,2)    =   MP_VA(:,2)    =   Ux  
        ! FN_VA(:,3)    =   MP_VA(:,3)    =   Uy  = 0d0
        ! FN_VA(:,4)    =   MP_VA(:,4)    =   Uz
        ! FN_VA(:,5)    =   MP_VA(:,5)    =   Ax
        ! FN_VA(:,6)    =   MP_VA(:,6)    =   Az
        !                   MP_VA(:,7)    =   weight

        II=0
        !! Wavemaker nodes
        DO IX=NODEID(-2)+1,NODEID(-3)
            II = II + 1
            F%MLPOI(II) = IX
        ENDDO

        !! Wall nodes
        DO IX=NODEID(-3)+1,NODEID(-1)
            IF(NODEID(IX).EQ.2) CYCLE !DONT INCLUDE BOTTOM NODES
            IF(COORX(IX) .LT. F%RLXLEN)THEN
                II = II + 1
                F%MLPOI(II) = IX
            ENDIF
        ENDDO

        !! Fluid nodes
        DO IX=1,NODEID(-2)
            IF(COORX(IX) .LT. F%RLXLEN)THEN
                II = II + 1
                F%MLPOI(II) = IX
            ENDIF
        ENDDO

        IF(II .NE. F%MLNP)THEN
            WRITE(8,'(" [ERR] ERROR IN SETCOUPLING")')
            STOP
        ENDIF

        WRITE(8,*)'[MSG] NUM OF OVERLAP POI',F%MLNP

        DO IX = 1, F%MLNP
            IY = F%MLPOI(IX)
            !Assuming left wall is vertical at X=0
            TMPR1 = COORX(IY) / F%RLXLEN 
            F%MLVAL(IX,7) = 1D0 - 3D0*TMPR1**2 + 2D0*TMPR1**3
            !MP_VA(IX,7)=1D0-DSIN(ATAN(1)*2D0*TMPR1)
        ENDDO
    END SUBROUTINE SETCOUPLING

    !!----------------------END SETCOUPLING----------------------!!

    !!-------------------------READFNPT--------------------------!!
    SUBROUTINE READFNPT(F,H0)
        IMPLICIT NONE

        CLASS(FNPTCPLTYP),INTENT(INOUT)::F
        REAL(KIND=8),INTENT(IN)::H0

        INTEGER(KIND=4)::IK,IX,IY
        REAL(KIND=8)::TMPR1

        ! FN_VA(:,1:6) and MP_VA(:,1:7) description
        ! FN_VA(:,1)    =   MP_VA(:,1)    =   pressure
        ! FN_VA(:,2)    =   MP_VA(:,2)    =   Ux  
        ! FN_VA(:,3)    =   MP_VA(:,3)    =   Uy  = 0d0
        ! FN_VA(:,4)    =   MP_VA(:,4)    =   Uz
        ! FN_VA(:,5)    =   MP_VA(:,5)    =   Ax
        ! FN_VA(:,6)    =   MP_VA(:,6)    =   Az
        !                   MP_VA(:,7)    =   weight

        IK=0
        DO IX = 1, F%NX
            DO IY = 1, F%NY
                IK = IK + 1
                READ(F%FILE, *)TMPR1, F%FNCO(IK,1:2), &
                F%FNVA(IK,1:2), F%FNVA(IK,4:6)
            ENDDO
            F%FNFS(IX,:)=F%FNCO(IK,:)
        ENDDO
        F%FNCO(:,1) = F%FNCO(:,1) - F%X0
        F%FNCO(:,2) = F%FNCO(:,2) + H0
        F%FNFS(:,1) = F%FNFS(:,1) - F%X0
        F%FNFS(:,2) = F%FNFS(:,2) + H0
        F%FNVA(:,1) = F%FNVA(:,1)*1000D0
        F%FNVA(:,3) = 0D0
    END SUBROUTINE READFNPT
    !!-----------------------END READFNPT------------------------!!

    !!------------------------FNPT_INTERP------------------------!!
    SUBROUTINE FNPT_INTERP(FN_NN, FN_CO, FN_VA, MP_NP, MP_POI, MP_VA, DDR)
        USE COMMONMOD
        USE CUDAFOR

        IMPLICIT NONE

        INTEGER(KIND=4),INTENT(IN)::FN_NN,MP_NP
        INTEGER(KIND=4),MANAGED,INTENT(IN)::MP_POI(MP_NP)
        REAL(KIND=8),MANAGED,INTENT(INOUT)::MP_VA(MP_NP,7)
        REAL(KIND=8),MANAGED,INTENT(IN)::FN_CO(FN_NN,2),FN_VA(FN_NN,6)
        REAL(KIND=8),MANAGED,INTENT(IN)::DDR(NODEID(0))

        INTEGER(KIND=4)::I,J,K,L,MP_I
        INTEGER(KIND=4)::NE_NP
        INTEGER(KIND=4)::NE_MA(FN_NN)
        REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
        REAL(KIND=8)::NE_VA(FN_NN,2)
        REAL(KIND=8)::NE_RAD=0.08D0

        DO MP_I=1,MP_NP
            NE_MA=0
            NE_NP=0     
            L=MP_POI(MP_I)   
            TMPR1=COORX(L,1)
            TMPR2=COORZ(L,1)
            TMPR7=0D0

            DO I=1,FN_NN
                TMPR3=DSQRT((FN_CO(I,1)-TMPR1)**2+(FN_CO(I,2)-TMPR2)**2)
                TMPR3=TMPR3/NE_RAD          
                IF(TMPR3.LE.1D0)THEN
                    NE_NP=NE_NP+1
                    NE_MA(NE_NP)=I
                    NE_VA(NE_NP,1)=TMPR3
                    IF(TMPR3.LE.0.5D0)THEN
                        NE_VA(NE_NP,2)=2D0/3D0-4D0*TMPR3**2+4D0*TMPR3**3
                        TMPR7=TMPR7+NE_VA(NE_NP,2)
                    ELSE
                        NE_VA(NE_NP,2)=4D0/3D0-4D0*TMPR3+4D0*TMPR3**2 - 4D0/3D0*TMPR3**3
                        TMPR7=TMPR7+NE_VA(NE_NP,2)
                    ENDIF
                ENDIF
            ENDDO

            NE_VA(1:NE_NP,2)=NE_VA(1:NE_NP,2)/TMPR7

            IF(NE_NP.LT.4)THEN
                WRITE(8,*)'[ERR] INCREASE FNPT NE RAD FOR',L,NODEID(L)
                WRITE(8,*)'[---] LOCATION',COORX(L,1),COORY(L,1),COORZ(L,1)
                WRITE(8,*)'[---] NUM NE',NE_NP
                STOP
            ENDIF

            TMPR1=0D0
            TMPR2=0D0
            TMPR3=0D0
            TMPR4=0D0
            TMPR5=0D0
            TMPR6=0D0
            DO I=1,NE_NP
                K=NE_MA(I)
                TMPR1=TMPR1+NE_VA(I,2)*FN_VA(K,1)
                TMPR2=TMPR2+NE_VA(I,2)*FN_VA(K,2)
                TMPR3=TMPR3+NE_VA(I,2)*FN_VA(K,3)
                TMPR4=TMPR4+NE_VA(I,2)*FN_VA(K,4)
                TMPR5=TMPR5+NE_VA(I,2)*FN_VA(K,5)
                TMPR6=TMPR6+NE_VA(I,2)*FN_VA(K,6)
            ENDDO
            MP_VA(MP_I,1)=TMPR1
            MP_VA(MP_I,2)=TMPR2
            MP_VA(MP_I,3)=TMPR3
            MP_VA(MP_I,4)=TMPR4
            MP_VA(MP_I,5)=TMPR5
            MP_VA(MP_I,6)=TMPR6        
        ENDDO
    END SUBROUTINE FNPT_INTERP
END MODULE FNPTCPLMOD
