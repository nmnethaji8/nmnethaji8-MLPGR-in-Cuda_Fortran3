SUBROUTINE COLLISION_3D_SHA(LNODE,NODN,DT,NODEID,COORX,COORY,COORZ,&
   DDR,SNX,SNY,SNZ)
   USE MLPGKINE
   USE NEIGHNODES
   IMPLICIT NONE

   INTEGER(KIND=4),INTENT(IN)::LNODE,NODN,NODEID(-7:NODN)
   REAL(KIND=8),INTENT(IN)::COORX(LNODE),COORY(LNODE),COORZ(LNODE)
   REAL(KIND=8),INTENT(IN)::DDR(NODN),DT
   REAL(KIND=8),INTENT(IN)::SNX(LNODE),SNY(LNODE),SNZ(LNODE)

   INTEGER,PARAMETER::MAXCOLNE=50
   INTEGER::I,I2,J,INOD,NWATER
   INTEGER::NN,ND(MAXCOLNE),NCOLL
   REAL(KIND=8)::VIJ(3),RIJ(3),EIJ(3),DR,DRSQ,VIJM,MI,MJ,DDRINOD,DV(3)
   REAL(KIND=8)::COLLRAD,COEFREST,COEFRAD1,COEFRAD2,MTOT
   LOGICAL::COLLWALL1,COLLWALL2
   TYPE(MLPGCONN)::ILINK

   WRITE(8,'(" [MSG] ENTERING COLLISION_3D_SHA")')

   COEFREST=1D0    !! COEFF OF RESTITUTION FOR COLLISION
   COEFRAD1=0.30D0 !! WATER WATER
   COEFRAD2=0.4D0 !! WATER WALL

   UX(:,3)=UX(:,1)
   UY(:,3)=UY(:,1)
   UZ(:,3)=UZ(:,1)
   NWATER=NODEID(-2)
   NCOLL=0

   DO INOD=1,NWATER

      NN=0
      COLLWALL1=.FALSE.
      ILINK=NLINK(INOD)
      MI=ROU(INOD)
      DDRINOD=DDR(INOD)

      DO I2=1,ILINK%I(0)
         I=ILINK%I(I2)

         RIJ(1)=COORX(INOD)-COORX(I)
         RIJ(2)=COORY(INOD)-COORY(I)
         RIJ(3)=COORZ(INOD)-COORZ(I)
         DR=DSQRT(RIJ(1)**2 + RIJ(2)**2 + RIJ(3)**2)

         IF(I.LE.NWATER)THEN
            COLLRAD=COEFRAD1*DDRINOD
            COLLWALL2=.FALSE.
         ELSE
            COLLRAD=COEFRAD2*DDRINOD
            COLLWALL2=.TRUE.
         ENDIF

         IF(DR.LT.COLLRAD)THEN
            VIJ(1)=UX(INOD,3)-UX(I,3)
            VIJ(2)=UY(INOD,3)-UY(I,3)
            VIJ(3)=UZ(INOD,3)-UZ(I,3)

            IF(COLLWALL2)THEN
               EIJ(1)=-SNX(I)
               EIJ(2)=-SNY(I)
               EIJ(3)=-SNZ(I)
            ELSE
               EIJ=RIJ/DR
            ENDIF

            VIJM=(VIJ(1)*EIJ(1))+(VIJ(2)*EIJ(2))+(VIJ(3)*EIJ(3))
            IF(VIJM.LT.0D0)THEN
               IF(-VIJM*DT.LT.1e-6)CYCLE
               IF(.NOT.COLLWALL1) COLLWALL1=COLLWALL2
               NN=NN+1
               IF(NN.GT.MAXCOLNE)THEN
                  WRITE(8,'("     [ERR] INCREASE MAXCOLNE FOR NODE")'),INOD
                  STOP
               ENDIF
               ND(NN)=I
            ENDIF
         ENDIF
      ENDDO

      IF(NN.EQ.0)CYCLE
      NCOLL=NCOLL+1

      DV=0D0
      MTOT=0D0
      IF(COLLWALL1) MI=0D0
      DO I2=1,NN
         I=ND(I2)
         IF(COLLWALL1)THEN
            IF(I.LE.NWATER)CYCLE
            MJ=ROU(I)
            VIJ(1)=UX(INOD,3)-UX(I,3)
            VIJ(2)=UY(INOD,3)-UY(I,3)
            VIJ(3)=UZ(INOD,3)-UZ(I,3)
            EIJ(1)=-SNX(I)
            EIJ(2)=-SNY(I)
            EIJ(3)=-SNZ(I)
            VIJM=(VIJ(1)*EIJ(1))+(VIJ(2)*EIJ(2))+(VIJ(3)*EIJ(3))

            VIJM=VIJM*MJ*(1D0+COEFREST)
            DV=DV+VIJM*EIJ
            MTOT=MTOT+MJ

         ELSE
            MJ=ROU(I)
            RIJ(1)=COORX(INOD)-COORX(I)
            RIJ(2)=COORY(INOD)-COORY(I)
            RIJ(3)=COORZ(INOD)-COORZ(I)
            DRSQ=(RIJ(1)**2 + RIJ(2)**2 + RIJ(3)**2)
            VIJ(1)=UX(INOD,3)-UX(I,3)
            VIJ(2)=UY(INOD,3)-UY(I,3)
            VIJ(3)=UZ(INOD,3)-UZ(I,3)
            VIJM=(VIJ(1)*RIJ(1))+(VIJ(2)*RIJ(2))+(VIJ(3)*RIJ(3))

            VIJM=VIJM/DRSQ*MJ*(1D0+COEFREST)
            DV=DV+VIJM*RIJ
            MTOT=MTOT+MJ
         ENDIF
      ENDDO

      IF(MTOT.LE.0D0)THEN
         WRITE(8,'("     [COL] ",2I10,L2)')INOD,NN,COLLWALL1
         WRITE(8,'("     [ERR] ERROR IN COLLISION, MTOT < 0",I10,F15.6)')&
            INOD,MTOT
      ENDIF
      DV=DV/(MTOT+MI)

      WRITE(8,'("     [COL] ",3I10,L2)')INOD,NODEID(INOD),NN,COLLWALL1
      WRITE(8,'("     [---] COR ",3F15.6)')COORX(INOD),&
         COORY(INOD),COORZ(INOD)
      WRITE(8,'("     [---] DV  ",3F15.6)')-DV
      DO I2=1,NN
         WRITE(8,'(A15,2I10)')"[---] NGH ",ND(I2),NODEID(ND(I2))
      ENDDO

      UX(INOD,1)=UX(INOD,1)-DV(1)
      UY(INOD,1)=UY(INOD,1)-DV(2)
      UZ(INOD,1)=UZ(INOD,1)-DV(3)
   ENDDO

   WRITE(8,'("     [COC] COLL COUNT ",I10)')NCOLL

   WRITE(8,'(" [MSG] EXITING COLLISION_3D_SHA")')
   WRITE(8,*)

END SUBROUTINE COLLISION_3D_SHA
