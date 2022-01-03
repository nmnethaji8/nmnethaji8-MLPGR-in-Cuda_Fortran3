PROGRAM THREED_BREAKINGWAVE
   USE COMMONMOD
   USE MLPGKINE
   USE DOMAINMOD
   USE NEIGHNODES
   USE NODELINKMOD
   USE PROBESMOD
   USE FNPTCPLMOD
   USE FLUXPLANEMOD
   USE VOLPLANEMOD
   USE CUDAFOR

   IMPLICIT NONE

   TYPE :: CPUTIMETYP
      INTEGER(KIND=8),MANAGED::TI(10)
      INTEGER(KIND=8)::SYSRATEINT
      REAL(KIND=8)::SYSRATE
   END TYPE CPUTIMETYP

   INTEGER(KIND=4)::IPRINT, I_PF, I_PF1, RESFREQ
   INTEGER(KIND=4)::NSTEPS, I_CAL_V, IFSI
   INTEGER(KIND=4)::ISTEP, ISTART, STEP0, IER, NNI, NNI2
   INTEGER(KIND=4)::I, J, II, IJ, IK, IX, IY
   INTEGER(KIND=4)::ISCYL, ISGHT, ICYL, IGHT

   REAL(KIND=8)::SCALE, VCOEFF, EPS_G
   REAL(KIND=8),MANAGED,ALLOCATABLE::FB(:),PTMP(:),P(:)
   REAL(KIND=8),MANAGED,ALLOCATABLE::UN(:,:),UM(:)
   !LOGICAL FL

   REAL(KIND=8),MANAGED::TEMPR(10)
   REAL(KIND=8),MANAGED,ALLOCATABLE::PRESS_DR(:)
   REAL(KIND=8),MANAGED,ALLOCATABLE::CSUXT1(:,:),CSUYT1(:,:),CSUZT1(:,:)
   REAL(KIND=8)::SPONGEX, DU, DV, DW
   REAL(KIND=8)::TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7

   CHARACTER(LEN=256)::MESHFILE, RESUMEFILE

   !INTEGER(KIND=4)::RNODE
   INTEGER(KIND=4)::REMESHFREQ
   INTEGER(KIND=4)::BNDNP,NODN,NPOI1,NTHR,MIRRNP
   INTEGER(KIND=4)::FSNOD1,FSNOD2,ERRSOL
   INTEGER(KIND=4),MANAGED,ALLOCATABLE::DZI(:,:)
   INTEGER(KIND=4),MANAGED,ALLOCATABLE::ERRTMP(:), MMTOIM(:)

   REAL(KIND=8),MANAGED,ALLOCATABLE::BNDXY(:,:), MIRRXY(:,:)
   INTEGER(KIND=4),MANAGED,ALLOCATABLE::BNDFIX(:,:)
   INTEGER(KIND=4),MANAGED,ALLOCATABLE::BNDFS(:)
   INTEGER(KIND=4)::NFS
   REAL(KIND=8),MANAGED,ALLOCATABLE::DDR(:)
   REAL(KIND=8),MANAGED,ALLOCATABLE::XFS(:),YFS(:),ZFS(:),ZFSTMP(:)
   REAL(KIND=8),MANAGED::DOMBL(3),DOMTR(3),DOMX(2),DOMY(2),DOMZ(2)
   REAL(KIND=8)::CYLX, CYLY, CYLR ,T1 ,T2
   INTEGER(KIND=cuda_count_kind) :: val,val2

   TYPE(FNPTCPLTYP),MANAGED,ALLOCATABLE::FP
   TYPE(NODELINKTYP),MANAGED,ALLOCATABLE::MLDOM,FSDOM,BOTDOM
   TYPE(DOMAINDATA),MANAGED,ALLOCATABLE::IM,CM,MM
   TYPE(CPUTIMETYP),MANAGED,ALLOCATABLE::CPUT
   TYPE(PROBETYP),MANAGED,ALLOCATABLE::WP,PP
   TYPE(FLUXPLANE),MANAGED,ALLOCATABLE,TARGET:: FLUXIN, FLUXOT
   TYPE(FLUXPLANE),MANAGED,POINTER:: F_PT
   TYPE(VOLPLANE),MANAGED,ALLOCATABLE::VOLPL

   LOGICAL,MANAGED,ALLOCATABLE::FSSRCH(:)
   LOGICAL::RESUMECHK=.FALSE.

   I=cudaDeviceGetLimit( val, cudaLimitMallocHeapSize )
   PRINT*,"cudaLimitMallocHeapSize",val
   val=val*(2**6)
   I = cudaDeviceSetLimit(cudaLimitMallocHeapSize,val)
   PRINT*,"cudaLimitMallocHeapSize",val

   OPEN(UNIT=8, FILE="mlpgTerOut.txt")
   OPEN(UNIT=9, FILE=" /home/vsriram/MATLAB-Drive/output2.txt")

   WRITE(8,*) '*****************************************************'
   WRITE(8,*) '3D WATER WAVE PROBLEM WITH MLPGR METHOD'
   WRITE(8,*) '*****************************************************'

   WRITE(*,*) '*****************************************************'
   WRITE(*,*) '3D WATER WAVE PROBLEM WITH MLPGR METHOD'
   WRITE(*,*) '*****************************************************'

   ! READ IN THE INPUT DATA
   ALLOCATE(WP, PP, FP)
   CALL INPUT(II, H0, DDL, SCALE, KW, MBAS, DT, TOTAL_TIME,  &
      IPRINT, I_PF, I_PF1, NSTEPS, &
      I_CAL_V, VCOEFF, I_WM, IFSI, NTHR, MESHFILE, DOMX, DOMY, DOMZ, &
      CYLX, CYLY, CYLR, SPONGEX, REMESHFREQ, RESFREQ, RESUMECHK, &
      RESUMEFILE, WP, PP, FP)

   CALL SETLNODE(II)

   CALL INITCOMMONMOD

   ALLOCATE( FB(LNODE), PTMP(LNODE), P(LNODE) )
   ALLOCATE( UN(LNODE,1:3), UM(LNODE) )
   ALLOCATE( CSUXT1(LNODE,1), CSUYT1(LNODE,1), CSUZT1(LNODE,1) )
   ALLOCATE( ERRTMP(LNODE), MMTOIM(LNODE) )

   !GENERATE THE NODES
   CALL NODEGRID(MESHFILE, FSNOD1, FSNOD2,DOMX, DOMY, DOMZ, CYLX, CYLY, CYLR)

   P = 0
   FB = 0
   PTMP = 0

   NODN=NODEID(0)       !THE TOTAL NUMBER

   WRITE(8,*)'[INF] FREE SURFACE REFERENCE NODES'
   IK=FSNOD1
   WRITE(8,'(I15,3F15.6)')IK,COORX(IK,1),COORY(IK,1),COORZ(IK,1)
   IK=FSNOD2
   WRITE(8,'(I15,3F15.6)')IK,COORX(IK,1),COORY(IK,1),COORZ(IK,1)

   ALLOCATE(DDR(NODEID(0)))
   DDR=DDL

   ALLOCATE(MLDOM,FSDOM,BOTDOM)
   CALL MLDOM%INITCELL(ICELLX,ICELLY,ICELLZ,500)
   CALL FSDOM%INITCELL(ICELLX,ICELLY,1,100)
   CALL BOTDOM%INITCELL(ICELLX,ICELLY,1,100)

   !! STORING BNDNODES FOR FORCED ADJUSTMENT IN NEWCOOR
   BNDNP=NODEID(0)-NODEID(-2)
   ALLOCATE(BNDXY(BNDNP,3),BNDFIX(3,0:BNDNP))
   ALLOCATE(BNDFS(0:BNDNP))
   !! BNDFIX DEFS (STORING BOTTOM NODES)
   !! BNDFIX(1,1:BNDFIX(1,0))=(/ BOTTOM POI /)
   !! BNDFIX(2,1:BNDFIX(2,0))=(/ NEAR SIDEWALL /)
   !! BNDFIX(3,1:BNDFIX(3,0))=(/ FAR SIDEWALL /)
   BNDFIX=0
   BNDFS=0
   TMPR7=0.0001D0

   BNDXY(1:BNDNP,1)=COORX(NODEID(-2)+1:NODEID(-2)+BNDNP,1)
   BNDXY(1:BNDNP,2)=COORY(NODEID(-2)+1:NODEID(-2)+BNDNP,1)
   BNDXY(1:BNDNP,3)=COORZ(NODEID(-2)+1:NODEID(-2)+BNDNP,1)

   DO II=1,BNDNP
      IJ=NODEID(-2)+II

      IF(NODEID(IJ).EQ.2)THEN
         BNDFIX(1,0)=BNDFIX(1,0)+1
         BNDFIX(1,BNDFIX(1,0))=IJ
      ENDIF

      IF(ABS(COORZ(IJ,1)-H0).LT.TMPR7)THEN
         IF((NODEID(IJ).EQ.1).OR.(NODEID(IJ).EQ.7))THEN
            BNDFS(0)=BNDFS(0)+1
            BNDFS(BNDFS(0))=IJ
         ENDIF
      ENDIF
   ENDDO

   WRITE(8,'(" [INF] BNDFIX : ",3I15)')BNDFIX(:,0)
   WRITE(8,'(" [---] ",3I15)')BNDFIX(:,1)
   WRITE(8,'(" [INF] BNDFS  : ",I15)')BNDFS(0)

   !! MIRROR NODES
   MIRRNP=NODEID(-1)-NODEID(-7)
   ALLOCATE(MIRRXY(MIRRNP,3))
   ISCYL=NODEID(-7)
   ISGHT=NODEID(-1)+NODEID(-7)-NODEID(-2)

   DO II=1,MIRRNP
      ICYL=ISCYL+II
      IGHT=ISGHT+II
      MIRRXY(II,1)=2*COORX(ICYL,1)-COORX(IGHT,1)
      MIRRXY(II,2)=2*COORY(ICYL,1)-COORY(IGHT,1)
      MIRRXY(II,3)=COORZ(IGHT,1)
   ENDDO
   WRITE(8,*)'[INF] MIRRNP ',MIRRNP
   WRITE(8,*)'[INF] IGHT ',IGHT
   WRITE(8,*)

   ALLOCATE(XFS(NODEID(-2)),YFS(NODEID(-2)),ZFS(NODEID(-2)))
   ALLOCATE(ZFSTMP(NODEID(-2)))
   ZFSTMP=0D0
   ALLOCATE(NLINK(NODEID(0)),STAT = IER)
   CALL ALLOCATEKINE(NODEID(0))

   !! STORING INITIAL MESH FOR REMESHING
   NPOI1=NODEID(0)
   ALLOCATE(IM,CM,MM)
   CALL IM%INITCOOR(NPOI1)
   IM%CX=COORX(1:NPOI1,1)
   IM%CY=COORY(1:NPOI1,1)
   IM%CZ=COORZ(1:NPOI1,1)

   !!---------------------------------------------------!!
   ! MM ONLY STORES MOVING NODES
   ! HENCE TYPE2 AND TYPE9 ARE EXCLUDED
   CALL MM%INITALL(NPOI1)
   IK = 0
   DO II = 1,NODEID(-1)
      IF(NODEID(II).LT.0)CYCLE
      IF(NODEID(II).EQ.2)CYCLE !BOTTOM NOT INCLUDED
      IF(NODEID(II).EQ.9)CYCLE !CYL NOT INCLUDED
      IK = IK+1
      MMTOIM(IK) = II
      MM%CX(IK) = IM%CX(II)
      MM%CY(IK) = IM%CY(II)
      MM%CZ(IK) = IM%CZ(II) !THOUGH THIS CHANGES WITH EACH REMESH
   ENDDO
   MM%NP = IK

   ! GETTING LOCAL DEP AT EACH OF THE MM% POINT FOR REMESH
   IK = NODEID(-4) - NODEID(-3)
   IX = NODEID(-3)+1
   IY = NODEID(-4)
   DOMBL(1)=MINVAL(IM%CX(IX:IY))
   DOMBL(2)=MINVAL(IM%CY(IX:IY))
   DOMBL(3)=0D0
   DOMTR(1)=MAXVAL(IM%CX(IX:IY))
   DOMTR(2)=MAXVAL(IM%CY(IX:IY))
   DOMTR(3)=0D0

   CALL BOTDOM%FILLCELL(IK, IM%CX(IX:IY), IM%CY(IX:IY), &
      ZFSTMP(1:IK), RCELL*DDL, DOMBL, DOMTR)

   NPOI1=MM%NP

   !NLMAXN=13
   CALL MLPG_GET_ETA(BOTDOM, IK, IM%CX(IX:IY), IM%CY(IX:IY), &
      IM%CZ(IX:IY), NPOI1, MM%CX(1:NPOI1), MM%CY(1:NPOI1), &
      MM%Z0(1:NPOI1), ERRTMP(1:NPOI1), DDL, 26, 1)
   WRITE(8,*),"[INF] LOCAL DEP ERRSUM ",SUM(ERRTMP(1:NPOI1))
   WRITE(8,*)
   !!---------------------------------------------------!!

   !INITIALISE THE WATER DENSITY AND V AND PRESSURE
   NODN=NODEID(0)

   CALL GIVENINITIALV_P(NODN,PTMP(1:NODN))
   P=PTMP

   WRITE(8,*)'DDL = ',DDL

   CALL ALLOCATESTORAGE(NODEID(-1))

   ALLOCATE(PRESS_DR(NODEID(-1)-NODEID(-2)))
   PRESS_DR=0D0

   ALLOCATE(FSSRCH(NODEID(0)))
   DO IX=1,NODEID(0)
      FSSRCH(IX)=.TRUE.
   ENDDO

   OPEN(WP%FILE, FILE='Export//WaveProbes.txt')
   OPEN(PP%FILE, FILE='Export//PressProbes.txt')
   OPEN(604,FILE='Export//CylindForce.txt')

   !-------------------READING FNPT--------------------!!
   NPOI1 = NODEID(0)
   CALL FP%SETCOUPLING( NPOI1, NODEID(-7:NPOI1), COORX(1:NPOI1,1) )
   OPEN( FP%FILE, FILE=TRIM(FP%FILENAME) )
   !!-----------------END READING FNPT------------------!!

   OPEN(1616,FILE='Export//CylPoi.txt')
   OPEN(1617,FILE='Export//BotDryPoi.txt')

   !! Crank Nicolson Setup
   CSUXT1=0D0
   CSUYT1=0D0
   CSUZT1=0D0

   ALLOCATE(CPUT)
   CALL SYSTEM_CLOCK(CPUT%TI(1),CPUT%SYSRATEINT)
   CPUT%SYSRATE=REAL(CPUT%SYSRATEINT,8)

   !! FLUXPLANE AND VOLPLANE SETTING
   ALLOCATE(FLUXIN, FLUXOT, VOLPL)
   CALL MLPGFLUXPLANE(DOMX, DOMY, DOMZ, DDL, &
      FLUXIN, FLUXOT, VOLPL)

   ISTART=0  !IF FRESH START
   STEP0 = ISTART
   ! TO RESUME FROM STORED STATE
   IF(RESUMECHK)THEN
      CALL READRESUME(RESUMEFILE, ISTEP, DDR, FSSRCH , P)

      !Reading the FNPT file till the resume ISTEP
      DO I = STEP0, ISTEP
         CALL FP%READFNPT( H0 )
      ENDDO

      STEP0 = ISTEP + 1
   ENDIF
   !WRITE(9,*),DOMBL(:),DOMTR(:),MM%CX(:),MM%CY(:),MM%CZ(:)

   !! LOOP ON TIME STEP START
   DO ISTEP = STEP0, NSTEPS

      IF ( (ISTEP.EQ.NSTEPS).AND. (ISTEP.NE.ISTART) )THEN
         CALL WRITERESUME((ISTEP-1), DDR, FSSRCH , P)
         EXIT
      ENDIF

      WRITE(8,*)
      WRITE(8,'(A)')"++++++ MLPGR FOR 3D BREAKING WAVE ++++++"
      WRITE(8,'(A,I10,A)')" +++++++ TIME STEP = ",ISTEP,"+++++++"
      WRITE(*,*)
      WRITE(*,'(A)')"++++++ MLPGR FOR 3D BREAKING WAVE ++++++"
      WRITE(*,'(A,I10,A)')" +++++++ TIME STEP = ",ISTEP," +++++++"

      CALL SYSTEM_CLOCK(CPUT%TI(3))

      NODN=NODEID(0)

      !!-------------------READING FNPT--------------------!!
      ! FN_VA(:,1:6) and MP_VA(:,1:7) description
      ! FN_VA(:,1)    =   MP_VA(:,1)    =   pressure
      ! FN_VA(:,2)    =   MP_VA(:,2)    =   Ux
      ! FN_VA(:,3)    =   MP_VA(:,3)    =   Uy  = 0d0
      ! FN_VA(:,4)    =   MP_VA(:,4)    =   Uz
      ! FN_VA(:,5)    =   MP_VA(:,5)    =   Ax
      ! FN_VA(:,6)    =   MP_VA(:,6)    =   Az
      !                   MP_VA(:,7)    =   weight

      CALL FP%READFNPT( H0 )
      !!-----------------END READING FNPT------------------!!

      !!---------------------REMESHING---------------------!!
      IF(REMESHFREQ.EQ.0) GOTO 501
      IF((ISTEP.GT.0).AND.(MOD(ISTEP,REMESHFREQ).EQ.0))THEN
         WRITE(8,'(" [REM]----REMESHING AT TIMESTEP  = ",I10,"----")')ISTEP

         DOMBL(1)=MINVAL(COORX(:,1))
         DOMBL(2)=MINVAL(COORY(:,1))
         DOMBL(3)=MINVAL(COORZ(:,1))
         DOMTR(1)=MAXVAL(COORX(:,1))
         DOMTR(2)=MAXVAL(COORY(:,1))
         DOMTR(3)=MAXVAL(COORZ(:,1))
         NPOI1=NODEID(0)

         CALL GEN_FNPT_3D(LNODE, NPOI1, NODEID(-7:NPOI1), COORX(:,1), &
            FP%NX, FP%NY, FP%NN, FP%FNCO, FP%FNFS, FP%FNVA, FP%DDL,  &
            DOMBL, DOMTR, NFS, XFS, YFS, ZFS, CM)

         DO I=1,NODEID(-2)
            IF(NODEID(I).EQ.4)THEN
               NFS=NFS+1
               XFS(NFS)=COORX(I,1)
               YFS(NFS)=COORY(I,1)
               ZFS(NFS)=COORZ(I,1)
            ENDIF
         ENDDO

         DO I=1,BNDFS(0)
            NFS=NFS+1
            XFS(NFS)=COORX(BNDFS(I),1)
            YFS(NFS)=COORY(BNDFS(I),1)
            ZFS(NFS)=COORZ(BNDFS(I),1)
         ENDDO

         DOMBL(1)=MINVAL(XFS)
         DOMBL(2)=MINVAL(YFS)
         DOMBL(3)=0D0
         DOMTR(1)=MAXVAL(XFS)
         DOMTR(2)=MAXVAL(YFS)
         DOMTR(3)=0D0

         CALL FSDOM%FILLCELL(NFS,XFS(1:NFS),YFS(1:NFS),ZFSTMP(1:NFS), &
            RCELL*DDL,DOMBL,DOMTR)

         NPOI1=MM%NP
         !NLMAXN=15
         CALL MLPG_GET_ETA(FSDOM, NFS, XFS(1:NFS), YFS(1:NFS), &
            ZFS(1:NFS), NPOI1, MM%CX(1:NPOI1), MM%CY(1:NPOI1), &
            MM%ET(1:NPOI1), ERRTMP(1:NPOI1), DDL, 30, 1)

         I=MINLOC(MM%ET(1:NPOI1),1)
         WRITE(8,'("     [---] MIN ETA, NID ",F15.6,I15)') &
            MM%ET(I),NODEID(MMTOIM(I))
         WRITE(8,'("     [---] MIN ETA  LOC ",2F15.6)') &
            MM%CX(I),MM%CY(I)

         DO I = 1, MM%NP
            ! MM%Z0 IS LOCAL DEPTH
            IX = MMTOIM(I)
            TMPR1 = MM%Z0(I)
            TMPR2 = IM%CZ(IX)
            TMPR3 = MM%ET(I)
            MM%CZ(I) = TMPR1 + (TMPR3-TMPR1)/(H0-TMPR1)*(TMPR2-TMPR1)
         ENDDO

         I=CM%NP+NODEID(-1)
         IF(I.GT.CM%NPS)THEN
            WRITE(8,'(" [ERR] ERR IN SIZE OF ARRAYS IN CM OBJECT P2")')
            STOP
         ENDIF
         DO I=1,NODEID(-1)
            IF(NWALLID(I,2).EQ.-10)CYCLE

            CM%NP=CM%NP+1
            CM%CX(CM%NP)=COORX(I,1)
            CM%CY(CM%NP)=COORY(I,1)
            CM%CZ(CM%NP)=COORZ(I,1)
            CM%PR(CM%NP)=P(I)
            CM%UX(CM%NP)=UX(I,1)
            CM%UY(CM%NP)=UY(I,1)
            CM%UZ(CM%NP)=UZ(I,1)
         ENDDO

         I=CM%NP
         DOMBL(1)=MINVAL(CM%CX(1:I))
         DOMBL(2)=MINVAL(CM%CY(1:I))
         DOMBL(3)=MINVAL(CM%CZ(1:I))
         DOMTR(1)=MAXVAL(CM%CX(1:I))
         DOMTR(2)=MAXVAL(CM%CY(1:I))
         DOMTR(3)=MAXVAL(CM%CZ(1:I))

         CALL MLDOM%FILLCELL(I,CM%CX(1:I),CM%CY(1:I),CM%CZ(1:I), &
            RCELL*DDL,DOMBL,DOMTR)

         I=CM%NP
         J=MM%NP

         !NLMAXN=61
         CALL MLPG_GET_UP(MLDOM,I,CM%CX(1:I),CM%CY(1:I),CM%CZ(1:I), &
            CM%UX(1:I),CM%UY(1:I),CM%UZ(1:I),CM%PR(1:I), &
            J,MM%CX(1:J),MM%CY(1:J),MM%CZ(1:J), &
            MM%UX(1:J),MM%UY(1:J),MM%UZ(1:J),MM%PR(1:J),DDL,62)
         
         DO I = 1, MM%NP
            IX = MMTOIM(I)
            COORX(IX,1)=MM%CX(I)
            COORY(IX,1)=MM%CY(I)
            COORZ(IX,1)=MM%CZ(I)
            UX(IX,1)=MM%UX(I)
            UY(IX,1)=MM%UY(I)
            UZ(IX,1)=MM%UZ(I)
            CSUXT1(IX,1)=MM%UX(I)
            CSUYT1(IX,1)=MM%UY(I)
            CSUZT1(IX,1)=MM%UZ(I)
            P(IX)=MM%PR(I)
         ENDDO

         CALL CM%ENDDOMAINDATA

         CALL U_BOUNDARY2(IFSI)

         ! PENDING
         ! CALL FOR BOTSLIPBC
         ! DONE AFTER JUDGEBOTTOM

         WRITE(8,'(" [REM]----REMESHING COMPLETED AT = ",I10,"----")') &
            ISTEP
         WRITE(8,*)

      ENDIF
501   CONTINUE
      !!-------------------END REMESHING-------------------!!

      DOMBL(1)=MINVAL(COORX(:,1))
      DOMBL(2)=MINVAL(COORY(:,1))
      DOMBL(3)=MINVAL(COORZ(:,1))
      DOMTR(1)=MAXVAL(COORX(:,1))
      DOMTR(2)=MAXVAL(COORY(:,1))
      DOMTR(3)=MAXVAL(COORZ(:,1))

      CALL SYSTEM_CLOCK(CPUT%TI(5))
      NPOI1=NODEID(0)

      CALL MLDOM%FILLCELL(NPOI1, COORX(1:NPOI1,1), &
         COORY(1:NPOI1,1), COORZ(1:NPOI1,1), RCELL*DDL, DOMBL, DOMTR)

      CALL NODELINK_3_SHA(MLDOM,LNODE,NPOI1,SCALE,DDL,DDR, &
         NODEID(-2:NPOI1),NWALLID,COORX(:,1),COORY(:,1),COORZ(:,1))

      CALL JUDGEBOTTOM(LNODE, NPOI1, NODEID(-7:NPOI1), NWALLID, &
         COORX(1:NPOI1,1), COORY(1:NPOI1,1), COORZ(1:NPOI1,1), &
         DDR(1:NPOI1), FSSRCH(1:NPOI1))

      !!---------------------REMESHING---------------------!!
      IF(REMESHFREQ.NE.0)THEN
         IF((ISTEP.GT.0) .AND. (MOD(ISTEP,REMESHFREQ).EQ.0))THEN

            ! THIS WAS NEEDED AFTER U_BOUNDARY2() IN THE REMESH CODE
            ! THAT HAD SET TYPE2 VEL TO 0
            ! THOUGH APPLICATION OF THIS NEEDED JUDGEBOTTOM
            ! WHICH INSTEAD ALSO NEEDED NEW MLDOM%FILLCELL
            ! AS THE OLD CALL WAS USED BY CM%
            CALL BOTSLIPBC(LNODE, NPOI1, NODEID(-7:NPOI1), NWALLID, &
               COORX(1:NPOI1,1), COORY(1:NPOI1,1), COORZ(1:NPOI1,1), &
               DDR(1:NPOI1), SNX, SNY, SNZ, SMX, SMY, SMZ, SSX, SSY, SSZ, &
               UX(1:NPOI1,1), UY(1:NPOI1,1), UZ(1:NPOI1,1))
         ENDIF
      ENDIF
      !!-------------------END REMESHING-------------------!!

      IF(ISTEP.EQ.ISTART)THEN
         CALL CHECK_NUMNE
         CALL GENERATEGHOST(DDR)

         CALL MLDOM%FILLCELL(NODEID(0),COORX(:,1),COORY(:,1), &
            COORZ(:,1),RCELL*DDL,DOMBL,DOMTR)

         CALL NODELINK_3_SHA(MLDOM,LNODE,NPOI1,SCALE,DDL,DDR, &
            NODEID(-2:NPOI1),NWALLID,COORX(:,1),COORY(:,1),COORZ(:,1))
         
         CALL JUDGEBOTTOM(LNODE, NPOI1, NODEID(-7:NPOI1), NWALLID, &
            COORX(1:NPOI1,1), COORY(1:NPOI1,1), COORZ(1:NPOI1,1), &
            DDR(1:NPOI1), FSSRCH(1:NPOI1))

         NPOI1=NODEID(0)
         CALL SRI2(LNODE,NPOI1,NODEID(-7:NPOI1),NWALLID, &
            COORX(:,2),COORY(:,2),COORZ(:,2), &
            DDR,FSSRCH,H0)
      ENDIF

      CALL SYSTEM_CLOCK(CPUT%TI(6))
      WRITE(8,'(" [TIM] TIME OF NODELINK_3_SHA IS ",F15.6)'), &
         1D0*(CPUT%TI(6)-CPUT%TI(5))/CPUT%SYSRATE
      WRITE(8,*)

      NODN=NODEID(-1)

      IF(ISTEP.NE.0)TOTAL_TIME=TOTAL_TIME+DT
      WRITE(8,*)'THE TIME STEP IS',DT
      WRITE(8,*)'TOTAL SIMULATING TIME IS',TOTAL_TIME

      !!-----------------INTERPOLATE FNPT------------------!!
      WRITE(8,*)'[MSG] FNPT INTERPOLATION START'

      CALL FNPT_INTERP( FP%NN, FP%FNCO, FP%FNVA, FP%MLNP,  &
         FP%MLPOI, FP%MLVAL, DDR)

      WRITE(8,*)'[MSG] FNPT INTERPOLATION STOP'
      !!---------------END INTERPOLATE FNPT----------------!!

      DU=0.D0
      DV=0.D0
      DW=0.D0

      DO I=1,NODEID(-2)

         !CALCULATE THE VISCOUS TERM
         !IF(I_CAL_V.EQ.1) CALL LAPLACIAN(I,DU,DV,DW)

         UX(I,2)=UX(I,1)+(VCOEFF*DU)*DT
         UY(I,2)=UY(I,1)+(VCOEFF*DV)*DT
         UZ(I,2)=UZ(I,1)+(VCOEFF*DW)*DT
         COORX(I,2)=COORX(I,1) !+0.5*(UX(I,2)+UX(I,1))*DT   !UX(I,2)*DT   !X COORDINATE VELOCITY
         COORY(I,2)=COORY(I,1) !+0.5*(UY(I,2)+UY(I,1))*DT   !UY(I,2)*DT   !Y COORDINATE VELOCITY
         COORZ(I,2)=COORZ(I,1) !+0.5*(UZ(I,2)+UZ(I,1))*DT   !UZ(I,2)*DT   !Y COORDINATE VELOCITY
         UN(I,1)=VCOEFF*DU*DT               !VISCOUS FORCE IN X, Y,Z DIRECTIONS
         UN(I,2)=VCOEFF*DV*DT
         UN(I,3)=VCOEFF*DW*DT
      ENDDO

      DO I = NODEID(-2)+1,NODEID(0)
         COORX(I,2) = COORX(I,1)
         COORY(I,2) = COORY(I,1)
         COORZ(I,2) = COORZ(I,1)
         UX(I,2)=UX(I,1) !+(VCOEFF*DU)*DT
         UY(I,2)=UY(I,1) !+(VCOEFF*DV)*DT
         UZ(I,2)=UZ(I,1) !+(VCOEFF*DW)*DT
      END DO

      NPOI1=NODEID(-1)
      CALL FIND_LAPLAC(NPOI1, UN(1:NPOI1,1:3))

      !! FNPT Coupling
      IF(I_WM.EQ.15)THEN
         WRITE(8,*)'[MSG] WAVEMAKER IS TYPE 15'

         DO NNI=NODEID(-2)+1,NODEID(-3)
            NNI2=NNI-NODEID(-2)
            IF(NODEID(NNI).NE.8)THEN
               WRITE(8,*)'[ERR] NOT A WAVEMAKER NODE',NNI
               STOP
            ENDIF

            COORX(NNI,2) = COORX(NNI,1)!+MP_VAT1(NNI2,2)*DT
            UX(NNI,1:2) = FP%MLVAL(NNI2,2)
            TANK_A(NNI,1) = FP%MLVAL(NNI2,5)

            COORZ(NNI,2) = COORZ(NNI,1)!+MP_VAT1(NNI2,4)*DT
            UZ(NNI,1:2) = FP%MLVAL(NNI2,4)
            TANK_A(NNI,3) = FP%MLVAL(NNI2,6)

            PRESS_DR(NNI2) = FP%MLVAL(NNI2,1) + (1000D0*9.81D0*(COORZ(NNI,2)-H0))
         ENDDO
      ENDIF

      ! SOLVE THE PARTICLE NUMBER DENSITY AND LAMDA VALUE
      !CALL PND_LAMDA1 !commented on 2021-04-07

      CALL GENERATEGHOST(DDR)

      CALL JUDGEFREESURFACE_SHA(DDR,FSSRCH,SPONGEX)

      WRITE(1616,'(" TIME = ",F15.6)')TOTAL_TIME
      !CALL SRI2(DDR,FSSRCH,H0)
      NPOI1=NODEID(0)
      CALL SRI2(LNODE,NPOI1,NODEID(-7:NPOI1),NWALLID, &
         COORX(:,2),COORY(:,2),COORZ(:,2), &
         DDR,FSSRCH,H0)

      NPOI1=NODEID(-1)
      CALL FIND_ACCN(NPOI1, UM(1:NPOI1))

      !!---------------------------------------------------!!

      !!---------------------------------------------------!!
      CALL SYSTEM_CLOCK(CPUT%TI(5))
      !BC3+MLPG (INNER PARTICLE FORMULATIONS)
      !CALL ASSEMATRIX_MLPG(FB,UM,ISTEP,WMA,PRESS_DR,NTHR)

      NPOI1=NODEID(0)
      !NLMAXN=200
      CALL ASSEMATRIX_MLPG_SHA(LNODE,NPOI1,NODEID(-2:NPOI1), &
         NWALLID,FB,PRESS_DR,GRA,H0,KW,COORX(:,2),COORY(:,2), &
         COORZ(:,2),SNX,SNY,SNZ,200,MBAS,I_WM)

      CALL FILL_MATRIX_SHA(NTHR,LNODE,NPOI1,NODEID(-2:NPOI1), &
         NWALLID,COORX(:,2),COORY(:,2), &
         COORZ(:,2),50,MBAS,KW,R1,DDR,DT,FB)

      CALL SYSTEM_CLOCK(CPUT%TI(6))
      WRITE(8,'(" [TIM] TOTAL TIME FILLING MATRIX IS ",F15.6)'), &
         1D0*(CPUT%TI(6)-CPUT%TI(5))/CPUT%SYSRATE

      CALL SYSTEM_CLOCK(CPUT%TI(5))
      DO I=1,NODEID(0)
         PTMP(I)=P(I)+((COORZ(I,2)-h0)*rou(i)*(-GRA))
      ENDDO

      EPS_G=1.E-10   !THE SECOND ERROR FOR GMRESWHOLEMA
      ERRSOL = 0
      NPOI1=NODEID(-1)

      CALL PRESSURE_SOLVER2(NPOI1,PTMP,FB,EPS_G,ERRSOL)
      CALL SYSTEM_CLOCK(CPUT%TI(6))
      WRITE(8,'(" [TIM] TIME OF PRESSURE_SOLVER IS ",F15.6)'), &
         1D0*(CPUT%TI(6)-CPUT%TI(5))/CPUT%SYSRATE
      WRITE(8,*)
  
      !CHANGE IT INTO REAL PRESSURE P=P-(Z-D)
      DO I=1,NODEID(0) 
         PTMP(I)=PTMP(I)-((COORZ(I,2)-h0)*rou(i)*(-GRA))
      ENDDO

      IF(ERRSOL.EQ.1) GOTO 201
      !x---------------------------------------------------x!

      !!---------------------------------------------------!!

      P = 0.d0

      NPOI1=NODEID(0)
      CALL GHOSTPART(LNODE,NPOI1,NODEID(-7:NPOI1),NWALLID, &
        COORX(:,2),COORY(:,2),COORZ(:,2),NLMAXN,MBAS,KW,R1,DDR,PTMP, &
        MIRRNP,MIRRXY)

      CALL PRESSURE_SMOOTH_SHA(LNODE,NPOI1,NODEID(-2:NPOI1), &
         NWALLID,COORX(:,2),COORY(:,2),COORZ(:,2),NLMAXN,MBAS,KW,R1, &
         DDR,PTMP,P)
      201 CONTINUE
   ENDDO

END PROGRAM THREED_BREAKINGWAVE
!!----------------------- END OF MAIN PROGRAM ---------------------!!

!!------------------------------INPUT------------------------------!!
SUBROUTINE INPUT(LNODETMP, H0, DDL, SCALE, KW, MBAS, DT,  &
   TOTAL_TIME, IPRINT, I_PF, I_PF1, &
   NSTEPS, I_CAL_V, VCOEFF, I_WM, IFSI, NTHR, MESHFILE,  &
   DOMX, DOMY, DOMZ, CYLX, CYLY, CYLR, SPONGEX, REMESHFREQ, &
   RESFREQ, RESUMECHK, RESUMEFILE, WP, PP, FP)
   USE PROBESMOD
   USE FNPTCPLMOD
   IMPLICIT NONE

   INTEGER(KIND=4),INTENT(OUT)::LNODETMP, KW, MBAS
   INTEGER(KIND=4),INTENT(OUT)::IPRINT, I_PF, I_PF1, NTHR
   INTEGER(KIND=4),INTENT(OUT)::NSTEPS, I_CAL_V, I_WM, IFSI
   INTEGER(KIND=4),INTENT(OUT)::REMESHFREQ, RESFREQ
   REAL(KIND=8),INTENT(OUT)::H0, DDL, SCALE, DT, TOTAL_TIME
   REAL(KIND=8),INTENT(OUT)::VCOEFF, SPONGEX
   REAL(KIND=8),INTENT(OUT)::DOMX(2), DOMY(2), DOMZ(2)
   REAL(KIND=8),INTENT(OUT)::CYLX, CYLY, CYLR
   CHARACTER(LEN=256),INTENT(OUT)::MESHFILE, RESUMEFILE
   TYPE(PROBETYP),INTENT(INOUT)::WP, PP
   TYPE(FNPTCPLTYP),INTENT(INOUT)::FP
   LOGICAL,INTENT(OUT)::RESUMECHK

   INTEGER(KIND=4)::I,NP
   CHARACTER(LEN=256)::TEXT1

   OPEN (UNIT=15, FILE='mlpgrInput.txt', STATUS='OLD')

   ! WATER DEPTH
   READ(15,*) TEXT1
   READ(15,*) H0

   ! READ DISTANCE BETWEEN NODES
   READ(15,*) TEXT1
   READ(15,*) DDL

   ! READ  SCALE TO DETERMINE THE DOMAIN OF INFLUENCE
   READ(15,*) TEXT1
   READ(15,*) SCALE

   ! READ  COEFFICIENT FOR GAUSS WEIGHT FUNCTION
   READ(15,*) TEXT1
   READ(15,*) KW

   ! READ  ORDER OF THE BASE FUNCTION
   READ(15,*) TEXT1
   READ(15,*) MBAS

   ! READ  THE TOTAL NUMBER OF TIME STEP AND TIME STEP
   READ(15,*) TEXT1
   READ(15,*) DT,TOTAL_TIME, NSTEPS

   ! READ RECORD DATA TIME,THE FREQUENCY OF PRINT OUT IN THE NUMBER OF TYPE STEPS
   ! RECORD THE TIME STEP, BEFORE THAT TIME, THE RECORD FREQUENCY IS EVRSTEP,AFTER
   ! THAT TIME, THE RECORD FREQUENCY IS EVRSTEP1
   READ(15,*) TEXT1
   READ(15,*) IPRINT,I_PF,I_PF1


   READ(15,*) TEXT1
   READ(15,*) RESFREQ
   READ(15,*) TEXT1
   READ(15,*) RESUMECHK
   READ(15,*) TEXT1
   READ(15,*) RESUMEFILE

   ! READ IF CALCULATE THE VICOUS, AND THE VISCOUS COEFFICIENT
   READ(15,*) TEXT1
   READ(15,*) I_CAL_V,VCOEFF

   ! THE KIND OF WAVEMAKER,AND FREQUENCY OF SMOOTH VELOCITY,OF SLOSH
   READ(15,*) TEXT1
   READ(15,*) I_WM

   ! IF ELASTIC STRUCTRE
   READ(15,*) TEXT1
   READ(15,*) IFSI

   READ(15,*) TEXT1
   READ(15,*) NTHR

   ! LNODE
   READ(15,*) TEXT1
   READ(15,*) LNODETMP

   READ(15,*) TEXT1
   READ(15,*) MESHFILE
   READ(15,*) TEXT1
   READ(15,*) DOMX(1), DOMY(1), DOMZ(1)
   READ(15,*) TEXT1
   READ(15,*) DOMX(2), DOMY(2), DOMZ(2)
   READ(15,*) TEXT1
   READ(15,*) CYLX, CYLY, CYLR

   READ(15,*) TEXT1
   READ(15,*) SPONGEX
   READ(15,*) TEXT1
   READ(15,*) REMESHFREQ

   !! WAVE PROBES
   READ(15,*) TEXT1
   READ(15,*) NP
   CALL WP%INITPROBE(NP,1,602)
   DO I=1,NP
      READ(15,*) WP%XYZ(I,1), WP%XYZ(I,2)
      WP%XYZ(I,3)=0D0
   ENDDO

   !! PRESSURE PROBES
   READ(15,*) TEXT1
   READ(15,*) NP
   CALL PP%INITPROBE(NP,1,603)

   DO I=1,NP
      READ(15,*) PP%XYZ(I,1), PP%XYZ(I,2), PP%XYZ(I,3)
   ENDDO

   !! FNPT Coupling Inputs
   READ(15,*) TEXT1
   READ(15,*) TEXT1
   READ(15,*) FP%FILENAME
   READ(15,*) TEXT1
   READ(15,*) FP%NX, FP%NY
   READ(15,*) TEXT1
   READ(15,*) FP%DDL
   READ(15,*) TEXT1
   READ(15,*) FP%X0
   READ(15,*) TEXT1
   READ(15,*) FP%RLXLEN
   FP%FILE=601
   FP%NN = FP%NX * FP%NY

   ! OUTPUT THE INPUT DATA TO CHECK IF CORRECT
   WRITE(8,'(1F12.6)') H0
   WRITE(8,'(1F12.6)') DDL
   WRITE(8,'(1F12.6)') SCALE
   WRITE(8,'(1I12)') KW
   WRITE(8,'(1I12)') MBAS
   WRITE(8,'(2F12.4,1I12)') DT,TOTAL_TIME, NSTEPS
   WRITE(8,'(3I12)') IPRINT,I_PF,I_PF1
   WRITE(8,'(1I12,1F12.6)') I_CAL_V,VCOEFF
   WRITE(8,'(2I12)') I_WM
   WRITE(8,'(I12)') IFSI
   WRITE(8,'(I12)') NTHR
   WRITE(8,'(A)') TRIM(MESHFILE)
   WRITE(8,'(3F12.6)') DOMX(1), DOMY(1), DOMZ(1)
   WRITE(8,'(3F12.6)') DOMX(2), DOMY(2), DOMZ(2)
   WRITE(8,'(3F12.6)') CYLX, CYLY, CYLR
   WRITE(8,'(F12.6)') SPONGEX
   WRITE(8,*)

   CLOSE(15)

   RETURN
END SUBROUTINE INPUT
!!----------------------------END INPUT----------------------------!!


!!--------------------------SETFLUXPLANE---------------------------!!
SUBROUTINE MLPGFLUXPLANE(DOMX, DOMY, DOMZ, MLDDL, &
   FLUXIN, FLUXOT, VOLPL)
   USE FLUXPLANEMOD
   USE VOLPLANEMOD
   IMPLICIT NONE

   TYPE(FLUXPLANE),INTENT(INOUT)::FLUXIN, FLUXOT
   TYPE(VOLPLANE),INTENT(INOUT)::VOLPL
   REAL(KIND=8),INTENT(IN)::DOMX(2), DOMY(2), DOMZ(2), MLDDL

   INTEGER(KIND=4)::I, J
   REAL(KIND=8)::TMPR1, TMPR2, BL(3), TR(3)


   TMPR1 = MLDDL    !! FLUXPLANE XY DISTANCE B/W NODES
   TMPR2 = MLDDL/2D0  !! FLUXPLANE Z DISTANCE B/W NODES
   I = FLOOR((DOMY(2) - DOMY(1))/TMPR1)+1  !NXY
   J = FLOOR((DOMZ(2) - DOMZ(1))/TMPR2)+1  !NZ


   !! FLUXIN
   TMPR1 = 3D0  !X POSITION OF THE PLANE
   BL(1) = TMPR1
   BL(2) = DOMY(1)
   BL(3) = DOMZ(1)
   TR(1) = TMPR1
   TR(2) = DOMY(2)
   TR(3) = DOMZ(2)
   CALL FLUXIN%SETFLUXPLANE(I, J, BL, TR)

   !! FLUXOT
   TMPR1 = 10D0  !X POSITION OF THE PLANE
   BL(1) = TMPR1
   BL(2) = DOMY(2)
   BL(3) = DOMZ(1)
   TR(1) = TMPR1
   TR(2) = DOMY(1)
   TR(3) = DOMZ(2)
   CALL FLUXOT%SETFLUXPLANE(I, J, BL, TR)

   !! VOLPL
   BL(1) = 3D0
   BL(2) = DOMY(2)
   TR(1) = 10D0
   TR(2) = DOMY(1)
   I = FLOOR( ABS(TR(1) - BL(1)) / MLDDL )+1
   J = FLOOR( ABS(TR(2) - BL(2)) / MLDDL )+1
   CALL VOLPL%SETVOLPLANE(I, J, BL(1:2), TR(1:2))


END SUBROUTINE MLPGFLUXPLANE
!!------------------------END SETFLUXPLANE-------------------------!!