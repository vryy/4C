c
c     Maintainer: Malte Neumann 
c                 neumann@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/neumann/ 
c                 0711 - 685-6121 
c
c     ---------------------------------------------------------------  
C=========================================================================
      SUBROUTINE SSPACE (A,ACOP,B,MAXA,R,EIGV,TT,W,
     1                   AR,BR,VEC,D,RTOLV,IND,NN,NNM,NWK,NWM,NC,
     2                   NSTA,ISOL,NNZ,EIGFOU,RFOU,NROOT,NITEM,NSMAX,
     3                   IFSS,ISUB,TOLEIG,TOLJAC,SHIFT,ISTLDL,INIT,
     4                   IPRINT)
C-------------------------------------------------------------------------
C     TOPIC : CONTROL PROGRAM FOR EIGENVALUE ANALYSIS 
C             TO SOLVE THE GENERALIZED PROBLEM (A-LAMDA*B)*V = 0
C             USING SUBSPACE ITERATION METHOD
C             RESULTS : EIGENVALUES AND EIGENVECTORS
C-------------------------------------------------------------------------
C     PROGRAMMED BY  HANS STEGMUELLER  , OCT 1988
C     ORIGINAL PROGRAM FROM BATHE/WILSON,
C     "NUMERICAL METHODS IN FINITE ELEMENT ANALYSIS", PAGE 509, (1976)
C------------------------------------------------------------------------
C     A      .... SYSTEM MATRIX TO BE MODIFIED IN SSPACE             (W) 
C     ACOP   .... SYSTEM MATRIX (ORIGINAL)                           (I)
C     B      .... MASS MATRIX  -OR-  GEOMETRIC STIFFNESS MATRIX      (I)
C     MAXA   .... POINTER VECTOR TO THE DIAGONAL ELEMENTS            (I)
C     R      .... EIGENVECTORS (ITERATIONVECTORS)                    (O)
C     EIGV   .... EIGENVALUES                                        (0)
C     TT     .... WORKING ARRAY                                      (W)
C     W      .... WORKING ARRAY                                      (W)
C     AR     .... SYSTEM MATRIX A REDUCED PROBLEM                    (W)
C     BR     .... SYSTEM MATRIX B REDUCED PROBLEM                    (W)
C     VEC    .... EIGENVECTORS OF REDUCED PROBLEM                    (W)
C     D      .... WORKING SPACE                                      (W)
C     RTOLV  .... RELATIVE TOLERANCES                                (O)
C     NN     .... REAL PROBLEM SIZE (NUMBER OF UNKNOWNS)             (I)
C     NNM    .... LENGTH OF POINTER ARRAY MAXA                       (I)
C     NWK    .... SIZE OF MATRIX A                                   (I)
C     NWM    .... SIZE OF MATRIX B                                   (I)
C     NC     .... SIZE OF REDUCED PROBLEM (NUM.OF.ITERATIONVECTORS)  (I)
C     NSTA   .... RESTART FLAG                                       (I)
C     NNZ    .... NUMBER OF NONZERO ELEMENTS IN SYSTEM MATRIX        (I)
C------------------------------------------------------------------------
C                 (W) ... PARAMETER ARE ONLY USED IN THIS PROGRAM
C                 (O) ... RESULTS OF THIS PROGRAM UNIT 
C------------------------------------------------------------------------
      IMPLICIT  REAL*8 (A-H,O-Z)
      INTEGER INFO
C       PARAMETER (ZERO=0.0D0,RIN=-0.99D+99)
      CHARACTER IFCTR*8,LINE*120
C       CHARACTER SNAME*6
      CHARACTER*4 OP
C      PARAMETER (SNAME='SSPACE')
      DIMENSION A(NWK),B(NWM),R(NN,NC),TT(NN),W(NN),EIGV(NC)
      DIMENSION D(NC),VEC(NC,NC),AR(NC,NC),BR(NC,NC),RTOLV(NC)
      DIMENSION ACOP(NWK),MAXA(NNM),IND(NC,3)
C      
      DIMENSION EIGFOU(NC),RFOU(NN,NC)
C
C-----------------------------------------------------------------------
      ZERO=0.0D0
      RIN=-0.99D+99
C--------------------------------------------- OPEN OUPUTFILE
      INFO=6
      IF (IPRINT.GT.1) THEN
      OPEN(INFO, FILE='eigen.doc')
      ENDIF
C--------------------------------------------- CHECK RANGE OF SEARCH     
C-----------------------------------------------------------------------
       IRANGE = 0
C-----------------------------------------------------------------------
      IF (IPRINT.GT.0) THEN
      IF (INIT.EQ.0) THEN
       LINE = 'E I G E N V A L U E   A N A L Y S I S'
       WRITE(INFO,*)LINE
       IF (ISUB.EQ.1) THEN
       WRITE(INFO,*) 'SUBSPACE ITERATION METHOD - JACOBI ROTATION '//
     +          'TO SOLVE SUBSPACE'
       ELSE IF (ISUB.EQ.2) THEN
         WRITE(INFO,*) 'SUBSPACE ITERATION METHOD - QZ-ALGORITHM '//
     +          'TO SOLVE SUBSPACE'
       ELSE
         WRITE (INFO,*) '... ILLEGAL PROGRAM PARAMETER ISUB'
         STOP
       ENDIF
       WRITE (INFO,2000) 'NUMBER OF REQUESTED EIGENVALUES : ',NROOT
       WRITE (INFO,2000) 'NUMBER OF ITERATION VECTORS     : ',NC   
       WRITE (INFO,2000) 'MAX. NUMBER OF ITERATIONS       : ',NITEM
       WRITE (INFO,2010) 'TOLERANCE FOR CONVERGENCE CHECK : ',TOLEIG
       WRITE (INFO,2010) 'RELATIVE SHIFT USED             : ',SHIFT
       IF (IRANGE .GT. 0) THEN
         WRITE (INFO,*)  'SEARCH IN A RANGE                         '
	 WRITE (INFO,2010) 'LOWER BOUND OF THE INTERVAL     : ',BOULO
	 WRITE (INFO,2010) 'UPPER BOUND OF THE INTERVAL     : ',BOUUP
       ENDIF
 2000  FORMAT (A,I12)
 2010  FORMAT (A,1P,E12.4)
      ENDIF
      ENDIF
C-----------------------------------------------------------------------
C                   SAVE INFORMATION AND SET NEW PARAMETER
C-----------------------------------------------------------------------
      ISH = 0    
C-----------------------------------------------------------------------
C                   SET STARTING ITERATION VECTORS
C-----------------------------------------------------------------------
      IF (NSTA.EQ.0) CALL SSPSTA (A,B,MAXA,R,W,TT,NN,NC,
     *                            NWK,NWM,NNM,ISOL,NNZ)
C-----------------------------------------------------------------------
C                   SET STARTING ITERATION VECTORS FOR A NEW SHIFT
C-----------------------------------------------------------------------
  555 CONTINUE
C-----------------------------------------------------------------------
C                   SHIFT EIGENSYSTEM - SET SHIFT FLAG
C-----------------------------------------------------------------------
      IF (SHIFT.NE.ZERO) THEN
        IF (ISOL.NE.4) THEN
          CALL MSHIFT (A,B,MAXA,NWK,NWM,NNM,NN,SHIFT)
        ENDIF
        IFSH=1
      ELSE
        IFSH=0
      ENDIF
C----------------------------------------------------------------------
C                   FACTORIZE MATRIX A INTO (L)*(D)*(L(T))
C----------------------------------------------------------------------
       ldl = 1
      IF (LDL.EQ.1) THEN
        ITIMSOL=0
        ILIN = 1
        NEQ = NN
       OP='TRI'
      CALL CALSOL (OP,ISOL,ITIMSOL,ILIN,NEQ,A,MAXA,NWK,NN,NNM,NSCH,TT)
      ENDIF
      NSCH0 = NSCH
      IF (NSCH0.GT.0) THEN
        WRITE (INFO,'(A,I3,A)') 'EIG001: SHIFTED PROBLEM HAS',NSCH0,
     +                          ' NEGATIVE DIAGONAL ELEMENTS'
      ENDIF
C-----------------------------------------------------------------------
C                   SAVE NSCH0 FOR FIRST STEP
C-----------------------------------------------------------------------
      IF (ISH.EQ.0 .AND. IRANGE .GT. 0) NSCHF = NSCH0
C      
C----------------------------------------STORE MATEIG-------------------
      IF (ISTLDL.GT.0) THEN
C        CALL VMCHAR ('EIGVAR:NAMLDL'  ,NAMLDL,'GET') ????????????
C        CALL AMDUP ('MATEIG',NAMLDL,IRC)
      ENDIF
C-----------------------------------------------------------------------
C      I T E R A T I O N   P R O G R A M                                 
C      - NROOT MAY BE INCREASED IF NEGATIVE EIGENVALUES ARE PRESENT
C      - ITERATION IS STOPPED IF * CONVERGENCE IS REACHED
C                                * ITERATION LIMIT IS REACHED
C                                * ALL EIGENVALUES REMAIN NEGATIVE
C                                  AFTER 4 ITERATIONS
C-----------------------------------------------------------------------
      CALL SSPITE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,IND,
     1             NN,NNM,NWK,NWM,NROOT,TOLEIG,TOLJAC,SHIFT,NC,NITEM,
     +             NSMAX,ISUB,IFSH,NERR,IFCTR,NITE,NSTA,
     +             NNEG,IST,NCONVR,IPRINT,ISOL,NNZ,INFO)
C-----------------------------------------------------------------------
C       >IST< CONTROLS FURTHER PROCESSING
C       IF EQ.0  ITERATION WAS SUCCESSFUL
C          EQ.2  ITERATION LIMIT REACHED
C          EQ.-1 ONLY NEGATIVE EIGENVALUES FOUND (STOP ANALYSIS)
C-----------------------------------------------------------------------
      IF (IST.EQ.99) GO TO 900
      IF (IST.EQ.2) THEN
      IF (IPRINT.GT.0) THEN
        LINE = 'EIG007: ITERATION LIMIT REACHED - '//
     +         'STURM CHECK WILL BE APPLIED ON CONVERGED EIGENVALUES'
        WRITE(INFO,*)LINE
      ENDIF
      ENDIF
C-----------------------------------------------------------------------
C       SHIFT CALCULATED EIGENVALUES BACK               
C       RESORT EIGENVALUES IN ASCENDING ORDER
C       DATA FOR EIGENVALUES, TOLERANCES, MODE SHAPES MUST BE MOVED !!
C----------------------------------------------------------------------- 
      IF (IFSH.GT.0) THEN
        CALL ESHIFT(EIGV,NC,1,SHIFT)  
        CALL MXSIND (EIGV, IND, NROOT,3,'A',IRC)
        CALL MXSDAT (EIGV, IND, NROOT, 1,1,IRC)
        CALL MXSDAT (RTOLV,IND, NROOT, 1,2,IRC)
        CALL MXSDAT (R    ,IND, NROOT,NN,2,IRC)
      ENDIF
C-----------------------------------------------------------------------
C       CHECK IF THE SOLUTION IS IN THE RANGE
C       AND FIND A NEW SHIFT-FACTOR IF NECESSARY
C-----------------------------------------------------------------------
      IF (ISH.NE.0) GOTO 555
C-----------------------------------------------------------------------
C       PUT SOLUTION ON EIGV AND R
C-----------------------------------------------------------------------      
C-----------------------------------------------------------------------
C       CALCULATE ERROR NORMS AND PRINT IF REQUESTED
C----------------------------------------------------------------------- 
      CALL SSPERR (ACOP,B,MAXA,R,EIGV,TT,W,D,NN,NWK,NWM,NNM,
     *             NC,NROOT,NNZ,ISOL)
C----------------------------------------------------------------------
C       I T E R A T I O N   PRINTOUT
C----------------------------------------------------------------------
C      CALL OSQPAG (NUMLIN)
C      IF (NUMLIN.LT.NROOT+6) CALL OSPAGE
C----------------------------------------------------------------------
      IF (IPRINT.GT.0) THEN
      WRITE (INFO,'(A,I2,A)') '>>>> FINAL RESULTS AFTER ',
     +                        NITE,' ITERATIONS'
      WRITE (INFO,'(A)') '  COUNT       EIGENVALUES     TOLERANCES'//
     +                                            '    ERROR-NORMS'
      DO 200 I=1,NROOT
        WRITE (INFO,'(I6,4X,1P,3E15.6)') I,EIGV(I),RTOLV(I),D(I)
  200 CONTINUE
      IF (IRANGE .GT. 0 .AND. NROOT.NE.NROOTS) THEN
        WRITE(INFO,'(A,I3,A)') ' THERE ARE ONLY ',NROOT,' EIGENVALUE'//
     *                       ' IN THE RANGE'
      ENDIF    
      ENDIF
C--------------------------------------------------------------------
C       PERFORM STURM SEQUENCE CHECK ON HIGHEST CONVERGED EIGENVALUE
C       USE ORIGINAL MATRIX A AND ADD SHIFUP*B
C---------------------------------------------------------------------
C       NCONVR ... NUMBER OF CONVERGED EIGENVALUES (EV)
C       NNEG   ... NUMBER OF NEGATIVE CONVERGED EIGENVALUES
C       NSCH0  ... NUMBER OF INITIAL NEGATIVE DIAGONAL ELEMENTS
C       NSCHP  ... NUMBER OF FINAL   NEGATIVE DIAGINAL ELEMENTS
C       CONDITION:
C                   NSCHP - NSCH0 = NCONVR - NNEG
C---------------------------------------------------------------------
C       COUNT FOR SHIFT FACTOR MUST INCLUDE CHECK ON TOLERANCE, BECAUSE
C       NOT ALL EIGENVALUES MAY BE CONVERGED (ITERATION LIMIT !)
C--------------------------------------------------------------------
      IF (NCONVR.GT.0) THEN
        IF (IRANGE .EQ. 0) THEN
	  SHIFUP = -1.0D20
          M=0
          DO 300 I=1,NROOT
            IF (RTOLV(I).LE.TOLEIG) THEN
              IF (EIGV(I).GT.SHIFUP) SHIFUP = EIGV(I)
              M=M+1
              IF (M.EQ.NCONVR) GO TO 310
            ENDIF
  300     CONTINUE
  310     CONTINUE
          SHIFUP = SHIFUP + ABS(SHIFUP)*0.01
	ELSE
	  SHIFUP = EIGV(NROOT) + ABS(EIGV(NROOT))*0.01
	ENDIF 
        CALL MXCR8  (ACOP,NWK,1,A)
        IF (ISOL.NE.4) THEN
          CALL MSHIFT (A,B,MAXA,NWK,NWM,NNM,NN,SHIFUP)
        ENDIF
        IF (LDL.EQ.1) THEN
        CALL CALSOL (OP,ISOL,ITIMSOL,ILIN,NEQ,A,MAXA,NWK,NN,NNM,NSCH,TT)
        ENDIF
        NSCHP = NSCH
C----------------------------------------------------------------------
        IF (IRANGE .EQ. 0) THEN
	  NMIS = (NSCHP - NSCH0) - (NCONVR - NNEG)
	ELSE
	  NMIS = (NSCHP - NSCHF) - NCONVR
	ENDIF
      IF (IPRINT.GT.0) THEN
        IF (NMIS.EQ.0) THEN
          LINE='CORRECT NUMBER OF EIGENVALUES VERIFIED BY '//
     +         'STURM SEQUENCE CHECK'
        WRITE(INFO,*)LINE
      ENDIF
        ELSE IF (NMIS.GT.0) THEN
      IF (IPRINT.GT.0) THEN
          WRITE (LINE,'(A,1PE12.5)') 'EIG002: VERIFY EIGENVALUES - '//
     +                               'SHIFT CHECK APPLIED AT ',SHIFUP
        WRITE(INFO,*)LINE
          WRITE (LINE,'(A,I2)') 'EIG002: NUMBER OF   M I S S I N G  '//
     +                          'EIGENVALUES = ',NMIS
        WRITE(INFO,*)LINE
      ENDIF
        ELSE IF (NMIS.LT.0) THEN
      IF (IPRINT.GT.0) THEN
          WRITE (LINE,'(A,1PE12.5)') 'EIG003: VERIFY EIGENVALUES - '//
     +                               'SHIFT CHECK APPLIED AT ',SHIFUP
        WRITE(INFO,*)LINE
          LINE = 'EIG003: M O R E  EIGENVALUES FOUND '//
     +           'THAN EXPECTED AFTER SHIFT'
        WRITE(INFO,*)LINE
          WRITE (LINE,'(A,I2)') '- EIGENVALUES FOUND        = ',NCONVR
        WRITE(INFO,*)LINE
          NUM = NSCHP - NSCH0 + NNEG
          WRITE (LINE,'(A,I2)') '- STURM CHECK EXPECTS ONLY = ',NUM
        WRITE(INFO,*)LINE
      ENDIF
        ENDIF
      ELSE
        LINE = 'EIG006: NO CONVEREGED POSITIVE EIGENVALUES FOUND '//
     +         ' - NO STURM SEQUENCE CHECK APPLIED'
        WRITE(INFO,*)LINE
      ENDIF
C------------------------------------------------------------------
 900  RETURN
      END
C=====================================================================
      SUBROUTINE SSPERR(A,B,MAXA,R,EIGV,TT,W,D,NN,NWK,NWM,NNM,
     *                  NC,NROOT,NNZ,ISOL)
C----------------------------------------------------------------------
C     TOPIC :   P R O G R A M                                                 
C              TO CALCULATE THE ERROR NORMS E(NN) OF                    
C   A(NN,NN) * R(NN,NROOT) -DIAG(EIGV(NROOT)) * B(NN,NN)*R(NN,NROOT)  
C                      = 0                                            
C            OUTPUT :                                                 
C       D(NROOT)  CONTAINS THE ERROR NORMS                            
C----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
C      CHARACTER SNAME*6
C      PARAMETER (SNAME='SSPERR',ZERO=0.0D0)
      DIMENSION A(NWK),B(NWM),R(NN,NC),EIGV(NC),TT(NN),W(NN),D(NC)
      INTEGER MAXA(NNM)
C----------------------------------------------------------------------
      ZERO=0.0D0
C----------------------------------------------------------------------
C
      DO 580 L=1,NROOT
             RT=EIGV(L)
             IF(ISOL.NE.4) THEN
               CALL XMULT (TT,A,R(1,L),MAXA,NN,NWK,NNM)
             ENDIF
             VNORM=0.
             DO 590 I=1,NN
                    VNORM=VNORM + TT(I)*TT(I)
 590         CONTINUE
             IF(ISOL.NE.4) THEN
               CALL XMULT (W,B,R(1,L),MAXA,NN,NWM,NNM)
             ENDIF
             WNORM=0.
             DO 600 I=1,NN
                    TT(I)=TT(I) - RT*W(I)
                    WNORM=WNORM + TT(I)*TT(I)
 600         CONTINUE
             VNORM=SQRT(VNORM)
             WNORM=SQRT(WNORM)
             IF (VNORM.EQ.ZERO) WRITE(*,*)'SSPERR',
     +    'ERROR NORM CANNOT BE CALCULATED BECAUSE OF VNORM=0' 
             D(L)=WNORM/VNORM
 580  CONTINUE
C
      RETURN
      END
C====================================================================
      SUBROUTINE SSPITE (A,B,MAXA,R,EIGV,TT,W,AR,BR,VEC,D,RTOLV,
     *                   IND,NN,NNM,NWK,NWM,NROOT,TOLEIG,TOLJAC,SHIFT,
     *                   NC,NITEM,NSMAX,ISUB,IFSH,NERR,IFCTR,NITE,NSTA,
     +                   NNEG,IST,NCONVR,IPRINT,ISOL,NNZ,INFO)
C----------------------------------------------------------------------
C      TOPIC : PERFORM THE SUBSPACE ITERATION LOOP
C              ISUB=1 :  GENERALIZED JACOBI ITERATION
C              ISUB=2 :  QZ - ALGORITHM
C-----------------------------------------------------------------------
C      PROGRAMMED BY HANS STEGMUELLER , OCT 1988  
C-------------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C      PARAMETER (ZERO=0.0D0)
      INTEGER INFO
      INTEGER MAXA (NNM),IND(NC,3)
      CHARACTER IFCTR*8,LINE*120
      CHARACTER*4 OP
      DIMENSION A(NWK),B(NWM),R(NN,NC),TT(NN),W(NN),EIGV(NC)
      DIMENSION D(NC),VEC(NC,NC),AR(NC,NC),BR(NC,NC),RTOLV(NC)
C-----------------------------------------------------------------------
      ZERO=0.0D0
C-----------------------------------------------------------------------
C       X ... PRINT EIGENVECTORS OF REDUCED PROBLEM
C       S ... PRINT SYSTEM MATRICES OF REDUCED PROBLEM
C       R ... PRINT RELATIVE TOLERANCES
C----------------------------------------------------------------------- 
      IFX=INDEX(IFCTR,'X')
      IFS=INDEX(IFCTR,'S')
      IFR=INDEX(IFCTR,'R')
C      IFE=IFX+IFS+IFR
      IFE=0
C-----------------------------------------------------------------------
      ICONV = 0
      NITE  = 0
      N1    = NC+1
C-----------------------------------------------------------------------
      DO 50 K=1,NC
        D(K) = ZERO
   50 CONTINUE
C--------------------------------------------------------------------
C          S T A R T   O F   I T E R A T I O N   L O O P
C--------------------------------------------------------------------
  100 NITE=NITE + 1
C
      IF (IPRINT.GT.0) THEN
      WRITE (INFO,'(A,I3,A)') '>>>> ITERATION  NO. : ',NITE,''
      IF (IFE.GT.0) WRITE (*,'(A,I2)') '>>> ITERATION NO. : ',NITE
      ENDIF
C
C---------------------------------------------------------------------
C         CALCULATE THE PROJECTION OF MATRIX A       
C---------------------------------------------------------------------
      DO 250 J=1,NC
        DO 200 K=1,NN
          TT(K)=R(K,J)
  200   CONTINUE
C        CALL CALSOL ('RED')
      OP='RED'
      CALL CALSOL (OP,ISOL,ITIMSOL,ILIN,NN,A,MAXA,NWK,NN,NNM,NSCH,TT)
        DO 220 I=J,NC
          ART=ZERO
          DO 210 K=1,NN
            ART=ART + R(K,I)*TT(K)
  210     CONTINUE
          AR(I,J)=ART
          AR(J,I)=ART
  220   CONTINUE
        DO 230 K=1,NN
          R(K,J)=TT(K)
  230   CONTINUE
  250 CONTINUE
C-------------------------------------------------------------------
C         CALCULATE THE PROJECTION OF MATRIX B
C-------------------------------------------------------------------
      DO 350 J=1,NC
        IF(ISOL.NE.4) THEN
          CALL XMULT (TT,B,R(1,J),MAXA,NN,NWM,NNM)
        ENDIF
        DO 320 I=J,NC
          BRT=ZERO
          DO 310 K=1,NN
            BRT=BRT + R(K,I)*TT(K)
  310     CONTINUE
          BR(I,J)=BRT
          BR(J,I)=BRT
  320   CONTINUE
        IF (ICONV.GT.0) GOTO 350
        DO 330 K=1,NN
          R(K,J)=TT(K)
  330   CONTINUE
  350 CONTINUE
C--------------------------------------------------------------------
      IF (IPRINT.GT.0) THEN
      IF (IFS.GT.0) THEN
        WRITE (INFO,'(A,I5)') 'ITERATION NO. :',NITE
        WRITE (INFO,'(A)') 'SYSTEM MATRICES OF REDUCED PROBLEM'
      ENDIF
      ENDIF
C---------------------------------------------------------------------
C     SOLVE FOR EIGENSYSTEM OF SUBSPACE OPERATORS
C--------------------------------------------------------------------
      IF (ISUB.NE.2) THEN
        CALL EIG_JACOBI (AR,BR,VEC,EIGV,W,NC,TOLJAC,SHIFT,IFSH,NSMAX)
      ELSE
        IERROR = 0
        CALL SSPEIG (AR,BR,VEC,EIGV,TT,W,IND,NC,IERROR)
      ENDIF
C---------------------------------------------------------------------
C     ARANGE EIGENVALUES WITH RESPECT TO ABSOLUTE VALUES
C     TO GET CONVERGENCE OF EIGENVALUES CLOSEST TO ZERO
C     MXSIND --> CREATE INDEX ARRAY (WITHOUT DATA MOVEMENT)
C     MXSDAT --> MOVE DATA INTO CORRECT POSITIONS
C---------------------------------------------------------------------
      CALL MXSIND (EIGV,IND,NC,3,'N',IRC)
      CALL MXSDAT (EIGV,IND,NC, 1,1,IRC)
      CALL MXSDAT (VEC ,IND,NC,NC,2,IRC)
C------------------------------------------------------------------
C     CALCULATE B TIMES APPROXIMATE EIGENVECTORS (ICONV.EQ.0)
C        -OR-   FINAL EIGENVECTOR APPROXIMATIONS (ICONV.GT.0)
C------------------------------------------------------------------
      DO 475 I=1,NN
        DO 450 J=1,NC
          TT(J) = R(I,J)
  450   CONTINUE
        DO 470 K=1,NC
          RT=ZERO
          DO 460 L=1,NC
            RT=RT + TT(L)*VEC(L,K)
  460     CONTINUE
          R(I,K)=RT
  470   CONTINUE
  475 CONTINUE
C-------------------------------------------------------------------
      IF (ICONV.GT.0) GOTO 900
C-------------------------------------------------------------------
C       EVALUATE TOLERANCES - COUNT FOR ABSOLUTE POSTIVE VALUES
C-------------------------------------------------------------------
      NNEG = 0
      NPOS = 0
      DO 480 I=1,NC
         DIF=ABS(EIGV(I)-D(I))
         RTOLV(I)=DIF/ABS(EIGV(I))
         EIGVAL = EIGV(I) + SHIFT
         IF (EIGVAL.GE.ZERO) NPOS=NPOS+1
  480 CONTINUE
C-------------------------------------------------------------------
C       ITERATION PRINTOUT 
C-------------------------------------------------------------------
      IF (IPRINT.GT.0) THEN
        WRITE (LINE,'(A)') '  COUNT       EIGENVALUES     TOLERANCES'
        WRITE(INFO,*)LINE
        DO 490 I=1,NC
          EV = EIGV(I)
          WRITE (INFO,'(1P,I6,4X,2E15.6)') I,EV,RTOLV(I)
  490   CONTINUE
      ENDIF
C---------------------------------------------------------------------
C       STOP ITERATION IF ALL EIGENVALUES ARE ABSOLUTLY NEGATIVE
C       RELATIVLY NEGATIVE WITH RESPECT TO SHIFT IS ALLOWED
C---------------------------------------------------------------------       
      IF (NPOS.EQ.0 .AND. NITE.GT.4) THEN
        IST = -1
        WRITE(*,*)'EIG005: ONLY NEGATIVE EIGENVALUES FOUND AFTER '//
     +         '4 ITERATIONS - SOLUTION STOPPED'
        GO TO 900
      ENDIF
C--------------------------------------------------------------------
C       C O N V E R G E N C E    C H E C K
C--------------------------------------------------------------------
C       NCONVR ... NUMBER OF CONVERGED EIGENVALUES
C       NNEG   ... NUMBER OF NEGATIVE CONVERGED EIGENVALUES
C--------------------------------------------------------------------
      NNEG   = 0
      NCONVR = 0
      DO 500 I=1,NROOT
        IF (RTOLV(I).GT.TOLEIG) GOTO 520
        IF (EIGV(I).LT.ZERO) NNEG = NNEG+1
        NCONVR = I
  500 CONTINUE
      ICONV=1
      IST = 0
      GOTO 100
C------------------------------------------------ CHECK ITERATION LIMIT
  520 CONTINUE
      IF (NITE.LT.NITEM) GOTO 530
      ICONV=2
      IST  =2
      GOTO 100
C------------------------------------------------------ NEXT ITERATION
  530 DO 540 I=1,NC
         D(I)=EIGV(I)
  540 CONTINUE
      GOTO 100
C---------------------------------------------------------------------
  900 RETURN
      END
C============================================================================
      SUBROUTINE SSPSTA (A,B,MAXA,R,W,TT,NN,NC,NWK,NWM,NNM,
     *                  ISOL,NNZ)
C =     P R O G R A M                                                =
C =          TO ESTABLISH STARTING ITERATION VECTORS FOR SUBSPACE    =
C =          ITERATION                                               =
C =                                                                  =
C =    WHERE :                                                       =
C =     A (NWK)  --->      VECTOR CONTAINING THE STIFFNESS MATRIX    =
C =     B (NWM)  --->      VECTOR CONTAINING THE MASS MATRIX         =
C =                     NWM EQ. NN    ..  LUMPED MASSES              =
C =                     NWM EQ. NWK   ..  CONSISTENT MASSES          =
C =     W (NN)   --->      WORKING VECTOR                            =
C =     TT(NN)   --->      WORKING VECTOR                            =
C =                     (CONTAINS THE EXCITED DEGREES OF FREEDOM ON  =
C =                      EXIT )                                      =
C =     MAXA (NWM) --->  ADRESSES OF THE DIAGONAL ELEMENTS IN A AND  =
C =                      B (CONS. MASSES)                            =
C =                                                                  =
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NWK),B(NWM),W(NN),TT(NN),R(NN,NC)
      INTEGER MAXA(NNM)
C
C
C     ESTABLISH STARTING ITERATION VECTORS
C
      ND=NN/NC
      IF (ISOL.NE.4) THEN
C      
        IF (NWM.GT.NN) GOTO4
        J=0
        DO 2 I=1,NN
           II=MAXA(I)
           R(I,1)=B(I)
           IF (B(I).GT.0.) J=J+1
           W(I)=B(I)/A(II)
 2      CONTINUE
        IF (NC.LE.J) GOTO 16
 4      J=0
        DO 10 I=1,NN
            II=MAXA(I)
            R(I,1)=B(II)
            IF (B(II).GT.0.) J=J+1
            W(I)=B(II)/A(II)
 10     CONTINUE
C
C
      ENDIF      
 16   DO 20 J=2,NC
      DO 20 I=1,NN
            R(I,J)=0.0D0
 20   CONTINUE
C
      L=NN-ND
      DO 30 J=2,NC
            RT=0.
            DO 40 I=1,L
                  IF (W(I).LT.RT) GOTO 40
                  RT=W(I)
                  IJ=I
 40         CONTINUE
            DO 50 I=L,NN
                  IF (W(I).LE.RT) GOTO 50
                  RT=W(I)
                  IJ=I
 50         CONTINUE
            TT(J)=REAL(IJ)
            W(IJ)=0.
            L=L-ND
            R(IJ,J)=1.
 30   CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MXCR8 (A,M,N,R)
C-----------------------------------------------------------------------        
C     KOPIERT A NACH R
C-----------------------------------------------------------------------        
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),R(1)
C
      MN=M*N
      DO 1 I=1,MN
      R(I)=A(I)
    1 CONTINUE
      RETURN
      END
C=======================================================================
C***********************************************************************
C   File:  @(#)solver_control.ftn	1.10   10/06/98   
C************************************************************************
      SUBROUTINE CALSOL(OP,ISOL,ITIMSOL,ILIN,NEQ,A,MAXA,NWK,NN,NNM,NSCH,
     1                   TT)
C      SUBROUTINE CALSOL(OP,ISOL,ITIMSOL,ILIN,NEQ,A,MAXA,NWK,NN,NNM,NSCH,
C     1                   TT)
C     +----------------------------------------------------------------+
C     |                                                                |
C     |                     INSTITUT FUER BAUSTATIK                    |
C     |                     UNIVERSITAET  STUTTGART                    |
C     |                                                                |
C     +----------------------------------------------------------------+
C     | MODUL CALSOL CONTAINS ALL PROGRAMS TO SOLVE THE  LIN. EQUATION |
C     |                                                                |
C     |                 MAT(NEQ,NEQ) * S(NEQ) = RHS(NEQ)               |
C     |                                                                |
C     | VERSION     : 1.3           PROGRAMMED BY : STEFAN     KIMMICH |
C     |                                             KAI-UWE BLETZINGER |
C     | DATE        : 29.05.1987                    HANS   STEGMUELLER |
C     |                                                                |
C     +----------------------------------------------------------------+
C     |  MODULE OPERATION PARAMETER (OP):                              |
C     |                                                                |
C     | 'INIT'   -->  INITIALIZE LOCAL DATA-MANAGEMENT FOR >CALSOL<    |
C     |               DEFINES ARRAYS MAT,RHS AND MAXA IF NEEDED        |
C     | 'TRI'    -->  REDUCTION OF SYSTEM MATRIX TO TRIANGULAR FORM    |
C     | 'RED'    -->  BACKSUBSTITUTION OF RIGHTHANDSIDE VECTOR         |
C     | 'ALL'    -->  EXECUTE REDUCTION AND BACKSUBSTITUTION TOGETHER  |
c     | 'CONT'   -->  redefines arrays mat,rhs and maxa if needed      |
C     +----------------------------------------------------------------+
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*4 OP
c      CHARACTER*8 OP*(*)
c      PARAMETER (ZERO=0.0D0)
c      PARAMETER (ISC=1)
C      SAVE IUMF
      DIMENSION A(NWK),RHS(1,1),MAXA(NNM),TT(NN)
      ICPRINT=6
      ZERO=0.0D0
      ISC=1
C---------------------------------------------------- CHECK COMMAND LINE
      IPR=0
      INIT=INDEX(OP,'INIT')
      IRED=INDEX(OP,'RED')
      IF (IRED.EQ.0) IRED=INDEX(OP,'ALL')
      ITRI=INDEX(OP,'TRI')
      IF (ITRI.EQ.0) ITRI=INDEX(OP,'ALL')
C - Abfrage linear zum Ausschalten der Dreieckszerlegung fuer it. Loeser
C ----------(nicht linear nicht noetig, da Zerlegung anders gesteuert !)
      IF (ISOL.EQ.2 .AND. ITRI.GT.0 .AND. ILIN.EQ.1) GOTO 1000
C
      ICONT=INDEX(OP,'CONT')
      IRHS=INDEX(OP,'RHS')
      ISYS=INDEX(OP,'SYS')
      IF (IRHS.GT.0 .OR. ISYS.GT.0) THEN
C         CALL FRLINE (OP,0)?????????????????????????????????????????
C         IF (IRHS.GT.0) CALL FRCHAR ('RHS',NAMRHS,8,1,IRC)
C         IF (ISYS.GT.0) CALL FRCHAR ('SYS',NAMMAT,8,1,IRC)
      ENDIF
C-----------------------------------------------------------------------
C   GET ALL REQUIRED GLOBAL VARIABLES BEFORE ACTIVATING BLOCK SOLVAR
C-----------------------------------------------------------------------
C------------------------------------------------------ INITIALIZE MODUL
C***********************************************************************
C       EXECUTION BLOCK : ISOL=1 : OPTIONS ARE TRI,RED,ALL 
C                         ISOL=2 : NO OPTIONS (ALL)
C                         ISOL=3 : OPTIONS ARE TRI,RED,ALL 
C                         ISOL=4/5 : OPTIONS ARE TRI,RED,ALL 
C                         ISOL=9 : NO OPTIONS (ALL)
C***********************************************************************
         IF (ISOL.EQ.1 .OR. ISOL.EQ.3 .OR. ISOL.EQ.4 .OR. ISOL.EQ.5)
     +   THEN
                           IACT=0
            IF (ITRI.GT.0) IACT=1
            IF (IRED.GT.0) IACT=2
            IF (ITRI.GT.0 .AND. IRED.GT.0) IACT=3
         ENDIF
         IF (ISOL.EQ.2) IACT=3
         IF (ISOL.EQ.9) IACT=3
         IF (IACT.EQ.0) WRITE(*,*)'CALSOL',OP
C-----------------------------------------------------------------------
C        GET NAME,DIMENSION AND ADDRESS OF SYSTEM MATRIX AND CHECK
C        IF ISYS>0 NAME WAS GIVEN BY COMMAND LINE
C-----------------------------------------------------------------------
         NEQMAT=NEQ
         IF (NEQMAT.LE.0) WRITE(*,*)'CALSOL','NEQMAT.LE.0 !'
C-----------------------------------------------------------------------
C        GET NAME,DIMENSION AND ADDRESS OF RHS AND CHECK
C        IF IRHS>0 NAME WAS GIVEN BY COMMAND LINE
C-----------------------------------------------------------------------
         IF (IACT.GT.1) THEN
            IF (IRHS.EQ.0) THEN
             NEQRHS= NN
             NUMRHS= 1
             N1RHS = 1
             N2RHS = 1
             NRHS  = 1
            ENDIF
         ELSE
            NEQRHS= 1
            NUMRHS= 1
            N1RHS = 1
            N2RHS = 1
            NRHS  = 1
         ENDIF
C----------------------------------------------- GET TYPE OF SOLVER USED
C
C****************************************** SOLVE LINEAR EQUATION SYSTEM
C************************************************** WITH SOLVER >COLSOL<
         IF (IACT.GT.1) THEN
               CALL COLSOL (A,TT,MAXA,NEQ,
     *                      NEQRHS,NUMRHS,NWK,NNM,N1RHS,N2RHS,IACT,
     *                      DET,ISC,NSCH,IPR,ICPRINT)
         ELSE
               CALL COLSOL (A,RHS,MAXA,NEQ,
     *                      NEQRHS,NUMRHS,NWK,NNM,N1RHS,N2RHS,IACT,
     *                      DET,ISC,NSCH,IPR,ICPRINT)
         ENDIF
C-------------------------------------------------------- SAVE VARIABLES
c               IF (IACT.EQ.1 .OR. IACT.EQ.3) THEN
c                 CALL VMGET ('DETINI',DET0)
c                 IF (DET0.EQ.ZERO) THEN
c                   DET0=DET
c                   CALL VMPUT ('DETINI',DET0)
c                 ENDIF
c                 DET = DET/DET0
c                 CALL VMPUT ('DETMAT',DET)
c                 CALL VMPUT ('ISC',ISC)
c                 CALL VMPUT ('NSCH',NSCH)
c               ENDIF
C
 1000 RETURN
      END
C=======================================================================
      SUBROUTINE XMULT (TT,B,RR,MAXA,NN,NWM,NNM)
C......................................................................
C
C     P R O G R A M
C          TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT
C
C......................................................................
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TT(NN),B(NWM),RR(NN),MAXA(NNM)
C
      IF (NWM.GT.NN) GOTO 20
      DO 10 I=1,NN
            TT(I)=B(I)*RR(I)
 10   CONTINUE
      RETURN
C
 20   DO 40 I=1,NN
            TT(I)=0.
 40   CONTINUE
      DO 100 I=1,NN
             KL=MAXA(I)
             KU=MAXA(I+1) - 1
             II=I + 1
             CC=RR(I)
             DO 110 KK=KL,KU
                    II=II - 1
                    TT(II)=TT(II) + B(KK)*CC
 110         CONTINUE
 100  CONTINUE
      IF (NN.EQ.1) THEN
        RETURN
      ENDIF
      DO 200 I=2,NN
             KL=MAXA(I) + 1
             KU=MAXA(I+1) - 1
             IF (KU-KL) 200,210,210
 210         II=I
             AA=0.
             DO 220 KK=KL,KU
                    II=II - 1
                    AA=AA + B(KK)*RR(II)
 220         CONTINUE
             TT(I)=TT(I) + AA
 200  CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MXSDAT (VEC,INDNUM,LVEC,NCOL,INDCOL,IRC)
C-----------------------------------------------------------------------        
C     TOPIC : REARRANGE SORTED DATA
C             DATA IS TO BE MOVED IN CORRECT POSITIONS
C-----------------------------------------------------------------------        
C     VEC     ...  VECTOR TO BE REARRANGED - MAY HAVE MULTIPLE COLUMNS
C     INDNUM  ...  INDNUM ARRAY CONTAINS POINTERS 
C     LVEC    ...  COLUMN LENGTH OF VEC AND INDNUM     
C     NCOL    ...  NUMBER OF COLUMNS OF VEC
C     INDCOL  ...  NUMBER OF COLUMNS OF INDNUM
C     IRC     ...  RETURN CODE (CURRENTLY NOT USED)
C-----------------------------------------------------------------------        
C     PROGRAMMED BY : HANS STEGMUELLER  NOV,1988
C-----------------------------------------------------------------------        
      INTEGER LVEC,NCOL,INDCOL,IRC
      REAL*8  VEC(NCOL,LVEC),VAL
      INTEGER INDNUM(LVEC,INDCOL)
      INTEGER I,K,M,N,SORT
C-----------------------------------------------------------------------        
      IRC = 1
C-----------------------------------------------------------------------        
C       RESTORE INDNUM COLUMN IF SECOND OR FURTHER REARRANGEMENT
C       IS TO BE DONE - AFTER FIRST RUN INDNUM DATA ARE DESTROYED
C-----------------------------------------------------------------------        
      IF (INDCOL.GT.1) THEN
        DO 10 K=1,LVEC
          INDNUM(K,1) = INDNUM(K,2)
   10   CONTINUE
      ENDIF
C-------------------------------------------------- START REARRANGE LOOP
      DO 500 M=1,100000
C--------------------------------------------------- CLEAR EXCHANGE FLAG
      SORT = 0
C----------------------------------------- LOOP OVER ALL VECTOR ELEMENTS
      DO 50 K=1,LVEC
C-----------------------------------------------------------------------        
C       IF INDNUM VALUE DOES NOT MATCH COUNTER VALUE 
C       DATA ARE NOT IN CORRECT ORDER AND MUST BE REARRANGED
C-----------------------------------------------------------------------        
        IF (INDNUM(K,1).NE.K) THEN
           I   = INDNUM(K,1)
           INDNUM(K,1) = INDNUM(I,1)
           INDNUM(I,1) = I
C---------------------------------------- REARRANGE DATA FOR ALL COLUMNS
           DO 40 N=1,NCOL
             VAL = VEC(N,K)
             VEC(N,K)   = VEC(N,I)
             VEC(N,I)   = VAL
   40      CONTINUE
C----------------------------------------------------- SET EXCHANGE FLAG
           SORT = 1
        ENDIF
   50 CONTINUE
C
      IF (SORT.EQ.0) GO TO 600
C
  500 CONTINUE
      WRITE (6,*) 'MXRARR: MAXIMUM SORT COUNTER REACHED - STOP'
      IRC = 1
  600 CONTINUE
C
  900 RETURN
      END
C=======================================================================
C***********************************************************************
C   File:  @(#)eigen_qz_solver.ftn	1.2   06/23/98   
C***********************************************************************
      SUBROUTINE SSPEIG (A,B,VEC,EIGV,BETA,ALFR,ITER,N,IERROR)
C   +---------------------------------------------------------------+
C   -    QZ-ALGORITHM FROM-
C   -    C.B.MOLER AND G.W.STEWART                                  -
C   -    REPORT NU. CNA-32                                          -
C   -    CENTER FOR NUMER. ANALYSIS                                 -
C   -    THE UNIVERSITY OF TEXAS AT AUSTIN :                        -
C   -    AN ALGORITHM FOR THE GENERALIZED MATRIX                    -
C   -    EIGENVALUE PROBLEM A*X=LAM.*B*X                            -
C   +---------------------------------------------------------------+
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(N,N),B(N,N),VEC(N,N),EIGV(N),ALFR(N),BETA(N)
      INTEGER ITER(N)
C   ------------------------------------------------------------------
C
      CALL EIGW(N,A,B,ALFR,BETA,VEC,ITER,IERROR)
      IF (IERROR.EQ.1) GOTO 999
C
      DO 110 I=1,N
       IF(ITER(I).NE.-1) GO TO 110
       WRITE(*,*)'SSPEIG'
  110 CONTINUE
      DO 200 I=1,N
      IF(ABS(BETA(I)).LT.1.0E-35) GO TO 190
      EIGV(I)=ALFR(I)/BETA(I)
      GO TO 200
  190 EIGV(I)=1.0E+34
  200 CONTINUE
      DO  350 I=1, N
      IF(EIGV(I).EQ.1.0E+34) GO TO 350
      DEIG=EIGV(I)
      DO  360 J=1, N
  360 VEC(J,I)=VEC(J,I)*DEIG
  350 CONTINUE
  999 RETURN
      END
C==================================================================
      SUBROUTINE MSHIFT (A,B,MAXA,NWK,NWM,NNM,NN,SHIFT)
C
C     P R O G R A M
C           TO SHIFT MATRIX A  :    A --->  A - SHIFT * B
C
       IMPLICIT REAL*8(A-H,O-Z)
       DIMENSION A(NWK),B(NWM)
       INTEGER MAXA(NNM)
          IF (NWM.GT.NN) GOTO 645
         DO 640 I=1,NN
             II=MAXA(I)
             A(II)=A(II) - B(I)*SHIFT
 640     CONTINUE
         RETURN
 645     DO 650 I=1,NWK
             A(I)=A(I) - B(I)*SHIFT
 650     CONTINUE
          RETURN
         END
C=======================================================================
      SUBROUTINE MXSIND (VEC,INDNUM,LVEC,INDCOL,TYP,IRC)
C-----------------------------------------------------------------------        
C     TOPIC : SORT DATA GIVEN IN VECTOR VEC 
C             NO DATA MOVEMENT IS DONE - CREATE ONLY INDNUM ARRAY
C-----------------------------------------------------------------------        
C     VEC     ... VECTOR TO BE SORTED
C     INDNUM  ... INDNUM ARRAY TO STORE POINTERS
C     LVEC    ... LENGTH OF VECTOR VEC
C     TYP     ... TYPE OF SORTING (A=ASCENDING)
C                                 D=DESCENDING
C                                 N=RELATIVE TO ZERO)
C     IRC     ...  RETURN CODE 
C-----------------------------------------------------------------------        
C     PROGRAMMED BY : HANS STEGMUELLER  NOV,1988
C-----------------------------------------------------------------------        
      INTEGER LVEC,INDCOL,IRC
      REAL*8  VEC(LVEC),START
      INTEGER INDNUM(LVEC,INDCOL)
      INTEGER POINT,COUNT,I,K,ISTART
      CHARACTER*(*) TYP,UPP*3,LOW*3,CTYP*1
      SAVE UPP,LOW
      DATA  UPP/'ADN'/ ,LOW/'ADN'/
C-----------------------------------------------------------------------        
      IRC = 0
C--------------------------------------------- CONVERT TYP TO UPPER CASE
      I = INDEX(LOW,TYP)
      IF (I.GT.0) THEN
        CTYP = UPP(I:I)
      ELSE
        CTYP = TYP
      ENDIF
C---------------------------------------------- INITIALIZE INDNUM ARRAYS
      DO 20 I=1,INDCOL
        DO 10 K=1,LVEC
         INDNUM(K,I) = 0
   10   CONTINUE
   20 CONTINUE
C-----------------------------------------------------------------------        
      COUNT = 0
C----------------------------------- START LOOP OVER ALL VECTOR ELEMENTS
      DO 200 K=1,LVEC
C---------------------------------------------- FIND NEXT UNSORTED VALUE
        DO 50 I=1,LVEC
          IF (INDNUM(I,1).EQ.0) THEN
             IF (CTYP.EQ.'N') THEN
               START  = ABS(VEC(I))
             ELSE
               START  = VEC(I)
             ENDIF
             POINT  = I
             ISTART = I+1
             GO TO 60
          ENDIF
   50   CONTINUE
   60   CONTINUE
C----------------------------------- START INNER LOOP TO FIND NEXT VALUE
        DO 100 I=ISTART,LVEC
          IF (INDNUM(I,1).GT.0) GO TO 100
          IF (CTYP.EQ.'N') THEN
            IF (ABS(VEC(I)).LT.START) THEN
               START = ABS(VEC(I))
               POINT = I
            ENDIF
          ELSE IF (CTYP.EQ.'A') THEN
            IF (VEC(I).LT.START) THEN
               START = VEC(I)
               POINT = I
            ENDIF
          ELSE IF (CTYP.EQ.'D') THEN
            IF (VEC(I).GT.START) THEN
               START = VEC(I)
               POINT = I
            ENDIF
          ELSE
            WRITE (6,*) 'MXSORT: INVALID SORT TYPE GIVEN - ',TYP
            IRC = 1
            GO TO 900
          ENDIF
  100   CONTINUE
C-----------------------------------------------------------------------        
C         IF INDNUM ARRAY CONTAINS MORE THAN ONE COLUMN
C         STORE SORT INFORMATION IN FIRST AND SECOND COLUMN 
C         SECOND COLUMN IS TO BE USED FOR REPEATED DATA 
C         REARRANGEMENTS BECAUSE INDNUM INFORMATION IS DESTROYED
C         DURING REARRANEGMENT PROCEDURE
C         3.COLUMN CONTAINS POINTERS TO ACCESS DATA WITH DIRECT
C         ADDRESSES (INDNUM(1) POINTS TO FIRST SORTED ELEMENT, ...)
C-----------------------------------------------------------------------        
                         COUNT           = COUNT + 1
                         INDNUM(POINT,1) = COUNT
        IF (INDCOL.GT.1) INDNUM(POINT,2) = COUNT        
        IF (INDCOL.GT.2) INDNUM(COUNT,3) = POINT
C
  200 CONTINUE            
C
  900 RETURN
      END
C=======================================================================
      SUBROUTINE EIGW (N,A,B,ALFR,BETA,X,ITER,IERROR)
C
C        PROGRAMMTEILE FUER LOESUNG BEI UNSYMMETRISCHEN
C        MATRIZEN A UND B ENTFERNT
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ITER(N)
      DIMENSION A(N,N),B(N,N),ALFR(N),BETA(N),X(N,N)
      DATA EPS/1.0D-13/
      CALL QZHES (N,N,A,B,X)
      CALL QZIT  (N,N,A,B,EPS,EPSA,EPSB,ITER,X)
      CALL QZVAL (N,N,A,B,EPSB,ALFR,BETA,X,IERROR)
      IF (IERROR.EQ.1) GOTO 999
      CALL QZVEC(N,N,A,B,EPSA,EPSB,ALFR,BETA,X)
  999 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE QZHES (ND,N,A,B,X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(ND,ND),B(ND,ND),X(ND,ND)
C
C  INITIALIZE X, USED TO SAVE TRANSFORMATIONS
C
      DO 3 I=1,N
         DO 2 J=1,N
            X(I,J) = 0.
    2    CONTINUE
         X(I,I) = 1.
    3 CONTINUE
C
C  REDUCE B TO UPPER TRIANGULAR
C
      NM1=N-1
      DO 100 L=1,NM1
         L1 = L+1
         S = 0.
         DO 20 I=L1,N
         S=S+ABS(B(I,L))
   20    CONTINUE
         IF (S.EQ.0.) GO TO 100
         S=S+ABS(B(L,L))
         R = 0.
         DO 25 I=L,N
            B(I,L) = B(I,L)/S
            R = R + B(I,L)**2
   25    CONTINUE
         R = SQRT(R)
         IF (B(L,L).LT.0.) R = -R
         B(L,L) = B(L,L) + R
         RHO = R*B(L,L)
         DO 50 J=L1,N
            T = 0.
            DO 30 I=L,N
               T = T + B(I,L)*B(I,J)
   30       CONTINUE
            T = -T/RHO
            DO 40 I=L,N
               B(I,J) = B(I,J) + T*B(I,L)
   40       CONTINUE
   50    CONTINUE
         DO 80 J=1,N
            T = 0.
            DO 60 I=L,N
               T = T + B(I,L)*A(I,J)
   60       CONTINUE
            T = -T/RHO
            DO 70 I=L,N
               A(I,J) = A(I,J) + T*B(I,L)
   70       CONTINUE
   80    CONTINUE
         B(L,L) = -S*R
         DO 90 I=L1,N
            B(I,L) = 0.
   90    CONTINUE
  100 CONTINUE
C
C  REDUCE A TO UPPER HESSENBERG, KEEP B TRIANGULAR
C
      IF (N.LE.2) GO TO 170
      NM2=N-2
      DO 160 K=1,NM2
         K1 = K+1
         NK1=N-K1
         DO 150 LB=1,NK1
            L = N-LB
            L1 = L+1
            CALL HSH2(A(L,K),A(L1,K),U1,U2,V1,V2)
            IF (U1.NE.1.) GO TO 125
            DO 110 J=K,N
               T = A(L,J) + U2*A(L1,J)
               A(L,J) = A(L,J) + T*V1
               A(L1,J) = A(L1,J) + T*V2
  110       CONTINUE
            A(L1,K) = 0.
            DO 120 J=L,N
               T = B(L,J) + U2*B(L1,J)
               B(L,J) = B(L,J) + T*V1
               B(L1,J) = B(L1,J) + T*V2
  120       CONTINUE
  125       CALL HSH2(B(L1,L1),B(L1,L),U1,U2,V1,V2)
            IF (U1.NE.1.) GO TO 150
            DO 130 I=1,L1
               T = B(I,L1) + U2*B(I,L)
               B(I,L1) = B(I,L1) + T*V1
               B(I,L) = B(I,L) + T*V2
  130       CONTINUE
            B(L1,L) = 0.
            DO 140 I=1,N
               T = A(I,L1) + U2*A(I,L)
               A(I,L1) = A(I,L1) + T*V1
               A(I,L) = A(I,L) + T*V2
  140       CONTINUE
            DO 145 I=1,N
               T = X(I,L1) + U2*X(I,L)
               X(I,L1) = X(I,L1) + T*V1
               X(I,L) = X(I,L) + T*V2
  145       CONTINUE
  150       CONTINUE
  160 CONTINUE
  170 CONTINUE
      RETURN
      END
      SUBROUTINE QZIT (ND,N,A,B,EPS,EPSA,EPSB,ITER,X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(ND,ND),B(ND,ND),X(ND,ND)
      DIMENSION ITER(N)
      LOGICAL MID
C
C  INITIALIZE ITER, COMPUTE EPSA,EPSB
C
      ANORM = 0.
      BNORM = 0.
      DO 185 I=1,N
         ITER(I) = 0
         ANI = 0.
         IF (I.NE.1) ANI = ABS(A(I,I-1))
         BNI = 0.
         DO 180 J=I,N
            ANI = ANI + ABS(A(I,J))
            BNI = BNI + ABS(B(I,J))
  180    CONTINUE
         IF (ANI.GT.ANORM) ANORM = ANI
         IF (BNI.GT.BNORM) BNORM = BNI
  185 CONTINUE
      IF(ANORM.EQ.0.) ANORM=EPS
      IF(BNORM.EQ.0.) BNORM=EPS
      EPSA = EPS*ANORM
      EPSB = EPS*BNORM
C
C  REDUCE A TO QUASI-TRIANGULAR, KEEP B TRIANGULAR
C
      M = N
  200 IF (M.LE.2) GO TO 390
C
C     CHECK FOR CONVERGENCE OR REDUCIBILITY
C
      DO 220 LB=1,M
         L = M+1-LB
         IF (L.EQ.1) GO TO 260
         IF (ABS(A(L,L-1)) .LE. EPSA) GO TO 230
  220 CONTINUE
  230 A(L,L-1) = 0.
      IF (L.LT.M-1) GO TO 260
      M = L-1
      GO TO 200
C
C     CHECK FOR SMALL TOP OF B
C
  260 IF (ABS(B(L,L)).GT.EPSB) GO TO 300
      B(L,L) = 0.
      L1 = L+1
      CALL HSH2(A(L,L),A(L1,L),U1,U2,V1,V2)
      IF (U1.NE.1.) GO TO 280
      DO 270 J=L,N
         T = A(L,J) + U2*A(L1,J)
         A(L,J) = A(L,J) + T*V1
         A(L1,J) = A(L1,J) + T*V2
         T = B(L,J) + U2*B(L1,J)
         B(L,J) = B(L,J) + T*V1
         B(L1,J) = B(L1,J) + T*V2
  270 CONTINUE
  280 L = L1
      GO TO 230
C
C     BEGIN ONE QZ STEP, ITERATION STRATEGY
C
  300 M1 = M - 1
      L1 = L + 1
      CONST = 0.75
      ITER(M) = ITER(M) + 1
      IF (ITER(M).EQ.1) GO TO 305
      IF (ABS(A(M,M-1)).LT.CONST*OLD1) GO TO 305
      IF (ABS(A(M-1,M-2)).LT.CONST*OLD2) GO TO 305
      IF (ITER(M).EQ.10) GO TO 310
      IF (ITER(M).GT.30) GO TO 380
C
C     ZEROTH COLUMN OF A
C
  305 B11 = B(L,L)
      B22 = B(L1,L1)
      IF (ABS(B22).LT.EPSB) B22 = EPSB
      B33 = B(M1,M1)
      IF (ABS(B33).LT.EPSB) B33 = EPSB
      B44 = B(M,M)
      IF (ABS(B44).LT.EPSB) B44 = EPSB
      A11 = A(L,L)/B11
      A12 = A(L,L1)/B22
      A21 = A(L1,L)/B11
      A22 = A(L1,L1)/B22
      A33 = A(M1,M1)/B33
      A34 = A(M1,M)/B44
      A43 = A(M,M1)/B33
      A44 = A(M,M)/B44
      B12 = B(L,L1)/B22
      B34 = B(M1,M)/B44
      A10 = ( (A33-A11)*(A44-A11) - A34*A43 + A43*B34*A11 )/A21
     1      + A12 - A11*B12
      A20 = (A22-A11-A21*B12) - (A33-A11) - (A44-A11) + A43*B34
      A30 = A(L+2,L1)/B22
      GO TO 315
C
C     AD HOC SHIFT
C
  310 A10 = 0.
      A20 = 1.
      A30 = 1.1605
C
  315 OLD1 = ABS(A(M,M-1))
      OLD2 = ABS(A(M-1,M-2))
C
C     BEGIN MAIN LOOP
C
      DO 360 K=L,M1
         MID = K.NE.M1
         K1 = K+1
         K2 = K+2
         K3 = K+3
         IF (K3.GT.M) K3 = M
         KM1 = K-1
         IF (KM1.LT.L) KM1 = L
         IF (K.EQ.L) CALL HSH3(A10,A20,A30,U1,U2,U3,V1,V2,V3)
         IF (K.GT.L.AND.K.LT.M1)
     1     CALL HSH3(A(K,KM1),A(K1,KM1),A(K2,KM1),U1,U2,U3,V1,V2,V3)
         IF (K.EQ.M1) CALL HSH2(A(K,KM1),A(K1,KM1),U1,U2,V1,V2)
         IF (U1.NE.1.) GO TO 325
         DO 320 J=KM1,N
            T = A(K,J) + U2*A(K1,J)
            IF (MID) T = T + U3*A(K2,J)
            A(K,J) = A(K,J) + T*V1
            A(K1,J) = A(K1,J) + T*V2
            IF (MID) A(K2,J) = A(K2,J) + T*V3
            T = B(K,J) + U2*B(K1,J)
            IF (MID) T = T + U3*B(K2,J)
            B(K,J) = B(K,J) + T*V1
            B(K1,J) = B(K1,J) + T*V2
            IF (MID) B(K2,J) = B(K2,J) + T*V3
  320    CONTINUE
         IF (K.EQ.L) GO TO 325
         A(K1,K-1) = 0.
         IF (MID) A(K2,K-1) = 0.
  325    IF (K.EQ.M1) GO TO 340
         CALL HSH3(B(K2,K2),B(K2,K1),B(K2,K),U1,U2,U3,V1,V2,V3)
         IF (U1.NE.1.) GO TO 340
         DO 330 I=1,K3
            T = A(I,K2) + U2*A(I,K1) + U3*A(I,K)
            A(I,K2) = A(I,K2) + T*V1
            A(I,K1) = A(I,K1) + T*V2
            A(I,K) = A(I,K) + T*V3
            T = B(I,K2) + U2*B(I,K1) + U3*B(I,K)
            B(I,K2) = B(I,K2) + T*V1
            B(I,K1) = B(I,K1) + T*V2
            B(I,K) = B(I,K) + T*V3
  330    CONTINUE
         B(K2,K) = 0.
         B(K2,K1) = 0.
         DO 335 I=1,N
            T = X(I,K2) + U2*X(I,K1) + U3*X(I,K)
            X(I,K2) = X(I,K2) + T*V1
            X(I,K1) = X(I,K1) + T*V2
            X(I,K) = X(I,K) + T*V3
  335    CONTINUE
  340    CALL HSH2(B(K1,K1),B(K1,K),U1,U2,V1,V2)
         IF (U1.NE.1.) GO TO 360
         DO 350 I=1,K3
            T = A(I,K1) + U2*A(I,K)
            A(I,K1) = A(I,K1) + T*V1
            A(I,K) = A(I,K) + T*V2
            T = B(I,K1) + U2*B(I,K)
            B(I,K1) = B(I,K1) + T*V1
            B(I,K) = B(I,K) + T*V2
  350    CONTINUE
         B(K1,K) = 0.
         DO 355 I=1,N
            T = X(I,K1) + U2*X(I,K)
            X(I,K1) = X(I,K1) + T*V1
            X(I,K) = X(I,K) + T*V2
  355    CONTINUE
  360    CONTINUE
C
C     END MAIN LOOP
C
      GO TO 200
C
C     END ONE QZ STEP
C
  380 DO 385 I=1,M
         ITER(I) = -1
  385 CONTINUE
  390 CONTINUE
      RETURN
      END
      SUBROUTINE QZVAL (ND,N,A,B,EPSB,ALFR,BETA,X,IERROR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C      PARAMETER (ZERO = 0.0D0)
      DIMENSION A(ND,ND),B(ND,ND),ALFR(N),BETA(N),X(ND,ND)
      LOGICAL FLIP
C
C  FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES
C
      ZERO = 0.0D0
      M = N
  400 CONTINUE
         IF (M.EQ.1) GO TO 410
         IF (A(M,M-1).NE.0.) GO TO 420
C
C        ONE-BY-ONE SUBMATRIX, ONE REAL ROOT
C
  410    ALFR(M) = A(M,M)
         IF(B(M,M).LT.0.0) ALFR(M)=-ALFR(M)
         BETA(M)=ABS(B(M,M))
         M = M-1
         GO TO 490
C
C        TWO-BY-TWO SUBMATRIX
C
  420    L = M-1
         IF (ABS(B(L,L)).GT.EPSB) GO TO 425
            B(L,L) = 0.
            CALL HSH2(A(L,L),A(M,L),U1,U2,V1,V2)
            GO TO 460
 425    IF (ABS(B(M,M)).GT.EPSB) GO TO 430
            B(M,M) = 0.
            CALL HSH2(A(M,M),A(M,L),U1,U2,V1,V2)
            BN = 0.
            GO TO 435
  430    AN = ABS(A(L,L))+ABS(A(L,M))+ABS(A(M,L))+ABS(A(M,M))
         BN = ABS(B(L,L))+ABS(B(L,M))+ABS(B(M,M))
         A11 = A(L,L)/AN
         A12 = A(L,M)/AN
         A21 = A(M,L)/AN
         A22 = A(M,M)/AN
         B11 = B(L,L)/BN
         B12 = B(L,M)/BN
         B22 = B(M,M)/BN
         E=A11/B11
         C=((A22-E*B22)/B22-(A21*B12)/(B11*B22))/2.
         D=C*C+(A21*(A12-E*B12))/(B11*B22)
         IF (D.LT.ZERO) THEN
           IERROR = 1 
           GOTO 999
         ELSE
           IERROR = 0
         ENDIF 
C
C        TWO REAL ROOTS
C        ZERO BOTH A(M,L) AND B(M,L)
C
         IF(C.GE.0.) E=E+(C+SQRT(D))
         IF(C.LT.0.) E=E+(C-SQRT(D))
         A11 = A11 - E*B11
         A12 = A12 - E*B12
         A22 = A22 - E*B22
         FLIP = (ABS(A11)+ABS(A12)).GE.(ABS(A21)+ABS(A22))
         IF (FLIP) CALL HSH2(A12,A11,U1,U2,V1,V2)
         IF (.NOT.FLIP) CALL HSH2(A22,A21,U1,U2,V1,V2)
  435    IF (U1.NE.1.) GO TO 450
         DO 440 I=1,M
            T = A(I,M) + U2*A(I,L)
            A(I,M) = A(I,M) + V1*T
            A(I,L) = A(I,L) + V2*T
            T = B(I,M) + U2*B(I,L)
            B(I,M) = B(I,M) + V1*T
            B(I,L) = B(I,L) + V2*T
  440    CONTINUE
         DO 445 I=1,N
            T = X(I,M) + U2*X(I,L)
            X(I,M) = X(I,M) + V1*T
            X(I,L) = X(I,L) + V2*T
  445    CONTINUE
  450    IF (BN.EQ.0.) GO TO 475
         FLIP = AN .GE. ABS(E)*BN
         IF (FLIP) CALL HSH2(B(L,L),B(M,L),U1,U2,V1,V2)
         IF (.NOT.FLIP) CALL HSH2(A(L,L),A(M,L),U1,U2,V1,V2)
  460    IF (U1.NE.1.) GO TO 475
         DO 470 J=L,N
            T = A(L,J) + U2*A(M,J)
            A(L,J) = A(L,J) + V1*T
            A(M,J) = A(M,J) + V2*T
            T = B(L,J) + U2*B(M,J)
            B(L,J) = B(L,J) + V1*T
            B(M,J) = B(M,J) + V2*T
  470    CONTINUE
  475    A(M,L) = 0.
         B(M,L) = 0.
         ALFR(L) = A(L,L)
         ALFR(M) = A(M,M)
         IF(B(L,L).LT.0.) ALFR(L)=-ALFR(L)
         IF(B(M,M).LT.0.) ALFR(M)=-ALFR(M)
         BETA(L) =ABS(B(L,L))
         BETA(M) = ABS(B(M,M))
         M = M-2
C
  490 IF (M.GT.0) GO TO 400
  999 CONTINUE
      RETURN
      END
      SUBROUTINE QZVEC (ND,N,A,B,EPSA,EPSB,ALFR,BETA,X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(ND,ND),B(ND,ND),ALFR(N),BETA(N),X(ND,ND)
      LOGICAL FLIP
C
C  FIND EIGENVECTORS OF QUASI-TRIANGULAR MATRICES
C  USE B FOR INTERMEDIATE STORAGE
C
C     DO 500 THRU 590 FOR M = N STEP (-1 OR -2) UNTIL 1
C
      M = N
  500 CONTINUE
C
C        REAL VECTOR
C
         ALFM = ALFR(M)
         BETM = BETA(M)
         B(M,M) = 1.
C
C        DO 510 THRU 540 FOR L = M-1 STEP (-1 OR -2) UNTIL 1
C
         L = M-1
         IF (L.EQ.0) GO TO 540
  510    CONTINUE
            L1 = L+1
            SL = 0.
            DO 515 J=L1,M
               SL = SL + (BETM*A(L,J)-ALFM*B(L,J))*B(J,M)
  515       CONTINUE
            IF (L.EQ.1) GO TO 520
            IF(BETM*A(L,L-1).NE.0.) GO TO 530
  520       D = BETM*A(L,L)-ALFM*B(L,L)
            IF (D.EQ.0.) D = (EPSA+EPSB)/2.
            B(L,M) = -SL/D
            L = L-1
            GO TO 540
C
  530       K = L-1
            SK = 0.
            DO 535 J=L1,M
               SK = SK + (BETM*A(K,J)-ALFM*B(K,J))*B(J,M)
  535       CONTINUE
            TKK = BETM*A(K,K) - ALFM*B(K,K)
            TKL = BETM*A(K,L) - ALFM*B(K,L)
            TLK = BETM*A(L,K)
            TLL = BETM*A(L,L) - ALFM*B(L,L)
            D = TKK*TLL - TKL*TLK
            IF (D.EQ.0.) D = (EPSA+EPSB)/2.
            B(L,M) = (TLK*SK - TKK*SL)/D
            FLIP = ABS(TKK) .GE. ABS(TLK)
            IF (FLIP) B(K,M) = -(SK + TKL*B(L,M))/TKK
            IF (.NOT.FLIP) B(K,M) = -(SL + TLL*B(L,M))/TLK
            L = L-2
  540    IF (L.GT.0) GO TO 510
         M = M-1
      IF (M.GT.0) GO TO 500
C
C  TRANSFORM TO ORIGINAL COORDINATE SYSTEM
C
      M = N
  600 CONTINUE
         DO 620 I=1,N
            S = 0.
            DO 610 J=1,M
               S = S  +X(I,J)*B(J,M)
  610       CONTINUE
            X(I,M) = S
  620    CONTINUE
         M = M-1
      IF (M.GT.0) GO TO 600
C
C  NORMALIZE SO THAT LARGEST COMPONENT = 1.
C
      M = N
  630 CONTINUE
         S = 0.
         DO 635 I=1,N
            R = ABS(X(I,M))
            IF (R.LT.S) GO TO 635
            S = R
            D = X(I,M)
  635    CONTINUE
         DO 640 I=1,N
            X(I,M) = X(I,M)/D
  640    CONTINUE
         M = M-1
      IF (M.GT.0) GO TO 630
  700 RETURN
      END
      SUBROUTINE HSH2 (A1,A2,U1,U2,V1,V2)
C
C  FINDS HOUSEHOLDER TRANSFORMATION THAT WILL ZERO A2
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (A2.EQ.0.) GO TO 10
      S = ABS(A1) + ABS(A2)
      U1 = A1/S
      U2 = A2/S
      R = SQRT(U1*U1+U2*U2)
      IF (U1.LT.0.) R = -R
      V1 = -(U1 + R)/R
      V2 = -U2/R
      U1 = 1.
      U2 = V2/V1
      RETURN
   10 U1 = 0.
      RETURN
      END
      SUBROUTINE HSH3 (A1,A2,A3,U1,U2,U3,V1,V2,V3)
C
C  FINDS HOUSEHOLDER TRANSFORMATION THAT WILL ZERO A2 AND A3
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IF (A2.EQ.0. .AND.A3.EQ.0.) GO TO 10
      S = ABS(A1) + ABS(A2) + ABS(A3)
      U1 = A1/S
      U2 = A2/S
      U3 = A3/S
      R = SQRT(U1*U1+U2*U2+U3*U3)
      IF (U1.LT.0.) R = -R
      V1 = -(U1 + R)/R
      V2 = -U2/R
      V3 = -U3/R
      U1 = 1.
      U2 = V2/V1
      U3 = V3/V1
      RETURN
   10 U1 = 0.
      RETURN
      END
      SUBROUTINE ESHIFT (E,NC,NR,SHIFT)
C ----------------------------------------------------------------------
C    RE-SHIFT EIGENVALUES
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION E(NR,NC)
      DO 100 I=1,NR
      DO 100 K=1,NC
      E(I,K)=E(I,K)+SHIFT
 100  CONTINUE
      RETURN
      END
      SUBROUTINE EIG_JACOBI (A,B,X,EIGV,D,N,TOLJAC,SHIFT,IFSH,NSMAX)
C    +-----------------------------------------------------------+
C    -    P R O G R A M                                          -
C    -       TO SOLVE THE GENERAL EIGENPROBLEM                   -
C    -       USING THE GENERALIZED JACOBI ITERATION METHOD       -
C    -                                                           -
C--
C--
C+-----------------------------------------------------------+
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N),B(N,N),X(N,N),EIGV(N),D(N)
C      CHARACTER TEXT*20
C
C     INITIALIZE EIGENVALUE AND EIGENVEKTOR MATRICES
C
      DO 10 I=1,N
       IF (A(I,I).GT.0. AND .B(I,I).GT.0.) GOTO 4
C         CALL EIGMSG (009,'EIG_JACOBI',TEXT,0)
    4 D(I) = A(I,I)/B(I,I)
   10 EIGV(I)=D(I)
      DO 30 I=1,N
        DO 20 J=1,N
   20   X(I,J)=0.
   30 X(I,I)=1.
C
      IF (N.NE.1) GOTO 39
        RETURN
   39 CONTINUE
C
C    INITIALIZE SWEEP COUNTER AND BEGIN ITERATION
C
      NSWEEP=0
      NR = N - 1
   40 NSWEEP = NSWEEP + 1
C
C    CHECK IF PRESENT OFF DIAGONAL ELEMENT IS LARGE ENOUGH TO REQUIRE
C      ZEROING
C
      EPS = (0.01**NSWEEP)**2
      DO 210 J=1,NR
       JJ=J+1
       DO 210 K=JJ,N
        EPTOLA=(A(J,K)/A(J,J))*(A(J,K)/A(K,K))
        EPTOLB=(B(J,K)/B(J,J))*(B(J,K)/B(K,K))
        IF ((EPTOLA.LT.EPS).AND.(EPTOLB.LT.EPS)) GOTO 210
C
C      IF ZEROING IS REQUIRED CALCULATE THE ROTATION MATRIX
C       ELEMENTS CA AND CG
C
       AKK=A(K,K)*B(J,K)-B(K,K)*A(J,K)
       AJJ=A(J,J)*B(J,K)-B(J,J)*A(J,K)
       AB=A(J,J)*B(K,K)-A(K,K)*B(J,J)
       SCALE = A(K,K)*B(K,K)
       ABCH=AB/SCALE
       AKKCH=AKK/SCALE
       AJJCH=AJJ/SCALE
       CHECK=(ABCH*ABCH+4.*AKKCH*AJJCH)/4.
       IF (CHECK) 50,60,60
   50  CONTINUE
       STOP
   60  SQCH=SCALE*SQRT(CHECK)
       D1=AB/2.+SQCH
       D2=AB/2.-SQCH
       DEN=D1
       IF (ABS(D2).GT.ABS(D1)) DEN=D2
       IF (DEN) 80,70,80
   70  CA=0.
       CG=-A(J,K)/A(K,K)
       GOTO 90
   80  CA=AKK/DEN
       CG=-AJJ/DEN
C
C      PERFORM THE GENERALIZED ROTATION TO ZERO THE PRESENT OFF-DIAGONAL
C
   90  IF (N-2) 100,190,100
  100  JP1=J+1
       JM1=J-1
       KP1=K+1
       KM1=K-1
       IF (JM1-1) 130,110,110
  110  DO 120 I=1,JM1
       AJ=A(I,J)
       BJ=B(I,J)
       AK=A(I,K)
       BK=B(I,K)
       A(I,J)=AJ+CG*AK
       B(I,J)=BJ+CG*BK
       A(I,K)=AK+CA*AJ
  120  B(I,K)=BK+CA*BJ
  130  IF (KP1-N) 140,140,160
  140  DO 150 I=KP1,N
       AJ=A(J,I)
       BJ=B(J,I)
       AK=A(K,I)
       BK=B(K,I)
       A(J,I)=AJ+CG*AK
       B(J,I)=BJ+CG*BK
       A(K,I)=AK+CA*AJ
  150  B(K,I)=BK+CA*BJ
  160  IF (JP1-KM1) 170,170,190
  170  DO 180 I=JP1,KM1
       AJ=A(J,I)
       BJ=B(J,I)
       AK=A(I,K)
       BK=B(I,K)
       A(J,I)=AJ+CG*AK
       B(J,I)=BJ+CG*BK
       A(I,K)=AK+CA*AJ
  180  B(I,K)=BK+CA*BJ
  190  AK=A(K,K)
       BK=B(K,K)
       A(K,K)=AK+2.*CA*A(J,K)+CA*CA*A(J,J)
       B(K,K)=BK+2.*CA*B(J,K)+CA*CA*B(J,J)
       A(J,J)=A(J,J)+2.*CG*A(J,K)+CG*CG*AK
       B(J,J)=B(J,J)+2.*CG*B(J,K)+CG*CG*BK
       A(J,K)=0.
       B(J,K)=0.
C
C      UPDATE THE EIGENVECTOR MATRIX AFTER EACH ROTATION
C
       DO 200 I=1,N
        XJ=X(I,J)
        XK=X(I,K)
        X(I,J)=XJ+CG*XK
  200   X(I,K)=XK+CA*XJ
  210   CONTINUE
C
C      UPDATE THE EIGENVALUES AFTER EACH SWEEP
C
       DO 220 I=1,N
         IF (A(I,I).GT.0. .AND. B(I,I). GT.0.) GOTO 220
C          CALL EIGMSG(009,'EIG_JACOBI',TEXT,0)
  220  EIGV(I)=A(I,I)/B(I,I)
C
C      CHECK FOR CONVERGENCE
C
  230  DO 240 I=1,N
       TOL=TOLJAC*D(I)
       DIF=ABS(EIGV(I)-D(I))
        IF (DIF.GT.TOL) GOTO 280
  240  CONTINUE
C
C      CHECK ALL OFF DIAGONAL ELEMENTS TO SEE IF ANOTHER SWEEP
C        IS REQUIRED
C
       EPS=TOLJAC**2
       DO 250 J=1,NR
        JJ=J+1
        DO 250 K=JJ,N
        EPSA=(A(J,K)/A(J,J))*(A(J,K)/A(K,K))
        EPSB=(B(J,K)/B(J,J))*(B(J,K)/B(K,K))
        IF( (EPSA.LT.EPS).AND.(EPSB.LT.EPS)) GOTO 250
       GOTO 280
  250  CONTINUE
C
C      FILL OUT BOTTOM TRIANGLE OF RESULTANT MATRICES AND SCALE EIGENVECT.
C
  255 DO 260 I=1,N
        DO 260 J=1,N
        A(J,I)=A(I,J)
  260 B(J,I)=B(I,J)
        DO 270 J=1,N
        BB=SQRT(B(J,J))
         DO 270 K=1,N
  270   X(K,J)=X(K,J)/BB
       RETURN
C
C      UPDATE D MATRIX AND START NEW SWEEP IF ALLOWED
C
  280  DO 290 I=1,N
  290  D(I)=EIGV(I)
       IF (NSWEEP.LT.NSMAX) GOTO 40
       GOTO 255
C
C
       END
