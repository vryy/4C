c
c     Maintainer: Andreas Lipka 
c                 lipka@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/lipka/ 
c                 0711 - 685-6575 
c
c     ---------------------------------------------------------------  
C-----------------------------------------------------------------------
      SUBROUTINE  FSDOC  (VAR,DF,DG,ETAI,ETHA,XDGO,RESU,RESL,VARUP,
     *                    VARLO,NUMVAR,BETA,ACCIT,DELTA,IPRINT )
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER   INFO
      DIMENSION VAR(NUMVAR),DF(NUMVAR),DG(NUMVAR),ETAI(NUMVAR)
      DIMENSION RESU(NUMVAR),RESL(NUMVAR),VARUP(NUMVAR),VARLO(NUMVAR)
C------------------------------------------------------ IN LINE FUNCTION
      POW(A,B) = EXP(B*LOG(A))
C
C-----------------------------------------------------------------------
      INFO=12
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      HALF=0.5D0
      PI=3.1415927D0
      PIHF=1.5707963D0
      UP=1.0D+08
      UL=1.0D-08
C
C------------------------------------------- OPEN OUPUTFILE
      IF (IPRINT.GT.0) THEN
      OPEN(INFO, FILE='ocalg.doc')
      ENDIF
C------------------------------------------- GET OPTIMIZATION PARAMETERS
       CONACT= 0.0
C------------------------------------------------------- INITIALIZE DATA
      ITRMAX=NUMVAR*10
C
C**********************Kontrollausgabe*********************************
      IF (IPRINT.GT.0) THEN
        WRITE(INFO,*)
        WRITE(INFO,*)'VAR,DF,DG,ETAI'
      DO 10 N=1,NUMVAR
        WRITE(INFO,'(4e14.7)')VAR(N),DF(N),DG(N),ETAI(N)
 10   CONTINUE
        WRITE(INFO,*)
        WRITE(INFO,*)'RESU,RESL,VARUP,VARLO'
      DO 20 N=1,NUMVAR
        WRITE(INFO,*)RESU(N),RESL(N),VARUP(N),VARLO(N)
 20   CONTINUE
      ENDIF
C **********************Kontrollausgabe*********************************
C
C------------------------------------ GET CURRENT UPPER AND LOWER BOUNDS
      DO NVAR=1,NUMVAR            
        V1 = VAR(NVAR)*(ONE-DELTA)
        V2 = VAR(NVAR)*(ONE+DELTA)
        VARUP(NVAR) = MIN(RESU(NVAR),MAX(V1,V2))
        VARLO(NVAR) = MAX(RESL(NVAR),MIN(V1,V2))
      ENDDO
C                                                                       
C------------------------------------------------- SHIFT UPDATE GRADIENT 
      SCLA=ZERO
C      
      IF(ETHA.NE.ZERO) THEN
        XUDF = ZERO
        XLDF = ZERO
C
        IF(XDGO.LT.UL) STOP "FSDOC1: ERROR: IOTA = 0 !"
C   
        DO NVAR=1,NUMVAR
          DGX  = -DG(NVAR)
          IF (ABS(DGX).GT.UL) THEN
            IF(ETHA.GT.ZERO) THEN
              DGF  = XDGO * (ETHA + (ONE-ETHA)*(DGX/XDGO))
            ELSE
              DGF  = DGX
            ENDIF
            FAC  = DF(NVAR)/DGF
            XUDF = MAX(FAC,XUDF)
            XLDF = MIN(FAC,XLDF)
          ENDIF
        ENDDO
C      
        IF(XUDF.GT.ZERO) THEN
C        
          SCLA= XUDF + UL 
C
      IF (IPRINT.GT.0) THEN
          WRITE(INFO,*) "" 
          WRITE(INFO,*) "FSDOC1: MODIFIKATION OF GRADIENTS IN FSDOC1"
          WRITE(INFO,*) "        XUDF:",XUDF
          WRITE(INFO,*) "        XLDF:",XLDF
      ENDIF
C
        ENDIF
C
      ENDIF             
C                                                                       
C                                                                       
C------------------------------------------- DETERMINE GRADIENT QUOTIENT 
      DO NVAR=1,NUMVAR
C
        IF (ETHA.NE.ZERO) THEN
C
          DGX  = -DG(NVAR)
          IF(ABS(DGX).GT.ZERO) THEN
            IF(ETHA.GT.ZERO) THEN
              DGF  = XDGO * (ETHA + (ONE-ETHA)*(DGX/XDGO))
            ELSE
              DGF  = DGX
            ENDIF
            FAC  = SCLA - DF(NVAR)/DGF
            IF(FAC.LT.ZERO) THEN
            IF (IPRINT.GT.0) THEN
              WRITE(INFO,*) "FSDOC1: ETHA > 0 , UPDATE RULE < 0"
            ENDIF
              FAC = ZERO
            ENDIF  
          ELSE
            IF(DF(NVAR).GT.ZERO) THEN
              FAC = VARLO(NVAR)/VAR(NVAR)
            ELSEIF(DF(NVAR).LT.ZERO) THEN
              FAC = VARUP(NVAR)/VAR(NVAR)
            ELSE
              FAC = ONE
            ENDIF
          ENDIF
C
        ELSE 
C
          IF(ABS(DG(NVAR)).GT.UL) THEN
            FAC = DF(NVAR) / DG(NVAR) 
            IF(FAC.LT.ZERO) THEN
            IF (IPRINT.GT.0) THEN
              WRITE(INFO,*) "FSDOC1: ETHA = 0 , UPDATE RULE < 0"
            ENDIF
              FAC = ZERO
            ENDIF
          ELSE
            IF(DF(NVAR).GT.ZERO) THEN
              FAC = VARLO(NVAR)/VAR(NVAR)
            ELSEIF(DF(NVAR).LT.ZERO) THEN
              FAC = VARUP(NVAR)/VAR(NVAR)
            ELSE
              FAC = ONE
            ENDIF  
          ENDIF
C
        ENDIF
C
        ETAI(NVAR) =  POW(FAC,BETA)
C
      ENDDO
C
C--------------------------------------------- DETERMINE MIN/MAX OF RLAM      
      RLAMIN =  ONE/UL
      RLAMAX = -ONE/UL
C      
      DO NVAR=1,NUMVAR
        FAC    = VAR(NVAR)*ETAI(NVAR)
        IF(ABS(FAC).GT.UL) THEN
          FAC1   = VARUP(NVAR)/FAC
          FAC2   = VARLO(NVAR)/FAC
        ELSE
          FAC1   = UP
          FAC2   = UP
        ENDIF
        RLAMIN = MIN(RLAMIN,FAC1,FAC2)
        RLAMAX = MAX(RLAMAX,FAC1,FAC2)
      ENDDO
C                     
C--------------- CHECK IF CONSTRAINT CAN BE SATISFIED BY RLAMIN / RLAMAX
      SRMIN = ZERO
      SRMAX = ZERO
C      
      DO NVAR=1,NUMVAR
        IF(ABS(DG(NVAR)).GT.UL) THEN
          X = RLAMIN*VAR(NVAR)*ETAI(NVAR)
          X = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
          SRMIN = SRMIN + DG(NVAR)*(X-VAR(NVAR))
          X = RLAMAX*VAR(NVAR)*ETAI(NVAR)
          X = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
          SRMAX = SRMAX + DG(NVAR)*(X-VAR(NVAR))
        ENDIF
      ENDDO
C      
      SRMIN = SRMIN + CONACT 
      SRMAX = SRMAX + CONACT
      SRQ   = SRMIN/SRMAX
C      
      IF(SRQ.GT.ZERO) THEN 
        WRITE(*,*) 'CONSTRAINTS CAN NOT BE SATISFIED'
        IF(SRQ.LT.ONE) THEN
          WRITE(*,*) '- ALL VARIABLES SET ON CURRENT LOWER BOUNDS'
          RLAM = RLAMIN
          GOTO 900
        ELSE
          WRITE(*,*) '- ALL VARIABLES SET ON CURRENT UPPER BOUNDS'
          RLAM = RLAMAX
          GOTO 900
        ENDIF
      ENDIF  
C
C-------------------------------------- INITIALIZE LOOP TO DETERMIN RLAM 
      ITR   = 0
      SRDIF = TWO*ACCIT
C
      DO WHILE (ABS(SRDIF).GT.ACCIT .AND. ITR.LT.ITRMAX)
        ITR   = ITR + 1
        SR    = ZERO
        RLAM  = HALF*(RLAMIN+RLAMAX)
        DO NVAR=1,NUMVAR
          IF(ABS(DG(NVAR)).GT.UL) THEN 
            X  = RLAM*VAR(NVAR)*ETAI(NVAR)
            X  = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
            SR = SR + DG(NVAR)*(X-VAR(NVAR))
          ENDIF 
        ENDDO
        SRDIF = SR + CONACT
        IF((SRDIF*SRMIN).GE.0 .AND. (SRDIF*SRMAX).LE.0 ) THEN
          RLAMIN = RLAM
          SRMIN  = SRDIF
        ELSEIF((SRDIF*SRMAX).GE.0 .AND. (SRDIF*SRMIN).LE.0 ) THEN
          RLAMAX = RLAM
          SRMAX  = SRDIF
        ELSE
          WRITE(*,*) 'ERROR IN FSDOC2: RLAM CAN NOT BE DETERMINED'
          STOP
        ENDIF
      ENDDO
C
C----------------------------------------------------- CHECK CONVERGENCE
      IF(ITR.GE.ITRMAX) THEN
        WRITE(*,*) 'CONSTRAINTS CAN NOT BE SATISFIED - MAXIMUM NUMBE'
     *               ,'R OF ITERATION REACHED'
        WRITE(*,*) 'ERROR IN CONSTRAINT :',SRDIF
      ENDIF
C
C---------------------------------------------------------- FINAL UPDATE  
  900 CONTINUE
      DO NVAR=1,NUMVAR
        IF(ABS(DG(NVAR)).GT.UL) THEN 
          X         = RLAM*VAR(NVAR)*ETAI(NVAR)
          VAR(NVAR) = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
        ELSE
          X         = VAR(NVAR)*ETAI(NVAR)
          VAR(NVAR) = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
        ENDIF 
      ENDDO
C
C------------------------------------------------------------------- END
C      
      RETURN
      END

