c
c     Maintainer: Malte Neumann 
c                 neumann@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/neumann/ 
c                 0711 - 685-6121 
c
c     ---------------------------------------------------------------  
C***********************************************************************
C   File:  @(#)solver_column_height.ftn	1.1   01/09/96   
C***********************************************************************
      SUBROUTINE COLSOL (A,V,MAXA,NN,NRR,NRC,NWA,NQM,NR1,NR2,KKK,
     *                   DET,ISC,NSCH,IPR,INFO)
C    +-----------------------------------------------------------------+
C    |  PROGRAM TO SOLVE  LINEAR EQUATION SYSTEM USING  COLUMN-HEIGHT  |
C    |  TECHNIQUE AND ONE-DIMENSIONAL STORED SYSTEM MATRIX.            |
C    +-----------------------------------------------------------------+
C    |  PARAMETERS :                                                   |
C    |                                                                 |
C    |  A      = SYSTEM MATRIX (ONE-DIMENSIONAL STORAGE)          (I)  |
C    |  V      = RIGHT-HAND-SIDE VECTOR AND SOLUTION              (I/O)|
C    |  MAXA   = ADDRESSES OF DIAGONAL ELEMENTS OF A              (I)  |
C    |  NN     = NUMBER OF EQUATIONS    (= NEQ)                   (I)  |
C    |  NRR    = NUMBER OF ROWS    OF RHS-ARRAY (= NEQRHS)        (I)  |
C    |  NRC    = NUMBER OF COLUMNS OF RHS-ARRAY (= NUMRHS)        (I)  |
C    |  NWA    = LENGTH OF ONE-DIMENSIONAL SYSTEM MATRIX          (I)  |
C    |  NQM    = LENGTH OF VECTOR  MAXA (= NEQ1)                  (I)  |
C    |  NR1    = FIRST VECTOR FOR SOLUTION                        (I)  |
C    |  NR2    = LAST  VECTOR FOR SOLUTION                        (I)  |
C    |  KKK    = FLAG TO CONTROL EXECUTION                        (I)  |
C    |           1 = DECOMPOSITION OF SYSTEM MATRIX                    |
C    |           2 = BACKSUBSTITUTION                                  |
C    |           3 = DECOMPOSITION AND BACKSUBSTITUTION                |
C    |  DET    = DETERMINANT STIFFNESS MATRIX                      (O) |
C    |  ISC    = NUMBER OF SCALINGS FOR DETERMINANT                (O) |
C    |  NSCH   = NUMBER OF NEGATIVE DIAGONAL ELEMENTS DURING FACT. (O) |
C    |  IPR    = TRACE PRINT PARAMETER PROCESSING EQUATIONS        (I) |
C    |  INFO   = UNIT NUMBER TO PRINT ERROR/INFORMATION MESSAGES   (I) |
C    +-----------------------------------------------------------------+
c      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT NONE
      REAL*8 ZERO,ONE
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c the fortran90 parameter statement causes the dde debugger to
c crash the globally defined structures, do not use!!!! m.gee 01/02
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER NN,NRR,NRC,NWA,NQM,NR1,NR2,KKK,ISC,NSCH,IPR
      INTEGER INFO,ICOMP,IRED,IERR,N,NDI,INC,KN,KL,KU,KH,K,IC
      INTEGER KLT,J,KI,ND,KK,L,KIL,KLTL,M,NRA
      REAL*8  DET,C,B,TT,T1,T2,T3,T4
      REAL*8 A(NWA),V(NRR,NRC)
      INTEGER MAXA(NQM)
      INTEGER PRTDEB
      DATA IRED/0/, ICOMP/0/ 
      SAVE ICOMP,IRED
C
C What's the use of these variables?
C
      T1=0.0
      T2=0.0
      T3=0.0
      T4=0.0
C
C
C       IPR = 0    ... NO TRACE PRINT DURING FACTORIZATION
C       IPR = NNN  ... PRINT MESSAGE AFTER PROCESSING NNN EQUATIONS
C
      zero=0.0
      one=1.0

      IF (IPR.GT.0) THEN
        PRTDEB=1
      ELSE
        PRTDEB=0
      ENDIF
C
      IF (KKK.EQ.1 .OR. KKK.EQ.3) THEN
C------------------------------------------------------------------
C        CHECK FOR ZERO DIAGONAL ELEMENTS ONLY FOR THE FIRST
C        CALL DURING ONE RUN OF PROGRAM
C------------------------------------------------------------------
        IF (ICOMP.EQ.0) THEN
         IERR=0
         IF (PRTDEB.GT.0) WRITE (INFO,'(/A)') 'DIAGONAL CHECK - START'
         DO 20 N=1,NN
           NDI=MAXA(N)
           IF (A(NDI).EQ.ZERO) THEN
             IF (IERR.EQ.0) THEN
               WRITE (INFO,'(A)') '... ZERO DIAGONAL ELEMENTS FOUND AT'
               IERR=1
             ENDIF
             WRITE (INFO,'(A,I10)') '... EQUATION',N
           ENDIF
   20    CONTINUE
         IF (PRTDEB.GT.0) WRITE (INFO,'(A)') 'DIAGONAL CHECK - END'
         IF (IERR.EQ.1) THEN
           WRITE (INFO,*) '... SINGULARITY DETECTED ON SYSTEM MATRIX'
           WRITE (INFO,*) '... PROGRAM STOP !!'
           STOP
         ENDIF
        ENDIF
C---------------------------------------------------------------------
        IF (PRTDEB.GT.0) THEN
          WRITE (INFO,'(/A,I6,A/)') 'PERFORM FACTORIZATION FOR',NN,
     +                           ' EQUATIONS'
          INC = IPR
        ENDIF
C-------------------------------------------- INITIALIZE CP-TIME COUNTER
C----------------------------------- INITIALIZE FACTORIZATION PARAMETERS
        DET=ONE
        ISC=0
        NSCH=0
C------------------------------- L*D*L(T) DECOMPOSITION OF SYSTEM MATRIX
        DO 140 N=1,NN
          KN=MAXA(N)
          KL=KN+1
          KU=MAXA(N+1)-1
          KH=KU-KL
          IF (KH) 110,90,50
   50     CONTINUE
          K=N-KH
          IC=0
          KLT=KU
          DO 80 J=1,KH
            IC=IC+1
            KLT=KLT-1
            KI=MAXA(K)
            ND=MAXA(K+1)-KI-1
            IF (ND.GT.0) THEN
              KK=MIN0(IC,ND)
              C=ZERO
              DO 70 L=1,KK
                KIL=KI+L
                KLTL=KLT+L
                C=C+A(KIL)*A(KLTL)
   70         CONTINUE
              A(KLT)=A(KLT)-C
             ENDIF
             K=K+1
   80     CONTINUE
   90     CONTINUE
          K=N
          B=ZERO
          DO 100 KK=KL,KU
            K=K-1
            KI=MAXA(K)
            C=A(KK)/A(KI)
            B=B+C*A(KK)
            A(KK)=C
  100     CONTINUE
          A(KN)=A(KN)-B
  110     CONTINUE
          CALL SOLSCL (DET,ISC,NSCH,A(KN))
          IF (PRTDEB.GT.0) THEN
            M=N/INC
            IF (M.GT.0) THEN
              WRITE (6,'(I6,A)') INC,' EQUATIONS PROCESSED'
              INC=INC+IPR
            ENDIF
          ENDIF
  140   CONTINUE
      ENDIF
C***********************************************************************
C***              REDUCTION OF RIGHT-HAND-SIDE AND BACKSUBSTITUTION  ***
C***********************************************************************
      IF (KKK.EQ.2 .OR. KKK.EQ.3) THEN
        IF (IRED.GT.0) PRTDEB=0
        DO 240 NRA = NR1,NR2
          IF (PRTDEB.GT.0) THEN
            WRITE (6,'(/A,I6,A)') 'PERFORM FORWARD REDUCTION FOR',NN,
     *                             ' EQUATIONS'
          ENDIF
          DO 180 N=1,NN
             KL=MAXA(N)+1
             KU=MAXA(N+1)-1
             IF (KU-KL.LT.0) GOTO 180
             K=N
             C=ZERO
             DO 170 KK=KL,KU
                K=K-1
                C=C+A(KK)*V(K,NRA)
  170        CONTINUE
             V(N,NRA)=V(N,NRA)-C
  180     CONTINUE
C
          DO 200 N=1,NN
             K=MAXA(N)
             V(N,NRA)=V(N,NRA)/A(K)
  200     CONTINUE
C----------------------------------
          IF (NN.EQ.1) GO TO 900
C----------------------------------
          IF (PRTDEB.GT.0) THEN
            WRITE (6,'(/A,I6,A)') 'PERFORM BACKSUBSTITUTION FOR ',NN,
     *                           ' EQUATIONS'
          ENDIF
          N=NN
          DO 230 L=2,NN
             KL=MAXA(N)+1
             KU=MAXA(N+1)-1
             IF (KU-KL.GE.0) THEN
                 K=N
                 DO 220 KK=KL,KU
                    K=K-1
                    V(K,NRA)=V(K,NRA)-A(KK)*V(N,NRA)
  220            CONTINUE
             ENDIF
             N=N-1
  230     CONTINUE
  240   CONTINUE
      ENDIF
C
  900 CONTINUE
C-----------------------------------------------------------------
      if (kkk.eq.1 .or. kkk.eq.3) then
       if (icomp.eq.0) then
        icomp=1
        tt = t2-t1
       endif
      endif
c
      if (kkk.eq.2 .or. kkk.eq.3) then
       if (ired.eq.0) then
        ired=1
        tt = t4-t3
       endif
      endif
c------------------------------------------------------------------
C
      RETURN
      END
C===================================================================
      SUBROUTINE SOLSCL (DET,ISC,NSCH,AD)
C
c      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none
C
C      PROGRAM TO SCALE DETERMINANT
C      MAX REAL NUMBER  1.0E+38
C
      real*8 det,ad,dmin,dmax,zero
      integer isc,nsch
c      PARAMETER (DMIN=1.0D-25,DMAX=1.0D+25,ZERO=0.0D0)
      dmin = 1.0D-25
      dmax = 1.0D+25
      zero = 0.0
C
      IF (ABS(DET).LE.DMIN) THEN
        DET=DET*DMAX
        ISC=ISC-1
      ELSE IF (ABS(DET).GE.DMAX) THEN
        DET=DET/DMAX
        ISC=ISC+1
      ENDIF
C--------------------
      IF (AD.LT.ZERO) NSCH=NSCH+1
      DET=DET*AD
C--------------------
      RETURN
      END
