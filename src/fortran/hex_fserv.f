c
c     Maintainer: Andreas Lipka 
c                 lipka@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/lipka/ 
c                 0711 - 685-6575 
c
c     ---------------------------------------------------------------  
C=======================================================================
      SUBROUTINE fortranpow (V,R,RE)
C-----------------------------------------------------------------------        
C     R = V^RE
C-----------------------------------------------------------------------        
      IMPLICIT REAL*8 (A-H,O-Z)
C
       R=V**RE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MXMAB (A,B,R,NZA,NSA,NSB)
C-----------------------------------------------------------------------        
C     R(I,J) = A(I,K)*B(K,J) ---  R = A*B
C-----------------------------------------------------------------------        
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),R(1)
C
       NN=NZA*NSB
       DO 5 I=1,NN
    5  R(I)=0.0
C
      IC0=0
      IB=0
      DO 30 I=1,NSB
      IA=0
      DO 20 J=1,NSA
      IB = IB+1
      S=B(IB)
      IC=IC0
      DO 10 K=1,NZA
      IC=IC+1
      IA=IA+1
   10 R(IC) = R(IC) + A(IA)*S
   20 CONTINUE
   30 IC0=IC0+NZA
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MXMATB (A,B,R,NZA,NSA,NSB,NULL)
C-----------------------------------------------------------------------        
C     R(I,J) = A(K,I)*B(K,J)  --- R = A(TRANS)*B
C-----------------------------------------------------------------------        
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1),R(1)
C
      IF (NULL.EQ.0) THEN
       NN=NSA*NSB
       DO 5 I=1,NN
    5  R(I)=0.0
      ENDIF
C
      IB0=0
      IC=0
      DO 30 I=1,NSB
      IA=0
      DO 20 K=1,NSA
      IC=IC+1
      S=0.0
      IB=IB0
      DO 10 J=1,NZA
      IA=IA+1
      IB=IB+1
   10 S=S + A(IA)*B(IB)
   20 R(IC) = R(IC) + S
   30 IB0=IB0+NZA
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MXMABT (A,M,N,L,B,R)
C-----------------------------------------------------------------------        
C     MATRIX MULTIPLICATION  A*B(TRANS)
C-----------------------------------------------------------------------        
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(1),B(1),R(1)
C
      DO 1 I=1,M
      IK=I-M
      DO 1 J=1,L
      HELP=0.0
      IA=I-M
      IB=J-L
      DO 2 K=1,N
      IA=IA+M
      IB=IB+L
    2 HELP=HELP+A(IA)*B(IB)
      IK=IK+M
      R(IK)=HELP
    1 CONTINUE
C
      RETURN
      END
C======================================================================
C
      SUBROUTINE SOLVEQ (A,B,NN,LL)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN,NN),B(NN,LL),ID(100)
C   
      DO 50 N=1,NN
 50   ID(N)=N
      DO 475 N=1,NN
      N1=N+1
      D=0.0D0
      DO 100 I=N,NN
      DO 100 J=N,NN
      IF(DABS(A(I,J))-D) 99,90,90
 90   D=DABS(A(I,J))
      II=I
      JJ=J
 99   CONTINUE
 100  CONTINUE
      DO 110 I=1,NN
      D=A(I,N)
      A(I,N) = A(I,JJ)
 110  A(I,JJ)= D
      I=ID(N)
      ID(N) = ID(JJ)
      ID(JJ)=I
      DO 120 J=N,NN
      D=A(N,J)
      A(N,J)=A(II,J)
 120  A(II,J)=D
      DO 130 L=1,LL
      D=B(N,L)
      B(N,L) = B(II,L)
 130  B(II,L)= D
      DO 150 L=1,LL
 150  B(N,L) = B(N,L)/A(N,N)
      IF(N-NN) 200,500,200
 200  CONTINUE
      DO 450 J=N1,NN
      IF(A(N,J)) 250,350,250
 250  A(N,J) = A(N,J)/A(N,N)
      DO 300 I=N1,NN
 300  A(I,J) = A(I,J)-A(I,N)*A(N,J)
 350  CONTINUE
      DO 400 L=1,LL
 400  B(J,L) = B(J,L)-A(J,N)*B(N,L)
 450  CONTINUE
 475  CONTINUE
 500  N1=N
      N=N-1
      IF(N) 700,700,550
 550  CONTINUE
      DO 600 L=1,LL
      DO 600 J=N1,NN
 600  B(N,L) = B(N,L)-A(N,J)*B(J,L)
      GO TO 500
 700  CONTINUE
      DO 950 N=1,NN
      DO 900 I=N,NN
      IF(ID(I)-N) 900,750,900
 750  CONTINUE
      DO 800 L=1,LL
      D=B(N,L)
      B(N,L) = B(I,L)
 800  B(I,L)=D
      GO TO 950
 900  CONTINUE
 950  ID(I) = ID(N)
      RETURN
      END
C
C======================================================================
C======================================================================
      SUBROUTINE C1JACB (A,V)
C     ******************************************************************
C     *                                                                *
C     *                      +---------------+                         *
C     *                      I  C 1 J A C B  I                         *
C     *                      +---------------+                         *
C     *                                                                *
C     * PROGRAM TO CALCULATE THE EIGENVALUES AND EIGENVECTORS BY       *
C     * THE JACOBI METHOD                                              *
C     ******************************************************************
C     * PARAMETER LIST:                                                *
C     * A      --> SYSTEM MATRIX (FOR EX. STRESS MATRIX)         (I/O) *
C     *            AFTER COMPUTATION IN DIAGONAL FORM (EIGENVALUES)    *
C     * V      --> PRINCIPAL DIRECTIONS (EIGENVECTORS)           (O)   *
C     ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      PARAMETER (ITM=200,ERR=1.D-9,ZERO=0.0D0,ONE=1.0D0,ZTOLP=1.D-9,
C     +TWO=2.0D0,FOUR=4.0D0,KKK=3)
      DIMENSION A(3,3),V(3,3)
C
      IT=0
C
C-----------------------------------------------------------------------
      ITM=200
      ERR=1.D-9
      ZERO=0.0D0
      ONE=1.0D0
      ZTOLP=1.D-9
      TWO=2.0D0
      FOUR=4.0D0
      KKK=3
C-------------------------------------------PUT A UNIT MATRIX IN ARRAY V
C
      DO 10 I=1,KKK
        V(I,I)=ONE
 10   CONTINUE
C
C---------------------------------FIND LARGEST OFF DIAGONAL COEFFICIENTS
C
 13   T=ZERO
      M=2
      DO 20 I=1,M
        J1=I+1
        DO 20 J=J1,KKK
        IF(ABS(A(I,J))-T) 20,20,2
 2      T=ABS(A(I,J))
        IR=I
        IC=J
 20   CONTINUE
      IF (IT) 5,4,5
C
C------------------------------------TAKE FIRST OFF DIAGONAL COEFFICIENT
C---------------------------------TIMES ERR AS COMPARISON VALUE FOR ZERO
C
 4    T1=T*ERR
 5    IF(T-T1) 999,999,6
C
C-------COMPUTE TAN(TA),SIN(S),AND COS(C) OF ROTATION ANGLE
C
 6    PS=A(IR,IR)-A(IC,IC)
      TA=(-PS+SQRT(PS*PS+FOUR*T*T))/(TWO*A(IR,IC))
      C=ONE/SQRT(ONE+TA*TA)
      S=C*TA
C
C-----------------------MULTIPLY ROTATION MATRIX TIMES V AND STORE IN V
C
      DO 50 I=1,KKK
        P=V(I,IR)
        V(I,IR)=C*P+S*V(I,IC)
        V(I,IC)=C*V(I,IC)-S*P
 50   CONTINUE
      I=1
 100  IF(I-IR) 7,200,7
C
C-------------------APPLY ORTHOGONAL TRANSFORMATION TO A AND STORE IN A
C
 7    P=A(I,IR)
      A(I,IR)=C*P+S*A(I,IC)
      A(I,IC)=C*A(I,IC)-S*P
      I=I+1
      GO TO 100
 200  I=IR+1
 300  IF(I-IC) 8,400,8
 8    P=A(IR,I)
      A(IR,I)=C*P+S*A(I,IC)
      A(I,IC)=C*A(I,IC)-S*P
      I=I+1
      GO TO 300
 400  I=IC+1
 500  IF(I-KKK) 9,9,600
 9    P=A(IR,I)
      A(IR,I)=C*P+S*A(IC,I)
      A(IC,I)=C*A(IC,I)-S*P
      I=I+1
      GO TO 500
 600  P=A(IR,IR)
      A(IR,IR)=C*C*P+TWO*C*S*A(IR,IC)+S*S*A(IC,IC)
      A(IC,IC)=C*C*A(IC,IC)+S*S*P-TWO*C*S*A(IR,IC)
      A(IR,IC)=ZERO
      IT=IT+1
      IF(IT-ITM) 13,13,999
C
 999  RETURN
      END  
C=======================================================================
      SUBROUTINE C1INV6(A,B,IRC)
C----------------------------------------------------------------------C
C!    INVERTIERUNG EINER 6x6 MATRIX
C         A   -->  INPUT
C         B   -->  INVERSE VON A 
C         IRC -->  RETURN CODE ( IRC = 0   OK
C                                IRC = 1   ERROR: MATRIX SINGULAER )  
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
C      PARAMETER (TOL=1.0D-20,N=6)
      DIMENSION A(6,6),B(6,6),MX(6),MY(6)

      TOL=1.0D-20
      N=6

      IRC=0                 

C     UMSPEICHERN DER MATRIX A AUF B
      DO 10 I=1,N
         DO 20 J=1,N
            B(I,J)=A(I,J)
   20    CONTINUE
   10 CONTINUE

C     VORBESETZEN DER PIVOT-VEKTOREN MX UND MY MIT NULL

      DO 30 I=1,N
         MX(I)=0
         MY(I)=0
   30 CONTINUE

C     BESTIMMUNG DES PIVOTELEMENTES

      DO 40 I=1,N
         PIVO=0.0
         DO 50 IX=1,N
            IF(MX(IX).EQ.0) THEN
               DO 60 IY=1,N
                  IF(MY(IY).EQ.0) THEN
                     IF(DABS(B(IX,IY)).GT.DABS(PIVO)) THEN
                        PIVO=B(IX,IY)
                        NX=IX
                        NY=IY
                     END IF
                  END IF      
   60          CONTINUE
            END IF
   50    CONTINUE
C        FALLS DAS PIVOTELEMENT NULL IST,IST DIE MATRIX SINGULAER
         IF(DABS(PIVO).LT.TOL) THEN
            IRC=1
            RETURN
         END IF
C        MERKEN DER INDIZES DES PIVOTELEMENTES
         MX(NX)=NY
         MY(NY)=NX
C        BERECHNUNG DER MATRIXELEMENTE
         HILF=1.0/PIVO
         DO 70 J=1,N
            IF(J.NE.NX) THEN
               FAKTOR=B(J,NY)*HILF
               DO 80 K=1,N
                  B(J,K)=B(J,K)-B(NX,K)*FAKTOR
   80          CONTINUE
               B(J,NY)=-FAKTOR
            END IF
   70    CONTINUE
         DO 90 K=1,N
            B(NX,K)=B(NX,K)*HILF
   90    CONTINUE
         B(NX,NY)=HILF
   40 CONTINUE

C     ZEILEN UND SPALTENVERTAUSCHUNG RUECKGAENGIG MACHEN

      DO 100 I=1,N-1
         DO 110 M=I,N
            IF(MX(M).EQ.I) GOTO 120
  110    CONTINUE
  120    J=M
         IF(J.NE.I) THEN
            DO 130 K=1,N
               H=B(I,K)
               B(I,K)=B(J,K)
               B(J,K)=H
  130       CONTINUE
            MX(J)=MX(I)
            MX(I)=I
         END IF
         DO 140 M=I,N
            IF(MY(M).EQ.I) GO TO 150
  140    CONTINUE
  150    J=M
         IF(J.NE.I) THEN
            DO 160 K=1,N
               H=B(K,I)
               B(K,I)=B(K,J)
               B(K,J)=H
  160       CONTINUE
            MY(J)=MY(I)
            MY(I)=I
         END IF
  100 CONTINUE

C     FEHLERBERECHNUNG
      S1=0.0
      S2=0.0
      DO 170 I=1,N
         DO 180 J=1,N
            H=0.0
            DO 190 K=1,N
               H=H+A(I,K)*B(K,J)
  190       CONTINUE
            IF(I.EQ.J) THEN
               S1=S1+DABS(H-1.0)
            ELSE
               S2=S2+DABS(H)
            END IF
  180    CONTINUE  
  170 CONTINUE
      
      RETURN
      END
C***********************************************************************
C   File:  @(#)s9_math_lib.ftn	1.4   02/02/01   
C***********************************************************************
      SUBROUTINE C1AB (A,B,R,NZA,NSA,NSB,NULL)                                  
C                                                                               
C       R(I,J) = A(I,K)*B(K,J) ---  R = A*B                                     
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                                 
      DIMENSION A(1),B(1),R(1)                                                  
C                                                                               
      IF (NULL.EQ.0) THEN                                                       
       NN=NZA*NSB                                                               
       DO 5 I=1,NN                                                              
    5  R(I)=0.0                                                                 
      ENDIF                                                                     
C                                                                               
      IC0=0                                                                     
      IB=0                                                                      
      DO 30 I=1,NSB                                                             
      IA=0                                                                      
      DO 20 J=1,NSA                                                             
      IB = IB+1                                                                 
      S=B(IB)                                                                   
      IC=IC0                                                                    
      DO 10 K=1,NZA                                                             
      IC=IC+1                                                                   
      IA=IA+1                                                                   
   10 R(IC) = R(IC) + A(IA)*S                                                   
   20 CONTINUE                                                                  
   30 IC0=IC0+NZA                                                               
C                                                                               
      RETURN                                                                    
      END                                                                       
C=======================================================================
      SUBROUTINE C1INVF (FN,FNI,DETF)
C----------------------------------------------------------------------
C     TOPIC : INVERSION OF 3X3 MATRIX 
C     FN  = 3 X 3 MATRIX IN VECTOR FORM TO BE INVERTED
C     FNI = INVERSE FROM FN
C     DETF= DETERMINANT FROM FNI 
C----------------------------------------------------------------------
      IMPLICIT  REAL*8 (A-H,O-Z)
C      PARAMETER (NDISD=9)
C      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
      DIMENSION  FN(9),FNI(9)

      NDISD=9
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      
      DUM = ( FN(1)*FN(2)*FN(3) +
     1        FN(5)*FN(7)*FN(8) +
     2        FN(9)*FN(4)*FN(6) -
     3        FN(8)*FN(2)*FN(9) -
     4        FN(6)*FN(7)*FN(1) -
     5        FN(3)*FN(4)*FN(5) )
      DETF=ONE/DUM
C  INVERSE FROM FN
      FNI(1) = DETF*(FN(2)*FN(3) - FN(7)*FN(6))  
      FNI(2) = DETF*(FN(1)*FN(3) - FN(8)*FN(9))  
      FNI(3) = DETF*(FN(1)*FN(2) - FN(5)*FN(4))  
      FNI(5) = DETF*(FN(9)*FN(6) - FN(5)*FN(3))  
      FNI(4) = DETF*(FN(7)*FN(8) - FN(4)*FN(3))  
      FNI(7) = DETF*(FN(9)*FN(4) - FN(1)*FN(7))  
      FNI(6) = DETF*(FN(5)*FN(8) - FN(1)*FN(6))  
      FNI(9) = DETF*(FN(5)*FN(7) - FN(9)*FN(2))  
      FNI(8) = DETF*(FN(4)*FN(6) - FN(2)*FN(8))  
C
C
      RETURN
      END
C
C=======================================================================
      SUBROUTINE C1PRDR (TENSOR,T1,T2,T3,FACTOR)
C     ******************************************************************
C     *                                                                *
C     *                       C 1 P R D R                              *
C     *                                                                *
C     * SUBROUTINE TO CALCULATE THE EIGENVALUES (T1,T2,T3) OF A SECOND *
C     * ORDER TENSOR "TENSOR"                                          *
C     ******************************************************************
C     *   PARAMETER LIST :                                             *
C     *   TENSOR   ===> TENSOR                                    (I/O)*
C     *   T1,T2,T3 ===> PRINCIPAL VALUES                           (I) *
C     *   FACTOR   ===> FACTOR THAT CLASSIFIES A SPECIFIC STATE    (I) *
C     ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER SNAME*6
      DIMENSION TENSOR(6),DEVIAT(6)
C      PARAMETER (SNAME='C1PRDR',ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
C      PARAMETER (THREE=3.0D0,TOL=1.E-5)
C
      SNAME='C1PRDR'
      ZERO=0.0D0
      ONE=1.0D0
      TWO=2.0D0
      THREE=3.0D0
      TOL=1.E-5
      
      FAC1 = - THREE*SQRT(THREE)/TWO
      TR = TENSOR(1)+TENSOR(2)+TENSOR(3)
C
      sh1 = TENSOR(4)/(tr/three)
      sh2 = TENSOR(5)/(tr/three)
      sh3 = TENSOR(6)/(tr/three) 
      if (sh1.lt.tol .and. sh2.lt.tol .and. sh3.lt.tol )then
        t1 = TENSOR(1)
        t2 = TENSOR(2)
        t3 = TENSOR(3)
        goto 900 
      endif
C
      DO 5 I=1,6
        DEVIAT(I) = TENSOR(I)
 5    CONTINUE
CC      CALL MXCR8 (TENSOR,6,1,DEVIAT)
      DO 10 I=1,3
        DEVIAT(I) = DEVIAT(I)-TR/THREE
 10   CONTINUE
C
      DJ2 =     (DEVIAT(1)**2 + DEVIAT(2)**2 + DEVIAT(3)**2 + 
     +     TWO *(DEVIAT(4)**2 + DEVIAT(5)**2 + DEVIAT(6)**2 ))/TWO
C
      DJ3 = DEVIAT(1)*DEVIAT(2)*DEVIAT(3)+DEVIAT(4)*DEVIAT(5)*DEVIAT(6)+
     +      DEVIAT(4)*DEVIAT(5)*DEVIAT(6)-DEVIAT(2)*DEVIAT(6)*DEVIAT(6)-
     +      DEVIAT(1)*DEVIAT(5)*DEVIAT(5)-DEVIAT(3)*DEVIAT(4)*DEVIAT(4)
C
      IF (DJ2.LE.TOL) THEN
       T1 = TENSOR(1)
       T2 = TENSOR(2)
       T3 = TENSOR(3)
       GOTO 900
      ENDIF 
C
      FAC1 = (FAC1*DJ3/(DJ2*SQRT(DJ2)))
      PI  = ACOS (-ONE)
C
      IF (ABS(FAC1).GT.ONE) THEN
        ANG = PI/TWO
      ELSE   
        ANG = ASIN (FAC1)
      ENDIF
C
      ANG = ANG/THREE 
C
      FAC = TWO*SQRT(DJ2/THREE)
C
      T1 = FAC*SIN(ANG+TWO*PI/THREE) + TR/THREE
      T2 = FAC*SIN(ANG) + TR/THREE
      T3 = FAC*SIN(ANG+TWO*TWO*PI/THREE) + TR/THREE
C
 900  FACTOR = (T1 + ABS(T1) + T2 + ABS(T2) + T3 + ABS(T3))/TWO
      DEN = ABS(T1)+ABS(T2)+ABS(T3)
      IF (DEN.EQ.ZERO) GOTO 950
      FACTOR = FACTOR/DEN
C
 950  CONTINUE      
C
      RETURN
      END
c=======================================================================
      subroutine c1inv3det (a,det)                                                 
c-----------------------------------------------------------------------        
c!    inversion of unsymmetric  3x3  matrix  a                                  
c-----------------------------------------------------------------------        
      implicit real*8 (a-h,o-z)                                                 
      dimension a(3,3)                                                          
c                                                                               
      b11 = a(1,1)                                                              
      b12 = a(1,2)                                                              
      b13 = a(1,3)                                                              
      b21 = a(2,1)                                                              
      b22 = a(2,2)                                                              
      b23 = a(2,3)                                                              
      b31 = a(3,1)                                                              
      b32 = a(3,2)                                                              
      b33 = a(3,3)                                                              
c                                                                               
      a(1,1) =   b22*b33 - b32*b23                                              
      a(2,1) = - b21*b33 + b31*b23                                              
      a(3,1) =   b21*b32 - b31*b22                                              
      a(1,2) = - b12*b33 + b32*b13                                              
      a(2,2) =   b11*b33 - b31*b13                                              
      a(3,2) = - b11*b32 + b31*b12                                              
      a(1,3) =   b12*b23 - b22*b13                                              
      a(2,3) = - b11*b23 + b21*b13                                              
      a(3,3) =   b11*b22 - b21*b12                                              
c                                                                               
      det = b11*a(1,1) + b12*a(2,1) + b13*a(3,1)                                
      do 10 i=1,3                                                               
       do 10 j=1,3                                                              
   10   a(i,j) = a(i,j)/det                                                     
c                                                                               
      return                                                                    
      end                                                                       
c=======================================================================
      subroutine c1inv3 (a)
c-----------------------------------------------------------------------
c     inversion of unsymmetric  3x3  matrix  a
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(3,3)
c
      b11 = a(1,1)
      b12 = a(1,2)
      b13 = a(1,3)
      b21 = a(2,1)
      b22 = a(2,2)
      b23 = a(2,3)
      b31 = a(3,1)
      b32 = a(3,2)
      b33 = a(3,3)
c
      a(1,1) =   b22*b33 - b32*b23
      a(2,1) = - b21*b33 + b31*b23
      a(3,1) =   b21*b32 - b31*b22
      a(1,2) = - b12*b33 + b32*b13
      a(2,2) =   b11*b33 - b31*b13
      a(3,2) = - b11*b32 + b31*b12
      a(1,3) =   b12*b23 - b22*b13
      a(2,3) = - b11*b23 + b21*b13
      a(3,3) =   b11*b22 - b21*b12
c
      det = b11*a(1,1) + b12*a(2,1) + b13*a(3,1)
      do 10 i=1,3
       do 10 j=1,3
   10   a(i,j) = a(i,j)/det
c
      return
      end
      

