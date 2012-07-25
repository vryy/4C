c
c     Maintainer: Michael Gee 
c                 gee@statik.uni-stuttgart.de 
c                 http://www.uni-stuttgart.de/ibs/members/gee/ 
c                 0711 - 685-6572 
c
c     ---------------------------------------------------------------  
C======================================================================
      SUBROUTINE S8JACB (A,V)
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
C      DIMENSION A(KKK,KKK),V(KKK,KKK)
      DIMENSION A(3,3),V(3,3)
C
      ITM=200
      ERR=1.D-12
      ZERO=0.0D0
      ONE=1.0D0
      ZTOLP=1.D-11
      TWO=2.0D0
      FOUR=4.0D0
      KKK=3
      IT=0
      T1 = 0.0
      IR = 0
      IC = 0
C
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
