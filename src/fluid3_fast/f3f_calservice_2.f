C--------------------------------------------------------------------------
C
C \file
C \brief service routines for fluid3_fast element
C
C
C <pre>
C Maintainer: Malte Neumann
C             neumann@statik.uni-stuttgart.de
C             http://www.uni-stuttgart.de/ibs/members/neumann/
C             0711 - 685-6121
C </pre>
C
C--------------------------------------------------------------------------



C--------------------------------------------------------------------------
C
C \brief routine to calculate velocities at integration point
C
C
C \param velint()()      real*8   (o) vel at int point
C \param funct()         real*8   (i) shape functions
C \param evel()()()      real*8   (i) element velocities
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fveli( velint, funct, evel, sizevec)

      implicit none

      integer sizevec(6)
      real*8 velint(sizevec(4),3)
      real*8 funct(sizevec(1))
      real*8 evel(sizevec(4),sizevec(1),3)

      integer i,j,k
      real*8 ZERO

      ZERO = 0.0

      do i=1,3
        do k=1,sizevec(5)
          velint(k,i)=ZERO
          do j=1,sizevec(2)
            velint(k,i) = velint(k,i)+funct(j)*evel(k,j,i)
          enddo
        enddo
      enddo


      return
      end subroutine






C--------------------------------------------------------------------------
C
C \brief routine to calculate velocity derivatives at integration point
C
C In this routine the derivatives of the velocity w.r.t x/y are calculated
C vderxy[0][2] = Ux,z
C
C \param vderxy()()()    real*8   (o) global vel derivatives
C \param derxy()()()     real*8   (i) global derivatives
C \param evel()()()      real*8   (i) element velocities
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fvder(vderxy, derxy, evel, sizevec)

      implicit none

      integer sizevec(6)
      real*8 vderxy(sizevec(4),3,3)
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 evel(sizevec(4),sizevec(1),3)

      integer i,j,k
      real*8 ZERO

      ZERO = 0.0

      do i=1,3
        do k=1,sizevec(5)
          vderxy(k,i,1)=ZERO
          vderxy(k,i,2)=ZERO
          vderxy(k,i,3)=ZERO
          do j=1,sizevec(2)
            vderxy(k,i,1)=vderxy(k,i,1)+derxy(k,j,i)*evel(k,j,1)
            vderxy(k,i,2)=vderxy(k,i,2)+derxy(k,j,i)*evel(k,j,2)
            vderxy(k,i,3)=vderxy(k,i,3)+derxy(k,j,i)*evel(k,j,3)
          enddo
        enddo
      enddo

      return
      end subroutine






C--------------------------------------------------------------------------
C
C \brief routine to calculate 2nd velocity derivatives at integration point
C
C In this routine the 2nd derivatives of the velocity
C w.r.t x/y/z are calculated
C    vderxy2[0][0] = Ux,xx
C    vderxy2[0][3] = Ux,xy
C    vderxy2[1][4] = Ux,xz
C    vderxy2[2][5] = Ux,yz
C
C \param vderxy2()()()   real*8   (o) 2nd global vel derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param evel()()()      real*8   (i) element velocities
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fvder2(vderxy2, derxy2, evel, sizevec)

      implicit none

      integer sizevec(6)
      real*8 vderxy2(sizevec(4),6,3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 evel(sizevec(4),sizevec(1),3)

      integer i,j,k
      real*8 ZERO

      ZERO = 0.0

      do i=1,6
        do k=1,sizevec(5)
          vderxy2(k,i,1)=ZERO
          vderxy2(k,i,2)=ZERO
          vderxy2(k,i,3)=ZERO
          do j=1,sizevec(2)
            vderxy2(k,i,1)=vderxy2(k,i,1)+derxy2(k,j,i)*evel(k,j,1)
            vderxy2(k,i,2)=vderxy2(k,i,2)+derxy2(k,j,i)*evel(k,j,2)
            vderxy2(k,i,3)=vderxy2(k,i,3)+derxy2(k,j,i)*evel(k,j,3)
          enddo
        enddo
      enddo

      return
      end subroutine



C--------------------------------------------------------------------------
C
C \brief convective velocities
C
C in this routine the convective velocity is calculated at the
C integration point:
C  u * grad(u)
C  e.g. 3D: COVx = Ux*Ux,x + Uy*Ux,y + Uz*Ux,z
C
C \param vderxy()()()    real*8   (i) global vel derivatives
C \param velint()()      real*8   (i) vel at int point
C \param covint()()      real*8   (o) conv. vel at int point
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fcovi( vderxy, velint, covint, sizevec)

      implicit none

      integer sizevec(6)
      real*8 vderxy(sizevec(4),3,3)
      real*8 velint(sizevec(4),3)
      real*8 covint(sizevec(4),3)

      integer i,j,k
      real*8 ZERO


      ZERO = 0.0

      do i=1,3
        do k=1,sizevec(5)
          covint(k,i)=ZERO
          do j=1,3
            covint(k,i)=covint(k,i)+velint(k,j)*vderxy(k,j,i)
          enddo
        enddo
      enddo

      return
      end subroutine






C--------------------------------------------------------------------------
C
C \brief routine to calculate pressure at integration point
C
C
C \param preint()        real*8   (i) pres at integration point
C \param funct()         real*8   (i) shape functions
C \param epren()()       real*8   (i) pres at time n
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fprei(preint, funct, epre, sizevec)

      implicit none

      integer sizevec(6)
      real*8 preint(sizevec(4))
      real*8 funct(sizevec(1))
      real*8 epre(sizevec(4),sizevec(1))

      integer j,k
      real*8 ZERO

      ZERO = 0.0

      do k=1,sizevec(5)
        preint(k) = ZERO
        do j=1,sizevec(2)
          preint(k) =preint(k) + funct(j) * epre(k,j)
        enddo
      enddo

      return
      end subroutine





C--------------------------------------------------------------------------
C
C \brief routine to calculate pressure derivatives at integration point
C
C In this routine derivatives of the pressure w.r.t x/y/z are calculated
C
C \param pderxy()()      real*8   (i) global pressure derivatives
C \param derxy()()()     real*8   (i) global derivatives
C \param epre()()        real*8   (i) pres at time n
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fpder(pderxy, derxy, epre, sizevec)

      implicit none

      integer sizevec(6)
      real*8 pderxy(sizevec(4),3)
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 epre(sizevec(4),sizevec(1))

      integer i,j,k

      do i=1,3
        do k=1,sizevec(5)
          pderxy(k,i) = 0.0
          do j=1,sizevec(2)
            pderxy(k,i)=pderxy(k,i)+derxy(k,j,i)*epre(k,j)
          enddo
        enddo
      enddo

      return
      end subroutine







C--------------------------------------------------------------------------
C
C \brief add estiff and emass
C
C
C \param estif()()()     real*8   (i/o) element stiffness matrix
C \param emass()()()     real*8   (i) element mass matrix
C \param thsl            real*8   (i)
C \param nis             integer  (i) flag
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fmast(estif, emass, thsl, nis, sizevec)

      implicit none

c     parameters
      integer sizevec(6)

      real*8   estif(sizevec(4),sizevec(3),sizevec(3))
      real*8   emass(sizevec(4),sizevec(3),sizevec(3))
      real*8   thsl
      integer  nis

c     some helpers
      integer  l,i,j
      integer  rows


      rows = sizevec(2) * 4

      if (nis.eq.0) then

        do l=1,sizevec(5)
          do i=1,rows
            do j=1,rows

              estif(l,i,j) = estif(l,i,j) * thsl + emass(l,i,j)

            enddo
          enddo
        enddo

      else

        do l=1,sizevec(5)
          do i=1,rows
            do j=1,rows

              estif(l,i,j) = estif(l,i,j) * thsl

            enddo
          enddo
        enddo

      endif


 999  return
      end



      SUBROUTINE f3finv6(B,IRC)
C----------------------------------------------------------------------C
C     INVERTIERUNG EINER 6x6 MATRIX
C         B   -->  INPUT & INVERSE VON A 
C         IRC -->  RETURN CODE ( IRC = 0   OK
C                                IRC = 1   ERROR: MATRIX SINGULAER )
C----------------------------------------------------------------------C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(6,6),MX(6),MY(6)

      TOL=1.0D-20
      N=6
      NX = 0
      NY = 0

      IRC=0

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
   60       CONTINUE
          END IF
   50   CONTINUE
C       FALLS DAS PIVOTELEMENT NULL IST,IST DIE MATRIX SINGULAER
        IF(DABS(PIVO).LT.TOL) THEN
          IRC=1
          RETURN
        END IF
C       MERKEN DER INDIZES DES PIVOTELEMENTES
        MX(NX)=NY
        MY(NY)=NX
C       BERECHNUNG DER MATRIXELEMENTE
        HILF=1.0/PIVO
        DO 70 J=1,N
          IF(J.NE.NX) THEN
            FAKTOR=B(J,NY)*HILF
            DO 80 K=1,N
              B(J,K)=B(J,K)-B(NX,K)*FAKTOR
   80       CONTINUE
            B(J,NY)=-FAKTOR
          END IF
   70   CONTINUE
        DO 90 K=1,N
          B(NX,K)=B(NX,K)*HILF
   90   CONTINUE
        B(NX,NY)=HILF
   40 CONTINUE

C     ZEILEN UND SPALTENVERTAUSCHUNG RUECKGAENGIG MACHEN
      DO 100 I=1,N-1
        DO 110 M=I,N
          IF(MX(M).EQ.I) GOTO 120
  110   CONTINUE
  120   J=M
        IF(J.NE.I) THEN
          DO 130 K=1,N
            H=B(I,K)
            B(I,K)=B(J,K)
            B(J,K)=H
  130     CONTINUE
          MX(J)=MX(I)
          MX(I)=I
        END IF
        DO 140 M=I,N
          IF(MY(M).EQ.I) GO TO 150
  140   CONTINUE
  150   J=M
        IF(J.NE.I) THEN
          DO 160 K=1,N
            H=B(K,I)
            B(K,I)=B(K,J)
            B(K,J)=H
  160     CONTINUE
          MY(J)=MY(I)
          MY(I)=I
        END IF
  100 CONTINUE


      RETURN
      END







      SUBROUTINE f3finv62(B,IRC,loop,loopl)

      implicit    none
      integer     loop
      integer     loopl
      real*8      B(loopl,32,32)
      integer     IRC

      integer     MX(loopl,6)
      integer     MY(loopl,6)
      integer     NX(loopl)
      integer     NY(loopl)
      real*8      FAKTOR
      real*8      HILF(loopl)
      real*8      PIVO(loopl)
      real*8      TOL
      integer     L,I,J,IX,IY,M
      real*8      H


      TOL=1.0D-20

      IRC=0


C     VORBESETZEN DER PIVOT-VEKTOREN MX UND MY MIT NULL
      DO 30 I=1,6
        DO 25 L=1,loop
          MX(L,I)=0
          MY(L,I)=0
   25   CONTINUE
   30 CONTINUE


C     BESTIMMUNG DES PIVOTELEMENTES
      DO 40 I=1,6

        DO 1000 L=1,loop
          PIVO(L)=0.0
 1000   continue

          DO 50 IX=1,6
            DO 60 IY=1,6
              DO 1100 L=1,loop

                IF(MX(L,IX).EQ.0.and.MY(L,IY).EQ.0.and.
     &              DABS(B(L,IX,IY)).GT.DABS(PIVO(L))) THEN
                  PIVO(L)=B(L,IX,IY)
                  NX(L)=IX
                  NY(L)=IY
                END IF

 1100         continue
   60       CONTINUE
   50     CONTINUE




        DO 1500 L=1,loop
C         MERKEN DER INDIZES DES PIVOTELEMENTES
          MX(L,NX(L))=NY(L)
          MY(L,NY(L))=NX(L)
C         BERECHNUNG DER MATRIXELEMENTE
          HILF(L)=1.0/PIVO(L)
 1500   continue


        DO 70 J=1,6

          DO 1600 L=1,loop

            IF(J.NE.NX(L)) THEN
              FAKTOR=B(L,J,NY(L))*HILF(L)

              B(L,J,1)=B(L,J,1)-B(L,NX(L),1)*FAKTOR
              B(L,J,2)=B(L,J,2)-B(L,NX(L),2)*FAKTOR
              B(L,J,3)=B(L,J,3)-B(L,NX(L),3)*FAKTOR
              B(L,J,4)=B(L,J,4)-B(L,NX(L),4)*FAKTOR
              B(L,J,5)=B(L,J,5)-B(L,NX(L),5)*FAKTOR
              B(L,J,6)=B(L,J,6)-B(L,NX(L),6)*FAKTOR

              B(L,J,NY(L))=-FAKTOR
            END IF

 1600     continue

   70   CONTINUE



        DO 1700 L=1,loop
          B(L,NX(L),1)=B(L,NX(L),1)*HILF(L)
          B(L,NX(L),2)=B(L,NX(L),2)*HILF(L)
          B(L,NX(L),3)=B(L,NX(L),3)*HILF(L)
          B(L,NX(L),4)=B(L,NX(L),4)*HILF(L)
          B(L,NX(L),5)=B(L,NX(L),5)*HILF(L)
          B(L,NX(L),6)=B(L,NX(L),6)*HILF(L)

          B(L,NX(L),NY(L))=HILF(L)
 1700   continue

   40 CONTINUE




C       ZEILEN UND SPALTENVERTAUSCHUNG RUECKGAENGIG MACHEN
      DO 100 I=1,5
        DO 110 M=I,6
          DO 2000 L=1,loop
            IF(MX(L,M).EQ.I.and.M.NE.I) then
              H=B(L,I,1)
              B(L,I,1)=B(L,M,1)
              B(L,M,1)=H
              H=B(L,I,2)
              B(L,I,2)=B(L,M,2)
              B(L,M,2)=H
              H=B(L,I,3)
              B(L,I,3)=B(L,M,3)
              B(L,M,3)=H
              H=B(L,I,4)
              B(L,I,4)=B(L,M,4)
              B(L,M,4)=H
              H=B(L,I,5)
              B(L,I,5)=B(L,M,5)
              B(L,M,5)=H
              H=B(L,I,6)
              B(L,I,6)=B(L,M,6)
              B(L,M,6)=H

              MX(L,M)=MX(L,I)
              MX(L,I)=I
            endif
 2000     continue
  110   continue
  100 continue



      DO 102 I=1,5
        DO 140 M=I,6
          DO 2002 L=1,loop
            IF(MY(L,M).EQ.I.and.M.NE.I) then
              H=B(L,1,I)
              B(L,1,I)=B(L,1,M)
              B(L,1,M)=H
              H=B(L,2,I)
              B(L,2,I)=B(L,2,M)
              B(L,2,M)=H
              H=B(L,3,I)
              B(L,3,I)=B(L,3,M)
              B(L,3,M)=H
              H=B(L,4,I)
              B(L,4,I)=B(L,4,M)
              B(L,4,M)=H
              H=B(L,5,I)
              B(L,5,I)=B(L,5,M)
              B(L,5,M)=H
              H=B(L,6,I)
              B(L,6,I)=B(L,6,M)
              B(L,6,M)=H

              MY(L,M)=MY(L,I)
              MY(L,I)=I
            endif
 2002     continue
  140   continue
  102 CONTINUE


      return

      end subroutine





C--------------------------------------------------------------------------
C
C \brief calculate stresse at integration point
C
C
C \param preint()        real*8   (i) pres at integration point
C \param vderxy()()()    real*8   (i) global vel derivatives
C \param sigint()()()    real*8   (i) stresses at GAUSS point
C \param visc            real*8   (i) vicosity
C \param iv              integer  (i) position
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fsigint( preint, vderxy, sigint, visc, iv, sizevec)

      implicit none

      integer sizevec(6)
      real*8 preint(sizevec(4))
      real*8 vderxy(sizevec(4),3,3)
      real*8 sigint(sizevec(4),sizevec(6),6)
      real*8 visc
      integer iv

      integer k

c              | Ux,x    Ux,y    Ux,z |
c              |                      |
c     vderxy = | Uy,x    Uy,y    Uy,z |
c              |                      |
c              | Uz,x    Uz,y    Uz,z |
c
c
c                  | Ux,x+Ux,x   Ux,y+Uy,x   Ux,z+Uz,x |
c                  |                                   |
c     eps(u) = 1/2 |             Uy,y+Uy,y   Uy,z+Uz,y |
c                  |                                   |
c                  |   symm.                 Uz,z+Uz,z |
c
c
c     SIGMA = -p_real*I + 2*nue * eps(u)

      do k=1,sizevec(5)
        sigint(k,iv,1) = -preint(k) + 2.0*visc * vderxy(k,1,1)
        sigint(k,iv,2) = -preint(k) + 2.0*visc * vderxy(k,2,2)
        sigint(k,iv,3) = -preint(k) + 2.0*visc * vderxy(k,3,3)
        sigint(k,iv,4) = (vderxy(k,2,1) + vderxy(k,1,2))*visc
        sigint(k,iv,5) = (vderxy(k,3,2) + vderxy(k,2,3))*visc
        sigint(k,iv,6) = (vderxy(k,3,1) + vderxy(k,1,3))*visc
      enddo



      return
      end subroutine

