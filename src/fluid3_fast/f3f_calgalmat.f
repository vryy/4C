C--------------------------------------------------------------------------
C
C \file
C \brief evaluate galerkin part of stiffness matrix
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




c     In this routine the galerkin part of matrix Kvv is calculated:
c     In this routine the galerkin part of matrix Kvp is calculated:
C--------------------------------------------------------------------------
C
C \brief evaluate galerkin part of Kvv and Kvp
C
C In this routine the galerkin part of matrix Kvv is calculated:
C
C EULER/ALE:
C     /
C    |  2 * nue * eps(v) : eps(u)   d_omega
C   /
C
C EULER/ALE:
C     /
C    |  v * u_old * grad(u)     d_omega
C   /
C
C     /
C    |  v * u * grad(u_old)     d_omega
C   /
C
C ALE:
C
C ale-convective velocity is split into
C  c = u - u_G
C    --> the known nonlinear term known from the EULER-case
C    --> a new term:
C
C     /
C  - |  v * u_G * grad(u)     d_omega
C   /
C
C    --> this is a linear term inside the domain and for explicit treatement
C        of  a free surface
C    --> for implicit treatement of the free surface this term is nonlinear
C        so it depends on the nonlinear iteration scheme if this term has to
C        be taken into account on the LHS!
C
C
C In this routine the galerkin part of matrix Kvp is calculated:
C
C     /
C    |  - div(v) * p     d_omega
C   /
C
C     /
C    | - q * div(u)      d_omega
C   /
C
C
C \param estif()()()     real*8   (o) element stiffness matrix
C \param velint()()      real*8   (i) vel at int point
C \param gridvint()()    real*8   (i) gridvel at int point
C \param vderxy()()()    real*8   (i) global vel derivatives
C \param funct()         real*8   (i) shape functions
C \param derxy()()()     real*8   (i) global derivatives
C \param fac()           real*8   (i) gauss integration factor
C \param paravec(2)      real*8   (i) some parameters
C \param flagvec(7)      integer  (i) some flags
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fcalgalk(estif, velint, gridvint, vderxy, funct,
     &     derxy, fac, paravec, flagvec, sizevec)

      implicit none

      real*8  paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 estif(sizevec(4),sizevec(3),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 gridvint(sizevec(4),3)
      real*8 vderxy(sizevec(4),3,3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 fac(sizevec(4))

      integer irow,icol,irn,icn,k
      real*8 c,aux

      integer posc
      real*8 aux2(3)


C--------------------------------------------------------------------------
C   Calculate full Galerkin part of matrix K:
c
C    /
C   |  2 * nue * eps(v) : eps(u)   d_omega
C  /
C--------------------------------------------------------------------------
      icol=1
      do icn=1,sizevec(2)
        irow=1
        do irn=1,sizevec(2)

          do k=1,sizevec(5)
            c=fac(k)*paravec(2)
            aux = derxy(k,irn,1)*derxy(k,icn,1)
     &      + derxy(k,irn,2)*derxy(k,icn,2)
     &      + derxy(k,irn,3)*derxy(k,icn,3)

            estif(k,icol,irow)  =estif(k,icol,irow)
     &                  +c*(aux+derxy(k,irn,1)*derxy(k,icn,1))
            estif(k,icol,irow+1)=estif(k,icol,irow+1)
     &                  +c*(    derxy(k,irn,1)*derxy(k,icn,2))
            estif(k,icol,irow+2)=estif(k,icol,irow+2)
     &                  +c*(    derxy(k,irn,1)*derxy(k,icn,3))

            estif(k,icol+1,irow)  =estif(k,icol+1,irow)
     &                  +c*(    derxy(k,irn,2)*derxy(k,icn,1))
            estif(k,icol+1,irow+1)=estif(k,icol+1,irow+1)
     &                  +c*(aux+derxy(k,irn,2)*derxy(k,icn,2))
            estif(k,icol+1,irow+2)=estif(k,icol+1,irow+2)
     &                  +c*(    derxy(k,irn,2)*derxy(k,icn,3))

            estif(k,icol+2,irow)  =estif(k,icol+2,irow)
     &                  +c*(    derxy(k,irn,3)*derxy(k,icn,1))
            estif(k,icol+2,irow+1)=estif(k,icol+2,irow+1)
     &                  +c*(    derxy(k,irn,3)*derxy(k,icn,2))
            estif(k,icol+2,irow+2)=estif(k,icol+2,irow+2)
     &                  +c*(aux+derxy(k,irn,3)*derxy(k,icn,3))

          enddo


          irow =irow + 4
        enddo
        icol =icol + 4
      enddo

C--------------------------------------------------------------------------
c   Calculate full Galerkin part of matrix Nc(u):
c
c    /
c   |  v * u_old * grad(u)     d_omega
c   /
C--------------------------------------------------------------------------

      icol=1
      do icn=1,sizevec(2)
        irow=1
        do irn=1,sizevec(2)

          do k=1,sizevec(5)
            aux = (velint(k,1)*derxy(k,icn,1)
     &            +  velint(k,2)*derxy(k,icn,2)
     &            +  velint(k,3)*derxy(k,icn,3))*funct(irn)*fac(k)
            estif(k,icol,irow)    =estif(k,icol,irow)  +aux
            estif(k,icol+1,irow+1)=estif(k,icol+1,irow+1)+aux
            estif(k,icol+2,irow+2)=estif(k,icol+2,irow+2)+aux
          enddo


          irow = irow + 4
        enddo
        icol = icol + 4
      enddo

C--------------------------------------------------------------------------
c   Calculate full Galerkin part of matrix Nr(u):
c
c    /
c   |  v * u * grad(u_old)     d_omega
c  /
C--------------------------------------------------------------------------
      if(flagvec(4).ne.0) then
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do k=1,sizevec(5)
              aux = funct(irn)*funct(icn)*fac(k)

              estif(k,icol,irow)  =estif(k,icol,irow)
     &             +aux*vderxy(k,1,1)
              estif(k,icol,irow+1)=estif(k,icol,irow+1)
     &             +aux*vderxy(k,1,2)
              estif(k,icol,irow+2)=estif(k,icol,irow+2)
     &             +aux*vderxy(k,1,3)

              estif(k,icol+1,irow)  =estif(k,icol+1,irow)
     &             +aux*vderxy(k,2,1)
              estif(k,icol+1,irow+1)=estif(k,icol+1,irow+1)
     &             +aux*vderxy(k,2,2)
              estif(k,icol+1,irow+2)=estif(k,icol+1,irow+2)
     &            +aux*vderxy(k,2,3)

              estif(k,icol+2,irow)  =estif(k,icol+2,irow)
     &             +aux*vderxy(k,3,1)
              estif(k,icol+2,irow+1)=estif(k,icol+2,irow+1)
     &             +aux*vderxy(k,3,2)
              estif(k,icol+2,irow+2)=estif(k,icol+2,irow+2)
     &             +aux*vderxy(k,3,3)
            enddo

            irow =irow + 4
          enddo
          icol =icol + 4
        enddo
      endif

C--------------------------------------------------------------------------
c   Calculate full Galerkin part due to split of ALE-convective velocity:
c     /
c  - |  v * u_G_old * grad(u)     d_omega
c   /
c   REMARK:
c    - this term is linear inside the domain and for explicit free surface
c    - for implicit free surface this term is nonlinear (Nc), since u_G is
c      unknown at the free surface. So it depends on the nonlinear iteration
c      scheme if this term has to be taken into account:
c        - Newton: YES
c        - fixpoint-like: YES
C--------------------------------------------------------------------------
c     evaluate only for ALE
      if(flagvec(7).ne.0) then
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do k=1,sizevec(5)
              aux = (gridvint(k,1)*derxy(k,icn,1)
     &            +  gridvint(k,2)*derxy(k,icn,2)
     &            +  gridvint(k,3)*derxy(k,icn,3))*funct(irn)*fac(k)
              estif(k,icol,irow)     = estif(k,icol,irow)     - aux
              estif(k,icol+1,irow+1) = estif(k,icol+1,irow+1) - aux
              estif(k,icol+2,irow+2) = estif(k,icol+2,irow+2) - aux
            enddo

            irow = irow + 4
          enddo
          icol = icol + 4
        enddo

      end if


C--------------------------------------------------------------------------
c   Calculate full Galerkin part of matrix Kvp:
c
c    /
c   |  - div(v) * p     d_omega
c  /
c
c   and matrix Kpv:
c
c    /
c   | - q * div(u)      d_omega
c  /
C--------------------------------------------------------------------------
      posc = 0
      do icol=1,sizevec(2)
        irow = 1

        posc = posc + 4
        do irn=1,sizevec(2)

          do k=1,sizevec(5)
            aux2(1) = funct(icol)*derxy(k,irn,1)*fac(k)
            aux2(2) = funct(icol)*derxy(k,irn,2)*fac(k)
            aux2(3) = funct(icol)*derxy(k,irn,3)*fac(k)

            estif(k,posc,irow)=estif(k,posc,irow)-aux2(1)
            estif(k,irow,posc)=estif(k,irow,posc)-aux2(1)

            estif(k,posc,irow+1)=estif(k,posc,irow+1)-aux2(2)
            estif(k,irow+1,posc)=estif(k,irow+1,posc)-aux2(2)

            estif(k,posc,irow+2)=estif(k,posc,irow+2)-aux2(3)
            estif(k,irow+2,posc)=estif(k,irow+2,posc)-aux2(3)
          enddo

          irow = irow + 4
        enddo
      enddo

      return
      end subroutine




c     In this routine the galerkin part of matrix Mvv is calculated:
C--------------------------------------------------------------------------
C
C \brief evaluate galerkin part of Mvv
C
C In this routine the galerkin part of matrix Mvv is calculated:
C
C     /
C    |  v * u    d_omega
C   /
C
C
C \param estif()()()     real*8   (o) element stiffness matrix
C \param funct()         real*8   (i) shape functions
C \param fac()           real*8   (i) gauss integration factor
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fcalgalm( estif, funct, fac, sizevec)

      implicit none

      integer sizevec(6)
      real*8 estif(sizevec(4),sizevec(3),sizevec(3))
      real*8 funct(sizevec(1))
      real*8 fac(sizevec(4))

      integer irow,icol,irn,icn,k
      integer nvdfe
      real*8 aux

      nvdfe = (3+1)*sizevec(2)

C--------------------------------------------------------------------------
c   Calculate full Galerkin part of matrix Mvv:
c
c    /
c   |  v * u    d_omega
c  /
C--------------------------------------------------------------------------
      icn=0
      do icol=1,nvdfe,4
        icn = icn+1
        irn=0
        do irow=1,nvdfe,4
          irn = irn+1

          do k=1,sizevec(5)
            aux = funct(icn)*funct(irn)*fac(k)
            estif(k,icol,irow)    =estif(k,icol,irow)  +aux
            estif(k,icol+1,irow+1)=estif(k,icol+1,irow+1)+aux
            estif(k,icol+2,irow+2)=estif(k,icol+2,irow+2)+aux
          enddo

        enddo
      enddo

      return
      end subroutine


