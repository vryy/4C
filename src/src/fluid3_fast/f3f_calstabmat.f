C--------------------------------------------------------------------------
C
C \file
C \brief stabilisation part of element stiffness matrix for fluid3_fast
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

c     In this routine the stabilisation part of matrix Kpv is calculated:c

c     In this routine the stabilisation part of matrix Kpp is calculated:c
C--------------------------------------------------------------------------
C
C \brief evaluate stabilisaton part of K
C
C In this routine the stabilisation part of matrix Kvv is calculated:
C
C EULER/ALE:
C     /
C    |  tau_c * div(v) * div(u)   d_omega
C   /
C
C EULER:
C
C     /
C    |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
C   /
C
C     /
C    |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
C   /
C
C ALE:
C
C     /
C    |  tau_mu * c * grad(v) * u_old * grad(u)   d_omega
C   /
C
C     /
C    |  tau_mu * c * grad(v) * u * grad(u_old)   d_omega
C   /
C
C     /
C  - |  tau_mu * c * grad(v) * u_G_old * grad(u)   d_omega
C   /
C
C EULER:
C     /
C    |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
C   /
C
C ALE:
C     /
C    |  -tau_mu * 2 * nue * c * grad(v) * div(eps(u))   d_omega
C   /
C
C EULER/ALE:
C     /
C    |  +/- tau_mp  * 4 * nue**2 * div(eps(v))  div(eps(u))d_omega
C   /
C
C EULER/ALE:
C     /
C    |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
C   /
C
C     /
C    |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
C   /
C
C ALE:
C     /
C    |  +/- tau_mp  * 2 * nue * div(eps(v)) * u_G_old * grad(u) d_omega
C   /
C
C
C In this routine the stabilisation part of matrix Kvp is calculated:
C
C EULER:
C     /
C    |  tau_mu * u_old * grad(v) * grad(p)   d_omega
C   /
C
C ALE:
C     /
C    |  tau_mu * c * grad(v) * grad(p)   d_omega
C   /
C
C EULER/ALE:
C     /
C    |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
C   /
C
C
C In this routine the stabilisation part of matrix Kpv is calculated:
C
C EULER/ALE:
C     /
C    |  - tau_mp * grad(q) * u_old * grad(u) d_omega
C   /
C
C     /
C    |  - tau_mp * grad(q) * u * grad(u_old) d_omega
C   /
C
C ALE:
C     /
C    |  tau_mp * grad(q) * u_G_old * grad(u) d_omega
C   /
C
C EULER/ALE:
C     /
C    |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
C   /
C
C
C In this routine the stabilisation part of matrix Kpp is calculated:
C
C     /
C    |  - tau_mp * grad(q) *grad(p) d_omega
C   /
C
C
C
C \param estif()()()     real*8   (o) element stiffness matrix
C \param velint()()      real*8   (i) vel at int point
C \param vel2int()()     real*8   (i) vel at int point
C \param gridvint()()    real*8   (i) gridvel at int point
C \param vderxy()()()    real*8   (i) global vel derivatives
C \param funct()         real*8   (i) shape functions
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param fac()           real*8   (i) gauss integration factor
C \param tau()()         real*8   (i) stabilisation parameter
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
      subroutine f3fcalstabk(estif,velint,vel2int,gridvint,vderxy,
     &      funct,derxy,derxy2,fac, tau, paravec, flagvec, sizevec)

      implicit none

      real*8 paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 estif(sizevec(4),sizevec(3),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 vel2int(sizevec(4),3)
      real*8 gridvint(sizevec(4),3)
      real*8 vderxy(sizevec(4),3,3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 fac(sizevec(4))
      real*8 tau(sizevec(4),3)

      integer irow,icol,irn,icn,j,k
      real*8 c,cc
      real*8 aux,auxr,auxc

      integer posc

      integer posr


      cc = 0.0

C--------------------------------------------------------------------------
C   Calculate continuity stabilisation part:
c
C    /
C   |  tau_c * div(v) * div(u)   d_omega
C  /
C--------------------------------------------------------------------------
      if(flagvec(1).ne.0) then
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              c = fac(j)*tau(j,3)
              estif(j,icol,irow)   = estif(j,icol,irow) +
     &    derxy(j,icn,1)*derxy(j,irn,1)*c
              estif(j,icol,irow+1) = estif(j,icol,irow+1) +
     &    derxy(j,icn,1)*derxy(j,irn,2)*c
              estif(j,icol,irow+2) = estif(j,icol,irow+2) +
     &    derxy(j,icn,1)*derxy(j,irn,3)*c

              estif(j,icol+1,irow)   = estif(j,icol+1,irow) +
     &    derxy(j,icn,2)*derxy(j,irn,1)*c
              estif(j,icol+1,irow+1) = estif(j,icol+1,irow+1) +
     &    derxy(j,icn,2)*derxy(j,irn,2)*c
              estif(j,icol+1,irow+2) = estif(j,icol+1,irow+2) +
     &    derxy(j,icn,2)*derxy(j,irn,3)*c

              estif(j,icol+2,irow)   = estif(j,icol+2,irow) +
     &    derxy(j,icn,3)*derxy(j,irn,1)*c
              estif(j,icol+2,irow+1) = estif(j,icol+2,irow+1) +
     &    derxy(j,icn,3)*derxy(j,irn,2)*c
              estif(j,icol+2,irow+2) = estif(j,icol+2,irow+2) +
     &    derxy(j,icn,3)*derxy(j,irn,3)*c
            enddo

            irow = irow + 4
          enddo
          icol = icol + 4
        enddo
      endif


c      calculate advection stabilisation part
      if(flagvec(2).ne.0) then

C--------------------------------------------------------------------------
c    Calculate advection stabilisation part Nc(u):
c
c EULER:
c     /
c    |  tau_mu * u_old * grad(v) * u_old * grad(u)   d_omega
c   /
c
c ALE:
c     /
c    |  tau_mu * c * grad(v) * u_old * grad(u)   d_omega
c   /
C--------------------------------------------------------------------------
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              cc=fac(j)*tau(j,1)
              auxc =(velint(j,1)*derxy(j,icn,1)
     &             + velint(j,2)*derxy(j,icn,2)
     &             + velint(j,3)*derxy(j,icn,3))*cc
              aux =(vel2int(j,1)*derxy(j,irn,1)
     &            + vel2int(j,2)*derxy(j,irn,2)
     &            + vel2int(j,3)*derxy(j,irn,3))*auxc
              estif(j,icol,irow)     = estif(j,icol,irow)     + aux
              estif(j,icol+1,irow+1) = estif(j,icol+1,irow+1) + aux
              estif(j,icol+2,irow+2) = estif(j,icol+2,irow+2) + aux
            enddo

            irow =irow + 4
          enddo
          icol = icol + 4
        enddo


C--------------------------------------------------------------------------
c    Calculate advection stabilisation part Nr(u):
c
c EULER:
c     /
c    |  tau_mu * u_old * grad(v) * u * grad(u_old)   d_omega
c   /
c
c ALE:
c     /
c    |  tau_mu * c * grad(v) * u * grad(u_old)   d_omega
c   /
C--------------------------------------------------------------------------
        if(flagvec(4).ne.0) then
          icol=1
          do icn=1,sizevec(2)

            irow=1
            do irn=1,sizevec(2)

              do j=1,sizevec(5)
                cc=fac(j)*tau(j,1)
                auxc = funct(icn)*cc
                aux = (vel2int(j,1)*derxy(j,irn,1)
     &               + vel2int(j,2)*derxy(j,irn,2)
     &               + vel2int(j,3)*derxy(j,irn,3))*auxc

                estif(j,icol,irow)   = estif(j,icol,irow)
     &               + aux*vderxy(j,1,1)
                estif(j,icol,irow+1) = estif(j,icol,irow+1)
     &               + aux*vderxy(j,1,2)
                estif(j,icol,irow+2) = estif(j,icol,irow+2)
     &               + aux*vderxy(j,1,3)

                estif(j,icol+1,irow)  =estif(j,icol+1,irow)
     &               + aux*vderxy(j,2,1)
                estif(j,icol+1,irow+1)=estif(j,icol+1,irow+1)
     &               + aux*vderxy(j,2,2)
                estif(j,icol+1,irow+2)=estif(j,icol+1,irow+2)
     &               + aux*vderxy(j,2,3)

                estif(j,icol+2,irow)  =estif(j,icol+2,irow)
     &               + aux*vderxy(j,3,1)
                estif(j,icol+2,irow+1)=estif(j,icol+2,irow+1)
     &               + aux*vderxy(j,3,2)
                estif(j,icol+2,irow+2)=estif(j,icol+2,irow+2)
     &               + aux*vderxy(j,3,3)
              enddo

              irow = irow + 4
            enddo
            icol = icol + 4
          enddo
        endif


C--------------------------------------------------------------------------
c    Calculate advection stabilisation part for ALE due to u_G:
c
c ALE:
c     /
c  - |  tau_mu * c * grad(v) * u_G_old * grad(u)   d_omega
c   /
c
c    REMARK:
c     - this term is linear inside the domain and for explicit free surface
c     - for implicit free surface this term is nonlinear (Nc), since u_G is
c       unknown at the free surface. So it depends on the nonlinear iteration
c       scheme if this term has to be taken into account:
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
              auxc = (gridvint(k,1)*derxy(k,icn,1)
     &              + gridvint(k,2)*derxy(k,icn,2)
     &              + gridvint(k,3)*derxy(k,icn,3))*cc
              aux = (vel2int(k,1)*derxy(k,icn,1)
     &            +  vel2int(k,2)*derxy(k,icn,2)
     &            +  vel2int(k,3)*derxy(k,icn,3))*auxc
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
c    Calculate advection stabilisation part for higher order elements:
c
c EULER:
c     /
c    |  -tau_mu * 2 * nue * u_old * grad(v) * div(eps(u))   d_omega
c   /
c
c ALE:
c     /
c    |  -tau_mu * 2 * nue * c * grad(v) * div(eps(u))   d_omega
c   /
C--------------------------------------------------------------------------
        if (flagvec(6).ne.0) then
          icol=1

          do icn=1,sizevec(2)
            irow=1

            do irn=1,sizevec(2)

              do j=1,sizevec(5)
                cc=fac(j)*tau(j,1)*paravec(2)
                auxc = derxy2(j,icn,1)+derxy2(j,icn,2)+derxy2(j,icn,3)
                aux =(vel2int(j,1)*derxy(j,irn,1)
     &              + vel2int(j,2)*derxy(j,irn,2)
     &              + vel2int(j,3)*derxy(j,irn,3))*cc

                estif(j,icol,irow)   = estif(j,icol,irow)  -
     &                  aux*(derxy2(j,icn,1)+ auxc)
                estif(j,icol,irow+1) = estif(j,icol,irow+1)-
     &                  aux* derxy2(j,icn,4)
                estif(j,icol,irow+2) = estif(j,icol,irow+2)-
     &                  aux* derxy2(j,icn,5)

                estif(j,icol+1,irow)   = estif(j,icol+1,irow) -
     &                  aux* derxy2(j,icn,4)
                estif(j,icol+1,irow+1) = estif(j,icol+1,irow+1)-
     &                  aux*(derxy2(j,icn,2)+ auxc)
                estif(j,icol+1,irow+2) = estif(j,icol+1,irow+2)-
     &                  aux* derxy2(j,icn,6)

                estif(j,icol+2,irow)   = estif(j,icol+2,irow)  -
     &                  aux* derxy2(j,icn,5)
                estif(j,icol+2,irow+1) = estif(j,icol+2,irow+1)-
     &                  aux* derxy2(j,icn,6)
                estif(j,icol+2,irow+2) = estif(j,icol+2,irow+2)-
     &                  aux*(derxy2(j,icn,3)+ auxc)
              enddo

              irow = irow + 4
            enddo
            icol = icol + 4
          enddo
        endif

      endif



      if((flagvec(6).ne.0).and.(flagvec(3).ne.0)) then

C--------------------------------------------------------------------------
c   Calculate viscous stabilisation part for higher order elements:
c
c    /
c   |  +/- tau_mp  * 4 * nue**2 * div(eps(v)) * div(eps(u))d_omega
c  /
C--------------------------------------------------------------------------
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              cc = fac(j) * tau(j,2) * paravec(2)*paravec(2)*paravec(1)
              auxc = derxy2(j,icn,1)+derxy2(j,icn,2)+derxy2(j,icn,3)
              auxr = derxy2(j,irn,1)+derxy2(j,irn,2)+derxy2(j,irn,3)

              estif(j,icol,irow)     = estif(j,icol,irow) +
     &   ((auxc + derxy2(j,icn,1))*(auxr + derxy2(j,irn,1))
     &          + derxy2(j,icn,4)*         derxy2(j,irn,4)
     &          + derxy2(j,icn,5)*         derxy2(j,irn,5) )*cc
              estif(j,icol,irow+1)   = estif(j,icol,irow+1) +
     &     ( derxy2(j,icn,4)*(auxr + derxy2(j,irn,1))
     &    + (auxc + derxy2(j,icn,2))*derxy2(j,irn,4)
     &            + derxy2(j,icn,6)* derxy2(j,irn,5) )*cc
              estif(j,icol,irow+2)   = estif(j,icol,irow+2) +
     &             ( derxy2(j,icn,5)*(auxr + derxy2(j,irn,1))
     &             + derxy2(j,icn,6)*derxy2(j,irn,4)
     &     + (auxc + derxy2(j,icn,3))*derxy2(j,irn,5) )*cc

              estif(j,icol+1,irow)   = estif(j,icol+1,irow) +
     &   ( (auxc + derxy2(j,icn,1))*derxy2(j,irn,4)
     &           + derxy2(j,icn,4)*(auxr + derxy2(j,irn,2))
     &           + derxy2(j,icn,5)*derxy2(j,irn,6) )*cc
              estif(j,icol+1,irow+1) = estif(j,icol+1,irow+1) +
     &            ( derxy2(j,icn,4)*derxy2(j,irn,4)
     &    + (auxc + derxy2(j,icn,2))*(auxr + derxy2(j,irn,2))
     &            + derxy2(j,icn,6)*derxy2(j,irn,6) )*cc
              estif(j,icol+1,irow+2) = estif(j,icol+1,irow+2) +
     &            ( derxy2(j,icn,5)*derxy2(j,irn,4)
     &            + derxy2(j,icn,6)*(auxr + derxy2(j,irn,2))
     &    + (auxc + derxy2(j,icn,3))*derxy2(j,irn,6) )*cc

              estif(j,icol+2,irow)  = estif(j,icol+2,irow) +
     &   ( (auxc + derxy2(j,icn,1))*derxy2(j,irn,5)
     &   + derxy2(j,icn,4)*derxy2(j,irn,6)
     &   + derxy2(j,icn,5)*(auxr + derxy2(j,irn,3)) )*cc
              estif(j,icol+2,irow+1) = estif(j,icol+2,irow+1) +
     &   ( derxy2(j,icn,4)*derxy2(j,irn,5)
     &    + (auxc + derxy2(j,icn,2))*derxy2(j,irn,6)
     &    + derxy2(j,icn,6)*(auxr + derxy2(j,irn,3)) )*cc
              estif(j,icol+2,irow+2) = estif(j,icol+2,irow+2) +
     &   ( derxy2(j,icn,5)*derxy2(j,irn,5)
     &  + derxy2(j,icn,6)*derxy2(j,irn,6)
     &  + (auxc+derxy2(j,irn,3))*(auxr+derxy2(j,irn,3)))*cc

            enddo

            irow = irow + 4
          enddo
          icol = icol + 4
        enddo


C--------------------------------------------------------------------------
c   Calculate viscous stabilisation part Nc(u) for higher order elements:
c
c    /
c   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u_old * grad(u) d_omega
c  /
C--------------------------------------------------------------------------
        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              cc = fac(j) * tau(j,2) * paravec(2) * paravec(1)
              aux =(velint(j,1)*derxy(j,icn,1)
     &              + velint(j,2)*derxy(j,icn,2)
     &              + velint(j,3)*derxy(j,icn,3))*cc
              auxr = derxy2(j,irn,1)+derxy2(j,irn,2)+derxy2(j,irn,3)
              estif(j,icol,irow)   =estif(j,icol,irow)   -
     &      (derxy2(j,irn,1)+auxr)*aux
              estif(j,icol,irow+1)   =estif(j,icol,irow+1) -
     &      derxy2(j,irn,4)*aux
              estif(j,icol,irow+2)   =estif(j,icol,irow+2) -
     &      derxy2(j,irn,5)*aux

              estif(j,icol+1,irow)   =estif(j,icol+1,irow)  -
     &      derxy2(j,irn,4)*aux
              estif(j,icol+1,irow+1) =estif(j,icol+1,irow+1)-
     &      (derxy2(j,irn,2)+auxr)*aux
              estif(j,icol+1,irow+2) =estif(j,icol+1,irow+2)-
     &      (derxy2(j,irn,6))*aux

              estif(j,icol+2,irow)   =estif(j,icol+2,irow)  -
     &      derxy2(j,irn,5)*aux
              estif(j,icol+2,irow+1) =estif(j,icol+2,irow+1)-
     &      (derxy2(j,irn,6))*aux
              estif(j,icol+2,irow+2) =estif(j,icol+2,irow+2)-
     &      (derxy2(j,irn,3)+auxr)*aux
            enddo

            irow = irow + 4
          enddo
          icol = icol + 4
        enddo


C--------------------------------------------------------------------------
c   Calculate viscous stabilisation part Nr(u) for higher order elements:
c
c    /
c   |  -/+ tau_mp  * 2 * nue * div(eps(v)) * u * grad(u_old) d_omega
c  /
C--------------------------------------------------------------------------
        if(flagvec(4).ne.0) then

          icol=1
          do icn=1,sizevec(2)
            irow=1
            do irn=1,sizevec(2)

              do j=1,sizevec(5)
                cc = fac(j) * tau(j,2) * paravec(2) * paravec(1)
                aux = funct(icn)*cc
                auxr = derxy2(j,irn,1)*derxy2(j,irn,2)*derxy2(j,irn,3)

                estif(j,icol,irow)  =  estif(j,icol,irow) -
     & (derxy2(j,irn,1)*vderxy(j,1,1)+derxy2(j,irn,4)*vderxy(j,1,2)
     &                       + derxy2(j,irn,5)*vderxy(j,1,3)
     &                       + auxr*vderxy(j,1,1))*aux
                estif(j,icol,irow+1) = estif(j,icol,irow+1) -
     &  (derxy2(j,irn,4)*vderxy(j,1,1)+derxy2(j,irn,2)*vderxy(j,1,2)
     &                       + derxy2(j,irn,6)*vderxy(j,1,3)
     &                       + auxr*vderxy(j,1,2))*aux
                estif(j,icol,irow+2) =  estif(j,icol,irow+2) -
     &  (derxy2(j,irn,5)*vderxy(j,1,1)+derxy2(j,irn,6)*vderxy(j,1,2)
     &               + derxy2(j,irn,3)*vderxy(j,1,3)
     &               + auxr*vderxy(j,1,3))*aux

                estif(j,icol+1,irow)  =  estif(j,icol+1,irow) -
     & (derxy2(j,irn,1)*vderxy(j,2,1)+derxy2(j,irn,4)*vderxy(j,2,2)
     &                 + derxy2(j,irn,5)*vderxy(j,2,3)
     &                 + auxr*vderxy(j,2,1))*aux
                estif(j,icol+1,irow+1) = estif(j,icol+1,irow+1) -
     &(derxy2(j,irn,4)*vderxy(j,2,1)+derxy2(j,irn,2)*vderxy(j,2,2)
     &                 + derxy2(j,irn,6)*vderxy(j,2,3)
     &                 + auxr*vderxy(j,2,2))*aux
                estif(j,icol+1,irow+2) = estif(j,icol+1,irow+2) -
     & (derxy2(j,irn,5)*vderxy(j,2,1)+derxy2(j,irn,6)*vderxy(j,2,2)
     &                 + derxy2(j,irn,3)*vderxy(j,2,3)
     &                 + auxr*vderxy(j,2,3))*aux

                estif(j,icol+2,irow)  =  estif(j,icol+2,irow) -
     &(derxy2(j,irn,1)*vderxy(j,3,1)+derxy2(j,irn,4)*vderxy(j,3,2)
     &                 + derxy2(j,irn,5)*vderxy(j,3,3)
     &                 + auxr*vderxy(j,3,1))*aux
                estif(j,icol+2,irow+1) =  estif(j,icol+2,irow+1) -
     &(derxy2(j,irn,4)*vderxy(j,3,1)+derxy2(j,irn,2)*vderxy(j,3,2)
     &                 + derxy2(j,irn,6)*vderxy(j,3,3)
     &                 + auxr*vderxy(j,3,2))*aux
                estif(j,icol+2,irow+2) =  estif(j,icol+2,irow+2) -
     &(derxy2(j,irn,5)*vderxy(j,3,1)+derxy2(j,irn,6)*vderxy(j,3,2)
     &                       + derxy2(j,irn,3)*vderxy(j,3,3)
     &                         + auxr*vderxy(j,3,3))*aux
              enddo

              irow = irow + 4
            enddo
            icol = irow + 4
          enddo
        endif


C--------------------------------------------------------------------------
c   Calculate viscous stabilisation part for ALE due to u_G:
C
c    /
c   |  +/- tau_mp  * 2 * nue * div(eps(v)) * u_G_old * grad(u) d_omega
c  /
C
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
     &            +  gridvint(k,3)*derxy(k,icn,3))*cc
              auxr = derxy2(k,irn,1)+derxy2(k,irn,2)+derxy2(k,irn,3)

              estif(k,icol,irow)   =estif(k,icol,irow)   -
     &      (derxy2(k,irn,1)+auxr)*aux
              estif(k,icol,irow+1)   =estif(k,icol,irow+1) -
     &      derxy2(k,irn,4)*aux
              estif(k,icol,irow+2)   =estif(k,icol,irow+2) -
     &      derxy2(k,irn,5)*aux

              estif(k,icol+1,irow)   =estif(k,icol+1,irow)  -
     &      derxy2(k,irn,4)*aux
              estif(k,icol+1,irow+1) =estif(k,icol+1,irow+1)-
     &      (derxy2(k,irn,2)+auxr)*aux
              estif(k,icol+1,irow+2) =estif(k,icol+1,irow+2)-
     &      (derxy2(k,irn,6))*aux

              estif(k,icol+2,irow)   =estif(k,icol+2,irow)  -
     &      derxy2(k,irn,5)*aux
              estif(k,icol+2,irow+1) =estif(k,icol+2,irow+1)-
     &      (derxy2(k,irn,6))*aux
              estif(k,icol+2,irow+2) =estif(k,icol+2,irow+2)-
     &      (derxy2(k,irn,3)+auxr)*aux

            enddo

            irow = irow + 4
          enddo
          icol = icol + 4
        enddo

      end if

      endif
c     endif higher order



C--------------------------------------------------------------------------
c  NOTATION:                                                              |
c    irow - row number in element matrix                                  |
c    icol - column number in element matrix                               |
c    irn  - row node: number of node considered for matrix-row            |
c    ird  - row dim.: number of spatial dimension at row node             |
c    posc - since there's only one full element stiffness matrix the      |
c           column number has to be changed!                              |
C--------------------------------------------------------------------------


C--------------------------------------------------------------------------
C    Calculate advection stabilisation:
C
C EULER:
C     /
C    |  tau_mu * u_old * grad(v) * grad(p)   d_omega
C   /
C
C ALE:
C     /
C    |  tau_mu * c * grad(v) * grad(p)   d_omega
C   /
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then

        posc = 0
        do icol=1,sizevec(2)
          irow=1
          posc= posc + 4
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              c = fac(j) * tau(j,1)
              aux = (vel2int(j,1)*derxy(j,irn,1)
     &             + vel2int(j,2)*derxy(j,irn,2)
     &             + vel2int(j,3)*derxy(j,irn,3))*c
              estif(j,posc,irow)  =estif(j,posc,irow)
     &               + derxy(j,icol,1)*aux
              estif(j,posc,irow+1)=estif(j,posc,irow+1)
     &               + derxy(j,icol,2)*aux
              estif(j,posc,irow+2)=estif(j,posc,irow+2)
     &               + derxy(j,icol,3)*aux
            enddo

            irow = irow + 4
          enddo
        enddo

      endif


C--------------------------------------------------------------------------
C   Calculate viscous stabilisation parts for higher order elements:
C
C    /
C   |  -/+ tau_mp * 2 * nue * div(eps(v)) * grad(p)  d_omega
C  /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        posc = 0
        do icol=1,sizevec(2)
          irow=1
          posc = posc + 4
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              c = fac(j) * tau(j,2) * paravec(2) * paravec(1)
              aux = derxy2(j,irn,1)+derxy2(j,irn,2)+derxy2(j,irn,3)
              estif(j,posc,irow)   =   estif(j,posc,irow)
     &                  -(derxy2(j,irn,1)*derxy(j,icol,1)
     &                  + derxy2(j,irn,4)*derxy(j,icol,2)
     &                  + derxy2(j,irn,5)*derxy(j,icol,3)
     &                  + aux*derxy(j,icol,1))*c
              estif(j,posc,irow+1) =   estif(j,posc,irow+1)
     &                  - (derxy2(j,irn,4)*derxy(j,icol,1)
     &                  + derxy2(j,irn,2)*derxy(j,icol,2)
     &                  + derxy2(j,irn,6)*derxy(j,icol,3)
     &                + aux*derxy(j,icol,2))*c
              estif(j,posc,irow+2) =   estif(j,posc,irow+2)
     &                  - (derxy2(j,irn,5)*derxy(j,icol,1)
     &                  + derxy2(j,irn,6)*derxy(j,icol,2)
     &                  + derxy2(j,irn,3)*derxy(j,icol,3)
     &              + aux*derxy(j,icol,3))*c
            enddo

            irow = irow + 4
          enddo
        enddo
      endif



C--------------------------------------------------------------------------
C  KPV
C--------------------------------------------------------------------------

C--------------------------------------------------------------------------
c   Calculate stabilisation part Nc(u):
c
c    /
c   |  - tau_mp * grad(q) * u_old * grad(u) d_omega
c  /
C--------------------------------------------------------------------------
      icol=1
      do icn=1,sizevec(2)
        posr = 0
        do irow=1,sizevec(2)
          posr = posr + 4

          do j=1,sizevec(5)
            c = fac(j) * tau(j,2)
            aux = (velint(j,1)*derxy(j,icn,1)
     &             + velint(j,2)*derxy(j,icn,2)
     &             + velint(j,3)*derxy(j,icn,3))*c
            estif(j,icol,posr)  =estif(j,icol,posr)
     &             -derxy(j,irow,1)*aux
            estif(j,icol+1,posr)=estif(j,icol+1,posr)
     &             -derxy(j,irow,2)*aux
            estif(j,icol+2,posr)=estif(j,icol+2,posr)
     &             -derxy(j,irow,3)*aux
          enddo

        enddo
        icol =icol + 4
      enddo


C--------------------------------------------------------------------------
c   Calculate stabilisation part Nr(u):
c
c    /
c   |  - tau_mp * grad(q) * u * grad(u_old) d_omega
c  /
C--------------------------------------------------------------------------
      if(flagvec(4).ne.0) then

        icol=1
        do icn=1,sizevec(2)
          posr = 0
          do irow=1,sizevec(2)
            posr = posr + 4

            do j=1,sizevec(5)
              c = fac(j) * tau(j,2)
              aux = funct(icn)*c
              estif(j,icol,posr)   = estif(j,icol,posr)
     &                  - aux*(derxy(j,irow,1)*vderxy(j,1,1)
     &                  + derxy(j,irow,2)*vderxy(j,1,2)
     &                  + derxy(j,irow,3)*vderxy(j,1,3))
              estif(j,icol+1,posr) = estif(j,icol+1,posr)
     &                  - aux*(derxy(j,irow,1)*vderxy(j,2,1)
     &                  + derxy(j,irow,2)*vderxy(j,2,2)
     &                  + derxy(j,irow,3)*vderxy(j,2,3))
              estif(j,icol+2,posr) = estif(j,icol+2,posr)
     &                  - aux*(derxy(j,irow,1)*vderxy(j,3,1)
     &                  + derxy(j,irow,2)*vderxy(j,3,2)
     &                  + derxy(j,irow,3)*vderxy(j,3,3))
            enddo

          enddo
          icol =icol + 4
        enddo
      endif


C--------------------------------------------------------------------------
c    Calculate stabilisation for ALE due to grid-velocity u_G:
c
c ALE:
c     /
c    |  tau_mp * grad(q) * u_G * grad(u) d_omega
c   /
c
c    REMARK:
c     - this term is linear inside the domain and for explicit free surface
c     - for implicit free surface this term is nonlinear (Nc), since u_G is
c       unknown at the free surface. So it depends on the nonlinear iteration
c       scheme if this term has to be taken into account:
c        - Newton: YES
c        - fixpoint-like: YES
C--------------------------------------------------------------------------
c     evaluate only for ALE
      if(flagvec(7).ne.0) then

        icol=1
        do icn=1,sizevec(2)
          posr = 0
          do irow=1,sizevec(2)
            posr = posr + 4

            do j=1,sizevec(5)
              c = fac(j) * tau(j,2)
              aux = (gridvint(j,1)*derxy(j,icn,1)
     &             + gridvint(j,2)*derxy(j,icn,2)
     &             + gridvint(j,3)*derxy(j,icn,3))*c
              estif(j,icol,posr)  = estif(j,icol,posr)
     &                            - derxy(j,irow,1)*aux
              estif(j,icol+1,posr)= estif(j,icol+1,posr)
     &                            - derxy(j,irow,2)*aux
              estif(j,icol+2,posr)= estif(j,icol+2,posr)
     &                            - derxy(j,irow,3)*aux
            enddo

          enddo
          icol =icol + 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate stabilisation part for higher order elements:
c
c    /
c   |  tau_mp * 2 * nue *grad(q) * div(eps(u)) d_omega
c  /
C--------------------------------------------------------------------------
      if(flagvec(6).ne.0) then

        icol=1
        do icn=1,sizevec(2)
          posr = 0
          do irow=1,sizevec(2)
            posr = posr + 4

            do j=1,sizevec(5)
              c = fac(j) * tau(j,2) * paravec(2)
              aux = derxy2(j,icn,1) + derxy2(j,icn,2) + derxy2(j,icn,3)
              estif(j,icol,posr)      = estif(j,icol,posr)
     &             +((derxy2(j,icn,1)+aux)*derxy(j,irow,1)
     &             + derxy2(j,icn,4)     *derxy(j,irow,2)
     &             + derxy2(j,icn,5)     *derxy(j,irow,3))*c

              estif(j,icol+1,posr) = estif(j,icol+1,posr)
     &       + ( derxy2(j,icn,4) *derxy(j,irow,1)
     &       +(derxy2(j,icn,2)+aux)*derxy(j,irow,2)
     &       + derxy2(j,icn,6)  *derxy(j,irow,3))*c

              estif(j,icol+2,posr) = estif(j,icol+2,posr)
     &             + (derxy2(j,icn,5)    *derxy(j,irow,1)
     &             + derxy2(j,icn,6)     *derxy(j,irow,2)
     &             +(derxy2(j,icn,3)+aux)*derxy(j,irow,3))*c

            enddo

          enddo
          icol =icol + 4
        enddo
      endif



C--------------------------------------------------------------------------
C Kpp
C--------------------------------------------------------------------------

C--------------------------------------------------------------------------
c   Calculate stabilisation part for matrix Kpp:
c
c    /
c   |  - tau_mp * grad(q) *grad(p) d_omega
c  /
C--------------------------------------------------------------------------
      posc = 0
      do icol=1,sizevec(2)
        posc = posc + 4
        posr = 0
        do irow=1,sizevec(2)
          posr = posr + 4

          do j=1,sizevec(5)
            c = fac(j) * tau(j,2)
            estif(j,posc,posr) = estif(j,posc,posr)
     &          -(derxy(j,irow,1)*derxy(j,icol,1)
     &          + derxy(j,irow,2)*derxy(j,icol,2)
     &          + derxy(j,irow,3)*derxy(j,icol,3))*c
          enddo

        enddo
      enddo

      return
      end subroutine







c     In this routine the stabilisation part of matrix Mvv and Mvp is calculated:
C--------------------------------------------------------------------------
C
C \brief evaluate stabilisaton part of M
C
C In this routine the stabilisation part of matrix Mvv is calculated:
C
C EULER:
C     /
C    |  -/+ tau_mu * u_old * grad(v) * u d_omega
C   /
C
C ALE:
C     /
C    |  -/+ tau_mu * c * grad(v) * u d_omega
C   /
C
C EULER/ALE:
C     /
C    |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
C   /
C
C In this routine the stabilisation part of matrix Mpv is calculated:
C
C     /
C    |  - tau_mp * grad(q) * u d_omega
C   /
C
C
C \param estif()()()     real*8   (o) element stiffness matrix
C \param velint()()      real*8   (i) vel at int point
C \param funct()         real*8   (i) shape functions
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param fac()           real*8   (i) gauss integration factor
C \param tau()()         real*8   (i) stabilisation parameter
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
      subroutine f3fcalstabm(estif , velint, funct, derxy, derxy2,
     &    fac, tau, paravec, flagvec, sizevec)

      implicit none

      real*8  paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 estif(sizevec(4),sizevec(3),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 fac(sizevec(4))
      real*8 tau(sizevec(4),3)

      integer irow,icol,irn,icn,j
      real*8 c,cc
      real*8 aux,auxc,auxr

      integer posr


C--------------------------------------------------------------------------
C    Calculate convection stabilisation part:
C
C EULER:
C     /
C    |  -/+ tau_mu * u_old * grad(v) * u d_omega
C   /
C
C ALE:
C     /
C    |  -/+ tau_mu * c * grad(v) * u d_omega
C   /
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then

        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              cc = fac(j) * tau(j,1)
              auxc = funct(icn)*cc

              aux = (velint(j,1)*derxy(j,irn,1)
     &             + velint(j,2)*derxy(j,irn,2)
     &             + velint(j,3)*derxy(j,irn,3) )*auxc

              estif(j,icol,irow)     =estif(j,icol,irow)    + aux
              estif(j,icol+1,irow+1) =estif(j,icol+1,irow+1)+ aux
              estif(j,icol+2,irow+2) =estif(j,icol+2,irow+2)+ aux
            enddo

            irow =irow + 4
          enddo
          icol =icol + 4
        enddo
      endif


C--------------------------------------------------------------------------
C   Calculate viscous stabilisation parts for higher order elements:
C
C    /
C   |  -/+ tau_mp * 2 * nue * div(eps(v)) * u  d_omega
C  /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        icol=1
        do icn=1,sizevec(2)
          irow=1
          do irn=1,sizevec(2)

            do j=1,sizevec(5)
              c = fac(j) * tau(j,2) * paravec(2) * paravec(1)
              aux = funct(icn)*c
              auxr = derxy2(j,irn,1) + derxy2(j,irn,2) + derxy2(j,irn,3)
              estif(j,icol,irow)    = estif(j,icol,irow)
     &             -(derxy2(j,irn,1) + auxr)*aux
              estif(j,icol,irow+1)  =estif(j,icol,irow+1)
     &             -derxy2(j,irn,4)*aux
              estif(j,icol,irow+2)  =estif(j,icol,irow+2)
     &             -derxy2(j,irn,5)*aux

              estif(j,icol+1,irow)  =estif(j,icol+1,irow)
     &             -derxy2(j,irn,4)*aux
              estif(j,icol+1,irow+1)=estif(j,icol+1,irow+1)
     &             -(derxy2(j,irn,2) + auxr)*aux
              estif(j,icol+1,irow+2)=estif(j,icol+1,irow+2)
     &             -derxy2(j,irn,6)*aux

              estif(j,icol+2,irow)  =estif(j,icol+2,irow)
     &             -derxy2(j,irn,5)*aux
              estif(j,icol+2,irow+1)=estif(j,icol+2,irow+1)
     &             -derxy2(j,irn,6)*aux
              estif(j,icol+2,irow+2)=estif(j,icol+2,irow+2)
     &             -(derxy2(j,irn,3) + auxr)*aux
            enddo

            irow =irow + 4
          enddo
          icol =icol + 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate stabilisation part for matrix Mpv:
c
c    /
c   |  - tau_mp * grad(q) * u d_omega
c  /
C--------------------------------------------------------------------------
      icol=1
      do icn=1,sizevec(2)
        posr = 0
        do irow=1,sizevec(2)
          posr = posr + 4

          do j=1,sizevec(5)
            c = fac(j) * tau(j,2)
            auxc = funct(icn)*c
            estif(j,icol,posr)  =estif(j,icol,posr)
     &             -derxy(j,irow,1)*auxc
            estif(j,icol+1,posr)=estif(j,icol+1,posr)
     &             -derxy(j,irow,2)*auxc
            estif(j,icol+2,posr)=estif(j,icol+2,posr)
     &             -derxy(j,irow,3)*auxc
          enddo

        enddo
        icol = icol + 4
      enddo

      return
      end subroutine

