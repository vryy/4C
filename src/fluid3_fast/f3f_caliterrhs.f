C--------------------------------------------------------------------------
C
C \file
C \brief iteration RHS for fluid3_fast element
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
C \brief iteration forces
C
C In this routine the galerkin part of the iteration forces for vel dofs
C is calculated:
C
C EULER:
C
C                    /
C    (+/-) THETA*dt |  v * u * grad(u)  d_omega
C                  /
C
C ALE:
C                    /
C    (+/-) THETA*dt |  v * c * grad(u)  d_omega
C                  /
C NOTE:
C   EULER:  covint = u*grad(u)
C   ALE:    covint = c*grad(u)
C
C
C In this routine the stabilisation part of the iteration forces for vel dofs
C is calculated:
C
C EULER:
C                    /
C    (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
C                  /
C
C                        /
C    (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
C                      /
C
C ALE:
C                    /
C    (+/-) THETA*dt |  tau_mu * c * grad(v) * c * grad(u)  d_omega
C                  /
C
C                        /
C    (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
C                      /
C
C NOTE:
C   EULER:  covint = u*grad(u)   velint = u
C   ALE:    covint = c*grad(u)   velint = alecovint (c)
C
C
C In this routine the stabilisation part of the iteration forces for pre
C dofs is calculated:
C
C EULER:
C                /
C    - THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
C              /
C
C ALE:
C                /
C    - THETA*dt |  tau_mp * grad(q) * c * grad(u)  d_omega
C              /
C
C
C \param eforce()()      real*8   (o) element iteration force
C \param covint()()      real*8   (i) conv. vel at int point
C \param velint()()      real*8   (i) vel at int point
C \param funct()         real*8   (i) shape functions
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param facsl()         real*8   (i)
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
      subroutine f3fcalif(eforce,covint,velint,funct,
     &  derxy,derxy2,facsl, tau, paravec, flagvec, sizevec)

      implicit none

      real*8  paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 eforce(sizevec(4),sizevec(3))
      real*8 covint(sizevec(4),3)
      real*8 velint(sizevec(4),3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 facsl(sizevec(4))
      real*8 tau(sizevec(4),3)

      integer inode,j
      integer irow
      real*8 cc
      real*8 aux
      real*8 fact(3)

      integer iel3

      iel3 = 4

C--------------------------------------------------------------------------
c   Calculate convective/convective stabilastion of iteration force vector:
c                   /
c   (+/-) THETA*dt |  tau_mu * u * grad(v) * u * grad(u)  d_omega
c    |            /
c    |
c    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then
        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fact(1) = tau(j,1)*covint(j,1)*facsl(j)
            fact(2) = tau(j,1)*covint(j,2)*facsl(j)
            fact(3) = tau(j,1)*covint(j,3)*facsl(j)

            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)

            eforce(j,irow)=eforce(j,irow)+aux*fact(1)
            eforce(j,irow+1)=eforce(j,irow+1)+aux*fact(2)
            eforce(j,irow+2)=eforce(j,irow+2)+aux*fact(3)
          enddo

          irow = irow + 4
        enddo

      endif
C--------------------------------------------------------------------------
c   Calculate convective/viscous stabilastion of iteration force vector:
c                       /
c   (-/+) -/+ THETA*dt |  tau_mp * 2*nue * div( eps(v) ) *  * grad(u)  d_omega
c    |                /
c    |
c    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            cc = facsl(j)*paravec(2)*tau(j,2)*paravec(1)
            aux = derxy2(j,inode,1)
     &          + derxy2(j,inode,2)
     &          + derxy2(j,inode,3)
            eforce(j,irow)   = eforce(j,irow) -
     &                  ((derxy2(j,inode,1)+aux)*covint(j,1)
     &                   + derxy2(j,inode,4)     *covint(j,2)
     &                   + derxy2(j,inode,5)     *covint(j,3))*cc
            eforce(j,irow+1) = eforce(j,irow+1) -
     &                  ((derxy2(j,inode,2)+aux)*covint(j,2)
     &                   + derxy2(j,inode,4)     *covint(j,1)
     &                   + derxy2(j,inode,6)     *covint(j,3))*cc
            eforce(j,irow+2) = eforce(j,irow+2) -
     &                  ((derxy2(j,inode,3)+aux)*covint(j,3)
     &                   + derxy2(j,inode,5)     *covint(j,1)
     &                   + derxy2(j,inode,6)     *covint(j,2))*cc
          enddo

          irow =irow+ 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate convective pressure stabilisation iteration force vector:
c                   /
c   (-/+) THETA*dt |  tau_mp * grad(q) * u * grad(u)  d_omega
c    |            /
c    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
C--------------------------------------------------------------------------

      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          fact(1) = covint(j,1)*tau(j,2)*facsl(j)
          fact(2) = covint(j,2)*tau(j,2)*facsl(j)
          fact(3) = covint(j,3)*tau(j,2)*facsl(j)
          eforce(j,iel3) = eforce(j,iel3)-(derxy(j,inode,1)*fact(1)
     &   + derxy(j,inode,2)*fact(2) + derxy(j,inode,3)*fact(3) )
        enddo

        iel3 = iel3 + 4
      enddo



C--------------------------------------------------------------------------
c galerkin part of iteration forces
C--------------------------------------------------------------------------


C--------------------------------------------------------------------------
c   Calculate convective forces of iteration force vector:
c                   /
c   (+/-) THETA*dt |  v * u * grad(u)  d_omega
c    |            /
c    |-> signs due to nonlin. iteration scheme (dynvar->sigma)
C--------------------------------------------------------------------------

      irow = 1
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          eforce(j,irow)   = eforce(j,irow)
     &                   + funct(inode)*covint(j,1)*facsl(j)
          eforce(j,irow+1) = eforce(j,irow+1)
     &                   + funct(inode)*covint(j,2)*facsl(j)
          eforce(j,irow+2) = eforce(j,irow+2)
     &                   + funct(inode)*covint(j,3)*facsl(j)
        enddo

        irow = irow + 4
      enddo

      return
      end subroutine

