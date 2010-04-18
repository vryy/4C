C--------------------------------------------------------------------------
C
C \file
C \brief external RHS for fluid3_fast element
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
C \brief galerkin part of external forces for vel dofs
C
C In this routine the galerkin part of the time forces for vel dofs
C is calculated:
C
C                    /
C    + (1-THETA)*dt |  v * b_old   d_omega
C                  /
C
C                /
C    + THETA*dt |  v * b   d_omega
C
C
C \param eforce()()      real*8   (o) element time force
C \param funct()         real*8   (i) shape functions
C \param edeadn()()      real*8   (i) dead load at n
C \param edeadng()()     real*8   (i) dead load at n+g
C \param fac()           real*8   (i) gauss integration factor
C \param facl            real*8   (i)
C \param facr            real*8   (i)
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fcalgalexf(eforce, funct, edeadn, edeadng,
     &          fac, facl, facr, sizevec)

      implicit none

      integer sizevec(6)
      real*8 eforce(sizevec(4),sizevec(3))
      real*8 funct(sizevec(1))
      real*8 edeadn(sizevec(4),3)
      real*8 edeadng(sizevec(4),3)
      real*8 fac(sizevec(4))
      real*8 facl
      real*8 facr

      real*8 facsl, facsr
      integer inode,irow,j



C--------------------------------------------------------------------------
C   Calculate galerkin part of external forces:
C
C                   /
C   + (1-THETA)*dt |  v * b_old   d_omega
C                 /
C
C               /
C   + THETA*dt |  v * b   d_omega
C             /
C
C--------------------------------------------------------------------------
      irow=1

      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          facsl = fac(j) * facl
          facsr = fac(j) * facr
          eforce(j,irow)=eforce(j,irow)+funct(inode)
     &        *(edeadn(j,1)*facsr+edeadng(j,1)*facsl)
          eforce(j,irow+1)=eforce(j,irow+1)+funct(inode)
     &        *(edeadn(j,2)*facsr+edeadng(j,2)*facsl)
          eforce(j,irow+2)=eforce(j,irow+2)+funct(inode)
     &        *(edeadn(j,3)*facsr+edeadng(j,3)*facsl)
        enddo

        irow=irow+4
      enddo

      return
      end subroutine





C--------------------------------------------------------------------------
C
C \brief stabilisation part of external forces for vel dofs
C
C In this routine the stabilisation part of the time forces for vel dofs
C is calculated:
C
C EULER:
C                      /
C    + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
C                    /
C
C ALE:
C                      /
C    + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
C                    /
C
C EULER/ALE:
C                        /
C    -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b^  d_omega
C                      /
C
C preassure dofs:
C ---------------
C                      /
C    - thetas(l,r)*dt  |  tau_mp * grad(q) * b^  d_omega
C                    /
C
C
C This routine is called twice with different values:
C 1. values are added to Iteration RHS (evaluation at n+1):
C     thetas(l,r) = THETA*dt
C     b^ = b = deadng
C 2. values are added to Time RHS (evaluation at n):
C    thetas(l,r) = (1-THETA)*dt
C    b^ = b_old = deadn
C
C \param eforce()()      real*8   (o) element time force
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param edeadn()()      real*8   (i) dead load at n
C \param velint()()      real*8   (i) vel at int point
C \param fac()           real*8   (i) gauss integration factor
C \param ths             real*8   (i)
C \param thp             real*8   (i)
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
      subroutine f3fcalstabexf(eforce,derxy,derxy2,edead,
     &          velint,fac,ths,thp,tau,
     &          paravec,flagvec,sizevec)

      implicit none

      real*8  paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 eforce(sizevec(4),sizevec(3))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 edead(sizevec(4),3)
      real*8 velint(sizevec(4),3)
      real*8 fac(sizevec(4))
      real*8 ths
      real*8 thp
      real*8 tau(sizevec(4),3)

      integer irow,inode,j
      real*8  aux
      real*8 c,fvts
      real*8 fact(3)

      integer iel3

c     set some factors
      c=ths


C--------------------------------------------------------------------------
C    Calculate external/convective stab-forces of time force vector:
C EULER:
C                      /
C    + thetas(l,r)*dt |  tau_mu * u * grad(v) * b^   d_omega
C                    /
C ALE:
C                      /
C    + thetas(l,r)*dt |  tau_mu * c * grad(v) * b^   d_omega
C                    /
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fact(1) = edead(j,1)*fac(j)*tau(j,1)*c
            fact(2) = edead(j,2)*fac(j)*tau(j,1)*c
            fact(3) = edead(j,3)*fac(j)*tau(j,1)*c

            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)

            eforce(j,irow) = eforce(j,irow) + aux*fact(1)
            eforce(j,irow+1) = eforce(j,irow+1) + aux*fact(2)
            eforce(j,irow+2) = eforce(j,irow+2) + aux*fact(3)
          enddo

          irow=irow+4
        enddo
      endif


C--------------------------------------------------------------------------
C   Calculate external/viscous stab-forces of time force vector:
C                       /
C   -/+ thetas(l,r)*dt |  tau_mp * 2*nue * div( eps(v) ) * b  d_omega
C                     /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        irow = 1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fvts = fac(j)*paravec(2)*tau(j,2)*paravec(1)*c
            aux = derxy2(j,1,inode)+derxy2(j,2,inode)+derxy2(j,3,inode)
            eforce(j,irow)   =  eforce(j,irow)
     &                  -( derxy2(j,inode,1)*edead(j,1)
     &                   + derxy2(j,inode,4)*edead(j,2)
     &      + derxy2(j,inode,5)*edead(j,3) + aux*edead(j,1))*fvts
            eforce(j,irow+1) = eforce(j,irow+1)
     &                  -( derxy2(j,inode,4)*edead(j,1)
     &                   + derxy2(j,inode,2)*edead(j,2)
     &      + derxy2(j,inode,6)*edead(j,3) + aux*edead(j,2))*fvts
            eforce(j,irow+2) = eforce(j,irow+2)
     &                  -( derxy2(j,inode,5)*edead(j,1)
     &                   + derxy2(j,inode,6)*edead(j,2)
     &      + derxy2(j,inode,3)*edead(j,3) + aux*edead(j,3))*fvts
          enddo

          irow = irow + 4
        enddo

      endif


C--------------------------------------------------------------------------
C     pre dofs
C--------------------------------------------------------------------------


c     set some factors
      iel3 = 4
      c=thp

C--------------------------------------------------------------------------
C   Calculate inertia/pressure stab forces of time force vector:
C                      /
C   - thetas(l,r)*dt  |  tau_mp * grad(q) * b  d_omega
C                    /
C--------------------------------------------------------------------------
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          fact(1) = edead(j,1)*tau(j,2)*fac(j)*c
          fact(2) = edead(j,2)*tau(j,2)*fac(j)*c
          fact(3) = edead(j,3)*tau(j,2)*fac(j)*c

          eforce(j,iel3)=eforce(j,iel3)-(derxy(j,inode,1)*fact(1)
     &    + derxy(j,inode,2)*fact(2)+derxy(j,inode,3)*fact(3))
        enddo

        iel3 =iel3 + 4
      enddo

      return
      end subroutine
