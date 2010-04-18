C--------------------------------------------------------------------------
C
C \file
C \brief time RHS for fluid3_fast element
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

c     In this routine the galerkin part of the time forces for vel dofs
c     is calculated:

c     In this routine the galerkin part of the time forces for pre dofs
c     is calculated:

c     In this routine the stabilisation part of the time forces for vel dofs
c     is calculated:

c     In this routine the stabilisation part of the time forces for vel dofs
c     is calculated:

C--------------------------------------------------------------------------
C
C \brief time forces
C
C In this routine the galerkin part of the time forces for vel dofs
C is calculated:
C
C EULER/ALE:
C       /
C    + |  v * u     d_omega
C     /
C
C EULER:
C                       /
C    (-) (1-THETA)*dt  |  v * u * grad(u)     d_omega
C                     /
C
C ALE:
C                       /
C    (-) (1-THETA)*dt  |  v * c * grad(u)     d_omega
C                     /
C
C EULER/ALE:
C                       /
C    (-) (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
C                     /
C
C                     /
C    + (1-THETA)*dt  |  div(v) * p  d_omega
C                   /
C
C In this routine the galerkin part of the time forces for pre dofs
C is calculated:
C
C                     /
C    + (1-THETA)*dt  |  q * div(u)  d_omega
C                   /
C
C In this routine the stabilisation part of the time forces for vel dofs
C is calculated:
C
C EULER:
C       /
C    + |  tau_mu * u * grad(v) * u  d_omega
C     /
C
C ALE:
C       /
C    + |  tau_mu * c * grad(v) * u  d_omega
C     /
C
C EULER/ALE:
C         /
C    -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
C       /
C
C EULER:
C                      /
C    (-) (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
C                    /
C
C ALE:
C                      /
C    (-) (1-THETA)*dt |  tau_mu * c * grad(v) * c * grad(u) d_omega
C                    /
C
C EULER:
C                      /
C    +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
C                    /
C
C ALE:
C                      /
C    +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
C                    /
C
C EULER:
C                    /
C    + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
C                  /
C
C ALE:
C                    /
C    + (1-THETA)*dt |  tau_mu * 2*nue * c *grad(v) * div( eps(u) )  d_omega
C                  /
C
C EULER/ALE:
C                      /
C    -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
C                    /
C
C        /
C   (-) |  tau_c * div(v) * div(u)  d_omega
C      /
C
C EULER:
C                     /
C   (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
C                   /
C
C ALE:
C                     /
C   (-) (1-THETA)*dt | tau_mu * c * grad(v) * grad(p)   d_omega
C                   /
C
C EULER/ALE:
C                     /
C   -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
C                   /
C
C In this routine the stabilisation part of the time forces for pre dofs
C is calculated:
C
C EULER/ALE:
C           /
C    (-)   |  tau_mp * grad(q) * u  d_omega
C         /
C
C EULER:
C                      /
C    + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
C                    /
C
C ALE:
C                      /
C    + (1-THETA)*dt   |  tau_mp * grad(q) * c * grad(u)  d_omega
C                    /
C
C EULER/ALE:
C                      /
C    + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
C                    /
C
C
C \param eforce()()      real*8   (o) element time force
C \param velint()()      real*8   (i) vel at int point
C \param vel2int()()     real*8   (i) vel at int point
C \param covint()()      real*8   (i) conv. vel at int point
C \param funct()         real*8   (i) shape functions
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (i) 2nd global derivatives
C \param vderxy()()()    real*8   (i) global vel derivatives
C \param vderxy2()()()   real*8   (i) 2nd global vel derivatives
C \param pderxy()()      real*8   (i) global pressure derivatives
C \param preint()        real*8   (i) pres at integration point
C \param fac()           real*8   (i) gauss integration factor
C \param ths             real*8   (i) factor
C \param thp             real*8   (i) factor
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
      subroutine f3fcaltf(eforce,velint,vel2int,covint,funct,
     &    derxy,derxy2,vderxy,vderxy2,pderxy,preint,fac,ths,thp,tau,
     &    paravec, flagvec, sizevec)

      implicit none

      real*8  paravec(2)
      integer flagvec(7)
      integer sizevec(6)

      real*8 eforce(sizevec(4),sizevec(3))
      real*8 velint(sizevec(4),3)
      real*8 vel2int(sizevec(4),3)
      real*8 covint(sizevec(4),3)
      real*8 funct(sizevec(1))
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 vderxy(sizevec(4),3,3)
      real*8 vderxy2(sizevec(4),6,3)
      real*8 pderxy(sizevec(4),3)
      real*8 preint(sizevec(4))
      real*8 fac(sizevec(4))
      real*8 ths
      real*8 thp
      real*8 tau(sizevec(4),3)

      integer irow,inode,k
      real*8 c
      real*8 aux
      real*8 facsr
      real*8 facpr
      real*8 fact(3)

      integer j
      integer iel3

      real*8 cc
      real*8 fvts,fvtsr,fvvtsr



C--------------------------------------------------------------------------
c   Calculate intertia forces of time force vector:
c
c      /
c   + |  v * u     d_omega
c    /
C--------------------------------------------------------------------------
      irow=1
      do inode=1,sizevec(2)

        do k=1,sizevec(4)
          fact(1) = vel2int(k,1)*fac(k)
          fact(2) = vel2int(k,2)*fac(k)
          fact(3) = vel2int(k,3)*fac(k)
          eforce(k,irow) =eforce(k,irow)+funct(inode)*fact(1)
          eforce(k,irow+1) =eforce(k,irow+1)+funct(inode)*fact(2)
          eforce(k,irow+2) =eforce(k,irow+2)+funct(inode)*fact(3)
        enddo

        irow=irow+4
      enddo


C--------------------------------------------------------------------------
c    Calculate convective forces of time force vector:
c
c EULER:
c                     /
c    - (1-THETA)*dt  |  v * u * grad(u)     d_omega
c                   /
c
c ALE:
c                     /
c    - (1-THETA)*dt  |  v * c * grad(u)     d_omega
c                   /
C--------------------------------------------------------------------------
      irow=1
      do inode=1,sizevec(2)

        do k=1,sizevec(5)
          facsr= fac(k) * ths
          eforce(k,irow)=eforce(k,irow)
     &          -funct(inode)*covint(k,1)*facsr
          eforce(k,irow+1)=eforce(k,irow+1)
     &          -funct(inode)*covint(k,2)*facsr
          eforce(k,irow+2)=eforce(k,irow+2)
     &          -funct(inode)*covint(k,3)*facsr
        enddo

        irow=irow+4
      enddo


C--------------------------------------------------------------------------
c   Calculate viscous forces of time force vector:
c
c                    /
c   - (1-THETA)*dt  |  2*nue * eps(v) : eps(u)  d_omega
c                  /
C--------------------------------------------------------------------------
      irow=1
      do inode=1,sizevec(2)

        do k=1,sizevec(5)
          c = fac(k) * ths * paravec(2)
          eforce(k,irow) =eforce(k,irow)-(derxy(k,inode,1)
     &          *(vderxy(k,1,1)+vderxy(k,1,1))*c)
          eforce(k,irow) =eforce(k,irow)-(derxy(k,inode,2)
     &          *(vderxy(k,2,1)+vderxy(k,1,2))*c)
          eforce(k,irow) =eforce(k,irow)-(derxy(k,inode,3)
     &          *(vderxy(k,3,1)+vderxy(k,1,3))*c)

          eforce(k,irow+1) =eforce(k,irow+1)-(derxy(k,inode,1)
     &          *(vderxy(k,1,2)+vderxy(k,2,1))*c)
          eforce(k,irow+1) =eforce(k,irow+1)-(derxy(k,inode,2)
     &          *(vderxy(k,2,2)+vderxy(k,2,2))*c)
          eforce(k,irow+1) =eforce(k,irow+1)-(derxy(k,inode,3)
     &          *(vderxy(k,3,2)+vderxy(k,2,3))*c)

          eforce(k,irow+2) =eforce(k,irow+2)-(derxy(k,inode,1)
     &          *(vderxy(k,1,3)+vderxy(k,3,1))*c)
          eforce(k,irow+2) =eforce(k,irow+2)-(derxy(k,inode,2)
     &          *(vderxy(k,2,3)+vderxy(k,3,2))*c)
          eforce(k,irow+2) =eforce(k,irow+2)-(derxy(k,inode,3)
     &          *(vderxy(k,3,3)+vderxy(k,3,3))*c)
        enddo

        irow=irow+4
      enddo


C--------------------------------------------------------------------------
c   Calculate pressure forces of time force vector:
c
c                    /
c   + (1-THETA)*dt  |  div(v) * p  d_omega
c                  /
C--------------------------------------------------------------------------
      irow=1
      do inode=1,sizevec(2)

        do k=1,sizevec(5)
          facpr = fac(k) * thp
          aux = preint(k) * facpr
          eforce(k,irow)=eforce(k,irow)+derxy(k,inode,1)*aux
          eforce(k,irow+1)=eforce(k,irow+1)+derxy(k,inode,2)*aux
          eforce(k,irow+2)=eforce(k,irow+2)+derxy(k,inode,3)*aux
        enddo

        irow=irow+4
      enddo




C--------------------------------------------------------------------------
c     galerkin part of time forces for pre dofs
C--------------------------------------------------------------------------

c     set some factors
      iel3=4

C--------------------------------------------------------------------------
c   Calculate continuity forces of time force vector:
c
c                    /
c   + (1-THETA)*dt  |  q * div(u)  d_omega
c                  /
C--------------------------------------------------------------------------
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          facsr = fac(j) * ths
          aux = facsr*(vderxy(j,1,1)+vderxy(j,2,2)+vderxy(j,3,3))
          eforce(j,iel3) =eforce(j,iel3)+funct(inode)*aux
        enddo

        iel3 =iel3 + 4
      enddo


C--------------------------------------------------------------------------
c\brief stabilisation part of time forces for vel dofs
C--------------------------------------------------------------------------

C--------------------------------------------------------------------------
c    Calculate inertia/convective stab-forces of time force vector:
c
c EULER:
c       /
c    + |  tau_mu * u * grad(v) * u  d_omega
c     /
c
c ALE:
c       /
c    + |  tau_mu * c * grad(v) * u  d_omega
c     /
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fact(1) = vel2int(j,1)*fac(j)*tau(j,1)
            fact(2) = vel2int(j,2)*fac(j)*tau(j,1)
            fact(3) = vel2int(j,3)*fac(j)*tau(j,1)
            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)

            eforce(j,irow)=eforce(j,irow)+aux*fact(1)
            eforce(j,irow+1)=eforce(j,irow+1)+aux*fact(2)
            eforce(j,irow+2)=eforce(j,irow+2)+aux*fact(3)
          enddo

          irow=irow+4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate inertia/viscous stab-forces of time force vector:
c
c        /
c   -/+ |  tau_mp * 2*nue * div( eps(v) ) * u  d_omega
c      /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        irow = 1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fvts = fac(j)*paravec(2)*tau(j,2)*paravec(1)
            aux = derxy2(j,inode,1)+derxy2(j,inode,2)+derxy2(j,inode,3)
            eforce(j,irow)   =  eforce(j,irow)  -
     &  ( derxy2(j,inode,1)*vel2int(j,1)+derxy2(j,inode,4)*vel2int(j,2)
     &  + derxy2(j,inode,5)*vel2int(j,3)+aux*vel2int(j,1))*fvts
            eforce(j,irow+1) = eforce(j,irow+1) -
     &  ( derxy2(j,inode,4)*vel2int(j,1)+derxy2(j,inode,2)*vel2int(j,2)
     &  + derxy2(j,inode,6)*vel2int(j,3)+aux*vel2int(j,2))*fvts
            eforce(j,irow+2) = eforce(j,irow+2) -
     &  ( derxy2(j,inode,5)*vel2int(j,1)+derxy2(j,inode,6)*vel2int(j,2)
     &  + derxy2(j,inode,3)*vel2int(j,3)+aux*vel2int(j,3))*fvts
          enddo

          irow = irow+4
        enddo

      endif


C--------------------------------------------------------------------------
c    Calculate convective/convective stab-forces of time force vector:
c
c EULER
c                    /
c    - (1-THETA)*dt |  tau_mu * u * grad(v) * u * grad(u) d_omega
c                  /
c
c ALE:
c                    /
c    - (1-THETA)*dt |  tau_mu * c * grad(v) * c * grad(u) d_omega
c                  /
c
C--------------------------------------------------------------------------
      if(flagvec(2).ne.0) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            facsr = fac(j) * ths
            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)
            fact(1) = tau(j,1)*covint(j,1)*facsr
            fact(2) = tau(j,1)*covint(j,2)*facsr
            fact(3) = tau(j,1)*covint(j,3)*facsr

            eforce(j,irow)   = eforce(j,irow)   - aux*fact(1)
            eforce(j,irow+1) = eforce(j,irow+1) - aux*fact(2)
            eforce(j,irow+2) = eforce(j,irow+2) - aux*fact(3)
          enddo

          irow=irow+4
        enddo

      endif

C--------------------------------------------------------------------------
c    Calculate convective/viscous stab-forces of time force vector:
c
c EULER:
c                      /
c    +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * u * grad(u)  d_omega
c                    /
c
c ALE:
c                      /
c    +/- (1-THETA)*dt |  tau_mp * 2*nue * div( eps(v) ) * c * grad(u)  d_omega
c                    /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then


        irow = 1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fvtsr =fac(j)*paravec(2)*tau(j,2)*paravec(1)*ths
            aux = derxy2(j,inode,1)+derxy2(j,inode,2)+derxy2(j,inode,3)
            eforce(j,irow)   = eforce(j,irow)
     &           +(derxy2(j,inode,1)*covint(j,1)
     &           + derxy2(j,inode,4)*covint(j,2)
     &           + derxy2(j,inode,5)*covint(j,3)+ aux*covint(j,1))*fvtsr
            eforce(j,irow+1) = eforce(j,irow+1)
     &           +(derxy2(j,inode,4)*covint(j,1)
     &           + derxy2(j,inode,2)*covint(j,2)
     &           + derxy2(j,inode,6)*covint(j,3)+ aux*covint(j,2))*fvtsr
            eforce(j,irow+2) = eforce(j,irow+2)
     &           +(derxy2(j,inode,5)*covint(j,1)
     &           + derxy2(j,inode,6)*covint(j,2)
     &           + derxy2(j,inode,3)*covint(j,3)+ aux*covint(j,3))*fvtsr
          enddo

          irow =irow + 4
        enddo

      endif

C--------------------------------------------------------------------------
c    Calculate viscous/convective stab-forces of time force vector:
c
c EULER:
c                    /
c    + (1-THETA)*dt |  tau_mu * 2*nue * u *grad(v) * div( eps(u) )  d_omega
c                  /
c
c ALE:
c                    /
c    + (1-THETA)*dt |  tau_mu * 2*nue * c *grad(v) * div( eps(u) )  d_omega
c                  /
C--------------------------------------------------------------------------
      if((flagvec(2).ne.0).and.(flagvec(6).ne.0)) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            cc = fac(j) * ths * paravec(2)*tau(j,1)

            fact(1) = (2.0*vderxy2(j,1,1) + vderxy2(j,2,1)
     &          + vderxy2(j,3,1)  + vderxy2(j,4,2) + vderxy2(j,5,3))*cc
            fact(2) = (2.0*vderxy2(j,2,2) + vderxy2(j,1,2)
     &          + vderxy2(j,3,2)  + vderxy2(j,4,1) + vderxy2(j,6,3))*cc
            fact(3) = (2.0*vderxy2(j,3,3) + vderxy2(j,1,3)
     &          + vderxy2(j,2,3)  + vderxy2(j,5,1) + vderxy2(j,6,2))*cc

            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)

            eforce(j,irow)   = eforce(j,irow) +  aux*fact(1)
            eforce(j,irow+1) = eforce(j,irow+1)+ aux*fact(2)
            eforce(j,irow+2) = eforce(j,irow+2)+ aux*fact(3)
          enddo

          irow =irow + 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate viscous/viscous stab-forces of time force vector:
c
c                     /
c   -/+ (1-THETA)*dt |  tau_mp * 4*nue^2 * div( eps(v) ) * div ( eps(u) )  d_omega
c                   /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then
        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            fvvtsr =fac(j)*paravec(2)*tau(j,2)*paravec(1)*ths*paravec(2)
            fact(1) = (2.0*vderxy2(j,1,1) + vderxy2(j,4,2)
     &      + vderxy2(j,5,3) + vderxy2(j,2,1) + vderxy2(j,3,1))*fvvtsr
            fact(2) = (2.0*vderxy2(j,2,2) + vderxy2(j,4,1)
     &      + vderxy2(j,6,3) + vderxy2(j,1,2) + vderxy2(j,3,2))*fvvtsr
            fact(3) = (2.0*vderxy2(j,3,3) + vderxy2(j,5,1)
     &      + vderxy2(j,6,2) + vderxy2(j,1,3) + vderxy2(j,2,3))*fvvtsr

            aux = derxy2(j,inode,1)+derxy2(j,inode,2)+derxy2(j,inode,3)
            eforce(j,irow)   = eforce(j,irow)
     &           -(derxy2(j,inode,1)*fact(1)
     &           + derxy2(j,inode,4)*fact(2)
     &           + derxy2(j,inode,5)*fact(3) + aux*fact(1))
            eforce(j,irow+1) = eforce(j,irow+1)
     &           -(derxy2(j,inode,4)*fact(1)
     &           + derxy2(j,inode,2)*fact(2)
     &           + derxy2(j,inode,6)*fact(3) + aux*fact(2))
            eforce(j,irow+2) = eforce(j,irow+2)
     &           -(derxy2(j,inode,5)*fact(1)
     &           + derxy2(j,inode,6)*fact(2)
     &           + derxy2(j,inode,3)*fact(3) + aux*fact(3))
          enddo

          irow = irow + 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate continuity stab-forces of time force vector:
c
c     /
c  - |  tau_c * div(v) * div(u)  d_omega
c   /
C--------------------------------------------------------------------------
      if(flagvec(1).ne.1) then
        irow = 1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            facsr = fac(j) * ths
            aux = tau(j,3) * facsr
     &       *(vderxy(j,1,1)+vderxy(j,2,2)+vderxy(j,3,3))
            eforce(j,irow)=eforce(j,irow)-derxy(j,inode,1)*aux
            eforce(j,irow+1)=eforce(j,irow+1)-derxy(j,inode,2)*aux
            eforce(j,irow+2)=eforce(j,irow+2)-derxy(j,inode,3)*aux
          enddo

          irow =irow + 4
        enddo

      endif


C--------------------------------------------------------------------------
c    Calculate pressure/convective stab-forces of time force vector:
c
c EULER:
c                     /
c   (-) (1-THETA)*dt | tau_mu * u * grad(v) * grad(p)   d_omega
c                   /
c
c ALE:
c                     /
c   (-) (1-THETA)*dt | tau_mu * c * grad(v) * grad(p)   d_omega
c                   /
C--------------------------------------------------------------------------
      if((flagvec(2).ne.0)) then
        irow = 1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            facpr = fac(j) * thp
            fact(1) = tau(j,1)*pderxy(j,1)*facpr
            fact(2) = tau(j,1)*pderxy(j,2)*facpr
            fact(3) = tau(j,1)*pderxy(j,3)*facpr
            aux = derxy(j,inode,1)*velint(j,1)
     &          + derxy(j,inode,2)*velint(j,2)
     &          + derxy(j,inode,3)*velint(j,3)

            eforce(j,irow) =eforce(j,irow)-aux*fact(1)
            eforce(j,irow+1) =eforce(j,irow+1)-aux*fact(2)
            eforce(j,irow+2) =eforce(j,irow+2)-aux*fact(3)
          enddo

          irow=irow+4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate pressure/viscous stab-forces of time force vector:
c
c                    /
c  -/+ (1-THETA)*dt | tau_mp * 2*nue * div( eps(v) ) * grad(p)    d_omega
c                  /
C--------------------------------------------------------------------------
      if((flagvec(3).ne.0).and.(flagvec(6).ne.0)) then

        irow=1
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            facpr = fac(j) * thp
            cc = facpr * paravec(2) * tau(j,2) * paravec(1)
            aux = derxy2(j,inode,1)+derxy2(j,inode,2)+derxy2(j,inode,3)
            eforce(j,irow)   =  eforce(j,irow)
     &                   +(derxy2(j,inode,1)*pderxy(j,1)
     &                   + derxy2(j,inode,4)*pderxy(j,2)
     &   + derxy2(j,inode,5)*pderxy(j,3) + aux*pderxy(j,1))*cc
            eforce(j,irow+1) =  eforce(j,irow+1)
     &                   +(derxy2(j,inode,4)*pderxy(j,1)
     &                   + derxy2(j,inode,2)*pderxy(j,2)
     &   + derxy2(j,inode,6)*pderxy(j,3) + aux*pderxy(j,2))*cc
            eforce(j,irow+2) =  eforce(j,irow+2)
     &                   +(derxy2(j,inode,5)*pderxy(j,1)
     &                   + derxy2(j,inode,6)*pderxy(j,2)
     &   + derxy2(j,inode,3)*pderxy(j,3) + aux*pderxy(j,3))*cc
          enddo

          irow =irow + 4
        enddo

      endif



C--------------------------------------------------------------------------
c     stabilisation part of time forces for pre dofs
C--------------------------------------------------------------------------

c     set some factors
      iel3 = 4

C--------------------------------------------------------------------------
c   Calculate inertia/pressure stab forces of time force vector:
c
c        /
c   -   |  tau_mp * grad(q) * u  d_omega
c      /
C--------------------------------------------------------------------------
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          fact(1) = velint(j,1)*tau(j,2)*fac(j)
          fact(2) = velint(j,2)*tau(j,2)*fac(j)
          fact(3) = velint(j,3)*tau(j,2)*fac(j)

          eforce(j,iel3)=eforce(j,iel3)-(derxy(j,inode,1)*fact(1)
     &    + derxy(j,inode,2)*fact(2) + derxy(j,inode,3)*fact(3))
        enddo

        iel3 = iel3 + 4
      enddo


C--------------------------------------------------------------------------
c    Calculate convective/pressure stab forces of time force vector:
c
c EULER:
c                      /
c    + (1-THETA)*dt   |  tau_mp * grad(q) * u * grad(u)  d_omega
c
c                    /
c ALE:
c                      /
c    + (1-THETA)*dt   |  tau_mp * grad(q) * c * grad(u)  d_omega
c                    /
C--------------------------------------------------------------------------
      iel3 = 4
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          facsr = fac(j) * ths
          fact(1) = covint(j,1)*tau(j,2)*facsr
          fact(2) = covint(j,2)*tau(j,2)*facsr
          fact(3) = covint(j,3)*tau(j,2)*facsr

          eforce(j,iel3)=eforce(j,iel3)+(derxy(j,inode,1)*fact(1)
     &    + derxy(j,inode,2)*fact(2) + derxy(j,inode,3)*fact(3))
        enddo

        iel3 = iel3 + 4
      enddo


C--------------------------------------------------------------------------
c   Calculate vsicous/pressure stab forces of time force vector:
c
c                     /
c   - (1-THETA)*dt   |  tau_mp * 2*nue * grad(q) * div( eps(u) )  d_omega
c                   /
C--------------------------------------------------------------------------
      if(flagvec(6).ne.0) then

        iel3 = 4
        do inode=1,sizevec(2)

          do j=1,sizevec(5)
            facsr = fac(j) * ths

            fact(1) = (2.0*vderxy2(j,1,1)+vderxy2(j,4,2)
     &           + vderxy2(j,5,3) + vderxy2(j,2,1)
     &           + vderxy2(j,3,1))*tau(j,2)*paravec(2)*facsr
            fact(2) = (2.0*vderxy2(j,2,2)+vderxy2(j,4,1)
     &           + vderxy2(j,6,3) + vderxy2(j,1,2)
     &           + vderxy2(j,3,2))*tau(j,2)*paravec(2)*facsr
            fact(3) = (2.0*vderxy2(j,3,3)+vderxy2(j,5,1)
     &           + vderxy2(j,6,2) + vderxy2(j,1,3)
     &           + vderxy2(j,2,3))*tau(j,2)*paravec(2)*facsr

            eforce(j,iel3) =eforce(j,iel3)-(derxy(j,inode,1)*fact(1)
     &    + derxy(j,inode,2)*fact(2)  + derxy(j,inode,3)*fact(3))
          enddo

          iel3 = iel3 + 4
        enddo

      endif


C--------------------------------------------------------------------------
c   Calculate pressure/pressure stab forces of time force vector:
c
c                     /
c   + (1-THETA)*dt   |  tau_mp * grad(q) * grad(p)  d_omega
c                   /
C--------------------------------------------------------------------------
      iel3 = 4
      do inode=1,sizevec(2)

        do j=1,sizevec(5)
          facpr = fac(j) * thp
          fact(1) = tau(j,2)*pderxy(j,1)*facpr
          fact(2) = tau(j,2)*pderxy(j,2)*facpr
          fact(3) = tau(j,2)*pderxy(j,3)*facpr

          eforce(j,iel3)=eforce(j,iel3)+(derxy(j,inode,1)*fact(1)
     &    + derxy(j,inode,2)*fact(2) + derxy(j,inode,3)*fact(3))
        enddo

        iel3 = iel3 + 4
      enddo


      return
      end subroutine

