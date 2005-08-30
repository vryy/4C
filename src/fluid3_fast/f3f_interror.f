C--------------------------------------------------------------------------
C
C \file
C \brief evaluate an error norm
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
C \brief evaluate an error norm at a gaus point
C
C
C \param elecord()()()   real*8   (i) nodal coordinates
C \param xyzint()()      real*8   (i) coordinates of gauss point
C \param evelng()()()    real*8   (i) vel at nodes
C \param velint()()      real*8   (i) vel at int point
C \param epren()()       real*8   (i) pre at nodes
C \param preint()        real*8   (i) pre at int point
C \param velerr          real*8   (i) norm of vel error
C \param preerr          real*8   (i) norm of pre error
C \param velsol          real*8   (i) norm of vel sol
C \param presol          real*8   (i) norm of pre sol
C \param norm            integer  (i) flag for the norm
C \param fac()           real*8   (i) some flags
C \param paravec(4)      real*8   (i) some flags
C \param sizevec(6)      integer  (i) some sizes
C
C \return integer
C
C \author mn
C \date   08/05
C
C--------------------------------------------------------------------------
      subroutine f3finterr(elecord, xyzint, evelng, velint, epren,
     &      preint, velerr, preerr, velsol, presol, norm,
     &      fac, paravec, sizevec)

      implicit none

      real*8  paravec(4)
      integer sizevec(6)

      real*8 elecord(sizevec(4),sizevec(1),3)
      real*8 xyzint(sizevec(4),3)
      real*8 evelng(sizevec(4),sizevec(1),3)
      real*8 velint(sizevec(4),3)
      real*8 epren(sizevec(4),sizevec(1))
      real*8 preint(sizevec(4))
      real*8 velerr
      real*8 preerr
      real*8 velsol
      real*8 presol
      integer norm
      real*8 fac(sizevec(4))


      integer l

      real*8 visc

      real*8 solx, soly, solz, solp
      real*8 difx, dify, difz, difp
      real*8 a,d,t,x1,x2,x3, pi
      real*8 init

      visc = paravec(2)
      init = paravec(3)

      pi   = 3.14159265359
      a    = pi/4.0
      d    = pi/2.0
      t    = paravec(1)




c     loop elements
      do l=1,sizevec(5)

        x1   = xyzint(l,1)
        x2   = xyzint(l,2)
        x3   = xyzint(l,3)



        if (init .eq. 8.0) then
          solp = -a*a/2.0*(exp(2.0*a*x1) + exp(2.0*a*x2) + exp(2.0*a*x3)
     &       + 2.0*sin(a*x1 + d*x2) * cos(a*x3 + d*x1) * exp(a*(x2+x3))
     &       + 2.0*sin(a*x2 + d*x3) * cos(a*x1 + d*x2) * exp(a*(x3+x1))
     &       + 2.0*sin(a*x3 + d*x1) * cos(a*x2 + d*x3) * exp(a*(x1+x2)))
     &       * exp(-2.0*visc*d*d*t)

          solx   = -a * ( exp(a*x1) * sin(a*x2 + d*x3)
     &           + exp(a*x3) * cos(a*x1 + d*x2) ) * exp(-visc*d*d*t)

          soly   = -a * ( exp(a*x2) * sin(a*x3 + d*x1)
     &           + exp(a*x1) * cos(a*x2 + d*x3) ) * exp(-visc*d*d*t)

          solz   = -a * ( exp(a*x3) * sin(a*x1 + d*x2)
     &           + exp(a*x2) * cos(a*x3 + d*x1) ) * exp(-visc*d*d*t)


          difx = velint(l,1) - solx
          dify = velint(l,2) - soly
          difz = velint(l,3) - solz
          difp = preint(l)   - solp
        else
          return 1
        endif


        if (norm .eq. 0) then
c         infinity error
          velerr = max( abs(difx) , velerr)
          velerr = max( abs(dify) , velerr)
          velerr = max( abs(difz) , velerr)
          preerr = max( abs(difp) , preerr)

          velsol = max( abs(solx) , velsol)
          velsol = max( abs(soly) , velsol)
          velsol = max( abs(solz) , velsol)
          presol = max( abs(solp) , presol)


          elseif (norm .eq. 1) then
c         L1 norm
          velerr = velerr + (abs(difx)+abs(dify)+abs(difz)) * fac(l)
          preerr = preerr + abs(difp) * fac(l)

          velsol = velsol + (abs(solx)+abs(soly)+abs(solz)) * fac(l)
          presol = presol + abs(solp) * fac(l)


          elseif (norm .eq. 2) then
c         L2 norm
          velerr = velerr + (difx*difx+dify*dify+difz*difz) * fac(l)
          preerr = preerr + difp*difp * fac(l)

          velsol = velsol + (solx*solx+soly*soly+solz*solz) * fac(l)
          presol = presol + solp*solp * fac(l)

        endif

      enddo



      return 0
      end subroutine


