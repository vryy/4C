C--------------------------------------------------------------------------
C
C \file
C \brief routine to evaluate iteration rhs on element base resulting from
C        mass and acceleration in generalised alpha time integration
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
C \brief calculates mass-acceleration-part of the rhs vector for
C        generalised-alpha time integration
C
C This routine performs the multiplication emass*acc(n) and writes its
C negative result on the elemental iteration force vector. This gives a
C part of the rhs required in the generalised alpha time integration scheme
C for fluids.
C
C \param emass()()()     real*8   (o) element mass matrix
C \param eaccn()()()     real*8   (i) nodal accelarations
C \param eiforce()()     real*8   (i) element iteration force
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fmassrhs(emass, eaccn, eiforce, sizevec)

      integer sizevec(6)
      real*8 emass(sizevec(4),sizevec(3),sizevec(3))
      real*8 eaccn(sizevec(4),sizevec(1),3)
      real*8 eiforce(sizevec(4),sizevec(3))

      integer i,j,k
      integer dim

c     perform multiplication
      dim=4

      do i=1,sizevec(2)
        do j=1,sizevec(2)
          do k=1,sizevec(5)
            eiforce(k,(i-1)*dim+1) = eiforce(k,(i-1)*dim+1)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+1) * eaccn(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+1) * eaccn(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+1) * eaccn(k,j,3))
            eiforce(k,(i-1)*dim+2) = eiforce(k,(i-1)*dim+2)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+2) * eaccn(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+2) * eaccn(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+2) * eaccn(k,j,3))
            eiforce(k,(i-1)*dim+3) = eiforce(k,(i-1)*dim+3)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+3) * eaccn(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+3) * eaccn(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+3 )* eaccn(k,j,3))
            eiforce(k,(i-1)*dim+4) = eiforce(k,(i-1)*dim+4)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+4) * eaccn(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+4) * eaccn(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+4) * eaccn(k,j,3))

          enddo
        enddo
      enddo


      return
      end subroutine

