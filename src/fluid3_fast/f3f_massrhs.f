C--------------------------------------------------------------------------
C
C \file
C \brief routine to evaluate rhs on element base
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
C This routine performs the multiplication emass*hist and writes its
C result on the elemental rhs force vector. This gives the right hand
C side for fluids.
C
C \param emass()()()     real*8   (o) element mass matrix
C \param ehist()()()     real*8   (i) nodal accelarations
C \param edeadng()()     real*8   (i) body force (constant in space)
C \param eiforce()()     real*8   (i) element iteration force
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fmassrhs(emass,ehist,edeadng,eiforce,hasext,
     &                      facl,sizevec)

      integer sizevec(6)
      integer hasext(sizevec(4))
      real*8 emass(sizevec(4),sizevec(3),sizevec(3))
      real*8 ehist(sizevec(4),sizevec(1),3)
      real*8 eiforce(sizevec(4),sizevec(3))
      real*8 edeadng(sizevec(4),3)
      real*8 facl

      integer i,j,k
      integer dim
      
      do i=1,sizevec(5)
        if(hasext(i).ne.0) then
          do j=1,sizevec(2)       ! count over nodes
            ehist(i,j,1) = ehist(i,j,1) + facl * edeadng(i,1)
            ehist(i,j,2) = ehist(i,j,2) + facl * edeadng(i,2)
            ehist(i,j,3) = ehist(i,j,3) + facl * edeadng(i,3)
          enddo
        else 
          hasext(i) = 1 ! this assures rhs assembly
        end if
      enddo

c     perform multiplication
      dim=4

      do i=1,sizevec(2)
        do j=1,sizevec(2)
          do k=1,sizevec(5)
            eiforce(k,(i-1)*dim+1) = eiforce(k,(i-1)*dim+1)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+1) * ehist(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+1) * ehist(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+1) * ehist(k,j,3))
            eiforce(k,(i-1)*dim+2) = eiforce(k,(i-1)*dim+2)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+2) * ehist(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+2) * ehist(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+2) * ehist(k,j,3))
            eiforce(k,(i-1)*dim+3) = eiforce(k,(i-1)*dim+3)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+3) * ehist(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+3) * ehist(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+3 )* ehist(k,j,3))
            eiforce(k,(i-1)*dim+4) = eiforce(k,(i-1)*dim+4)
     &              +(emass(k,(j-1)*dim+1,(i-1)*dim+4) * ehist(k,j,1)
     &              + emass(k,(j-1)*dim+2,(i-1)*dim+4) * ehist(k,j,2)
     &              + emass(k,(j-1)*dim+3,(i-1)*dim+4) * ehist(k,j,3))

          enddo
        enddo
      enddo


      return
      end subroutine

