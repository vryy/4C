C--------------------------------------------------------------------------
C
C \file
C \brief contains 'faddmsr' which adds a fortran estiff into a msr
C        sparse matrix
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
C \brief assembling of one element matrix
C
C
C \param b      DOUBLE  (o) explanation
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine faddmsr(estif,lm,owner,invupdate,bindx,invbindx,
     &                   val,
     &                   myrank,nprocs,numeqtotal,numeq,
     &                   numnp,nd,aloopl,loopl,ele )

      implicit none

c     some helpers
      integer  i,j,k,ii,jj

      integer  iscouple
      integer  iowner,iiind
      integer  ind

c     parameters
      integer  numnp
      integer  nd
      integer  aloopl
      integer  loopl
      integer  ele

      integer  myrank
      integer  nprocs
      integer  numeqtotal
      integer  numeq

      real*8   estif(loopl,numnp,numnp)

      integer  lm(numnp)
      integer  owner(numnp)

      integer  invupdate(numeqtotal)
      integer  bindx(*)
      integer  invbindx(*)

      real*8   val(*)



c      /* now start looping the dofs */
c      /* loop over i (the element row) */
      iscouple  = 0
      iowner    = myrank

      do i=1,nd


        ii = lm(i)

c       /* check for boundary condition */
        if (ii.gt.numeqtotal) then
          goto 100
        end if


        iiind = invupdate(ii)
CC        if (iiind==-1) dserror("dof ii not found on this proc");


c       /* NO initialization: only values that are set below will be used for this row!! */
c       /* fill invbindx */
        do k=bindx(iiind+0)+1,bindx(iiind+1)
          invbindx(bindx(k)+1) = k
        end do

c       /* loop over j (the element column) */
c       /* This is the full unsymmetric version ! */
        do j=1,nd


          jj = lm(j)


c         /* check for boundary condition */
          if (jj.gt.numeqtotal) then
            goto 80
          end if


c         /* do main-diagonal entry */
          if (ii.eq.jj) then
            val(iiind) = val(iiind) + estif(ele+1,j,i)

          else
c         /* do off-diagonal entry in row ii */
            ind = invbindx(jj)
            val(ind) = val(ind) + estif(ele+1,j,i)
          end if



   80     continue
        end do ! /* end loop over j */


  100   continue
      end do !  /* end loop over i */


 999  return
      end subroutine





C--------------------------------------------------------------------------
C
C \brief assembling of one element matrix for parallel case
C
C
C \param b      DOUBLE  (o) explanation
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine faddmsrp(estif,lm,owner,invupdate,bindx,invbindx,
     &                   val,
     &                   myrank,nprocs,numeqtotal,numeq,
     &                   numnp,nd,aloopl,loopl,ele )

      implicit none

c     some helpers
      integer  i,j,k,ii,jj

      integer  iscouple
      integer  iowner,iiind
      integer  ind

c     parameters
      integer  numnp
      integer  nd
      integer  aloopl
      integer  loopl
      integer  ele

      integer  myrank
      integer  nprocs
      integer  numeqtotal
      integer  numeq

      real*8   estif(loopl,numnp,numnp)

      integer  lm(numnp)
      integer  owner(numnp)

      integer  invupdate(numeqtotal)
      integer  bindx(*)
      integer  invbindx(*)

      real*8   val(*)
CC      integer  nsend
CC      integer  isend(*)
CC      real*8   dsend(*)



c      /* now start looping the dofs */
c      /* loop over i (the element row) */
      iscouple  = 0
      iowner    = myrank

      do i=1,nd


        ii = lm(i)

C       ifdef PARALLEL
C       loop only my own rows
        if (owner(i).ne.myrank) goto 100
C       endif

c       /* check for boundary condition */
        if (ii.gt.numeqtotal) then
          goto 100
        end if

CCc       /* check for coupling condition */
CC        #ifdef PARALLEL
CC        if (ncdofs)
CC            {
CC              iscouple = 0;
CC              iowner    = -1;
CC
CC              for (k=0; k<ncdofs; k++)
CC              {
CC                if (ii==cdofs[k][0])
CC                {
CC                  iscouple=1;
CC                  for (i=1; i<=nprocs; i++)
CC                  {
CC                    if (cdofs[k][i]==2)
CC                    {
CC                      iowner=i-1;
CC                      break;
CC                    }
CC                  }
CC                }
CC              }
CC
CC            }
CC#endif


        iiind = invupdate(ii)
CC        if (iiind==-1) dserror("dof ii not found on this proc");


c       /* NO initialization: only values that are set below will be used for this row!! */
c       /* fill invbindx */
        do k=bindx(iiind+0)+1,bindx(iiind+1)
          invbindx(bindx(k)+1) = k
        end do

c       /* loop over j (the element column) */
c       /* This is the full unsymmetric version ! */
        do j=1,nd


          jj = lm(j)


c         /* check for boundary condition */
          if (jj.gt.numeqtotal) then
            goto 80
          end if

c         /* either not a coupled dof or I am master owner */
          if (iscouple.eq.0 .or. iowner.eq.myrank) then

c           /* do main-diagonal entry */
            if (ii.eq.jj) then
              val(iiind) = val(iiind) + estif(ele+1,j,i)

            else
c           /* do off-diagonal entry in row ii */
              ind = invbindx(jj)
              val(ind) = val(ind) + estif(ele+1,j,i)
            end if

c         /* a coupled dof and I am slave owner */
          else
C            add_msr_sendbuff(ii,jj,i,j,iowner,isend,dsend,estif,nsend);
          end if


   80     continue
        end do ! /* end loop over j */


  100   continue
      end do !  /* end loop over i */


 999  return
      end subroutine
