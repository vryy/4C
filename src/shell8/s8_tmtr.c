/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates all metrics                                m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_tmtr(DOUBLE   **x,
                DOUBLE   **a3,
                DOUBLE     e3,
                DOUBLE   **gkov,
                DOUBLE   **gkon,
                DOUBLE   **gmkov,
                DOUBLE   **gmkon,
                DOUBLE    *det,
                DOUBLE    *funct,
                DOUBLE   **deriv,
                INT        iel,
                DOUBLE     condfac,
                INT        flag)
{
INT    i,j,k,idim,ialpha,inode;
INT    i1=0,i2=0,i3=0,i4=0,i5=0;
DOUBLE zeta,det_dummy;
#ifdef DEBUG
dstrc_enter("s8_tmtr");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- make decision about what to do */
if (flag==0)/* "all" */
{
   i1=1;
   i2=0;
   i3=1;
   i4=1;
   i5=1;
   goto start;
}
else if (flag==1) i1=1;/* "3" */
else if (flag==2) i2=1;/* "det" */
else if (flag==3) i3=1;/* "kon_bas" */
else if (flag==4) i4=1;/* "kov_met" */
else if (flag==5) i5=1;/* "kon_met" */
if (i2==1) {i1=1; }
if (i3==1) {i1=1; i2=0;}
if (i4==1) {i1=1; }
if (i5==1) {i1=1; i4=1;}
start:
/*---------------------------------------------------- sdc-conditioning */
zeta = e3 / condfac;
/*------------------------------------ interpolation of kovariant g1,g2 */
for (ialpha=0; ialpha<2; ialpha++)
{
   for (idim=0; idim<3; idim++)
   {
      gkov[idim][ialpha]=0.0;
      for (inode=0; inode<iel; inode++)
      {
         gkov[idim][ialpha] +=

            deriv[ialpha][inode] * (x[idim][inode]+zeta*a3[idim][inode]);
      }
   }
}
/*------------------------------------------------- interpolation of g3 */
if (i1==1)
{
   for (idim=0; idim<3; idim++)
   {
      gkov[idim][2]=0.0;
      for (inode=0; inode<iel; inode++)
      {
         gkov[idim][2] +=

            funct[inode] * a3[idim][inode];
      }
   }
}
/*--------------------------------------------------------- det of gkov */
if (i2==1)
{
   dserror("i2 not yet implemented ");
}
/*--------------- kontravariant basis vectors g1,g2,g3 (inverse of kov) */
if (i3==1)
{
   math_array_copy(gkov,3,3,gkon);
   math_inv3(gkon,det);
   math_tran(gkon,3);
}
/*--------------------------------------------- kovariant metrik tensor */
if (i4==1)
{
   for (i=0; i<3; i++)
   {
      for (j=i; j<3; j++)
      {
         gmkov[i][j]=0.0;
         for (k=0; k<3; k++)
         gmkov[i][j] += gkov[k][i]*gkov[k][j];
      }
   }
         gmkov[1][0] = gmkov[0][1];
         gmkov[2][0] = gmkov[0][2];
         gmkov[2][1] = gmkov[1][2];
}
/*----------------------------------------- kontravariant metrik tensor */
if (i5==1)
{
    math_array_copy(gmkov,3,3,gmkon);
    math_inv3(gmkon,&det_dummy);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_tmtr */
#endif
