/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates all metrics                                m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8mtr(DOUBLE   **x,
              DOUBLE   **a3,
              DOUBLE     e3,
              DOUBLE   **gkov,
              DOUBLE   **gkon,
              DOUBLE   **gmkov,
              DOUBLE   **gmkon,
              DOUBLE    *det,
              DOUBLE    *funct,
              DOUBLE   **deriv,
              DOUBLE    *hte,
              INT        iel,
              DOUBLE     condfac,
              char       string[])
{
INT    i,j,k,idim,ialpha,inode;
INT    i1=0,i2=0,i3=0,i4=0,i5=0;
DOUBLE zeta,det_dummy;
DOUBLE h2;
DOUBLE a[3],b[3],c[3];
#ifdef DEBUG
dstrc_enter("s8_tmtr");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------- make decision about what to do */
if (strstr(string,"all")!=NULL)/* "all" */
{
   i1=1;
   i2=0;
   i3=1;
   i4=1;
   i5=1;
   goto start;
}
if (strstr(string,"3")!=NULL) i1=1;
if (strstr(string,"det")!=NULL) i2=1;
if (strstr(string,"kon_bas")!=NULL) i3=1;
if (strstr(string,"kov_met")!=NULL) i4=1;
if (strstr(string,"kon_met")!=NULL) i5=1;
if (i2) {i1=1; }
if (i3) {i1=1; i2=0;}
if (i4) {i1=1; }
if (i5) {i1=1; i4=1;}
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
         h2 = hte[inode]/2.0;

         gkov[idim][ialpha] +=

            deriv[ialpha][inode] * (x[idim][inode]+h2*zeta*a3[idim][inode]);
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
         h2 = hte[inode]/2.0;

         gkov[idim][2] +=

            funct[inode] * h2*a3[idim][inode];
      }
   }
}
/*--------------------------------------------------------- det of gkov */
if (i2==1)
{
   for (i=0; i<3; i++)
   {
      a[i] = gkov[i][0];
      b[i] = gkov[i][1];
      c[i] = gkov[i][2];
   }
   math_sppr(det,a,b,c);
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
