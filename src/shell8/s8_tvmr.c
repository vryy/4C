#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates all metrics                                m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_tvmr(double   **x,
                double   **a3,
                double   **akov,
                double   **akon,
                double   **amkov,
                double   **amkon,
                double    *det,
                double    *funct,
                double   **deriv,
                int        iel,
                double   **a3kvp,
                int        flag)
{
int    i,j,k,idim,ialpha,inode;
int    i1=0,i2=0,i3=0,i4=0,i5=0,i6=0;
double det_dummy;
#ifdef DEBUG 
dstrc_enter("s8_tvmr");
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
   i6=1;
   goto start;
}
else if (flag==1) i1=1;/* "3" */
else if (flag==2) i2=1;/* "det" */
else if (flag==3) i3=1;/* "kon_bas" */
else if (flag==4) i4=1;/* "kov_met" */
else if (flag==5) i5=1;/* "kon_met" */
else if (flag==6) i6=1;/* "a3p" */
if (i2==1) {i1=1; }
if (i3==1) {i1=1; i2=0;}
if (i4==1) {i1=1; }
if (i5==1) {i1=1; i4=1;}
start:
/*------------------------------------ interpolation of kovariant a1,a2 */
for (ialpha=0; ialpha<2; ialpha++)
{
   for (idim=0; idim<3; idim++)
   {
      akov[idim][ialpha]=0.0;
      for (inode=0; inode<iel; inode++)
      {
         akov[idim][ialpha] += 
         
            deriv[ialpha][inode] * x[idim][inode];
      }
   }
}
/*------------------------------------------------- interpolation of a3 */
if (i1==1)
{
   for (idim=0; idim<3; idim++)
   {
      akov[idim][2]=0.0;
      for (inode=0; inode<iel; inode++)
      {
         akov[idim][2] +=
         
            funct[inode] * a3[idim][inode];
      }
   }
}
/*--------------------------------------------------------- det of akov */
if (i2==1)
{
   dserror("i2 not yet implemented ");
}
/*--------------- kontravariant basis vectors g1,g2,g3 (inverse of kov) */
if (i3==1)
{
   math_array_copy(akov,3,3,akon);
   math_inv3(akon,det);
   math_tran(akon,3);
}
/*--------------------------------------------- kovariant metrik tensor */
if (i4==1)
{
   for (i=0; i<3; i++)
   {
      for (j=i; j<3; j++)
      {
         amkov[i][j]=0.0;
         for (k=0; k<3; k++) 
         amkov[i][j] += akov[k][i]*akov[k][j];
      }
   }   
         amkov[1][0] = amkov[0][1];
         amkov[2][0] = amkov[0][2];
         amkov[2][1] = amkov[1][2];
}
/*----------------------------------------- kontravariant metrik tensor */
if (i5==1)
{
    math_array_copy(amkov,3,3,amkon);
    math_inv3(amkon,&det_dummy);
}
/*------------------------------------------- partial derivatives of a3 */
if (i6==1)
{
   for (ialpha=0; ialpha<2; ialpha++)
   {
      for (idim=0; idim<3; idim++)
      {
         a3kvp[idim][ialpha]=0.0;
         for (inode=0; inode<iel; inode++)
         {
            a3kvp[idim][ialpha] += deriv[ialpha][inode]*a3[idim][inode];
         }
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tvmr */
#endif
