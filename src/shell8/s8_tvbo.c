#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | B-Operator for compatible strains                      m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_tvbo(double      e1,
                double      e2,
                double    **bop,
                double     *funct,
                double    **deriv,
                int         iel,
                int         numdf,
                double    **akov,
                double    **a3kvp,
                int         nsansq)
{
int inode,node_start;
double a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z;
double a31x,a31y,a31z,a32x,a32y,a32z;
double pk,pk1,pk2;
#ifdef DEBUG 
dstrc_enter("s8_tvbo");
#endif
/*----------------------------------------------------------------------*/
      a1x=akov[0][0];
      a1y=akov[1][0];
      a1z=akov[2][0];
      a2x=akov[0][1];
      a2y=akov[1][1];
      a2z=akov[2][1];
      a3x=akov[0][2];
      a3y=akov[1][2];
      a3z=akov[2][2];
      a31x=a3kvp[0][0];
      a31y=a3kvp[1][0];
      a31z=a3kvp[2][0];
      a32x=a3kvp[0][1];
      a32y=a3kvp[1][1];
      a32z=a3kvp[2][1];
/*----------------------------------------------------------------------*/
for (inode=0; inode<iel; inode++)
{
   pk  = funct[inode];
   pk1 = deriv[0][inode];
   pk2 = deriv[1][inode];
   
   node_start = inode*numdf;

if (nsansq != 6)   
{
   bop[0][node_start+0]= pk1*a1x;
   bop[0][node_start+1]= pk1*a1y;
   bop[0][node_start+2]= pk1*a1z;
   bop[0][node_start+3]= 0.0;
   bop[0][node_start+4]= 0.0;
   bop[0][node_start+5]= 0.0;

   bop[1][node_start+0]= pk2*a1x + pk1*a2x; 
   bop[1][node_start+1]= pk2*a1y + pk1*a2y;
   bop[1][node_start+2]= pk2*a1z + pk1*a2z;
   bop[1][node_start+3]= 0.0;
   bop[1][node_start+4]= 0.0;
   bop[1][node_start+5]= 0.0;
}
if (nsansq == 0)
{   
   bop[2][node_start+0]= pk1*a3x;
   bop[2][node_start+1]= pk1*a3y;
   bop[2][node_start+2]= pk1*a3z;
   bop[2][node_start+3]= pk *a1x;
   bop[2][node_start+4]= pk *a1y;
   bop[2][node_start+5]= pk *a1z;
}
if (nsansq != 6)
{   
   bop[3][node_start+0]= pk2*a2x;
   bop[3][node_start+1]= pk2*a2y;
   bop[3][node_start+2]= pk2*a2z;
   bop[3][node_start+3]= 0.0;
   bop[3][node_start+4]= 0.0;
   bop[3][node_start+5]= 0.0;
}   
if (nsansq == 0)
{   
   bop[4][node_start+0]= pk2*a3x;
   bop[4][node_start+1]= pk2*a3y;
   bop[4][node_start+2]= pk2*a3z;
   bop[4][node_start+3]= pk *a2x;
   bop[4][node_start+4]= pk *a2y;
   bop[4][node_start+5]= pk *a2z;
}   
   bop[5][node_start+0]= 0.0;
   bop[5][node_start+1]= 0.0;
   bop[5][node_start+2]= 0.0;
   bop[5][node_start+3]= pk *a3x;
   bop[5][node_start+4]= pk *a3y;
   bop[5][node_start+5]= pk *a3z;
   
   bop[6][node_start+0]= pk1*a31x;
   bop[6][node_start+1]= pk1*a31y;
   bop[6][node_start+2]= pk1*a31z;
   bop[6][node_start+3]= pk1*a1x;
   bop[6][node_start+4]= pk1*a1y;
   bop[6][node_start+5]= pk1*a1z;
   
   bop[7][node_start+0]= pk1*a32x + pk2*a31x;
   bop[7][node_start+1]= pk1*a32y + pk2*a31y;
   bop[7][node_start+2]= pk1*a32z + pk2*a31z;
   bop[7][node_start+3]= pk1* a2x + pk2* a1x;
   bop[7][node_start+4]= pk1* a2y + pk2* a1y;
   bop[7][node_start+5]= pk1* a2z + pk2* a1z;
   
   bop[8][node_start+0]= 0.0;
   bop[8][node_start+1]= 0.0;
   bop[8][node_start+2]= 0.0;
   bop[8][node_start+3]= pk *a31x + pk1*a3x;
   bop[8][node_start+4]= pk *a31y + pk1*a3y;
   bop[8][node_start+5]= pk *a31z + pk1*a3z;
   
   bop[9][node_start+0]= pk2*a32x;
   bop[9][node_start+1]= pk2*a32y;
   bop[9][node_start+2]= pk2*a32z;
   bop[9][node_start+3]= pk2*a2x;
   bop[9][node_start+4]= pk2*a2y;
   bop[9][node_start+5]= pk2*a2z;
   
   bop[10][node_start+0]= 0.0;
   bop[10][node_start+1]= 0.0;
   bop[10][node_start+2]= 0.0;
   bop[10][node_start+3]= pk *a32x + pk2*a3x;
   bop[10][node_start+4]= pk *a32y + pk2*a3y;
   bop[10][node_start+5]= pk *a32z + pk2*a3z;
   
   bop[11][node_start+0]= 0.0;
   bop[11][node_start+1]= 0.0;
   bop[11][node_start+2]= 0.0;
   bop[11][node_start+3]= 0.0;
   bop[11][node_start+4]= 0.0;
   bop[11][node_start+5]= 0.0;
   
   
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tvbo */
