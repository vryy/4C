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
 |                                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*
C.......................................................................
C!    GEAENDERTE METRIK DES VERF. SCHALENRAUMS INFOLGE ENHANCED STRAIN .
C.......................................................................
*/
void s8_tvkg(DOUBLE **estif,
             DOUBLE  *stress_r,
             DOUBLE  *funct,
             DOUBLE **deriv,
             INT      numdf,
             INT      iel,
             DOUBLE   weight,
             DOUBLE   e1,
             DOUBLE   e2)
{
INT     inode,jnode;
INT     i_indiz,j_indiz;
DOUBLE  pi;
DOUBLE  pj;
DOUBLE  d11;
DOUBLE  d12;
DOUBLE  d21;
DOUBLE  d22;

DOUBLE  pd1ij;
DOUBLE  pd1ji;
DOUBLE  pd2ij;
DOUBLE  pd2ji;

DOUBLE  xn;
DOUBLE  xm;
DOUBLE  yu;
DOUBLE  yo;
DOUBLE  yy;
DOUBLE  z;

DOUBLE  sn11;
DOUBLE  sn21;
DOUBLE  sn31;
DOUBLE  sn22;
DOUBLE  sn32;
DOUBLE  sn33;
DOUBLE  sm11;
DOUBLE  sm21;
DOUBLE  sm31;
DOUBLE  sm22;
DOUBLE  sm32;

#ifdef DEBUG
dstrc_enter("s8_tvkg");
#endif
/*----------------------------------------------------------------------*/
sn11 = stress_r[0];
sn21 = stress_r[1];
sn31 = stress_r[2];
sn22 = stress_r[3];
sn32 = stress_r[4];
sn33 = stress_r[5];
sm11 = stress_r[6];
sm21 = stress_r[7];
sm31 = stress_r[8];
sm22 = stress_r[9];
sm32 = stress_r[10];
/*----------------------------------------------------------------------*/
for (inode=0; inode<iel; inode++)
{
   for (jnode=0; jnode<=inode; jnode++)
   {
      pi = funct[inode];
      pj = funct[jnode];

      d11 = deriv[0][inode] * deriv[0][jnode];
      d12 = deriv[0][inode] * deriv[1][jnode];
      d21 = deriv[1][inode] * deriv[0][jnode];
      d22 = deriv[1][inode] * deriv[1][jnode];

      pd1ij = deriv[0][inode] * pj;
      pd1ji = deriv[0][jnode] * pi;
      pd2ij = deriv[1][inode] * pj;
      pd2ji = deriv[1][jnode] * pi;

      xn = (sn11*d11 + sn21*(d12+d21) + sn22*d22) * weight;

      xm = (sm11*d11 + sm21*(d12+d21) + sm22*d22) * weight;

      yu = (sn31*pd1ji + sn32*pd2ji) * weight;

      yo = (sn31*pd1ij + sn32*pd2ij) * weight;

      yy = (sm31*(pd1ij+pd1ji) + sm32*(pd2ij+pd2ji)) * weight;

      z  = pi*pj*sn33*weight;

      i_indiz = inode*numdf;
      j_indiz = jnode*numdf;

      estif[inode*numdf+0][jnode*numdf+0] += xn;
      estif[inode*numdf+1][jnode*numdf+1] += xn;
      estif[inode*numdf+2][jnode*numdf+2] += xn;

      estif[inode*numdf+3][jnode*numdf+0] += (xm+yu);
      estif[inode*numdf+4][jnode*numdf+1] += (xm+yu);
      estif[inode*numdf+5][jnode*numdf+2] += (xm+yu);

      estif[inode*numdf+0][jnode*numdf+3] += (xm+yo);
      estif[inode*numdf+1][jnode*numdf+4] += (xm+yo);
      estif[inode*numdf+2][jnode*numdf+5] += (xm+yo);

      estif[inode*numdf+3][jnode*numdf+3] += (yy+z);
      estif[inode*numdf+4][jnode*numdf+4] += (yy+z);
      estif[inode*numdf+5][jnode*numdf+5] += (yy+z);

      if (inode!=jnode)
      {
         estif[jnode*numdf+0][inode*numdf+0] += xn;
         estif[jnode*numdf+1][inode*numdf+1] += xn;
         estif[jnode*numdf+2][inode*numdf+2] += xn;

         estif[jnode*numdf+0][inode*numdf+3] += (xm+yu);
         estif[jnode*numdf+1][inode*numdf+4] += (xm+yu);
         estif[jnode*numdf+2][inode*numdf+5] += (xm+yu);

         estif[jnode*numdf+3][inode*numdf+0] += (xm+yo);
         estif[jnode*numdf+4][inode*numdf+1] += (xm+yo);
         estif[jnode*numdf+5][inode*numdf+2] += (xm+yo);

         estif[jnode*numdf+3][inode*numdf+3] += (yy+z);
         estif[jnode*numdf+4][inode*numdf+4] += (yy+z);
         estif[jnode*numdf+5][inode*numdf+5] += (yy+z);
      }

   } /* end loop over jnode */
} /* end loop over inode */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_tvkg */
#endif
