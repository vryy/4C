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
 | integrate material law and stresses in thickness direction of shell  |
 |                                                      m.gee 06/01     |
 *----------------------------------------------------------------------*/
void s8_tvma(DOUBLE **D, DOUBLE **C, DOUBLE *stress, DOUBLE *stress_r,
                DOUBLE e3, DOUBLE fact, DOUBLE condfac)
{
INT i,i6,j,j6;
DOUBLE zeta;
DOUBLE stress_fact, C_fact;
#ifdef DEBUG
dstrc_enter("s8_tvma");
#endif
/*----------------------------------------------------------------------*/
zeta = e3/condfac;
for (i=0; i<6; i++)
{
   i6=i+6;
   stress_fact   = stress[i]*fact;
   stress_r[i]  += stress_fact;
   stress_r[i6] +=  stress_fact * zeta;
   for (j=0; j<6; j++)
   {
      j6 = j+6;
      C_fact     = C[i][j]*fact;
      D[i][j]   += C_fact;
      D[i6][j]  += C_fact*zeta;
      D[i6][j6] += C_fact*zeta*zeta;
   }
}
/*-------------------------------------------------------- symmetrize D */
for (i=0; i<12; i++)
{
   for (j=i+1; j<12; j++)
   {
      D[i][j]=D[j][i];
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_tvma */
#endif

