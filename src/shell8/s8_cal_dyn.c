#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 | calculation of gaussian points mass with update of displacements     |
 *----------------------------------------------------------------------*/
void s8_tmas(double *funct, double *thick, double **emass, int iel, int numdf,
             double facv, double facw, double facvw)
{
int           i,j,k;
double        he,hehe;
double        helpf;
double        help;
#ifdef DEBUG 
dstrc_enter("s8_tmas");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------- half element thickness at gaussian point */
he  = 0.0;
for (i=0; i<iel; i++) he += thick[i] * funct[i];
he   /= 2.0;
hehe  = he*he;
/*---------------------------------------------------- make mass matrix */
for (i=0; i<iel; i++)
for (j=0; j<iel; j++)
{
    helpf = funct[i] * funct[j];

    help = facv * helpf;
    for (k=0; k<3; k++)
    {
       emass[j*numdf+k][i*numdf+k] += help;
    }

    help = facw * helpf * hehe;
    for (k=3; k<6; k++)
    {
       emass[j*numdf+k][i*numdf+k] += help;
    }    
    
    if (ABS(facvw) >= EPS14)
    {
       help = facvw * helpf * he;
       emass[j*numdf+3][i*numdf+0] += help;
       emass[j*numdf+4][i*numdf+1] += help;
       emass[j*numdf+5][i*numdf+2] += help;
       emass[j*numdf+0][i*numdf+3] += help;
       emass[j*numdf+1][i*numdf+4] += help;
       emass[j*numdf+2][i*numdf+5] += help;
    }
    
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tmas */

#endif
