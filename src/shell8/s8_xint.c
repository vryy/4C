#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | SUPERPOSITION OF WEIGHTED FUNCTION VALUES              m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_xint(double *result, double *values, double *funct, int iel)
{
int              i;
#ifdef DEBUG 
dstrc_enter("s8_xint");
#endif
/*----------------------------------------------------------------------*/
*result=0.0;
for (i=0; i<iel; i++) *result += (funct[i]*values[i]);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_xint */
#endif
