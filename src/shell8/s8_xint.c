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
 | SUPERPOSITION OF WEIGHTED FUNCTION VALUES              m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_xint(DOUBLE *result, DOUBLE *values, DOUBLE *funct, INT iel)
{
INT              i;
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
