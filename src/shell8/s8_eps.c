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
 | calc green lagrange strains from shell metrics         m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_eps(DOUBLE *strain,DOUBLE **gmkovc, DOUBLE **gmkovr)
{
#ifdef DEBUG
dstrc_enter("s8_eps");
#endif
/*----------------------------------------------------------------------*/
strain[0] = 0.5*(gmkovc[0][0] - gmkovr[0][0]);
strain[1] = 0.5*(gmkovc[0][1] - gmkovr[0][1]);
strain[2] = 0.5*(gmkovc[0][2] - gmkovr[0][2]);
strain[3] = 0.5*(gmkovc[1][1] - gmkovr[1][1]);
strain[4] = 0.5*(gmkovc[1][2] - gmkovr[1][2]);
strain[5] = 0.5*(gmkovc[2][2] - gmkovr[2][2]);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_eps */
#endif
