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
 |                                                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
/*
C.......................................................................
C!    GEAENDERTE METRIK DES VERF. SCHALENRAUMS INFOLGE ENHANCED STRAIN .
C.......................................................................
*/
void s8_vthv(DOUBLE **gmkovc,
             DOUBLE **gmkonc,
             DOUBLE  *epsh,
             DOUBLE  *detc,
             DOUBLE   e3,
             DOUBLE   condfac)
{
DOUBLE det_dummy;
DOUBLE zeta;

#ifdef DEBUG
dstrc_enter("s8_vthv");
#endif
/*----------------------------------------------------------------------*/
zeta = e3 / condfac;
/*----------------------------------------------------------------------*/
gmkovc[0][0] = gmkovc[0][0] + 2.0 * (epsh[0]+zeta*epsh[6]);
gmkovc[1][0] = gmkovc[1][0] +       (epsh[1]+zeta*epsh[7]);
gmkovc[2][0] = gmkovc[2][0] +       (epsh[2]+zeta*epsh[8]);
gmkovc[1][1] = gmkovc[1][1] + 2.0 * (epsh[3]+zeta*epsh[9]);
gmkovc[2][1] = gmkovc[2][1] +       (epsh[4]+zeta*epsh[10]);
gmkovc[2][2] = gmkovc[2][2] + 2.0 * (epsh[5]+zeta*epsh[11]);
gmkovc[0][2] = gmkovc[2][0];
gmkovc[1][2] = gmkovc[2][1];
gmkovc[0][1] = gmkovc[1][0];
/*----------------------------------------------------------------------*/
math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = -det_dummy;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_vthv */
#endif
