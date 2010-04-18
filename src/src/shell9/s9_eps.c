/*!----------------------------------------------------------------------
\file
\brief contains the routine
 - s9_eps: which calculates the green lagrange strains from the kovariant
           shell9 metrics


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calc green lagrange strains from shell metrics

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the green lagrange strains from the kovariant
shell9 metrics in current and reference configuration
</pre>
\param  DOUBLE *strain   (o) green lagrange strains
\param  DOUBLE **gmkovc  (i) kovariant metric in current configuration
\param  DOUBLE **gmkovr  (i) kovariant metric in referent configuration

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_call_mat()   [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_eps(DOUBLE strain[6],DOUBLE **gmkovc, DOUBLE **gmkovr)
{
#ifdef DEBUG
dstrc_enter("s9_eps");
#endif
/*----------------------------------------------------------------------*/
strain[0] = 0.5*(gmkovc[0][0] - gmkovr[0][0]);
strain[1] = 0.5*(gmkovc[0][1] - gmkovr[0][1]);
strain[2] = 0.5*(gmkovc[1][1] - gmkovr[1][1]);
strain[3] = 0.5*(gmkovc[0][2] - gmkovr[0][2]);
strain[4] = 0.5*(gmkovc[1][2] - gmkovr[1][2]);
strain[5] = 0.5*(gmkovc[2][2] - gmkovr[2][2]);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_eps */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
#endif
