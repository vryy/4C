/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'saxiintg' which gives the coordinates and
weight factors for numerical integration of a axisymmetric shell element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*-----------------------------------------------------------------------*/
#ifdef D_AXISHELL
#include "../headers/standardtypes.h"
#include "axishell.h"

/*! 
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief coordinates and weight factors for numerical integration 

<pre>                                                              mn 05/03 
This routine  gives the coordinates and weight factors for numerical
integration of a axisymmetric shell element. Gauss-Integration with 5 GP 
from 0 to 1 is used here.

</pre>
\param *data   SAXI_DATA  (o)   structure containing the coordinates and weighting factors

\warning There is nothing special to this routine
\return void                                               
\sa calling:   ---; 
    called by: saxiinit, saxi_static_ke, saxi_cal_fext, saxi_cal_stress;

*----------------------------------------------------------------------*/
void saxiintg(SAXI_DATA   *data)
{
#ifdef DEBUG 
dstrc_enter("saxiintg");
#endif

/* gauss sampling points in xsi-direction */
data->xgr[0] = 0.046910077030668;  
data->xgr[1] = 0.2307653449471585;   
data->xgr[2] = 0.5;
data->xgr[3] = 1.0-data->xgr[1];
data->xgr[4] = 1.0-data->xgr[0];

/* weighting factors */  
data->wgt[0] = 0.236926885056189 / 2.0; 
data->wgt[1] = 0.478628670499366 / 2.0; 
data->wgt[2] = 0.568888888888889 / 2.0; 
data->wgt[3] = data->wgt[1];
data->wgt[4] = data->wgt[0];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of saxiintg */
/*----------------------------------------------------------------------*/
#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
