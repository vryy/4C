/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_intg' which gives the coordinates and
weight factors for numerical integration of a 2D ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief coordinates and weight factors for numerical integration 

<pre>                                                              mn 06/02 
This routine  gives the coordinates and weight factors for numerical
integration of a 2D ale element.

</pre>
\param *ele    ELEMENT    (i)   the element
\param *data   ALE2_DATA  (o)   structure containing the coordinates and weighting factors

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_intg(const ELEMENT   *ele,
              ALE2_DATA        *data)
{
INT i, k;

DOUBLE zero  = 0.0;
DOUBLE one   = 1.0;
DOUBLE two   = 2.0;
DOUBLE three = 3.0;

DOUBLE  q12, q13, q16, q23;
DOUBLE  xgr[13][8],xgs[13][8],wgtt[13][8];
static DOUBLE xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("ale2_intg");
#endif
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    Q U A D                 ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
switch(ele->e.ale2->nGP[0])/* direction r */
{
case 2:
  data->xgpr[0] = -0.5773502691896;
  data->xgpr[1] =  0.5773502691896;
  
  data->wgtr[0] = 1.0            ;
  data->wgtr[1] = 1.0            ;
break;
case 1:
  data->xgpr[0] = 0.0;
  
  data->wgtr[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale2_intg");
break;
}
switch(ele->e.ale2->nGP[1])/* direction s */
{
case 2:
  data->xgps[0] = -0.5773502691896;
  data->xgps[1] =  0.5773502691896;
  
  data->wgts[0] = 1.0            ;
  data->wgts[1] = 1.0            ;
break;
case 1:
  data->xgps[0] = 0.0;
  
  data->wgts[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale2_intg");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2_intg */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
