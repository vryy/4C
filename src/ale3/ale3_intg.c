/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_intg' which gives the coordinates and
weight factors for numerical integration of a 3D ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief coordinates and weight factors for numerical integration 

<pre>                                                              mn 06/02 
This routine  gives the coordinates and weight factors for numerical
integration of a 3D ale element.

</pre>
\param *ele    ELEMENT    (i)   the element
\param *data   ALE3_DATA  (o)   structure containing the coordinates and weighting factors

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale3_static_ke

*----------------------------------------------------------------------*/
void ale3_intg(const ELEMENT   *ele,
              ALE3_DATA        *data)
{
int i, k;

double zero  = 0.0;
double one   = 1.0;
double two   = 2.0;
double three = 3.0;

double  q12, q13, q16, q23;
double  xgr[13][8],xgs[13][8],wgtt[13][8];
static double xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("ale3_intg");
#endif
/*----------------------------------------------------------------------*  
 |     INTEGRATION PARAMETERS FOR    H E X A H E D R A L     ELEMENTS   |
 |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
 |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
 *----------------------------------------------------------------------*/       
switch(ele->e.ale3->nGP[0])/* direction r */
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
   dserror("unknown number of gaussian points in ale3_intg");
break;
}
switch(ele->e.ale3->nGP[1])/* direction s */
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
   dserror("unknown number of gaussian points in ale3_intg");
break;
}
switch(ele->e.ale3->nGP[2])/* direction t */
{
case 2:
  data->xgpt[0] = -0.5773502691896;
  data->xgpt[1] =  0.5773502691896;
  
  data->wgtt[0] = 1.0            ;
  data->wgtt[1] = 1.0            ;
break;
case 1:
  data->xgpt[0] = 0.0;
  
  data->wgtt[0] = 2.0;
break;
default:
   dserror("unknown number of gaussian points in ale3_intg");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale3_intg */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
