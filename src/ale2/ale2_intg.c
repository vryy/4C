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

DOUBLE  q12 = 1.0/2.0;
DOUBLE  q13 = 1.0/3.0;
DOUBLE  q16 = 1.0/6.0;
DOUBLE  q23 = 2.0/3.0;
DOUBLE  xgr[13][8],xgs[13][8],wgtt[13][8];
static DOUBLE xg[6][6],wgt[6][6];
#ifdef DEBUG 
dstrc_enter("ale2_intg");
#endif

switch (ele->distyp)
{
    case quad4:
    case quad8:
    case quad9:
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
        break; /* end case quad489*/
    case tri3:
    case tri6:
        /*----------------------------------------------------------------------*  
          |     INTEGRATION PARAMETERS FOR    T R I                   ELEMENTS   |
          |     GAUSS SAMPLING POINTS  AT     R/S-COORDINATES     RESPECTIVELY   |
          |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
         *----------------------------------------------------------------------*/       
        switch(ele->e.ale2->nGP[0])
        {
            case 1:
                data->xgpr[0] = q13;
                data->xgps[0] = q13;
                data->wgtr[0] = q12;
                break;
            case 3:
                data->xgpr[0] = q12;
                data->xgpr[1] = q12;
                data->xgpr[2] = ZERO;
                data->xgps[0] = ZERO;
                data->xgps[1] = q12;
                data->xgps[2] = q12;
                data->wgtr[0] = q16;
                data->wgtr[1] = q16;
                data->wgtr[2] = q16;
                break;
            default:
                dserror("unknown number of gaussian points");
                break;
        }
        break; /* end case tri36 */
    default:
        dserror("unknown typ of discretisation");
        break; /* end default */
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
