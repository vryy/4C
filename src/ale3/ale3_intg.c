/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_intg' which gives the coordinates and
weight factors for numerical integration of a 3D ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

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
void ale3_intg(
    const ELEMENT   *ele,
    ALE3_DATA       *data
    )
{

  DOUBLE  q14, q16, q124;
  DOUBLE  palpha,pbeta;

#ifdef DEBUG
  dstrc_enter("ale3_intg");
#endif

  switch(ele->distyp)
  {
    case hex8:
    case hex20:
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
      break; /* end case hex820 */
    case tet4:
    case tet10:
      q14 = 1.0/4.0;
      q16 = 1.0/6.0;
      q124= 1.0/24.0;
      palpha = (5.0+3.0*sqrt(5.0))/20.0;
      pbeta  = (5.0-sqrt(5.0))/20.0;

      /*----------------------------------------------------------------------*
        |     INTEGRATION PARAMETERS FOR    T E T R A H E D R A L   ELEMENTS   |
        |     GAUSS SAMPLING POINTS  AT     R/S/T-COORDINATES   RESPECTIVELY   |
        |                            AND    CORRESPONDING WEIGHTING  FACTORS   |
       *----------------------------------------------------------------------*/

      /*----------------------------------------------------------------------*
        |    GAUSS INTEGRATION         1 SAMPLING POINT, DEG.OF PRECISION 1    |
        |                              CASE 0                                  |
       *----------------------------------------------------------------------*/
      switch(ele->e.ale3->nGP[0])/* direction r */
      {
        case 1:
          data->xgpr[0]    =  q14 ;
          data->xgps[0]    =  q14 ;
          data->xgpt[0]    =  q14 ;
          data->wgtr[0]    =  q16 ;
          break;
          /*----------------------------------------------------------------*
            | GAUSS INTEGRATION    4 SAMPLING POINTS, DEG.OF PRECISION 2    |
            |                      CASE 1                                   |
           *----------------------------------------------------------------*/
        case 4:
          data->xgpr[0]    =    pbeta ;
          data->xgpr[1]    =    palpha;
          data->xgpr[2]    =    pbeta ;
          data->xgpr[3]    =    pbeta ;
          data->xgps[0]    =    pbeta ;
          data->xgps[1]    =    pbeta ;
          data->xgps[2]    =    palpha;
          data->xgps[3]    =    pbeta ;
          data->xgpt[0]    =    pbeta ;
          data->xgpt[1]    =    pbeta ;
          data->xgpt[2]    =    pbeta ;
          data->xgpt[3]    =    palpha;
          data->wgtr[0]    =    q124  ;
          data->wgtr[1]    =    q124  ;
          data->wgtr[2]    =    q124  ;
          data->wgtr[3]    =    q124  ;
          break;
        default:
          dserror("unknown number of gausian points");
          break;
      }

      break; /* end if tet410 */
    default:
      dserror("unknown typ of discretisation");
      break;
  } /* end switch */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ale3_intg */


#endif
/*! @} (documentation module close)*/
