/*!----------------------------------------------------------------------
\file
\brief ls2_intg.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/
static DOUBLE     q12 = 1.0/2.0;
static DOUBLE     q13 = 1.0/3.0;



/*!----------------------------------------------------------------------
\brief construct integration data for 2D level set element

<pre>                                                            irhan 05/04
construct integration data for 2D level set element
</pre>

*----------------------------------------------------------------------*/
void ls2_intg(
  LS2_INTG_DATA   *data
  )
{
  INT        i,j;
  DOUBLE     DUMMY = 0.0;

#ifdef DEBUG
  dstrc_enter("ls2_intg");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      data->wgtq[i][j] = 0.0;
      data->wgtt[i][j] = 0.0;
      data->xgq[i][j] = 0.0;
      data->xgtr[i][j] = 0.0;
      data->xgts[i][j] = 0.0;
    }
  }

  /* for quadrilateral element  */
  /* weights */
  data->wgtq[0][0] = 2.0000000000000;
  data->wgtq[0][1] = 1.0000000000000;
  data->wgtq[1][1] = 1.0000000000000;
  data->wgtq[0][2] = 0.5555555555555;
  data->wgtq[1][2] = 0.8888888888888;
  data->wgtq[2][2] = 0.5555555555555;
  /* Gauss points */
  data->xgq[0][0] =  0.0000000000000;
  data->xgq[0][1] = -0.5773502691896;
  data->xgq[1][1] =  0.5773502691896;
  data->xgq[0][2] = -0.7745966692414;
  data->xgq[1][2] =  0.0000000000000;
  data->xgq[2][2] =  0.7745966692414;

  /* for triangular element */
  /* weights */
  data->wgtt[0][0] = 0.5000000000000;
  data->wgtt[0][1] = DUMMY;
  data->wgtt[1][1] = DUMMY;
  data->wgtt[0][2] = 0.1666666666667;
  data->wgtt[1][2] = 0.1666666666667;
  data->wgtt[2][2] = 0.1666666666667;
  /* Gauss points */
  data->xgtr[0][0] = 0.3333333333333;
  data->xgtr[0][1] = DUMMY;
  data->xgtr[1][1] = DUMMY;
  data->xgtr[0][2] = 0.6666666666667;
  data->xgtr[1][2] = 0.1666666666667;
  data->xgtr[2][2] = 0.1666666666667;

  data->xgts[0][0] = 0.3333333333333;
  data->xgts[0][1] = DUMMY;
  data->xgts[1][1] = DUMMY;
  data->xgts[0][2] = 0.1666666666667;
  data->xgts[1][2] = 0.6666666666667;
  data->xgts[2][2] = 0.1666666666667;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of ls2_intg */
/*! @} (documentation module close)*/
#endif
