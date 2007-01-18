/*======================================================================*/
/*!
\file
\brief Deformation and displacement gradients

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>
*/
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"


/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Calculate material deformation gradient and displacement gradient
       at Gauss point

\param   enod         INT     (i)  number of element nodes
\param   deriv[][]    DOUBLE  (i)  natural derivatives of shape functions
\param   xji[][]      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   disgrdv[]    DOUBLE  (o)  displacement gradient vect at (r,s,t)
\param   defgrd[][]   DOUBLE  (o)  deformation gradient tensor at (r,s,t)

\return void

\author bborn
\date 12/06
*/
void so3_def_grad(INT enod,
                  DOUBLE edis[MAXNOD_SOLID3][NDIM_SOLID3],
                  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE disgrdv[NUMDFGR_SOLID3],
                  DOUBLE defgrd[NDIM_SOLID3][NDIM_SOLID3])
{
  /*--------------------------------------------------------------------*/
  INT inod;  /* current node */
  INT idfg;
  DOUBLE N_X, N_Y, N_Z;  /* derivative w.r. to reference coords */
/*   DOUBLE defgrdm[NUMDFGR_SOLID3][NUMSTR_SOLID3];  def. grad. matrix */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_def_grad");
#endif

  /*--------------------------------------------------------------------*/
  /* initialise vector to zero */
  for (idfg=0; idfg<NUMDFGR_SOLID3; idfg++)
  {
    disgrdv[idfg] = 0.0;
  }

  /*--------------------------------------------------------------------*/
  /* construct displacement gradient in vector */
  /* loop over element nodes */
  for (inod=0; inod<enod; inod++)
  {
    /* derivatives of shape function with respect to material XYZ-frame */
    N_X = xji[0][0]*deriv[inod][0]
        + xji[0][1]*deriv[inod][1]
        + xji[0][2]*deriv[inod][2];
    N_Y = xji[1][0]*deriv[inod][0]
        + xji[1][1]*deriv[inod][1]
        + xji[1][2]*deriv[inod][2];
    N_Z = xji[2][0]*deriv[inod][0]
        + xji[2][1]*deriv[inod][1]
        + xji[2][2]*deriv[inod][2];

    /* add contribution of shape function at current node */
    disgrdv[0] += N_X * edis[inod][0];  /* u_{1,1} = d u_1/d X_1 */
    disgrdv[1] += N_Y * edis[inod][1];  /* u_{2,2} = d u_2/d X_2 */
    disgrdv[2] += N_Z * edis[inod][2];  /* u_{3,3} = d u_3/d X_3 */
    disgrdv[3] += N_Y * edis[inod][0];  /* u_{1,2} = d u_1/d X_2 */
    disgrdv[4] += N_X * edis[inod][1];  /* u_{2,1} = d u_2/d X_1 */
    disgrdv[5] += N_Z * edis[inod][1];  /* u_{2,3} = d u_2/d X_3 */
    disgrdv[6] += N_Y * edis[inod][2];  /* u_{3,2} = d u_3/d X_2 */
    disgrdv[7] += N_X * edis[inod][2];  /* u_{3,1} = d u_3/d X_1 */
    disgrdv[8] += N_Z * edis[inod][0];  /* u_{1,3} = d u_1/d X_3 */
  }

  /*--------------------------------------------------------------------*/
  /* set deformation gradient matrix */
  /* [ 1+u_{1,1}                       |   u_{1,2}               u_{1,3} ]
   * [            1+u_{2,2}            |   u_{2,1}    u_{2,3}            ]
   * [                       1+u_{3,3} |              u_{3,2}    u_{3,1} ]
   * [ ---------  ---------  ---------   ---------  ---------  --------- ]
   * [              u_{1,2}            | 1+u_{1,1}    u_{1,3}            ]
   * [   u_{2,1}                       | 1+u_{2,2}               u_{2,3} ]
   * [ ---------  ---------  ---------   ---------  ---------  --------- ]
   * [                         u_{2,3} |            1+u_{2,2}    u_{2,1} ]
   * [              u_{3,2}            |   u_{3,1}  1+u_{3,3}            ]
   * [ ---------  ---------  ---------   ---------  ---------  --------- ]
   * [   u_{3,1}                       |   u_{3,2}             1+u_{3,3} ]
   * [                         u_{1,3} |              u_{1,2}  1+u_{1,1} ]
   */
/*   defgrdm[0][0] = 1.0 + disgrd[0]; */
/*   defgrdm[0][3] =       disgrd[3]; */
/*   defgrdm[0][5] =       disgrd[8]; */
/*   defgrdm[1][1] = 1.0 + disgrd[1]; */
/*   defgrdm[1][3] =       disgrd[4]; */
/*   defgrdm[1][4] =       disgrd[5]; */
/*   defgrdm[2][2] = 1.0 + disgrd[2]; */
/*   defgrdm[2][4] =       disgrd[6]; */
/*   defgrdm[2][5] =       disgrd[7]; */
/*   defgrdm[3][1] =       disgrd[3]; */
/*   defgrdm[3][3] = 1.0 + disgrd[0]; */
/*   defgrdm[3][4] =       disgrd[8]; */
/*   defgrdm[4][0] =       disgrd[4]; */
/*   defgrdm[4][3] = 1.0 + disgrd[1]; */
/*   defgrdm[4][5] =       disgrd[5]; */
/*   defgrdm[5][2] =       disgrd[5]; */
/*   defgrdm[5][4] = 1.0 + disgrd[1]; */
/*   defgrdm[5][5] =       disgrd[4]; */
/*   defgrdm[6][1] =       disgrd[6]; */
/*   defgrdm[6][3] =       disgrd[7]; */
/*   defgrdm[6][4] = 1.0 + disgrd[2]; */
/*   defgrdm[7][0] =       disgrd[7]; */
/*   defgrdm[7][3] =       disgrd[6]; */
/*   defgrdm[7][5] = 1.0 + disgrd[2]; */
/*   defgrdm[8][2] =       disgrd[8]; */
/*   defgrdm[8][4] =       disgrd[3]; */
/*   defgrdm[8][5] = 1.0 + disgrd[0]; */

  /*--------------------------------------------------------------------*/
  /* (material) deformation gradient tensor */
  /*        d x   d (X+u)   [ 1+u_{1,1}     u_{1,2}     u_{1,3} ]
   *    F = --- = ------- = [   u_{2,1}   1+u_{2,2}     u_{2,3} ]
   *        d X     d X     [   u_{3,1}     u_{3,2}   1+u_{3,3} ]
   */
  defgrd[0][0] = 1.0 + disgrdv[0];
  defgrd[0][1] =       disgrdv[3];
  defgrd[0][2] =       disgrdv[8];
  defgrd[1][0] =       disgrdv[4];
  defgrd[1][1] = 1.0 + disgrdv[1];
  defgrd[1][2] =       disgrdv[5];
  defgrd[2][0] =       disgrdv[7];
  defgrd[2][1] =       disgrdv[6];
  defgrd[2][2] = 1.0 + disgrdv[2];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_grad */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
