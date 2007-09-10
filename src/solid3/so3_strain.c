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
#ifndef CCADISCRET
#ifdef D_SOLID3


/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"


/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/

/*======================================================================*/
/*!
\brief Calculate linear strains at Gauss point

\param   enod       INT     (i)  number of element nodes
\param   **deriv    DOUBLE  (i)  natural derivatives of shape functions
\param   **xji      DOUBLE  (i)  Inverse Jacobi matrix at (r,s,t)
\param   *grdisv    DOUBLE  (o)  displacement gradient at (r,s,t)

\return void

\author bborn
\date 12/06
*/
void so3_strain_lin(ELEMENT *ele,
                    DOUBLE disgrdv[NUMDFGR_SOLID3],
                    DOUBLE strain[NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
  DOUBLE u11, u12, u13, u21, u22, u23, u31, u32, u33;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_strain_lin");
#endif

  /*--------------------------------------------------------------------*/
  /* temporary displacement gradient components */
  u11 = disgrdv[0];
  u22 = disgrdv[1];
  u33 = disgrdv[2];
  u12 = disgrdv[3];
  u21 = disgrdv[4];
  u23 = disgrdv[5];
  u32 = disgrdv[6];
  u31 = disgrdv[7];
  u13 = disgrdv[8];

  /*--------------------------------------------------------------------*/
  /* Linear strain measure */
  strain[0] = u11;  /* eps_{11} = u_{1,1} */
  strain[1] = u22;  /* eps_{22} = u_{2,2} */
  strain[2] = u33;  /* eps_{33} = u_{3,3} */
  strain[3] = u12 + u21;  /* 2*eps_{12} = u_{1,2}+u_{2,1} */
  strain[4] = u23 + u32;  /* 2*eps_{23} = u_{2,3}+u_{3,2} */
  strain[5] = u31 + u13;  /* 2*eps_{31} = u_{3,1}+u_{1,3} */


  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_strain_lin */


/*======================================================================*/
/*!
\brief Calculate Green-Lagrange strains at Gauss point

\param   ele        ELEMENT* (i)  current element
\param   disgrdv    DOUBLE[] (i)  displacement gradient vector
\param   strain     DOUBLE[] (o)  Green-Lagrange strain vector
\return void

\author bborn
\date 12/06
*/
void so3_strain_gl(ELEMENT *ele,
                   DOUBLE disgrdv[NUMDFGR_SOLID3],
                   DOUBLE strain[NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
  DOUBLE u11, u12, u13, u21, u22, u23, u31, u32, u33;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_strain_gl");
#endif

  /*--------------------------------------------------------------------*/
  /* temporary displacement gradient components */
  u11 = disgrdv[0];
  u22 = disgrdv[1];
  u33 = disgrdv[2];
  u12 = disgrdv[3];
  u21 = disgrdv[4];
  u23 = disgrdv[5];
  u32 = disgrdv[6];
  u31 = disgrdv[7];
  u13 = disgrdv[8];

  /*--------------------------------------------------------------------*/
  /* Green-Lagrange strain measure */
  /* E_{11} = u_{1,1} + 1/2*(u_{1,1}^2 + u_{2,1}^2 + u_{3,1}^2) */
  strain[0] = u11 + 0.5*(u11*u11 + u21*u21 + u31*u31);
  /* E_{22} = u_{2,2} + 1/2*(u_{1,2}^2 + u_{2,2}^2 + u_{3,2}^2) */
  strain[1] = u22 + 0.5*(u12*u12 + u22*u22 + u32*u32);
  /* E_{33} = u_{3,3} + 1/2*(u_{1,3}^2 + u_{2,3}^2 + u_{3,3}^2) */
  strain[2] = u33 + 0.5*(u13*u13 + u23*u23 + u33*u33);
  /* 2*E_{12} = u_{1,2}+u_{2,1}
   *          + u_{1,1}*u_{1,2} + u_{2,1}*u_{2,2} + u_{3,1}*u_{3,2} */
  strain[3] = u12 + u21 + u11*u12 + u21*u22 + u31*u32;
  /* 2*E_{23} = u_{2,3}+u_{3,2}
   *          + u_{1,2}*u_{1,3} + u_{2,2}*u_{2,3} + u_{3,2}*u_{3,3} */
  strain[4] = u23 + u32 + u32*u33 + u22*u23 + u12*u13;
  /* 2*E_{31} = u_{3,1}+u_{1,3}
   *          + u_{1,3}*u_{1,1} + u_{2,3}*u_{2,1} + u_{3,3}*u_{3,1} */
  strain[5] = u31 + u13 + u31*u33 + u21*u23 + u11*u13;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_grad */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
#endif
