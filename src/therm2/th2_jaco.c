/*======================================================================*/
/*!
\file
\brief contains the routine 'w1_jaco' which calculates the operator
       matrix at point r,s for a wall element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/
#ifndef CCADISCRET
#ifdef D_THERM2

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*/
/*!
\brief Calculate operator matrix at point r,s

\param   **deriv     DOUBLE   (o)   natural derivatives of shape functions
                                    point r,s
\param   **xjm       DOUBLE   (o)   Jacobi matrix at r,s
                                    (isoparametric concept)
\param   *det        DOUBLE   (o)   Jacobi determinant at r,s
\param   *ele        ELEMENT  (i)   pointer to current element
\param   iel         INT      (i)   number of element nodes
\return void

\author bborn
\date 03/06
*/
void th2_jaco(DOUBLE **deriv,
              DOUBLE **xjm,
              DOUBLE *det,
              ELEMENT *ele,
              INT iel)
{
  /*--------------------------------------------------------------------*/
  INT k;  /* counter */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_jaco");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize Jacobian */
  xjm[0][0] = 0.0 ;
  xjm[0][1] = 0.0 ;
  xjm[1][0] = 0.0 ;
  xjm[1][1] = 0.0 ;
  /*--------------------------------------------------------------------*/
  /* Jacobian at point r,s */
  /* loop all element nodes */
  for (k=0; k<iel; k++)
  {
    xjm[0][0] += deriv[0][k] * ele->node[k]->x[0] ;
    xjm[0][1] += deriv[0][k] * ele->node[k]->x[1] ;
    xjm[1][0] += deriv[1][k] * ele->node[k]->x[0] ;
    xjm[1][1] += deriv[1][k] * ele->node[k]->x[1] ;
  }
  /* determinant of Jacobian */
  *det = xjm[0][0]*xjm[1][1] - xjm[1][0]*xjm[0][1];
  /*--------------------------------------------------------------------*/
  /* check negative Jacobian determinant */
  if (*det < 0.0)
  {
    dserror("NEGATIVE JACOBIAN DETERMINANT");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_jaco */


/*======================================================================*/
#endif /*end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
#endif
