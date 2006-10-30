/*======================================================================*/
/*!
\file
\brief Metrics of THERM3 element, i.e. metric conversion between physical
       and parameter space. 

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 09/06
*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/

/*======================================================================*/
/*!
\brief Jacobian matrix at coordinate (r,s,t) (e.g. Gauss point) by
       applying isoparametric concept

\param ele      ELEMENT*    (i)   the element
\param enod     INT         (i)   number of element nodes
\param deriv    DOUBLE*     (i)   derivatives of the shape functions
\param flag     INT         (i)   flag for inverse of Jacobian
                                     1 : inverse Jacobian
                                     else : don't inverse Jacobian
\param xjm      DOUBLE**    (o)   the Jacobian matrix
\param det      DOUBLE*     (o)   determinant of the Jacobian matrix
\param xji      DOUBLE**    (o)   inverse of Jacobian
\return void

\author bborn
\date 09/06
*/
void th3_metr_jaco(ELEMENT *ele,
                   INT      enod,
                   DOUBLE   deriv[MAXNOD_THERM3][NDIM_THERM3],
                   INT      flag,
                   DOUBLE   xjm[NDIM_THERM3][NDIM_THERM3],
                   DOUBLE  *det,
                   DOUBLE   xji[NDIM_THERM3][NDIM_THERM3])
{
  INT i, j, k;
  DOUBLE nodxyz;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_metr_jaco");
#endif
  
  /*--------------------------------------------------------------------*/
  /* initialise Jacobian to zero */
  for (i=0; i<NDIM_THERM3; i++)
  {
    for (j=0; j<NDIM_THERM3; j++)
    {
      xjm[i][j] = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
  /* loop all element nodes / shape functions and dimensions */
  /* The Jacobian matrix looks like:
   *         [ x_,r  y_,r  z_,r ]
   *     J = [ x_,s  y_,s  z_,s ]
   *         [ x_,t  y_,t  z_,t ]
   * This is a usual form used in the finite element environment,
   * however, it is the transpose of the common Jacobian
   *              d x_i
   *     J_ij^T = ------
   *              d r_j
   * This form is for instance applied in the material
   * deformation tensors  F  in geometrically non-linear mechanics.
   * In finite elements  J  and not  J^T  is used, because we want to
   * multiply on the _left_ with its inverse  J^{-1}
   *     [ N_,x ]^k           [ N_,r ]^k
   *     [ N_,y ]   =  J^{-1} [ N_,s ]
   *     [ N_,z ]             [ N_,t ]
   * with node index  k=1,...,elenod, shape functions  N
   */
  for (k=0; k<enod; k++)
  {
    for (i=0; i<NDIM_THERM3; i++)
    {
      nodxyz = ele->node[k]->x[i];
      for (j=0; j<NDIM_THERM3; j++)
      {
        xjm[j][i] += deriv[k][j] * nodxyz;
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* determinat */
  *det = xjm[0][0] * xjm[1][1] * xjm[2][2]  /* Sarrus' rule */
       + xjm[0][1] * xjm[1][2] * xjm[2][0]  /* [Bronstein, */
       + xjm[0][2] * xjm[1][0] * xjm[2][1]  /*  Taschenbuch der */
       - xjm[0][2] * xjm[1][1] * xjm[2][0]  /*  der Mathematik] */
       - xjm[0][0] * xjm[1][2] * xjm[2][1]
       - xjm[0][1] * xjm[1][0] * xjm[2][2];
  /* check if negative Jacobian */
  if (*det < 0.0)
  {
     dserror("NEGATIVE JACOBIAN DETERMINANT");
  }

  /*--------------------------------------------------------------------*/
  /* inverse of Jacobian */
  if (flag == 1)
  {
    xji[0][0] = xjm[1][1]*xjm[2][2] - xjm[1][2]*xjm[2][1];
    xji[0][1] = xjm[0][2]*xjm[2][1] - xjm[0][1]*xjm[2][2];
    xji[0][2] = xjm[0][1]*xjm[1][2] - xjm[0][2]*xjm[1][1];
    xji[1][0] = xjm[1][2]*xjm[2][0] - xjm[1][0]*xjm[2][2];
    xji[1][1] = xjm[0][0]*xjm[2][2] - xjm[0][2]*xjm[2][0];
    xji[1][2] = xjm[0][2]*xjm[1][0] - xjm[0][0]*xjm[1][2];
    xji[2][0] = xjm[1][0]*xjm[2][1] - xjm[1][1]*xjm[2][0];
    xji[2][1] = xjm[0][1]*xjm[2][0] - xjm[0][0]*xjm[2][1];
    xji[2][2] = xjm[0][0]*xjm[1][1] - xjm[0][1]*xjm[1][0];
    for (i=0; i<NDIM_THERM3; i++)
    {
      for (j=0; j<NDIM_THERM3; j++)
      {
        xji[i][j] = xji[i][j] / (*det);
      }  /* end of for */
    }  /* end of for */
  }  /* end of if */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_metr_jaco */


/*======================================================================*/
/*!
\brief Determine metric at Gauss point by reducing the Jacobian
       matrix 

\param ele     ELEMENT*     (i)  pointer to current element
\param nelelod INT          (i)  number of element nodes
\param deriv   DOUBLE**     (i)  derivatives of shape fct at Gauss point
\param sidredm DOUBLE**     (i)  matrix for dimension reduction
\param metr    DOUBLE*      (o)  metric
\return void

\author bborn
\date 09/06
*/
void th3_metr_surf(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3], 
                   DOUBLE sidredm[DIMSID_THERM3][NDIM_THERM3],
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian (dummy) */
  DOUBLE gamt[DIMSID_THERM3][NDIM_THERM3];  /* differential */
  DOUBLE metm[DIMSID_THERM3][DIMSID_THERM3];  /* metric matrix */
  INT idim, jdim, kdimsid, jdimsid;  /* dimension indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_metr_surf");
#endif

  /*--------------------------------------------------------------------*/
  /* get Jacobian matrix (xji is flagged off) */
  th3_metr_jaco(ele, nelenod, deriv, 0, xjm, &det, xji);

  /*--------------------------------------------------------------------*/
  /* build gamma tensor */
  /* Assume this integration is applied for a surface of the volumetric
   * element for which  t=0  . The surface is a plane in parameter space
   * with the Cartesian coordinates (r,s). The surface in 3dim
   * physical space is generally not a plane but a hypersuface.
   * The physical surface is described by
   *     { [-1,1]x[-1,1]  --->  R^3
   *     {     (r,s)      +-->  (x,y,z)
   * which assumes a surface of hexadron element.
   * If the physical coordinates  (x,y,z)  are varied with respect
   * to the parameters  (r,s)  the  gamma  tensor (a Jacobian) is
   * obtained
   *             [ x_,r  x_,s ]|
   *     gamma = [ y_,r  y_,s ]|
   *             [ z_,r  z_,s ]|(r,s)
   * The  gamma  tensor transforms locally (at (r,s)) parameter 
   * to physical coordinates.
   * Here the transposed of  gamma  is determined:
   *                               [ x_,r  y_,r  z_,r ]|
   *    gamma^T(r,s) = [ 1  0  0 ] [ x_,s  y_,s  z_,s ]|
   *                   [ 0  1  0 ] [ x_,t  y_,t  z_,t ]|(r,s)
   *                  =   sidredm            J
   * with  sidredm  the side reduction matrix and  J  the Jacobian
   * matrix (xjm) defined in FE-fashion (see th3_metr_jaco).
   * sidredm  depends on the surface parameter space (r,s). 
   */
  for (kdimsid=0; kdimsid<DIMSID_THERM3; kdimsid++)
  {
    for (idim=0; idim<NDIM_THERM3; idim++)
    {
      for (jdim=0; jdim<NDIM_THERM3; jdim++)
      {
        gamt[kdimsid][jdim] = sidredm[kdimsid][idim] * xjm[idim][jdim];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* build metric tensor */
  /* The metric tensor  g  (i.e. metm) is achieved by
   *    g = gamma^T gamma
   *      = [ x_,r^2 + y_,r^2 + z_,r^2   x_,r*x_,s + y_,r*y_,s + y_,r*y_,s ]
   *        [          sym               (x_,s)^2 + (y_,s)^2 + (z_,s)^2    ]
   */
  /* Tensor  g  might be swapped, but swapping does not affect its 
   * determinant,
   * instead of [ m_00  m_01 ] swapped means [ m_11 m_10 ]
   *            [ m_10  m_11 ]               [ m_01 m_00 ]
   */
  for (kdimsid=0; kdimsid<DIMSID_THERM3; kdimsid++)
  {
    for (jdimsid=0; jdim<DIMSID_THERM3; jdim++)
    {
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
        metm[kdimsid][jdimsid] = gamt[kdimsid][idim] * gamt[jdimsid][idim];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* determinat of metric tensor = squared scalar density */
  /* This result is equivalent to
   *    metr = || gamma_r x gamma_s ||
   * in which
   *    gamma_r = [ x_,r  y_,r  z_,r ]^T
   *    gamma_s = [ x_,s  y_,s  z_,s ]^T
   * and || () || is the Euclidean norm.
   */
  *metr = sqrt( metm[0][0] * metm[1][1] - metm[0][1] * metm[1][0] );

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_metr_surf */


/*======================================================================*/
/*!
\brief Determine metric at Gauss point by reducing the Jacobian
       matrix for element edge

\param ele     ELEMENT*     (i)  pointer to current element
\param nelelod INT          (i)  number of element nodes
\param deriv   DOUBLE**     (i)  derivatives of shape fct at Gauss point
\param sired   DOUBLE**     (i)  matrix for dimension reduction
\param metr    DOUBLE*      (o)  metric
\return void

\author bborn
\date 09/06
*/
void th3_metr_line(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3], 
                   DOUBLE linredv[NDIM_THERM3],
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian (dummy) */
  DOUBLE gamt[NDIM_THERM3];  /* differential */
  DOUBLE mets;  /* metric scalar */
  INT idim, jdim;  /* dimension indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_metr_line");
#endif

  /*--------------------------------------------------------------------*/
  /* get Jacobian matrix (xji is not written) */
  th3_metr_jaco(ele, nelenod, deriv, 0, xjm, &det, xji);

  /*--------------------------------------------------------------------*/
  /* build gamma matrix */
  for (idim=0; idim<NDIM_THERM3; idim++)
  {
    for (jdim=0; jdim<NDIM_THERM3; jdim++)
    {
      gamt[jdim] = linredv[idim] * xjm[idim][jdim];
    }
  }

  /*--------------------------------------------------------------------*/
  /* build squared metric scalar by inner product */
  mets = 0.0;
  for (idim=0; idim<NDIM_THERM3; idim++)
  {
    mets = mets + gamt[idim] * gamt[idim];
  }

  /*--------------------------------------------------------------------*/
  /* metric */
  *metr = sqrt(mets);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_metr_line */

/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/

