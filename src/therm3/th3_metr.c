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
                   DOUBLE **deriv,
                   INT      flag,
                   DOUBLE **xjm,
                   DOUBLE  *det,
                   DOUBLE **xji)
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
    for (j=0; i<NDIM_THERM3; i++)
    {
      xjm[i][j] = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
  /* loop all element nodes / shape functions and dimensions */
  /* inod = ele->numnp; */
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
\param sired   DOUBLE**     (i)  matrix for dimension reduction
\param metr    DOUBLE*      (o)  metric
\return void

\author bborn
\date 09/06
*/
void th3_metr_surf(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE **deriv, 
                   DOUBLE **sidred,
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian matrix */
  DOUBLE gamt[DIMSID_THERM3][NDIM_THERM3];  /* differential */
  DOUBLE metm[DIMSID_THERM3][DIMSID_THERM3];  /* metric matrix */
  INT idim, jdim, kdimsid;  /* dimension indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_metr_surf");
#endif

  /*--------------------------------------------------------------------*/
  /* get Jacobian matrix (xji is flagged off)
   * The Jacobian matrix looks like:
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
   *     [ N_,x ]^k           [ N_,r ]
   *     [ N_,y ]   =  J^{-1} [ N_,s ]
   *     [ N_,z ]             [ N_,t ]
   */
  th3_metr_jaco(ele, nelenod, deriv, 0, xjm, &det, xji);

  /*--------------------------------------------------------------------*/
  /* build gamma tensor */
  for (kdimsid=0; kdimsid<DIMSID_THERM3; kdimsid++)
  {
    for (idim=0; idim<NDIM_THERM3; idim++)
    {
      for (jdim=0; jdim<NDIM_THERM3; jdim++)
      {
        gamt[kdimsid][jdim] = sidred[kdimsid][idim] * xjm[idim][jdim];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* build metric tensor --- might be swapped, but it does not affect
   * its determinant
   * instead of [ m_00  m_01 ] swapped means [ m_11 m_10 ]
   *            [ m_10  m_11 ]               [ m_01 m_00 ] */
  for (kdimsid=0; kdimsid<DIMSID_THERM3; kdimsid++)
  {
    for (jdim=0; jdim<DIMSID_THERM3; jdim++)
    {
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
        metm[kdimsid][jdim] = gamt[kdimsid][idim] * gamt[jdim][idim];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* determinat of metric tensor -- metric */
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
                   DOUBLE **deriv, 
                   DOUBLE *linredv,
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian matrix */
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
  /* build metric scalar */
  for (idim=0; idim<NDIM_THERM3; idim++)
  {
    mets = gamt[idim] * gamt[idim];
  }

  /*--------------------------------------------------------------------*/
  /* determinat of metric tensor -- metric */
  *metr = mets;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_metr_line */

/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/

