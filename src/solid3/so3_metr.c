/*======================================================================*/
/*!
\file
\brief Metrics of SOLID3 element, i.e. metric conversion between physical
       and parameter space. 

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/frenzel
	    089-289-15240
</pre>

\author mf
\date 10/06
*/
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
\brief Jacobian matrix at coordinate (r,s,t) (e.g. Gauss point) by
       applying isoparametric concept

\param ele      ELEMENT*    (i)   the element
\param enod     INT         (i)   number of element nodes
\param deriv    DOUBLE[]    (i)   derivatives of the shape functions
\param flag     INT         (i)   flag for inverse of Jacobian
                                     1 : inverse Jacobian
                                     else : don't inverse Jacobian
\param xjm      DOUBLE[][]  (o)   the Jacobian matrix
\param det      DOUBLE*     (o)   determinant of the Jacobian matrix
\param xji      DOUBLE[][]  (o)   inverse of Jacobian
\return void

\author mf
\date 10/06
*/
void so3_metr_jaco(ELEMENT *ele,
                   INT      enod,
                   DOUBLE   ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE   deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                   INT      flag,
                   DOUBLE   xjm[NDIM_SOLID3][NDIM_SOLID3],
                   DOUBLE  *det,
                   DOUBLE   xji[NDIM_SOLID3][NDIM_SOLID3])
{
  INT i, j, k;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_metr_jaco");
#endif
  
  /*--------------------------------------------------------------------*/
  /* initialise Jacobian to zero */
  for (i=0; i<NDIM_SOLID3; i++)
  {
    for (j=0; j<NDIM_SOLID3; j++)
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
    for (i=0; i<NDIM_SOLID3; i++)
    {
      for (j=0; j<NDIM_SOLID3; j++)
      {
         xjm[j][i] += deriv[k][j] * ex[k][i];
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* determinat */
  /* determinat is calculated by Sarrus' rule 
   * [Bronstein, Taschenbuch der Mathematik */
  *det = xjm[0][0] * xjm[1][1] * xjm[2][2]
       + xjm[0][1] * xjm[1][2] * xjm[2][0]
       + xjm[0][2] * xjm[1][0] * xjm[2][1]
       - xjm[0][2] * xjm[1][1] * xjm[2][0]
       - xjm[0][0] * xjm[1][2] * xjm[2][1]
       - xjm[0][1] * xjm[1][0] * xjm[2][2];
  /* check if negative Jacobian */
  if (*det == 0.0)
  {
    dserror("ZERO JACOBIAN DETERMINANT");
  }
  else if (*det < 0.0)
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
    for (i=0; i<NDIM_SOLID3; i++)
    {
      for (j=0; j<NDIM_SOLID3; j++)
      {
        xji[i][j] = xji[i][j] / (*det);
      }
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_metr_jaco */


/*======================================================================*/
/*!
\brief Determine metric at Gauss point by reducing the Jacobian
       matrix 
       Determine physical orientation of parametric base vectors 
       of surface at Gauss point
       This can be either the material/initial/reference configuration
       or the spatial/deformed/current configuration depending on
       the kind of node coordinates.

\param ele     ELEMENT*     (i)  pointer to current element
\param nelelod INT          (i)  number of element nodes
\param ex      DOUBLE[][]   (i)  element node coordinates
\param deriv   DOUBLE**     (i)  derivatives of shape fct at Gauss point
\param sired   DOUBLE**     (i)  matrix for dimension reduction
\param gamtt   DOUBLE[][]   (o)  physically oriented parametric base vectors
                                 (optional)
\param metr    DOUBLE*      (o)  metric
\return void

\author mf
\date 10/06
*/
void so3_metr_surf(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE sidredm[DIMSID_SOLID3][NDIM_SOLID3],
                   DOUBLE gamtt[DIMSID_SOLID3][NDIM_SOLID3],
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];  /* inverse Jacobian matrix */
  DOUBLE gamt[DIMSID_SOLID3][NDIM_SOLID3];  /* differential */
  DOUBLE metm[DIMSID_SOLID3][DIMSID_SOLID3];  /* metric matrix */
  INT idim, jdim, kdimsid, jdimsid;  /* dimension indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_metr_surf");
#endif

  /*--------------------------------------------------------------------*/
  /* get Jacobian matrix xjm (xji is not computed) */
  so3_metr_jaco(ele, nelenod, ex, deriv, 0, xjm, &det, NULL);

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
   * matrix (xjm) defined in FE-fashion (see so3_metr_jaco).
   * sidredm  depends on the surface parameter space (r,s). 
   */
  /* set to zero */
  memset(gamt, 0, sizeof(gamt));
  /* matrix product */
  for (kdimsid=0; kdimsid<DIMSID_SOLID3; kdimsid++)
  {
    for (idim=0; idim<NDIM_SOLID3; idim++)
    {
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        gamt[kdimsid][jdim] += sidredm[kdimsid][idim] * xjm[idim][jdim];
      }
    }
  }
  /* set output */
  if (gamtt != NULL)
  {
    for (kdimsid=0; kdimsid<DIMSID_SOLID3; kdimsid++)
    {
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        gamtt[kdimsid][jdim] = gamt[kdimsid][jdim];
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
  /* set to zero */
  memset(metm, 0, sizeof(metm));
  /* matrix product */
  for (kdimsid=0; kdimsid<DIMSID_SOLID3; kdimsid++)
  {
    for (jdimsid=0; jdimsid<DIMSID_SOLID3; jdimsid++)
    {
      for (idim=0; idim<NDIM_SOLID3; idim++)
      {
        metm[kdimsid][jdimsid] += gamt[kdimsid][idim] * gamt[jdimsid][idim];
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
}  /* end of so3_metr_surf */


/*======================================================================*/
/*!
\brief Determine metric at Gauss point by reducing the Jacobian
       matrix for element edge

\param ele     ELEMENT*     (i)  pointer to current element
\param nelelod INT          (i)  number of element nodes
\param ex      DOUBLE[][]   (i)  element node coordinates
\param deriv   DOUBLE[][]   (i)  derivatives of shape fct at Gauss point
\param sired   DOUBLE[][]   (i)  matrix for dimension reduction
\param metr    DOUBLE*      (o)  metric
\return void

\author mf
\date 10/06
*/
void so3_metr_line(ELEMENT *ele, 
                   INT nelenod, 
                   DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                   DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3], 
                   DOUBLE linredv[NDIM_SOLID3],
                   DOUBLE *metr)
{
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];  /* Jacobian matrix */
  DOUBLE det;  /* determinat of Jacobian matrix */
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];  /* inverse Jacobian matrix */
  DOUBLE gamt[NDIM_SOLID3];  /* differential */
  DOUBLE mets;  /* metric scalar */
  INT idim, jdim;  /* dimension indices */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_metr_line");
#endif

  /*--------------------------------------------------------------------*/
  /* get Jacobian matrix (xji is not written) */
  so3_metr_jaco(ele, nelenod, ex, deriv, 0, xjm, &det, xji);

  /*--------------------------------------------------------------------*/
  /* build gamma matrix */
  /* set to zero */
  memset(gamt, 0, sizeof(gamt));
  /* matrix product */
  for (idim=0; idim<NDIM_SOLID3; idim++)
  {
    for (jdim=0; jdim<NDIM_SOLID3; jdim++)
    {
      gamt[jdim] += linredv[idim] * xjm[idim][jdim];
    }
  }

  /*--------------------------------------------------------------------*/
  /* build squared metric scalar by inner product */
  mets = 0.0;
  for (idim=0; idim<NDIM_SOLID3; idim++)
  {
    mets += gamt[idim] * gamt[idim];
  }

  /*--------------------------------------------------------------------*/
  /* determinat of metric tensor -- metric */
  *metr = sqrt(mets);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_metr_line */


/*======================================================================*/
/*!
\brief Determine rotational content in Jacobian matrix

\param xjm     DOUBLE[][]   (i)  Jacobian matrix J
\param xrm     DOUBLE[][]   (o)  rotational content in J
\param xrvm    DOUBLE[][]   (o)  rotational content in J prepared
                                 such that symmetric tensors
                                 stored vectorially can be
                                 rotated
\param         DOUBLE[][]   (o)  inverse rotational content in J
                                 prepared such that symmetric tensors
                                 stored vectorially can be
                                 rotated
\return void

\author bborn
\date 01/07
*/
void so3_metr_rot(DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE xrm[NDIM_SOLID3][NDIM_SOLID3],
                  DOUBLE xrvm[NUMSTR_SOLID3][NUMSTR_SOLID3],
                  DOUBLE xrvi[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  DOUBLE r00, r01, r02, r10, r11, r12, r20, r21, r22;  /* components
                                                        * of intermediate
                                                        * rotation matrix */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_metr_rot");
#endif

  /*--------------------------------------------------------------------*/
  /* get rotational `content' R in isoparametric Jacobian
   * WARNING: Jacobian is stored in classical FE-fashion, 
   *          i.e. it is transposed if viewed as `deformation gradient'
   *          describing the metric mapping of parameter to  material
   *          frame.
   *          Thus rotation matrix is transposed
   *          (and left and right stretch are switched) */
  so3_tns3_plrdcmp(xjm, xrm, NULL, NULL);

  /*--------------------------------------------------------------------*/
  /* rotation matrix R in components (accounting for transposedness) */
  r00 = xrm[0][0];  r10 = xrm[0][1];  r20 = xrm[0][2];
  r01 = xrm[1][0];  r11 = xrm[1][1];  r21 = xrm[1][2];
  r02 = xrm[2][0];  r12 = xrm[2][1];  r22 = xrm[2][2];

  /*--------------------------------------------------------------------*/
  /* Rotational component of J is prepared such that 
   * a symmetric 2-tensor A referring to the parameter space (r,s,t)
   *         [ A_11  A_12  A_31 ]
   *     A = [ A_12  A_22  A_23 ]
   *         [ A_31  A_23  A_33 ]
   * denoted vectorially Av, i.e.
   *     Av^T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]
   * can be computed by matrix-vector product
   *     (A)_{XYZ} = R^T . (A)_{rst} . R
   * becomes
   *     (Av)_{XYZ} = xrvm . (Av)_{rst} */
  xrvm[0][0] = r00*r00;
  xrvm[0][1] = r10*r10;
  xrvm[0][2] = r20*r20;
  xrvm[0][3] = 2.0*r00*r10;
  xrvm[0][4] = 2.0*r10*r20;
  xrvm[0][5] = 2.0*r00*r20;
  xrvm[1][0] = r01*r01;
  xrvm[1][1] = r11*r11;
  xrvm[1][2] = r21*r21;
  xrvm[1][3] = 2.0*r01*r11;
  xrvm[1][4] = 2.0*r11*r21;
  xrvm[1][5] = 2.0*r01*r21;
  xrvm[2][0] = r02*r02;
  xrvm[2][1] = r12*r12;
  xrvm[2][2] = r22*r22;
  xrvm[2][3] = 2.0*r02*r12;
  xrvm[2][4] = 2.0*r12*r22;
  xrvm[2][5] = 2.0*r02*r22;
  /* ~~~ */
  xrvm[3][0] = r00*r01;
  xrvm[3][1] = r10*r11;
  xrvm[3][2] = r20*r21;
  xrvm[3][3] = r00*r11 + r10*r01;
  xrvm[3][4] = r10*r21 + r20*r11;
  xrvm[3][5] = r00*r21 + r20*r01;
  xrvm[4][0] = r01*r02;
  xrvm[4][1] = r11*r12;
  xrvm[4][2] = r21*r22;
  xrvm[4][3] = r01*r12 + r11*r02;
  xrvm[4][4] = r11*r22 + r21*r12;
  xrvm[4][5] = r01*r22 + r21*r02;
  xrvm[5][0] = r00*r02;
  xrvm[5][1] = r10*r12;
  xrvm[5][2] = r20*r22;
  xrvm[5][3] = r00*r12 + r10*r02;
  xrvm[5][4] = r10*r22 + r20*r12;
  xrvm[5][5] = r00*r22 + r20*r02;

  /*--------------------------------------------------------------------*/
  /* The rotational component of J^{-1} is R^{-1} in which R^{-1}
   * is the rotation matrix computed for the inverse map. The inverse
   * of R is R^T.
   * Thus we achieve for a 2-tensor B referring to the material space (XYZ)
   *         [ B_11  B_12  B_13 ]
   *     B = [ B_21  B_22  B_23 ]
   *         [ B_31  B_32  B_33 ]
   * rotated in parameter space directions
   *     (B)_{rst} = R^{-T} . (B)_{XYZ} . R^{-1}
   * introducing R^{-1} = R^T = xrm results in
   *     (B)_{rst} = R . (B)_{XYZ} . R^T */
  /* R^{-1} = R^T is set component by component */
  r00 = xrm[0][0];  r01 = xrm[0][1];  r02 = xrm[0][2];
  r10 = xrm[1][0];  r11 = xrm[1][1];  r12 = xrm[1][2];
  r20 = xrm[2][0];  r21 = xrm[2][1];  r22 = xrm[2][2];

  /*--------------------------------------------------------------------*/
  /* Rotational component R^{-1} of J^{-1} is prepared such that 
   * the mapping 
   *     (B)_{rst} = R^{-T} . (B)_{XYZ} . R^{-1}
   * of a symmetric 2-tensor B referring to the material 
   * frame (X,Y,Z), 
   *         [ B_11  B_12  B_31 ]
   *     B = [ B_12  B_22  B_23 ]
   *         [ B_31  B_23  B_33 ]
   * is achieved by a matrix-vector product
   *     (Bv)_{XYZ} = xrvi . (Bv)_{XYZ}
   * with vectorially denoted Bv, i.e.
   *     Bv^T = [ B_11  B_22  B_33  B_12  B_23  B_31 ] */
  xrvi[0][0] = r00*r00;
  xrvi[0][1] = r10*r10;
  xrvi[0][2] = r20*r20;
  xrvi[0][3] = 2.0*r00*r10;
  xrvi[0][4] = 2.0*r10*r20;
  xrvi[0][5] = 2.0*r00*r20;
  xrvi[1][0] = r01*r01;
  xrvi[1][1] = r11*r11;
  xrvi[1][2] = r21*r21;
  xrvi[1][3] = 2.0*r01*r11;
  xrvi[1][4] = 2.0*r11*r21;
  xrvi[1][5] = 2.0*r01*r21;
  xrvi[2][0] = r02*r02;
  xrvi[2][1] = r12*r12;
  xrvi[2][2] = r22*r22;
  xrvi[2][3] = 2.0*r02*r12;
  xrvi[2][4] = 2.0*r12*r22;
  xrvi[2][5] = 2.0*r02*r22;
  /* ~~~ */
  xrvi[3][0] = r00*r01;
  xrvi[3][1] = r10*r11;
  xrvi[3][2] = r20*r21;
  xrvi[3][3] = r00*r11 + r10*r01;
  xrvi[3][4] = r10*r21 + r20*r11;
  xrvi[3][5] = r00*r21 + r20*r01;
  xrvi[4][0] = r01*r02;
  xrvi[4][1] = r11*r12;
  xrvi[4][2] = r21*r22;
  xrvi[4][3] = r01*r12 + r11*r02;
  xrvi[4][4] = r11*r22 + r21*r12;
  xrvi[4][5] = r01*r22 + r21*r02;
  xrvi[5][0] = r00*r02;
  xrvi[5][1] = r10*r12;
  xrvi[5][2] = r20*r22;
  xrvi[5][3] = r00*r12 + r10*r02;
  xrvi[5][4] = r10*r22 + r20*r12;
  xrvi[5][5] = r00*r22 + r20*r02;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/

