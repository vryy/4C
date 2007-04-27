/*======================================================================*/
/*!
\file
\brief Iterative update of internal variables

<pre>
Maintainer: Burkhard Bornemann 
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author mf
\date 10/06
*/
extern GENPROB genprob;



/*======================================================================*/
/*!
\brief Iterative update of element variables

The element stiffness matrix, ie the tangent operator, is determined
for the linear, 3dim elasticity

\param   *container     CONTAINER        (i)  see container.h
\param   *ele           ELEMENT          (i)  pointer to current element
\param   *gpshade       SO3_GPSHAPEDERIV (i)  Gauss point coords 
\param   *mat           MATERIAL         (i)
\param   *eforc_global  ARRAY            (o)  global vector for internal 
                                              forces (initialized!)
\param   *estif_global  ARRAY            (o)  element stiffness matrix
\param   *emass_global  ARRAY            (o)  element mass matrix
\return void

\author mf
\date 10/06
*/
void so3_iv_upd(const CONTAINER* container,
                ELEMENT* ele,
                const SO3_GPSHAPEDERIV* gpshade,
                const MATERIAL* mat)
{
  /* general variables/constants */
  INT nelenod;  /* numnp of this element */
  INT neledof;  /* total number of element DOFs */
  DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3];  /* material coord. of element */
  DOUBLE edis[MAXNOD_SOLID3][NDIM_SOLID3];  /* cur. element displacements */

  /* integration */
  INT ngp;  /* total number of Gauss points in domain */
  INT igp;  /* current total index of Gauss point */
  DOUBLE fac;  /* integration factors */
  INT inod;  /* node index */
  NODE *actnode;  /* pointer to current node */
  INT jdim; /* dimension index */

  /* quantities at Gauss point */
  SO3_GEODEFSTR gds;  /* isoparametric Jacobian, deformation grad, etc at
                       * Gauss point */

  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("so3_iv_upd");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;
  neledof = NUMDOF_SOLID3 * nelenod;
  for (inod=0; inod<nelenod; inod++)
  {
    actnode = ele->node[inod];
    for (jdim=0; jdim<NDIM_SOLID3; jdim++)
    {
      ex[inod][jdim] = actnode->x[jdim];
      /* THE FOLLOWING HARD-CODED `0' IS SHARED AMONG ALL STRUCTURE ELEMENTS
       * AND SOLUTION TECHNIQUES (BOTH STATICS, DYNAMICS AND FSI).
       * IT ACCESSES THE CURRENT NODAL DISPLACEMENTS STORED IN sol ARRAY.
       * THIS HARD-CODED `0' SHOULD BE REPLACED BY A SOFT-CODED VERSION.
       * NEW SOFT-CODED INDEX FOR OLD DISCRETISATION ==> array_position.h
       * NEW SOFT-CODED INDEX FOR NEW DISCRETISATION ==> TO BE ANNOUNCED
       * ==> IN FEM WE TRUST. */
      edis[inod][jdim] = actnode->sol.a.da[0][jdim];
    }
  }
  ngp = gpshade->gptot;  /* total number of Gauss points in domain */

  /*--------------------------------------------------------------------*/
  /* integration loop */
  for (igp=0; igp<ngp; igp++)
  {
    /*------------------------------------------------------------------*/
    /* Gauss point */
    gds.igp = igp;
    gds.gpc[0] = gpshade->gpco[igp][0];  /* r-coordinate */
    gds.gpc[1] = gpshade->gpco[igp][1];  /* s-coordinate */
    gds.gpc[2] = gpshade->gpco[igp][2];  /* t-coordinate */
    fac = gpshade->gpwg[igp];  /* weight */
    /*------------------------------------------------------------------*/
    /* compute Jacobian matrix, its determinant and inverse */
    so3_metr_jaco(ele, nelenod, ex, gpshade->gpderiv[igp], 1, 
                  gds.xjm, &(gds.xjdet), gds.xji);
    /*------------------------------------------------------------------*/
    /* integration (quadrature) factor */
    fac = fac * gds.xjdet;
    /*------------------------------------------------------------------*/
    /* deformation tensor and displacement gradient */
    so3_def_grad(nelenod, edis, gpshade->gpderiv[igp], gds.xji, 
                 gds.disgrdv, gds.defgrd);
    /*------------------------------------------------------------------*/
    /* Linear/Green-Lagrange strain vector */
    if (ele->e.so3->kintype == so3_geo_lin)
    {
      so3_strain_lin(ele, gds.disgrdv, gds.stnengv);
    }
    else if (ele->e.so3->kintype == so3_total_lagr)
    {
      so3_strain_gl(ele, gds.disgrdv, gds.stnglv);
    }
    else
    {
      dserror("Cannot digest chosen type of spatial kinematic\n");
    }
    /*------------------------------------------------------------------*/
    /* calculate B-operator
     * bop differs depending on geometrically linearity/non-linearity */
    so3_bop(ele, nelenod, gpshade->gpderiv[igp], gds.xji, 
            gds.defgrd, gds.bopn, gds.bop);
    /*------------------------------------------------------------------*/
    /* call material law and update MIVs */
    so3_mat_mivupd(container, ele, mat, igp, &gds);
  }

  /*--------------------------------------------------------------------*/
  /* local co-system */
  dsassert(ele->locsys == locsys_no,
           "locsys not implemented for this element!\n");

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_iv_upd() */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
