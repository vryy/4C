/*======================================================================*/
/*!
\file
\brief Update of internal variables

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

/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 03/06
*/
#ifdef D_TSI
extern FIELD *field;
#endif

/*======================================================================*/
/*!
\brief Test whether an update of internal variables is required
\author bborn
\date 05/07
*/
void so3_iv_updreq(ELEMENT* ele,
                   MATERIAL* actmat,
                   INT* updreq)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_iv_updreq");
#endif

  /*--------------------------------------------------------------------*/
  /* initialise */
  *updreq = 0;

  /*--------------------------------------------------------------------*/
  /* check material internal variables */
  so3_mat_mivupdreq(ele, actmat, updreq);

  /*--------------------------------------------------------------------*/
  /* Possibly, other update requests for internal variables will have to be
   * checked here, e.g., EAS parameters. */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Calculate linear stiffness matrix of SOLID3 element

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
void so3_iv_upditer(CONTAINER* container,
                    ELEMENT* ele,
                    SO3_GPSHAPEDERIV* gpshade,
                    MATERIAL* mat)
{
  /* general variables/constants */
  const INT nelenod = ele->numnp;  /* number of element nodes */
  const INT neledof = NUMDOF_SOLID3 * nelenod;  /* total element DOFs */

  /* indices in NODE sol* fiels */
  INT idisn;
  INT idisres;
#ifdef D_TSI
  if (genprob.probtyp == prb_tsi)
  {
    const ARRAY_POSITION_SOL* isol
      = &(field[genprob.numsf].dis[container->disnum_s].ipos.isol);
    idisn = isol->disn;
    const ARRAY_POSITION_SOLRES* isolres
      = &(field[genprob.numsf].dis[container->disnum_s].ipos.isolres);
    idisres = isolres->disres;
  }
  else
#endif
  {
    /* The following hard-coded `0' is shared among all structure elements
     * and solution techniques (both statics, dynamics and FSI).
     * It accesses the current NODE displacements stored in its sol array.
     * This hard-coded `0' should be replaced by a soft-coded version. */
    idisn = 0;
    idisres = 0;
  }

  DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3];  /* material coord. of element */
  DOUBLE edis[MAXNOD_SOLID3][NDIM_SOLID3];  /* cur. element displacements */
  DOUBLE edisii[MAXDOF_SOLID3];  /* iterative/residual displacements */

  /* integration */
  const INT ngp = gpshade->gptot;  /* total number of Gauss points in domain */
  INT igp;  /* current total index of Gauss point */
  DOUBLE fac;  /* integration factors */
  INT inod;  /* node index */
  NODE* actnode;  /* pointer to current node */
  INT jdim; /* dimension index */

  /* quantities at Gauss point */
  SO3_GEODEFSTR gds;  /* isoparametric Jacobian, deformation grad, etc at
                       * Gauss point */
  DOUBLE epsii[NUMSTR_SOLID3];  /* iterative/residual total strain */


  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("so3_iv_upditer");
#endif

  /*--------------------------------------------------------------------*/
  /* element coordinates and displacements */
  for (inod=0; inod<nelenod; inod++)
  {
    actnode = ele->node[inod];
    for (jdim=0; jdim<NDIM_SOLID3; jdim++)
    {
      /* nodal coordinates */
      ex[inod][jdim] = actnode->x[jdim];
      /* displacement at current step \inc\D_{n+1}^<i> */
      edis[inod][jdim] = actnode->sol.a.da[idisn][jdim];
      /* build iterative/residual displacements \iinc\D_{n+1}^<i> */
      edisii[inod*NDIM_SOLID3+jdim] 
        = actnode->sol_residual.a.da[idisres][jdim];
    }
  }
  

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
    /* iterative strains */
    {
      INT istr;
      for (istr=0; istr<NUMSTR_SOLID3; istr++)
      {
        DOUBLE rcsum = 0.0;
        INT jdof;
        for (jdof=0; jdof<neledof; jdof++)
        {
          rcsum += gds.bop[istr][jdof] * edisii[jdof];
        }
        epsii[istr] = rcsum;
      }
    }
    /*------------------------------------------------------------------*/
    /* call material law and update MIVs */
    so3_mat_mivupditer(container, ele, mat, igp, &gds, epsii);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Incremental update of element variables

\param   *container     CONTAINER        (i)  see container.h
\param   *ele           ELEMENT          (i)  pointer to current element
\param   *mat           MATERIAL         (i)  element material
\return void

\author mf
\date 05/07
*/
void so3_iv_updincr(const CONTAINER* container,
                ELEMENT* ele,
                const MATERIAL* mat)
{

  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("so3_iv_updincr");
#endif

  /*--------------------------------------------------------------------*/
  /* material internal variables */
  so3_mat_mivupdincr(container, ele, mat);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_iv_upd() */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/




