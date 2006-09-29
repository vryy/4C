/*======================================================================*/
/*!
\file
\brief (tangent) stiffness matrix, 
       mass matrix,
       internal forces of SOLID3 element

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15240
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
\brief Calculate linear stiffness matrix of SOLID3 element

The element stiffness matrix, ie the tangent operator, is determined
for the linear, 3dim elasticity

\param   *ele           ELEMENT     (i)  pointer to current element
\param   *data          SO3_DATA    (i)  common element data
\param   *estif_global  ARRAY       (o)  element stiffness matrix
\param   *emass_global  ARRAY       (o)  element mass matrix
\param   *force         DOUBLE      (o)  global vector for internal 
                                           forces (initialized!)
\return void

\author mf
\date 10/06
*/
void so3_lin_stiff(const ELEMENT *ele,
                   const SO3_DATA *data,
                         MATERIAL *mat,
                         ARRAY *estif_global,
                         ARRAY *emass_global,
                         DOUBLE *force)
{
  /* general variables/constants */
  INT nelenod;  /* numnp of this element */
  INT neledof;  /* total number of element DOFs */

  /* integration */
  INT igpr, igps, igpt;  /* Gauss point indices */
  INT gpnumr, gpnums, gpnumt;  /* Gauss point numbers */
  INT gpintcr, gpintcs, gpintct;  /* Gauss point integration cases */
  INT ip;  /* current total index of Gauss point */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpcr;  /* r-coord current GP */
  DOUBLE gpcs;  /* s-coord current GP */
  DOUBLE gpct;  /* t-coord current GP */

  /* quantities at Gauss point */
  DOUBLE shape[MAXNOD_SOLID3];  /* shape functions */
  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3];   /* shape fct. derivatives */
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];  /* inverse Jacobian matrix */
  DOUBLE bop[NDIM_SOLID3][NUMDOF_SOLID3*MAXNOD_SOLID3]; /* B-operator */
  DOUBLE cmat[NUMSTSS_SOLID3][NUMSTRN_SOLID3];  /* material matrix */
  DOUBLE stress[NUMSTSS_SOLID3];  /* stress */

  /* convenience */
  DOUBLE **estif;  /* element stiffness matrix */

  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("so3_lin_stiff");
#endif

  /*--------------------------------------------------------------------*/
  /* element matrix fields are reused for every element, thus
  /* have to be reinitialized to zero */
  amzero(estif_global);  /* element tangent matrix */
  estif = estif_global->a.da;
  amzero(emass_global); /* element mass matrix */
  emass = emass_global->a.da;

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;
  neledof = NUMDOF_SOLID3 * nelenod;

  /*--------------------------------------------------------------------*/
  /* Gauss integraton data */
  switch (ele->distyp)
  {
      /* hexahedra elements */
      case hex8: case hex20: case hex27:
          gpnumr = ele->e.so3->gpnum[0];
          gpintcr = ele->e.so3->gpintc[0];
          gpnums = ele->e.so3->gpnum[1];
          gpintcs = ele->e.so3->gpintc[1];
          gpnumt = ele->e.so3->gpnum[2];
          gpintct = ele->e.so3->gpintc[2];
          break;
      /* tetrahedra elements */
      /* tets are not simply rst-oriented and have just one GP-set nr */
      case tet4: case tet10:
          gpnumr = 1;
          gpnums = 1;
          gpnumt = ele->e.so3->gpnum[0];
          gpintcr = 1;
          gpintcs = 1;
          gpintct = ele->e.so3->gpintc[0];
          break;
      default:
          dserror("ele->distyp unknown!");
  }  /* end of switch(ele->distyp) */

  /*--------------------------------------------------------------------*/
  /* integration loops */
  for (igpr=0; igpr<gpnumr; igpr++)
  {
    for (igps=0; igps<gpnums; igps++)
    {
      for (igpt=0; igpt<gpnumt; igpt++)
      {
        /*--------------------------------------------------------------*/
        /* obtain current Gauss coordinates and weights */
        switch (ele->distyp)
        {
            /* hexahedra */
            case hex10: case hex20: case hex27:
              gpcr = data.ghlc[gpintcr][igpr];  /* r-coordinate */
              gpcs = data.ghlc[gpintcs][igps];  /* s-coordinate */
              gpct = data.ghlc[gpintct][igpt];  /* t-coordinate */
              fac = data.ghlw[gpintcr][igpr]  /* weight */
                  * data.ghlw[gpintcs][igps]
                  * data.ghlw[gpintct][igpt];
              break;
            /* tetrahedra */
            case tet4: case tet10:
              gpcr = data.gtdcr[gpintct][igpt];  /* r-coordinate */
              gpcs = data.gtdcs[gpintct][igpt];  /* s-coordinate */
              gpct = data.gtdct[gpintct][igpt];  /* t-coordinate */
              fac = data.gtdw[gpintct][igpt];  /* weight */
              break;
            default:
              dserror("ele->distyp unknown!");
        }  /* end of switch (ele->distyp) */
        /*--------------------------------------------------------------*/
        /* shape functions and their derivatives */
        so3_shape_deriv(ele->distyp, gpcr, gpcs, gpct, 1, shape, deriv);
        /*--------------------------------------------------------------*/
        /* compute Jacobian matrix, its determinant and inverse */
        so3_metr_jaco(ele, nelenod, deriv, 1, xjm, &det, xji);
        /*--------------------------------------------------------------*/
        /* integration (quadrature) factor */
        fac = fac * det;
        /*--------------------------------------------------------------*/
        /* calculate linear B-operator */
        so3_lin_bop(nelenod, deriv, xji, bop);
        /*--------------------------------------------------------------*/
        /* call material law */
        ip = (igpr+1) * (igps+1) * (igpt+1)  /* total Gauss point index */
        so3_mat_sel(ele, mat, bop, ip, stress, cmat);
        /*--------------------------------------------------------------*/
        /* element linear stiffness matrix estif add contribution at Gauss point
        so3_lin_stiff_bcb(neledof, bop, cmat, fac, estif);
        /*--------------------------------------------------------------*/
        /* element internal force from integration of stresses
        if (force == 1)
        {
          so3_lin_fint(stress, fac, bop, neledof, force);
        }
      }  /* end of for */
    }  /* end of for */
  }  /* end of for */

  /*--------------------------------------------------------------------*/
  /* local co-system */
  dsassert(ele->locsys == locsys_no,
           "locsys not implemented for this element!\n");

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_lin_stiff(...) */


/*======================================================================*/
/*!
\brief Classic B^T C B operation for one Gauss point
       Adds contribution of this Gauss point to element tangent

\param   neledof    INT     (i)    number of element DOFs
\param **bop        DOUBLE  (i)    B-operator
\param **cmat       DOUBLE  (i)    constitutive matrix
\param   fac        DOUBLE  (i)    integration factor of current GP
\param **tmat       DOUBLE  (io)   element stiffness matrix increment of
                                     current Gauss point
\return void

\author mf
\date 10/06
*/
void so2_lin_stiff_bcb(const INT       neledof,
                       const DOUBLE  **bop,
                       const DOUBLE  **cmat,
                       const DOUBLE    fac,
                       DOUBLE  **tmat)
{
  INT i, j, k, l, m;  /* counters */
  DOUBLE bopcmati, tkl;  /* intermediate sums */
  DOUBLE bopcmat[NUMSTRN_SOLID3];  /* bopcmat_ki = bop_kj * cmat_ji */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_lin_stiff_bcb");
#endif

  /*--------------------------------------------------------------------*/
  for (k=0; k<neledof; k++)
  {
    /* bopcmat_ki = bop_kj * cmat_ji */
    for (i=0; i<NUMSTRN_SOLID3; i++)
    {
      bopcmati = 0.0;
      for (j=0; j<NUMSTSS_SOLID3; j++)
      {
        bopcmati = bopcmati + bop[j][k]*cmat[j][i]*fac;
      }
      bopcmat[i] = bopcmati;
    }
    /* tmat_kl = bopcmat_ki * bop_il */
    for (l=0; l<neledof; l++)
    {
      tkl = 0.0;
      for (i=0; i<NUMSTRN_SOLID3; i++)
      {
        tkl = tkl + bopcmat[i]*bop[i][l];
      }
      tmat[k][l] = tmat[k][l] + tkl;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
   dstrc_exit();
#endif
   return;
} /* end of so3_lin_stiff_bcb(...) */




/*======================================================================*/
/*!
\brief Evaluate element nodal forces

The element nodal forces 'fie' are incremented by the contribution
of the current Gauss point

\param     neledof   INT      (i)    number of element nodes
\param   **bop       DOUBLE   (i)    B-operator for current GP
\param    *stress    DOUBLE   (i)    stress for current GP
\param     fac       DOUBLE   (i)    Gauss quadrature factor mult. etc.
\param    *intfor    DOUBLE   (io)   element internal force
\return void

\author mf
\date 03/06
*/
void so3_lin_fint(const INT      neledof,
                  const DOUBLE **bop,
                  const DOUBLE  *stress,
                  const DOUBLE   fac,
                        DOUBLE  *intfor)
{
  /*--------------------------------------------------------------------*/
  INT i, k;  /* counters */
  DOUBLE intfork;  /* stress multiplied by 'fac' */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_lin_fint");
#endif
  for (k=0; k<neledof; k++)
  {
    intfork = 0.0;
    for (i=0, i<NUMSTSS_SOLID3; i++)
    {
      intfork = intfork + bop[i][k]*stress[i]*fac;
    }
    intfor[k] += intfork;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_lin_fint(...) */


/*======================================================================*/
#endif /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
