/*======================================================================*/
/*!
\file
\brief Tangent (stiffness) matrix, 
       capacity (mass) matrix,
       nodal heat fluxes (internal forces) of THERM3 element

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*----------------------------------------------------------------------*/
/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/06
*/
extern GENPROB genprob;



/*======================================================================*/
/*!
\brief Calculate tangent matrix of linear heat conduction THERM3
       element

The element stiffness matrix, ie the tangent operator, is determined
for the linear, 3dim heat conduction

\param   *ele           ELEMENT     (i)  pointer to current element
\param   *data          TH3_DATA    (i)  common element data
\param   *estif_global  ARRAY       (o)  element tangent matrix
\param   *emass_global  ARRAY       (o)  element mass matrix
\param   *force         DOUBLE      (o)  global vector for internal 
                                           forces (initialized!)
\return void

\author bborn
\date 09/06
*/
void th3_lin_tang(CONTAINER *cont,
                  ELEMENT *ele,
                  TH3_DATA *data,
                  MATERIAL *mat,
                  ARRAY *estif_global,
                  ARRAY *emass_global,
                  DOUBLE *force)
{
  /* general variables/constants */
  INT nelenod;  /* numnp of this element */
  INT neledof;  /* total number of element DOFs */
  DIS_TYP distyp;  /* type of discretisation */

  /* integration */
  INT igpr, igps, igpt;  /* Gauss point indices */
  INT gpnumr, gpnums, gpnumt;  /* Gauss point numbers */
  INT gpintcr, gpintcs, gpintct;  /* Gauss point integration cases */
  INT igp;  /* current total index of Gauss point */
  DOUBLE fac;  /* integration factors */
  DOUBLE gpcr;  /* r-coord current GP */
  DOUBLE gpcs;  /* s-coord current GP */
  DOUBLE gpct;  /* t-coord current GP */

  /* quantities at Gauss point */
  DOUBLE shape[MAXNOD_THERM3];  /* shape functions */
  DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3];   /* shape fct. derivatives */
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];  /* inverse Jacobian matrix */
  DOUBLE bop[NUMTMGR_THERM3][NUMDOF_THERM3*MAXNOD_THERM3]; /* B-operator */
  DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3];  /* conductivity matrix */
  DOUBLE hflux[NUMHFLX_THERM3];  /* heat flux */

  /* convenience */
  DOUBLE **estif;  /* element stiffness matrix */

  /*--------------------------------------------------------------------*/
  /* start */
#ifdef DEBUG
  dstrc_enter("th3_lin_tang");
#endif

  /*--------------------------------------------------------------------*/
  /* some of the fields have to be reinitialized to zero */
  amzero(estif_global);  /* element tangent matrix */
  estif = estif_global->a.da;
  /* element capacity (mass) matrix is not dealt with --- yet */
  if (emass_global != NULL)
  {
    amzero(emass_global);
  }

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;
  neledof = NUMDOF_THERM3 * nelenod;
  distyp = ele->distyp;

  /*--------------------------------------------------------------------*/
  /* Gauss integraton data */
  switch (distyp)
  {
    /* hexahedra elements */
    case hex8: case hex20: case hex27:
      gpnumr = ele->e.th3->gpnum[0];
      gpintcr = ele->e.th3->gpintc[0];
      gpnums = ele->e.th3->gpnum[1];
      gpintcs = ele->e.th3->gpintc[1];
      gpnumt = ele->e.th3->gpnum[2];
      gpintct = ele->e.th3->gpintc[2];
      break;
    /* tetrahedra elements */
    case tet4: case tet10:
      gpnumr = 1;
      gpnums = 1;
      gpnumt = ele->e.th3->gpnum[0];
      gpintcr = 1;
      gpintcs = 1;
      gpintct = ele->e.th3->gpintc[0];
      break;
    default:
      dserror("distyp unknown!");
  }  /* end of switch(distyp) */

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
        switch (distyp)
        {
          /* hexahedra */
          case hex8: case hex20: case hex27:
            gpcr = data->ghlc[gpintcr][igpr];  /* r-coordinate */
            gpcs = data->ghlc[gpintcs][igps];  /* s-coordinate */
            gpct = data->ghlc[gpintct][igpt];  /* t-coordinate */
            fac = data->ghlw[gpintcr][igpr]  /* weight */
                * data->ghlw[gpintcs][igps]
                * data->ghlw[gpintct][igpt];
            break;
          /* tetrahedra */
          case tet4: case tet10:
            gpcr = data->gtdc[gpintct][igpt][0];  /* r-coordinate */
            gpcs = data->gtdc[gpintct][igpt][1];  /* s-coordinate */
            gpct = data->gtdc[gpintct][igpt][2];  /* t-coordinate */
            fac = data->gtdw[gpintct][igpt];  /* weight */
            break;
          default:
            dserror("distyp unknown!");
        }  /* end of switch (distyp) */
        /*--------------------------------------------------------------*/
        /* shape functions and their derivatives */
        th3_shape_deriv(distyp, gpcr, gpcs, gpct, 1, shape, deriv);
        /*--------------------------------------------------------------*/
        /* compute Jacobian matrix, its determinant and inverse */
        th3_metr_jaco(ele, nelenod, deriv, 1, xjm, &det, xji);
        /*--------------------------------------------------------------*/
        /* integration (quadrature) factor */
        fac = fac * det;
        /*--------------------------------------------------------------*/
        /* calculate B-operator */
        th3_bop(nelenod, deriv, xji, bop);
        /*--------------------------------------------------------------*/
        /* call material law */
        th3_mat_sel(cont, ele, mat, bop, igp, hflux, cmat);
        /*--------------------------------------------------------------*/
        /* element tangent matrix estif add contribution at Gauss point
         * (like element stiffness matrix) */
        th3_lin_bcb(neledof, bop, cmat, fac, estif);
        /*--------------------------------------------------------------*/
        /* element nodal heat flux from integration of heat fluxes
         * (like element internal forces) */
        if (force != NULL)
        {
          th3_lin_fint(neledof, bop, hflux, fac, force);
        }
        /*--------------------------------------------------------------*/
        /* increment total Gauss point index */
        igp++;
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
} /* end of th3_lin_tang(...) */


/*======================================================================*/
/*!
\brief Classic B^T C B operation for one Gauss point
       Adds contribution of this Gauss point to element tangent

\param   neledof    INT     (i)    number of element DOFs
\param **bop        DOUBLE  (i)    B-operator
\param **cmat       DOUBLE  (i)    constitutive matrix
\param   fac        DOUBLE  (i)    integration factor of current GP
\param **tmat       DOUBLE  (io)   element tangent matrix increment of
                                     current Gauss point
\return void

\author bborn
\date 03/06
*/
void th3_lin_bcb(INT       neledof,
                 DOUBLE    bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                 DOUBLE    cmat[NUMHFLX_THERM3][NUMTMGR_THERM3],
                 DOUBLE    fac,
                 DOUBLE  **tmat)
{
  INT i, j, k, l;  /* counters */
  DOUBLE bopcmati, tkl;  /* intermediate sums */
  DOUBLE bopcmat[NUMTMGR_THERM3];  /* bopcmat_ki = bop_kj * cmat_ji */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_lin_bcb");
#endif

  /*--------------------------------------------------------------------*/
  for (k=0; k<neledof; k++)
  {
    /* bopcmat_ki = bop_kj * cmat_ji */
    for (i=0; i<NUMTMGR_THERM3; i++)
    {
      bopcmati = 0.0;
      for (j=0; j<NUMHFLX_THERM3; j++)
      {
        bopcmati = bopcmati + bop[j][k]*cmat[j][i]*fac;
      }
      bopcmat[i] = bopcmati;
    }
    /* tmat_kl = bopcmat_ki * bop_il */
    for (l=0; l<neledof; l++)
    {
      tkl = 0.0;
      for (i=0; i<NUMTMGR_THERM3; i++)
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
} /* end of th3_lin_bcb(...) */



/*======================================================================*/
/*!
\brief Evaluate element nodal forces

The element nodal forces 'fie' are incremented by the contribution
of the current Gauss point

\param     neledof   INT      (i)    number of element nodes
\param   **bop       DOUBLE   (i)    B-operator for current GP
\param    *hflux     DOUBLE   (i)    heat flux for current GP
\param     fac       DOUBLE   (i)    Gauss quadrature factor mult. etc.
\param    *intfor    DOUBLE   (io)   element internal force
\return void

\author bborn
\date 03/06
*/
void th3_lin_fint(INT      neledof,
                  DOUBLE   bop[NDIM_THERM3][NUMDOF_THERM3*MAXNOD_THERM3],
                  DOUBLE   hflux[NUMHFLX_THERM3],
                  DOUBLE   fac,
                  DOUBLE  *intfor)
{
  /*--------------------------------------------------------------------*/
  INT i, k;  /* counters */
  DOUBLE intfork;  /* heatflux multiplied by 'fac' */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_lin_fint");
#endif

  /*--------------------------------------------------------------------*/
  /* geometrically linear */
  for (k=0; k<neledof; k++)
  {
    intfork = 0.0;
    for (i=0; i<NUMHFLX_THERM3; i++)
    {
      intfork = intfork + bop[i][k]*hflux[i]*fac;
    }
    intfor[k] += intfork;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_lin_fint(...) */


/*======================================================================*/
#endif /* end of #ifdef D_THERM3 */
/*! @} (documentation module close) */
