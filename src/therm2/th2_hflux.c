/*======================================================================*/
/*!
\file
\brief contains the routine 'w1_cal_stress' which evaluates the element
       stresses for 2D isoparametric degenerated element
       contains the routine 'w1_mami' which evaluates the principal
       stresses and directions for 2D isoparametric degenerated element
       contains the routine 'w1rsn' which returns R/S coordinates of
       gauss integration points 2D isoparametric degenerated element
       contains the routine 'w1recs' which extrapolates from gauss
       points for rectangles

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
#include "../headers/standardtypes.h"
#include "therm2.h"


/*! \addtogroup THERM2 */
/*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable GENPROB genprob is defined in global_control.c

\author bborn
\date 03/06
*/
extern struct _GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Locally globals (also known as dynamic statics)

\author bborn
\date 03/06
*/
static INT allocated = 0;  /* flag storing allocation status */
static ARRAY    cmat_a;      /* material tensor */
static DOUBLE **cmat;
static ARRAY    shape_a;  /* shape functions */
static DOUBLE  *shape;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    xjm0_a;    /* jacobian matrix at r,s=0*/
static DOUBLE **xjm0;
static ARRAY    xji_a;    /* inverse of jacobian matrix */
static DOUBLE **xji;
static ARRAY    hflux_a;      /* dummy matrix for saving stresses at gauss point*/
static DOUBLE  *hflux;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;


/*======================================================================*/
/*!
\brief Allocate element heat fluxes and  working arrays

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 03/06
*/
void th2_hflux_init(PARTITION *actpart)
{
  PARTITION *tpart;
  INT i, j;  /* counters */
  ELEMENT *actele;  /* pointer to current element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_hflux_init");
#endif

  /*--------------------------------------------------------------------*/
  /* allocate heat flux arrays stored at each element */
  tpart = actpart;  /*&(actpart[genprob.numtf]);*/
  /* loop over all discretisations of partition thermal field */
  for (j=0; j<tpart->ndis; j++)
  {
    /* loop over all elements of current discretisation */
    for (i=0; i<tpart->pdis[j].numele; i++)
    {
      /* set current element */
      actele = tpart->pdis[j].element[i];
      /* check wether THERM2 element */
      if (actele->eltyp == el_therm2)
      {
        /* allocate heat flux arrays per element */
        am4def("hflux_gp", &(actele->e.th2->hflux_gp),
               1, 3*NUMHFLX_THERM2, MAXGAUSS, 0, "D3");
        am4def("flux_nd",&(actele->e.th2->hflux_nd),
               1, 3*NUMHFLX_THERM2, MAXNOD, 0, "D3");
      }  /* end of if (actele->eltyp == el_therm2) */
    }  /* end of for (i=0; i<actpart->pdis[j].numele; i++) */
  }  /* end of for (j=0; j<actpart->ndis; j++) */

  /*--------------------------------------------------------------------*/
  /* allocate working arrays */
  if (allocated == 0)
  {
    shape = amdef("shape", &shape_a, MAXNOD_THERM2, 1, "DV");
    deriv = amdef("deriv", &deriv_a, NDIM_THERM2, MAXNOD_THERM2, "DA");
    cmat = amdef("cmat", &cmat_a, NUMTMGR_THERM2, NUMTMGR_THERM2, "DA");
    xjm = amdef("xjm", &xjm_a, NDIM_THERM2, NDIM_THERM2, "DA");
    xjm0 = amdef("xjm0", &xjm0_a, NDIM_THERM2, NDIM_THERM2, "DA");
    xji = amdef("xji", &xji_a, NDIM_THERM2, NDIM_THERM2, "DA");
    hflux = amdef("hflux", &hflux_a, NUMHFLX_THERM2, 1,"DV");
    bop = amdef("bop", &bop_a,
                NUMTMGR_THERM2, NUMDOF_THERM2*MAXNOD_THERM2, "DA");
    /* locally allocatable arrays/vectors have been allocated */
    allocated = 1;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_hflux_init */


/*======================================================================*/
/*!
\brief Deallocate locally globals

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 03/06
*/
void th2_hflux_final(PARTITION *actpart)
{
  PARTITION *tpart;
  ELEMENT *actele;
  INT i, j;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_hlux_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate heat flux arrays */
  tpart = actpart; /*&(actpart[genprob.numtf]); */
  /* loop over all discretisations of partition */
  for (j=0; j<tpart->ndis; j++)
  {
    /* loop over all elements of current discretisation */
    for (i=0; i<tpart->pdis[j].numele; i++)
    {
      /* set current element */
      actele = tpart->pdis[j].element[i];
      /* check wether THERM2 element */
      if (actele->eltyp == el_therm2)
      {
        /* deallocate heat flux arrays at element */
        am4del(&(actele->e.th2->hflux_gp));
        am4del(&(actele->e.th2->hflux_nd));
      }  /* end of if (actele->eltyp == el_therm2) */
    }  /* end of for (i=0; i<actpart->pdis[j].numele; i++) */
  }  /* end of for (j=0; j<actpart->ndis; j++) */

  /*--------------------------------------------------------------------*/
  if (allocated == 1)
  {
    amdel(&shape_a);
    amdel(&deriv_a);
    amdel(&cmat_a);
    amdel(&xjm_a);
    amdel(&xjm0_a);
    amdel(&xji_a);
    amdel(&hflux_a);
    amdel(&bop_a);
    /* reset alloacation status */
    allocated = 0;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_hflux_final */


/*======================================================================*/
/*!
\brief Evaluate element heat flux

2-D isoparametric degenerated element

\author bborn
\date 03/06
*/
void th2_hflux_cal(ELEMENT *ele,
                   TH2_DATA *data,
                   MATERIAL *mat,
                   INT kstep)  /* number of current load step */
{

  /* constants */
  const INT numdf = NUMDOF_THERM2;
  INT i;  /* some loopers */
  INT nir = 0;  /* num GP in r/s/t direction */
  INT nis = 0;  /* num GP in r/s/t direction */
  INT lr, ls;  /* loopers over GP */
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* total number of element DOFs */
  INT ip;  /* total index Gauss point */
  INT intc = 0;  /* "integration case" for tri-element */
  INT newval = 1;  /* controls evaluation of new stresses    */

  DOUBLE fac;  /* a factor */
  DOUBLE e1 = 0.0;  /* GP r-coord */
  DOUBLE e2 = 0.0;  /* GP s-coord */
  DOUBLE rs[NDIM_THERM2];  /* natural coordinate a point */

  DOUBLE hfluxnd;  /* heat flux component (temporary variable) */

  DOUBLE det;  /* Jacobi determinants */
  INT inode;  /* element nodal index */
  INT npoint;  /* total number of Gauss points in element domain */

  /*--------------------------------------------------------------------*/
  /* start */
  /* Working arrays (locally globals) MUST be initialised! */
#ifdef DEBUG
  dstrc_enter("th2_hflux_cal");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;
  neledof = numdf * nelenod;
  /*--------------------------------------------------------------------*/
  /* integration parameters */
  switch (ele->distyp)
  {
    case quad4: case quad8: case quad9:  /* --> quad - element */
      nir = ele->e.th2->nGP[0];
      nis = ele->e.th2->nGP[1];
      break;
    case tri3: case tri6: /* --> tri - element */
      nir  = 1;
      nis  = ele->e.th2->nGP[0];
      intc = ele->e.th2->gpintc;
      break;
    default:
      dserror("ele->distyp unknown!");
  }  /* end switch(ele->distyp) */
  /* total number of Gauss points */
  npoint  = nir*nis;

  /*====================================================================*/
  /* integration loops */
  ip = 0;  /* init total Gauss point counter */
  for (lr=0; lr<nir; lr++)
  {
    for (ls=0; ls<nis; ls++)
    {
      /*----------------------------------------------------------------*/
      /* get values of  shape functions and their derivatives */
      switch(ele->distyp)
      {
        /* --> quad - element */
        case quad4: case quad8: case quad9:
          e1   = data->gqlc[lr][nir-1];  /* r-ccord of cur. GP */
          e2   = data->gqlc[ls][nis-1];  /* s-ccord of cur. GP */
          break;
        /* --> tri - element */
        case tri3: case tri6:
          e1   = data->gtdcr[ls][intc];
          e2   = data->gtdcs[ls][intc];
      break;
      default:
         dserror("ele->distyp unknown!");
      } /* end switch(ele->distyp) */

      /*----------------------------------------------------------------*/
      /* shape functions and their derivatives (flagged by last 1) */
      th2_shape_deriv(shape, deriv, e1, e2, ele->distyp, 1);
      /*----------------------------------------------------------------*/
      /* compute Jacobian matrix and its determinant */
      th2_jaco(deriv, xjm, &det, ele, nelenod);
      /*----------------------------------------------------------------*/
      /* integration (quadrature) factor */
      fac = det;
      /*----------------------------------------------------------------*/
      /* calculate operator B */
      amzero(&bop_a);
      th2_bop(bop, deriv, xjm, det, nelenod);
      /*----------------------------------------------------------------*/
      /* call material law
       * to obtain heat fluxes 'hflux' globally xy-oriented */
      newval = 1;  /* Flag to calculate stresses */
      th2_mat_sel(ele, mat, bop, ip, hflux, cmat);
      /*----------------------------------------------------------------*/
      /* transform globally xy-oriented heat fluxes
       * to local rs-orientation */
/*       switch(ele->e.th2->hfluxtype) */
/*       { */
/*         case  th2_hflux_gprs: case th2_hflux_ndrs: */
/*           w1_tram(xjm, transm, transmi, work); */
/*           w1_lss(hflux, transmi, transm, it); */
/*           break; */
/*         default: */
/*           break; */
/*       } */
      /*----------------------------------------------------------------*/
      /* store heat flux of current Gauss point */
      for (i=0; i<NUMHFLX_THERM2; i++)
      {
	ele->e.th2->hflux_gp.a.d3[kstep][i][ip] = hflux[i];
      }
      /*----------------------------------------------------------------*/
      /* absolute heat flux plus orientation angle */
      th2_hflux_steep(hflux,
                      &(ele->e.th2->hflux_gp.a.d3[kstep][3][ip]),
                      &(ele->e.th2->hflux_gp.a.d3[kstep][4][ip]),
                      &(ele->e.th2->hflux_gp.a.d3[kstep][5][ip]));
      /*----------------------------------------------------------------*/
      /* increment total GP counter */
      ip++;
    }  /* end of for (ls=0; ls<nis; ls++) */
  }  /* end of for (lr=0; lr<nir; lr++) */


  /*--------------------------------------------------------------------*/
  /* Heat fluxes at Gauss points are extrapolated to element nodes. */
  /* loop all element nodes */
  for (inode=0; inode<nelenod; inode++)
  {
    /*------------------------------------------------------------------*/
    /* get local coordinates of node ==> rs */
    th2_cfg_noders(ele, inode, rs);
    /*------------------------------------------------------------------*/
    /* extrapolate values now */
    /* heat flux in x-direction : q_1 */
    th2_hflux_extrpol(ele, data, npoint,
                      &(ele->e.th2->hflux_gp.a.d3[kstep][0][0]), rs,
                      &(hfluxnd));
    ele->e.th2->hflux_nd.a.d3[kstep][0][inode] = hfluxnd;
    /* heat flux in y-direction : q_2 */
    th2_hflux_extrpol(ele, data, npoint,
                      &(ele->e.th2->hflux_gp.a.d3[kstep][1][0]), rs,
                      &(hfluxnd));
    ele->e.th2->hflux_nd.a.d3[kstep][1][inode] = hfluxnd;
    /* heat flux in z-direction : q_3 */
    th2_hflux_extrpol(ele, data, npoint,
                      &(ele->e.th2->hflux_gp.a.d3[kstep][2][0]), rs,
                      &(hfluxnd));
    ele->e.th2->hflux_nd.a.d3[kstep][2][inode] = hfluxnd;

    /* absolute heat flux and its direction at element nodes */
    hflux[0] = ele->e.th2->hflux_nd.a.d3[kstep][0][inode];
    hflux[1] = ele->e.th2->hflux_nd.a.d3[kstep][1][inode];
    hflux[2] = ele->e.th2->hflux_nd.a.d3[kstep][2][inode];
    th2_hflux_steep(hflux,
                    &(ele->e.th2->hflux_nd.a.d3[kstep][3][inode]),
                    &(ele->e.th2->hflux_nd.a.d3[kstep][4][inode]),
                    &(ele->e.th2->hflux_nd.a.d3[kstep][5][inode]));
  }  /* end of for (inode=0; inode<nelenod; inode++) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_hflux_cal */


/*======================================================================*/
/*!
\brief Modulus of heat flux and its direction of steepest increase/decrease

\param *hflux       DOUBLE   (i)   heat flux (at GP)
\param *hfluxmod    DOUBLE   (o)   modulus of heat flux
\param *hfluxang    DOUBLE   (o)   angle of direction of steepest
                                     increase/decrease
                                     angle in [-pi/2,pi/2] with respect
                                     to x-axis
\param *dum         DOUBLE   (o)   dummy for future purposes
\return void

\author bborn
\date 03/06
*/
void th2_hflux_steep(DOUBLE *hflux,
                     DOUBLE *hfluxmod,
                     DOUBLE *hfluxang,
                     DOUBLE *dum)

{
  INT ihf;  /* count heat flux components */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_hflux_steep");
#endif

  /*--------------------------------------------------------------------*/
  /* absolute value of heat flux (modulus) */
  *hfluxmod = 0.0;
  for (ihf=0; ihf<NUMHFLX_THERM2; ihf++)
  {
    *hfluxmod += hflux[ihf]*hflux[ihf];
  }
  *hfluxmod = sqrt(*hfluxmod);

  /*--------------------------------------------------------------------*/
  /* direction of heat flux */
  if (hflux[0] != 0.0)
  {
    if (hflux[0] > 0.0)
    {
      *hfluxang = atan(hflux[1]/hflux[0]);
    }
    else
    {
      if (hflux[1] >= 0.0)
      {
        *hfluxang = PI + atan(hflux[1]/hflux[0]);
      }
      else
      {
        *hfluxang = -PI + atan(hflux[1]/hflux[0]);
      }
    }
  }
  else
  {
    if (hflux[1] > 0.0)
    {
      *hfluxang = PI/2.0;
    }
    else if (hflux[1] < 0.0)
    {
      *hfluxang = -PI/2.0;
    }
    else
    {
      *hfluxang = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
  /* unused dummy */
  *dum = 0.0;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th2_hflux_steep */


/*======================================================================*/
/*!
\brief Extrapolate heat fluxes at Gauss points to element nodes

These procedure does not provide more accurate heat fluxes as the direct
computation, but maybe saves a bit of computing time.

\param  *ele     ELEMENT     (i)   pointer to active element
\param  *data    TH2_DATA    (i)   pointer to THERM2 data (GPs coords etc)
\param  ngauss   INT         (i)   total number of Gauss points
\param  *hfluxgp DOUBLE      (i)   heat flux at Gauss points
\param  *rs      DOUBLE      (i)   element node natural coordinates
\param  *hfluxnd DOUBLE      (o)   extrapolated heat flux at node
\return void

\author bborn
\date 03/06
*/
void th2_hflux_extrpol(ELEMENT *ele,
                       TH2_DATA *data,
                       INT ngauss,
                       DOUBLE *hfluxgp,
                       DOUBLE *rs,
                       DOUBLE *hfluxnd)
{
  INT i, j;  /* loop counter */
  INT lr, ls;  /* Gauss point loop counter */
  INT nir=0, nis=0;  /* Gauss point numbers in r-/s-direction */
  INT intc;  /* Gauss piont integration set for tri elements */
  INT igauss;  /* total Gauss point counter */
  DOUBLE grlr;  /* r-coord of current Gauss point */
  DOUBLE gsls;  /* s-ccord of current Gauss point */
  DOUBLE gri;  /* Lagrange stage */
  DOUBLE gsj;  /* Lagrange stage */
  DOUBLE funinc;  /* increment due to contribution of cur. GP */

  /*--------------------------------------------------------------------*/
  /* start up */
#ifdef DEBUG
  dstrc_enter("th2_hflux_extrpol");
#endif

  /*--------------------------------------------------------------------*/
  /* number of Gauss points in each direction */
  switch (ele->distyp)
  {
    /* --> quad - element */
    case quad4: case quad8: case quad9:
      nir = ele->e.th2->nGP[0];
      nis = ele->e.th2->nGP[1];
      break;
    /* --> tri - element */
    case tri3: case tri6:
      nir  = 1;
      nis  = ele->e.th2->nGP[0];
      intc = ele->e.th2->gpintc;
      break;
    default:
      dserror("ele->distyp unknown!");
  }  /* end switch(ele->distyp) */
  /* check total number of Gauss points */
  if (nir*nis != ngauss)
  {
    dserror("Total number of Gauss points do not match!");
  }

  /*--------------------------------------------------------------------*/
  /* extrapolation */
  switch (ele->distyp)
  {
    /*------------------------------------------------------------------*/
    /* quad - element */
    case quad4: case quad8: case quad9:
      /* initialise extrapolated value */
      *hfluxnd = 0.0;
      /* initialise total Gauss point index/counter */
      igauss = 0;
      /* loop Gauss points in r-direction */
      for (lr=0; lr<nir; lr++)
      {
        /* r-coordinate of current Gauss point */
        grlr = data->gqlc[lr][nir-1];
        /* loop Gauss points in s-direction */
        for (ls=0; ls<nis; ls++)
        {
          /* s-coordinate of current Gauss point */
          gsls = data->gqlc[ls][nis-1];
          /* determine func increment due to extrapolation function
           * at Gauss point */
          /* initialise increment/contribution of current Gauss point */
          /* IMPORTANT: The calling sequence of the nir*nis Gauss points
           *            must be indentically as in the superroutine
           *            th2_hflux_cal, as the heat fluxes are stored
           *            in a vector
           */
          funinc = hfluxgp[igauss];
          /* build bi-directional Lagrange polynomials at Gauss points */
          /* l_{lr,ls}(r,s) = sig_{lr,ls}
           *                * l_{lr}^{nir-1}(r)
           *                * l_{ls}^{nis-1}(s)
           * with l_{lr}^{nir-1}(r) = \prod_{i=0,i!=lr}^{nir-1}
           *                            (r - r_i)/(r_lr - r_i)
           * and  l_{ls}^{nis-1}(s) = \prod_{j=0,j!=ls}^{nis-1}
           *                            (s - s_i)/(s_lr - s_i)
           */
          /* in r-direction : l_{lr}^{nir-1}(r) */
          for (i=0; i<nir; i++)
          {
            if (i != lr)
            {
              gri = data->gqlc[i][nir-1];
              funinc = funinc * (rs[0] - gri)/(grlr - gri);
            }
          }
          /* in s-direction : l_{ls}^{nis-1}(s) */
          for (j=0; j<nis; j++)
          {
            if (j != ls)
            {
              gsj = data->gqlc[j][nis-1];
              funinc = funinc * (rs[1] - gsj)/(gsls - gsj);
            }
          }
          /* add increment of current Gauss point */
          *hfluxnd = *hfluxnd + funinc;
          /* increment total Gauss point counter */
          igauss++;
        }  /* end of for (ls=0; ls<nis; ls++) */
      }  /* end of for (lr=0; lr<nir; lr++) */
      break;
    /*------------------------------------------------------------------*/
    /* tri - element */
    case tri3: case tri6:
      dserror("Extrapolation is not available to tri elements!");
      break;
    /*------------------------------------------------------------------*/
    /* catch unfitting elements */
    default:
      dserror("Discretisation type is impossible!");
  }  /* end of switch (ele->distyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return ;
}  /* end of th2_hflux_extrpol */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM2 */
/*! @} (documentation module close)*/
#endif
