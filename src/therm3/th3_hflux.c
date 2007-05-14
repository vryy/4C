/*======================================================================*/
/*!
\file
\brief Evaluate the element heat fluxes 

The element heat flux is evaluated at certain points in the element:
       (1) at the Gauss points
       (2) at the element nodes.
The element heat flux is calculated in different co-ordinate bases:
       (1) in global (X,Y,Z) components
       (2) in local (r,s,t) components
       (3) modulus and angles enclosed with X,Y,Z-axes

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 10/06
*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*----------------------------------------------------------------------*/
/*! 
\addtogroup THERM3
@{ (documentation module open)
*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

global variable GENPROB genprob is defined in global_control.c 

\author bborn
\date 03/06
*/
extern GENPROB genprob;

/*----------------------------------------------------------------------*/
/*!
\brief Fields
       vector of numfld FIELDs, defined in global_control.c
\author bborn
\date 05/07
*/
extern FIELD* field;

/*----------------------------------------------------------------------*/
/*!
\brief Locally globals (also known as dynamic statics)

\author bborn
\date 03/06
*/
static INT allocated = 0;  /* flag storing allocation status */
static ARRAY hflux_a;  /* temporary vector for heat flux at gauss point */
static DOUBLE *hflux;


/*======================================================================*/
/*!
\brief Allocate element heat fluxes and  working arrays

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 03/06
*/
void th3_hflux_init(PARTITION *actpart)
{
  PARTITION *tpart;
  INT i, j;  /* counters */
  ELEMENT *actele;  /* pointer to current element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_hflux_init");
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
      /* check if THERM3 element */
      if (actele->eltyp == el_therm3)
      {
        /* allocate heat flux arrays per element */
        am4def("hflux_gp_xyz", &(actele->e.th3->hflux_gp_xyz), 
               1, MAXGAUSS, NUMHFLX_THERM3, 0, "D3");
        am4def("hflux_gp_rst", &(actele->e.th3->hflux_gp_rst), 
               1, MAXGAUSS, NUMHFLX_THERM3, 0, "D3");
        am4def("hflux_gp_123", &(actele->e.th3->hflux_gp_123), 
               1, MAXGAUSS, 1+NUMHFLX_THERM3, 0, "D3");
        am4def("hflux_nd_xyz", &(actele->e.th3->hflux_nd_xyz), 
               1, MAXNOD_THERM3, NUMHFLX_THERM3, 0, "D3");
        am4def("hflux_nd_rst", &(actele->e.th3->hflux_nd_rst), 
               1, MAXNOD_THERM3, NUMHFLX_THERM3, 0, "D3");
        am4def("hflux_nd_123", &(actele->e.th3->hflux_nd_123), 
               1, MAXNOD_THERM3, 1+NUMHFLX_THERM3, 0, "D3");
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
  /* allocate working arrays */
  if (allocated == 0)
  {
    hflux = amdef("hflux", &hflux_a, NUMHFLX_THERM3, 1,"DV");
    /* locally allocatable arrays/vectors have been allocated */
    allocated = 1;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_hflux_init */


/*======================================================================*/
/*!
\brief Deallocate local globals

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 10/06
*/
void th3_hflux_final(PARTITION *actpart)
{
  PARTITION *tpart;
  ELEMENT *actele;
  INT i, j;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_hlux_final");
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
      /* check wether THERM3 element */
      if (actele->eltyp == el_therm3)
      {
        /* deallocate heat flux arrays at element */
        am4del(&(actele->e.th3->hflux_gp_xyz));
        am4del(&(actele->e.th3->hflux_gp_rst));
        am4del(&(actele->e.th3->hflux_gp_123));
        am4del(&(actele->e.th3->hflux_nd_xyz));
        am4del(&(actele->e.th3->hflux_nd_rst));
        am4del(&(actele->e.th3->hflux_nd_123));
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
  if (allocated == 1)
  {
    amdel(&hflux_a);
    /* reset allocation status */
    allocated = 0;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_hflux_final */


/*======================================================================*/
/*!
\brief Evaluate element heat flux

\author bborn
\date 10/06
*/
void th3_hflux_cal(CONTAINER *cont,
                   ELEMENT *ele,
                   TH3_DATA *data,
                   MATERIAL *mat)
{
  /* locator */
#ifdef D_TSI
  const ARRAY_POSITION_SOL* isol 
    = &(field[genprob.numtf].dis[cont->disnum_t].ipos.isol);
  const INT itemn = isol->temn;  /* curr. temperature index */
#else
  const INT itemn = 0;  /* curr. temperature index */
#endif
  const INT place = 0;

  INT nelenod = ele->numnp;  /* number of element nodes */
  INT neledof = NUMDOF_THERM3 * nelenod;  /* total number of element DOFs */
  DIS_TYP distyp = ele->distyp;  /* type of discretisation */
  const INT tdof = 0;
  DOUBLE ex[MAXNOD_THERM3][NDIM_THERM3];  /* material coord. of element */
  DOUBLE etem[MAXDOF_THERM3];  /* curr. element temperature vector */
  
  INT gpnumr = 0;  /* Gauss points in r-direction */
  INT gpnums = 0;  /* Gauss points in s-direction */
  INT gpnumt = 0;  /* Gauss points in t-direction */
  INT gpintcr = 0;  /* GP integration case in r-direction */
  INT gpintcs = 0;  /* GP integration case in s-direction */
  INT gpintct = 0;  /* GP integration case in t-direction */
  INT gpnum = 0;  /* total number of Gauss points in element domain */
  INT igpr, igps, igpt;  /* Gauss point index in r-, s- and t-direction */
  INT igp;  /* total index Gauss point */
  DOUBLE fac;  /* a factor */
  DOUBLE gpcr = 0.0;  /* GP r-coord */
  DOUBLE gpcs = 0.0;  /* GP s-coord */
  DOUBLE gpct = 0.0;  /* GP t-coord */
  INT inod;  /* nodal index */
  NODE* actnode;  /* pointer to current node */
  INT jdim;  /* dimension index */

  DOUBLE det;  /* Jacobi determinants */
  DOUBLE rst[NDIM_THERM3];  /* natural coordinate a point */
  DOUBLE shape[MAXNOD_THERM3];
  DOUBLE deriv[MAXNOD_THERM3][NDIM_THERM3];
  DOUBLE xjm[NDIM_THERM3][NDIM_THERM3];
  DOUBLE xji[NDIM_THERM3][NDIM_THERM3];
  DOUBLE bop[NUMTMGR_THERM3][NUMDOF_THERM3*MAXNOD_THERM3];
  DOUBLE cmat[NUMHFLX_THERM3][NUMTMGR_THERM3];
  DOUBLE tmgr[NUMTMGR_THERM3];  /* temp. gradient */
  DOUBLE hflux[NUMHFLX_THERM3];  /* heat flux */

  INT ihflx;  /* heat flux index */
  INT inode;  /* element nodal index */

  /*--------------------------------------------------------------------*/
  /* start */
  /* Working arrays (locally globals) MUST be initialised! */ 
#ifdef DEBUG
  dstrc_enter("th3_hflux_cal");
#endif

  /*--------------------------------------------------------------------*/
  /* number of Gauss points */
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
      /* total number of Gauss points */
      gpnum = gpnumr * gpnums * gpnumt;
      break;
    /* tetrahedra elements */
    case tet4: case tet10:
      gpnumr = 1;
      gpnums = 1;
      gpnumt = ele->e.th3->gpnum[0];
      gpintcr = 1;
      gpintcs = 1;
      gpintct = ele->e.th3->gpintc[0];
      /* total number of Gauss points */
      gpnum = gpnumt;
      break;
    default:
      dserror("distyp unknown!");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
  /* element geometry and temperature */
  for (inod=0; inod<nelenod; inod++)
  {
    actnode = ele->node[inod];
    for (jdim=0; jdim<NDIM_THERM3; jdim++)
    {
      ex[inod][jdim] = actnode->x[jdim];
    }
    etem[inod] = actnode->sol.a.da[itemn][tdof];  /* inod==idof */
  }

  /*--------------------------------------------------------------------*/
  /* heat flux at every Gauss point */
  /* init total Gauss point counter */
  igp = 0;
  /* parse all Gauss points and evaluate stress */
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
        /* factor */
        fac = det;
        /*--------------------------------------------------------------*/
        /* calculate B-operator */
        th3_bop(nelenod, deriv, xji, bop);
        /*--------------------------------------------------------------*/
        /* temperature gradient */
        th3_lin_temgrad(ele, bop, etem, tmgr);
        /*--------------------------------------------------------------*/
        /* call material law */
        th3_mat_sel(cont, ele, mat, igp, tmgr, hflux, cmat);
        /*--------------------------------------------------------------*/
        /* store heat flux of current Gauss point */
        for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
        {
          ele->e.th3->hflux_gp_xyz.a.d3[place][igp][ihflx] = hflux[ihflx];
        }
        /*--------------------------------------------------------------*/
        /* store heat flux in parameter space co-ordinates (r,s,t) */
        th3_hflux_rst(xjm, hflux, igp, ele->e.th3->hflux_gp_rst.a.d3[place]);
        /*--------------------------------------------------------------*/
        /* store modulus and angles at current Gauss point */
        th3_hflux_modang(hflux, igp, ele->e.th3->hflux_gp_123.a.d3[place]);
        /*--------------------------------------------------------------*/
        /* increment total Gauss point counter */
        igp++;
      }  /* end of for */
    }  /* end of for */
  }  /* end of for */

  /*--------------------------------------------------------------------*/
  /* Heat fluxes at Gauss points are extrapolated to element nodes. */
  /* loop all element nodes */
  for (inode=0; inode<nelenod; inode++)
  {
    /*------------------------------------------------------------------*/
    /* get local coordinates of node ==> rs */
    th3_cfg_noderst(ele, data, inode, rst);
    /*------------------------------------------------------------------*/
    /* extrapolate values now */
    /* heat flux in x-direction : q_1 */
    th3_hflux_extrpol(ele, data, gpnum, 
                      ele->e.th3->hflux_gp_xyz.a.d3[place], rst,
                      hflux);
    for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
    {
      ele->e.th3->hflux_nd_xyz.a.d3[place][inode][ihflx] = hflux[ihflx];
    }
    /* absolute heat flux and its direction at element nodes */
    th3_hflux_modang(hflux, inode, ele->e.th3->hflux_nd_123.a.d3[place]);
  }  /* end for */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_hflux_cal */


/*======================================================================*/
/*!
\brief Convert heat flux with respect to (X,Y,Z) to (r,s,t)

\param **xjm         DOUBLE   (i)   Jacobi matrix at GP
\param  *hflux       DOUBLE   (i)   heat flux (at GP)
\param   igp         INT      (i)   index of current Gauss point
\param **hflux123    DOUBLE   (o)   modulus and angles at GP
\return void

\author bborn
\date 10/06
*/
void th3_hflux_rst(DOUBLE xjm[NDIM_THERM3][NDIM_THERM3],
                   DOUBLE *hflux,
                   INT igp,
                   DOUBLE **hfluxrst)
{
  INT idim, jdim;
  DOUBLE hfluc;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_hflux_rst");
#endif

  /*--------------------------------------------------------------------*/
  dsassert(NDIM_THERM3==NUMHFLX_THERM3,
           "Warning dimension of heat flux vector differs to problem "
           "dimension\n");

  /*--------------------------------------------------------------------*/
  /* get locally defined heat flux */
  for (idim=0; idim<NDIM_THERM3; idim++)
  {
    hfluc = 0.0;
    for (jdim=0; jdim<NDIM_THERM3; jdim++)
    {
      /* we have to use the 'transposed Jacobian' here,
       * referred to the isoparametric Jacobian notation
       * commonly used, cf. th3_metr.c */
      hfluc = hfluc + xjm[jdim][idim] * hflux[jdim];
    }
    hfluxrst[igp][idim] = hfluc;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void void th3_hflux_rst */

/*======================================================================*/
/*!
\brief Modulus of heat flux and its direction measured in angles with
       respect to the global X,Y,Z-axes

\param  *hflux       DOUBLE   (i)   heat flux (at GP)
\param   igp         INT      (i)   index of current Gauss point
\param **hflux123    DOUBLE   (o)   modulus and angles at GP
\return void

\author bborn
\date 10/06
*/
void th3_hflux_modang(DOUBLE *hflux,
                      INT igp,
                      DOUBLE **hflux123)

{
  INT ihf;  /* count heat flux components */
  DOUBLE hfluxmod;  /* heat flux modulus */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_hflux_modang");
#endif

  /*--------------------------------------------------------------------*/
  /* absolute value of heat flux (modulus) */
  hfluxmod = 0.0;
  for (ihf=0; ihf<NUMHFLX_THERM3; ihf++)
  {
    hfluxmod += hflux[ihf]*hflux[ihf];
  }
  hfluxmod = sqrt(hfluxmod);
  hflux123[igp][NUMHFLX_THERM3] = hfluxmod;  /* store in ELEMENT array */

  /*--------------------------------------------------------------------*/
  /* direction of heat flux with respect to X,Y,Z-axes */
  if (hfluxmod != 0.0)
  {
    for (ihf=0; ihf<NUMHFLX_THERM3; ihf++)
    {
      hflux123[igp][ihf] = acos(hflux[ihf]/hfluxmod);
    }
  }
  else
  {
    for (ihf=0; ihf<NUMHFLX_THERM3; ihf++)
    {
      hflux123[igp][ihf] = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_hflux_modang */


/*======================================================================*/
/*!
\brief Extrapolate heat fluxes at Gauss points to element nodes

These procedure does not provide more accurate heat fluxes as the direct
computation, but maybe saves a bit of computing time.

\param  *ele      ELEMENT     (i)   pointer to active element
\param  *data     TH3_DATA    (i)   pointer to THERM3 data (GPs coords etc)
\param  ngauss    INT         (i)   total number of Gauss points
\param  **hfluxgp DOUBLE      (i)   heat flux at Gauss points
\param  *rst      DOUBLE      (i)   element node natural coordinates
\param  *hfluxnd  DOUBLE      (o)   extrapolated heat flux at node
\return void

\author bborn
\date 03/06
*/
void th3_hflux_extrpol(ELEMENT *ele,
                       TH3_DATA *data,
                       INT ngauss,
                       DOUBLE **hfluxgp,
                       DOUBLE *rst,
                       DOUBLE *hfluxnd)
{
  INT gpnumr = 0;  /* Gauss points in r-direction */
  INT gpnums = 0;  /* Gauss points in s-direction */
  INT gpnumt = 0;  /* Gauss points in t-direction */
  INT gpintcr = 0;  /* GP integration case in r-direction */
  INT gpintcs = 0;  /* GP integration case in s-direction */
  INT gpintct = 0;  /* GP integration case in t-direction */
  INT gpnum = 0;  /* total number of Gauss points in element domain */
  INT igpr, igps, igpt;  /* Gauss point loop counter */
  INT igauss;  /* total Gauss point counter */
  DOUBLE gpcr, gpcs, gpct;  /* GP (r,s,t) co-ordinate */

  INT ihflx;  /* heat flux index */

  INT i, j, k;  /* loop counter */
  DOUBLE gri, gsj, gtk;  /* Lagrange stages */
  DOUBLE funinc[NUMHFLX_THERM3];  /* increment due to contribution of cur. GP */

  /*--------------------------------------------------------------------*/
  /* start up */
#ifdef DEBUG
  dstrc_enter("th3_hflux_extrpol");
#endif

  /*--------------------------------------------------------------------*/
  /* number of Gauss points */
  switch (ele->distyp)
  {
    /* hexahedra elements */
    case hex8: case hex20: case hex27:
      gpnumr = ele->e.th3->gpnum[0];
      gpintcr = ele->e.th3->gpintc[0];
      gpnums = ele->e.th3->gpnum[1];
      gpintcs = ele->e.th3->gpintc[1];
      gpnumt = ele->e.th3->gpnum[2];
      gpintct = ele->e.th3->gpintc[2];
      /* total number of Gauss points */
      gpnum = gpnumr * gpnums * gpnumt;
      break;
    /* tetrahedra elements */
    case tet4: case tet10:
      gpnumr = 1;
      gpnums = 1;
      gpnumt = ele->e.th3->gpnum[0];
      gpintcr = 1;
      gpintcs = 1;
      gpintct = ele->e.th3->gpintc[0];
      /* total number of Gauss points */
      gpnum = gpnumt;
      break;
    default:
      dserror("distyp unknown!");
  }  /* end switch */
  /*--------------------------------------------------------------------*/
  /* check total number of Gauss points */
  dsassert(gpnum==ngauss, "Total number of Gauss points do not match!");

  /*--------------------------------------------------------------------*/
  /* initialise nodal heat flux vector */
  for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
  {
    hfluxnd[ihflx] = 0.0;
  }

  /*--------------------------------------------------------------------*/
  /* extrapolation */
  switch (ele->distyp)
  {
    /*------------------------------------------------------------------*/
    /* hexahedron element */
    case hex8: case hex20: case hex27:
      /* initialise total Gauss point index/counter */
      igauss = 0;
      /* loop Gauss points in r-direction */
      for (igpr=0; igpr<gpnumr; igpr++)
      {
        /* r-coordinate of current Gauss point */
        gpcr = data->ghlc[gpintcr][igpr];
        /* loop Gauss points in s-direction */
        for (igps=0; igps<gpnums; igps++)
        {
          /* s-coordinate of current Gauss point */
          gpcs = data->ghlc[gpintcs][igps];
          /* loop Gauss points in t-direction */
          for (igpt=0; igpt<gpnumt; igpt++)
          {
            /* t-coordinate of current Gauss point */
            gpct = data->ghlc[gpintct][igpt];
            /* determine func increment due to extrapolation function
             * at Gauss point */
            /* initialise increment/contribution of current Gauss point */
            /* IMPORTANT: The calling sequence of the nir*nis Gauss points
             *            must be indentically as in the superroutine
             *            th2_hflux_cal, as the heat fluxes are stored
             *            in a vector
             */
            for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
            {
              funinc[ihflx] = hfluxgp[igauss][ihflx];
            }
            /* build tri-directional Lagrange polynomials at Gauss points */
            /* l_{igpr,igps,igpt}(r,s,t) = hflux_{igpr,igps,igpt}
             *                           * l_{igpr}^{gpintcr}(r)
             *                           * l_{igps}^{gpintcs}(s)
             *                           * l_{igps}^{gpintct}(t)
             * with 
             * l_{igpr}^{gpintcr}(r) 
             *      = \prod_{i=0,i!=igpr}^{gpintcr} (r - r_i)/(r_igpr - r_i)
             * l_{igps}^{gpintcs}(s) 
             *      = \prod_{i=0,i!=igps}^{gpintcs} (s - s_i)/(s_igps - s_i)
             * l_{igpt}^{gpintct}(t) 
             *      = \prod_{i=0,i!=igpt}^{gpintct} (t - t_i)/(t_igpt - t_i)
             */
            /* in r-direction : l_{igpr}^{gpintcr}(r) */
            for (i=0; i<gpnumr; i++)
            {
              if (i != igpr)
              {
                gri = data->ghlc[gpintcr][i];
                for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
                {
                  funinc[ihflx] = funinc[ihflx] 
                    * (rst[0] - gri)/(gpcr - gri);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* in s-direction : l_{igps}^{gpintcs}(s) */
            for (j=0; j<gpnums; j++)
            {
              if (j != igps)
              {
                gsj = data->ghlc[gpintcs][j];
                for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
                {
                  funinc[ihflx] = funinc[ihflx] 
                    * (rst[1] - gsj)/(gpcs - gsj);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* in t-direction : l_{igpt}^{gpintct}(t) */
            for (k=0; k<gpnumt; k++)
            {
              if (k != igpt)
              {
                gtk = data->ghlc[gpintct][k];
                for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
                {
                  funinc[ihflx] = funinc[ihflx] 
                    * (rst[2] - gtk)/(gpct - gtk);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* add increment of current Gauss point */
            for (ihflx=0; ihflx<NUMHFLX_THERM3; ihflx++)
            {
              hfluxnd[ihflx] = hfluxnd[ihflx] + funinc[ihflx];
            }
            /* increment total Gauss point counter */
            igauss++;
          }  /* end for */
        }  /* end for */
      }  /* end for */
      break;
    /*------------------------------------------------------------------*/
    /* tetrahedron elements */
    case tet4: case tet10:
      dserror("Extrapolation is not available to tet elements!");
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
  return;
}  /* end void th3_hflux_extrpol */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
/*! @} (documentation module close)*/
