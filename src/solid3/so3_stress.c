/*======================================================================*/
/*!
\file
\brief Evaluate the element stresses

The element stresses are evaluated at certain points in the element:
       (1) at the Gauss points
       (2) at the element nodes.
The element stresses are calculated in different co-ordinate bases:
       (1) in global (X,Y,Z) components
       (2) in local (r,s,t) components
       (3) principal components and angles enclosed with X,Y,Z-axes

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 12/06
*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*----------------------------------------------------------------------*/
/*! 
\addtogroup SOLID3
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
\brief filewise globals (also known as "dynamic statics")

\author bborn
\date 03/06
*/
static ARRAY stress_a;  /* temporary vector for stresses at gauss point */
static DOUBLE *stress;


/*======================================================================*/
/*!
\brief Allocate element heat fluxes and  working arrays

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 03/06
*/
void so3_stress_init(PARTITION *actpart)
{
  PARTITION *tpart;
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT *actele;  /* pointer to current element */
  SOLID3 *actso3;  /* pointer to current SOLID3 element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_stress_init");
#endif
  
  /*--------------------------------------------------------------------*/
  /* allocate heat flux arrays stored at each element */
  tpart = actpart;  /*&(actpart[genprob.numtf]);*/
  /* loop over all discretisations of partition thermal field */
  for (jdis=0; jdis<tpart->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<tpart->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = tpart->pdis[jdis].element[iele];
      /* check if SOLID3 element */
      if (actele->eltyp == el_solid3)
      {
        /* set pointer to SOLID3 */
        actso3 = actele->e.so3;
        /* allocate heat flux arrays per element */
        am4def("stress_gpxyz", &(actso3->stress_gpxyz), 
               1, MAXGAUSS_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_gprst", &(actso3->stress_gprst), 
               1, MAXGAUSS_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_gp123", &(actso3->stress_gp123), 
               1, MAXGAUSS_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_ndxyz", &(actso3->stress_ndxyz), 
               1, MAXNOD_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_ndrst", &(actso3->stress_ndrst), 
               1, MAXNOD_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_nd123", &(actso3->stress_nd123), 
               1, MAXNOD_SOLID3, NUMSTR_SOLID3, 0, "D3");
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
  /* allocate working arrays */
  if (stress == NULL)
  {
    stress = amdef("stress", &(stress_a), NUMSTR_SOLID3, 1,"DV");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_stress_init */


/*======================================================================*/
/*!
\brief Deallocate local globals

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 10/06
*/
void so3_stress_final(PARTITION *actpart)
{
  PARTITION *tpart;
  ELEMENT *actele;
  INT iele, jdis;
  SOLID3 *actso3;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_hlux_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate heat flux arrays */
  tpart = actpart; /*&(actpart[genprob.numtf]); */
  /* loop over all discretisations of partition */
  for (jdis=0; jdis<tpart->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<tpart->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = tpart->pdis[jdis].element[iele];
      /* check wether THERM2 element */
      if (actele->eltyp == el_solid3)
      {
        /* set current SOLID3 element */
        actso3 = actele->e.so3;
        /* deallocate heat flux arrays at element */
        am4del(&(actso3->stress_gpxyz));
        am4del(&(actso3->stress_gprst));
        am4del(&(actso3->stress_gp123));
        am4del(&(actso3->stress_ndxyz));
        am4del(&(actso3->stress_ndrst));
        am4del(&(actso3->stress_nd123));
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
  if (stress)
  {
    amdel(&stress_a);
    stress = NULL;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_stress_final */


/*======================================================================*/
/*!
\brief Evaluate element stresses 

\author bborn
\date 10/06
*/
void so3_stress_cal(CONTAINER *cont,
                    ELEMENT *ele,
                    SO3_DATA *data,
                    MATERIAL *mat)
{
  const INT place = 0;
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* total number of element DOFs */
  DIS_TYP distyp;  /* type of discretisation */
  
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

  DOUBLE det;  /* Jacobi determinants */
  DOUBLE rst[NDIM_SOLID3];  /* natural coordinate a point */
  DOUBLE shape[MAXNOD_SOLID3];
  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3];
  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];
  DOUBLE bop[NUMTMGR_SOLID3][NUMDOF_SOLID3*MAXNOD_SOLID3];
  DOUBLE cmat[NUMHFLX_SOLID3][NUMTMGR_SOLID3];

  INT ihflx;  /* heat flux index */
  INT inode;  /* element nodal index */

  /*--------------------------------------------------------------------*/
  /* start */
  /* Working arrays (locally globals) MUST be initialised! */ 
#ifdef DEBUG
  dstrc_enter("so3_stress_cal");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  nelenod = ele->numnp;
  neledof = NUMDOF_SOLID3 * nelenod;
  distyp = ele->distyp;

  /*--------------------------------------------------------------------*/
  /* number of Gauss points */
  switch (distyp)
  {
    /* hexahedra elements */
    case hex8: case hex20: case hex27:
      gpnumr = ele->e.so3->gpnum[0];
      gpintcr = ele->e.so3->gpintc[0];
      gpnums = ele->e.so3->gpnum[1];
      gpintcs = ele->e.so3->gpintc[1];
      gpnumt = ele->e.so3->gpnum[2];
      gpintct = ele->e.so3->gpintc[2];
      /* total number of Gauss points */
      gpnum = gpnumr * gpnums * gpnumt;
      break;
    /* tetrahedra elements */
    case tet4: case tet10:
      gpnumr = 1;
      gpnums = 1;
      gpnumt = ele->e.so3->gpnum[0];
      gpintcr = 1;
      gpintcs = 1;
      gpintct = ele->e.so3->gpintc[0];
      /* total number of Gauss points */
      gpnum = gpnumt;
      break;
    default:
      dserror("distyp unknown!");
  }  /* end switch */

  /*--------------------------------------------------------------------*/
  /* stress at every Gauss point */
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
        so3_shape_deriv(distyp, gpcr, gpcs, gpct, 1, shape, deriv);
        /*--------------------------------------------------------------*/
        /* compute Jacobian matrix, its determinant and inverse */
        so3_metr_jaco(ele, nelenod, deriv, 1, xjm, &det, xji);
        /*--------------------------------------------------------------*/
        /* factor */
        fac = det;
        /*--------------------------------------------------------------*/
        /* calculate B-operator */
        so3_bop(nelenod, deriv, xji, bop);
        /*--------------------------------------------------------------*/
        /* call material law */
        so3_mat_sel(cont, ele, mat, bop, igp, stress, cmat);
        /*--------------------------------------------------------------*/
        /* store heat flux of current Gauss point */
        for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
        {
          ele->e.so3->stress_gp_xyz.a.d3[place][igp][ihflx] = stress[ihflx];
        }
        /*--------------------------------------------------------------*/
        /* store heat flux in parameter space co-ordinates (r,s,t) */
        so3_stress_rst(xjm, stress, igp, ele->e.so3->stress_gp_rst.a.d3[place]);
        /*--------------------------------------------------------------*/
        /* store modulus and angles at current Gauss point */
        so3_stress_modang(stress, igp, ele->e.so3->stress_gp_123.a.d3[place]);
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
    so3_cfg_noderst(ele, data, inode, rst);
    /*------------------------------------------------------------------*/
    /* extrapolate values now */
    /* heat flux in x-direction : q_1 */
    so3_stress_extrpol(ele, data, gpnum, 
                      ele->e.so3->stress_gp_xyz.a.d3[place], rst,
                      stress);
    for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
    {
      ele->e.so3->stress_nd_xyz.a.d3[place][ihflx][inode] = stress[ihflx];
    }
    /* absolute heat flux and its direction at element nodes */
    so3_stress_modang(stress, inode, ele->e.so3->stress_nd_123.a.d3[place]);
  }  /* end for */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_stress_cal */


/*======================================================================*/
/*!
\brief Convert heat flux with respect to (X,Y,Z) to (r,s,t)

\param **xjm         DOUBLE   (i)   Jacobi matrix at GP
\param  *stress       DOUBLE   (i)   heat flux (at GP)
\param   igp         INT      (i)   index of current Gauss point
\param **stress123    DOUBLE   (o)   modulus and angles at GP
\return void

\author bborn
\date 10/06
*/
void so3_stress_rst(DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3],
                   DOUBLE *stress,
                   INT igp,
                   DOUBLE **stressrst)
{
  INT idim, jdim;
  DOUBLE hfluc;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_stress_rst");
#endif

  /*--------------------------------------------------------------------*/
  dsassert(NDIM_SOLID3 == NUMSTR_SOLID3,
           "Warning dimension of heat flux vector differs to problem "
           "dimension\n");

  /*--------------------------------------------------------------------*/
  /* get locally defined heat flux */
  for (idim=0; idim<NDIM_SOLID3; idim++)
  {
    hfluc = 0.0;
    for (jdim=0; jdim<NDIM_SOLID3; jdim++)
    {
      /* we have to use the 'transposed Jacobian' here,
       * referred to the isoparametric Jacobian notation
       * commonly used, cf. so3_metr.c */
      hfluc = hfluc + xjm[jdim][idim] * stress[jdim];
    }
    stressrst[igp][idim] = hfluc;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void void so3_stress_rst */

/*======================================================================*/
/*!
\brief Modulus of heat flux and its direction measured in angles with
       respect to the global X,Y,Z-axes

\param  *stress       DOUBLE   (i)   heat flux (at GP)
\param   igp         INT      (i)   index of current Gauss point
\param **stress123    DOUBLE   (o)   modulus and angles at GP
\return void

\author bborn
\date 10/06
*/
void so3_stress_modang(DOUBLE *stress,
                      INT igp,
                      DOUBLE **stress123)

{
  INT ihf;  /* count heat flux components */
  DOUBLE stressmod;  /* heat flux modulus */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_stress_modang");
#endif

  /*--------------------------------------------------------------------*/
  /* absolute value of heat flux (modulus) */
  stressmod = 0.0;
  for (ihf=0; ihf<NUMHFLX_SOLID3; ihf++)
  {
    stressmod += stress[ihf]*stress[ihf];
  }
  stressmod = sqrt(stressmod);
  stress123[igp][NUMHFLX_SOLID3] = stressmod;  /* store in ELEMENT array */

  /*--------------------------------------------------------------------*/
  /* direction of heat flux with respect to X,Y,Z-axes */
  if (stressmod != 0.0)
  {
    for (ihf=0; ihf<NUMHFLX_SOLID3; ihf++)
    {
      stress123[igp][ihf] = acos(stress[ihf]/stressmod);
    }
  }
  else
  {
    for (ihf=0; ihf<NUMHFLX_SOLID3; ihf++)
    {
      stress123[igp][ihf] = 0.0;
    }
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_stress_modang */


/*======================================================================*/
/*!
\brief Extrapolate heat fluxes at Gauss points to element nodes

These procedure does not provide more accurate heat fluxes as the direct
computation, but maybe saves a bit of computing time.

\param  *ele      ELEMENT     (i)   pointer to active element
\param  *data     SO3_DATA    (i)   pointer to THERM2 data (GPs coords etc)
\param  ngauss    INT         (i)   total number of Gauss points
\param  **stressgp DOUBLE      (i)   heat flux at Gauss points
\param  *rst      DOUBLE      (i)   element node natural coordinates
\param  *stressnd  DOUBLE      (o)   extrapolated heat flux at node
\return void

\author bborn
\date 03/06
*/
void so3_stress_extrpol(ELEMENT *ele,
                       SO3_DATA *data,
                       INT ngauss,
                       DOUBLE **stressgp,
                       DOUBLE *rst,
                       DOUBLE *stressnd)
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
  DOUBLE funinc[NUMHFLX_SOLID3];  /* increment due to contribution of cur. GP */

  /*--------------------------------------------------------------------*/
  /* start up */
#ifdef DEBUG
  dstrc_enter("so3_stress_extrpol");
#endif

  /*--------------------------------------------------------------------*/
  /* number of Gauss points */
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
      /* total number of Gauss points */
      gpnum = gpnumr * gpnums * gpnumt;
      break;
    /* tetrahedra elements */
    case tet4: case tet10:
      gpnumr = 1;
      gpnums = 1;
      gpnumt = ele->e.so3->gpnum[0];
      gpintcr = 1;
      gpintcs = 1;
      gpintct = ele->e.so3->gpintc[0];
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
  for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
  {
    stressnd[ihflx] = 0.0;
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
             *            th2_stress_cal, as the heat fluxes are stored
             *            in a vector
             */
            for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
            {
              funinc[ihflx] = stressgp[ihflx][igauss];
            }
            /* build tri-directional Lagrange polynomials at Gauss points */
            /* l_{igpr,igps,igpt}(r,s,t) = stress_{igpr,igps,igpt}
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
                for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
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
                for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
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
                for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
                {
                  funinc[ihflx] = funinc[ihflx] 
                    * (rst[2] - gtk)/(gpct - gtk);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* add increment of current Gauss point */
            for (ihflx=0; ihflx<NUMHFLX_SOLID3; ihflx++)
            {
              stressnd[ihflx] = stressnd[ihflx] + funinc[ihflx];
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
}  /* end void so3_stress_extrpol */


/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close)*/
