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
               1, MAXGAUSS_SOLID3, 4*NDIM_SOLID3, 0, "D3");
        am4def("stress_ndxyz", &(actso3->stress_ndxyz), 
               1, MAXNOD_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_ndrst", &(actso3->stress_ndrst), 
               1, MAXNOD_SOLID3, NUMSTR_SOLID3, 0, "D3");
        am4def("stress_nd123", &(actso3->stress_nd123), 
               1, MAXNOD_SOLID3, 4*NDIM_SOLID3, 0, "D3");
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

\param  cont      CONTAINER*         (i)    see container.h
\param  ele       ELE*               (io)   current element
\param  data      SO3_DATA*          (i)    common SOLID3 data
\param  gpshade   SO3_GPSHAPEDERIV*  (i)    Gauss point coords
                                            and shape functions

\param  mat       MATERIAL*          (i)    
\return void

\author bborn
\date 01/07
*/
void so3_stress(CONTAINER *cont,
                ELEMENT *ele,
                SO3_DATA *data,
                SO3_GPSHAPEDERIV *gpshade,
                INT imyrank,
                MATERIAL *mat)
{
  const INT place = 0;
  INT nelenod;  /* number of element nodes */
  INT neledof;  /* total number of element DOFs */
  DIS_TYP distyp;  /* type of discretisation */
  DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3];
  DOUBLE edis[MAXNOD_SOLID3][NDIM_SOLID3];
  INT gpnum = 0;  /* total number of Gauss points in element domain */

  INT igp;  /* total index Gauss point */
  INT inod;  /* element nodal index */
  INT jdim;  /* dimension index */
  INT istr;  /* stress/strain index */
  NODE *actnode;  /* current node */

  DOUBLE fac;  /* a factor */
  DOUBLE gpcr = 0.0;  /* GP r-coord */
  DOUBLE gpcs = 0.0;  /* GP s-coord */
  DOUBLE gpct = 0.0;  /* GP t-coord */
  DOUBLE rst[NDIM_SOLID3];  /* natural coordinate a point */
  SO3_GEODEFSTR gds;  /* isoparametric Jacobian, deformation grad, etc at
                       * Gauss point */
  DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3];
  DOUBLE stress[NUMSTR_SOLID3];
  DOUBLE bopn[MAXDOF_SOLID3][NUMDOF_SOLID3];
  DOUBLE bop[NUMSTR_SOLID3][MAXDOF_SOLID3];



  /*--------------------------------------------------------------------*/
  /* start */
  /* Working arrays (locally globals) MUST be initialised! */ 
#ifdef DEBUG
  dstrc_enter("so3_stress");
#endif

  /*--------------------------------------------------------------------*/
  /* element properties */
  distyp = ele->distyp;
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
  gpnum = gpshade->gptot;  /* total number of Gauss points in domain */


  /*--------------------------------------------------------------------*/
  /* stress at every Gauss point */
  /* parse all Gauss points and evaluate stress */
  for (igp=0; igp<gpnum; igp++)
  {
    /*------------------------------------------------------------------*/
    /* Gauss point */
    gpcr = gpshade->gpco[igp][0];  /* r-coordinate */
    gpcs = gpshade->gpco[igp][1];  /* s-coordinate */
    gpct = gpshade->gpco[igp][2];  /* t-coordinate */
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
    /* calculate B-operator */
    so3_bop(ele, nelenod, gpshade->gpderiv[igp], gds.xji, gds.defgrd,
            bopn, bop);
    /*------------------------------------------------------------------*/
    /* call material law ==> 2nd PK-stresses and constitutive matrix */
    so3_mat_sel(ele, mat, igp, &gds, stress, cmat);
    /*------------------------------------------------------------------*/
    /* store stress of current Gauss point */
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      ele->e.so3->stress_gpxyz.a.d3[place][igp][istr] = stress[istr];
    }
    /*------------------------------------------------------------------*/
    /* construct rotational component of inverse FE-Jacobian */
    so3_metr_rot(gds.xjm, gds.xrm, gds.xrvm, gds.xrvi);
    /*------------------------------------------------------------------*/
    /* store heat flux in parameter space co-ordinates (r,s,t) */
    so3_stress_rst(gds.xrvi, stress, 
                   ele->e.so3->stress_gprst.a.d3[place][igp]);
    /*------------------------------------------------------------------*/
    /* store principle and direction angles at current Gauss point */
    so3_stress_123(stress, ele->e.so3->stress_gp123.a.d3[place][igp]);
  }

  /*--------------------------------------------------------------------*/
  /* Stresses at Gauss points are extrapolated to element nodes. */
  /* loop all element nodes */
  for (inod=0; inod<nelenod; inod++)
  {
    /*------------------------------------------------------------------*/
    /* get local coordinates of node ==> rst */
    so3_cfg_noderst(ele, data, inod, rst);
    /*------------------------------------------------------------------*/
    /* extrapolate values now */
    so3_stress_extrpol(ele, data, gpnum, 
                       ele->e.so3->stress_gpxyz.a.d3[place], rst,
                       stress);
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      ele->e.so3->stress_ndxyz.a.d3[place][inod][istr] = stress[istr];
    }
    /* store principle and direction angles at current Gauss point */
    so3_stress_123(stress, ele->e.so3->stress_nd123.a.d3[place][inod]);
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_stress_cal */


/*======================================================================*/
/*!
\brief Redirect stress vector with respect to (X,Y,Z) to (r,s,t)

\param   xrvi        DOUBLE[][]  (i)   spec. rotation matrix (rst)->(XYZ)
                                       for symmtric 2-tensors stored
                                       vectorially
\param   stress      DOUBLE[]    (i)   stress (at GP)
\param   stressrst   DOUBLE**    (o)   modulus and angles at GP
\return void

\author bborn
\date 01/07
*/
void so3_stress_rst(DOUBLE xrvi[NUMSTR_SOLID3][NUMSTR_SOLID3],
                    DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE *stressrst)
{
  INT istr, jstr;  /* indices */
  DOUBLE sc;  /* intermediate stress component */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_stress_rst");
#endif

  /*--------------------------------------------------------------------*/
  /* get locally oriented stress vector */
  for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
  {
    sc = 0.0;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      sc += xrvi[jstr][istr] * stress[istr];
    }
    stressrst[jstr] = sc;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void void so3_stress_rst */

/*======================================================================*/
/*!
\brief Principial stresses and principial directions at current
       Gauss point

\param  *stress       DOUBLE   (i)   stress vector
\param  *stress123    DOUBLE   (o)   principial stresses
\return void

\author bborn
\date 01/07
*/
void so3_stress_123(DOUBLE stress[NUMSTR_SOLID3],
                    DOUBLE *stress123)

{
  DOUBLE stresst[NDIM_SOLID3][NDIM_SOLID3];  /* stress tensor */
  DOUBLE ew[NDIM_SOLID3];  /* principial stresses */
  INT err;  /* error 0=PASSED, 1=ERROR */

#ifdef DEBUG
  dstrc_enter("so3_stress_123");
#endif

  /* stress tensor */
  so3_tns3_v2tsym(stress, stresst);

  /* principial stresses */
  so3_tns3_spcdcmp(stresst, &err, ew);

  /* write principial stresses -- if found */
  if (err == 0)
  {
    stress123[0] = ew[0];  /* 1st principal stress */
    stress123[1] = ew[1];  /* 2nd principal stress */
    stress123[2] = ew[2];  /* 3rd principal stress */
    stress123[3] = -1.0;  /* 1st stress 1st direction */
    stress123[4] = -1.0;  /* 1st stress 2nd direction */
    stress123[5] = -1.0;  /* 1st stress 3rd direction */
    stress123[6] = -1.0;  /* 2st stress 1st direction */
    stress123[7] = -1.0;  /* 2st stress 2nd direction */
    stress123[8] = -1.0;  /* 2st stress 3rd direction */
    stress123[9] = -1.0;  /* 3rd stress 1st direction */
    stress123[10] = -1.0;  /* 3rd stress 2nd direction */
    stress123[11] = -1.0;  /* 3rd stress 3rd direction */
  }
  else
  {
    stress123[0] = -1.0;  /* 1st principal stress */
    stress123[1] = -1.0;  /* 2nd principal stress */
    stress123[2] = -1.0;  /* 3rd principal stress */
    stress123[3] = -1.0;  /* 1st stress 1st direction */
    stress123[4] = -1.0;  /* 1st stress 2nd direction */
    stress123[5] = -1.0;  /* 1st stress 3rd direction */
    stress123[6] = -1.0;  /* 2st stress 1st direction */
    stress123[7] = -1.0;  /* 2st stress 2nd direction */
    stress123[8] = -1.0;  /* 2st stress 3rd direction */
    stress123[9] = -1.0;  /* 3rd stress 1st direction */
    stress123[10] = -1.0;  /* 3rd stress 2nd direction */
    stress123[11] = -1.0;  /* 3rd stress 3rd direction */
  }

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

\param  *ele       ELEMENT     (i)   pointer to active element
\param  *data      SO3_DATA    (i)   pointer to THERM2 data (GPs coords etc)
\param  ngauss     INT         (i)   total number of Gauss points
\param  **stressgp DOUBLE      (i)   heat flux at Gauss points
\param  rst[]      DOUBLE      (i)   element node natural coordinates
\param  stressnd[] DOUBLE      (o)   extrapolated heat flux at node
\return void

\author bborn
\date 03/06
*/
void so3_stress_extrpol(ELEMENT *ele,
                        SO3_DATA *data,
                        INT ngauss,
                        DOUBLE **stressgp,
                        DOUBLE rst[NDIM_SOLID3],
                        DOUBLE stressnd[NUMSTR_SOLID3])
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

  INT istr;  /* stress index */

  INT i, j, k;  /* loop counter */
  DOUBLE gri, gsj, gtk;  /* Lagrange stages */
  DOUBLE funinc[NUMSTR_SOLID3];  /* increment due to contribution of cur. GP */

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
  for (istr=0; istr<NUMSTR_SOLID3; istr++)
  {
    stressnd[istr] = 0.0;
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
            for (istr=0; istr<NUMSTR_SOLID3; istr++)
            {
              funinc[istr] = stressgp[igauss][istr];
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
                for (istr=0; istr<NUMSTR_SOLID3; istr++)
                {
                  funinc[istr] *= (rst[0] - gri)/(gpcr - gri);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* in s-direction : l_{igps}^{gpintcs}(s) */
            for (j=0; j<gpnums; j++)
            {
              if (j != igps)
              {
                gsj = data->ghlc[gpintcs][j];
                for (istr=0; istr<NUMSTR_SOLID3; istr++)
                {
                  funinc[istr] *= (rst[1] - gsj)/(gpcs - gsj);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* in t-direction : l_{igpt}^{gpintct}(t) */
            for (k=0; k<gpnumt; k++)
            {
              if (k != igpt)
              {
                gtk = data->ghlc[gpintct][k];
                for (istr=0; istr<NUMSTR_SOLID3; istr++)
                {
                  funinc[istr] *= (rst[2] - gtk)/(gpct - gtk);
                }  /* end for */
              }  /* end if */
            }  /* end for */
            /* add increment of current Gauss point */
            for (istr=0; istr<NUMSTR_SOLID3; istr++)
            {
              stressnd[istr] += funinc[istr];
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
