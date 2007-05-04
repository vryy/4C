/*======================================================================*/
/*!
\file
\brief De St. Venant--Kirchhoff's linear elastic material
       (linear elastic in a material sense)
       (or Hooke's material in geometrically linear set-up)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/bornemann
	    089-289-15237
</pre>

\author bborn
\date 04/07
*/
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/07
*/
extern GENPROB genprob;

/*======================================================================*/
/*!
\brief Apply St.Venant--Kirchhoff's material
\author bborn
\date 04/07
*/
void so3_mat_stvenant_sel(CONTAINER* container,
                          ELEMENT* ele,
                          MATERIAL* mat,
                          INT ip,
                          SO3_GEODEFSTR* gds,
                          DOUBLE stress[NUMSTR_SOLID3],
                          DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_stvenant_sel");
#endif

  /* Young's modulus (modulus of elasticity */
  DOUBLE Emod = mat->m.stvenant->youngs;
  /* Poisson's ratio */
  DOUBLE nu = mat->m.stvenant->possionratio;
  /* INT istrn, istss, istr, jstr; */
  DOUBLE strain[NUMSTR_SOLID3];  /* strain vector */

  /*--------------------------------------------------------------------*/
  /* isotropic elasticity tensor C in matrix notion */
  /*                       [ 1-nu     nu     nu |          0    0    0 ]
   *                       [        1-nu     nu |          0    0    0 ]
   *           E           [               1-nu |          0    0    0 ]
   *   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
   *       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
   *                       [                    |      (1-2*nu)/2    0 ]
   *                       [ symmetric          |           (1-2*nu)/2 ]
   */
#ifdef D_MAT
  {
    /* We need to create a dynamically allocated array to hold
     * elasticity matrix. The reason is due to C being incapable of
     * treating statically and dynamically allocated multi-dimensional
     * arrays. In short: mat_el_iso(...) requires a DOUBLE**, whereas
     * our cmat is DOUBLE[][]. According to, e.g., Section 2.5 in
     * http://www.lysator.liu.se/c/c-faq/c-2.html
     * these types cannot be reliably recasted in each other. */
    /* make dynamically allocated array */
    ARRAY ccmat_a;
    DOUBLE **ccmat;
    ccmat = amdef("elasticitymatrix", &ccmat_a,
                  NUMSTR_SOLID3, NUMSTR_SOLID3, "DA");
    /* retrieve elasticity matrix */
    mat_el_iso(Emod, nu, ccmat);
    /* copy arrays */
    INT istr, jstr;
    for (istr=0; istr<NUMSTR_SOLID3; istr++)
    {
      for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
      {
        cmat[istr][jstr] = ccmat[istr][jstr];
      }
    }
    /* deallocate */
    amdel(&ccmat_a);
  }
#else
  {
    DOUBLE mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
    /* constitutive matrix */
    /* set the whole thing to zero */
    memset(cmat, 0, NUMSTR_SOLID3*NUMSTR_SOLID3*sizeof(DOUBLE));
/*     INT istr, jstr; */
/*     for (istr=0; istr<NUMSTR_SOLID3; istr++) */
/*     { */
/*       for (jstr=0; jstr<NUMSTR_SOLID3; jstr++) */
/*       { */
/*         cmat[istr][jstr] = 0.0; */
/*       } */
/*     } */
    /* write non-zero components */
    cmat[0][0] = mfac*(1.0-nu);
    cmat[0][1] = mfac*nu;
    cmat[0][2] = mfac*nu;
    cmat[1][0] = mfac*nu;
    cmat[1][1] = mfac*(1.0-nu);
    cmat[1][2] = mfac*nu;
    cmat[2][0] = mfac*nu;
    cmat[2][1] = mfac*nu;
    cmat[2][2] = mfac*(1.0-nu);
    /* ~~~ */
    cmat[3][3] = mfac*0.5*(1.0-2.0*nu);
    cmat[4][4] = mfac*0.5*(1.0-2.0*nu);
    cmat[5][5] = mfac*0.5*(1.0-2.0*nu);
  }
#endif
  /*--------------------------------------------------------------*/
  /* set local strain vector */
  if (ele->e.so3->kintype == so3_geo_lin)
  {
    /* linear (engineering) strain vector */
    INT istrn;
    for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
    {
      strain[istrn] = gds->stnengv[istrn];
    }
  } 
  else if (ele->e.so3->kintype == so3_total_lagr)
  {
    /* Green-Lagrange strain vector */
    INT istrn;
    for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
    {
      strain[istrn] = gds->stnglv[istrn];
    }
  }
  else
  {
    dserror("Cannot digest chosen type of spatial kinematic\n");
  }
  /*--------------------------------------------------------------*/
  /* strain due to thermal expansion */
#ifdef D_TSI
  if (genprob.probtyp == prb_tsi)
  {
    DOUBLE tem;  /* temperature */
    /* temperature at Gauss point */
    so3_tsi_temper(container, ele,
                   gds->gpc[0], gds->gpc[1], gds->gpc[2], 
                   &tem);
    /* coefficient of linear thermal expansion */
    DOUBLE thermexpans = mat->m.stvenant->thermexpans;

    /* thermal strain vector */
    strain[0] -= thermexpans * tem;  /* E_xx */
    strain[1] -= thermexpans * tem;  /* E_yy */
    strain[2] -= thermexpans * tem;  /* E_zz */
    /* do nothing with E_xy, E_yz and E_zx */
  }
#endif
  /*--------------------------------------------------------------*/
  /* compute stress vector */
  INT istss;
  for (istss=0; istss<NUMSTR_SOLID3; istss++)
  {
    INT istrn;
    DOUBLE stresssum = 0.0;  /* intermediate row * column sum */
    for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
    {
      stresssum += cmat[istss][istrn] * strain[istrn];
    }
    stress[istss] = stresssum;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

#endif
