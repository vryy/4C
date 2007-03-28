/*======================================================================*/
/*!
\file
\brief Select proper material law

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
#ifdef D_MAT
#include "../materials/mat_prototypes.h"
#endif

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*/
/*!
\brief General problem data

\author bborn
\date 03/07
*/
extern GENPROB genprob;


/*======================================================================*/
/*!
\brief Select proper material law

\param *container CONTAINER      (i)   container
\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void

\author bborn
\date 03/07
*/
void so3_mat_sel(CONTAINER *container,
                 ELEMENT *ele,
                 MATERIAL *mat,
                 INT ip,
                 SO3_GEODEFSTR *gds,
                 DOUBLE stress[NUMSTR_SOLID3],
                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_mat_sel");
#endif

  /*====================================================================*/
  /* ===> central material routines
   *
   *      These materials are supposed to be connected to the 
   *      existant (or new?) central material routines.
   *      Right now, only the simple St.Venant-Kirchhoff material
   *      is included to test the element.
   *
   * ===> central material routines
   */

  /*====================================================================*/
  /* the material law (it's a material world!) */
  switch (mat->mattyp)
  {
    /*------------------------------------------------------------------*/
    /* St. Venant-Kirchhoff material */
    case m_stvenant:
      {
        /* Young's modulus (modulus of elasticity */
        DOUBLE Emod = mat->m.stvenant->youngs;
        /* Poisson's ratio */
        DOUBLE nu = mat->m.stvenant->possionratio;
        INT istrn, istss, istr, jstr;
        DOUBLE strain[NUMSTR_SOLID3];  /* strain vector */
        /*--------------------------------------------------------------*/
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
          /* We need to create a dynmically allocated array to hold
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
          for (istr=0; istr<NUMSTR_SOLID3; istr++)
          {
            for (jstr=0; jstr<NUMSTR_SOLID3; jstr++)
            {
              cmat[istr][jstr] = 0.0;
            }
          }
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
          for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
          {
            strain[istrn] = gds->stnengv[istrn];
          }
        } 
        else if (ele->e.so3->kintype == so3_total_lagr)
        {
          /* Green-Lagrange strain vector */
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
#ifdef D_THERM3
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
#endif
        /*--------------------------------------------------------------*/
        /* compute stress vector */
        for (istss=0; istss<NUMSTR_SOLID3; istss++)
        {
          DOUBLE stresssum = 0.0;  /* intermediate row * column sum */
          for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
          {
            stresssum += cmat[istss][istrn] * strain[istrn];
          }
          stress[istss] = stresssum;
        }
      }
      break;
    default:
      dserror("Type of material law is not applicable");
      break;
  }  /* end of switch (mat->mattyp) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of so3_mat_sel */

/*======================================================================*/
/*!
\brief get density out of material law

\param  mat       MATERIAL*   (i)   material data
\param  density   DOUBLE*     (o)   density value
\return void

\author bborn
\date 01/07
*/
void so3_mat_density(MATERIAL *mat, 
                     DOUBLE *density)
{

#ifdef DEBUG
  dstrc_enter("so3_mat_density");
#endif

  /* switch material type */
  switch(mat->mattyp)
  {
    /* ST.VENANT-KIRCHHOFF-MATERIAL */
    case m_stvenant:
      *density = mat->m.stvenant->density;
      break;
    /* kompressible neo-hooke */
    case m_neohooke:
      *density = mat->m.neohooke->density;
      break;
    /* porous linear elastic */
    case m_stvenpor:
      *density = mat->m.stvenpor->density;
      break;
    /* hyperelastic polyconvex material */
    case m_hyper_polyconvex:
      *density = mat ->m.hyper_polyconvex->density;
      break;
    /* */
    default:
      dserror("Density of chosen material is not defined!");
      break;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
