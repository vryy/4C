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

/*!
\addtogroup SOLID3
*//*! @{ (documentation module open)*/


/*======================================================================*/
/*!
\brief Select proper material law

\param *ele       ELEMENT        (i)   pointer to current element
\param *mat       MATERIAL       (i)   pointer to current material
\param ip         INT            (i)   current Gauss point index
\param *gds       SO3_GEODEFSTR  (i)   geom. & def. data at Gauss point
\param stress[]   DOUBLE         (o)   linear(Biot)/2.Piola-Kirchhoff stress
\param cmat[][]   DOUBLE         (o)   constitutive matrix
\return void

\author mf
\date 10/06
*/
void so3_mat_sel(ELEMENT *ele,
                 MATERIAL *mat,
                 INT ip,
                 SO3_GEODEFSTR *gds,
                 DOUBLE stress[NUMSTR_SOLID3],
                 DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3])
{
  INT imat, jmat, istrn, istss, inode;  /* counters */
  DOUBLE Emod, nu, mfac, strainsum, stresssum;
  DOUBLE strain[NUMSTR_SOLID3];  /* strain vector */
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
    case struct_stvenant:
      /*----------------------------------------------------------------*/
      /* isotropic elasticity tensor C in matrix notion */
      /*                       [ 1-nu     nu     nu |          0    0    0 ]
       *                       [        1-nu     nu |          0    0    0 ]
       *           E           [               1-nu |          0    0    0 ]
       *   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
       *       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
       *                       [                    |      (1-2*nu)/2    0 ]
       *                       [ symmetric          |           (1-2*nu)/2 ]
       */
      Emod = mat->youngs;  /* Young's modulus (modulus of elasticity */
      nu = mat->possionratio;  /* Poisson's ratio */
      mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
      /* constitutive matrix */
      /* set the whole thing to zero */
      memset(cmat, 0, sizeof(cmat));
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
      /*----------------------------------------------------------------*/
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
      /*----------------------------------------------------------------*/
      /* compute stress vector */
      for (istss=0; istss<NUMSTR_SOLID3; istss++)
      {
        stresssum = 0.0;
        for (istrn=0; istrn<NUMSTR_SOLID3; istrn++)
        {
          stresssum += cmat[istss][istrn] * strain[istrn];
        }
        stress[istress] = stresssum;
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
#endif  /* end of #ifdef D_SOLID3 */
/*! @} (documentation module close) */
