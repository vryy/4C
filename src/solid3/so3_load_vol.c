/*======================================================================*/
/*!
\file
\brief Spatial integration of loads (ie body forces/traction) applied
       to element domain (volume)

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
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
\date 05/07
*/
extern GENPROB genprob;



/*======================================================================*/
/*!
\brief Spatial integration of 
       body force in element domain (volume) [force/volume]
\param  ele        ELEMENT*           (i)   pointer to current element
\param  gpshade    SO3_GPSHAPEDERIV*  (i)   Gauss point coords etc
\param  ex         DOUBLE[][]         (i)   material coord. of element
\param  timen      DOUBLE             (i)   curr. load factor/curr. time
\param  gvol       GVOL*              (io)  geometry volume of element
\param  eload      DOUBLE[][]         (io)  element load
\author bborn
\date 04/07
*/
void so3_load_vol_int(ELEMENT* ele,
                      SO3_GPSHAPEDERIV* gpshade,
                      DOUBLE ex[MAXNOD_SOLID3][NDIM_SOLID3],
                      const DOUBLE timen,
                      GVOL* gvol,
                      DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  const INT nelenod = ele->numnp;  /* number of element nodes */

  const INT curve = gvol->neum->curve;  /* load curve */
  DOUBLE cfac;  /* curve factor */

  const INT ngp = gpshade->gptot;  /* total number of Gauss points */
  INT jgp;  /* Gauss point index */
  INT fac;  /* integration factor */

  DOUBLE xjm[NDIM_SOLID3][NDIM_SOLID3];  /* Jacobian matrix */
  DOUBLE det;  /* Jacobi determinant */
  DOUBLE xji[NDIM_SOLID3][NDIM_SOLID3];  /* inverse Jacobian matrix */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_vol_int");
#endif

  /*--------------------------------------------------------------------*/
  /* load-step-dependent/time-step-dependent scale */
  if (curve < 0)
  {
    cfac = 1.0;
  }
  else
  {
    if (genprob.timetyp == time_static)
    {
      /* not implemented, could be something like: */
      /* dyn_facfromcurve(curve, timen, &(cfac)); */
      cfac = 1.0;  /* remove this */
    }
    else if (genprob.timetyp == time_dynamic)
    {
      dyn_facfromcurve(curve, timen, &(cfac));
    }
  }

  /*--------------------------------------------------------------------*/
  /* integration loop */
  for (jgp=0; jgp<ngp; jgp++)
  {
    /* initialise intgration factor */
    fac = cfac * gpshade->gpwg[jgp];  /* Curve factor * Gauss weight */
    /* compute Jacobian matrix, its determinant
     * inverse Jacobian is not calculated */
    so3_metr_jaco(ele, nelenod, ex, gpshade->gpderiv[jgp], 0,
                  xjm, &det, xji);
    /* integration (quadrature) factor */
    fac *= det;
    /* volume-load  ==> eload modified */
    so3_load_vol_val(ele, nelenod, gpshade->gpshape[jgp], fac, eload);
  }  /* end of for */
  /*--------------------------------------------------------------------*/
  /* the volume load of this element has been done,
   * -> switch off the volume load condition
   * will be switched on again at beginning of next time step */
  gvol->neum = NULL;

  /*--------------------------------------------------------------------*/
  /* finish */
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief Determine load due to body force on element domain (volume)

\param    ele      ELEMENT*     (i)    actual element
\param    nelenod  INT          (i)    number of element nodes
\param    shape    DOUBLE[]     (i)    shape function at Gauss point
\param    fac      DOUBLE       (i)    integration factor
\param    eload    DOUBLE[][]   (io)   element load vector contribution
\return void

\author bborn
\date 01/07
*/
void so3_load_vol_val(ELEMENT *ele,
                      INT nelenod,
                      DOUBLE shape[MAXNOD_SOLID3],
                      DOUBLE fac,
                      DOUBLE eload[MAXNOD_SOLID3][NUMDOF_SOLID3])
{
  DOUBLE source[NUMDOF_SOLID3];  /* body force [force/vol] */
  INT idof, inode;  /* loopers (i = loaddirection x or y)(j=node) */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_load_vol_val");
#endif
  /*--------------------------------------------------------------------*/
  /* distinguish load type */
  switch(ele->g.gvol->neum->neum_type)
  {
    /*------------------------------------------------------------------*/
    /* uniform (density-proportional) dead load */
    case neum_dead:
      for (idof=0; idof<NUMDOF_SOLID3; idof++)
      {
        source[idof] = ele->g.gvol->neum->neum_val.a.dv[idof];
      }
      /* add load vector component to element load vector */
      for (inode=0; inode<nelenod; inode++)
      {
        for (idof=0; idof<NUMDOF_SOLID3; idof++)
        {
          eload[inode][idof] += shape[inode] * source[idof] * fac;
        }
      }
      break;
    /*------------------------------------------------------------------*/
    default:
      dserror("load case unknown");
      break;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of so3_load_vol_val(...) */

/*======================================================================*/
#endif  /*end of #ifdef D_SOLID3 */
