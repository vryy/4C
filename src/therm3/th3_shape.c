/*======================================================================*/
/*!
\file
\brief Calculate shape functions and their derivatives with respect to
       the parameter space

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 09/06
*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"

/*!
\addtogroup THERM3
*//*! @{ (documentation module open)*/

/*======================================================================*/
/*!
\brief Shape functions and their natural derivatives at point (r,s,t)

\param  typ        DIS_TYP   (i)    discretisation type
\param  r          DOUBLE    (i)    r-coord of (r,s,t)
\param  s          DOUBLE    (i)    s-coord of (r,s,t)
\param  t          DOUBLE    (i)    t-coord of (r,s,t)
\param  option     INT       (i)    option
                                    ==0 : only shape functions
                                    ==1 : shape functions and derivatives
\param  *shape     DOUBLE    (o)    shape function value at (r,s)
\param  **deriv    DOUBLE    (o)    shape function derivative at (r,s)

\return void

\author bborn
\date 09/06
*/
void th3_shape_deriv(DIS_TYP     typ,
                     DOUBLE      r,
                     DOUBLE      s,
                     DOUBLE      t,
                     INT         option,
                     DOUBLE      shape[MAXNOD_THERM3],
                     DOUBLE      deriv[MAXNOD_THERM3][NDIM_THERM3])
{
  DOUBLE rp, sp, tp, rm, sm, tm, rrm, ssm, ttm;
  DOUBLE u;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_shape_deriv");
#endif

  /*--------------------------------------------------------------------*/
  /* a few parameters depending on coord of current point (r,s) */
  rm  = 1.0 - r;
  rp  = 1.0 + r;
  sm  = 1.0 - s;
  sp  = 1.0 + s;
  tm  = 1.0 - t;
  tp  = 1.0 + t;
  rrm = 1.0 - r*r;
  ssm = 1.0 - s*s;
  ttm = 1.0 - t*t;
  /*--------------------------------------------------------------------*/
  /* switch according to discretisation type */
  switch (typ)
  {
    /*==================================================================*/
    /* hexahedral elements */
    /* linear interpolation */
    case hex8:
      shape[0] = 0.125*rm*sm*tm;  /* shape fct. assoc. to node 0 */
      shape[1] = 0.125*rp*sm*tm;
      shape[2] = 0.125*rp*sp*tm;
      shape[3] = 0.125*rm*sp*tm;
      shape[4] = 0.125*rm*sm*tp;
      shape[5] = 0.125*rp*sm*tp;
      shape[6] = 0.125*rp*sp*tp;
      shape[7] = 0.125*rm*sp*tp;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* shape fct. assoc. to node 0 */
        deriv[0][0] = -0.125*sm*tm;  /* differentiated in r */
        deriv[0][1] = -0.125*rm*tm;  /* differentiated in s */
        deriv[0][2] = -0.125*rm*sm;  /* differentiated in t */
        /* shape fct. assoc. to node 1 */
        deriv[1][0] = 0.125*sm*tm;
        deriv[1][1] = -0.125*rp*tm;
        deriv[1][2] = -0.125*rp*sm;
        /* shape fct. assoc. to node 2 */
        deriv[2][0] = 0.125*sp*tm;
        deriv[2][1] = 0.125*rp*tm;
        deriv[2][2] = -0.125*rp*sp;
        /* shape fct. assoc. to node 3 */
        deriv[3][0] = -0.125*sp*tm;
        deriv[3][1] = 0.125*rm*tm;
        deriv[3][2] = -0.125*rm*sp;
        /* shape fct. assoc. to node 4 */
        deriv[4][0] = -0.125*sm*tp;
        deriv[4][1] = -0.125*rm*tp;
        deriv[4][2] = 0.125*rm*sm;
        /* shape fct. assoc. to node 5 */
        deriv[5][0] = 0.125*sm*tp;
        deriv[5][1] = -0.125*rp*tp;
        deriv[5][2] = 0.125*rp*sp;
        /* shape fct. assoc. to node 6 */
        deriv[6][0] = 0.125*sp*tp;
        deriv[6][1] = 0.125*rp*tp;
        deriv[6][2] = 0.125*rp*sp;
        /* shape fct. assoc. to node 7 */
        deriv[7][0] = -0.125*sp*tp;
        deriv[7][1] = 0.125*rm*tp;
        deriv[7][2] = 0.125*rm*sp;
      }
      break;
    /* quadratic interpolation without central nodes (serendipity) */
    case hex20:
      /* shape functions associated to vertex nodes k=1,...,8
       * N^k = 1/8 (1 + r^k r) (1 + s^k s) (1 + t^k k) 
       *           (r^k r + s^k s + t^k t - 2)
       * with r^k,s^k,t^k = -1,+1
       * [Zienkiewicz, Methode der Finiten Elemente, Hanser, 1975]
       * However, here the slightly different notation is used
       * N^k = 1/8 (1 + r^k r) (1 + s^k s) (1 + t^k k) 
       *           ( (1 + r^k r) + (1 + s^k s) + (1 + t^k t) - 2 - 3) */
      shape[0] = 0.125*rm*sm*tm*(rm+sm+tm-5.0);  
      shape[1] = 0.125*rp*sm*tm*(rp+sm+tm-5.0);
      shape[2] = 0.125*rp*sp*tm*(rp+sp+tm-5.0);
      shape[3] = 0.125*rm*sp*tm*(rm+sp+tm-5.0);
      shape[4] = 0.125*rm*sm*tp*(rm+sm+tp-5.0);
      shape[5] = 0.125*rp*sm*tp*(rp+sm+tp-5.0);
      shape[6] = 0.125*rp*sp*tp*(rp+sp+tp-5.0);
      shape[7] = 0.125*rm*sp*tp*(rm+sp+tp-5.0);
      /* shape functions associated to middle nodes on edges k=8,...,19
       * N^k = 1/4 (1 - r r) (1 + s^k s) (1 + t^k t)
       * with r^k=0, s^k,t^k = -1,+1
       * analogously for s^k,t^k=0 
       * [Zienkiewicz, Methode der Finiten Elemente, Hanser, 1975] */
      shape[8] = 0.25*rrm*sm*tm;
      shape[9] = 0.25*rp*ssm*tm;
      shape[10] = 0.25*rrm*sp*tm;
      shape[11] = 0.25*rm*ssm*tm;
      shape[16] = 0.25*rrm*sm*tp;
      shape[17] = 0.25*rp*ssm*tp;
      shape[18] = 0.25*rrm*sp*tp;
      shape[19] = 0.25*rm*ssm*tp;
      shape[12] = 0.25*rm*sm*ttm;
      shape[13] = 0.25*rp*sm*ttm;
      shape[14] = 0.25*rp*sp*ttm;
      shape[15] = 0.25*rm*sp*ttm;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* corners */
        deriv[0][0] = -0.125*sm*rm*(2.0*rm+sm+tm-5.0); 
        deriv[0][1] = -0.125*rm*tm*(rm+2.0*sm+tm-5.0); 
        deriv[0][2] = -0.125*rm*sm*(rm+sm+2.0*tm-5.0);
        deriv[1][0] = 0.125*sm*tm*(2.0*rp+sm+tm-5.0);
        deriv[1][1] = -0.125*rp*tm*(rp+2.0*sm+tm-5.0); 
        deriv[1][2] = -0.125*rp*sm*(rp+sm+2.0*tm-5.0);
        deriv[2][0] = 0.125*sp*tm*(2.0*rp+sp+tm-5.0);
        deriv[2][1] = 0.125*rp*tm*(rp+2.0*sp+tm-5.0);
        deriv[2][2] = -0.125*rp*sp*(rp+sp+2.0*tm-5.0);
        deriv[3][0] = -0.125*sp*tm*(2.0*rm+sp+tm-5.0);
        deriv[3][1] = 0.125*rm*tm*(rm+2.0*sp+tm-5.0);
        deriv[3][2] = -0.125*rm*sp*(rm+sp+2.0*tm-5.0);
        deriv[4][0] = -0.125*sm*tp*(2.0*rm+sm+tp-5.0);
        deriv[4][1] = -0.125*rm*tp*(rm+2.0*sm+tp-5.0);
        deriv[4][2] = 0.125*rm*sm*(rm+sm+2.0*tp-5.0);
        deriv[5][0] = 0.125*sm*tp*(2.0*rp+sm+tp-5.0);
        deriv[5][1] = -0.125*rp*tp*(rp+2.0*sm+tp-5.0);
        deriv[5][2] = 0.125*rp*sm*(rp+sm+2.0*tp-5.0);
        deriv[6][0] = 0.125*sp*tp*(2.0*rp+sp+tp-5.0);
        deriv[6][1] = 0.125*rp*tp*(rp+2.0*sp+tp-5.0);
        deriv[6][2] = 0.125*rp*sp*(rp+sp+2.0*tp-5.0);
        deriv[7][0] = -0.125*sp*tp*(2.0*rm+sp+tp-5.0);
        deriv[7][1] = 0.125*rm*tp*(rm+2.0*sp+tp-5.0);
        deriv[7][2] = 0.125*rm*sp*(rm+sp+2.0*tp-5.0);
        /* centres in (t=-1)-plane */
        deriv[8][0] = 0.5*(rm-1.0)*sm*tm;
        deriv[8][1] = -0.25*rrm*tm;
        deriv[8][2] = -0.25*rrm*sm;
        deriv[9][0] = 0.25*ssm*tm;
        deriv[9][1] = 0.5*rp*(sm-1.0)*tm;
        deriv[9][2] = -0.25*rp*ssm;
        deriv[10][0] = 0.5*(rm-1.0)*sp*tm;
        deriv[10][1] = 0.25*rrm*tm;
        deriv[10][2] = -0.25*rrm*sp;
        deriv[11][0] = -0.25*ssm*tm;
        deriv[11][1] = 0.5*rm*(sm-1.0)*tm;
        deriv[11][2] = -0.25*rm*ssm;
        /* centres in (t=1)-plane */
        deriv[16][0] = 0.5*(rm-1.0)*sm*tp;
        deriv[16][1] = -0.25*rrm*tp;
        deriv[16][2] = 0.25*rrm*sm;
        deriv[17][0] = 0.25*ssm*tp;
        deriv[17][1] = 0.5*rp*(sm-1.0)*tp;
        deriv[17][2] = 0.25*rp*ssm;
        deriv[18][0] = 0.5*(rm-1.0)*sp*tp;
        deriv[18][1] = 0.25*rrm*tp;
        deriv[18][2] = 0.25*rrm*sp;
        deriv[19][0] = -0.25*ssm*tp;
        deriv[19][1] = 0.5*rm*(sm-1.0)*tp;
        deriv[19][2] = 0.25*rm*ssm;
        /* centres in (t=0)-plane */
        deriv[12][0] = -0.25*sm*ttm;
        deriv[12][1] = -0.25*rm*ttm;
        deriv[12][2] = 0.5*rm*sm*(tm-1.0);
        deriv[13][0] = 0.25*sm*ttm;
        deriv[13][1] = -0.25*rp*ttm;
        deriv[13][2] = 0.5*rp*sm*(tm-1.0);
        deriv[14][0] = 0.25*sp*ttm;
        deriv[14][1] = 0.25*rp*ttm;
        deriv[14][2] = 0.5*rp*sp*(tm-1.0);
        deriv[15][0] = -0.25*sp*ttm;
        deriv[15][1] = 0.25*rm*ttm;
        deriv[15][2] = 0.5*rm*sp*(tm-1.0);
      }
      break;
    /* quadratic interpolation with central node (Lagrangian) */
    case hex27:
      dserror("hex27 shape functions are not implemented");
      /* optionally include derivatives */
      if (option == 1)
      {
        dserror("Derivatives of hex27 shape functions are not implemented");
      }
      break;
    /*==================================================================*/
    /* tetrahedral elements */
    /* linear interpolation */
    case tet4:
      shape[0] = r;
      shape[1] = s;
      shape[2] = t;
      shape[3] = 1.0 - r - s - t;
      /* optionally include derivatives */
      if (option == 1)
      {
        deriv[0][0] = 1.0;
        deriv[0][1] = 0.0;
        deriv[0][2] = 0.0;
        deriv[1][0] = 0.0;
        deriv[1][1] = 1.0;
        deriv[1][2] = 0.0;
        deriv[2][0] = 0.0;
        deriv[2][1] = 0.0;
        deriv[2][2] = 1.0;
        deriv[3][0] = -1.0;
        deriv[3][1] = -1.0;
        deriv[3][2] = -1.0;
      }
      break;
    /* Quadratic interpolation */
    case tet10:
      u = 1.0-r-s-t;
      dserror("tet10 shape functions are not implemented");
      /* optionally include derivatives */
      if (option == 1)
      {
        dserror("Derivatives of tet10 shape functions are not implemented");
      }
      break;
    /* catch erroneous discretisation types */
    default:
      dserror("unknown typ of interpolation");
      break;
  } /* end of switch typ */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_shape_deriv */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */



