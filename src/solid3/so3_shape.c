/*======================================================================*/
/*!
\file
\brief Calculate shape functions and their derivatives with respect to
       the parameter space

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089-289-15240
</pre>

\author mf
\date 10/06
*/
#ifndef CCADISCRET
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
\brief Shape functions and their natural derivatives at point (r,s,t)

\param  typ        DIS_TYP   (i)    discretisation type
\param  r          DOUBLE    (i)    r-coord of (r,s,t)
\param  s          DOUBLE    (i)    s-coord of (r,s,t)
\param  t          DOUBLE    (i)    t-coord of (r,s,t)
\param  option     INT       (i)    option
                                    ==0 : only shape functions
                                    ==1 : shape functions and derivatives
\param  shape[]    DOUBLE    (o)    shape function value at (r,s)
\param  deriv[][]  DOUBLE    (o)    shape function derivative at (r,s)

\return void

\author mf
\date 10/06
*/
void so3_shape_deriv(DIS_TYP typ,
                     DOUBLE r,
                     DOUBLE s,
                     DOUBLE t,
                     INT option,
                     DOUBLE shape[MAXNOD_SOLID3],
                     DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3])
{
/*----------------------------------------------------------------------*/
/* a few parameters depending on coord of current point (r,s,t) */
  DOUBLE rp, sp, tp;
  DOUBLE rm, sm, tm;
  DOUBLE rrm, ssm, ttm;
  DOUBLE rrp, ssp, ttp;
  DOUBLE u, r4, s4, t4, u4;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_deriv");
#endif

  /*--------------------------------------------------------------------*/
  /* switch according to discretisation type */
  switch (typ)
  {
    /*==================================================================*/
    /* hexahedral elements */
    /* linear interpolation */
    case hex8:
      /* auxiliary variables */
      rm  = 1.0 - r;
      rp  = 1.0 + r;
      sm  = 1.0 - s;
      sp  = 1.0 + s;
      tm  = 1.0 - t;
      tp  = 1.0 + t;
      /* shape functions at (r,s,t) */
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
        deriv[5][2] = 0.125*rp*sm;
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
    /* serendipity quadratic interpolation */
    case hex20:
      /* auxiliary variables */
      rm  = 1.0 - r;
      rp  = 1.0 + r;
      sm  = 1.0 - s;
      sp  = 1.0 + s;
      tm  = 1.0 - t;
      tp  = 1.0 + t;
      rrm = 1.0 - r*r;
      ssm = 1.0 - s*s;
      ttm = 1.0 - t*t;
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
      shape[12] = 0.25*rm*sm*ttm;
      shape[13] = 0.25*rp*sm*ttm;
      shape[14] = 0.25*rp*sp*ttm;
      shape[15] = 0.25*rm*sp*ttm;
      shape[16] = 0.25*rrm*sm*tp;
      shape[17] = 0.25*rp*ssm*tp;
      shape[18] = 0.25*rrm*sp*tp;
      shape[19] = 0.25*rm*ssm*tp;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* corners */
        deriv[0][0] = -0.125*sm*tm*(2.0*rm+sm+tm-5.0);
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
      }
      break;
    /* Lagrange quadratic interpolation */
    case hex27:
      /* auxiliary variables */
      rm  = 1.0 - r;
      rp  = 1.0 + r;
      sm  = 1.0 - s;
      sp  = 1.0 + s;
      tm  = 1.0 - t;
      tp  = 1.0 + t;
      rrm = 1.0 - 2.0*r;
      rrp = 1.0 + 2.0*r;
      ssm = 1.0 - 2.0*s;
      ssp = 1.0 + 2.0*s;
      ttm = 1.0 - 2.0*t;
      ttp = 1.0 + 2.0*t;
      /* corner nodes */
      shape[0] = -0.125* r*rm * s*sm * t*tm;
      shape[1] = 0.125 * r*rp * s*sm * t*tm;
      shape[2] = -0.125 * r*rp * s*sp * t*tm;
      shape[3] = 0.125 * r*rm * s*sp * t*tm;
      shape[4] = 0.125 * r*rm * s*sm * t*tp;
      shape[5] = -0.125 * r*rp * s*sm * t*tp;
      shape[6] = 0.125 * r*rp * s*sp * t*tp;
      shape[7] = -0.125 * r*rm * s*sp * t*tp;
      /* edge centre nodes */
      shape[8] = 0.25 * rm*rp * s*sm * t*tm;
      shape[9] = -0.25 * r*rp * sm*sp * t*tm;
      shape[10] = -0.25 * rm*rp * s*sp * t*tm;
      shape[11] = 0.25 * r*rm * sm*sp * t*tm;
      shape[12] = 0.25 * r*rm * s*sm * tm*tp;
      shape[13] = -0.25 * r*rp * s*sm * tm*tp;
      shape[14] = 0.25 * r*rp * s*sp * tm*tp;
      shape[15] = -0.25 * r*rm * s*sp * tm*tp;
      shape[16] = -0.25 * rm*rp * s*sm * t*tp;
      shape[17] = 0.25 * r*rp * sm*sp * t*tp;
      shape[18] = 0.25 * rm*rp * s*sp * t*tp;
      shape[19] = -0.25 * r*rm * sm*sp * t*tp;
      /* side centre nodes */
      shape[20] = -0.5 * rm*rp * sm*sp * t*tm;
      shape[21] = -0.5 * rm*rp * s*sm * tm*tp;
      shape[22] = 0.5 * r*rp * sm*sp * tm*tp;
      shape[23] = 0.5 * rm*rp * s*sp * tm*tp;
      shape[24] = -0.5 * r*rm * sm*sp * tm*tp;
      shape[25] = 0.5 * rm*rp * sm*sp * t*tp;
      /* volume centre node */
      shape[26] = rm*rp * sm*sp * tm*tp;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* shape fct. assoc. to node 0 */
        deriv[0][0] = -0.125 * rrm * s*sm * t*tm;
        deriv[0][1] = -0.125 * r*rm * ssm * t*tm;
        deriv[0][2] = -0.125 * r*rm * s*sm * ttm;
        /* shape fct. assoc. to node 1 */
        deriv[1][0] = 0.125 * rrp * s*sm * t*tm;
        deriv[1][1] = 0.125 * r*rp * ssm * t*tm;
        deriv[1][2] = 0.125 * r*rp * s*sm * ttm;
        /* shape fct. assoc. to node 2 */
        deriv[2][0] = -0.125 * rrp * s*sp * t*tm;
        deriv[2][1] = -0.125 * r*rp * ssp * t*tm;
        deriv[2][2] = -0.125 * r*rp * s*sp * ttm;
        /* shape fct. assoc. to node 3 */
        deriv[3][0] = 0.125 * rrm * s*sp * t*tm;
        deriv[3][1] = 0.125 * r*rm * ssp * t*tm;
        deriv[3][2] = 0.125 * r*rm * s*sp * ttm;
        /* shape fct. assoc. to node 4 */
        deriv[4][0] = 0.125 * rrm * s*sm * t*tp;
        deriv[4][1] = 0.125 * r*rm * ssm * t*tp;
        deriv[4][2] = 0.125 * r*rm * s*sm * ttp;
        /* shape fct. assoc. to node 5 */
        deriv[5][0] = -0.125 * rrp * s*sm * t*tp;
        deriv[5][1] = -0.125 * r*rp * ssm * t*tp;
        deriv[5][2] = -0.125 * r*rp * s*sm * ttp;
        /* shape fct. assoc. to node 6 */
        deriv[6][0] = 0.125 * rrp * s*sp * t*tp;
        deriv[6][1] = 0.125 * r*rp * ssp * t*tp;
        deriv[6][2] = 0.125 * r*rp * s*sp * ttp;
        /* shape fct. assoc. to node 7 */
        deriv[7][0] = -0.125 * rrm * s*sp * t*tp;
        deriv[7][1] = -0.125 * r*rm * ssp * t*tp;
        deriv[7][2] = -0.125 * r*rm * s*sp * ttp;
        /* shape fct. assoc. to node 8 */
        deriv[8][0] = -0.5 * r * s*sm * t*tm;
        deriv[8][1] = 0.25 * rm*rp * ssm * t*tm;
        deriv[8][2] = 0.25 * rm*rp * s*sm * ttm;
        /* shape fct. assoc. to node 9 */
        deriv[9][0] = -0.25 * rrp * sm*sp * t*tm;
        deriv[9][1] = 0.5 * r*rp * s * t*tm;
        deriv[9][2] = -0.25 * r*rp * sm*sp * ttm;
        /* shape fct. assoc. to node 10 */
        deriv[10][0] = 0.5 * r * s*sp * t*tm;
        deriv[10][1] = -0.25 * rm*rp * ssp * t*tm;
        deriv[10][2] = -0.25 * rm*rp * s*sp * ttm;
        /* shape fct. assoc. to node 11 */
        deriv[11][0] = 0.25 * rrm * sm*sp * t*tm;
        deriv[11][1] = -0.5 * r*rm * s * t*tm;
        deriv[11][2] = 0.25 * r*rm * sm*sp * ttm;
        /* shape fct. assoc. to node 12 */
        deriv[12][0] = 0.25 * rrm * s*sm * tm*tp;
        deriv[12][1] = 0.25 * r*rm * ssm * tm*tp;
        deriv[12][2] = -0.5 * r*rm * s*sm * t;
        /* shape fct. assoc. to node 13 */
        deriv[13][0] = -0.25 * rrp * s*sm * tm*tp;
        deriv[13][1] = -0.25 * r*rp * ssm * tm*tp;
        deriv[13][2] = 0.5 * r*rp * s*sm * t;
        /* shape fct. assoc. to node 14 */
        deriv[14][0] = 0.25 * rrp * s*sp * tm*tp;
        deriv[14][1] = 0.25 * r*rp * ssp * tm*tp;
        deriv[14][2] = -0.5 * r*rp * s*sp * t;
        /* shape fct. assoc. to node 15 */
        deriv[15][0] = -0.25 * rrm * s*sp * tm*tp;
        deriv[15][1] = -0.25 * r*rm * ssp * tm*tp;
        deriv[15][2] = 0.5 * r*rm * s*sp * t;
        /* shape fct. assoc. to node 16 */
        deriv[16][0] = 0.5 * r * s*sm * t*tp;
        deriv[16][1] = -0.25 * rm*rp * ssm * t*tp;
        deriv[16][2] = -0.25 * rm*rp * s*sm * ttp;
        /* shape fct. assoc. to node 17 */
        deriv[17][0] = 0.25 * rrp * sm*sp * t*tp;
        deriv[17][1] = -0.5 * r*rp * s * t*tp;
        deriv[17][2] = 0.25 * r*rp * sm*sp * ttp;
        /* shape fct. assoc. to node 18 */
        deriv[18][0] = -0.5 * r * s*sp * t*tp;
        deriv[18][1] = 0.25 * rm*rp * ssp * t*tp;
        deriv[18][2] = 0.25 * rm*rp * s*sp * ttp;
        /* shape fct. assoc. to node 19 */
        deriv[19][0] = -0.25 * rrm * sm*sp * t*tp;
        deriv[19][1] = 0.5 * r*rm * s * t*tp;
        deriv[19][2] = -0.25 * r*rm * sm*sp * ttp;
        /* shape fct. assoc. to node 20 */
        deriv[20][0] = r * sm*sp * t*tm;
        deriv[20][1] = rm*rp * s * t*tm;
        deriv[20][2] = -0.5 * rm*rp * sm*sp * ttm;
        /* shape fct. assoc. to node 21 */
        deriv[21][0] = r * s*sm * tm*tp;
        deriv[21][1] = -0.5 * rm*rp * ssm * tm*tp;
        deriv[21][2] = rm*rp * s*sm * t;
        /* shape fct. assoc. to node 22 */
        deriv[22][0] = 0.5 * rrp * sm*sp * tm*tp;
        deriv[22][1] = -r*rp * s * tm*tp;
        deriv[22][2] = -r*rp * sm*sp * t;
        /* shape fct. assoc. to node 23 */
        deriv[23][0] = -r * s*sp * tm*tp;
        deriv[23][1] = 0.5 * rm*rp * ssp * tm*tp;
        deriv[23][2] = -rm*rp * s*sp * t;
        /* shape fct. assoc. to node 24 */
        deriv[24][0] = -0.5 * rrm * sm*sp * tm*tp;
        deriv[24][1] = r*rm * s * tm*tp;
        deriv[24][2] = r*rm * sm*sp * t;
        /* shape fct. assoc. to node 25 */
        deriv[25][0] = -r * sm*sp * t*tp;
        deriv[25][1] = -rm*rp * s * t*tp;
        deriv[25][2] = 0.5 * rm*rp * sm*sp * ttp;
        /* shape fct. assoc. to node 26 */
        deriv[26][0] = -2.0 * r * sm*sp * tm*tp;
        deriv[26][1] = -2.0 * rm*rp * s * tm*tp;
        deriv[26][2] = -2.0 * rm*rp * sm*sp * t;
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
      /* auxiliary variables */
      u = 1.0 - r - s - t;
      r4 = 4.0*r;
      s4 = 4.0*s;
      t4 = 4.0*t;
      u4 = 4.0*u;
      /* shape functions at (r,s,t)
       * [T.J.R. Hughes, The finite element method, Dover, 2000] */
      /* corner nodes */
      shape[0] = r*(2.0*r - 1.0);
      shape[1] = s*(2.0*s - 1.0);
      shape[2] = t*(2.0*t - 1.0);
      shape[3] = u*(2.0*u - 1.0);
      /* middle nodes */
      shape[4] = r4*s;
      shape[5] = s4*t;
      shape[6] = r4*t;
      shape[7] = r4*u;
      shape[8] = s4*u;
      shape[9] = t4*u;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* shape fct. assoc. to node 0 */
        deriv[0][0] = r4 - 1.0;  /* differentiated with respect to r */
        deriv[0][1] = 0.0;       /* differentiated with respect to s */
        deriv[0][2] = 0.0;       /* differentiated with respect to t */
        /* shape fct. assoc. to node 1 */
        deriv[1][0] = 0.0;
        deriv[1][1] = s4 - 1.0;
        deriv[1][2] = 0.0;
        /* shape fct. assoc. to node 2 */
        deriv[2][0] = 0.0;
        deriv[2][1] = 0.0;
        deriv[2][2] = t4 - 1.0;
        /* shape fct. assoc. to node 3 */
        deriv[3][0] = 1.0 - u4;
        deriv[3][1] = 1.0 - u4;
        deriv[3][2] = 1.0 - u4;
        /* shape fct. assoc. to node 4 */
        deriv[4][0] = s4;
        deriv[4][1] = r4;
        deriv[4][2] = 0.0;
        /* shape fct. assoc. to node 5 */
        deriv[5][0] = 0.0;
        deriv[5][1] = t4;
        deriv[5][2] = s4;
        /* shape fct. assoc. to node 6 */
        deriv[6][0] = t4;
        deriv[6][1] = 0.0;
        deriv[6][2] = r4;
        /* shape fct. assoc. to node 7 */
        deriv[7][0] = u4 - r4;
        deriv[7][1] = -r4;
        deriv[7][2] = -r4;
        /* shape fct. assoc. to node 8 */
        deriv[8][0] = -s4;
        deriv[8][1] = u4 - s4;
        deriv[8][2] = -s4;
        /* shape fct. assoc. to node 9 */
        deriv[9][0] = -t4;
        deriv[9][1] = -t4;
        deriv[9][2] = u4 - t4;
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
} /* end of so3_shape_deriv */


/*======================================================================*/
/*!
\brief Initialise arrays for Gauss point coordinate, weights, shape
       functions + their derivatives

\param so3_gpshade SO3_GPSHAPEDERIV* (o)    coord & weights
                                            & shape functions
                                            & derivatives

\return void

\author bborn
\date 12/06
*/
void so3_shape_gpshade_init(SO3_GPSHAPEDERIV* so3_gpshade)
{
  INT idim;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_gpshade_init");
#endif

  so3_gpshade->distyp = dis_none;
  for (idim=0; idim<NDIM_SOLID3; idim++)
  {
    so3_gpshade->gpintc[idim] = 0;
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Shape functions and their natural derivatives at point (r,s,t)

\param  ele         ELEMENT*          (i)    pointer to current element
\param  data        SO3_DATA*         (i)    constant Gauss point data
\param  so3_gpshade SO3_GPSHAPEDERIV* (o)    coord & weights
                                             & shape functions
                                             & derivatives

\return void

\author bborn
\date 12/06
*/
void so3_shape_gpshade(ELEMENT *ele,
                       SO3_DATA *data,
                       SO3_GPSHAPEDERIV *so3_gpshade)
{
  INT calc_gpshade;  /* operation flag */
  INT gpnumr=0, gpnums=0, gpnumt=0;  /* auxiliar number of Gauss points */
  INT gpintcr=0, gpintcs=0, gpintct=0;  /* auxiliar integration case */
  INT igpr, igps, igpt;  /* directional Gauss point index */
  INT igp;  /* total Gauss point index (in domain) */
  INT idim;  /* dimension index */
  DOUBLE gpcr=1.0, gpcs=1.0, gpct=1.0;  /* Gauss point coordinate */
  DOUBLE fac=1.0;  /* Gauss weight */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_gpshade");
#endif

  /*--------------------------------------------------------------------*/
  /* check if current element discretisation + Gauss integration case
   * is identical to stored combination */
  calc_gpshade = 0;
  switch (ele->distyp)
  {
    /* hexahedron elements */
    case hex8: case hex20: case hex27:
      if (ele->distyp != so3_gpshade->distyp)
      {
        calc_gpshade = 1;
      }
      else
      {
        for (idim=0; idim<NDIM_SOLID3; idim++)
        {
          if (ele->e.so3->gpintc[idim] != so3_gpshade->gpintc[idim])
          {
            calc_gpshade = 1;
            break;
          }
        }
      }
      break;
    /* tetrahedron elements */
    case tet4: case tet10:
      if (ele->distyp != so3_gpshade->distyp)
      {
        calc_gpshade = 1;
      }
      else
      {
        if (ele->e.so3->gpintc[0] != so3_gpshade->gpintc[0])
        {
          calc_gpshade = 1;
        }
      }
      break;
    default:
      dserror("Discretisation type is not available\n");
      break;
  }

  /*--------------------------------------------------------------------*/
  /* determine shape functions at Gauss points */
  if (calc_gpshade)
  {
    /* set discretisation type of stored data */
    so3_gpshade->distyp = ele->distyp;
    /* select new Gauss point set */
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
        so3_gpshade->gpintc[0] = gpintcr;
        so3_gpshade->gpintc[1] = gpintcs;
        so3_gpshade->gpintc[2] = gpintct;
        so3_gpshade->gptot = gpnumr*gpnums*gpnumt;
        break;
      /* tetrahedra elements */
      /* tets are not simply rst-oriented and have just one GP-set nr */
      case tet4: case tet10:
        gpnumr = 1;
        gpnums = 1;
        gpnumt = ele->e.so3->gpnum[0];
        gpintcr = 1;
        gpintcs = 1;
        gpintct = ele->e.so3->gpintc[0];
        so3_gpshade->gpintc[0] = gpintct;
        so3_gpshade->gptot = gpnumt;
        break;
      default:
        dserror("ele->distyp unknown!");
    }
    /* initialise total domain GP index */
    igp = 0;
    /* walk along every Gauss point */
    /* WARNING: Do not change this ordering unless you change it
     *          also in so3_stress_extrpol() */
    for (igpr=0; igpr<gpnumr; igpr++)
    {
      for (igps=0; igps<gpnums; igps++)
      {
        for (igpt=0; igpt<gpnumt; igpt++)
        {
          /* retrieve Gauss point coordinate and weight */
          switch (ele->distyp)
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
              dserror("ele->distyp unknown!");
              break;
          }
          /* set coordinate and weight in total array */
          so3_gpshade->gpco[igp][0] = gpcr;
          so3_gpshade->gpco[igp][1] = gpcs;
          so3_gpshade->gpco[igp][2] = gpct;
          so3_gpshade->gpwg[igp] = fac;
          /* shape functions and their derivatives at igp */
          so3_shape_deriv(ele->distyp, gpcr, gpcs, gpct, 1,
                          so3_gpshade->gpshape[igp],
                          so3_gpshade->gpderiv[igp]);
          /* increment absolute Gauss point index */
          igp++;
        }
      }
    }
    /* verify total number of Gauss points */
    dsassert(so3_gpshade->gptot == igp,
             "Broken total Gauss point number\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/*!
\brief Test shape functions and their derivatives at element nodes

\param  data        SO3_DATA*         (i)    constant Gauss point data
\return void

\author bborn
\date 01/07
*/
#ifdef TEST_SOLID3
void so3_shape_test(SO3_DATA *data)
{
  FILE *filetest;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_test");
#endif

  /*--------------------------------------------------------------------*/
  /* open test file */
  filetest = fopen("so3_shape_test.testout", "w");

  /*--------------------------------------------------------------------*/
  /* hex8 */
  if (MAXNOD_SOLID3 >= 8)
  {
    /* test every shape function node by node */
    so3_shape_test_shp(data, hex8, "HEX8", 8, filetest);
    /* test every shape-function derivatives
     * by summing over nodes */
    so3_shape_test_drv(data, hex8, "HEX8", 8, filetest);
  }
  /* hex20 */
  if (MAXNOD_SOLID3 >= 20)
  {
    so3_shape_test_shp(data, hex20, "HEX20", 20, filetest);
    so3_shape_test_drv(data, hex20, "HEX20", 20, filetest);
  }
  /* hex27 */
  if (MAXNOD_SOLID3 >= 27)
  {
    so3_shape_test_shp(data, hex27, "HEX27", 27, filetest);
    so3_shape_test_drv(data, hex27, "HEX27", 27, filetest);
  }
  /* tet4 */
  if (MAXNOD_SOLID3 >= 4)
  {
    so3_shape_test_shp(data, tet4, "TET4", 4, filetest);
    so3_shape_test_drv(data, tet4, "TET4", 4, filetest);
  }
  /* tet10 */
  if (MAXNOD_SOLID3 >= 10)
  {
    so3_shape_test_shp(data, tet10, "TET10", 10, filetest);
    so3_shape_test_drv(data, tet10, "TET10", 10, filetest);
  }

  /*--------------------------------------------------------------------*/
  fclose(filetest);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif /* end #ifdef TEST_SOLID3 */


/*======================================================================*/
/*!
\brief Evaluate every shape function at every node
       We expect to get 1 if shape function corresponds
       to node otheriwse we should get 0

\param  data        SO3_DATA*         (i)    constant Gauss point data
\param  dis         DIS_TYP           (i)    discretisation type
\param  text        CHAR*             (i)    header line to print
\param  nelenod     INT               (i)    num. elem. nodes
\param  filetest    FILE*             (i/o)  output file
\return void

\author bborn
\date 01/07
*/
#ifdef TEST_SOLID3
void so3_shape_test_shp(SO3_DATA *data,
                        DIS_TYP dis,
                        CHAR *text,
                        INT nelenod,
                        FILE *filetest)
{
  INT inod, jnod, idim;  /* indices */
  DOUBLE rst[NDIM_SOLID3];  /* parameter coords of a node */
  DOUBLE shape[MAXNOD_SOLID3]; /* shape functions eval. at a node */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_test_shp");
#endif

  /*--------------------------------------------------------------------*/
  fprintf(filetest, "%s\n", text);
  for (inod=0; inod<nelenod; inod++)
  {
    /* get rst-coordinate of node */
    if ( (dis == hex8) || (dis == hex20) || (dis == hex27) )
    {
      for (idim=0; idim<NDIM_SOLID3; idim++)
      {
        rst[idim] = data->nodhrst[inod][idim];
      }
    }
    else if ( (dis == tet4) || (dis == tet10) )
    {
      for (idim=0; idim<NDIM_SOLID3; idim++)
      {
        rst[idim] = data->nodtrst[inod][idim];
      }
    }
    /* print node index plus its rst-coord */
    fprintf(filetest, "%d:(%.1f,%.1f,%.1f):  ",
            inod, rst[0], rst[1], rst[2]);
    /* get shape functions at current rst-triple
     * should return shape[inod] = 1
     * and           shape[jnod] = 0 for every jnod != inod */
    so3_shape_deriv(dis, rst[0], rst[1], rst[2], 0, shape, NULL);
    for (jnod=0; jnod<nelenod; jnod++)
    {
      fprintf(filetest, "%d:%d  ", jnod, (INT)shape[jnod]);
    }
    fprintf(filetest, "\n");
  }
  fprintf(filetest, "\n");

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif /* end #ifdef TEST_SOLID3 */


/*======================================================================*/
/*!
\brief Evaluate every derivatives of shape functions at every node
       We expect to get 0 for any point in any direction

\param  data        SO3_DATA*         (i)    constant Gauss point data
\param  dis         DIS_TYP           (i)    discretisation type
\param  text        CHAR*             (i)    header line to print
\param  nelenod     INT               (i)    num. elem. nodes
\param  filetest    FILE*             (i/o)  output file
\return void

\author bborn
\date 01/07
*/
#ifdef TEST_SOLID3
void so3_shape_test_drv(SO3_DATA *data,
                        DIS_TYP dis,
                        CHAR *text,
                        INT nelenod,
                        FILE *filetest)
{
  INT inod, jnod, idim, jdim;
  DOUBLE rst[NDIM_SOLID3];
  DOUBLE der[NDIM_SOLID3];
  DOUBLE shape[MAXNOD_SOLID3];
  DOUBLE deriv[MAXNOD_SOLID3][NDIM_SOLID3];

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_shape_test_drv");
#endif

  /*--------------------------------------------------------------------*/
  fprintf(filetest, "%s\n", text);
  for (inod=0; inod<nelenod; inod++)
  {
    if ( (dis == hex8) || (dis == hex20) || (dis == hex27) )
    {
      for (idim=0; idim<NDIM_SOLID3; idim++)
      {
        rst[idim] = data->nodhrst[inod][idim];
      }
    }
    else if ( (dis == tet4) || (dis == tet10) )
    {
      for (idim=0; idim<NDIM_SOLID3; idim++)
      {
        rst[idim] = data->nodtrst[inod][idim];
      }
    }
    fprintf(filetest, "%d:(%.1f,%.1f,%.1f):  ",
            inod, rst[0], rst[1], rst[2]);
    so3_shape_deriv(dis, rst[0], rst[1], rst[2], 1, shape, deriv);
    der[0] = 0.0;
    der[1] = 0.0;
    der[2] = 0.0;
    for (jnod=0; jnod<nelenod; jnod++)
    {
      for (jdim=0; jdim<NDIM_SOLID3; jdim++)
      {
        der[jdim] += deriv[jnod][jdim];
      }
    }
    for (jdim=0; jdim<NDIM_SOLID3; jdim++)
    {
      fprintf(filetest, "%f  ", der[jdim]);
    }
    fprintf(filetest, "\n");
  }
  fprintf(filetest, "\n");

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif /* end #ifdef TEST_SOLID3 */

/*======================================================================*/
#endif  /* end of #ifdef D_SOLID3 */



#endif
