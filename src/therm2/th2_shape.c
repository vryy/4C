/*======================================================================*/
/*!
\file
\brief contains the routine 'w1_funct_deriv' which calculates the shape
       functions and derivatives for a wall element

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
/* headers */
#include "../headers/standardtypes.h"
#include "therm2.h"

/*!
\addtogroup THERM2
*//*! @{ (documentation module open)*/

/*======================================================================*/
/*!
\brief Shape functions and their natural derivatives at point (r,s)

\param  *shape     DOUBLE    (o)    shape function value at (r,s)
\param  **deriv    DOUBLE    (o)    shape function derivative at (r,s)
\param  r          DOUBLE    (i)    r-coord of (r,s)
\param  s          DOUBLE    (i)    s-coord of (r,s)
\param  typ        DIS_TYP   (i)    discretisation type
\param  option     INT       (i)    option
                                    ==0 : only shape functions
                                    ==1 : shape functions and derivatives
\return void

\author bborn
\date 03/06
*/
void th2_shape_deriv(DOUBLE     *shape,
                     DOUBLE    **deriv,
                     DOUBLE      r,
                     DOUBLE      s,
                     DIS_TYP     typ,
                     INT         option)
{
  INT            i, ii;
  DOUBLE         rr,ss,rp,rm,sp,sm,r2,s2,t;
  DOUBLE         rh,sh,rs,rhp,rhm,shp,shm;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th2_shape_deriv");
#endif

  /*--------------------------------------------------------------------*/
  /* a few parameters depending on coord of current point (r,s) */
  rr = r*r;
  ss = s*s;
  rp = 1.0+r;
  rm = 1.0-r;
  sp = 1.0+s;
  sm = 1.0-s;
  r2 = 1.0-rr;
  s2 = 1.0-ss;
  /*--------------------------------------------------------------------*/
  /* switch according to discretisation type */
  switch (typ)
  {
    /*------------------------------------------------------------------*/
    /* quadrilateral elements */
    /* linear interpolation */
    case quad4:
      shape[0] = 0.25*rp*sp;
      shape[1] = 0.25*rm*sp;
      shape[2] = 0.25*rm*sm;
      shape[3] = 0.25*rp*sm;
      /* optionally include derivatives */
      if (option == 1)
      {
        deriv[0][0]= 0.25*sp;
        deriv[0][1]=-0.25*sp;
        deriv[0][2]=-0.25*sm;
        deriv[0][3]= 0.25*sm;
        deriv[1][0]= 0.25*rp;
        deriv[1][1]= 0.25*rm;
        deriv[1][2]=-0.25*rm;
        deriv[1][3]=-0.25*rp;
      }
      break;
    /* quadratic interpolation without central node (serendipity) */
    case quad8:
      shape[0] = 0.25*rp*sp;
      shape[1] = 0.25*rm*sp;
      shape[2] = 0.25*rm*sm;
      shape[3] = 0.25*rp*sm;
      shape[4] = 0.5*r2*sp;
      shape[5] = 0.5*rm*s2;
      shape[6] = 0.5*r2*sm;
      shape[7] = 0.5*rp*s2;
      shape[0] = shape[0] - 0.5*(shape[4] + shape[7]);
      /* optionally include derivatives */
      if (option == 1)
      {
        deriv[0][0]= 0.25*sp;
        deriv[0][1]=-0.25*sp;
        deriv[0][2]=-0.25*sm;
        deriv[0][3]= 0.25*sm;
        deriv[1][0]= 0.25*rp;
        deriv[1][1]= 0.25*rm;
        deriv[1][2]=-0.25*rm;
        deriv[1][3]=-0.25*rp;
        deriv[0][4]=-1.0*r*sp;
        deriv[0][5]=-0.5*  s2;
        deriv[0][6]=-1.0*r*sm;
        deriv[0][7]= 0.5*  s2;
        deriv[1][4]= 0.5*r2  ;
        deriv[1][5]=-1.0*rm*s;
        deriv[1][6]=-0.5*r2  ;
        deriv[1][7]=-1.0*rp*s;

        deriv[0][0]=deriv[0][0] - 0.5*(deriv[0][4] + deriv[0][7]);
        deriv[1][0]=deriv[1][0] - 0.5*(deriv[1][4] + deriv[1][7]);
      }
      for (i=1; i<=3; i++)
      {
        ii = i + 3;
        shape[i]=shape[i] - 0.5*(shape[ii] + shape[ii+1]);
        /* optionally include derivatives */
        if (option == 1)
        {
          deriv[0][i]=deriv[0][i] - 0.5*(deriv[0][ii] + deriv[0][ii+1]);
          deriv[1][i]=deriv[1][i] - 0.5*(deriv[1][ii] + deriv[1][ii+1]);
        }
      }
      break;
    /* quadratic interpolation with central node (Lagrangian) */
    case quad9:
      rh  = 0.5*r;
      sh  = 0.5*s;
      rs  = rh*sh;
      rhp = r+0.5;
      rhm = r-0.5;
      shp = s+0.5;
      shm = s-0.5;
      shape[0] = rs*rp*sp;
      shape[1] =-rs*rm*sp;
      shape[2] = rs*rm*sm;
      shape[3] =-rs*rp*sm;
      shape[4] = sh*sp*r2;
      shape[5] =-rh*rm*s2;
      shape[6] =-sh*sm*r2;
      shape[7] = rh*rp*s2;
      shape[8] = r2*s2;
      /* optionally include derivatives */
      if (option == 1)
      {
        deriv[0][0]= rhp*sh*sp;
        deriv[0][1]= rhm*sh*sp;
        deriv[0][2]=-rhm*sh*sm;
        deriv[0][3]=-rhp*sh*sm;
        deriv[0][4]=-2.0*r*sh*sp;
        deriv[0][5]= rhm*s2;
        deriv[0][6]= 2.0*r*sh*sm;
        deriv[0][7]= rhp*s2;
        deriv[0][8]=-2.0*r*s2;
        deriv[1][0]= shp*rh*rp;
        deriv[1][1]=-shp*rh*rm;
        deriv[1][2]=-shm*rh*rm;
        deriv[1][3]= shm*rh*rp;
        deriv[1][4]= shp*r2;
        deriv[1][5]= 2.0*s*rh*rm;
        deriv[1][6]= shm*r2;
        deriv[1][7]=-2.0*s*rh*rp;
        deriv[1][8]=-2.0*s*r2;
      }
      break;
    /*------------------------------------------------------------------*/
    /* triangular elements */
    /* linear interpolation */
    case tri3:
      shape[0]=1.0-r-s;
      shape[1]=r;
      shape[2]=s;
      /* optionally include derivatives */
      if (option == 1)
      {
        deriv[0][0]=-1.0;
        deriv[1][0]=-1.0;
        deriv[0][1]= 1.0;
        deriv[1][1]=0.0;
        deriv[0][2]=0.0;
        deriv[1][2]= 1.0;
      }
      break;
    /* Quadratic interpolation */
    case tri6:
      t = 1.0-r-s;
      shape[0] = t*(2.0*t-1.0);
      shape[1] = r*(2.0*r-1.0);
      shape[2] = s*(2.0*s-1.0);
      shape[3] = 4.0*r*t;
      shape[4] = 4.0*r*s;
      shape[5] = 4.0*s*t;
      /* optionally include derivatives */
      if (option == 1)
      {
        /* first natural derivative of shape[0] with respect to r */
        deriv[0][0] = -4.0*t + 1.0;
        /* first natural derivative of shape[0] with respect to s */
        deriv[1][0] = -4.0*t + 1.0;
        deriv[0][1] = 4.0*r - 1.0;
        deriv[1][1] = 0.0;
        deriv[0][2] = 0.0;
        deriv[1][2] = 4.0*s - 1.0;
        deriv[0][3] = 4.0*t - 4.0*r;
        deriv[1][3] = -4.0*r;
        deriv[0][4] = 4.0*s;
        deriv[1][4] = 4.0*r;
        deriv[0][5] = -4.0*s;
        deriv[1][5] = 4.0*t - 4.0*s;
      } /* end if (option == 1) */
      break;
    /* catch errenous discretisation types */
    default:
      dserror("unknown typ of interpolation");
      break;
  } /* end of switch typ */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th2_shape_deriv */


/*======================================================================*/
#endif  /* end of #ifdef D_THERM2 */



#endif
