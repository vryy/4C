/*!----------------------------------------------------------------------
\file
\brief xfem_calfuncderiv.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM
#include "../headers/standardtypes.h"
#include "xfem_prototypes.h"
/*!
\addtogroup XFEM
*//*! @{ (documentation module open)*/
static DOUBLE Q14 = ONE/FOUR;



/*!----------------------------------------------------------------------
\brief shape function and derivative evaluation for enriched fluid element

<pre>                                                            irhan 05/04
shape function and derivative evaluation for enriched fluid element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_funct(
  DOUBLE     *funct,
  DOUBLE    **deriv,
  DOUBLE    **deriv2,
  DOUBLE      r,
  DOUBLE      s,
  DIS_TYP     typ,
  DOUBLE     *lset01,
  INT         iel,
  INT         is_elcut
  )
{
  INT        i;
  DOUBLE     lset;
  DOUBLE     rp,rm,sp,sm;


#ifdef DEBUG
  dstrc_enter("xfem_f2_funct");
#endif
/*----------------------------------------------------------------------*/

  switch (typ)
  {
      case tri3:
        /* evaluate shape functions */
        funct[0] = r;
        funct[1] = s;
        funct[2] = ONE-r-s;

        /* evaluate first derivatives */
        deriv[0][0] =  ONE;
        deriv[1][0] =  ZERO;
        deriv[0][1] =  ZERO;
        deriv[1][1] =  ONE;
        deriv[0][2] = -ONE;
        deriv[1][2] = -ONE;
        break;
      case quad4:
        /*
         * partI => shape functions
         */
        rp = ONE+r;
        rm = ONE-r;
        sp = ONE+s;
        sm = ONE-s;

        funct[0] = Q14*rp*sp;
        funct[1] = Q14*rm*sp;
        funct[2] = Q14*rm*sm;
        funct[3] = Q14*rp*sm;
        /* first derivative evaluation */
        deriv[0][0] = Q14*sp;
        deriv[1][0] = Q14*rp;

        deriv[0][1] = -Q14*sp;
        deriv[1][1] = Q14*rm;

        deriv[0][2] = -Q14*sm;
        deriv[1][2] = -Q14*rm;

        deriv[0][3] = Q14*sm;
        deriv[1][3] = -Q14*rp;
        /* second derivative evaluation */
        deriv2[0][0] =  ZERO;
        deriv2[1][0] =  ZERO;
        deriv2[2][0] =  Q14;

        deriv2[0][1] =  ZERO;
        deriv2[1][1] =  ZERO;
        deriv2[2][1] = -Q14;

        deriv2[0][2] =  ZERO;
        deriv2[1][2] =  ZERO;
        deriv2[2][2] =  Q14;

        deriv2[0][3] =  ZERO;
        deriv2[1][3] =  ZERO;
        deriv2[2][3] = -Q14;
        break;
      default:
        dserror("distyp unknown");
  }

  /* check */
  if (is_elcut==1 || is_elcut==2)
  {
    /*
     * partII => enriched shape functions
     */

    /* compute distance function at Gauss point */
    lset = 0.0;
    for (i=0; i<iel; i++)
    {
      lset += funct[i]*lset01[i];
    }
    /* compute enriched shape function at Gauss point */
    for (i=0; i<iel; i++)
    {
      funct[iel+i] = funct[i]*(fabs(lset)-fabs(lset01[i]));
    }
  }
  else
  {
    /* compute enriched shape function at Gauss point */
    for (i=0; i<iel; i++)
    {
      funct[iel+i] = 0.0;
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_funct */



/*!----------------------------------------------------------------------
\brief shape function evaluation for enriched fluid element

<pre>                                                            irhan 05/04
shape function evaluation for enriched fluid element (used by level set element)
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_funct1(
  DOUBLE     *funct,
  DOUBLE      r,
  DOUBLE      s,
  DIS_TYP     typ,
  INT         iel,
  DOUBLE     *lset01,
  INT         is_elcut
  )
{
  INT        i;
  DOUBLE     lset;

#ifdef DEBUG
  dstrc_enter("xfem_f2_funct1");
#endif
/*----------------------------------------------------------------------*/

  switch (typ)
  {
      case tri3:
      case quad4:
        /*
         * partI => standard
         */
        if (is_elcut==1 || is_elcut==2)
        {
          /*
           * partII => enriched
           */

          /* compute distance function at Gauss point */
          lset = 0.0;
          for (i=0; i<iel; i++)
          {
            lset += funct[i]*lset01[i];
          }
          /* compute enriched shape function at Gauss point */
          for (i=0; i<iel; i++)
          {
            funct[iel+i] = funct[i]*(fabs(lset)-fabs(lset01[i]));
          }
        }
        else
        {
          /* compute enriched shape function at Gauss point */
          for (i=0; i<iel; i++)
          {
            funct[iel+i] = 0.0;
          }
        }
        break;
/*----------------------------------------------------------------------*/
      default:
        dserror("distyp unknown");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_funct1 */




/*!----------------------------------------------------------------------
\brief first global derivative evaluation for enriched fluid2 element

<pre>                                                            irhan 05/04
global derivative evaluation for enriched fluid2 element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_derxy(
  DOUBLE    **derxy,
  DOUBLE    **deriv,
  DOUBLE    **xjm,
  DOUBLE      det,
  INT         iel,
  DOUBLE     *funct,
  DOUBLE     *lset01,
  INT         is_elcut
  )
{
  INT        i,j;
  DOUBLE     lset,dlset[2];
  DOUBLE     sign;

#ifdef DEBUG
  dstrc_enter("xfem_f2_derxy");
#endif
/*----------------------------------------------------------------------*/

  /*
   * partI => standard
   */

  /* initialize */
  for(i=0; i<iel; i++)
  {
    derxy[0][i] = ZERO;
    derxy[1][i] = ZERO;
  }

  for (i=0; i<iel; i++)
  {
    derxy[0][i] += ( xjm[1][1]*deriv[0][i]  +
                     (-xjm[0][1]*deriv[1][i]))/det;
    derxy[1][i] +=((-xjm[1][0]*deriv[0][i]) +
                   xjm[0][0]*deriv[1][i])/det;
  }

  if (is_elcut==1 || is_elcut==2)
  {
    /*
     * partII => enriched
     */

    /* compute N*phi */
    lset = 0.0;
    for (i=0; i<iel; i++)
    {
      lset += funct[i]*lset01[i];
    }
    sign = lset/fabs(lset);
    /* compute gradN*phi */
    for (i=0; i<2; i++)
    {
      dlset[i] = 0.0;
      for (j=0; j<iel; j++)
      {
        dlset[i] += derxy[i][j]*lset01[j];
      }
    }
    /*
     * compute enriched global derivatives of shape functions
     * at Gauss point
     */
    for (i=0; i<2; i++)
    {
      for (j=0; j<iel; j++)
      {
        derxy[i][iel+j] = derxy[i][j]*(fabs(lset)-fabs(lset01[j])) +
          funct[j]*(sign*dlset[i]);
      }
    }
  }
  else
  {
    /*
     * compute enriched global derivatives of shape functions
     * at Gauss point
     */
    for (i=0; i<2; i++)
    {
      for (j=0; j<iel; j++)
      {
        derxy[i][iel+j] = 0.0;
      }
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_derxy */



/*!----------------------------------------------------------------------
\brief second global derivative evaluation for enriched fluid2 element

<pre>                                                            irhan 05/04
second global derivative evaluation for enriched fluid2 element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_derxy2(
  DOUBLE     **xyze,
  DOUBLE     **xjm,
  DOUBLE     **bi,
  DOUBLE     **xder2,
  DOUBLE     **derxy,
  DOUBLE     **derxy2,
  DOUBLE     **deriv2,
  INT          iel,
  DOUBLE     *funct,
  DOUBLE     *lset01,
  INT         is_elcut
  )
{
  INT        i,j;
  DOUBLE     lset,dlset[2],d2lset[3];
  DOUBLE     sign;
  INT        ind1,ind2;
  INT        index[2][3] = {
                            {1, 2, 1},
                            {1, 2, 2}
                           };


  DOUBLE x00,x01,x02,x10,x11,x12,x20,x21,x22;
  DOUBLE det,dum;
  DOUBLE r0,r1,r2;

#ifdef DEBUG
  dstrc_enter("xfem_f2_der2xy");
#endif
/*----------------------------------------------------------------------*/

  /*
   * partI => standard
   */

  /* calculate elements of jacobian_bar matrix */
  x00 = xjm[0][0]*xjm[0][0];
  x01 = xjm[0][1]*xjm[0][1];
  x02 = TWO*xjm[0][0]*xjm[0][1];
  x10 = xjm[1][0]*xjm[1][0];
  x11 = xjm[1][1]*xjm[1][1];
  x12 = TWO*xjm[1][0]*xjm[1][1];
  x20 = xjm[0][0]*xjm[1][0];
  x21 = xjm[0][1]*xjm[1][1];
  x22 = xjm[0][0]*xjm[1][1] + xjm[0][1]*xjm[1][0];

  /* inverse of jacobian_bar matrix */
  det =   x00*x11*x22 + x01*x12*x20 + x10*x21*x02
        - x20*x11*x02 - x00*x12*x21 - x01*x10*x22 ;
  dum = ONE/det;
  bi[0][0] =   dum*(x11*x22 - x12*x21);
  bi[1][0] =  -dum*(x10*x22 - x20*x12);
  bi[2][0] =   dum*(x10*x21 - x20*x11);
  bi[0][1] =  -dum*(x01*x22 - x21*x02);
  bi[1][1] =   dum*(x00*x22 - x02*x20);
  bi[2][1] =  -dum*(x00*x21 - x20*x01);
  bi[0][2] =   dum*(x01*x12 - x11*x02);
  bi[1][2] =  -dum*(x00*x12 - x10*x02);
  bi[2][2] =   dum*(x00*x11 - x01*x10);

  /* initialise */
  for (i=0; i<3; i++)
  {
    xder2[i][0]=ZERO;
    xder2[i][1]=ZERO;
  }
  for (i=0;i<iel;i++)
  {
    derxy2[0][i]=ZERO;
    derxy2[1][i]=ZERO;
    derxy2[2][i]=ZERO;
  }

  /* determine 2nd derivatives of coord.-functions */
  for (i=0; i<iel; i++)
  {
    xder2[0][0] += deriv2[0][i] * xyze[0][i];
    xder2[1][0] += deriv2[1][i] * xyze[0][i];
    xder2[2][0] += deriv2[2][i] * xyze[0][i];
    xder2[0][1] += deriv2[0][i] * xyze[1][i];
    xder2[1][1] += deriv2[1][i] * xyze[1][i];
    xder2[2][1] += deriv2[2][i] * xyze[1][i];
  }

  /* calculate second global derivatives */
  for (i=0; i<iel; i++)
  {
    r0 = deriv2[0][i] - xder2[0][0]*derxy[0][i] - xder2[0][1]*derxy[1][i];
    r1 = deriv2[1][i] - xder2[1][0]*derxy[0][i] - xder2[1][1]*derxy[1][i];
    r2 = deriv2[2][i] - xder2[2][0]*derxy[0][i] - xder2[2][1]*derxy[1][i];

    derxy2[0][i] += bi[0][0]*r0 + bi[0][1]*r1 + bi[0][2]*r2;
    derxy2[1][i] += bi[1][0]*r0 + bi[1][1]*r1 + bi[1][2]*r2;
    derxy2[2][i] += bi[2][0]*r0 + bi[2][1]*r1 + bi[2][2]*r2;
  }

  if (is_elcut==1 || is_elcut==2)
  {
    /*
     * partII => enriched
     */

    /* compute distance function at Gauss point */
    lset = 0.0;
    for (i=0; i<iel; i++)
    {
      lset += funct[i]*lset01[i];
    }
    sign = lset/fabs(lset);
    /* compute gradN*phi */
    for (i=0; i<2; i++)
    {
      dlset[i] = 0.0;
      for (j=0; j<iel; j++)
      {
        dlset[i] += derxy[i][j]*lset01[j];
      }
    }
    /* compute gradgradN*phi */
    for (i=0; i<3; i++)
    {
      d2lset[i] = 0.0;
      for (j=0; j<iel; j++)
      {
        d2lset[i] += derxy2[i][j]*lset01[j];
      }
    }
    /*
     * compute enriched second global derivatives of shape functions
     * at Gauss point
     */
    for (i=0; i<3; i++)
    {
      for (j=0; j<iel; j++)
      {
        ind1 = index[0][i]-1;
        ind2 = index[1][i]-1;
        derxy2[i][j+iel] =
          derxy2[i][j]*(fabs(lset)-fabs(lset01[j])) +
          derxy[ind1][j]*(sign*dlset[ind2]) +
          derxy[ind2][j]*(sign*dlset[ind1]) +
          funct[j]*(sign*d2lset[i]);
      }
    }
  }
  else
  {
    /*
     * compute enriched second global derivatives of shape functions
     * at Gauss point
     */
    for (i=0; i<3; i++)
      for (j=0; j<iel; j++)
        derxy2[i][j+iel] = 0.0;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_der2xy */



/*!----------------------------------------------------------------------
\brief evaluation of velocity components (at a Gauss point) for enriched fluid2 element

<pre>                                                            irhan 05/04
 evaluation of velocity components (at a Gauss point) for enriched fluid2 element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_veli(
  DOUBLE  *velint,
  DOUBLE  *funct,
  DOUBLE **evel,
  INT      iel
  )
{
  INT     i,j;

#ifdef DEBUG
  dstrc_enter("xfem_f2_veli");
#endif
/*---------------------------------------------------------------------*/

  for (i=0; i<2; i++)
  {
    velint[i] = 0.0;
    for (j=0; j<TWO*iel; j++)
    {
      velint[i] += funct[j]*evel[i][j];
    }
  }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_veli */



/*!----------------------------------------------------------------------
\brief evaluation of the global first derivatives of the velocity
components (at a Gauss point) for enriched fluid2 element

<pre>                                                            irhan 05/04
evaluation of the global first derivatives of the velocity
components (at a Gauss point) for enriched fluid2 element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_vder(
  DOUBLE **vderxy,
  DOUBLE **derxy,
  DOUBLE **evel,
  INT      iel
  )
{
  INT     i,j;

#ifdef DEBUG
  dstrc_enter("xfem_f2_vder");
#endif
/*---------------------------------------------------------------------*/

  /* loop */
  for (i=0; i<2; i++)
  {
    vderxy[0][i] = ZERO;
    vderxy[1][i] = ZERO;
    for (j=0; j<TWO*iel; j++)
    {
      vderxy[0][i] += derxy[i][j]*evel[0][j];
      vderxy[1][i] += derxy[i][j]*evel[1][j];
    }
  }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of xfem_f2_vder */



/*!----------------------------------------------------------------------
\brief evaluation of the global second derivatives of the velocity
components (at a Gauss point) for enriched fluid2 element

<pre>                                                            irhan 05/04
evaluation of the global second derivatives of the velocity
components (at a Gauss point) for enriched fluid2 element
</pre>

*----------------------------------------------------------------------*/
void xfem_f2_vder2(
  DOUBLE **vderxy2,
  DOUBLE **derxy2,
  DOUBLE **evel,
  INT      iel
  )
{
  INT     i,j;

#ifdef DEBUG
  dstrc_enter("xfem_f2_vderxy2");
#endif
/*---------------------------------------------------------------------*/

  for (i=0; i<3; i++)
  {
    vderxy2[0][i] = ZERO;
    vderxy2[1][i] = ZERO;
    for (j=0; j<TWO*iel; j++)
    {
      vderxy2[0][i] += derxy2[i][j]*evel[0][j];
      vderxy2[1][i] += derxy2[i][j]*evel[1][j];
    }
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of xfem_f2_vderxy2 */
/*! @} (documentation module close)*/
#endif
