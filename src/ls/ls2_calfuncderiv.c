/*!----------------------------------------------------------------------
\file
\brief ls2_calfuncderiv.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/
static DOUBLE Q14 = ONE/FOUR;



/*!----------------------------------------------------------------------
\brief shape function and derivative evaluation for 2D level set element

<pre>                                                            irhan 05/04
</pre>

*----------------------------------------------------------------------*/
void ls2_funct(
  DOUBLE     *funct,
  DOUBLE    **deriv,
  DOUBLE      r,
  DOUBLE      s,
  DIS_TYP     typ
  )
{
  DOUBLE     rp,rm,sp,sm;

#ifdef DEBUG
  dstrc_enter("ls2_funct");
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
        rp=ONE+r;
        rm=ONE-r;
        sp=ONE+s;
        sm=ONE-s;
        /* evaluate shape functions */
        funct[0]=Q14*rp*sp;
        funct[1]=Q14*rm*sp;
        funct[2]=Q14*rm*sm;
        funct[3]=Q14*rp*sm;
        /* evaluate first derivatives */
   	deriv[0][0]= Q14*sp;
        deriv[1][0]= Q14*rp;

        deriv[0][1]=-Q14*sp;
        deriv[1][1]= Q14*rm;

        deriv[0][2]=-Q14*sm;
        deriv[1][2]=-Q14*rm;

        deriv[0][3]= Q14*sm;
        deriv[1][3]=-Q14*rp;
        break;
      default:
        dserror("distyp unknown");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls2_funct */



/*!----------------------------------------------------------------------
\brief computation of Jacobian and its determinant

<pre>                                                            irhan 05/04
computation of Jacobian and its determinant
</pre>

*----------------------------------------------------------------------*/
void ls2_jaco(
  DOUBLE    **xyze,
  DOUBLE     *funct,
  DOUBLE    **deriv,
  DOUBLE    **xjm,
  DOUBLE     *det,
  INT         iel,
  ELEMENT    *ele
  )
{
  INT     i;

#ifdef DEBUG
  dstrc_enter("ls2_jaco");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  xjm[0][0] = ZERO ;
  xjm[0][1] = ZERO ;
  xjm[1][0] = ZERO ;
  xjm[1][1] = ZERO ;
  /* Jacobian */
  for (i=0; i<iel; i++)
  {
    xjm[0][0] += deriv[0][i]*xyze[0][i];
    xjm[0][1] += deriv[0][i]*xyze[1][i];
    xjm[1][0] += deriv[1][i]*xyze[0][i];
    xjm[1][1] += deriv[1][i]*xyze[1][i];
  }
  /* determinant of Jacobian */
  *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

  if(*det<ZERO)
  {
    printf("\n");
    printf("GLOBAL ELEMENT %i\n",ele->Id);
    printf("NEGATIVE JACOBIAN DETERMINANT\n");
    printf("STOPPING ...");
#ifdef PARALLEL
    dserror("not regulary!\n");
#endif
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls2_jaco */


/*!----------------------------------------------------------------------
\brief global derivative evaluation for 2D level set element

<pre>                                                            irhan 05/04
global derivative evaluation for 2D level set element
</pre>

*----------------------------------------------------------------------*/
void ls2_gder(
  DOUBLE   **derxy,
  DOUBLE   **deriv,
  DOUBLE   **xjm,
  DOUBLE     det,
  INT        iel
  )
{
  INT     i;

#ifdef DEBUG
  dstrc_enter("ls2_gder");
#endif
/*----------------------------------------------------------------------*/

  /* initialize */
  for(i=0; i<iel; i++)
  {
    derxy[0][i]=ZERO;
    derxy[1][i]=ZERO;
  }
  /* construct */
  for (i=0; i<iel; i++)
  {
    derxy[0][i] =(  xjm[1][1]*deriv[0][i]  + (-xjm[0][1]*deriv[1][i]))/det;
    derxy[1][i] =((-xjm[1][0]*deriv[0][i]) +   xjm[0][0]*deriv[1][i] )/det;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls2_gder */
/*! @} (documentation module close)*/
#endif
