/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_funct_deriv' which calculates the shape
       functions and derivatives for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*
 | shape functions and derivatives                               al 9/01|
 *----------------------------------------------------------------------*/
void w1_funct_deriv(DOUBLE     *funct,
                    DOUBLE    **deriv,
                    DOUBLE      r,
                    DOUBLE      s,
                    DIS_TYP     typ,
                    INT         option)
{
INT            i, ii;
const DOUBLE   q12 = ONE/TWO;
const DOUBLE   q14 = ONE/FOUR;
DOUBLE         rr,ss,rp,rm,sp,sm,r2,s2,t;
DOUBLE         rh,sh,rs,rhp,rhm,shp,shm;
#ifdef DEBUG
dstrc_enter("w1_funct_deriv");
#endif
/*----------------------------------------------------------------------*/
/* if option ==0 only funtion evaluation, if option==1 also derivatives */
/*----------------------------------------------------------------------*/
rr = r*r;
ss = s*s;
rp = 1.0+r;
rm = 1.0-r;
sp = 1.0+s;
sm = 1.0-s;
r2 = 1.0-rr;
s2 = 1.0-ss;
switch(typ)
{
/*------------------------------------------------ rectangular elements */
case quad4:/*------------------------- linear rectangular interpolation */
   funct[0] = q14*rp*sp;
   funct[1] = q14*rm*sp;
   funct[2] = q14*rm*sm;
   funct[3] = q14*rp*sm;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= q14*sp;
      deriv[0][1]=-q14*sp;
      deriv[0][2]=-q14*sm;
      deriv[0][3]= q14*sm;
      deriv[1][0]= q14*rp;
      deriv[1][1]= q14*rm;
      deriv[1][2]=-q14*rm;
      deriv[1][3]=-q14*rp;
   }
break;
case quad8:/*------------- quadratic interpolation without central node */
   funct[0] = q14*rp*sp;
   funct[1] = q14*rm*sp;
   funct[2] = q14*rm*sm;
   funct[3] = q14*rp*sm;
   funct[4] = q12*r2*sp;
   funct[5] = q12*rm*s2;
   funct[6] = q12*r2*sm;
   funct[7] = q12*rp*s2;
   funct[0] = funct[0] - q12*(funct[4] + funct[7]);
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= q14*sp;
      deriv[0][1]=-q14*sp;
      deriv[0][2]=-q14*sm;
      deriv[0][3]= q14*sm;
      deriv[1][0]= q14*rp;
      deriv[1][1]= q14*rm;
      deriv[1][2]=-q14*rm;
      deriv[1][3]=-q14*rp;
      deriv[0][4]=-ONE*r*sp;
      deriv[0][5]=-q12*  s2;
      deriv[0][6]=-ONE*r*sm;
      deriv[0][7]= q12*  s2;
      deriv[1][4]= q12*r2  ;
      deriv[1][5]=-ONE*rm*s;
      deriv[1][6]=-q12*r2  ;
      deriv[1][7]=-ONE*rp*s;

      deriv[0][0]=deriv[0][0] - q12*(deriv[0][4] + deriv[0][7]);
      deriv[1][0]=deriv[1][0] - q12*(deriv[1][4] + deriv[1][7]);
   }
   for (i=1; i<=3; i++)
   {
      ii=i + 3;
      funct[i]=funct[i] - q12*(funct[ii] + funct[ii+1]);
      if (option==1)              /*--- check for derivative evaluation ---*/
      {
          deriv[0][i]=deriv[0][i] - q12*(deriv[0][ii] + deriv[0][ii+1]);
          deriv[1][i]=deriv[1][i] - q12*(deriv[1][ii] + deriv[1][ii+1]);
      }
   }
break;
case quad9:/*---------------- quadratic interpolation with central node */
   rh  = q12*r;
   sh  = q12*s;
   rs  = rh*sh;
   rhp = r+q12;
   rhm = r-q12;
   shp = s+q12;
   shm = s-q12;
   funct[0] = rs*rp*sp;
   funct[1] =-rs*rm*sp;
   funct[2] = rs*rm*sm;
   funct[3] =-rs*rp*sm;
   funct[4] = sh*sp*r2;
   funct[5] =-rh*rm*s2;
   funct[6] =-sh*sm*r2;
   funct[7] = rh*rp*s2;
   funct[8] = r2*s2;
   if (option==1)
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
/*------------------------------------------------- triangular elements */
case tri3: /* LINEAR shape functions and their natural derivatives -----*/
/*----------------------------------------------------------------------*/
   funct[0]=ONE-r-s;
   funct[1]=r;
   funct[2]=s;

   if(option==1) /* --> first derivative evaluation */
   {
      deriv[0][0]=-ONE;
      deriv[1][0]=-ONE;
      deriv[0][1]= ONE;
      deriv[1][1]=ZERO;
      deriv[0][2]=ZERO;
      deriv[1][2]= ONE;
   } /* endif (option==1) */
break;
/*-------------------------------------------------------------------------*/
case tri6: /* Quadratic shape functions and their natural derivatives -----*/
    t = ONE-r-s;

    funct[0] = t*(TWO*t-ONE);
    funct[1] = r*(TWO*r-ONE);
    funct[2] = s*(TWO*s-ONE);
    funct[3] = FOUR*r*t;
    funct[4] = FOUR*r*s;
    funct[5] = FOUR*s*t;
    
    if (option == 1) /* --> first derivative evaluation */
    {
        /* first natural derivative of funct[0] with respect to r */
        deriv[0][0] = -FOUR*t + ONE;
        /* first natural derivative of funct[0] with respect to s */
        deriv[1][0] = -FOUR*t + ONE;
        deriv[0][1] = FOUR*r - ONE;
        deriv[1][1] = ZERO;
        deriv[0][2] = ZERO;
        deriv[1][2] = FOUR*s - ONE;
        deriv[0][3] = FOUR*t - FOUR*r;
        deriv[1][3] = -FOUR*r;
        deriv[0][4] = FOUR*s;
        deriv[1][4] = FOUR*r;
        deriv[0][5] = -FOUR*s;
        deriv[1][5] = FOUR*t - FOUR*s;
    } /* end if (option==1) */
break;
default:
   dserror("unknown typ of interpolation");
break;
} /* end of switch typ */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_funct_deriv */

/*----------------------------------------------------------------------*
 | shape functions and derivatives                        genk 10/02    |
 | for a degenerated 2D isoparametric element                           |
 *----------------------------------------------------------------------*/
void w1_degfuncderiv(DOUBLE     *funct,
                      DOUBLE    **deriv,
                      DOUBLE      r,
                      DIS_TYP     typ,
                      INT         option)
{
const DOUBLE   q12 = ONE/TWO;
DOUBLE         rr,rp,rm,r2;
#ifdef DEBUG
dstrc_enter("w1_edgefuncderiv");
#endif
/*----------------------------------------------------------------------*/
/* if option ==0 only funtion evaluation, if option==1 also derivatives */
/*----------------------------------------------------------------------*/
rr = r*r;
rp = ONE+r;
rm = ONE-r;
r2 = 1.0-rr;

switch(typ)
{
/*------------------------------------------------ rectangular elements */
case quad4: case tri3:/*-------------- linear rectangular interpolation */
   funct[0] = q12*rm;
   funct[1] = q12*rp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -q12;
      deriv[0][1]=  q12;
   }
break;
case quad8: case quad9: case tri6:/*---------- quadratic interpolation  */
/*   dserror("quadratic interpolation for deg. funcs not checked yet!\n");*/
   funct[0] = -q12*r*rm;
   funct[1] = r2;
   funct[2] = q12*r*rp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -q12+r;
      deriv[0][1]= -TWO*r;
      deriv[0][2]=  q12+r;
   }
break;
default:
   dserror("unknown typ of interpolation");
break;
} /* end of switch typ */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_funct_deriv */

#endif



