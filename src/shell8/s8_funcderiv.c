#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | shape functions and derivatives                            m.gee 6/01|
 *----------------------------------------------------------------------*/
void s8_funct_deriv(double     *funct, 
                       double    **deriv, 
                       double      r, 
                       double      s,
                       DIS_TYP     typ,
                       int         option)
{
int            i;
const double   q12 = 1.0/2.0;
const double   q14 = 1.0/4.0;
const double   q16 = 1.0/6.0;
double         rr,ss,rp,rm,sp,sm,r2,s2;
double         rh,sh,rs,rhp,rhm,shp,shm;
#ifdef DEBUG 
dstrc_enter("s8_funct_deriv");
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
   if (option==1)
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
case quad8:
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
case tri3:
break;
case tri6:
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
} /* end of s8_funct_deriv */
#endif





