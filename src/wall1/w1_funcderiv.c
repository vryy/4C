#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*----------------------------------------------------------------------*
 | shape functions and derivatives                               al 9/01|
 *----------------------------------------------------------------------*/
void w1_funct_deriv(double     *funct, 
                    double    **deriv, 
                    double      r, 
                    double      s,
                    DIS_TYP   typ,
                    int         option)
{
int            i, ii;
const double   q12 = 1.0/2.0;
const double   q14 = 1.0/4.0;
const double   q16 = 1.0/6.0;
double         rr,ss,rp,rm,sp,sm,r2,s2;
double         rh,sh,rs,rhp,rhm,shp,shm;
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
      deriv[0][4]=-1.0*r*sp;
      deriv[0][5]=-q12*  s2;
      deriv[0][6]=-1.0*r*sm;
      deriv[0][7]= q12*  s2;
      deriv[1][4]= q12*r2  ;
      deriv[1][5]=-1.0*rm*s;
      deriv[1][6]=-q12*r2  ;
      deriv[1][7]=-1.0*rp*s;
      
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
} /* end of w1_funct_deriv */
#endif



