/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | shape functions and derivatives                            m.gee 6/01|
 *----------------------------------------------------------------------*/
void s8_funct_deriv(
    DOUBLE     *funct,
    DOUBLE    **deriv,
    DOUBLE      r,
    DOUBLE      s,
    DIS_TYP     typ,
    INT         option)
{

  const DOUBLE   q12 = 1.0/2.0;
  const DOUBLE   q14 = 1.0/4.0;
  DOUBLE         rr,ss,rp,rm,sp,sm,r2,s2;
  DOUBLE         rh,sh,rs,rhp,rhm,shp,shm;

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
      funct[0] = -q14*(1-r)*(1-s)*(1+r+s);
      funct[1] = -q14*(1+r)*(1-s)*(1-r+s);
      funct[2] = -q14*(1+r)*(1+s)*(1-r-s);
      funct[3] = -q14*(1-r)*(1+s)*(1+r-s);
      funct[4] =  q12*(1-r*r)*(1-s);
      funct[5] =  q12*(1+r)*(1-s*s);
      funct[6] =  q12*(1-r*r)*(1+s);
      funct[7] =  q12*(1-r)*(1-s*s);
      if (option==1)
      {
        deriv[0][0]=  q14*(1-s)*(2*r+s);
        deriv[0][1]=  q14*(1-s)*(2*r-s);
        deriv[0][2]=  q14*(1+s)*(2*r+s);
        deriv[0][3]=  q14*(1+s)*(2*r-s);
        deriv[0][4]= -r*(1-s);
        deriv[0][5]=  q12*(1-s*s);
        deriv[0][6]= -r*(1+s);
        deriv[0][7]= -q12*(1-s*s);

        deriv[1][0]=  q14*(1-r)*(r+2*s);
        deriv[1][1]=  q14*(1+r)*(-r+2*s);
        deriv[1][2]=  q14*(1+r)*(r+2*s);
        deriv[1][3]=  q14*(1-r)*(-r+2*s);
        deriv[1][4]= -q12*(1-r*r);
        deriv[1][5]= -s*(1+r);
        deriv[1][6]=  q12*(1-r*r);
        deriv[1][7]= -s*(1-r);
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


    /* triangular elements */
    case tri3:
      funct[0]=1-r-s;
      funct[1]=r;
      funct[2]=s;
      if (option==1)
      {
        deriv[0][0]=  -1.0;
        deriv[0][1]=  1.0;
        deriv[0][2]=  0.0;

        deriv[1][0]=  -1.0;
        deriv[1][1]=  0.0;
        deriv[1][2]=  1.0;
      }
      break;

    case tri6:
      funct[0]=(1-2*r-2*s)*(1-r-s);
      funct[1]=2*r*r-r;
      funct[2]=2*s*s-s;
      funct[3]=4*(r-r*r-r*s);
      funct[4]=4*r*s;
      funct[5]=4*(s-s*s-s*r);
      if (option==1)
      {
        deriv[0][0]= -3.0+4.0*r+4.0*s;
        deriv[0][1]= 4.0*r-1.0;
        deriv[0][2]= 0.0;
        deriv[0][3]= 4.0*(1-2.0*r-s);
        deriv[0][4]= 4.0*s;
        deriv[0][5]= -4.0*s;

        deriv[1][0]= -3.0+4.0*r+4.0*s;
        deriv[1][1]= 0.0;
        deriv[1][2]= 4.0*s-1.0;
        deriv[1][3]= -4.0*r;
        deriv[1][4]= 4.0*r;
        deriv[1][5]= 4.0*(1.0-2.0*s-r);
      }
/*      dserror("tri3 not possible for shell8!!"); */
      break;

    default:
      dserror("unknown typ of interpolation");
      break;
  } /* end of switch typ */


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of s8_funct_deriv */
#endif





