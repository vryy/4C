/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_funct_deriv' which calculates the 
shape functions for a 3d ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"


/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the shape functions and derivatives

<pre>                                                              mn 06/02 
This routine calcuates the shape functions and their derivatives at a
point r,s,t for an 3D-ale-element.

</pre>
\param *funct  double  (o)   shape functions
\param **deriv double  (o)   the derivatives of the shape functions
\param r       double  (i)   r coordinate
\param s       double  (i)   s coordinate
\param t       double  (i)   t coordinate
\param typ     DIS_TYP (i)   type of dicretization
\param option  int     (i)   option == 0 : only functions,
                             option == 1 : also derivatives

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale3_static_ke

*----------------------------------------------------------------------*/
void ale3_funct_deriv(double     *funct, 
                       double    **deriv, 
                       double      r, 
                       double      s,
                       double      t,
                       DIS_TYP     typ,
                       int         option)
{
int            i, ii;
const double   q12 = 1.0/2.0;
const double   q14 = 1.0/4.0;
const double   q16 = 1.0/6.0;
const double   q18 = 1.0/8.0;
const double   q64 = 1.0/64.0;
const double   q964= 9.0/64.0;
double         rr,ss,tt,rp,sp,tp,rm,sm,tm,rrm,ssm,ttm;
#ifdef DEBUG 
dstrc_enter("ale3_funct_deriv");
#endif
/*----------------------------------------------------------------------*/
/* if option ==0 only funtion evaluation, if option==1 also derivatives */
/*----------------------------------------------------------------------*/
rp  = 1.0+r;
rm  = 1.0-r;
sp  = 1.0+s;
sm  = 1.0-s;
tp  = 1.0+t;
tm  = 1.0-t;
rrm = 1.0-r*r;
ssm = 1.0-s*s;
ttm = 1.0-t*t;

switch(typ)
{
/*----------------------------------------------------------------------*
 |  L I N E A r     sHAPE FUNCtIONs AND tHEIr NAtUrAL DErIVAtIVEs       |
 |  sErENDIPItY     ( 8-NODED ELEMENt )                                 |
 *---------------------------------------------------------------------*/
case hex8:
   funct[0]=q18*rp*sm*tm;
   funct[1]=q18*rp*sp*tm;
   funct[2]=q18*rm*sp*tm;
   funct[3]=q18*rm*sm*tm;
   funct[4]=q18*rp*sm*tp;
   funct[5]=q18*rp*sp*tp;
   funct[6]=q18*rm*sp*tp;
   funct[7]=q18*rm*sm*tp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= q18*sm*tm  ;
      deriv[0][1]= q18*sp*tm  ;
      deriv[0][2]=-deriv[0][1];
      deriv[0][3]=-deriv[0][0];
      deriv[0][4]= q18*sm*tp  ;
      deriv[0][5]= q18*sp*tp  ;
      deriv[0][6]=-deriv[0][5];
      deriv[0][7]=-deriv[0][4];
      deriv[1][0]=-q18*tm*rp  ;
      deriv[1][1]=-deriv[1][0];
      deriv[1][2]= q18*tm*rm  ;
      deriv[1][3]=-deriv[1][2];
      deriv[1][4]=-q18*tp*rp  ;
      deriv[1][5]=-deriv[1][4];
      deriv[1][6]= q18*tp*rm  ;
      deriv[1][7]=-deriv[1][6];
      deriv[2][0]=-q18*rp*sm  ;
      deriv[2][1]=-q18*rp*sp  ;
      deriv[2][2]=-q18*rm*sp  ;
      deriv[2][3]=-q18*rm*sm  ;
      deriv[2][4]=-deriv[2][0];
      deriv[2][5]=-deriv[2][1];
      deriv[2][6]=-deriv[2][2];
      deriv[2][7]=-deriv[2][3];
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
} /* end of ale3_funct_deriv */
#endif
/*! @} (documentation module close)*/
