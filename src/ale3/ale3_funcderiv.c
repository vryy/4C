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
point r,s,t for 3D hex and tet ale-elements.

</pre>
\param *funct  DOUBLE  (o)   shape functions
\param **deriv DOUBLE  (o)   the derivatives of the shape functions
\param r       DOUBLE  (i)   r coordinate
\param s       DOUBLE  (i)   s coordinate
\param t       DOUBLE  (i)   t coordinate
\param typ     DIS_TYP (i)   type of dicretization
\param option  INT     (i)   option == 0 : only functions,
                             option == 1 : also derivatives

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale3_static_ke

*----------------------------------------------------------------------*/
void ale3_funct_deriv(DOUBLE     *funct, 
                       DOUBLE    **deriv, 
                       DOUBLE      r, 
                       DOUBLE      s,
                       DOUBLE      t,
                       DIS_TYP     typ,
                       INT         option)
{
INT            i, ii;
const DOUBLE   q12 = 1.0/2.0;
const DOUBLE   q14 = 1.0/4.0;
const DOUBLE   q16 = 1.0/6.0;
const DOUBLE   q18 = 1.0/8.0;
const DOUBLE   q64 = 1.0/64.0;
const DOUBLE   q964= 9.0/64.0;
DOUBLE         rr,ss,tt,rp,sp,tp,rm,sm,tm,rrm,ssm,ttm;
DOUBLE         t1,t2,t3,t4;
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
case tet4: /* LINEAR shape functions and their natural derivatives -----*/
/*--------------------------------------------------- form basic values */
   t1=r;
   t2=s;
   t3=t;
   t4=ONE-r-s-t;
      
   funct[0]= t1;
   funct[1]= t2;
   funct[2]= t3;
   funct[3]= t4;
   
   if(option==1) /* --> first derivative evaluation */
   {         
      deriv[0][0]= ONE;
      deriv[0][1]= ZERO;
      deriv[0][2]= ZERO;
      deriv[0][3]=-ONE;

      deriv[1][0]= ZERO;
      deriv[1][1]= ONE;
      deriv[1][2]= ZERO;
      deriv[1][3]=-ONE;

      deriv[2][0]= ZERO;
      deriv[2][1]= ZERO;
      deriv[2][2]= ONE;
      deriv[2][3]=-ONE;
   } /* endif (option==1) */  
break;
   
case tet10: /*  QUADRATIC shape functions and their natural derivatives */
			 
   dserror("shape functions for tet10 not yet implemented \n"); 
/*--------------------------------------------------- form basic values */
/*   t1=r;
   t2=s;
   t3=t;
   t4=ONE-r-s-t;
   
   funct[0] =  ;
   funct[1] =  ;
   funct[2] = ;
   funct[3] = ;
   funct[4] = ;
   funct[5] = ;
   funct[6] = ;
   funct[7] = ;
   funct[8] = ;
   funct[9] = ;


   if(option==1) /* --> first derivative evaluation */
/*   {
/*      deriv[0][0] = ;
      deriv[1][0] = ;
      deriv[2][0] = ;

      deriv[0][1] = ;
      deriv[1][1] = ;
      deriv[2][1] = ;
		
      deriv[0][2] = ;
      deriv[1][2] = ;
      deriv[2][2] = ;

      deriv[0][3] = ;
      deriv[1][3] = ;
      deriv[2][3] = ;

      deriv[0][4] = ;
      deriv[1][4] = ;
      deriv[2][4] = ;

      deriv[0][5] = ;
      deriv[1][5] = ;
      deriv[2][5] = ;

      deriv[0][6] = ;
      deriv[1][6] = ;
      deriv[2][6] = ;

      deriv[0][7] = ;
      deriv[1][7] = ;
      deriv[2][7] = ;

      deriv[0][8] = ;
      deriv[1][8] = ;
      deriv[2][8] = ;

      deriv[0][9] = ;
      deriv[1][9] = ;
      deriv[2][9] = ;
   }*/
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
