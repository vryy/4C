/*!----------------------------------------------------------------------
\file
\brief contains the routine 'if_funcderiv' which calculates ansatzfunctions,
       its derivatives with respect to xi and the jacobi-determinant

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h" 

/*! 
\addtogroup INTERF
*/
/*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculates ansatzfunctions,
        its derivatives with respect to xi and the jacobi-determinant
<pre>                                                              ah 05/03 
This routine calculates ansatzfunctions,
             its derivatives with respect to xi and the jacobi-determinant

</pre>
\param   e1        DOUBLE  (I)   xi - coordinate of actual gaussian point
\param   typ       DIS_TYP (I)   quad4 or quad8
\param  *x_mid     DOUBLE  (I)   x-coordinates of fictive nodes on midline of element
\param  *y_mid     DOUBLE  (I)   y-coordinates of fictive nodes on midline of element
\param   b_parabel DOUBLE  (I)   y = a + b*x + c*x^2 
\param   c_parabel DOUBLE  (I)   y = a + b*x + c*x^2 
\param  *funct     DOUBLE  (O)   shape functions for [u] 
\param   co        DOUBLE  (O)   cosinus of angle bet x-dir and orient. of IF-ele
\param   si        DOUBLE  (O)   sinus of angle bet x-dir and orient. of IF-ele
\param  *det       DOUBLE  (O)   determinants of jacobian matrix  
\return void                                               

*----------------------------------------------------------------------*/
void if_funcderiv(DOUBLE  e1,
                  DIS_TYP typ,
                  DOUBLE *x_mid,
                  DOUBLE *y_mid,
                  DOUBLE  b_parabel,
                  DOUBLE  c_parabel,
                  DOUBLE *funct,
                  DOUBLE *co,
                  DOUBLE *si,
                  DOUBLE *det)
{
DOUBLE x_GP;
DOUBLE dy_dx;
DOUBLE dx_dxi;
DOUBLE deltax,deltay;
DOUBLE deriv[3];
DOUBLE beta;
DOUBLE alpha=0.0;

#ifdef DEBUG 
dstrc_enter("if_funcderiv");
#endif
/*----------------------------------------------------------------------*/
switch(typ)
{
case quad4:

   /*----------------------------------------------- ansatzfunctions ---*/
   funct[0] = (ONE/TWO) * (1 - e1);
   funct[1] = (ONE/TWO) * (1 + e1);
   deltax = x_mid[1] - x_mid[0];
   deltay = y_mid[1] - y_mid[0];
   /*-------------------------------------- Jacobidet = ds/dxi = L/2 ---*/
   *det = sqrt(deltax*deltax + deltay*deltay)/TWO;
   /*--------------------------------------- cosinus and sinus alpha ---*/
   beta = atan( deltay / deltax ); /* angle between X and xi direction in [0,..pi/4] */
   if(deltay >=0 && deltax >0)           alpha = beta;
   else if (deltay > 0 && deltax ==0)    alpha =  PI/TWO;
   else if (deltax <0)                   alpha = beta + PI;
   else if (deltay < 0 && deltax ==0)    alpha =  (3*PI)/TWO;
   else if (deltay < 0 && deltax >0)     alpha =  beta + TWO*PI;
   
break;
/*-----------------------------------------------------------------------*/
case quad8:

   /*---------------------------- ansatzfunctions and its derivatives ---*/
   funct[0] = (ONE/TWO) * (e1*e1 - e1);
   funct[1] = (ONE/TWO) * (e1*e1 + e1);
   funct[2] = (ONE - e1*e1);
   deriv[0] = e1 - ONE/TWO;
   deriv[1] = e1 + ONE/TWO;
   deriv[2] = - TWO*e1;
   /*-------------------------------------------- x-coordinate of GP ---*/
   x_GP   = funct[0]*x_mid[0]+funct[1]*x_mid[1]+funct[2]*x_mid[2];
   /*------------------------------------------------------ f'=dy/dx ---*/
   dy_dx  = b_parabel + 2*c_parabel*x_GP;
   /*-------------------------------------------------------- dx/dxi ---*/
   dx_dxi = deriv[0]*x_mid[0]+deriv[1]*x_mid[1]+deriv[2]*x_mid[2];
   /*---------------------- Jacobidet = ds/dxi=(1 + f'^2)^0,5*dx/dxi ---*/
   *det   = sqrt(1+dy_dx*dy_dx)*FABS(dx_dxi);
   /*--------------------------------------- cosinus and sinus alpha ---*/
   beta   = atan(dy_dx);
   if(dy_dx >=0 && dx_dxi >0)           alpha = beta;
   else if (dy_dx > 0 && dx_dxi ==0)    alpha =  PI/TWO;
   else if (dx_dxi <0)                  alpha = beta + PI;
   else if (dy_dx < 0 && dx_dxi ==0)    alpha =  (3*PI)/TWO;
   else if (dy_dx < 0 && dx_dxi >0)     alpha =  beta + TWO*PI;
   
  
break;
default:
   dserror("discretisation unknown for Interface");
break;
}

*co = cos(alpha);
*si = sin(alpha);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of if_funcderiv */
/*----------------------------------------------------------------------*/


#endif /*D_INTERF*/
/*! @} (documentation module close)*/
