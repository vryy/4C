/*!----------------------------------------------------------------------
\file
\brief contains the routine 'if_funcderiv' which calculates ansatzfunctions,
       its derivatives with respect to xi and the jacobi-determinant

*----------------------------------------------------------------------*/
#ifdef D_INTERF
#include "../headers/standardtypes.h"
#include "interf.h"
#include "interf_prototypes.h" 

/*! 
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculates ansatzfunctions,
        its derivatives with respect to xi and the jacobi-determinant
<pre>                                                              mn 05/03 
This routine calculates ansatzfunctions,
             its derivatives with respect to xi and the jacobi-determinant

</pre>
\param *ele  ELEMENT  (o)   the element

\warning buffer[50] is not needed locally
\return void                                               
\sa caling:    ---; 
    called by: inp_struct_field()

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
INT  i;
DOUBLE ell;
DOUBLE x_GP;
DOUBLE dy_dx;
DOUBLE dx_dxi;
DOUBLE deltax,deltay;
DOUBLE deriv[3];
DOUBLE beta,alpha;

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



/*----------------------------------------------------------------------*/
void if_func_bope(DOUBLE   e1,
                  DOUBLE  *x_mid,
                  DOUBLE  *y_mid,
                  DOUBLE  *functe,
                  DOUBLE  *coe,
                  DOUBLE  *sie,
                  DOUBLE  *dete,
                  INT      flag,
                  DOUBLE **bope)
{
INT  i;
DOUBLE ell;
DOUBLE x_GP;
DOUBLE dy_dx;
DOUBLE dx_dxi;
DOUBLE deltax,deltay,deltay_help;
DOUBLE beta,alpha;
DOUBLE derive[4];

#ifdef DEBUG 
dstrc_enter("if_func_bope");
#endif
/*----------------------------------------------------------------------*/
deltax = x_mid[1] - x_mid[0];
deltay = y_mid[1] - y_mid[0];
/*-------------------------------------- Jacobidet = ds/dxi = L/2 ---*/
*dete = sqrt(deltax*deltax + deltay*deltay)/TWO;
/*--------------------------------------- cosinus and sinus alpha ---*/
beta = atan( deltay / deltax ); /* angle between X and xi direction in [0,..pi/4] */

if(flag==1)
{
 functe[0] = (ONE/TWO) *(ONE/TWO) * (1 - e1);
 functe[1] = (ONE/TWO) *(ONE/TWO) * (1 + e1);
 functe[2] = (ONE/TWO) *(ONE/TWO) * (1 + e1);
 functe[3] = (ONE/TWO) *(ONE/TWO) * (1 - e1);

 derive[0] = -(ONE/TWO) *(ONE/TWO);
 derive[1] =  (ONE/TWO) *(ONE/TWO);
 derive[2] =  (ONE/TWO) *(ONE/TWO);
 derive[3] = -(ONE/TWO) *(ONE/TWO);
 
}
else if(flag==2)
{
 functe[0] = (ONE/TWO) *(ONE/TWO) * (1 - e1);
 functe[1] = (ONE/TWO) *(ONE/TWO) * (1 - e1);
 functe[2] = (ONE/TWO) *(ONE/TWO) * (1 + e1);
 functe[3] = (ONE/TWO) *(ONE/TWO) * (1 + e1);

 derive[0] = -(ONE/TWO) *(ONE/TWO);
 derive[1] = -(ONE/TWO) *(ONE/TWO);
 derive[2] =  (ONE/TWO) *(ONE/TWO);
 derive[3] =  (ONE/TWO) *(ONE/TWO);
}
if(deltax !=0)
{
 bope[0][0]= derive[0]* (TWO/deltax);
 bope[0][1]= derive[1]* (TWO/deltax);
 bope[0][2]= derive[2]* (TWO/deltax);
 bope[0][3]= derive[3]* (TWO/deltax);
}
else if(deltax ==0)
{
 bope[0][0]= 0.0;
 bope[0][1]= 0.0;
 bope[0][2]= 0.0;
 bope[0][3]= 0.0;
}

if(deltay !=0)
{

 bope[1][0]= derive[0]* (TWO/deltay);
 bope[1][1]= derive[1]* (TWO/deltay);
 bope[1][2]= derive[2]* (TWO/deltay);
 bope[1][3]= derive[3]* (TWO/deltay);

}
else if(deltay ==0)
{
 deltay_help=1.0E-50;
 bope[1][0]= derive[0]* (TWO/deltay_help);
 bope[1][1]= derive[1]* (TWO/deltay_help);
 bope[1][2]= derive[2]* (TWO/deltay_help);
 bope[1][3]= derive[3]* (TWO/deltay_help);
}


 if(deltay >=0 && deltax >0)           alpha = beta;
 else if (deltay > 0 && deltax ==0)    alpha =  PI/TWO;
 else if (deltax <0)                   alpha = beta + PI;
 else if (deltay < 0 && deltax ==0)    alpha =  (3*PI)/TWO;
 else if (deltay < 0 && deltax >0)     alpha =  beta + TWO*PI;


*coe = cos(alpha);
*sie = sin(alpha);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of if_func_bope */
/*----------------------------------------------------------------------*/







#endif /*D_INTERF*/
/*! @} (documentation module close)*/
