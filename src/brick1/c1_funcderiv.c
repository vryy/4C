/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1_funct_deriv' which calclate 
       shape functions and derivatives for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief compute shape functions and derivatives

<pre>                                                              al 06/02
This routine calcuates shape functions and derivatives for an 3D-hex-element.

</pre>
\param    funct   DOUBLE*  (o)   shape functions  
\param    deriv   DOUBLE** (o)   the derivatives of the shape functions  
\param        r   DOUBLE   (i)   r coordinate           
\param        s   INT      (i)   s coordinate           
\param        t   INT      (i)   t coordinate   
\param      typ   INT      (i)   type of dicretization   
\param   option   INT      (i)   option == 0 : only functions  
                                option == 1 : also derivatives

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_funct_deriv(double     *funct, 
                    double    **deriv, 
                    double          r, 
                    double          s,
                    double          t,
                    int           typ,
                    int         option)
{
const double   q12 = 1.0/2.0;
const double   q14 = 1.0/4.0;
const double   q16 = 1.0/6.0;
const double   q18 = 1.0/8.0;
const double   q64 = 1.0/64.0;
const double   q964= 9.0/64.0;
double         rp,sp,tp,rm,sm,tm,rrm,ssm,ttm;
#ifdef DEBUG 
dstrc_enter("c1_funct_deriv");
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
 |  l i n e a r     shape functions and their natural derivatives       |
 |  serendipity     ( 8-noded element )                                 |
 *---------------------------------------------------------------------*/
case 8:
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
/*----------------------------------------------------------------------*
 |   q u a d r a t i c  shape functions and their natural derivatives   |
 |   serendipity        ( 20-noded element )                            |
 *---------------------------------------------------------------------*/
case 20:
   funct[0]  = q18*rp*sm*tm*(rp+sm+tm-5.0);
   funct[1]  = q18*rp*sp*tm*(rp+sp+tm-5.0);
   funct[2]  = q18*rm*sp*tm*(rm+sp+tm-5.0);
   funct[3]  = q18*rm*sm*tm*(rm+sm+tm-5.0);
   funct[4]  = q18*rp*sm*tp*(rp+sm+tp-5.0);
   funct[5]  = q18*rp*sp*tp*(rp+sp+tp-5.0);
   funct[6]  = q18*rm*sp*tp*(rm+sp+tp-5.0);
   funct[7]  = q18*rm*sm*tp*(rm+sm+tp-5.0);
   funct[8]  = q14*rp*ssm*tm;
   funct[9]  = q14*rrm*sp*tm;
   funct[10] = q14*rm*ssm*tm;
   funct[11] = q14*rrm*sm*tm;
   funct[12] = q14*rp*ssm*tp;
   funct[13] = q14*rrm*sp*tp;
   funct[14] = q14*rm*ssm*tp;
   funct[15] = q14*rrm*sm*tp;
   funct[16] = q14*rp*sm*ttm;
   funct[17] = q14*rp*sp*ttm;
   funct[18] = q14*rm*sp*ttm;
   funct[19] = q14*rm*sm*ttm;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
     deriv[0][0] = q18*sm*tm*(2.0*rp+sm+tm-5.0); 
     deriv[0][1] = q18*sp*tm*(2.0*rp+sp+tm-5.0); 
     deriv[0][2] =-q18*sp*tm*(2.0*rm+sp+tm-5.0); 
     deriv[0][3] =-q18*sm*tm*(2.0*rm+sm+tm-5.0); 
     deriv[0][4] = q18*sm*tp*(2.0*rp+sm+tp-5.0); 
     deriv[0][5] = q18*sp*tp*(2.0*rp+sp+tp-5.0); 
     deriv[0][6] =-q18*sp*tp*(2.0*rm+sp+tp-5.0); 
     deriv[0][7] =-q18*sm*tp*(2.0*rm+sm+tp-5.0); 
     deriv[0][8] = q14*ssm*tm                  ; 
     deriv[0][9] =-q12*r*sp*tm                 ; 
     deriv[0][10]=-deriv[0][8]                 ; 
     deriv[0][11]=-q12*r*sm*tm                 ; 
     deriv[0][12]= q14*ssm*tp                  ; 
     deriv[0][13]=-q12*r*sp*tp                 ; 
     deriv[0][14]=-deriv[0][12]                ; 
     deriv[0][15]=-q12*r*sm*tp                 ; 
     deriv[0][16]= q14*sm*ttm                  ; 
     deriv[0][17]= q14*sp*ttm                  ; 
     deriv[0][18]=-deriv[0][17]                ; 
     deriv[0][19]=-deriv[0][16]                ; 
     deriv[1][0] =-q18*tm*rp*(2.0*sm+tm+rp-5.0); 
     deriv[1][1] = q18*tm*rp*(2.0*sp+tm+rp-5.0); 
     deriv[1][2] = q18*tm*rm*(2.0*sp+tm+rm-5.0); 
     deriv[1][3] =-q18*tm*rm*(2.0*sm+tm+rm-5.0); 
     deriv[1][4] =-q18*tp*rp*(2.0*sm+tp+rp-5.0); 
     deriv[1][5] = q18*tp*rp*(2.0*sp+tp+rp-5.0); 
     deriv[1][6] = q18*tp*rm*(2.0*sp+tp+rm-5.0); 
     deriv[1][7] =-q18*tp*rm*(2.0*sm+tp+rm-5.0); 
     deriv[1][8] =-q12*s*tm*rp                 ; 
     deriv[1][9] = q14*rrm*tm                  ; 
     deriv[1][10]=-q12*s*tm*rm                 ; 
     deriv[1][11]=-deriv[1][9]                 ; 
     deriv[1][12]=-q12*s*tp*rp                 ; 
     deriv[1][13]= q14*rrm*tp                  ; 
     deriv[1][14]=-q12*s*tp*rm                 ; 
     deriv[1][15]=-deriv[1][13]                ; 
     deriv[1][16]=-q14*ttm*rp                  ; 
     deriv[1][17]=-deriv[1][16]                ; 
     deriv[1][18]= q14*ttm*rm                  ; 
     deriv[1][19]=-deriv[1][18]                ; 
     deriv[2][0] =-q18*rp*sm*(2.0*tm+rp+sm-5.0); 
     deriv[2][1] =-q18*rp*sp*(2.0*tm+rp+sp-5.0); 
     deriv[2][2] =-q18*rm*sp*(2.0*tm+rm+sp-5.0); 
     deriv[2][3] =-q18*rm*sm*(2.0*tm+rm+sm-5.0); 
     deriv[2][4] = q18*rp*sm*(2.0*tp+rp+sm-5.0); 
     deriv[2][5] = q18*rp*sp*(2.0*tp+rp+sp-5.0); 
     deriv[2][6] = q18*rm*sp*(2.0*tp+rm+sp-5.0); 
     deriv[2][7] = q18*rm*sm*(2.0*tp+rm+sm-5.0); 
     deriv[2][8] =-q14*ssm*rp                  ; 
     deriv[2][9] =-q14*rrm*sp                  ; 
     deriv[2][10]=-q14*ssm*rm                  ; 
     deriv[2][11]=-q14*rrm*sm                  ; 
     deriv[2][12]=-deriv[2][8]                 ; 
     deriv[2][13]=-deriv[2][9]                 ; 
     deriv[2][14]=-deriv[2][10]                ; 
     deriv[2][15]=-deriv[2][11]                ; 
     deriv[2][16]=-q12*t*rp*sm                 ; 
     deriv[2][17]=-q12*t*rp*sp                 ; 
     deriv[2][18]=-q12*t*rm*sp                 ; 
     deriv[2][19]=-q12*t*rm*sm                 ; 
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
} /* end of c1_funct_deriv */
#endif
/*! @} (documentation module close)*/

































































