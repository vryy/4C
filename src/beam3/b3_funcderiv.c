/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_funct_deriv' which calculates the values
of the shape functions and the correspondent derivatives at point r

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM3
#include "../headers/standardtypes.h"
#include "beam3.h"
#include "beam3_prototypes.h"

/*!
\addtogroup BEAM3
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the values of the shape functions and the correspondent
derivatives at point r

<pre>                                                              fh 10/02
This routine calculates the values of the shape functions and the
correspondent derivatives if necessary at the actual gauss point r

</pre>
\param *funct    DOUBLE    (o)  value for shape functions at point r
\param **deriv   DOUBLE    (o)  value for derivatives at point r
\param r         DOUBLE    (i)  actual gauss point coordinate
\param typ       DIS_TYP   (i)  LIN2, LIN3 or LIN4
\param option    INT       (i)  flag for derivatives (1) or not (0)


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_funct_deriv(DOUBLE     *funct,
                    DOUBLE    **deriv,
                    DOUBLE      r,
                    DIS_TYP     typ,
                    INT         option)
{
const DOUBLE   q12 = 1.0/2.0;
const DOUBLE   q116 = 1.0/16.0;
DOUBLE         r2,r3,rp,rm;     /* some terms for shape functions */

#ifdef DEBUG
dstrc_enter("b3_funct_deriv");
#endif
/*----------------------------------------------------------------------*/
/* if option ==0 only function evaluation, if option==1 also derivatives*/
/*----------------------------------------------------------------------*/
r2 = r*r;
r3 = r*r*r;
rp = 1.0+r;
rm = 1.0-r;
/*------------------------------------------------------- beam elements */
/*--------------------------------- linear beam interpolation ----------*/
switch(typ)
{
case line2:/*-----------------------linear interpolation-----------------*/
   funct[0] = q12*rm;
   funct[1] = q12*rp;
   if (option==1)  	    /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -q12;
      deriv[0][1]= +q12;
   }
break;
case line3:/*------------------------quadratic interpolation ------------*/
   funct[0] = q12*rm-q12*(1.0-r2);
   funct[1] = q12*rp-q12*(1.0-r2);
   funct[2] = 1.0-r2;
   if (option==1)           /*--- check for derivative evaluation ---*/
   {
      deriv[0][0] = -q12 + r;
      deriv[0][1] = +q12 + r;
      deriv[0][2] = -2.0*r;
   }
break;
default:
break;
} /* end of switch typ */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_funct_deriv */
#endif
/*! @} (documentation module close)*/
