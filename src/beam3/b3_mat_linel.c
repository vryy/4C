/*!----------------------------------------------------------------------
\file
\brief contains the routines 'b3_mat_linel' which calculates the
constitutive matrix for a beam element

<pre>
Maintainer: Frank Huber
            huber@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/huber/
            0771 - 685-6120
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
\brief calculates the constitutive matrix for a beam

<pre>                                                              fh 10/02
This routine calculates the linear elastic constitutive matrix for a
spatial 1d-Timoshenko-beam element

</pre>
\param *ele      ELEMENT   (i)  actual element
\param ym        DOUBLE    (i)  Youngs Modulus
\param pv        DOUBLE    (i)  Poissons ratio
\param **d       DOUBLE    (o)  Constitutive matrix


\warning There is nothing special in this routine
\return void
\sa calling:   ---;
    called by: b3_call_mat()

*----------------------------------------------------------------------*/
void b3_mat_linel(ELEMENT *ele,
                  DOUBLE   ym,
                  DOUBLE   pv,
		  DOUBLE **d)
{
/* Values for cross section and material */
DOUBLE e, g, k, a, iy, iz, it, h, gay, gaz, gas;
INT    ike;

#ifdef DEBUG
dstrc_enter("b3_mat_linel");
#endif
/*----------get value for shear correction factor k=1/gs ----------------*/
ike=ele->e.b3->ike;
k = ele->e.b3->gs;
k = 1./k;

e=ym;
g=ym/(2*(1+pv));
if (ike==2)
{
   a = ele->e.b3->area;
   iy= ele->e.b3->iuu;
   iz= ele->e.b3->ivv;
   it= ele->e.b3->it;
   h = ele->e.b3->length;
   gas=k*g*a;

   /*----------Residual Bending Flexibility: see Hughes p. 378-------------*/
/*   gay=1./(1./gas + h*h/(12.*e*iy));
   gaz=1./(1./gas + h*h/(12.*e*iz));*/
   gay=gas;
   gaz=gas;

   d[0][0]=e*a;
   d[1][1]=gay;
   d[2][2]=gaz;
   d[3][3]=g*it;
   d[4][4]=e*iy;
   d[5][5]=e*iz;
}
else
{
   d[0][0]=e;
   d[1][1]=g*k;
   d[2][2]=g*k;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_mat_linel */
#endif
/*! @} (documentation module close)*/
