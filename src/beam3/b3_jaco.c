/*!----------------------------------------------------------------------
\file
\brief contains the routine 'b3_jaco' which calculates the Jacobian Matrix
at point r,s,t and the corresponding determinant

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
\brief calculates the Jacobian Matrix at point r,s,t and the corresponding
determinant

<pre>                                                              fh 02/03
This routine calculates the Jacobian Matrix at point r,s,t and the
corresponding determinant of the actual Timoshenko beam element

</pre>
\param *funct    DOUBLE    (i)  values of shape functions at actual GP
\param **deriv   DOUBLE    (i)  values of derivatives at actual GP
\param **xjm     DOUBLE    (o)  Jacobian Matrix
\param **ijm     DOUBLE    (o)  Inverse of Jacobian Matrix
\param **V       DOUBLE    (i)  Vs, Vt directions of local vectors s,t
\param *det      DOUBLE    (o)  determinant of Jacobian Matrix
\param s         DOUBLE    (i)  coordinate in local s-direction
\param t         DOUBLE    (i)  coordinate in local t-direction
\param *ele      ELEMENT   (i)  actual element
\param iel       INT       (i)  number of element nodes


\warning This routine is not used yet.
\return void
\sa calling:   ---;
    called by: b3_cal_ele()

*----------------------------------------------------------------------*/
void b3_jaco(DOUBLE     *funct,
             DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE    **ijm,
	     DOUBLE    **V,
	     DOUBLE     *det,
             DOUBLE      s,
	     DOUBLE      t,
	     ELEMENT    *ele,
             INT         iel)
{
/*----------------------------------------------------------------------*/
INT i,j,k;        /* some loopers */
DOUBLE b,a;       /* a = element height, b = element width */
DOUBLE N=0.;      /* value for shape function at actual node */
DOUBLE dN;        /* value for shape function derivative at actual node */
NODE *actnode;    /* actual element node */
#ifdef DEBUG
dstrc_enter("b3_jaco");
#endif
/*---------------------------------------------- initialize jacobian ---*/
/*      CALL MXIR8 (XJM,NUMDFE,NUMDFE,ZERO) */
/*---------------------------------- determine jacobian at point r,s ---*/
a = ele->e.b3->height;
b = ele->e.b3->width;
xjm[0][0] = 0.0;
xjm[0][1] = 0.0;
xjm[0][2] = 0.0;
xjm[1][0] = 0.0;
xjm[1][1] = 0.0;
xjm[1][2] = 0.0;
xjm[2][0] = 0.0;
xjm[2][1] = 0.0;
xjm[2][2] = 0.0;
/* constant vector Vs, Vt over the element (no curvature inside element */
for (k=0; k<iel; k++)
{
     actnode = ele->node[k];
     dN=deriv[0][k];
     xjm[0][0] += dN * actnode->x[0] + 0.5*t*a*dN*V[2][0] + 0.5*s*b*dN*V[1][0];
     xjm[0][1] += dN * actnode->x[1] + 0.5*t*a*dN*V[2][1] + 0.5*s*b*dN*V[1][1];
     xjm[0][2] += dN * actnode->x[2] + 0.5*t*a*dN*V[2][2] + 0.5*s*b*dN*V[1][2];
     N += funct[k];
}
     xjm[1][0] = 0.5*b*N*V[1][0];
     xjm[1][1] = 0.5*b*N*V[1][1];
     xjm[1][2] = 0.5*b*N*V[1][2];
     xjm[2][0] = 0.5*a*N*V[2][0];
     xjm[2][1] = 0.5*a*N*V[2][1];
     xjm[2][2] = 0.5*a*N*V[2][2];

for (i=0; i<3; i++)
{  for (j=0; j<3; j++) ijm[i][j]=xjm[i][j];
}
/*-----calculation of inverse and determinant of jacobian matrix -------*/
math_inv3(ijm,det);
if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of b3_jaco */
#endif
/*! @} (documentation module close)*/
