/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_hourglass' which calculates the 
hourglass stabilization matrix for a 2D ale element

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"
/*----------------------------------------------------------------------*
 |                                                         mn 06/02     |
 | vector of material laws                                              |
 | defined in global_control.c
*----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the additional stiffness matrix for hourglass stabilization

<pre>                                                              mn 06/02 
This routine calcuates the additional stiffness matrix for hourglass
stabilization for a 2D element.
  
see also:                                                             
   T. Belytschko and L.P. Bindeman:                                      
   Assumed strain stabilization of the 4-node quadrilateral with 1-point 
   quadrature for nonlinear problems.                                    
   Comp. Meth. Appl. Mech. Eng. 88 (1991) p. 311-340.                    
        
</pre>
\param *ele  ELEMENT  (i)   the element
\param **s   DOUBLE   (i/o) (i) the one point quadrature matrix
                            (o) the complete, stabilized stiffness matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_hourglass(ELEMENT  *ele,
		    DOUBLE  **s) 
{

MATERIAL      *actmat;

INT            i, j;

DOUBLE         h[4];
DOUBLE         x[4];
DOUBLE         y[4];
DOUBLE         y24,y13,x24,x13;

DOUBLE         hxx,hyy,hxy;

DOUBLE         c1,c2,c3;

DOUBLE         ee,nu,mu;

DOUBLE         bx[4];
DOUBLE         by[4];

DOUBLE         hx,hy;

DOUBLE         gam[4];
DOUBLE         gg[4][4];

DOUBLE         aa;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_hourglass");
#endif
/*----------------------------------------------------------------------*
 | material data: ee .. youngs modulus                                  |
 |                nu .. poisson ratio                                   |
 |                mu .. shear modulus                                   |
 *----------------------------------------------------------------------*/
actmat = &(mat[ele->mat-1]);
ee = actmat->m.stvenant->youngs;
nu = actmat->m.stvenant->possionratio;

mu = ee / (2*(1+nu));


/*----------------------------------------------------------------------*
 | special hourglass mode vector                                        |
 *----------------------------------------------------------------------*/
h[0] = 1.0;
h[1] = -1.0;
h[2] = 1.0;
h[3] = -1.0;

hx = 0.0;
hy = 0.0;


/*----------------------------------------------------------------------*
 | constants for ASOI                                                   |
 *----------------------------------------------------------------------*/
c1 = 4.0 * mu;
c2 = 0.0;
c3 = -4.0 * mu;



for(i=0; i<4; i++)
{
  x[i] = ele->node[i]->x[0];
  y[i] = ele->node[i]->x[1];

  hx = hx + h[i] * x[i];
  hy = hy + h[i] * y[i];
}

y24 = y[1] - y[3];
y13 = y[0] - y[2];
x24 = x[1] - x[3];
x13 = x[0] - x[2];


/*----------------------------------------------------------------------*
 | element area                                                         |
 *----------------------------------------------------------------------*/
aa = 0.5*(x13*y24 - x24*y13);


/*----------------------------------------------------------------------*
 | B-operator elements                                                  |
 *----------------------------------------------------------------------*/
bx[0] = + 1/(2*aa) * y24;
bx[1] = - 1/(2*aa) * y13;
bx[2] = - 1/(2*aa) * y24;
bx[3] = + 1/(2*aa) * y13;

by[0] = - 1/(2*aa) * x24;
by[1] = + 1/(2*aa) * x13;
by[2] = + 1/(2*aa) * x24;
by[3] = - 1/(2*aa) * x13;

hxx = 2/(3*aa) * (y13*y13 + y24*y24);
hyy = 2/(3*aa) * (x13*x13 + x24*x24);
hxy = -2/(3*aa) * (x13*y13 + x24*y24);


/*----------------------------------------------------------------------*
 | general hourglass mode vector                                        |
 *----------------------------------------------------------------------*/
for(i=0; i<4; i++)
{
  gam[i] = 0.25 * ( h[i] - hx * bx[i] - hy * by[i]);
}


/*----------------------------------------------------------------------*
 | generating the stabilization matrix and adding it to the one-point   |
 | quadrature matrix                                                    |
 | rearranging order of dofs                                            |
 *----------------------------------------------------------------------*/
for(i=0; i<4; i++)
{
  for(j=0; j<4; j++)
  {
    gg[i][j] = gam[i] * gam[j];
    s[2*i][2*j]  = s[2*i][2*j] + (c1*hxx + c2*hyy) * gg[i][j];
    s[2*i+1][2*j]  = s[2*i+1][2*j] + (c3*hxy) * gg[i][j];
    s[2*i][2*j+1]  = s[2*i][2*j+1] + (c3*hxy) * gg[i][j];
    s[2*i+1][2*j+1]  = s[2*i+1][2*j+1] + (c1*hyy + c2*hxx) * gg[i][j];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2_hourglass */
#endif
/*! @} (documentation module close)*/
