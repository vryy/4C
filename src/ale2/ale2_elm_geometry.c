/*!----------------------------------------------------------------------
\file
\brief contains the routines to serve the element stiffness matrices

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

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  area linear quadrilaterals

<pre>                                                              ck 01/03
This routine evaluates the area of a 4 noded quad ale element

</pre>
\param **xyz    DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y

\warning There is nothing special to this routine
\return DOUBLE element area
\sa calling:
             called by: ale2_static_ke_prestress(),
	                ale2_static_ke_step2()

*----------------------------------------------------------------------*/
DOUBLE ale2_el_area(DOUBLE **xyz)
{
DOUBLE a, b, c;  /* geometrical values */
DOUBLE el_area;  /* element area */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_el_area");
#endif
/*----------------------------------------------------------------------*/
a = (xyz[0][0]-xyz[1][0])*(xyz[0][0]-xyz[1][0])
   +(xyz[0][1]-xyz[1][1])*(xyz[0][1]-xyz[1][1]); /* line 0-1 squared */
b = (xyz[1][0]-xyz[2][0])*(xyz[1][0]-xyz[2][0])
   +(xyz[1][1]-xyz[2][1])*(xyz[1][1]-xyz[2][1]); /* line 1-2 squared */
c = (xyz[2][0]-xyz[0][0])*(xyz[2][0]-xyz[0][0])
   +(xyz[2][1]-xyz[0][1])*(xyz[2][1]-xyz[0][1]); /* diag 2-0 squared */
el_area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
a = (xyz[2][0]-xyz[3][0])*(xyz[2][0]-xyz[3][0])
   +(xyz[2][1]-xyz[3][1])*(xyz[2][1]-xyz[3][1]); /* line 2-3 squared */
b = (xyz[3][0]-xyz[0][0])*(xyz[3][0]-xyz[0][0])
   +(xyz[3][1]-xyz[0][1])*(xyz[3][1]-xyz[0][1]); /* line 3-0 squared */
/*-------------------------------------------- evaluate element area ---*/
el_area += 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return el_area;
}


/*!----------------------------------------------------------------------
\brief  area linear triangles

<pre>                                                              ck 01/03
This routine evaluates the area of a 3 noded triangular ale element

</pre>
\param **xyz    DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y
\param   i      INT       (i)   node i
\param   j      INT       (i)   node j
\param   k      INT       (i)   node k

\warning There is nothing special to this routine
\return DOUBLE element area
\sa calling:
             called by: ale2_torsional()

*----------------------------------------------------------------------*/
DOUBLE ale2_area_tria(DOUBLE **xyz, INT i, INT j, INT k)
{
DOUBLE a, b, c;  /* geometrical values */
DOUBLE el_area;  /* element area */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_area_tria");
#endif
/*----------------------------------------------------------------------*/
a = (xyz[i][0]-xyz[j][0])*(xyz[i][0]-xyz[j][0])
   +(xyz[i][1]-xyz[j][1])*(xyz[i][1]-xyz[j][1]); /* line i-j squared */
b = (xyz[j][0]-xyz[k][0])*(xyz[j][0]-xyz[k][0])
   +(xyz[j][1]-xyz[k][1])*(xyz[j][1]-xyz[k][1]); /* line j-k squared */
c = (xyz[k][0]-xyz[i][0])*(xyz[k][0]-xyz[i][0])
   +(xyz[k][1]-xyz[i][1])*(xyz[k][1]-xyz[i][1]); /* line k-i squared */
el_area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
#ifdef DEBUG
dstrc_exit();
#endif

return el_area;
}



/*!----------------------------------------------------------------------
\brief  determine the length and direction of a line

<pre>                                                             ck 06/03
This routine evaluates the length and direction (sin and cos) of a line
between two points i and j.

</pre>
\param **xyz    DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y
\param   i      INT       (i)   node i
\param   j      INT       (i)   node j
\param   k      INT       (i)   node k

\warning There is nothing special to this routine
\return void
\sa calling:
             called by: ale2_static_ke_spring()

*----------------------------------------------------------------------*/
void edge_geometry(INT      i,
                   INT      j,
         	   DOUBLE **xyz,
	           DOUBLE  *length,
	           DOUBLE  *sin_alpha,
	           DOUBLE  *cos_alpha)
{
DOUBLE   delta_x, delta_y;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("edge_length");
#endif
/*---------------------------------------------- x- and y-difference ---*/
delta_x = xyz[j][0]-xyz[i][0];
delta_y = xyz[j][1]-xyz[i][1];
/*------------------------------- determine distance between i and j ---*/
*length = sqrt( delta_x * delta_x
              + delta_y * delta_y );
if (*length < EPS14) dserror("edge or diagonal of element has zero length");
/*--------------------------------------- determine direction of i-j ---*/
*sin_alpha = delta_y / *length;
*cos_alpha = delta_x / *length;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}

/*!----------------------------------------------------------------------
\brief determines the torsional spring stiffness of a triangle i, j, k

<pre>                                                             ck 06/03
This routine evaluates the torsional springs attached to the vertices of
a triangle i, j, k following Farhat et al. in 'Torsional springs for
two-dimensional dynamic unstructured fluid meshes' in Comput. Mech.
Engrg. 163 (1998) 231-245
It determines the 6 x 6 matrix which corresponds to the displacement
degrees of freedom u and v of the nodes.

</pre>
\param   i         INT       (i)   node i
\param   j         INT       (i)   node j
\param   k         INT       (i)   node k
\param **xyz       DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y
\param **k_torsion DOUBLE    (o)   torsional stiffness matrix (6x6)
\param   init      INT       (i)   initialisation flag

\warning There is nothing special to this routine
\return void
\sa calling:
             called by: ale2_tors_spring_quad4(),
	                ale2_tors_spring_tri3()

*----------------------------------------------------------------------*/
void ale2_torsional(INT      i,
                    INT      j,
		    INT      k,
		    DOUBLE **xyz,
		    DOUBLE **k_torsion,
		    INT      init)
/*
                          k
                           *
			  / \ l_jk
    y,v ^	   l_ki  /   \
	|		/     \
	--->	     i *-------* j
	  x,u	          l_ij
*/
{
DOUBLE x_ij, x_jk, x_ki;  /* x-differences between nodes */
DOUBLE y_ij, y_jk, y_ki;  /* y-differences between nodes */
DOUBLE l_ij, l_jk, l_ki;  /* side lengths */
DOUBLE a_ij, a_jk, a_ki;  /* auxiliary values same as in Farhat et al. */
DOUBLE b_ij, b_jk, b_ki;  /*                  - " -                    */
DOUBLE area;              /* area of the triangle */
static ARRAY    R_a;      /* rotation matrix same as in Farhat et al. */
static DOUBLE **R;
static ARRAY    C_a;      /* torsion stiffness matrix same as in Farhat et al. */
static DOUBLE **C;
static ARRAY    A_a;      /* auxiliary array of intermediate results */
static DOUBLE **A;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_torsional");
#endif
/*--------------------------------------------------- initialisation ---*/
if (init==1)
{
   R = amdef("R", &R_a, 3, 6, "DA");
   C = amdef("C", &C_a, 3, 3, "DA");
   A = amdef("A", &A_a, 6, 3, "DA");
   goto end;
}
/*----------------------------------------------------------------------*/
amzero(&C_a);
/*--------------------------------- determine basic geometric values ---*/
x_ij = xyz[j][0] - xyz[i][0];
x_jk = xyz[k][0] - xyz[j][0];
x_ki = xyz[i][0] - xyz[k][0];
y_ij = xyz[j][1] - xyz[i][1];
y_jk = xyz[k][1] - xyz[j][1];
y_ki = xyz[i][1] - xyz[k][1];

l_ij = sqrt( x_ij*x_ij + y_ij*y_ij );
l_jk = sqrt( x_jk*x_jk + y_jk*y_jk );
l_ki = sqrt( x_ki*x_ki + y_ki*y_ki );
/*----------------------------------------------- check edge lengths ---*/
if (l_ij < EPS14) dserror("edge or diagonal of element has zero length");
if (l_jk < EPS14) dserror("edge or diagonal of element has zero length");
if (l_ki < EPS14) dserror("edge or diagonal of element has zero length");
/*-------------------------------------------- fill auxiliary values ---*/
a_ij = x_ij / (l_ij*l_ij);
a_jk = x_jk / (l_jk*l_jk);
a_ki = x_ki / (l_ki*l_ki);
b_ij = y_ij / (l_ij*l_ij);
b_jk = y_jk / (l_jk*l_jk);
b_ki = y_ki / (l_ki*l_ki);
/*--------------------------------------------------- determine area ---*/
area = ale2_area_tria(xyz,i,j,k);
/*---------------------------------- determine torsional stiffnesses ---*/
C[0][0] = l_ij*l_ij * l_ki*l_ki / (4.0*area*area);
C[1][1] = l_ij*l_ij * l_jk*l_jk / (4.0*area*area);
C[2][2] = l_ki*l_ki * l_jk*l_jk / (4.0*area*area);
/*--------------------------------------- fill transformation matrix ---*/
R[0][0] = - b_ki - b_ij;
R[0][1] = a_ij + a_ki;
R[0][2] = b_ij;
R[0][3] = - a_ij;
R[0][4] = b_ki;
R[0][5] = - a_ki;

R[1][0] = b_ij;
R[1][1] = - a_ij;
R[1][2] = - b_ij - b_jk;
R[1][3] = a_jk + a_ij;
R[1][4] = b_jk;
R[1][5] = - a_jk;

R[2][0] = b_ki;
R[2][1] = - a_ki;
R[2][2] = b_jk;
R[2][3] = - a_jk;
R[2][4] = - b_jk - b_ki;
R[2][5] = a_ki + a_jk;
/*----------------------------------- perform matrix multiplications ---*/
math_mattrnmatdense(A,R,C,6,3,3,0,1.0);         /* A = R^t * C */
math_matmatdense(k_torsion,A,R,6,3,6,0,1.0);   /* stiff = A * R */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}

/*!----------------------------------------------------------------------
\brief adds the torsional spring stiffness part to a quad 4 ale element

<pre>                                                             ck 06/03
divides a quad 4 ale element into 4 subtriangles which contain each three
of the four elemental nodes, determines the torsional springs for these
subtriangles and adds these values to the appropriate places in the elemental
stiffness matrix.

</pre>
\param **estif     DOUBLE    (o)   element stiffness matrix
\param **xyz       DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y
\param   init      INT       (i)   initialisation flag

\warning There is nothing special to this routine
\return void
\sa calling: ale2_torsional()
             called by: ale2_static_ke_spring()

*----------------------------------------------------------------------*/
void ale2_tors_spring_quad4(DOUBLE **estif, DOUBLE **xyz, INT init)
{
INT i, j;                 /* counters */

static ARRAY    k_tria_a; /* rotational stiffness matrix of one triangle */
static DOUBLE **k_tria;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_tors_spring_quad4");
#endif
/*--------------------------------------------------- initialisation ---*/
if (init==1)
{
   k_tria = amdef("k_tria", &k_tria_a, 6, 6, "DA");
   ale2_torsional(0,0,0,NULL,NULL,1);
   goto end;
}
/*--- pass all nodes and determine the triangle defined by the node i and
adjunct nodes... ---*/
ale2_torsional(0,1,2,xyz,k_tria,0);
/*---------- ...sort everything into the element stiffness matrix... ---*/
for (i=0; i<6; i++)
   for (j=0; j<6; j++)
      estif[i][j] += k_tria[i][j];
/*--------------------------------- ...repeat for second triangle... ---*/
ale2_torsional(1,2,3,xyz,k_tria,0);
for (i=2; i<8; i++)
   for (j=2; j<8; j++)
      estif[i][j] += k_tria[i-2][j-2];
/*------------------------------------------ ...and for the third... ---*/
ale2_torsional(2,3,0,xyz,k_tria,0);
for (i=4; i<8; i++)
   for (j=4; j<8; j++)
      estif[i][j] += k_tria[i-4][j-4];
for (i=0; i<2; i++)
   for (j=0; j<2; j++)
      estif[i][j] += k_tria[i+4][j+4];
for (i=4; i<8; i++)
{
   for (j=0; j<2; j++)
   {
      estif[i][j] += k_tria[i-4][j+4];
      estif[j][i] += k_tria[i-4][j+4];
   }
}
/*------------------------------- ...and eventually for a forth time ---*/
ale2_torsional(3,0,1,xyz,k_tria,0);
for (i=0; i<4; i++)
   for (j=0; j<4; j++)
      estif[i][j] += k_tria[i+2][j+2];
for (i=6; i<8; i++)
   for (j=6; j<8; j++)
      estif[i][j] += k_tria[i-6][j-6];
for (i=6; i<8; i++)
{
   for (j=0; j<4; j++)
   {
      estif[i][j] += k_tria[i-6][j+2];
      estif[j][i] += k_tria[i-6][j+2];
   }
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
}



/*!----------------------------------------------------------------------
\brief adds the torsional spring stiffness part to a quad 4 ale element

<pre>                                                             ck 06/03
adds the torisonal stiffness part to the elemental stiffness of a
linear triangle ale element

</pre>

\param **estif     DOUBLE    (o)   element stiffness matrix
\param **xyz       DOUBLE    (i)   elemental coordinates 1st index: node
                                                      2dn index: x or y
\param   init      INT       (i)   initialisation flag

\warning There is nothing special to this routine
\return void
\sa calling: ale2_torsional()
             called by: ale2_static_ke_spring()

*----------------------------------------------------------------------*/
void ale2_tors_spring_tri3(DOUBLE **estif, DOUBLE **xyz, INT init)
{
INT i, j;      /* counters */

static ARRAY    k_tria_a; /* rotational stiffness matrix of one triangle */
static DOUBLE **k_tria;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("ale2_tors_spring_tri3");
#endif
/*--------------------------------------------------- initialisation ---*/
if (init==1)
{
   k_tria = amdef("k_tria", &k_tria_a, 6, 6, "DA");
   ale2_torsional(0,0,0,NULL,NULL,1);
   goto end;
}
/*-------------------------------- evaluate torsional stiffness part ---*/
ale2_torsional(0,1,2,xyz,k_tria,0);
/*-------------------------- add everything to the element stiffness ---*/
for (i=0; i<6; i++)
   for (j=0; j<6; j++)
      estif[i][j] += k_tria[i][j];
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return; /* end of ale2_tors_spring_tri3 */
}

/*!----------------------------------------------------------------------
\brief calculate global derivatives at point r,s

<pre>                                                            ck 06/03
This routine calcuates global derivatives of the shape functions at the
given point r,s for an 2D-ale-element.

</pre>
\param **deriv_xy  DOUBLE  (o)   global derivatives of the shape functions
\param **deriv     DOUBLE  (i)   the derivatives of the shape functions
\param **xjm       DOUBLE  (i)   the Jacobian matrix
\param det         DOUBLE  (i)   the determinant of the Jacobian matrix
\param iel         INT     (i)   number of nodes per element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: ale2_static_ke_laplace()

*----------------------------------------------------------------------*/
void ale2_deriv_xy(DOUBLE    **deriv_xy,
                   DOUBLE    **deriv,
                   DOUBLE    **xjm,
                   DOUBLE      det,
                   INT         iel)
{
/*----------------------------------------------------------------------*/
INT inode;
DOUBLE dum;
DOUBLE xji[2][2];
#ifdef DEBUG
dstrc_enter("ale2_bop");
#endif
/*---------------------------------------------- inverse of jacobian ---*/
dum = 1.0/det;
xji[0][0] = xjm[1][1]* dum;
xji[0][1] =-xjm[0][1]* dum;
xji[1][0] =-xjm[1][0]* dum;
xji[1][1] = xjm[0][0]* dum;
/*------------------------------------------- get global derivatives ---*/
for (inode=0; inode<iel; inode++)
{
  deriv_xy[0][inode] += xji[0][0] * deriv[0][inode]
                     +  xji[0][1] * deriv[1][inode];
  deriv_xy[1][inode] += xji[1][0] * deriv[0][inode]
                     +  xji[1][1] * deriv[1][inode];
} /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of ale2_deriv_xy */
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
#endif
