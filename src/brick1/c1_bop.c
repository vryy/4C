/*!----------------------------------------------------------------------
\file
\brief contains the routines which calculate the operator 
       matrix for a 3D hex element

*----------------------------------------------------------------------*/
#ifdef D_BRICK1

#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculate operator matrix at point r,s,t                                            

<pre>                                                              al 06/01
This routine calcuates the operator matrix b at the given point r,s,t
for an 3D-hex-element.
</pre>
\param bop       DOUBLE**  (o)   the calculated operator matrix
\param bn        DOUBLE**  (o)   the calculated operator without zero values ?
\param deriv     DOUBLE**  (i)   the derivatives of the shape functions
\param xjm       DOUBLE**  (i)   the Jacobian matrix
\param det       DOUBLE    (i)   the determinant of the Jacobian matrix
\param iel       INT       (i)   number of nodes per element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_bop(DOUBLE    **b,
            DOUBLE    **bn,
            DOUBLE    **deriv,
            DOUBLE    **xjm,
            DOUBLE      det,
            INT         iel)
{
/*----------------------------------------------------------------------*/
INT i,node_start;
DOUBLE dum;
DOUBLE x1r, x2r, x3r, x1s, x2s, x3s, x1t, x2t, x3t;
DOUBLE xi11, xi12, xi13, xi21, xi22, xi23, xi31, xi32, xi33;
DOUBLE hr, hs, ht;
DOUBLE h1, h2, h3;
#ifdef DEBUG 
dstrc_enter("c1_bop");
#endif
/*---------------------------------------------- inverse of jacobian ---*/
  x1r = xjm[0][0];
  x2r = xjm[0][1];
  x3r = xjm[0][2];
  x1s = xjm[1][0];
  x2s = xjm[1][1];
  x3s = xjm[1][2];
  x1t = xjm[2][0];
  x2t = xjm[2][1];
  x3t = xjm[2][2];
 
  dum=1.0/det;

  xi11=dum*(x2s*x3t - x2t*x3s);
  xi12=dum*(x3r*x2t - x2r*x3t);
  xi13=dum*(x2r*x3s - x3r*x2s);
  xi21=dum*(x3s*x1t - x3t*x1s);
  xi22=dum*(x1r*x3t - x3r*x1t);
  xi23=dum*(x3r*x1s - x1r*x3s);
  xi31=dum*(x1s*x2t - x1t*x2s);
  xi32=dum*(x2r*x1t - x1r*x2t);
  xi33=dum*(x1r*x2s - x2r*x1s);
/*----------------------------- get operator b of global derivatives ---*/
  for (i=0; i<iel; i++)         
  {
    node_start = i*3;
  
    hr   = deriv[0][i];
    hs   = deriv[1][i];
    ht   = deriv[2][i];

    h1 = xi11*hr + xi12*hs + xi13*ht;
    h2 = xi21*hr + xi22*hs + xi23*ht;
    h3 = xi31*hr + xi32*hs + xi33*ht;
    
    bn[0][i] = h1;
    bn[1][i] = h2;
    bn[2][i] = h3;

    b[0][node_start+0] = h1 ;
    b[0][node_start+1] = 0.0;
    b[0][node_start+2] = 0.0;
    b[1][node_start+0] = 0.0;
    b[1][node_start+1] = h2 ;
    b[1][node_start+2] = 0.0;
    b[2][node_start+0] = 0.0;
    b[2][node_start+1] = 0.0;
    b[2][node_start+2] = h3 ;
    b[3][node_start+0] = h2 ;
    b[3][node_start+1] = h1 ;
    b[3][node_start+2] = 0.0;
    b[4][node_start+0] = 0.0;
    b[4][node_start+1] = h3 ;
    b[4][node_start+2] = h2 ;
    b[5][node_start+0] = h3 ;
    b[5][node_start+1] = 0.0;
    b[5][node_start+2] = h1 ;
  } /* end of loop over nodes */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of c1_bop */

/*!----------------------------------------------------------------------
\brief include initial displacements to derivative operator                                            

<pre>                                                              al 06/01
This routine includes initial displacements to derivative operator
for an 3D-hex-element.
</pre>
\param bop       DOUBLE**  (o)   the calculated operator matrix
\param disd      DOUBLE*   (i)   displacement derivatives
\param iel       INT       (i)   number of nodes per element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_cint()

*----------------------------------------------------------------------*/
void c1_bdis( DOUBLE   **bop, /* b-operator matrix                      */
              DOUBLE   *disd, /* displacement derivatives               */
              INT        iel) /* number of nodes at actual element      */
{
/*----------------------------------------------------------------------*/
INT node_start, inode;
DOUBLE rl11, rl12, rl13, rl21, rl22, rl23, rl31, rl32, rl33;
DOUBLE h1, h2, h3;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_bdis");
#endif
/*----------------------------------------------------------------------*/
  rl11 = disd[0];
  rl12 = disd[3];
  rl13 = disd[7];
  rl21 = disd[4];
  rl22 = disd[1];
  rl23 = disd[5];
  rl31 = disd[8];
  rl32 = disd[6];
  rl33 = disd[2];

  for (inode=0; inode<iel; inode++)
  {
    node_start = inode*3;

   h1 = bop[0][node_start+0];
   h2 = bop[1][node_start+1];
   h3 = bop[2][node_start+2];
   bop[0][node_start+0]  += rl11*h1;
   bop[0][node_start+1]  += rl21*h1;
   bop[0][node_start+2]  += rl31*h1;
   bop[1][node_start+0]  += rl12*h2;
   bop[1][node_start+1]  += rl22*h2;
   bop[1][node_start+2]  += rl32*h2;
   bop[2][node_start+0]  += rl13*h3;
   bop[2][node_start+1]  += rl23*h3;
   bop[2][node_start+2]  += rl33*h3;
   bop[3][node_start+0]  += rl11*h2 + rl12*h1;
   bop[3][node_start+1]  += rl21*h2 + rl22*h1;
   bop[3][node_start+2]  += rl31*h2 + rl32*h1;
   bop[4][node_start+0]  += rl12*h3 + rl13*h2;
   bop[4][node_start+1]  += rl22*h3 + rl23*h2;
   bop[4][node_start+2]  += rl32*h3 + rl33*h2;
   bop[5][node_start+0]  += rl11*h3 + rl13*h1;
   bop[5][node_start+1]  += rl21*h3 + rl23*h1;
   bop[5][node_start+2]  += rl31*h3 + rl33*h1;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of c1_bdis */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
