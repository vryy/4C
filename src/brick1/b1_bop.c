#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void b1_bop(double    **b,
            double    **deriv,
            double    **xjm,
            double      det,
            int         iel)
{
/*----------------------------------------------------------------------*/
int i,node_start;
double dum;
double x1r, x2r, x3r, x1s, x2s, x3s, x1t, x2t, x3t;
double xi11, xi12, xi13, xi21, xi22, xi23, xi31, xi32, xi33;
double hr, hs, ht;
double h1, h2, h3;
#ifdef DEBUG 
dstrc_enter("b1_bop");
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
} /* end of b1_bop */
/*----------------------------------------------------------------------*/
