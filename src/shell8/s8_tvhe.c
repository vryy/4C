/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  calculates metrics (geom. nonlinear)                  m.gee 6/01    |
 *----------------------------------------------------------------------*/
/*
C.......................................................................
C!    METRIK DES SCHALENRAUMS 
c!    gmkovc wird hier erneut berechnet, um die Vernachlaessigung der
c!    in e3 quadratischen Anteile zu beruecksichtigen, d.h.
c!    gmkovc_ij ungleich gkovc_i*gkovc_j                                         .
C.......................................................................
*/
void s8_tvhe(DOUBLE **gmkovr,
             DOUBLE **gmkovc,
             DOUBLE **gmkonr,
             DOUBLE **gmkonc,
             DOUBLE **gkovr,
             DOUBLE **gkovc,
             DOUBLE  *detr,
             DOUBLE  *detc,
             DOUBLE **amkovc,
             DOUBLE **amkovr,
             DOUBLE **akovc,
             DOUBLE **akovr,
             DOUBLE **a3kvpc,
             DOUBLE **a3kvpr,
             DOUBLE   e3,
             DOUBLE   condfac)
{
INT i;
DOUBLE b11c=0.0;
DOUBLE b12c=0.0;
DOUBLE b21c=0.0;
DOUBLE b22c=0.0;
DOUBLE b31c=0.0;
DOUBLE b32c=0.0;

DOUBLE b11r=0.0;
DOUBLE b12r=0.0;
DOUBLE b21r=0.0;
DOUBLE b22r=0.0;
DOUBLE b31r=0.0;
DOUBLE b32r=0.0;

DOUBLE det_dummy;

DOUBLE zeta;

#ifdef DEBUG 
dstrc_enter("s8_tvhe");
#endif
/*----------------------------------------------------------------------*/
zeta = e3 / condfac;
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++) b11c += akovc[i][0]*a3kvpc[i][0];
for (i=0; i<3; i++) b12c += akovc[i][0]*a3kvpc[i][1];
for (i=0; i<3; i++) b21c += akovc[i][1]*a3kvpc[i][0];
for (i=0; i<3; i++) b22c += akovc[i][1]*a3kvpc[i][1];
for (i=0; i<3; i++) b31c += akovc[i][2]*a3kvpc[i][0];
for (i=0; i<3; i++) b32c += akovc[i][2]*a3kvpc[i][1];

for (i=0; i<3; i++) b11r += akovr[i][0]*a3kvpr[i][0];
for (i=0; i<3; i++) b12r += akovr[i][0]*a3kvpr[i][1];
for (i=0; i<3; i++) b21r += akovr[i][1]*a3kvpr[i][0];
for (i=0; i<3; i++) b22r += akovr[i][1]*a3kvpr[i][1];
for (i=0; i<3; i++) b31r += akovr[i][2]*a3kvpr[i][0];
for (i=0; i<3; i++) b32r += akovr[i][2]*a3kvpr[i][1];

gmkovc[0][0] = gmkovr[0][0] + (amkovc[0][0]-amkovr[0][0]) + zeta*2.0*(b11c-b11r);
gmkovc[1][1] = gmkovr[1][1] + (amkovc[1][1]-amkovr[1][1]) + zeta*2.0*(b22c-b22r);
gmkovc[2][2] = gmkovr[2][2] + (amkovc[2][2]-amkovr[2][2]);
gmkovc[0][1] = gmkovr[0][1] + (amkovc[0][1]-amkovr[0][1]) + zeta*(b21c+b12c-b21r-b12r);
gmkovc[0][2] = gmkovr[0][2] + (amkovc[0][2]-amkovr[0][2]) + zeta*(b31c-b31r);
gmkovc[1][2] = gmkovr[1][2] + (amkovc[1][2]-amkovr[1][2]) + zeta*(b32c-b32r);
gmkovc[2][0] = gmkovc[0][2]; 
gmkovc[2][1] = gmkovc[1][2]; 
gmkovc[1][0] = gmkovc[0][1]; 
/*----------------------------------------------------------------------*/
math_array_copy(gmkovr,3,3,gmkonr);
math_inv3(gmkonr,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detr = sqrt(det_dummy);

math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tvhe */



/*----------------------------------------------------------------------*
 |  calculates metrics (geom. linear)                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
/*
C.......................................................................
C!    METRIK DES SCHALENRAUMS 
c!    gmkovc wird hier erneut berechnet, um die Vernachlaessigung der
c!    in e3 quadratischen Anteile zu beruecksichtigen, d.h.
c!    gmkovc_ij ungleich gkovc_i*gkovc_j                                         .
C.......................................................................
*/
void s8_tvhe_lin(DOUBLE **gmkovr,
                 DOUBLE **gmkovc,
                 DOUBLE **gmkonr,
                 DOUBLE **gmkonc,
                 DOUBLE **gkovr,
                 DOUBLE **gkovc,
                 DOUBLE  *detr,
                 DOUBLE  *detc,
                 DOUBLE **amkovc,
                 DOUBLE **amkovr,
                 DOUBLE **akovc,
                 DOUBLE **akovr,
                 DOUBLE **a3kvpc,
                 DOUBLE **a3kvpr,
                 DOUBLE   e3)
{
INT i,j,k;
DOUBLE b11c=0.0;
DOUBLE b12c=0.0;
DOUBLE b21c=0.0;
DOUBLE b22c=0.0;
DOUBLE b31c=0.0;
DOUBLE b32c=0.0;

DOUBLE b11r=0.0;
DOUBLE b12r=0.0;
DOUBLE b21r=0.0;
DOUBLE b22r=0.0;
DOUBLE b31r=0.0;
DOUBLE b32r=0.0;

DOUBLE heps[3][3];
DOUBLE det_dummy;


#ifdef DEBUG 
dstrc_enter("s8_tvhe_lin");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   heps[i][j]=0.0;
   for (k=0; k<3; k++)  heps[i][j] += gkovc[k][i] * gkovr[k][j];
}
/*----------------------------------------------------------------------*/
gmkovc[0][0] = 2.0 * heps[0][0] - gmkovr[0][0];
gmkovc[1][1] = 2.0 * heps[1][1] - gmkovr[1][1];
gmkovc[2][2] = 2.0 * heps[2][2] - gmkovr[2][2];

gmkovc[0][1] = heps[0][1]+heps[1][0] - gmkovr[0][1];
gmkovc[0][2] = heps[0][2]+heps[2][0] - gmkovr[0][2];
gmkovc[1][2] = heps[1][2]+heps[2][1] - gmkovr[1][2];

gmkovc[1][0] = gmkovc[0][1];
gmkovc[2][0] = gmkovc[0][2];
gmkovc[2][1] = gmkovc[1][2];
/*----------------------------------------------------------------------*/
math_array_copy(gmkovr,3,3,gmkonr);
math_inv3(gmkonr,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detr = sqrt(det_dummy);

math_array_copy(gmkovc,3,3,gmkonc);
math_inv3(gmkonc,&det_dummy);
if (det_dummy <= 0.0) det_dummy = EPS8;
*detc = sqrt(det_dummy);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tvhe_lin */
#endif
