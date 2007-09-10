/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | collocation points, shape functions and derivatives    m.gee 2/02    |
 | and all necessary basis vectors and metrics at these points          |
 *----------------------------------------------------------------------*/
void s8_ans_colloqpoints(INT nsansq,INT iel,INT ans,DIS_TYP distyp,
                        DOUBLE xr1[],DOUBLE xs1[],DOUBLE xr2[],DOUBLE xs2[],
                        DOUBLE *funct1q[],DOUBLE **deriv1q[],
                        DOUBLE *funct2q[],DOUBLE **deriv2q[],
                        DOUBLE **xrefe,DOUBLE **a3r,DOUBLE **xcure,DOUBLE **a3c,
                        DOUBLE **akovr1q[] ,DOUBLE **akonr1q[],
                        DOUBLE **amkovr1q[],DOUBLE **amkonr1q[],
                        DOUBLE **a3kvpr1q[],
                        DOUBLE **akovc1q[] ,DOUBLE **akonc1q[],
                        DOUBLE **amkovc1q[],DOUBLE **amkonc1q[],
                        DOUBLE **a3kvpc1q[],
                        DOUBLE **akovr2q[] ,DOUBLE **akonr2q[],
                        DOUBLE **amkovr2q[],DOUBLE **amkonr2q[],
                        DOUBLE **a3kvpr2q[],
                        DOUBLE **akovc2q[] ,DOUBLE **akonc2q[],
                        DOUBLE **amkovc2q[],DOUBLE **amkonc2q[],
                        DOUBLE **a3kvpc2q[],
                        DOUBLE *detr, DOUBLE *detc)
{
INT i;
#ifdef DEBUG
dstrc_enter("s8_ans_colloqpoints");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------- get coordinates of collocation points */
s8_ans_colloqcoords(xr1,xs1,xr2,xs2,iel,ans);
/*--------------------------- loop over all collocation points (2 or 6) */
for (i=0; i<nsansq; i++)
{
  s8_funct_deriv(funct1q[i],deriv1q[i],xr1[i],xs1[i],distyp,1);

  s8_tvmr(xrefe,a3r,akovr1q[i],akonr1q[i],amkovr1q[i],amkonr1q[i],detr,
          funct1q[i],deriv1q[i],iel,a3kvpr1q[i],0);

  s8_tvmr(xcure,a3c,akovc1q[i],akonc1q[i],amkovc1q[i],amkonc1q[i],detc,
          funct1q[i],deriv1q[i],iel,a3kvpc1q[i],0);


  s8_funct_deriv(funct2q[i],deriv2q[i],xr2[i],xs2[i],distyp,1);

  s8_tvmr(xrefe,a3r,akovr2q[i],akonr2q[i],amkovr2q[i],amkonr2q[i],detr,
          funct2q[i],deriv2q[i],iel,a3kvpr2q[i],0);

  s8_tvmr(xcure,a3c,akovc2q[i],akonc2q[i],amkovc2q[i],amkonc2q[i],detc,
          funct2q[i],deriv2q[i],iel,a3kvpc2q[i],0);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_ans_colloqpoints */

/*----------------------------------------------------------------------*
 | collocation points, shape functions and derivatives    m.gee 2/02    |
 *----------------------------------------------------------------------*/
void s8_ans_colloqcoords(DOUBLE xqr1[], DOUBLE xqs1[],
                        DOUBLE xqr2[], DOUBLE xqs2[],
                        INT iel, INT ans)
{
DOUBLE rthreei;
#ifdef DEBUG
dstrc_enter("s8_ans_colloq");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------- ans==1 -> ans for querschublocking */
if (ans==1)
{
   if (iel==4)/* 4-noded element */
   {
      xqr1[0] =  0.0;     xqs1[0] = -1.0; /* point ( 0.0/-1.0) */
      xqr1[1] =  0.0;     xqs1[1] =  1.0; /* point ( 0.0/ 1.0) */

      xqr2[0] = -1.0;     xqs2[0] =  0.0; /* point (-1.0/ 0.0) */
      xqr2[1] =  1.0;     xqs2[1] =  0.0; /* point ( 1.0/ 0.0) */
   }
   else if (iel==9)/* 9-noded element */
   {
      rthreei = 1.0 / (sqrt(3.0));
      xqr1[0] = -rthreei; xqs1[0] = -1.0;
      xqr1[1] = -rthreei; xqs1[1] =  0.0;
      xqr1[2] = -rthreei; xqs1[2] =  1.0;
      xqr1[3] =  rthreei; xqs1[3] = -1.0;
      xqr1[4] =  rthreei; xqs1[4] =  0.0;
      xqr1[5] =  rthreei; xqs1[5] =  1.0;

      xqr2[0] = -1.0;     xqs2[0] = -rthreei;
      xqr2[1] = 0.0;     xqs2[1] =  -rthreei;
      xqr2[2] =  1.0;     xqs2[2] = -rthreei;
      xqr2[3] = -1.0;     xqs2[3] =  rthreei;
      xqr2[4] =  0.0;     xqs2[4] =  rthreei;
      xqr2[5] =  1.0;     xqs2[5] =  rthreei;
   }
} /* end of (ans==1) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_ans_colloq */


/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 | shape functions for assumed strains (querschub)                      |
 *----------------------------------------------------------------------*/
void s8_ansq_funct(DOUBLE frq[], DOUBLE fsq[], DOUBLE r, DOUBLE s,
                  INT iel, INT nsansq)
{
DOUBLE rthreei;
DOUBLE pr[3],ps[3];
DOUBLE qr[2],qs[2];
DOUBLE rr[3],rs[3];
#ifdef DEBUG
dstrc_enter("s8_ansq_funct");
#endif
/*----------------------------------------------------------------------*/
if (iel==4)
{
   frq[0] = 0.5 * (1.0 - s);
   frq[1] = 0.5 * (1.0 + s);

   fsq[0] = 0.5 * (1.0 - r);
   fsq[1] = 0.5 * (1.0 + r);
}
if (iel==9)
{
   rthreei = 1.0 / (sqrt(3.0));

   pr[0] = -0.5 * s * (1.0-s);
   pr[1] =  (1.0-s) * (1.0+s);
   pr[2] =  0.5 * s * (1.0+s);

   qr[0] =  0.5 * (1.0-r/rthreei);
   qr[1] =  0.5 * (1.0+r/rthreei);

   rr[0] =  1.0/6.0 - 0.5 * s;
   rr[1] =  2.0/3.0;
   rr[2] =  1.0/6.0 + 0.5 * s;

   ps[0] = -0.5 * r * (1.0-r);
   ps[1] =  (1.0-r) * (1.0+r);
   ps[2] =  0.5 * r * (1.0+r);

   qs[0] =  0.5 * (1.0-s/rthreei);
   qs[1] =  0.5 * (1.0+s/rthreei);

   rs[0] =  1.0/6.0 - 0.5 * r;
   rs[1] =  2.0/3.0;
   rs[2] =  1.0/6.0 + 0.5 * r;

   frq[0] = pr[0] * qr[0];
   frq[1] = pr[1] * qr[0];
   frq[2] = pr[2] * qr[0];
   frq[3] = pr[0] * qr[1];
   frq[4] = pr[1] * qr[1];
   frq[5] = pr[2] * qr[1];

   fsq[0] = ps[0] * qs[0];
   fsq[1] = ps[1] * qs[0];
   fsq[2] = ps[2] * qs[0];
   fsq[3] = ps[0] * qs[1];
   fsq[4] = ps[1] * qs[1];
   fsq[5] = ps[2] * qs[1];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_ansq_funct */


/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 |                                                                      |
 *----------------------------------------------------------------------*/
void s8_ans_bbar_q(DOUBLE **bop, DOUBLE frq[], DOUBLE fsq[],
                  DOUBLE  *funct1q[],  DOUBLE  *funct2q[],
                  DOUBLE **deriv1q[],  DOUBLE **deriv2q[],
                  DOUBLE **akovc1q[],  DOUBLE **akovc2q[],
                  DOUBLE **a3kvpc1q[], DOUBLE **a3kvpc2q[],
                  INT iel, INT numdf, INT nsansq)
{
INT inode, node_start;
INT isamp;
DOUBLE a1x1;
DOUBLE a1y1;
DOUBLE a1z1;
DOUBLE a3x1;
DOUBLE a3y1;
DOUBLE a3z1;
DOUBLE a2x2;
DOUBLE a2y2;
DOUBLE a2z2;
DOUBLE a3x2;
DOUBLE a3y2;
DOUBLE a3z2;
DOUBLE a31x1;
DOUBLE a31y1;
DOUBLE a31z1;
DOUBLE a32x2;
DOUBLE a32y2;
DOUBLE a32z2;
DOUBLE p1k;
DOUBLE p2k;
DOUBLE pk1;
DOUBLE pk2;
DOUBLE fris;
DOUBLE fsis;
#ifdef DEBUG
dstrc_enter("s8_ans_bbar_q");
#endif
/*----------------------------------------------------------------------*/
for (inode=0; inode<iel; inode++)
{
   node_start = inode*numdf;

   bop[2][node_start+0]= 0.0;
   bop[2][node_start+1]= 0.0;
   bop[2][node_start+2]= 0.0;
   bop[2][node_start+3]= 0.0;
   bop[2][node_start+4]= 0.0;
   bop[2][node_start+5]= 0.0;

   bop[4][node_start+0]= 0.0;
   bop[4][node_start+1]= 0.0;
   bop[4][node_start+2]= 0.0;
   bop[4][node_start+3]= 0.0;
   bop[4][node_start+4]= 0.0;
   bop[4][node_start+5]= 0.0;

   for (isamp=0; isamp<nsansq; isamp++)
   {
      a1x1 = akovc1q[isamp][0][0];
      a1y1 = akovc1q[isamp][1][0];
      a1z1 = akovc1q[isamp][2][0];
      a3x1 = akovc1q[isamp][0][2];
      a3y1 = akovc1q[isamp][1][2];
      a3z1 = akovc1q[isamp][2][2];

      a2x2 = akovc2q[isamp][0][1];
      a2y2 = akovc2q[isamp][1][1];
      a2z2 = akovc2q[isamp][2][1];
      a3x2 = akovc2q[isamp][0][2];
      a3y2 = akovc2q[isamp][1][2];
      a3z2 = akovc2q[isamp][2][2];

      a31x1 = a3kvpc1q[isamp][0][0];
      a31y1 = a3kvpc1q[isamp][1][0];
      a31z1 = a3kvpc1q[isamp][2][0];

      a32x2 = a3kvpc2q[isamp][0][1];
      a32y2 = a3kvpc2q[isamp][1][1];
      a32z2 = a3kvpc2q[isamp][2][1];

      p1k = funct1q[isamp][inode];
      p2k = funct2q[isamp][inode];

      pk1 = deriv1q[isamp][0][inode];
      pk2 = deriv2q[isamp][1][inode];

      fris = frq[isamp];
      fsis = fsq[isamp];
/*--------------------------------------------------E13(CONST)-------- */
      bop[2][node_start+0]+= pk1*a3x1*fris;
      bop[2][node_start+1]+= pk1*a3y1*fris;
      bop[2][node_start+2]+= pk1*a3z1*fris;
      bop[2][node_start+3]+= p1k*a1x1*fris;
      bop[2][node_start+4]+= p1k*a1y1*fris;
      bop[2][node_start+5]+= p1k*a1z1*fris;
/*--------------------------------------------------E23(CONST)-------- */
      bop[4][node_start+0]+= pk2*a3x2*fsis;
      bop[4][node_start+1]+= pk2*a3y2*fsis;
      bop[4][node_start+2]+= pk2*a3z2*fsis;
      bop[4][node_start+3]+= p2k*a2x2*fsis;
      bop[4][node_start+4]+= p2k*a2y2*fsis;
      bop[4][node_start+5]+= p2k*a2z2*fsis;

   } /* end of (isamp=0; isamp<nsansq; isamp++) */
} /* end of (inode=0; inode<iel; inode++) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_ans_bbar_q */



/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 | modifications to metrics of shell body due to ans for querschub      |
 *----------------------------------------------------------------------*/
void s8_ans_tvhe_q(DOUBLE **gmkovr,DOUBLE **gmkovc,DOUBLE **gmkonr,DOUBLE **gmkonc,
                  DOUBLE **gkovr,DOUBLE **gkovc,DOUBLE **amkovc,DOUBLE **amkovr,
                  DOUBLE **akovc,DOUBLE **akovr,DOUBLE **a3kvpc,DOUBLE **a3kvpr,
                  DOUBLE *detr,   DOUBLE *detc,
                  DOUBLE **amkovr1q[], DOUBLE **amkovc1q[],
                  DOUBLE **akovr1q[] , DOUBLE **akovc1q[] ,
                  DOUBLE **a3kvpr1q[], DOUBLE **a3kvpc1q[],
                  DOUBLE **amkovr2q[], DOUBLE **amkovc2q[],
                  DOUBLE **akovr2q[] , DOUBLE **akovc2q[] ,
                  DOUBLE **a3kvpr2q[], DOUBLE **a3kvpc2q[],
                  DOUBLE frq[], DOUBLE fsq[], DOUBLE e3, INT nansq, INT iel,
                  DOUBLE condfac)
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
dstrc_enter("s8_ans_tvhe_q");
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
/*----------------------------------------------------------------------*/
gmkovc[0][0] = gmkovr[0][0] + (amkovc[0][0]-amkovr[0][0]) + zeta * 2.0 * (b11c-b11r);
gmkovc[1][1] = gmkovr[1][1] + (amkovc[1][1]-amkovr[1][1]) + zeta * 2.0 * (b22c-b22r);
gmkovc[2][2] = gmkovr[2][2] + (amkovc[2][2]-amkovr[2][2]);
gmkovc[0][1] = gmkovr[0][1] + (amkovc[0][1]-amkovr[0][1]) + zeta * (b21c+b12c-b21r-b12r);
gmkovc[0][2] = gmkovr[0][2]                               + zeta * (b31c-b31r);
gmkovc[1][2] = gmkovr[1][2]                               + zeta * (b32c-b32r);
gmkovc[2][0] = gmkovc[0][2];
gmkovc[2][1] = gmkovc[1][2];
gmkovc[1][0] = gmkovc[0][1];
/*----------------------------------------------------------------------*/
for (i=0; i<nansq; i++)
{
   gmkovc[0][2] += (amkovc1q[i][0][2]-amkovr1q[i][0][2]) * frq[i];
   gmkovc[1][2] += (amkovc2q[i][1][2]-amkovr2q[i][1][2]) * fsq[i];
}
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
} /* end of s8_ans_tvhe_q */



/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 | geometric stiffness matrix kg with ans for querschub                 |
 *----------------------------------------------------------------------*/
void s8_ans_tvkg(DOUBLE **estif,DOUBLE *stress_r,DOUBLE *funct,DOUBLE **deriv,
                INT numdf,INT iel,DOUBLE weight,DOUBLE e1,DOUBLE e2,
                DOUBLE frq[], DOUBLE fsq[],DOUBLE *funct1q[],DOUBLE  *funct2q[],
                DOUBLE **deriv1q[], DOUBLE **deriv2q[],INT ansq, INT nsansq)
{
INT     i,inode,jnode;
INT     i_indiz,j_indiz;
DOUBLE  pi;
DOUBLE  pj;
DOUBLE  d11;
DOUBLE  d12;
DOUBLE  d21;
DOUBLE  d22;

DOUBLE  pd1ij;
DOUBLE  pd1ji;
DOUBLE  pd2ij;
DOUBLE  pd2ji;

DOUBLE  xn;
DOUBLE  xm;
DOUBLE  yu;
DOUBLE  yo;
DOUBLE  yy;
DOUBLE  z;

DOUBLE  sn11;
DOUBLE  sn21;
DOUBLE  sn31;
DOUBLE  sn22;
DOUBLE  sn32;
DOUBLE  sn33;
DOUBLE  sm11;
DOUBLE  sm21;
DOUBLE  sm31;
DOUBLE  sm22;
DOUBLE  sm32;

#ifdef DEBUG
dstrc_enter("s8_ans_tvkg");
#endif
/*----------------------------------------------------------------------*/
sn11 = stress_r[0];
sn21 = stress_r[1];
sn31 = stress_r[2];
sn22 = stress_r[3];
sn32 = stress_r[4];
sn33 = stress_r[5];
sm11 = stress_r[6];
sm21 = stress_r[7];
sm31 = stress_r[8];
sm22 = stress_r[9];
sm32 = stress_r[10];
/*----------------------------------------------------------------------*/
for (inode=0; inode<iel; inode++)
{
   for (jnode=0; jnode<=inode; jnode++)
   {
      pi = funct[inode];
      pj = funct[jnode];

      d11 = deriv[0][inode] * deriv[0][jnode];
      d12 = deriv[0][inode] * deriv[1][jnode];
      d21 = deriv[1][inode] * deriv[0][jnode];
      d22 = deriv[1][inode] * deriv[1][jnode];

      xn = (sn11*d11 + sn21*(d12+d21) + sn22*d22) * weight;
      xm = (sm11*d11 + sm21*(d12+d21) + sm22*d22) * weight;

      /*----------------------------------------- no ans for querschub */
      if (ansq==0)
      {
         pd1ij = deriv[0][inode] * pj;
         pd1ji = deriv[0][jnode] * pi;
         pd2ij = deriv[1][inode] * pj;
         pd2ji = deriv[1][jnode] * pi;
         yu = (sn31*pd1ji + sn32*pd2ji) * weight;
         yo = (sn31*pd1ij + sn32*pd2ij) * weight;
      }
      else/*----------------------------------------- ans for querschub */
      {
         yu=0.0;
         yo=0.0;
         for (i=0; i<nsansq; i++)
         {
            pd1ij = deriv1q[i][0][inode] * funct1q[i][jnode] * frq[i];
            pd1ji = deriv1q[i][0][jnode] * funct1q[i][inode] * frq[i];
            pd2ij = deriv2q[i][1][inode] * funct2q[i][jnode] * fsq[i];
            pd2ji = deriv2q[i][1][jnode] * funct2q[i][inode] * fsq[i];

            yu += (sn31*pd1ji + sn32*pd2ji) * weight;
            yo += (sn31*pd1ij + sn32*pd2ij) * weight;
         }
      }
      /*---------------------------------------------------------------*/

      /*---------------- linear part of querschub is always unmodified */
      pd1ij = deriv[0][inode] * pj;
      pd1ji = deriv[0][jnode] * pi;
      pd2ij = deriv[1][inode] * pj;
      pd2ji = deriv[1][jnode] * pi;

      yy = (sm31*(pd1ij+pd1ji) + sm32*(pd2ij+pd2ji)) * weight;
      z  = pi*pj*sn33*weight;

      i_indiz = inode*numdf;
      j_indiz = jnode*numdf;

      estif[inode*numdf+0][jnode*numdf+0] += xn;
      estif[inode*numdf+1][jnode*numdf+1] += xn;
      estif[inode*numdf+2][jnode*numdf+2] += xn;

      estif[inode*numdf+3][jnode*numdf+0] += (xm+yu);
      estif[inode*numdf+4][jnode*numdf+1] += (xm+yu);
      estif[inode*numdf+5][jnode*numdf+2] += (xm+yu);

      estif[inode*numdf+0][jnode*numdf+3] += (xm+yo);
      estif[inode*numdf+1][jnode*numdf+4] += (xm+yo);
      estif[inode*numdf+2][jnode*numdf+5] += (xm+yo);

      estif[inode*numdf+3][jnode*numdf+3] += (yy+z);
      estif[inode*numdf+4][jnode*numdf+4] += (yy+z);
      estif[inode*numdf+5][jnode*numdf+5] += (yy+z);

      if (inode!=jnode)
      {
         estif[jnode*numdf+0][inode*numdf+0] += xn;
         estif[jnode*numdf+1][inode*numdf+1] += xn;
         estif[jnode*numdf+2][inode*numdf+2] += xn;

         estif[jnode*numdf+0][inode*numdf+3] += (xm+yu);
         estif[jnode*numdf+1][inode*numdf+4] += (xm+yu);
         estif[jnode*numdf+2][inode*numdf+5] += (xm+yu);

         estif[jnode*numdf+3][inode*numdf+0] += (xm+yo);
         estif[jnode*numdf+4][inode*numdf+1] += (xm+yo);
         estif[jnode*numdf+5][inode*numdf+2] += (xm+yo);

         estif[jnode*numdf+3][inode*numdf+3] += (yy+z);
         estif[jnode*numdf+4][inode*numdf+4] += (yy+z);
         estif[jnode*numdf+5][inode*numdf+5] += (yy+z);
      }

   } /* end loop over jnode */
} /* end loop over inode */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8_ans_tvkg */
#endif
#endif
