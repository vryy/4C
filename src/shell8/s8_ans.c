#include "../headers/standardtypes.h"
#include "shell8.h"

/*----------------------------------------------------------------------*
 | collocation points, shape functions and derivatives    m.gee 2/02    |
 | and all necessary basis vectors and metrics at these points          |
 *----------------------------------------------------------------------*/
int s8_ans_colloqpoints(int nsansq,int iel,int ans,DIS_TYP distyp,
                        double xr1[],double xs1[],double xr2[],double xs2[],
                        double *funct1q[],double **deriv1q[], 
                        double *funct2q[],double **deriv2q[],
                        double **xrefe,double **a3r,double **xcure,double **a3c,
                        double **akovr1q[] ,double **akonr1q[],
                        double **amkovr1q[],double **amkonr1q[],
                        double **a3kvpr1q[],
                        double **akovc1q[] ,double **akonc1q[],
                        double **amkovc1q[],double **amkonc1q[],
                        double **a3kvpc1q[],
                        double **akovr2q[] ,double **akonr2q[],
                        double **amkovr2q[],double **amkonr2q[],
                        double **a3kvpr2q[],
                        double **akovc2q[] ,double **akonc2q[],
                        double **amkovc2q[],double **amkonc2q[],
                        double **a3kvpc2q[],
                        double *detr, double *detc)
{
int i;
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
int s8_ans_colloqcoords(double xqr1[], double xqs1[],
                        double xqr2[], double xqs2[],
                        int iel, int ans)
{
double rthreei;
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
   if (iel==9)/* 9-noded element */
   {
      rthreei = 1.0 / (sqrt(3.0));
      xqr1[0] = -rthreei; xqs1[0] = -1.0; 
      xqr1[1] =  rthreei; xqs1[1] = -1.0; 
      xqr1[2] = -rthreei; xqs1[2] =  0.0; 
      xqr1[3] =  rthreei; xqs1[3] =  0.0; 
      xqr1[4] = -rthreei; xqs1[4] =  1.0; 
      xqr1[5] =  rthreei; xqs1[5] =  1.0; 

      xqr2[0] = -1.0;     xqs2[0] = -rthreei; 
      xqr2[1] = -1.0;     xqs2[1] =  rthreei; 
      xqr2[2] =  0.0;     xqs2[2] = -rthreei; 
      xqr2[3] =  0.0;     xqs2[3] =  rthreei; 
      xqr2[4] =  1.0;     xqs2[4] = -rthreei; 
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
int s8_ansq_funct(double frq[], double fsq[], double r, double s,
                  int iel, int nsansq)
{
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_ansq_funct */


/*----------------------------------------------------------------------*
 |                                                        m.gee 2/02    |
 |                       |
 *----------------------------------------------------------------------*/
int s8_ans_bbar_q(double **bop, double frq[], double fsq[],
                  double  *funct1q[],  double  *funct2q[],
                  double **deriv1q[],  double **deriv2q[],
                  double **akovr1q[],  double **akovr2q[],
                  double **a3kvpr1q[], double **a3kvpr2q[],
                  int iel, int numdf, int nsansq)
{
int inode, node_start;
int isamp;
double a1x1; 
double a1y1;
double a1z1;
double a3x1;
double a3y1;
double a3z1;
double a2x2;
double a2y2;
double a2z2;
double a3x2;
double a3y2;
double a3z2;
double a31x1;
double a31y1;
double a31z1;
double a32x2;
double a32y2;
double a32z2;
double p1k;
double p2k;
double pk1;
double pk2;
double fris;
double fsis;
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
      a1x1 = akovr1q[isamp][0][0];
      a1y1 = akovr1q[isamp][1][0];
      a1z1 = akovr1q[isamp][2][0];
      a3x1 = akovr1q[isamp][0][2];
      a3y1 = akovr1q[isamp][1][2];
      a3z1 = akovr1q[isamp][2][2];
      
      a2x2 = akovr2q[isamp][0][1];
      a2y2 = akovr2q[isamp][1][1];
      a2z2 = akovr2q[isamp][2][1];
      a3x2 = akovr2q[isamp][0][2];
      a3y2 = akovr2q[isamp][1][2];
      a3z2 = akovr2q[isamp][2][2];

      a31x1 = a3kvpr1q[isamp][0][0];
      a31y1 = a3kvpr1q[isamp][1][0];
      a31z1 = a3kvpr1q[isamp][2][0];

      a32x2 = a3kvpr2q[isamp][0][1];
      a32y2 = a3kvpr2q[isamp][1][1];
      a32z2 = a3kvpr2q[isamp][2][1];
      
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
 | modifications to metrics of shell ody due to ans for querschub       |
 *----------------------------------------------------------------------*/
int s8_ans_tvhe_q(double **gmkovr,double **gmkovc,double **gmkonr,double **gmkonc,
                  double *detr,   double *detc,
                  double **amkovr1q[], double **amkovc1q[], 
                  double **akovr1q[] , double **akovc1q[] ,
                  double **a3kvpr1q[], double **a3kvpc1q[],
                  double **amkovr2q[], double **amkovc2q[], 
                  double **akovr2q[] , double **akovc2q[] ,
                  double **a3kvpr2q[], double **a3kvpc2q[],
                  double frq[], double fsq[], double e3, int nansq, int iel)
{
int i;
#ifdef DEBUG 
dstrc_enter("s8_ans_tvhe_q");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<nansq; i++)
{
  /* gmkovc[0][2] +=    */
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_ans_tvhe_q */
