/*!----------------------------------------------------------------------
\file
\brief contains the routines ...
 - 's9_ans_colloqpoints': Gets the collocation points for Querschub-ANS,
             the shape functions and derivatives and calculates all
             necessary basis vectors an metrics at these points. For shell9
             the metrics have to be calculated for every kinematic layer
 - 's9_ans_colloqcoords': Gets the coordinates of the collocation points
 - 's9_ansq_funct': Shape functions for Querschub-ANS
 - 's9_ans_bbar_q': Calculates the values for alpha13 and alpha23 of the
             B-Operator Matrix -> Querschub-ANS
 - 's9_ans_tvhe_q': Modifies the metrics due to Querschub-ANS
 - 's9_ans_tvkg': Calculates the geometric stiffness Matrix due to Querschub-ANS


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief collocation points, shape functions and derivatives

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine calculates the coordinates of the collocation points, gets
the shape functions and derivatives at these points and calculates all
necessary basis vectors and metrics at these points
</pre>
\param  INT       nsansq       (i) number of collocation points
\param  INT       ans          (i) type of ANS (1: Querschub)
\param  DOUBLE    xr1[],...    (o) natural coordinates of collocation points
\param  DOUBLE   *funct1q[],...(o) shape functions at collocation points
\param  DOUBLE  **deriv1q[],...(o) shape functions derivatives at collocation points

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_ans_colloqpoints(INT       nsansq,     INT       iel,
                         INT       ans,        DIS_TYP   distyp,
                         DOUBLE    xr1[],      DOUBLE    xs1[],
                         DOUBLE    xr2[],      DOUBLE    xs2[],
                         DOUBLE   *funct1q[],  DOUBLE  **deriv1q[],
                         DOUBLE   *funct2q[],  DOUBLE  **deriv2q[],
                         DOUBLE  **xrefe,      DOUBLE ***a3r,
                         DOUBLE  **xcure,      DOUBLE ***a3c,
                         DOUBLE ***akovr1q[],  DOUBLE ***akonr1q[],
                         DOUBLE ***amkovr1q[], DOUBLE ***amkonr1q[],
                         DOUBLE ***a3kvpr1q[],
                         DOUBLE ***akovc1q[] , DOUBLE ***akonc1q[],
                         DOUBLE ***amkovc1q[], DOUBLE ***amkonc1q[],
                         DOUBLE ***a3kvpc1q[],
                         DOUBLE ***akovr2q[] , DOUBLE ***akonr2q[],
                         DOUBLE ***amkovr2q[], DOUBLE ***amkonr2q[],
                         DOUBLE ***a3kvpr2q[],
                         DOUBLE ***akovc2q[] , DOUBLE ***akonc2q[],
                         DOUBLE ***amkovc2q[], DOUBLE ***amkonc2q[],
                         DOUBLE ***a3kvpc2q[],
                         DOUBLE  **akovh,
                         DOUBLE  **akonh,
                         DOUBLE  **amkovh,
                         DOUBLE  **amkonh,
                         INT       num_klay)
{
INT    i;
DOUBLE det_dummy;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ans_colloqpoints");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------- get coordinates of collocation points */
s9_ans_colloqcoords(xr1,xs1,xr2,xs2,iel,ans);
/*--------------------------- loop over all collocation points (2 or 6) */
for (i=0; i<nsansq; i++)
{
  s9_funct_deriv(funct1q[i],deriv1q[i],xr1[i],xs1[i],distyp,1);

  s9_tvmr(xrefe,a3r,akovr1q[i],akonr1q[i],amkovr1q[i],amkonr1q[i],
          akovh,akonh,amkovh,amkonh,
          &det_dummy,funct1q[i],deriv1q[i],iel,a3kvpr1q[i],num_klay);

  s9_tvmr(xcure,a3c,akovc1q[i],akonc1q[i],amkovc1q[i],amkonc1q[i],
          akovh,akonh,amkovh,amkonh,
          &det_dummy,funct1q[i],deriv1q[i],iel,a3kvpc1q[i],num_klay);


  s9_funct_deriv(funct2q[i],deriv2q[i],xr2[i],xs2[i],distyp,1);

  s9_tvmr(xrefe,a3r,akovr2q[i],akonr2q[i],amkovr2q[i],amkonr2q[i],
          akovh,akonh,amkovh,amkonh,
          &det_dummy,funct2q[i],deriv2q[i],iel,a3kvpr2q[i],num_klay);

  s9_tvmr(xcure,a3c,akovc2q[i],akonc2q[i],amkovc2q[i],amkonc2q[i],
          akovh,akonh,amkovh,amkonh,
          &det_dummy,funct2q[i],deriv2q[i],iel,a3kvpc2q[i],num_klay);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_ans_colloqpoints */

/*!----------------------------------------------------------------------
\brief get natural coordinates of collocation points

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine calculates the natural coordinates of the collocation points
</pre>
\param  DOUBLE   xqr1[],...(o) natural coordinates of collocation points
\param  INT      iel       (i) number of nodes to this element
\param  INT      ans       (i) type of ANS (1: Querschub)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9_ans_colloqpoints()   [s9_ans.c]

*----------------------------------------------------------------------*/
void s9_ans_colloqcoords(DOUBLE xqr1[], DOUBLE xqs1[],
                         DOUBLE xqr2[], DOUBLE xqs2[],
                         INT    iel,    INT    ans)
{
DOUBLE rthreei;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ans_colloq");
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
   if (iel==8)/* 8-noded element */
   {
      dserror("ANS Q for QUAD8 not implemented 's9_ans_colloqcoords'  ");
   }
   if (iel==9)/* 9-noded element */
   {
      rthreei = 1.0 / (sqrt(3.0));
      xqr1[0] = -rthreei; xqs1[0] = -1.0;
      xqr1[1] = -rthreei; xqs1[1] =  0.0;
      xqr1[2] = -rthreei; xqs1[2] =  1.0;
      xqr1[3] =  rthreei; xqs1[3] = -1.0;
      xqr1[4] =  rthreei; xqs1[4] =  0.0;
      xqr1[5] =  rthreei; xqs1[5] =  1.0;

      xqr2[0] = -1.0;     xqs2[0] = -rthreei;
      xqr2[1] =  0.0;     xqs2[1] = -rthreei;
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
} /* end of s9_ans_colloq */


/*!----------------------------------------------------------------------
\brief get the shape functions for ans (Querschub)

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine calculates the shape functions for ANS (Querschub) at the
natural coordinates r,s
</pre>
\param  DOUBLE   frq[],... (o) shape functions for ANS
\param  INT      iel       (i) number of nodes to this element
\param  INT      nsansq    (i) number of collocation points

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_ansq_funct(DOUBLE frq[], DOUBLE fsq[], DOUBLE r, DOUBLE s,
                   INT    iel)
{
DOUBLE rthreei;
DOUBLE pr[3],ps[3];
DOUBLE qr[2],qs[2];
DOUBLE rr[3],rs[3];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ansq_funct");
#endif
/*----------------------------------------------------------------------*/
if (iel==4)
{
   frq[0] = 0.5 * (1.0 - s);
   frq[1] = 0.5 * (1.0 + s);

   fsq[0] = 0.5 * (1.0 - r);
   fsq[1] = 0.5 * (1.0 + r);
}
if (iel==8)/* 8-noded element */
{
   dserror("ANS Q for QUAD8 not implemented   's9_ansq_funct'   ");
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
} /* end of s9_ansq_funct */


/*!----------------------------------------------------------------------
\brief gets the values for alpha13/alpha23 (Querschub) of the B-Operator
 matrix

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine gets the values for alpha13/alpha23 (Querschub) of the
B-Operator matrix
</pre>
\param  DOUBLE  **bop         (i/o) B-Operator to be modified (alpha13,alpha23)
\param  DOUBLE    frq[],...    (i)  shape functions for ANS
\param  DOUBLE   *funct1q[],...(i)  shape functions at collocation points
\param  DOUBLE  **deriv1q[],...(i)  shape functions derivatives at collocation points
\param  INT       num_klay     (i)  number of kinematic layers to this element
\param  INT       klay         (i)  actual kinematic layer of which is looped
\param  INT       nsansq       (i)  number of collocation points

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_ans_bbar_q(DOUBLE  **bop,        DOUBLE    frq[],      DOUBLE fsq[],
                   DOUBLE   *funct1q[],  DOUBLE   *funct2q[],
                   DOUBLE  **deriv1q[],  DOUBLE  **deriv2q[],
                   DOUBLE ***akovc1q[],  DOUBLE ***akovc2q[],
                   DOUBLE ***a3kvpc1q[], DOUBLE ***a3kvpc2q[],
                   INT       iel,        INT       numdf,
                   INT       num_klay,   INT       klay,
                   DOUBLE    condfac,    INT       nsansq)
{
INT    jlay,inode, node_start;
INT    isamp;
DOUBLE a1x1, a1y1, a1z1;
DOUBLE a2x2, a2y2, a2z2;

DOUBLE a3xi1, a3yi1, a3zi1;
DOUBLE a3xj1, a3yj1, a3zj1;

DOUBLE a3xi2, a3yi2, a3zi2;
DOUBLE a3xj2, a3yj2, a3zj2;

DOUBLE a31xj1, a31yj1, a31zj1;
DOUBLE a32xj2, a32yj2, a32zj2;

DOUBLE p1k;
DOUBLE p2k;
DOUBLE pk1;
DOUBLE pk2;
DOUBLE fris;
DOUBLE fsis;

INT    idof, jdof;
DOUBLE fac,  fac1;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ans_bbar_q");
#endif
/*----------------------------------------------------------------------*/
  idof = 3 + (klay) * 3;
/*---------- loop over all kinematic layers ---------------------------*/
for (jlay=0; jlay<num_klay; jlay++)
{
  /* is jlay on trajectory of reference layer to klay (=ilay) ? */
  fac1 = s9notr(num_klay,klay,jlay);
  if (fac1 == 0.0) continue;

  /* independent part of xsi */
  fac = s9ksi(num_klay,klay,jlay,condfac);
  jdof = 3 + (jlay) * 3;

/*---------- loop over all nodes ---------------------------------------*/
  for (inode=0; inode<iel; inode++)
  {
     node_start = inode*numdf;

     for (isamp=0; isamp<nsansq; isamp++)
     {
        /*values on klay -> actual kinematic layer*/
        a3xi1 = akovc1q[isamp][0][2][klay];
        a3yi1 = akovc1q[isamp][1][2][klay];
        a3zi1 = akovc1q[isamp][2][2][klay];

        a3xi2 = akovc2q[isamp][0][2][klay];
        a3yi2 = akovc2q[isamp][1][2][klay];
        a3zi2 = akovc2q[isamp][2][2][klay];

        /*values on jlay -> layers to be looped (continuity matrix)*/
        a1x1 = akovc1q[isamp][0][0][jlay];
        a1y1 = akovc1q[isamp][1][0][jlay];
        a1z1 = akovc1q[isamp][2][0][jlay];

        a3xj1 = akovc1q[isamp][0][2][jlay];
        a3yj1 = akovc1q[isamp][1][2][jlay];
        a3zj1 = akovc1q[isamp][2][2][jlay];

        a2x2 = akovc2q[isamp][0][1][jlay];
        a2y2 = akovc2q[isamp][1][1][jlay];
        a2z2 = akovc2q[isamp][2][1][jlay];

        a3xj2 = akovc2q[isamp][0][2][jlay];
        a3yj2 = akovc2q[isamp][1][2][jlay];
        a3zj2 = akovc2q[isamp][2][2][jlay];

        a31xj1 = a3kvpc1q[isamp][0][0][jlay];
        a31yj1 = a3kvpc1q[isamp][1][0][jlay];
        a31zj1 = a3kvpc1q[isamp][2][0][jlay];

        a32xj2 = a3kvpc2q[isamp][0][1][jlay];
        a32yj2 = a3kvpc2q[isamp][1][1][jlay];
        a32zj2 = a3kvpc2q[isamp][2][1][jlay];

        p1k = funct1q[isamp][inode];
        p2k = funct2q[isamp][inode];

        pk1 = deriv1q[isamp][0][inode];
        pk2 = deriv2q[isamp][1][inode];

        fris = frq[isamp];
        fsis = fsq[isamp];
  /*-------------------------------------------------------------------- */
        if (klay == jlay)
        {
  /*--------------------------------------------------E13(CONST)-------- */
          bop[3][node_start+0]+= pk1*a3xi1*fris;
          bop[3][node_start+1]+= pk1*a3yi1*fris;
          bop[3][node_start+2]+= pk1*a3zi1*fris;

          bop[3][node_start+idof+0]+= p1k*a1x1*fris;
          bop[3][node_start+idof+1]+= p1k*a1y1*fris;
          bop[3][node_start+idof+2]+= p1k*a1z1*fris;
  /*--------------------------------------------------E23(CONST)-------- */
          bop[4][node_start+0]+= pk2*a3xi2*fsis;
          bop[4][node_start+1]+= pk2*a3yi2*fsis;
          bop[4][node_start+2]+= pk2*a3zi2*fsis;

          bop[4][node_start+idof+0]+= p2k*a2x2*fsis;
          bop[4][node_start+idof+1]+= p2k*a2y2*fsis;
          bop[4][node_start+idof+2]+= p2k*a2z2*fsis;
        }
        else if (klay != jlay)
        {
  /*--------------------------------------------------E13(CONST)-------- */
          bop[3][node_start+jdof+0]+= fac * pk1*a3xi1*fris;
          bop[3][node_start+jdof+1]+= fac * pk1*a3yi1*fris;
          bop[3][node_start+jdof+2]+= fac * pk1*a3zi1*fris;

          bop[3][node_start+idof+0]+= fac * p1k*a31xj1*fris;
          bop[3][node_start+idof+1]+= fac * p1k*a31yj1*fris;
          bop[3][node_start+idof+2]+= fac * p1k*a31zj1*fris;
  /*--------------------------------------------------E23(CONST)-------- */
          bop[4][node_start+jdof+0]+= fac * pk2*a3xi2*fsis;
          bop[4][node_start+jdof+1]+= fac * pk2*a3yi2*fsis;
          bop[4][node_start+jdof+2]+= fac * pk2*a3zi2*fsis;

          bop[4][node_start+idof+0]+= fac * p2k*a32xj2*fsis;
          bop[4][node_start+idof+1]+= fac * p2k*a32yj2*fsis;
          bop[4][node_start+idof+2]+= fac * p2k*a32zj2*fsis;
        }

     } /* end of (isamp=0; isamp<nsansq; isamp++) */
  } /* end of (inode=0; inode<iel; inode++) */
} /* end of loop over all kinematic layers*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_ans_bbar_q */



/*!----------------------------------------------------------------------
\brief modifies the metrics of shell body due to ans for querschub

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine modifies the metrics 'gmkovc' and 'gmkonc' of shell body
due to ans for querschub.
</pre>
\param  DOUBLE  **gmkovc       (i/o) kovariant metric in current configuration
\param  DOUBLE  **gmkonc       (i/o) kontravariant metric in current configuration
\param  DOUBLE    e3            (i)  local coordinate zeta (theta3) of act. mlay in act. klay
\param  DOUBLE   *klayhgt       (i)  hight of all kinematic layers in % of total thickness of shell
\param  DOUBLE   *mlayhgt       (i)  hight of all material layers in % of total thickness of actual kin. layer
\param  INT       num_klay      (i)  number of kinematic layers to this element
\param  INT       num_mlay      (i)  number of material layers to this kinematic layer
\param  INT       klay          (i)  actual kinematic layer of which is looped
\param  INT       mlay          (i)  actual material layer of within the actual kinematic layer

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_ans_tvhe_q(DOUBLE  **gmkovr,     DOUBLE  **gmkovc,
                   DOUBLE  **gmkonr,     DOUBLE  **gmkonc,
                   DOUBLE ***amkovc,     DOUBLE ***amkovr,
                   DOUBLE ***akovc,      DOUBLE ***akovr,
                   DOUBLE ***a3kvpc,     DOUBLE ***a3kvpr,
                   DOUBLE   *detr,       DOUBLE   *detc,
                   DOUBLE ***amkovr1q[], DOUBLE ***amkovc1q[],
                   DOUBLE ***amkovr2q[], DOUBLE ***amkovc2q[],
                   DOUBLE    frq[],      DOUBLE    fsq[],
                   DOUBLE    e3,         INT       nansq,
                   DOUBLE    h,           /* total thickness of this element */
                   DOUBLE   *klayhgt,     /* hight of kin layer in % of total thickness of shell */
                   DOUBLE   *mlayhgt,     /* hight of mat layer in % of adjacent kin layer */
                   INT       num_klay,    /* number of kin layers to this element */
                   INT       klay,        /* actual kin layer */
                   INT       mlay,        /* actual mat layer of this kin layer */
                   DOUBLE    condfac)
{
INT i,jlay;
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

DOUBLE deltah, h_mlay, h_kl;
DOUBLE zeta_kl,zeta,det_dummy;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ans_tvhe_q");
#endif
/*----------------------------------------------------------------------*/
gmkovc[0][0] = gmkovr[0][0] + (amkovc[0][0][klay]-amkovr[0][0][klay]);
gmkovc[1][1] = gmkovr[1][1] + (amkovc[1][1][klay]-amkovr[1][1][klay]);
gmkovc[2][2] = gmkovr[2][2] + (amkovc[2][2][klay]-amkovr[2][2][klay]);
gmkovc[0][1] = gmkovr[0][1] + (amkovc[0][1][klay]-amkovr[0][1][klay]);
gmkovc[0][2] = gmkovr[0][2]; /* + (amkovc[0][2][klay]-amkovr[0][2][klay])*/
gmkovc[1][2] = gmkovr[1][2]; /* + (amkovc[1][2][klay]-amkovr[1][2][klay])*/
/*----------------------------------------------------------------------*/
/*- calculate zeta_kl of kinematic layer due to local zeta_ml of material layer -> old s9tmtr --*/
h_kl   = (klayhgt[klay]/100.)*h;              /* absolute hight of the actual kinematic layer */
deltah = (mlayhgt[mlay]/100.)*h_kl;           /* absolute hight of the actual material layer */
h_mlay   = 0.0;
for (i=0; i<=mlay; i++)                        /* sum of the absolute hights of the material layers  */
{                                             /* within the actual kinematic layer up to the actual */
   h_mlay += (mlayhgt[i]/100.)*h_kl;          /* material layer */
}
zeta_kl = -1. + (-deltah*(1.-e3)+2.*h_mlay)/h_kl;  /* equals the theta3 value of the act. kinematic layer, see equation (5.45) in Dis. Braun */

/*--- loop over all kinematic layers to get the continuity coeficients (Pic. 7.4 on p. 116 in Dis. Braun)-> old s9mmtr ------*/
for (jlay=0; jlay<num_klay; jlay++)
{
   zeta = zeta_kl;
   zeta = s9con(zeta,num_klay,klay,jlay,condfac);
   if (zeta != 0.0)  /*---- interpolation -------- g1,g2 (kov.) -------*/
   {
      b11c=0.0;
      b12c=0.0;
      b21c=0.0;
      b22c=0.0;
      b31c=0.0;
      b32c=0.0;

      b11r=0.0;
      b12r=0.0;
      b21r=0.0;
      b22r=0.0;
      b31r=0.0;
      b32r=0.0;

      for (i=0; i<3; i++) b11c += akovc[i][0][0]   *a3kvpc[i][0][jlay];
      for (i=0; i<3; i++) b12c += akovc[i][0][0]   *a3kvpc[i][1][jlay];
      for (i=0; i<3; i++) b21c += akovc[i][1][0]   *a3kvpc[i][0][jlay];
      for (i=0; i<3; i++) b22c += akovc[i][1][0]   *a3kvpc[i][1][jlay];
      for (i=0; i<3; i++) b31c += akovc[i][2][klay]*a3kvpc[i][0][jlay];
      for (i=0; i<3; i++) b32c += akovc[i][2][klay]*a3kvpc[i][1][jlay];

      for (i=0; i<3; i++) b11r += akovr[i][0][0]   *a3kvpr[i][0][jlay];
      for (i=0; i<3; i++) b12r += akovr[i][0][0]   *a3kvpr[i][1][jlay];
      for (i=0; i<3; i++) b21r += akovr[i][1][0]   *a3kvpr[i][0][jlay];
      for (i=0; i<3; i++) b22r += akovr[i][1][0]   *a3kvpr[i][1][jlay];
      for (i=0; i<3; i++) b31r += akovr[i][2][klay]*a3kvpr[i][0][jlay];
      for (i=0; i<3; i++) b32r += akovr[i][2][klay]*a3kvpr[i][1][jlay];

      gmkovc[0][0] = gmkovc[0][0] + zeta*2.0*(b11c-b11r);
      gmkovc[1][1] = gmkovc[1][1] + zeta*2.0*(b22c-b22r);
      gmkovc[2][2] = gmkovc[2][2] ;
      gmkovc[0][1] = gmkovc[0][1] + zeta*(b21c+b12c-b21r-b12r);
      gmkovc[0][2] = gmkovc[0][2] + zeta*(b31c-b31r);
      gmkovc[1][2] = gmkovc[1][2] + zeta*(b32c-b32r);
      /*----------------------------------------------------------------------*/
      for (i=0; i<nansq; i++)
      {
         gmkovc[0][2] += (amkovc1q[i][0][2][jlay]-amkovr1q[i][0][2][jlay]) * frq[i];
         gmkovc[1][2] += (amkovc2q[i][1][2][jlay]-amkovr2q[i][1][2][jlay]) * fsq[i];
      }
   }
} /*================================ end loop over all layers ========*/

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
} /* end of s9_ans_tvhe_q */



/*!----------------------------------------------------------------------
\brief geometric stiffness matrix kg with ans for querschub

<pre>                     m.gee 2/02              modified by    sh 02/03
This routine calculates the geometric stiffness matrix kg if the ans
method for querschub is used
</pre>
\param  DOUBLE  **estif        (i/o) element stiffness metrix to be modified
\param  DOUBLE  *stress_r       (i)  stress resultants
\param  DOUBLE  *funct          (i)  shape functions at integration point
\param  DOUBLE **deriv          (i)  shape function derivatives at integration point
\param  DOUBLE   weight         (i)  weight at integration point
\param  DOUBLE   frq[],...      (i)  shape functions for ANS
\param  DOUBLE  *funct1q[],...  (i)  shape functions at collocation points
\param  DOUBLE **deriv1q[],...  (i)  shape function derivatives at collocation points
\param  INT      ansq           (i)  flag if querschub ans
\param  INT      nsansq         (i)  number of collocation points
\param  INT      num_klay       (i)  number of kinematic layers to this element
\param  INT      klay           (i)  actual kinematic layer of which is looped

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_ans_tvkg(DOUBLE **estif,     DOUBLE  *stress_r,
                 DOUBLE  *funct,     DOUBLE **deriv,
                 INT      numdf,     INT      iel,
                 DOUBLE   weight,
                 DOUBLE   frq[],     DOUBLE   fsq[],
                 DOUBLE  *funct1q[], DOUBLE  *funct2q[],
                 DOUBLE **deriv1q[], DOUBLE **deriv2q[],
                 INT      ansq,      INT       nsansq,
                 INT      klay,                       /* actual kin layer */
                 INT      num_klay)                   /* number of kin layers to this element */
{
INT     i,inode,jnode;
INT     idof,jdof;
INT     jlay;
DOUBLE  fac1;

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
DOUBLE  y1,y2,y;
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

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_ans_tvkg");
#endif
/*----------------------------------------------------------------------*/
sn11 = stress_r[0];
sn21 = stress_r[1];
sn22 = stress_r[2];
sn31 = stress_r[3];
sn32 = stress_r[4];
sn33 = stress_r[5];
sm11 = stress_r[6];
sm21 = stress_r[7];
sm22 = stress_r[8];
sm31 = stress_r[9];
sm32 = stress_r[10];
/*----------------------------------------------------------------------*/
idof = 3 + (klay) * 3;
/*----------------------------------------------------------------------*/
/*---------- loop over all kinematic layers ---------------------------*/
for (jlay=0; jlay<num_klay; jlay++)
{
  /* is jlay on trajectory of reference layer to klay (=ilay) ? */
  fac1 = s9notr(num_klay,klay,jlay);
  if (fac1 == 0.0) continue;
  jdof = 3 + (jlay) * 3;

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

        /*-------------------------------------E11,E12,E22(CONST)*/
        xn = (sn11*d11 + sn21*(d12+d21) + sn22*d22) * weight;
        /*---------------------------------------E11,E12,E22(LIN)*/
        xm = (sm11*d11 + sm21*(d12+d21) + sm22*d22) * weight;


        /*----------------------------------------- no ans for querschub */
        if (ansq==0)
        {
           pd1ij = deriv[0][inode] * pj;
           pd1ji = deriv[0][jnode] * pi;
           pd2ij = deriv[1][inode] * pj;
           pd2ji = deriv[1][jnode] * pi;
           /*-----------------------------------------E13,E23(CONST)*/
           yu = (sn31*pd1ji + sn32*pd2ji) * weight;
           yo = (sn31*pd1ij + sn32*pd2ij) * weight;
        }
        else/*----------------------------------------- ans for querschub */
        {
           yu = 0.0;
           yo = 0.0;
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

        /*-------------------------------------------E13,E23(LIN)*/
        y1 = (sm31*pd1ij + sm32*pd2ij) * weight;
        y2 = (sm31*pd1ji + sm32*pd2ji) * weight;
        y  = y1 + y2;
        /*---------------------------------------------E33(CONST)*/
        z  = pi*pj*sn33*weight;

        if (klay == jlay)
        {
           estif[inode*numdf+0][jnode*numdf+0] += xn;
           estif[inode*numdf+1][jnode*numdf+1] += xn;
           estif[inode*numdf+2][jnode*numdf+2] += xn;

           estif[inode*numdf+idof+0][jnode*numdf+0] += (xm+yu);
           estif[inode*numdf+idof+1][jnode*numdf+1] += (xm+yu);
           estif[inode*numdf+idof+2][jnode*numdf+2] += (xm+yu);

           estif[inode*numdf+0][jnode*numdf+idof+0] += (xm+yo);
           estif[inode*numdf+1][jnode*numdf+idof+1] += (xm+yo);
           estif[inode*numdf+2][jnode*numdf+idof+2] += (xm+yo);

           estif[inode*numdf+idof+0][jnode*numdf+idof+0] += (y+z);
           estif[inode*numdf+idof+1][jnode*numdf+idof+1] += (y+z);
           estif[inode*numdf+idof+2][jnode*numdf+idof+2] += (y+z);

           /*make symmetric*/
           if (inode!=jnode)
           {
             estif[jnode*numdf+0][inode*numdf+0] += xn;
             estif[jnode*numdf+1][inode*numdf+1] += xn;
             estif[jnode*numdf+2][inode*numdf+2] += xn;

             estif[jnode*numdf+0][inode*numdf+idof+0] += (xm+yu);
             estif[jnode*numdf+1][inode*numdf+idof+1] += (xm+yu);
             estif[jnode*numdf+2][inode*numdf+idof+2] += (xm+yu);

             estif[jnode*numdf+idof+0][inode*numdf+0] += (xm+yo);
             estif[jnode*numdf+idof+1][inode*numdf+1] += (xm+yo);
             estif[jnode*numdf+idof+2][inode*numdf+2] += (xm+yo);

             estif[jnode*numdf+idof+0][inode*numdf+idof+0] += (y+z);
             estif[jnode*numdf+idof+1][inode*numdf+idof+1] += (y+z);
             estif[jnode*numdf+idof+2][inode*numdf+idof+2] += (y+z);
           }
        }
        else if (klay != jlay)
        {
           estif[jnode*numdf+jdof+0][inode*numdf+0] += xn;
           estif[jnode*numdf+jdof+1][inode*numdf+1] += xn;
           estif[jnode*numdf+jdof+2][inode*numdf+2] += xn;

           estif[jnode*numdf+0][inode*numdf+jdof+0] += xn;
           estif[jnode*numdf+1][inode*numdf+jdof+1] += xn;
           estif[jnode*numdf+2][inode*numdf+jdof+2] += xn;

           estif[jnode*numdf+jdof+0][inode*numdf+idof+0] += yo;
           estif[jnode*numdf+jdof+1][inode*numdf+idof+1] += yo;
           estif[jnode*numdf+jdof+2][inode*numdf+idof+2] += yo;

           estif[jnode*numdf+idof+0][inode*numdf+jdof+0] += yu;
           estif[jnode*numdf+idof+1][inode*numdf+jdof+1] += yu;
           estif[jnode*numdf+idof+2][inode*numdf+jdof+2] += yu;

           /*make symmetric*/
           if (inode!=jnode)
           {
             estif[inode*numdf+0][jnode*numdf+jdof+0] += xn;
             estif[inode*numdf+1][jnode*numdf+jdof+1] += xn;
             estif[inode*numdf+2][jnode*numdf+jdof+2] += xn;

             estif[inode*numdf+jdof+0][jnode*numdf+0] += xn;
             estif[inode*numdf+jdof+1][jnode*numdf+1] += xn;
             estif[inode*numdf+jdof+2][jnode*numdf+2] += xn;

             estif[inode*numdf+idof+0][jnode*numdf+jdof+0] += yo;
             estif[inode*numdf+idof+1][jnode*numdf+jdof+1] += yo;
             estif[inode*numdf+idof+2][jnode*numdf+jdof+2] += yo;

             estif[inode*numdf+jdof+0][jnode*numdf+idof+0] += yu;
             estif[inode*numdf+jdof+1][jnode*numdf+idof+1] += yu;
             estif[inode*numdf+jdof+2][jnode*numdf+idof+2] += yu;
           }
        }

     } /* end loop over jnode */
  } /* end loop over inode */

} /* end of loop over all kinematic layers*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_ans_tvkg */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
