/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_cal_stress' which evaluates the element
       stresses for 2D isoparametric degenerated element
       contains the routine 'w1_mami' which evaluates the principal
       stresses and directions for 2D isoparametric degenerated element
       contains the routine 'w1rsn' which returns R/S coordinates of
       gauss integration points 2D isoparametric degenerated element
       contains the routine 'w1recs' which extrapolates from gauss
       points for rectangles

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | evaluate element stresses                                 al 9/01    |
 | 2-D isoparametric degenerated element                                |
 *----------------------------------------------------------------------*/
void w1_cal_stress(ELEMENT   *ele,
                   W1_DATA   *data,
                   MATERIAL  *mat,
                   INT	      kstep,  /* number of current load step */
		       INT        init)
{
INT                 i,j,k,l;            /* some loopers */
INT                 nir=0;          /* num GP in r/s/t direction */
INT                 nis=0;          /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
INT                 ip;
INT                 intc=0;      /* "integration case" for tri-element */
INT	              it=0;      /* flag for transformation global/local   */
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 1;/* controls evaluation of new stresses    */
const INT           numdf  = 2;
const INT           numeps = 3;
const INT           numstr = 4;

DOUBLE              fac;
DOUBLE              e1=0.0;            /*GP-coords*/
DOUBLE              e2=0.0;            /*GP-coords*/
DOUBLE              facr=1.0;        /* weights at GP */
DOUBLE              facs=1.0;        /* weights at GP */
DOUBLE              dum;
DOUBLE deltad[8], help[4];
DOUBLE knninv[4][4], knc[8][4];

static ARRAY    D_a;      /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;  /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;  /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;    /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    xjm0_a;    /* jacobian matrix at r,s=0*/
static DOUBLE **xjm0;
static ARRAY    xji_a;    /* inverse of jacobian matrix */
static DOUBLE **xji;
static ARRAY    F_a;      /* dummy matrix for saving stresses at gauss point*/
static DOUBLE  *F;
static ARRAY    bop_a;    /* B-operator */
static DOUBLE **bop;
static ARRAY    gop_a;    /* incomp_modes: G-operator */
static DOUBLE  *gop;
static ARRAY    alpha_a;  /* incomp-modes: internal dof */
static DOUBLE  *alpha;
static ARRAY    spar_a;   /* function parameters for extrapolation */
static DOUBLE **spar;
static ARRAY    transm_a; /* transformation matrix sig(loc)-->sig(glob) */
static DOUBLE **transm;
static ARRAY    transmi_a;/* inverse transformation matrix sig(glob)-->sig(loc) */
static DOUBLE **transmi;
static ARRAY    work_a;
static DOUBLE **work;     /* working array */

DOUBLE det,det0;
INT node, npoint;
DOUBLE fv, r, s;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_cal_stress");
#endif
/*--------------------------------------------- allocate working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_WALL1 ,"DA");
D         = amdef("D"      ,&D_a   ,6,6             ,"DA");
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf     ,"DA");
xjm0      = amdef("xjm0"   ,&xjm0_a,numdf,numdf     ,"DA");
xji       = amdef("xji"    ,&xji_a ,numdf,numdf     ,"DA");
F         = amdef("F"      ,&F_a,4,1     ,"DV");
bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_WALL1),"DA");
gop       = amdef("gop"  ,&gop_a ,4,1,"DV");
alpha     = amdef("alpha"  ,&alpha_a ,4,1,"DV");
transm    = amdef("transm"  ,&transm_a ,numstr,numstr,"DA");
transmi   = amdef("transmi"  ,&transmi_a ,numstr,numstr,"DA");
spar      = amdef("spar"  ,&spar_a ,numstr,16,"DA");
work      = amdef("work"  ,&work_a ,4,4,"DA");
goto end;
}
/*--------------------------------------------- deallocate working arrays */
if (init==-1)
{
amdel(&funct_a);
amdel(&deriv_a);
amdel(&D_a);
amdel(&xjm_a);
amdel(&xjm0_a);
amdel(&xji_a);
amdel(&F_a);
amdel(&bop_a);
amdel(&gop_a);
amdel(&alpha_a);
amdel(&transm_a);
amdel(&transmi_a);
amdel(&spar_a);
amdel(&work_a);
goto end;
}
/*------------------------------------------- integration parameters ---*/
w1intg(ele,data,1);
/*------------------------------------------- integration parameters ---*/
iel     = ele->numnp;
nd      = numdf * iel;
/*------- get integraton data ---------------------------------------- */
switch (ele->distyp)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   nir = ele->e.w1->nGP[0];
   nis = ele->e.w1->nGP[1];
break;
case tri3: case tri6: /* --> tri - element */
   nir  = ele->e.w1->nGP[0];
   nis  = 1;
   intc = ele->e.w1->nGP[1]-1;
break;
default:
   dserror("ele->distyp unknown!");
} /* end switch(ele->distyp) */
npoint  = nir*nis;

if(ele->e.w1->modeltype == incomp_mode)
{
 /*------------------ shape functions and their derivatives at r,s=0 ---*/
 w1_funct_deriv(funct,deriv,0,0,ele->distyp,1);
 /*-------------------------------- compute jacobian matrix at r,s=0 ---*/
 w1_jaco (deriv,xjm0,&det0,ele,iel);
 amzero(&alpha_a);
 if(mat->mattyp != m_stvenant)
 {
  for (i=0;i<4;i++)
  alpha[i] = ele->e.w1->elewa[0].imodewa[0].alpha[i];
 }
 else
 {
  for (i=0; i<4; i++)
  {
    deltad[2*i]   = ele->node[i]->sol.a.da[0][0];
    deltad[2*i+1] = ele->node[i]->sol.a.da[0][1];
    for(j=0;j<4;j++)
      knninv[i][j]=ele->e.w1->elewa[0].imodewa[0].knninv[i][j];
    for(j=0;j<8;j++)
      knc[j][i]=ele->e.w1->elewa[0].imodewa[0].knc[j][i];
  }
/*-------------------------------------------- evaluate alpha ---*/
  for (i=0; i<4; i++)
  {
    dum=0.0;
    for (k=0; k<4; k++)
    {
     help[k] = 0.0;
     for (l=0; l<8; l++)
     {
       help[k] += knc[l][k] * deltad[l];
     }
     /*-   help[k] += fintn[k];--*/
     dum += knninv[i][k] * help[k];
    }
    alpha[i] -= dum ;
  }
 }
}
/*================================================ integration loops ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
   for (ls=0; ls<nis; ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(ele->distyp)
      {
      case quad4: case quad8: case quad9:  /* --> quad - element */
       e1   = data->xgrr[lr];
       e2   = data->xgss[ls];
      break;
      case tri3: case tri6:  /* --> tri - element */
	 e1   = data->txgr[lr][intc];
	 e2   = data->txgs[lr][intc];
      break;
      default:
         dserror("ele->distyp unknown!");
      } /* end switch(ele->distyp) */
      ip++;
       /*------------------------- shape functions and their derivatives */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------------------------------------ compute jacobian matrix ---*/
      w1_jaco (deriv,xjm,&det,ele,iel);
      fac = facr * facs * det;
      /*--------------------------------------- calculate operator B ---*/
      amzero(&bop_a);
      w1_bop(bop,deriv,xjm,det,iel);
      /*--------------------------------------- calculate operator G ---*/
      if(ele->e.w1->modeltype == incomp_mode)
      {
        amzero(&gop_a);
        w1_gop(gop,xjm0,det0,det,e1,e2);
      }
     /*------------------------------------------ call material law ---*/
/*-----------------------------------------------------fh 06/02---*/
      newval=1; /* Flag to calculate stresses */
      w1_call_mat(ele, mat,ele->e.w1->wtype, bop,gop,alpha, ip, F,D, istore,newval);

      /* transformation of global stresses into local stresses and vice versa */

      switch(ele->e.w1->stresstyp){
      case w1_rs:
      w1_tram(xjm,transm,transmi,work);
      w1_lss(F,transmi,transm,it);
      break;
      default:
      break;
      }
      for (i=0; i<4; i++)
      {
	ele->e.w1->stress_GP.a.d3[kstep][i][ip]= F[i];
      }


      w1_mami(F, &ele->e.w1->stress_GP.a.d3[kstep][4][ip],
                 &ele->e.w1->stress_GP.a.d3[kstep][5][ip],
                 &ele->e.w1->stress_GP.a.d3[kstep][6][ip]);

   }/*============================================= end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
for (node=1; node<=iel; node++)
{
/*----------------------------------- get local coordinates of nodes ---*/
  r = w1rsn (node,1,iel) ;
  s = w1rsn (node,2,iel) ;
/*------------------------------------------- extrapolate values now ---*/
  /* stress in xx-direction (sig_11) */
  w1recs (&fv,r,s,&ele->e.w1->stress_GP.a.d3[kstep][0][0],spar[0],npoint,node);
  ele->e.w1->stress_ND.a.d3[kstep][0][node-1] = fv;

  /* stress in yy-direction (sig_22) */
  w1recs (&fv,r,s,&ele->e.w1->stress_GP.a.d3[kstep][1][0],spar[1],npoint,node);
  ele->e.w1->stress_ND.a.d3[kstep][1][node-1] = fv;

  /* stress in xy-direction (sig_12) */
  w1recs (&fv,r,s,&ele->e.w1->stress_GP.a.d3[kstep][2][0],spar[2],npoint,node);
  ele->e.w1->stress_ND.a.d3[kstep][2][node-1] = fv;

  /* stress in zz-direction (sig_33), sometimes set to zero */
  w1recs (&fv,r,s,&ele->e.w1->stress_GP.a.d3[kstep][3][0],spar[3],npoint,node);
  ele->e.w1->stress_ND.a.d3[kstep][3][node-1] = fv;

  F[0] = ele->e.w1->stress_ND.a.d3[kstep][0][node-1];
  F[1] = ele->e.w1->stress_ND.a.d3[kstep][1][node-1];
  F[2] = ele->e.w1->stress_ND.a.d3[kstep][2][node-1];
  F[3] = ele->e.w1->stress_ND.a.d3[kstep][3][node-1];

  w1_mami(F, &ele->e.w1->stress_ND.a.d3[kstep][4][node-1],
             &ele->e.w1->stress_ND.a.d3[kstep][5][node-1],
             &ele->e.w1->stress_ND.a.d3[kstep][6][node-1]);
}
/*----------------------------------------------------------------*/

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_cal_stress */
/*----------------------------------------------------------------------*
 |                                                           al 9/01    |
 |      principal stresses and directions                               |
 |      2D isoparametric degenerated element                            |
 |----------------------------------------------------------------------|
 |      angle of principal stresses:                                    |
 |                                               2*TAU XY               |
 |                           TAN (2*ALFA) = -------------------         |
 |                                          (SIGMA X - SIGMA Y)         |
 *----------------------------------------------------------------------*/
void w1_mami(DOUBLE *stress,
             DOUBLE *fps, /* FIRST  PRINCIPAL STRESS       */
             DOUBLE *sps, /* SECOND PRINCIPAL STRESS       */
             DOUBLE *aps) /* ANGLE OF PRINCIPAL DIRECTION  */
{
/*----------------------------------------------------------------------*/
DOUBLE 		ag;
static DOUBLE 	rad=0.;
char		jobz[1];
char		uplo[1];
INT		n=2;
INT		lda=2;
DOUBLE		A[4];
DOUBLE		W[2];
INT		lwork=5;
DOUBLE		work[5];
INT		info=0;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_mami");
#endif
/*----------------------------------------------------------------------*/
jobz[0]='V';
uplo[0]='L';

A[0]=stress[0];
A[1]=stress[2];
A[2]=stress[2];
A[3]=stress[1];

/*------calculate eigenvalues in ascending order and corresponding------*/
/*------orthonormal eigenvectors with LAPACK FORTRAN ROUTINE------------*/

/* We've got a symmetric 2x2 matrix. We can do this by hand. */
/*
dsyev(jobz,
      uplo,
      &n,
      &(A[0]),
      &lda,
      &(W[0]),
      &(work[0]),
      &lwork,
      &info);
*/

/* eigenvalues */
work[0] = sqrt(A[0]*A[0] - 2*A[0]*A[3] + 4*A[1]*A[1] + A[3]*A[3]);
W[0] = .5*(A[0] + A[3] - work[0]);
W[1] = .5*(A[0] + A[3] + work[0]);

/* eigenvectors */
/* first component of both eigenvectors */
work[1] = -.5/A[1]*(-A[0] + A[3] + work[0]);
work[2] = -.5/A[1]*(-A[0] + A[3] - work[0]);

/* lengths for normalization */
work[3] = sqrt(work[1]*work[1] + 1);
work[4] = sqrt(work[2]*work[2] + 1);

/* fill result in A. fortran style: colunm wise */
A[0] = work[1]/work[3];
A[1] = 1./work[3];
A[2] = work[2]/work[4];
A[3] = 1./work[4];

  if (rad == 0.)
  {
    rad = atan(1.) / 45.;
  }
ag=acos(A[2]) /rad;

if (A[3]<0.)
  {
    ag=180-ag;
  }

/*----------------------------------------- store principal stresses ---*/
*fps=W[1];
*sps=W[0];
*aps=ag;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_mami */

/*----------------------------------------------------------------------*
 |                                                           al 9/01    |
 |      returns R/S coordinates of gauss integration points             |
 |      2D isoparametric degenerated element                            |
 |----------------------------------------------------------------------|
 |      xr489 ---> R,S values for rectangles with 4,8,9 nodes           |
 |      xr16  ---> R,S values for rectangles with 12,16 nodes           |
 |      xt36  ---> R,S values for triangle   with 3,6   nodes           |
 |      xt10  ---> R,S values for triangle   with 10    nodes           |
 *----------------------------------------------------------------------*/
DOUBLE w1rsn (INT node, /* number of actual integration point */
              INT  irs, /* r/s identifier */
              INT  iel) /* number of nodes at actual element */
{
  static DOUBLE xr489[18] = { 1.,1.,-1.,1.,-1.,-1.,1.,
          -1.,0.,1.,-1.,0.,0.,-1.,1.,0.,0.,0. };
  static DOUBLE xr16[32] = { 1.,1.,-1.,1.,-1.,-1.,1.,
          -1.,.3333333333333333,1.,-.3333333333333333,1.,-1.,
          .3333333333333333,-1.,-.3333333333333333,-.3333333333333333,-1.,
          .3333333333333333,-1.,1.,-.3333333333333333,1.,.3333333333333333,
          .3333333333333333,.3333333333333333,-.3333333333333333,
          .3333333333333333,-.3333333333333333,-.3333333333333333,
          .3333333333333333,-.3333333333333333 };
  static DOUBLE xt36[12] = { 0.,0.,1.,0.,0.,1.,.5,0.,.5,.5,0.,.5 };
  static DOUBLE xt10[20] = { 0.,0.,1.,0.,0.,1.,
          .3333333333333333,0.,.6666666666666667,0.,.6666666666666667,
          .3333333333333333,.3333333333333333,.6666666666666667,0.,
          .6666666666666667,0.,.3333333333333333,.3333333333333333,
          .3333333333333333 };
  static INT label[16] = { 5,5,3,1,5,3,5,1,1,4,5,2,5,5,5,2 };

  /* System generated locals */
  DOUBLE ret_val;

  /* Local variables */
  static INT jump;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1rsn");
#endif
/*---------------------------------- labels for element nodal points ---*/
    jump = label[iel - 1];
    switch (jump) {
	case 1:  goto L1;
	case 2:  goto L2;
	case 3:  goto L3;
	case 4:  goto L4;
	case 5:  goto L5;
    }
/*----------------------------- set r/s coordinates for nodal points ---*/

/* although bitwise left-shift <<, ie node << 1, appears pretty cool, it
   translates simply to 2*node, yours bb 8/05 */
L1:
    ret_val = xr489[irs + (node << 1) - 3];
    goto L100;
L2:
    ret_val = xr16[irs + (node << 1) - 3];
    goto L100;
L3:
    ret_val = xt36[irs + (node << 1) - 3];
    goto L100;
L4:
    ret_val = xt10[irs + (node << 1) - 3];
    goto L100;
L5:
    dserror("illegal number of nodes given");

L100:

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ret_val;
} /* w1rsn */
/*----------------------------------------------------------------------*
 |                                                           al 9/01    |
 |      extrapolation form gauss points for rectangles                  |
 *----------------------------------------------------------------------*/
void w1recs(DOUBLE *funval,/* function value  */
            DOUBLE  r,     /* r/s coordinates */
            DOUBLE  s,
            DOUBLE *fval,  /* function values at gauss points */
            DOUBLE *fpar,  /* function parameters */
            INT     igauss,/* number of gauss points */
            INT     icode) /* ==1 initialize function parameters */
                           /* > 1            function evaluation */
{
  static INT i;
  static DOUBLE f1, f2, f3, f4, f5,
                p1, p2, p3, p4, p5, p6, p7, p8, p9,
                r3, s3,
                p10, p11, p12, p13, p14, p15, p16,
                ra, rb, rq, rr, rs, ss,
                qa2, qb2, ra3, rb3, raa, rab, rbb, qab, qba;
/*----------------------------------------------------------------------*/
  /* Parameter adjustments */
  --fpar;
  --fval;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1recs");
#endif
/*------------------------------------- evaluate function parameters ---*/
    if (icode == 1) {
	if (igauss >= 1) {
	    p1 = fval[1];
	    p2 = fval[2];
	}
	if (igauss >= 3) {
	    p3 = fval[3];
        }
	if (igauss >= 4) {
	    p2 = fval[2];
	    p3 = fval[3];
	    p4 = fval[4];
	}
	if (igauss >= 9) {
	    p5 = fval[5];
	    p6 = fval[6];
	    p7 = fval[7];
	    p8 = fval[8];
	    p9 = fval[9];
	}
	if (igauss >= 16) {
	    p10 = fval[10];
	    p11 = fval[11];
	    p12 = fval[12];
	    p13 = fval[13];
	    p14 = fval[14];
	    p15 = fval[15];
	    p16 = fval[16];
	}
	if (igauss == 1) {
/*----------------------------------------- constant stress function ---*/
	    fpar[1] = p1;
	} else if (igauss == 2) {
/*------------------------ constant stress function for TRIANGLES ------*/
	    fpar[1] = p1;
	    fpar[2] = p2;
        } else if (igauss == 3) {
/*-------------------------- linear stress function for TRIANGLES ------*/
	    fpar[1] = p1;
	    fpar[2] = p2;
	    fpar[3] = p3;
	} else if (igauss == 4) {
/*------------------------------------------- linear stress function ---*/
	    f1 = .25;
	    f2 = f1 * sqrt(3.);
	    f3 = f1 * 3.;
	    fpar[1] = f2 * (-p1 + p3 - p2 + p4);
	    fpar[2] = f3 * (p1 - p3 - p2 + p4);
	    fpar[3] = f2 * (-p1 - p3 + p2 + p4);
	    fpar[4] = f1 * (p1 + p3 + p2 + p4);
	} else if (igauss == 9) {
/*---------------------------------------- quadratic stress function ---*/
	    f3 = .83333333333333337;
	    f4 = f3 / 2.;
	    f1 = f3 * f3;
	    f2 = f1 * sqrt(.59999999999999998);
	    f5 = f3 * sqrt(.59999999999999998);
	    fpar[1] = f1 * (p1 + p7 + p3 + p9 - (p4 + p2 + p8 + p6) * 2. + p5
		    * 4.);
	    fpar[2] = f2 * (-p1 - p7 + p3 + p9 + (p4 - p6) * 2.);
	    fpar[3] = f2 * (-p1 + p7 - p3 + p9 + (p2 - p8) * 2.);
	    fpar[4] = f3 * (p2 + p8 - p5 * 2.);
	    fpar[5] = f4 * (p1 - p7 - p3 + p9);
	    fpar[6] = f3 * (p4 + p6 - p5 * 2.);
	    fpar[7] = f5 * (-p2 + p8);
	    fpar[8] = f5 * (-p4 + p6);
	    fpar[9] = p5;
	} else if (igauss == 16) {
/*-------------------------------------------- cubic stress function ---*/
	    ra = sqrt(1.2) * 2.;
	    rb = sqrt((3. - ra) / 7.);
	    ra = sqrt((ra + 3.) / 7.);
	    rab = ra * rb;
	    raa = ra * ra;
	    rbb = rb * rb;
	    ra3 = ra * raa;
	    rb3 = rb * rbb;
	    qab = ra / rb;
	    qba = rb / ra;
	    qa2 = qab * qab;
	    qb2 = qba * qba;
	    rq = .63802083333333337;
	    fpar[1] = (p1 - p13 - p4 + p16) / raa + (p6 - p10 - p7 + p11) /
		    rbb + (-p5 + p9 - p2 + p14 + p3 - p15 + p8 - p12) / rab;
	    fpar[2] = (-p1 + p13 + p2 - p14 + p3 - p15 - p4 + p16) / ra + (p5
		    - p9 - p6 + p10 - p7 + p11 + p8 - p12) / rb;
	    fpar[3] = (-p1 + p5 + p9 - p13 + p4 - p8 - p12 + p16) / ra + (p2
		    - p6 - p10 + p14 - p3 + p7 + p11 - p15) / rb;
	    fpar[4] = (-p6 + p10 + p7 - p11) * qa2 + (-p1 + p13 + p4 - p16) *
		    qb2 + (p2 - p14 - p3 + p15) * qab + (p5 - p9 - p8 + p12) *
		     qba;
	    fpar[5] = p1 - p5 - p9 + p13 - p2 + p6 + p10 - p14 - p3 + p7 +
		    p11 - p15 + p4 - p8 - p12 + p16;
	    fpar[6] = (-p6 + p10 + p7 - p11) * qa2 + (-p1 + p13 + p4 - p16) *
		    qb2 + (p2 - p14 - p3 + p15) * qba + (p5 - p9 - p8 + p12) *
		     qab;
	    fpar[7] = (p6 - p10 + p7 - p11) * qab * ra + (p1 - p13 + p4 - p16)
		     * qba * rb + (-p2 + p14 - p3 + p15) * ra + (-p5 + p9 -
		    p8 + p12) * rb;
	    fpar[8] = (-p2 + p6 + p10 - p14 + p3 - p7 - p11 + p15) * qab * ra
		    + (p1 - p5 - p9 + p13 - p4 + p8 + p12 - p16) * qba * rb;
	    fpar[9] = (-p5 + p9 + p6 - p10 + p7 - p11 - p8 + p12) * qab * ra
		    + (p1 - p13 - p2 + p14 - p3 + p15 + p4 - p16) * qba * rb;
	    fpar[10] = (p6 + p10 - p7 - p11) * qab * ra + (p1 + p13 - p4 -
		    p16) * qba * rb + (-p2 - p14 + p3 + p15) * rb + (-p5 - p9
		    + p8 + p12) * ra;
	    fpar[11] = (p2 - p6 - p10 + p14 + p3 - p7 - p11 + p15) * raa + (
		    -p1 + p5 + p9 - p13 - p4 + p8 + p12 - p16) * rbb;
	    fpar[12] = (p6 - p10 - p7 + p11) * raa * qa2 + (p1 - p13 - p4 +
		    p16) * rbb * qb2 + (-p5 + p9 - p2 + p14 + p3 - p15 + p8 -
		    p12) * rab;
	    fpar[13] = (p5 + p9 - p6 - p10 - p7 - p11 + p8 + p12) * raa + (
		    -p1 - p13 + p2 + p14 + p3 + p15 - p4 - p16) * rbb;
	    fpar[14] = (-p6 + p10 - p7 + p11) * qab * ra3 + (-p1 + p13 - p4 +
		    p16) * qba * rb3 + (p2 - p14 + p3 - p15) * rab * rb + (p5
		    - p9 + p8 - p12) * rab * ra;
	    fpar[15] = (-p6 - p10 + p7 + p11) * qab * ra3 + (-p1 - p13 + p4 +
		    p16) * qba * rb3 + (p2 + p14 - p3 - p15) * rab * ra + (p5
		    + p9 - p8 - p12) * rab * rb;
	    fpar[16] = (p6 + p10 + p7 + p11) * raa * raa + (p1 + p13 + p4 +
		    p16) * rbb * rbb + (-p5 - p9 - p2 - p14 - p3 - p15 - p8 -
		    p12) * raa * rbb;
	    for (i = 1; i <= 16; ++i) {
		fpar[i] *= rq;
	    }
	} else {
        dserror("TOO MANY GAUSS POINTS FOR EXTRAPOLATION");
	}
    }
/*---------------------------------- calculate function value F(r,s) ---*/
    if (igauss == 1) {
/*---------------------------------- constant function extrapolation ---*/
	*funval = fpar[1];
    } else if (igauss == 2) {
/*----------------- constant function extrapolation for TRIANGLES ------*/
	p1 = fpar[1];
	p2 = fpar[2];
	*funval = (p1+p2)/2.0;
    } else if (igauss == 3) {
/*------------------- linear function extrapolation for TRIANGLES ------*/
        /* assumes (ele->w1.e->nGP[1] == 1) type of Gauss-points */
	p1 = fpar[1];
	p2 = fpar[2];
	p3 = fpar[3];
	*funval = p2 + 2.0*(p2-p1)*(s-0.5) + 2.0*(p2-p3)*(r-0.5);
    } else if (igauss == 4) {
/*------------------------------------ linear function extrapolation ---*/
	p1 = fpar[1];
	p2 = fpar[2];
	p3 = fpar[3];
	p4 = fpar[4];
	*funval = p1 * r + p2 * r * s + p3 * s + p4;
/*--------------------------------- quadratic function extrapolation ---*/
    } else if (igauss == 9) {
	p1 = fpar[1];
	p2 = fpar[2];
	p3 = fpar[3];
	p4 = fpar[4];
	p5 = fpar[5];
	p6 = fpar[6];
	p7 = fpar[7];
	p8 = fpar[8];
	p9 = fpar[9];
	rs = r * s;
	*funval = p1 * rs * rs + p2 * rs * r + p3 * rs * s + p4 * r * r +
		p5 * rs + p6 * s * s + p7 * r + p8 * s + p9;
/*------------------------------------- cubic function extrapolation ---*/
    } else if (igauss == 16) {
	p1 = fpar[1];
	p2 = fpar[2];
	p3 = fpar[3];
	p4 = fpar[4];
	p5 = fpar[5];
	p6 = fpar[6];
	p7 = fpar[7];
	p8 = fpar[8];
	p9 = fpar[9];
	p10 = fpar[10];
	p11 = fpar[11];
	p12 = fpar[12];
	p13 = fpar[13];
	p14 = fpar[14];
	p15 = fpar[15];
	p16 = fpar[16];

	rs = r * s;
	rr = r * r;
	r3 = r * rr;
	ss = s * s;
	s3 = s * ss;
	*funval = p1 * r3 * s3 + p2 * r3 * ss + p3 * rr * s3 + p4 * r3 * s +
		p5 * rr * ss + p6 * r * s3 + p7 * r3 + p8 * rr * s + p9 *
		r * ss + p10 * s3 + p11 * rr + p12 * rs + p13 * ss + p14 * r
		+ p15 * s + p16;
    } else {
        dserror("TOO MANY GAUSS POINTS FOR EXTRAPOLATION");
    }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
    return ;
} /* w1recs */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
#endif
