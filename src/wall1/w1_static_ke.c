/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1static_ke' which forms the linear stiffness
       ke for wall1 element
 contains the routine 'w1fi' which evaluates the element forces

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | integration of linear stiffness ke for wall1 element      al 9/01    |
 *----------------------------------------------------------------------*/
void w1static_ke(ELEMENT   *ele,
                 W1_DATA   *data,
                 MATERIAL  *mat,
                 ARRAY     *estif_global,
                 ARRAY     *emass_global, /* global vector for mass */
                 DOUBLE    *force,  /* global vector for internal forces (initialized!) */
                 INT        init)
{
INT                 i,a,b;            /* some loopers */
INT                 nir=0;          /* num GP in r/s direction */
INT                 nis=0;          /* num GP in r/s direction */
INT                 lr, ls;           /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
INT                 ip;
INT                 intc;      /* "integration case" for tri-element     */
INT                 lanz, maxreb;
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 0;/* controls evaluation of new stresses    */
const INT           numdf  = 2;
const INT           numeps = 3;

DOUBLE              fac,facm;         /* integration factors            */
DOUBLE              stifac;           /* thickness                      */
DOUBLE              e1=0.0;            /* GP-coords                      */
DOUBLE              e2=0.0;            /* GP-coords                      */
DOUBLE              facr=0.0;        /* weights at GP                  */
DOUBLE              facs=0.0;        /* weights at GP                  */
DOUBLE              density;          /* for calculation of mass matrix */
INT                 imass;            /* flag for calc of mass matrix   */

static ARRAY    D_a;         /* material tensor */
static DOUBLE **D;
static ARRAY    funct_a;     /* shape functions */
static DOUBLE  *funct;
static ARRAY    deriv_a;     /* derivatives of shape functions */
static DOUBLE **deriv;
static ARRAY    xjm_a;       /* jacobian matrix */
static DOUBLE **xjm;
static ARRAY    xjm0_a;      /* jacobian matrix at r = s = 0 */
static DOUBLE **xjm0;
static ARRAY    bop_a;       /* B-operator */
static DOUBLE **bop;
static ARRAY    gop_a;       /* incomp modes: strain = Bop*d + Gop*alpha */
static DOUBLE  *gop;
static ARRAY    alpha_a;     /* incomp modes: alpha = vector of element internal dof*/
static DOUBLE  *alpha;
static ARRAY    knc_a;       /* incomp modes: knc = BT C G (8*4)     */
static DOUBLE **knc;
static ARRAY    knn_a;       /* incomp modes: knn = GT C G (4*4)     */
static DOUBLE **knn;
static ARRAY    knninv_a;    /* incomp modes:knninv = inverse(knn)   */
static DOUBLE **knninv;
static ARRAY    deltak_a;    /* incomp modes:deltak=knc*knninv*kcn   */
static DOUBLE **deltak;
static ARRAY    fintn_a;     /* incomp modes: fintn = GT sig (4*1)   */
static DOUBLE  *fintn;
static ARRAY    deltaf_a;    /* incomp modes:deltaf=knc*knninv*fintn */
static DOUBLE  *deltaf;

static DOUBLE **estif;       /* element stiffness matrix ke */
static DOUBLE **emass;       /* mass matrix */

DOUBLE F[4];                 /* stress */

DOUBLE det,det0;             /* Jacobi-det, Jacobi-det at r=s=0 */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1static_ke");
#endif
/*------------------------------------------------- some working arrays */
istore = 0;

if (init==1)
{
funct     = amdef("funct"  ,&funct_a ,MAXNOD_WALL1,1 ,"DV");
deriv     = amdef("deriv"  ,&deriv_a ,2,MAXNOD_WALL1 ,"DA");
D         = amdef("D"      ,&D_a     ,6,6             ,"DA");
xjm       = amdef("xjm"    ,&xjm_a   ,numdf,numdf     ,"DA");
xjm0      = amdef("xjm0"   ,&xjm0_a  ,numdf,numdf     ,"DA");
bop       = amdef("bop"    ,&bop_a   ,numeps,(numdf*MAXNOD_WALL1),"DA");
gop       = amdef("gop"    ,&gop_a   ,4,1,"DV");
alpha     = amdef("alpha"  ,&alpha_a ,4,1,"DV");
knc       = amdef("knc"    ,&knc_a   ,8,4,"DA");
knn       = amdef("knn"    ,&knn_a   ,4,4,"DA");
knninv    = amdef("knninv" ,&knninv_a,4,4,"DA");
deltak    = amdef("deltak" ,&deltak_a,8,8,"DA");
fintn     = amdef("fintn"  ,&fintn_a ,4,1,"DV");
deltaf    = amdef("deltaf" ,&deltaf_a,8,1,"DV");

amzero(&D_a);
goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
   amdel(&funct_a);
   amdel(&deriv_a);
   amdel(&D_a);
   amdel(&xjm_a);
   amdel(&xjm0_a);
   amdel(&bop_a);
   amdel(&gop_a);
   amdel(&alpha_a);
   amdel(&knc_a);
   amdel(&knn_a);
   amdel(&knninv_a);
   amdel(&deltak_a);
   amdel(&deltaf_a);

   goto end;
}

else if(init==2)
{
  istore = 1;
}
/*------------------------------------------- integration parameters ---*/
w1intg(ele,data,1);
/*------------------------------------ check calculation of mass matrix */
if (emass_global)
{
   imass = 1;
   amzero(emass_global);
   emass = emass_global->a.da;
   w1_getdensity(mat,&density);
}
else
{
   imass   = 0;
   emass   = NULL;
   density = 0.0;
}
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
iel     = ele->numnp;
nd      = numdf * iel;
/*----------------------------------------------------------------------*/
stifac = ele->e.w1->thick;
/*----------------------------------------------------------------------*/
if(ele->e.w1->modeltype == incomp_mode)
{
 if(iel != 4)
 {
 dserror("methode of incompatible modes only for 4-node element");
 }
 if(ele->e.w1->wtype != plane_stress)
 {
 dserror("methode of incompatible modes only for Plane Stress");
 }
 /*------------------ shape functions and their derivatives at r,s=0 ---*/
 w1_funct_deriv(funct,deriv,0,0,ele->distyp,1);
 /*-------------------------------- compute jacobian matrix at r,s=0 ---*/
 w1_jaco (deriv,xjm0,&det0,ele,iel);
 /*---------------------------------------------------------------------*/
 amzero(&alpha_a);
 if(mat->mattyp != m_stvenant)
 {
  /*-------------------------------------- update internal dof alpha ---*/
  w1_updalpha(alpha,ele,knc,knninv,fintn,istore);
  amzero(&fintn_a);
  amzero(&deltaf_a);
 }
 /*---------------------------------------------------------------------*/
 amzero(&knc_a);
 amzero(&knn_a);
 amzero(&deltak_a);
}
/*------------------------------- loop concrete reinforcement steel ----*/
if(mat->mattyp==m_pl_epc) maxreb = mat->m.pl_epc->maxreb;
else                      maxreb = 0;

for (lanz=0; lanz<maxreb+1; lanz++)
{
/*------- get integraton data ---------------------------------------- */
switch (ele->distyp)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   nir = ele->e.w1->nGP[0];
   nis = ele->e.w1->nGP[1];
break;
case tri3: case tri6:  /* --> tri - element */
   nir  = ele->e.w1->nGP[0];
   nis  = 1;
   intc = ele->e.w1->nGP[1]-1;
break;
default:
   dserror("ele->distyp unknown! in 'w1_statik_ke.c' ");
} /* end switch(ele->distyp) */
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
       facr = data->wgtr[lr];
       e2   = data->xgss[ls];
       facs = data->wgts[ls];
      break;
      case tri3: case tri6:  /* --> tri - element */
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
      break;
      default:
         dserror("ele->distyp unknown!");
      } /* end switch(ele->distyp) */
      ip++;
      /*------------------------- shape functions and their derivatives */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------------------------------------ compute jacobian matrix ---*/
      w1_jaco (deriv,xjm,&det,ele,iel);
      /*--------------------------- thickness rebar  (input rebar %) ---*/
      if(lanz>0)
      {
        stifac = mat->m.pl_epc->reb_area[lanz-1];
      }
      /*------------------------------------ integration factor  -------*/
      fac = facr * facs * det * stifac;
      /*------------------------------ compute mass matrix if imass=1---*/
      if (imass == 1)
      {
       facm = fac * density;
       for (a=0; a<iel; a++)
       {
        for (b=0; b<iel; b++)
        {
         emass[2*a][2*b]     += facm * funct[a] * funct[b]; /* a,b even */
         emass[2*a+1][2*b+1] += facm * funct[a] * funct[b]; /* a,b odd  */
        }
       }
      }
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
      if(lanz==0)
      {
        w1_call_mat(ele,mat,ele->e.w1->wtype,bop,gop,alpha,ip, F, D,istore,newval);
      }
      else
      {
        w1_mat_rebar(ele,mat,bop,NULL,NULL,xjm,F,D,ip,lanz,istore);
      }
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*-------------------------------- element stiffness matrix ke ---*/
       w1_keku(estif,bop,D,fac,nd,numeps);
      /*------------- additional stiffness parts due to incom. modes ---*/
       if(ele->e.w1->modeltype == incomp_mode)
       {
        w1_knc(knc,bop,gop,D,fac);
        w1_knn(knn,gop,D,fac);
       }
      /*--------------- nodal forces fi from integration of stresses ---*/
       if (force)
       {
        w1fi (F,fac,bop,nd,force);
        if(ele->e.w1->modeltype == incomp_mode)
        {
         w1_fintn(F,fac,gop,fintn);
        }
      }
     }
   }/*============================================= end of loop over ls */
}/*================================================ end of loop over lr */
}/*------------------------------ loop concrete reinforcement steel ----*/

/*----------------- static condensation of internal dof for inc. modes--*/
if(ele->e.w1->modeltype == incomp_mode  && istore==0)
{
 w1_inverse_matrix(4,knn,knninv);   /*- knn is destroyed afterwards! ---*/
 if (force)    w1_stat_cond(knninv,knc,deltak,fintn,deltaf,ele);
 else          w1_stat_cond(knninv,knc,deltak,NULL,NULL,ele);
 amadd(estif_global,&deltak_a,-1.0,0);
 if (force)
 {
  for(i=0;i<8;i++)  force[i] -= deltaf[i];
 }
}
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1static_ke */
/*----------------------------------------------------------------------*
 | evaluates element forces                              al    9/01     |
 *----------------------------------------------------------------------*/
void w1fi( DOUBLE  *F,
           DOUBLE   fac,
           DOUBLE **bop,
           INT      nd,
           DOUBLE  *fie)
{
/*----------------------------------------------------------------------*/
INT i,j;
DOUBLE tau11, tau12, tau22, tau33;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1fi");
#endif
/*----------------------------------------------------------------------*/
  tau11 = F[0]*fac;
  tau22 = F[1]*fac;
  tau12 = F[2]*fac;
  tau33 = F[3]*fac;
/*----------------------------- updated lagrange or geometric linear ---*/
  for (j=1; j<nd; j+=2)
  {
    i=j-1;
    fie[i]+=bop[0][i]*tau11 + bop[2][i]*tau12;
    fie[j]+=bop[1][j]*tau22 + bop[2][j]*tau12;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1fi */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

#endif /*D_WALL1*/
