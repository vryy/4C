/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1static_ke' which forms the linear stiffness
       ke for wall1 element
 contains the routine 'w1fi' which evaluates the element forces

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
                 double    *force,  /* global vector for internal forces (initialized!) */
                 int        init)
{
int                 i,j,k,a,b;        /* some loopers */
int                 nir,nis;          /* num GP in r/s/t direction */
int                 lr, ls;           /* loopers over GP */
int                 iel;              /* numnp to this element */
int                 dof;
int                 nd;
int                 ip;
int                 lanz, maxreb;
int                 istore = 0;/* controls storing of new stresses to wa */
int                 newval = 0;/* controls evaluation of new stresses    */
const int           numdf  = 2;
const int           numeps = 3;

double              fac,facm;         /* integration factors            */
double              stifac;           /* thickness                      */
double              e1,e2;            /* GP-coords                      */
double              facr,facs;        /* weights at GP                  */
double              weight;
double              density;          /* for calculation of mass matrix */       
int                 imass;            /* flag for calc of mass matrix   */

static ARRAY    D_a;         /* material tensor */     
static double **D;         
static ARRAY    funct_a;     /* shape functions */    
static double  *funct;     
static ARRAY    deriv_a;     /* derivatives of shape functions */   
static double **deriv;     
static ARRAY    xjm_a;       /* jacobian matrix */     
static double **xjm;         
static ARRAY    xjm0_a;      /* jacobian matrix at r = s = 0 */     
static double **xjm0;         
static ARRAY    bop_a;       /* B-operator */   
static double **bop;       
static ARRAY    gop_a;       /* incomp modes: strain = Bop*d + Gop*alpha */   
static double  *gop;       
static ARRAY    alpha_a;     /* incomp modes: alpha = vector of element internal dof*/   
static double  *alpha;       
static ARRAY    knc_a;       /* incomp modes: knc = BT C G (8*4)     */   
static double **knc;
static ARRAY    knn_a;       /* incomp modes: knn = GT C G (4*4)     */   
static double **knn;
static ARRAY    knninv_a;    /* incomp modes:knninv = inverse(knn)   */   
static double **knninv;
static ARRAY    deltak_a;    /* incomp modes:deltak=knc*knninv*kcn   */   
static double **deltak;
static ARRAY    fintn_a;     /* incomp modes: fintn = GT sig (4*1)   */   
static double  *fintn;
static ARRAY    deltaf_a;    /* incomp modes:deltaf=knc*knninv*fintn */   
static double  *deltaf;
       
static double **estif;       /* element stiffness matrix ke */
static double **emass;       /* mass matrix */

double F[4];                 /* stress */
double fie[18];              /* internal force */

double det,det0;             /* Jacobi-det, Jacobi-det at r=s=0 */

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
for (i=0; i<18; i++) fie[i] = 0.0;
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
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
 w1_jaco (funct,deriv,xjm0,&det0,ele,iel);                         
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
/*================================================ integration loops ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
   /*=============================== gaussian point and weight at it ===*/
   e1   = data->xgrr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      ip++;
      /*============================ gaussian point and weight at it ===*/
      e2   = data->xgss[ls];
      facs = data->wgts[ls];
      /*------------------------- shape functions and their derivatives */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*------------------------------------ compute jacobian matrix ---*/       
      w1_jaco (funct,deriv,xjm,&det,ele,iel);                         
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
        w1_call_mat(ele,mat,ele->e.w1->wtype,bop,gop,alpha,xjm,ip, F, D,istore,newval);
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
void w1fi( double  *F,
           double   fac,
           double **bop,
           int      nd,
           double  *fie)
{
/*----------------------------------------------------------------------*/
int i,j;
double tau11, tau12, tau21, tau22, tau33;
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
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
