#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | integration of nonlinear stiffness ke                                |
   and internal forces for shell8 element                m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8static_keug(ELEMENT   *ele,                         /* the element structure */
                   S8_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix      (NOT initialized!) */
                   double    *force,                      /* for internal forces (initialized!) */
                   int        kstep,                      /* actual step in nonlinear analysis */
                   int        init)                       /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*/
/* if force==NULL no internal forces are calculated                     */
/*----------------------------------------------------------------------*/
{
int                 i,j,k,l;                                /* some loopers */
double              sum;

int                 imass;                                  /* flag for calculating mass matrix */
double              density;                                /* density of the material */
double              facv,facw,facvw;                        /* variables for mass integration */
double             *thick;

int                 dof;
int                 nir,nis,nit;                            /* num GP in r/s/t direction */
int                 lr, ls, lt;                             /* loopers over GP */
int                 iel;                                    /* numnp to this element */
int                 nd;                                     /* ndofs to this element (=numdf*numnp) */
const int           numdf=NUMDOF_SHELL8;                    /* ndofs per node to this element */

double              e1,e2,e3;                               /*GP-coords*/
double              fac,facr,facs,fact;                     /* weights at GP */
double              xnu;                                    /* value of shell shifter */
double              weight;

double              condfac;                                /* sdc conditioning factor */
double              h2;                                     /* half nodal height */

double              detr0;                                  /* jacobian determinant in ref conf. at mid point*/
double              detc0;                                  /* jacobian determinant in cur conf. at mid point*/
double              detr;                                   /* jacobian determinant in ref conf. at gp */
double              detc;                                   /* jacobian determinant in cur conf. at gp */

double              h[3];                                   /* working array */
double              da;                                     /* area on mid surface */

double              stress[6], stress_r[12];                /* stress and stress resultants */
double              strain[6];                              /* strains */
double             *intforce;

static ARRAY        C_a;         static double **C;         /* material tensor */
static ARRAY        D_a;         static double **D;         /* material tensor integrated in thickness direction */

static ARRAY        a3r_a;       static double **a3r;       /* a3 in reference config (lenght h/2) */
static ARRAY        a3c_a;       static double **a3c;       /* a3 in current   config (lenght h/2 + disp) */

static ARRAY        a3kvpr_a;    static double **a3kvpr;    /* partiel derivatives of normal vector ref.config. */
static ARRAY        a3kvpc_a;    static double **a3kvpc;    /* partiel derivatives of normal vector cur.config. */

static ARRAY        xrefe_a;     static double **xrefe;     /* coords of midsurface in ref config */
static ARRAY        xcure_a;     static double **xcure;     /* coords of midsurface in cur condig */
                                        double **a3ref;     /* elements directors (lenght 1) */

static ARRAY        funct_a;     static double  *funct;     /* shape functions */
static ARRAY        deriv_a;     static double **deriv;     /* derivatives of shape functions */

/* mid surface basis vectors and metric tensors */
static ARRAY        akovr_a;     static double **akovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        akonr_a;     static double **akonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        amkovr_a;    static double **amkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        amkonr_a;    static double **amkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        akovc_a;     static double **akovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        akonc_a;     static double **akonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        amkovc_a;    static double **amkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        amkonc_a;    static double **amkonc;    /* kontravar.--------------"------------ current.config. */

static ARRAY        akovr0_a;    static double **akovr0;    /* kovariant basis vectors at mid point ref.config. */
static ARRAY        akonr0_a;    static double **akonr0;    /* kontravar.--------------"----------- ref.config. */
static ARRAY        amkovr0_a;   static double **amkovr0;   /* kovaraiant metric tensor at mid point ref.config. */
static ARRAY        amkonr0_a;   static double **amkonr0;   /* kontravar.--------------"------------ ref.config. */

static ARRAY        akovc0_a;    static double **akovc0;    /* kovariant basis vectors at mid point current.config. */
static ARRAY        akonc0_a;    static double **akonc0;    /* kontravar.--------------"----------- current.config. */
static ARRAY        amkovc0_a;   static double **amkovc0;   /* kovaraiant metric tensor at mid point current.config. */
static ARRAY        amkonc0_a;   static double **amkonc0;   /* kontravar.--------------"------------ current.config. */

/* shell body basis vectors and metric tensors */
static ARRAY        gkovr_a;     static double **gkovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        gkonr_a;     static double **gkonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        gmkovr_a;    static double **gmkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        gmkonr_a;    static double **gmkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        gkovc_a;     static double **gkovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        gkonc_a;     static double **gkonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        gmkovc_a;    static double **gmkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        gmkonc_a;    static double **gmkonc;    /* kontravar.--------------"------------ current.config. */

static ARRAY        bop_a;       static double **bop;       /* B-Operator for compatible strains */
                                 static double **estif;     /* element stiffness matrix ke */
                                 static double **emass;     /* element mass matrix */
static ARRAY        work_a;      static double **work;      /* working array to do Bt*D*B */

/* arrays for eas */
int                 nhyb;                                   /* scnd dim of P */
double              epsh[12];                               /* transformed eas strains */
double             *alfa;
double            **oldDtildinv;
double            **oldLt;
double             *oldRtild;
static ARRAY        P_a;         static double **P;         /* eas matrix M */
static ARRAY        transP_a;    static double **transP;    /* eas matrix M */
static ARRAY        T_a;         static double **T;         /* transformation matrix for eas */
static ARRAY        workeas_a;   static double **workeas;   /* eas working array */
static ARRAY        workeas2_a;  static double **workeas2;   /* eas working array */
static ARRAY        Lt_a;        static double **Lt;        /* eas matrix L transposed */
static ARRAY        Dtild_a;     static double **Dtild;     /* eas matrix Dtilde */
static ARRAY        Dtildinv_a;  static double **Dtildinv;  /* inverse of eas matrix Dtilde */
static ARRAY        Rtild_a;     static double  *Rtild;     /* eas part of internal forces */
static ARRAY        eashelp_a;   static double  *eashelp;   /* working vector for eas */


/* arrays for ANS */
int                 ansq;
const int           nsansmax=6;
int                 nsansq;                                  /* number of sampling points for ans */
double              xr1[6];                                  /* coordinates of collocation points for ANS */
double              xs1[6];
double              xr2[6];
double              xs2[6];                     
double              frq[6];
double              fsq[6];

static ARRAY       funct1q_a[6];  static double  *funct1q[6];    /* shape functions at collocation points */
static ARRAY       deriv1q_a[6];  static double **deriv1q[6];    /* derivation of these shape functions */

static ARRAY       akovr1q_a[6];  static double **akovr1q[6];
static ARRAY       akonr1q_a[6];  static double **akonr1q[6];
static ARRAY       amkovr1q_a[6]; static double **amkovr1q[6];
static ARRAY       amkonr1q_a[6]; static double **amkonr1q[6];
static ARRAY       a3kvpr1q_a[6]; static double **a3kvpr1q[6];

static ARRAY       akovc1q_a[6];  static double **akovc1q[6];
static ARRAY       akonc1q_a[6];  static double **akonc1q[6];
static ARRAY       amkovc1q_a[6]; static double **amkovc1q[6];
static ARRAY       amkonc1q_a[6]; static double **amkonc1q[6];
static ARRAY       a3kvpc1q_a[6]; static double **a3kvpc1q[6];


static ARRAY       funct2q_a[6];  static double  *funct2q[6];    /* shape functions at collocation points */
static ARRAY       deriv2q_a[6];  static double **deriv2q[6];    /* derivation of these shape functions */

static ARRAY       akovr2q_a[6];  static double **akovr2q[6];
static ARRAY       akonr2q_a[6];  static double **akonr2q[6];
static ARRAY       amkovr2q_a[6]; static double **amkovr2q[6];
static ARRAY       amkonr2q_a[6]; static double **amkonr2q[6];
static ARRAY       a3kvpr2q_a[6]; static double **a3kvpr2q[6];

static ARRAY       akovc2q_a[6];  static double **akovc2q[6];
static ARRAY       akonc2q_a[6];  static double **akonc2q[6];
static ARRAY       amkovc2q_a[6]; static double **amkovc2q[6];
static ARRAY       amkonc2q_a[6]; static double **amkonc2q[6];
static ARRAY       a3kvpc2q_a[6]; static double **a3kvpc2q[6];

#ifdef DEBUG 
dstrc_enter("s8static_keug");
#endif
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
if (init==1)
{
iel = MAXNOD_SHELL8; /* maximum number of nodes for this type of shell */

xrefe     = amdef("xrefe"  ,&xrefe_a,3,MAXNOD_SHELL8,"DA");       
xcure     = amdef("xcure"  ,&xcure_a,3,MAXNOD_SHELL8,"DA");       
a3r       = amdef("a3r"    ,&a3r_a,3,MAXNOD_SHELL8,"DA");         
a3c       = amdef("a3c"    ,&a3c_a,3,MAXNOD_SHELL8,"DA");         

a3kvpr    = amdef("a3kvpr" ,&a3kvpr_a,3,2,"DA");        
a3kvpc    = amdef("a3kvpc" ,&a3kvpc_a,3,2,"DA");        

funct     = amdef("funct"  ,&funct_a,MAXNOD_SHELL8,1,"DV");       
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL8,"DA");       

akovr0    = amdef("akovr0" ,&akovr0_a,3,3,"DA");        
akonr0    = amdef("akonr0" ,&akonr0_a,3,3,"DA");        
amkovr0   = amdef("amkovr0",&amkovr0_a,3,3,"DA");       
amkonr0   = amdef("amkonr0",&amkonr0_a,3,3,"DA");       

akovc0    = amdef("akovc0" ,&akovc0_a,3,3,"DA");        
akonc0    = amdef("akonc0" ,&akonc0_a,3,3,"DA");        
amkovc0   = amdef("amkovc0",&amkovc0_a,3,3,"DA");       
amkonc0   = amdef("amkonc0",&amkonc0_a,3,3,"DA");       

akovr     = amdef("akovr"  ,&akovr_a,3,3,"DA");         
akonr     = amdef("akonr"  ,&akonr_a,3,3,"DA");         
amkovr    = amdef("amkovr" ,&amkovr_a,3,3,"DA");        
amkonr    = amdef("amkonr" ,&amkonr_a,3,3,"DA");        

akovc     = amdef("akovc"  ,&akovc_a,3,3,"DA");         
akonc     = amdef("akonc"  ,&akonc_a,3,3,"DA");         
amkovc    = amdef("amkovc" ,&amkovc_a,3,3,"DA");        
amkonc    = amdef("amkonc" ,&amkonc_a,3,3,"DA");        

gkovr     = amdef("gkovr"  ,&gkovr_a,3,3,"DA");         
gkonr     = amdef("gkonr"  ,&gkonr_a,3,3,"DA");         
gmkovr    = amdef("gmkovr" ,&gmkovr_a,3,3,"DA");        
gmkonr    = amdef("gmkonr" ,&gmkonr_a,3,3,"DA");        

gkovc     = amdef("gkovc"  ,&gkovc_a,3,3,"DA");         
gkonc     = amdef("gkonc"  ,&gkonc_a,3,3,"DA");         
gmkovc    = amdef("gmkovc" ,&gmkovc_a,3,3,"DA");        
gmkonc    = amdef("gmkonc" ,&gmkonc_a,3,3,"DA");        

bop       = amdef("bop"    ,&bop_a ,12,(numdf*MAXNOD_SHELL8),"DA");
C         = amdef("C"      ,&C_a   ,6 ,6                    ,"DA");             
D         = amdef("D"      ,&D_a   ,12,12                   ,"DA");           
work      = amdef("work"   ,&work_a,12,(MAXNOD_SHELL8*numdf),"DA");  

/* for eas */
P         = amdef("P"      ,&P_a       ,12                   ,MAXHYB_SHELL8        ,"DA");                 
transP    = amdef("transP" ,&transP_a  ,12                   ,MAXHYB_SHELL8        ,"DA");       
T         = amdef("T"      ,&T_a       ,12                   ,12                   ,"DA");
workeas   = amdef("workeas",&workeas_a ,MAXHYB_SHELL8        ,(MAXNOD_SHELL8*numdf),"DA");
workeas2  = amdef("workeas",&workeas2_a,(MAXNOD_SHELL8*numdf),MAXHYB_SHELL8        ,"DA");
Lt        = amdef("Lt"     ,&Lt_a      ,MAXHYB_SHELL8        ,(MAXNOD_SHELL8*numdf),"DA");
Dtild     = amdef("Dtild"  ,&Dtild_a   ,MAXHYB_SHELL8        ,MAXHYB_SHELL8        ,"DA");
Dtildinv  = amdef("Dtildi" ,&Dtildinv_a,MAXHYB_SHELL8        ,MAXHYB_SHELL8        ,"DA");
Rtild     = amdef("Rtild"  ,&Rtild_a   ,MAXHYB_SHELL8        ,1                    ,"DV");
eashelp   = amdef("eashelp",&eashelp_a ,MAXHYB_SHELL8        ,1                    ,"DV");

/* for ans */
for (i=0; i<nsansmax; i++) funct1q[i]  = amdef("funct1q",&(funct1q_a[i]),MAXNOD_SHELL8,1,"DV");
for (i=0; i<nsansmax; i++) deriv1q[i]  = amdef("deriv1q",&(deriv1q_a[i]),2,MAXNOD_SHELL8,"DA");
for (i=0; i<nsansmax; i++) akovr1q[i]  = amdef("akovr1q" ,&(akovr1q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) akonr1q[i]  = amdef("akonr1q" ,&(akonr1q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) amkovr1q[i] = amdef("amkovr1q",&(amkovr1q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) amkonr1q[i] = amdef("amkonr1q",&(amkonr1q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) a3kvpr1q[i] = amdef("a3kvpr1q",&(a3kvpr1q_a[i]),3,2,"DA");
for (i=0; i<nsansmax; i++) akovc1q[i]  = amdef("akovc1q" ,&(akovc1q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) akonc1q[i]  = amdef("akonc1q" ,&(akonc1q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) amkovc1q[i] = amdef("amkovc1q",&(amkovc1q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) amkonc1q[i] = amdef("amkonc1q",&(amkonc1q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) a3kvpc1q[i] = amdef("a3kvpc1q",&(a3kvpc1q_a[i]),3,2,"DA");

for (i=0; i<nsansmax; i++) funct2q[i]  = amdef("funct2q",&(funct2q_a[i]),MAXNOD_SHELL8,1,"DV");
for (i=0; i<nsansmax; i++) deriv2q[i]  = amdef("deriv2q",&(deriv2q_a[i]),2,MAXNOD_SHELL8,"DA");
for (i=0; i<nsansmax; i++) akovr2q[i]  = amdef("akovr2q" ,&(akovr2q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) akonr2q[i]  = amdef("akonr2q" ,&(akonr2q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) amkovr2q[i] = amdef("amkovr2q",&(amkovr2q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) amkonr2q[i] = amdef("amkonr2q",&(amkonr2q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) a3kvpr2q[i] = amdef("a3kvpr2q",&(a3kvpr2q_a[i]),3,2,"DA");
for (i=0; i<nsansmax; i++) akovc2q[i]  = amdef("akovc2q" ,&(akovc2q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) akonc2q[i]  = amdef("akonc2q" ,&(akonc2q_a[i]) ,3,3,"DA");
for (i=0; i<nsansmax; i++) amkovc2q[i] = amdef("amkovc2q",&(amkovc2q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) amkonc2q[i] = amdef("amkonc2q",&(amkonc2q_a[i]),3,3,"DA");
for (i=0; i<nsansmax; i++) a3kvpc2q[i] = amdef("a3kvpc2q",&(a3kvpc2q_a[i]),3,2,"DA");

goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase                                                         */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
amdel(&xrefe_a);
amdel(&xcure_a);
amdel(&a3r_a);
amdel(&a3c_a);

amdel(&a3kvpr_a);
amdel(&a3kvpc_a);

amdel(&funct_a);
amdel(&deriv_a);

amdel(&akovr0_a);   
amdel(&akonr0_a);   
amdel(&amkovr0_a);  
amdel(&amkonr0_a);  

amdel(&akovc0_a);   
amdel(&akonc0_a);   
amdel(&amkovc0_a);  
amdel(&amkonc0_a);  

amdel(&akovr_a);   
amdel(&akonr_a);   
amdel(&amkovr_a);  
amdel(&amkonr_a);  

amdel(&akovc_a);   
amdel(&akonc_a);   
amdel(&amkovc_a);  
amdel(&amkonc_a);

amdel(&gkovr_a);   
amdel(&gkonr_a);   
amdel(&gmkovr_a);  
amdel(&gmkonr_a);  

amdel(&gkovc_a);   
amdel(&gkonc_a);   
amdel(&gmkovc_a);  
amdel(&gmkonc_a);

amdel(&bop_a);
amdel(&C_a);
amdel(&D_a);
amdel(&work_a);

amdel(&P_a);
amdel(&transP_a);
amdel(&T_a);
amdel(&workeas_a);
amdel(&workeas2_a);
amdel(&Lt_a);
amdel(&Dtild_a);
amdel(&Dtildinv_a);
amdel(&Rtild_a);
amdel(&eashelp_a);

/* uninit von ans feldern fehlt noch */

goto end;  
}
/*----------------------------------------------------------------------*/
/* calculation phase                                                    */
/*----------------------------------------------------------------------*/
/*-------------------------------------------- init the gaussian points */
s8intg(ele,data,0);
/*------------------------------------------------------- check for eas */
nhyb=ele->e.s8->nhyb;
if (nhyb>0) 
{
   amzero(&P_a);
   amzero(&T_a);
   amzero(&Lt_a);
   amzero(&Dtild_a);
   amzero(&Rtild_a);
   /*--------------------------------------- update of eas strains alfa */
   /*----------------------------------------- resize alfa if necessary */
   if (kstep >= ele->e.s8->alfa.fdim)
   amredef(&(ele->e.s8->alfa),ele->e.s8->alfa.fdim+5,nhyb,"DA");
   /*--------------------------------------- set pointer to actual alfa */
   alfa = ele->e.s8->alfa.a.da[kstep];
   /*-------------------- set pointer to storage of old Dtildinv and Lt */
   oldDtildinv = ele->e.s8->Dtildinv.a.da;
   oldLt       = ele->e.s8->Lt.a.da;
   oldRtild    = ele->e.s8->Rtilde.a.dv;
   /*---------------- set number of displacement - dofs to this element */
   iel = ele->numnp;
   nd  = iel*NUMDOF_SHELL8;
   /*------------------- make multiplication eashelp = Lt * disp[kstep] */
   for (i=0; i<nhyb; i++)
   {
      sum=0.0;
      for (j=0; j<iel; j++)
      for (k=0; k<numdf; k++)
      {
         l = j*numdf+k;
         sum += oldLt[i][l]*ele->node[j]->sol_residual.a.da[0][k];
      }
      eashelp[i] = sum;
   }
   /*---------------------------------------- add old Rtilde to eashelp */
   for (i=0; i<nhyb; i++) eashelp[i] += oldRtild[i];
   /*----------------- make multiplication alfa -= olDtildinv * eashelp */
   math_matvecdense(alfa,oldDtildinv,eashelp,nhyb,nhyb,1,-1.0);

}
/*------------------------------------ check calculation of mass matrix */
if (emass_global) 
{
   imass = 1;
   amzero(emass_global);
   emass = emass_global->a.da;
   s8_getdensity(mat,&density);
   thick = ele->e.s8->thick_node.a.dv;
} 
else 
{
   imass   = 0;
   emass   = NULL;
   density = 0.0; 
   thick   = NULL;
}   
/*----------------- some of the fields have to be reinitialized to zero */
amzero(&D_a);
amzero(estif_global);
estif     = estif_global->a.da;
/*----------------------------------------------- integrationsparameter */
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
nit     = ele->e.s8->nGP[2];
iel     = ele->numnp;
nd      = iel*NUMDOF_SHELL8;
condfac = ele->e.s8->sdc;
a3ref   = ele->e.s8->a3ref.a.da;
          amzero(&(ele->e.s8->intforce));
intforce= ele->e.s8->intforce.a.dv;
/*----------------------------------------------------- geometry update */
for (k=0; k<iel; k++)
{
   h2 = ele->e.s8->thick_node.a.dv[k];
   h2 /= 2.0;
   h2 *= condfac;
   
   a3r[0][k] = a3ref[0][k] * h2;
   a3r[1][k] = a3ref[1][k] * h2;
   a3r[2][k] = a3ref[2][k] * h2;

   xrefe[0][k] = ele->node[k]->x[0];
   xrefe[1][k] = ele->node[k]->x[1];
   xrefe[2][k] = ele->node[k]->x[2];

   xcure[0][k] = xrefe[0][k] + ele->node[k]->sol.a.da[0][0];
   xcure[1][k] = xrefe[1][k] + ele->node[k]->sol.a.da[0][1];
   xcure[2][k] = xrefe[2][k] + ele->node[k]->sol.a.da[0][2];
 
   a3c[0][k] = a3r[0][k]     + ele->node[k]->sol.a.da[0][3];
   a3c[1][k] = a3r[1][k]     + ele->node[k]->sol.a.da[0][4];
   a3c[2][k] = a3r[2][k]     + ele->node[k]->sol.a.da[0][5];
}
/*============ metric and shape functions at collocation points (ANS=1) */
/*----------------------------------------------------- 4-noded element */
if (ele->e.s8->ans==1 || ele->e.s8->ans==3)/*------------ querschub_ans */
{
                 ansq=1;
   if (iel==4) nsansq=2;
   if (iel==9) nsansq=6;
   s8_ans_colloqpoints(nsansq,iel,ele->e.s8->ans,ele->distyp,
                       xr1,xs1,xr2,xs2,
                       funct1q,deriv1q,funct2q,deriv2q,
                       xrefe,a3r,xcure,a3c,
                       akovr1q,akonr1q,amkovr1q,amkonr1q,a3kvpr1q,
                       akovc1q,akonc1q,amkovc1q,amkonc1q,a3kvpc1q,
                       akovr2q,akonr2q,amkovr2q,amkonr2q,a3kvpr2q,
                       akovc2q,akonc2q,amkovc2q,amkonc2q,a3kvpc2q,
                       &detr,&detc);
}
else
{
   ansq=0;
   nsansq=0;
}
/*=============================== metric of element mid point (for eas) */
if (nhyb>0)
{
   s8_funct_deriv(funct,deriv,0.0,0.0,ele->distyp,1);

   s8_tmtr(xrefe,a3r,0.0,akovr0,akonr0,amkovr0,amkonr0,&detr0,
              funct,deriv,iel,condfac,0);

   s8_tmtr(xcure,a3c,0.0,akovc0,akonc0,amkovc0,amkonc0,&detc0,
              funct,deriv,iel,condfac,0);
}
/*=================================================== integration loops */
for (lr=0; lr<nir; lr++)
{
   /*================================== gaussian point and weight at it */
   e1   = data->xgpr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      /*=============================== gaussian point and weight at it */
      e2   = data->xgps[ls];
      facs = data->wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      s8_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*----------------------------- shape functions for querschub-ans */
      if (ansq==1) s8_ansq_funct(frq,fsq,e1,e2,iel,nsansq);
      /*-------- init mid surface material tensor and stress resultants */
      amzero(&D_a);
      for (i=0; i<12; i++) stress_r[i]=0.0;
      /*------------------------------------ init mass matrix variables */
      if (imass)
      {
         facv  = 0.0;
         facw  = 0.0;
         facvw = 0.0;
      }
      /*------------------------------------- metrics at gaussian point */
      s8_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,&detr,
                 funct,deriv,iel,a3kvpr,0);
      s8_tvmr(xcure,a3c,akovc,akonc,amkovc,amkonc,&detc,
                 funct,deriv,iel,a3kvpc,0);
      /*------------------------- make h as cross product in ref config.
                                       to get area da on shell mid surf */
      h[0] = akovr[1][0]*akovr[2][1] - akovr[2][0]*akovr[1][1];
      h[1] = akovr[2][0]*akovr[0][1] - akovr[0][0]*akovr[2][1];
      h[2] = akovr[0][0]*akovr[1][1] - akovr[1][0]*akovr[0][1];
      /*------------------------------------- make director unit lenght 
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);
      /*--------------------------------------- make eas if switched on */
      if (nhyb>0)
      {
         /*------------------- make shape functions for incomp. strains */
         s8_eas(nhyb,e1,e2,iel,ele->e.s8->eas,P);
         /*-------------------- transform basis of Eij to Gausian point */
         s8_transeas(P,transP,T,akovr,akonr0,detr,detr0,nhyb);
         /*------------------------ transform strains to Gaussian point */
         math_matvecdense(epsh,transP,alfa,12,nhyb,0,1.0);
      }
      /*------------------------ make B-operator for compatible strains */
      s8_tvbo(e1,e2,bop,funct,deriv,iel,numdf,akovc,a3kvpc,nsansq);
      /*-------------------------------------- modifications due to ans */
      if (ansq==1)
      s8_ans_bbar_q(bop,frq,fsq,funct1q,funct2q,deriv1q,deriv2q,
                    akovc1q,akovc2q,a3kvpc1q,a3kvpc2q,iel,numdf,nsansq);
      /*============================== loop GP in thickness direction t */
      for (lt=0; lt<nit; lt++)
      {
         /*---------------------------- gaussian point and weight at it */
         e3   = data->xgpt[lt];
         fact = data->wgtt[lt];
         /*-------------------- basis vectors and metrics at shell body */ 
         s8_tmtr(xrefe,a3r,e3,gkovr,gkonr,gmkovr,gmkonr,&detr,
                    funct,deriv,iel,condfac,0);

         s8_tmtr(xcure,a3c,e3,gkovc,gkonc,gmkovc,gmkonc,&detc,
                    funct,deriv,iel,condfac,0);
         /*--------------------------------- metric at gp in shell body */     
         if (ansq==0)
         s8_tvhe(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,&detr,&detc,
                 amkovc,amkovr,akovc,akovr,a3kvpc,a3kvpr,e3);     
         /*- modifications to metric of shell body due to querschub-ans */
         else
         s8_ans_tvhe_q(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,amkovc,amkovr,
                       akovc,akovr,a3kvpc,a3kvpr,&detr,&detc,
                       amkovr1q,amkovc1q,akovr1q,akovc1q,a3kvpr1q,a3kvpc1q,
                       amkovr2q,amkovc2q,akovr2q,akovc2q,a3kvpr2q,a3kvpc2q,
                       frq,fsq,e3,nsansq,iel);
         /*----------- calc shell shifter and put it in the weight fact */
         xnu   = (1.0/condfac)*(detr/da);
         fact *= xnu; 
         /*----------------------- change to current metrics due to eas */
         if (nhyb>0) s8_vthv(gmkovc,gmkonc,epsh,&detc,e3);
         /*------------------------------------------ call material law */
         s8_tmat(ele,mat,stress,strain,C,gmkovc,gmkonc,gmkovr,gmkonr,
                    gkovc,gkonc,gkovr,gkonr,detc,detr,e3,0);
         /*---------------- do thickness integration of material tensor */           
         s8_tvma(D,C,stress,stress_r,e3,fact,condfac);
         /*-------------------------- mass matrix thickness integration */
         if (imass)
         {
            facv  += data->wgtt[lt] * detr;
            facw  += data->wgtt[lt] * detr * DSQR(e3);
            facvw += data->wgtt[lt] * detr * e3;
         }
      }/*========================================== end of loop over lt */
      /*------------ product of all weights and jacobian of mid surface */
      weight = facr*facs*da;
      /*----------------------------------- elastic stiffness matrix ke */
      s8_BtDB(estif,bop,D,iel,numdf,weight,work);
      /*--------------------------------- geometric stiffness matrix kg */
      if (ansq==0)
      s8_tvkg(estif,stress_r,funct,deriv,numdf,iel,weight,e1,e2);
      else
      s8_ans_tvkg(estif,stress_r,funct,deriv,numdf,iel,weight,e1,e2,
                  frq,fsq,funct1q,funct2q,deriv1q,deriv2q,ansq,nsansq);
      /*-------------------------------- calculation of internal forces */
      if (force) s8_intforce(intforce,stress_r,bop,iel,numdf,12,weight);
      /*------------- mass matrix : gaussian point on shell mid surface */
      if (imass)
      {
         fac    = facr * facs * density;
         facv  *= fac;
         facw  *= fac;
         facvw *= fac;
         s8_tmas(funct,thick,emass,iel,numdf,facv,facw,facvw);
      }
      /*----------------------------------- integration of eas matrices */
      if (nhyb>0)
      {
         /*=============================================================*/
         /*  Ltrans(nhyb,nd) = Mtrans(nhyb,12) * D(12,12) * B(12,nd)
         /*=============================================================*/
         /*----------------------------------------------------- DB=D*B */
         math_matmatdense(workeas,D,bop,12,12,nd,0,0.0);
         /*----------------------------------- Ltransposed = Mt * D * B */
         math_mattrnmatdense(Lt,transP,workeas,nhyb,12,nd,1,weight);
         /*=============================================================*/
         /*  Dtilde(nhyb,nhyb) = Mtrans(nhyb,12) * D(12,12) * M(12,nhyb)
         /*=============================================================*/
         /*----------------------------------------------------DM = D*M */
         math_matmatdense(workeas,D,transP,12,12,nhyb,0,0.0);
         /*-------------------------------------------- Dtilde = Mt*D*M */
         math_mattrnmatdense(Dtild,transP,workeas,nhyb,12,nhyb,1,weight);
         /*=============================================================*/
         /*  Rtilde(nhyb) = Mtrans(nhyb,12) * Forces(12)
         /*=============================================================*/
         /*------------------------- eas part of internal forces Rtilde */
         math_mattrnvecdense(Rtild,transP,stress_r,nhyb,12,1,weight);
      }
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */
/*----------------- make modifications to stiffness matrices due to eas */
if (nhyb>0)
{
   /*------------------------------------ make inverse of matrix Dtilde */
   amcopy(&Dtild_a,&Dtildinv_a);
   math_sym_inv(Dtildinv,nhyb);
   /*===================================================================*/
   /* estif(nd,nd) = estif(nd,nd) - Ltrans(nhyb,nd) * Dtilde^-1(nhyb,nhyb) * L(nd,nhyb)
   /*===================================================================*/
   /*------------------------------------------- make Ltrans * Dtildinv */
   math_mattrnmatdense(workeas2,Lt,Dtildinv,nd,nhyb,nhyb,0,0.0);
   /*---------------------------------- make estif -= Lt * Dtildinv * L */
   math_matmatdense(estif,workeas2,Lt,nd,nhyb,nd,1,-1.0);
   /*===================================================================*/
   /* R(nd) = R(nd) - Ltrans(nhyb,nd) * Dtilde^-1(nhyb,nhyb) * Rtilde(nhyb)
   /*===================================================================*/
   /*--------------------------- make intforce -= Lt * Dtildinv * Rtild */
   math_matvecdense(intforce,workeas2,Rtild,nd,nhyb,1,-1.0);
   /*------------------------------------------ put Dtildinv to storage */
   for (i=0; i<nhyb; i++)
   for (j=0; j<nhyb; j++) ele->e.s8->Dtildinv.a.da[i][j] = Dtildinv[i][j];
   /*------------------------------------------------ put Lt to storage */
   for (i=0; i<nhyb; i++)
   for (j=0; j<nd; j++) ele->e.s8->Lt.a.da[i][j] = Lt[i][j];
   /*-------------------------------------------- put Rtilde to storage */
   for (i=0; i<nhyb; i++) ele->e.s8->Rtilde.a.dv[i] = Rtild[i];
   /*-------------------------------------------------------------- end */
}
/*- add internal forces to global vector, if a global vector was passed */
/*                                                      to this routine */
for (i=0; i<ele->numnp; i++)
{
   for (j=0; j<ele->node[i]->numdf; j++)
   {
      force[i*numdf+j] += intforce[i*numdf+j];
   }
}
/*----------------------------------------------------------------- end */
/*-------------------------------------- just to test symmetry of estif */
/*
for (i=0; i<nd; i++)
{
   for (j=0; j<nd; j++)
   {
      if (i==j); 
      else
      {
         sum = emass[i][j] - emass[j][i];
         sum = FABS(sum);
         if (sum > EPS11)
         {
            printf("emass[%d][%d]=%f emass[%d][%d]=%f diff=%20.10E\n",
            i,j,emass[i][j],
            j,i,emass[j][i],
            sum);
         }
      }
   }
}*/
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8static_keug */
