/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9static_keug: which caluclates the element stiffness matrix 
                 (geom. nonliner) elementtechnology like ANS for 
                 'Querschub' and EAS for a different locking type 
                 effects are implemented


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
\brief vector of multilayer material law

<pre>                                                            sh 10/02
This structure struct _MULTIMAT  *multimat is defined in global_control.c
and the type is in materials.h                                                  
It holds all information about the layered section data
</pre>
*----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;



/*!----------------------------------------------------------------------
\brief integration of linear stiffness ke for shell9 element                                      

<pre>                     m.gee 6/01              modified by    sh 02/03
This routine performs the integration of the linear element stiffness
matrix ke for a shell9 element
</pre>
\param  ELEMENT   *ele          (i)  the element structure
\param  S9_DATA   *data         (i)  element integration data
\param  MATERIAL  *mat          (i)  the material structure
\param  ARRAY     *estif_global (o)  element stiffness matrix (NOT initialized!)
\param  ARRAY     *emass_global (o)  element mass matrix (NOT initialized!)
\param  INT        kintyp       (i)  kintyp=0: geo_lin; =1: upd_lagr; =2: tot_lagr 
\param  DOUBLE    *force       (i/o) global vector for internal forces (initialized!)
\param  INT        kstep        (i)  actual step in nonlinear analysis
\param  INT        init         (i)  init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9static_keug(ELEMENT   *ele,                        /* the element structure */
                   S9_DATA   *data,                       /* element integration data */
                   MATERIAL  *mat,                        /* the material structure */
                   ARRAY     *estif_global,               /* element stiffness matrix (NOT initialized!) */
                   ARRAY     *emass_global,               /* element mass matrix (NOT initialized!) */
                   INT        kintyp,                     /* typ of kinematic formulation */
                   DOUBLE    *force,                      /* global vector for internal forces (initialized!) */
                   INT        kstep,                      /* actual step in nonlinear analysis */
                   INT        init)                       /* init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase */
/*----------------------------------------------------------------------*/
/* if force==NULL no internal forces are calculated                     */
/*----------------------------------------------------------------------*/
{
INT                 i,j,k,l,kl,ml;                          /* some loopers */
DOUBLE              sum;
INT                 ngauss;

INT                 nir,nis,nit;                            /* num GP in r/s/t direction */
INT                 lr, ls, lt;                             /* loopers over GP */
INT                 iel;                                    /* numnp to this element */
INT                 nd;                                     /* ndofs to this element (=numdf*numnp) */

INT                 istore = 0;                             /* controls storing of new stresses to wa */
INT                 newval = 0;                             /* controls evaluation of new stresses    */

INT                 ip;                                     /* actual integration point */
INT                 actlay;                                 /* actual layer */
INT                 num_mlay;                               /* number of material layers to actual kinematic layer */  
INT                 num_klay;                               /* number of kinematic layers to this element*/
INT                 numdf;                                  /* ndofs per node to this element */
INT                 numdof_shell9; 
DOUBLE             *klayhgt;                                /* hight of kinematic layer in percent of total thicknes of element*/
DOUBLE             *mlayhgt;                                /* hight of material layer in percent of adjacent kinematic layer*/
MULTIMAT           *actmultimat;                            /* material of actual material layer */
INT                 rot_axis;                               /* rotation axis of laminat (1=x-; 2=y-; 3=z-axis) */
DOUBLE              phi;                                    /* angle of rotation about rot_axis */

DOUBLE              e1,e2,e3;                               /*GP-coords*/
DOUBLE              facr,facs,fact;                         /* weights at GP */
DOUBLE              xnu;                                    /* value of shell shifter */
DOUBLE              weight;

DOUBLE              condfac;                                /* sdc conditioning factor */
DOUBLE              h2;                                     /* half nodal height */
DOUBLE              hgt;                                    /* element thickness */
DOUBLE              hte[MAXNOD_SHELL9];                     /* element thickness at nodal points */

DOUBLE              detr;                                   /* jacobian determinant in ref conf. at gp */
DOUBLE              detc;                                   /* jacobian determinant in cur conf. at gp */

DOUBLE              h[3];                                   /* working array */
DOUBLE              da;                                     /* area on mid surface */

DOUBLE              stress[6], stress_r[12];                /* stress and stress resultants */
DOUBLE              strain[6];                              /* strains */
DOUBLE             *intforce;

static ARRAY        C_a;         static DOUBLE **C;         /* material tensor */
static ARRAY        D_a;         static DOUBLE **D;         /* material tensor integrated in thickness direction */

static ARRAY4D      a3r_a;       static DOUBLE ***a3r;      /* a3 in reference config -> for each kinematic layer */
static ARRAY4D      a3c_a;       static DOUBLE ***a3c;      /* a3 in current   config (a3r + disp) */

static ARRAY4D      a3kvpr_a;    static DOUBLE ***a3kvpr;   /* partiel derivatives of normal vector ref.config. */
static ARRAY4D      a3kvpc_a;    static DOUBLE ***a3kvpc;   /* partiel derivatives of normal vector cur.config. */

static ARRAY        xrefe_a;     static DOUBLE **xrefe;     /* coords of midsurface in ref config */
static ARRAY        xcure_a;     static DOUBLE **xcure;     /* coords of midsurface in cur condig */
                                        DOUBLE **a3ref;     /* elements directors (lenght 1) */

static ARRAY        funct_a;     static DOUBLE  *funct;     /* shape functions */
static ARRAY        deriv_a;     static DOUBLE **deriv;     /* derivatives of shape functions */

/* mid surface basis vectors and metric tensors */
static ARRAY4D      akovr_a;     static DOUBLE ***akovr;    /* kovariant basis vectors at Int point ref.config. */
static ARRAY4D      akonr_a;     static DOUBLE ***akonr;    /* kontravar.--------------"----------- ref.config. */
static ARRAY4D      amkovr_a;    static DOUBLE ***amkovr;   /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY4D      amkonr_a;    static DOUBLE ***amkonr;   /* kontravar.--------------"------------ ref.config. */

static ARRAY4D      akovc_a;     static DOUBLE ***akovc;    /* kovariant basis vectors at Int point current.config. */
static ARRAY4D      akonc_a;     static DOUBLE ***akonc;    /* kontravar.--------------"----------- current.config. */
static ARRAY4D      amkovc_a;    static DOUBLE ***amkovc;   /* kovaraiant metric tensor at Int point current.config. */
static ARRAY4D      amkonc_a;    static DOUBLE ***amkonc;   /* kontravar.--------------"------------ current.config. */

/* mid surface basis vectors and metric tensors -> help for s9_tvmr.c */
static ARRAY        akovh_a;     static DOUBLE **akovh;     
static ARRAY        akonh_a;     static DOUBLE **akonh;     
static ARRAY        amkovh_a;    static DOUBLE **amkovh;    
static ARRAY        amkonh_a;    static DOUBLE **amkonh;    

/* shell body basis vectors and metric tensors */
static ARRAY        gkovr_a;     static DOUBLE **gkovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        gkonr_a;     static DOUBLE **gkonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        gmkovr_a;    static DOUBLE **gmkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        gmkonr_a;    static DOUBLE **gmkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        gkovc_a;     static DOUBLE **gkovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        gkonc_a;     static DOUBLE **gkonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        gmkovc_a;    static DOUBLE **gmkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        gmkonc_a;    static DOUBLE **gmkonc;    /* kontravar.--------------"------------ current.config. */

static ARRAY        bop_a;       static DOUBLE **bop;       /* B-Operator for compatible strains */
                                 static DOUBLE **estif;     /* element stiffness matrix ke */
                               /*  static DOUBLE **emass;*/     /* element mass matrix */
static ARRAY        work_a;      static DOUBLE **work;      /* working array to do Bt*D*B */

/* arrays for eas */
DOUBLE              detr0[MAXKLAY_SHELL9];                  /* jacobian determinant in ref conf. at mid point*/
INT                 nhyb;                                   /* scnd dim of P */
DOUBLE              epsh[12];                               /* transformed eas strains  */
DOUBLE              mlhgt_eas[1];                           /* hight of material layer in percent of adjacent kinematic layer*/

static ARRAY        alfa_a[MAXKLAY_SHELL9]; static DOUBLE  *alfa[MAXKLAY_SHELL9];    /* alfa-values for each kinematic layer */
static ARRAY        eashelp_a;              static DOUBLE  *eashelp;                 /* working vector for eas */

DOUBLE            **oldDtildinv;
DOUBLE            **oldL;
DOUBLE             *oldRtild;

static ARRAY        P_a;         static DOUBLE **P;         /* eas matrix M */
static ARRAY        transP_a;    static DOUBLE **transP;    /* eas matrix M */
static ARRAY        T_a;         static DOUBLE **T;         /* transformation matrix for eas */
static ARRAY        workeas_a;   static DOUBLE **workeas;   /* eas working array */
static ARRAY        workeas2_a;  static DOUBLE **workeas2;  /* eas working array */

static ARRAY   LKl_a[MAXKLAY_SHELL9];        static DOUBLE **L_kl[MAXKLAY_SHELL9];        /* eas matrix L                 -> one kinematic Layer*/
static ARRAY   LtKl_a[MAXKLAY_SHELL9];       static DOUBLE **Lt_kl[MAXKLAY_SHELL9];       /* eas matrix L transposed      -> one kinematic Layer*/
static ARRAY   DtildKl_a[MAXKLAY_SHELL9];    static DOUBLE **Dtild_kl[MAXKLAY_SHELL9];    /* eas matrix Dtilde            -> one kinematic Layer*/
static ARRAY   DtildinvKl_a[MAXKLAY_SHELL9]; static DOUBLE **Dtildinv_kl[MAXKLAY_SHELL9]; /* eas matrix Dtilde inv        -> one kinematic Layer*/
static ARRAY   RtildKl_a[MAXKLAY_SHELL9];    static DOUBLE  *Rtild_kl[MAXKLAY_SHELL9];    /* eas part of internal forces  -> one kinematic Layer*/


static ARRAY   akovr0_a[MAXKLAY_SHELL9];   static DOUBLE **akovr0[MAXKLAY_SHELL9];  /* kovariant basis vectors at mid point ref.config. -> vor each kinematic layer */          
static ARRAY   akonr0_a[MAXKLAY_SHELL9];   static DOUBLE **akonr0[MAXKLAY_SHELL9];  /* kontravar.--------------"----------- ref.config. -> vor each kinematic layer */          
static ARRAY   amkovr0_a[MAXKLAY_SHELL9];  static DOUBLE **amkovr0[MAXKLAY_SHELL9]; /* kovaraiant metric tensor at mid point ref.config. -> vor each kinematic layer */          
static ARRAY   amkonr0_a[MAXKLAY_SHELL9];  static DOUBLE **amkonr0[MAXKLAY_SHELL9]; /* kontravar.--------------"------------ ref.config. -> vor each kinematic layer */          

/* arrays for ANS */
INT                 ansq;
const INT           nsansmax=6;
INT                 nsansq;                                  /* number of sampling points for ans */
DOUBLE              xr1[6];                                  /* coordinates of collocation points for ANS */
DOUBLE              xs1[6];
DOUBLE              xr2[6];
DOUBLE              xs2[6];                     
DOUBLE              frq[6];
DOUBLE              fsq[6];

/*a metric for each kinematic layer*/
static ARRAY       funct1q_a[6];  static DOUBLE   *funct1q[6];    /* shape functions at collocation points */
static ARRAY       deriv1q_a[6];  static DOUBLE  **deriv1q[6];    /* derivation of these shape functions */

static ARRAY4D     akovr1q_a[6];  static DOUBLE ***akovr1q[6];
static ARRAY4D     akonr1q_a[6];  static DOUBLE ***akonr1q[6];
static ARRAY4D     amkovr1q_a[6]; static DOUBLE ***amkovr1q[6];
static ARRAY4D     amkonr1q_a[6]; static DOUBLE ***amkonr1q[6];
static ARRAY4D     a3kvpr1q_a[6]; static DOUBLE ***a3kvpr1q[6];

static ARRAY4D     akovc1q_a[6];  static DOUBLE ***akovc1q[6];
static ARRAY4D     akonc1q_a[6];  static DOUBLE ***akonc1q[6];
static ARRAY4D     amkovc1q_a[6]; static DOUBLE ***amkovc1q[6];
static ARRAY4D     amkonc1q_a[6]; static DOUBLE ***amkonc1q[6];
static ARRAY4D     a3kvpc1q_a[6]; static DOUBLE ***a3kvpc1q[6];


static ARRAY       funct2q_a[6];  static DOUBLE   *funct2q[6];    /* shape functions at collocation points */
static ARRAY       deriv2q_a[6];  static DOUBLE  **deriv2q[6];    /* derivation of these shape functions */

static ARRAY4D     akovr2q_a[6];  static DOUBLE ***akovr2q[6];
static ARRAY4D     akonr2q_a[6];  static DOUBLE ***akonr2q[6];
static ARRAY4D     amkovr2q_a[6]; static DOUBLE ***amkovr2q[6];
static ARRAY4D     amkonr2q_a[6]; static DOUBLE ***amkonr2q[6];
static ARRAY4D     a3kvpr2q_a[6]; static DOUBLE ***a3kvpr2q[6];

static ARRAY4D     akovc2q_a[6];  static DOUBLE ***akovc2q[6];
static ARRAY4D     akonc2q_a[6];  static DOUBLE ***akonc2q[6];
static ARRAY4D     amkovc2q_a[6]; static DOUBLE ***amkovc2q[6];
static ARRAY4D     amkonc2q_a[6]; static DOUBLE ***amkonc2q[6];
static ARRAY4D     a3kvpc2q_a[6]; static DOUBLE ***a3kvpc2q[6];


#ifdef DEBUG 
dstrc_enter("s9static_keug");
#endif
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
if (init==1)
{
iel = MAXNOD_SHELL9; /* maximum number of nodes for this type of shell */

xrefe     = amdef("xrefe"  ,&xrefe_a,3,MAXNOD_SHELL9,"DA");       
xcure     = amdef("xcure"  ,&xcure_a,3,MAXNOD_SHELL9,"DA");       
a3r       = am4def("a3r"    ,&a3r_a,3,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         
a3c       = am4def("a3c"    ,&a3c_a,3,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         

a3kvpr    = am4def("a3kvpr" ,&a3kvpr_a,3,2,MAXKLAY_SHELL9,0,"D3");         
a3kvpc    = am4def("a3kvpc" ,&a3kvpc_a,3,2,MAXKLAY_SHELL9,0,"D3");        

funct     = amdef("funct"  ,&funct_a,MAXNOD_SHELL9,1,"DV");       
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL9,"DA");       

akovr     = am4def("akovr"  ,&akovr_a,3,3,MAXKLAY_SHELL9,0,"D3");         
akonr     = am4def("akonr"  ,&akonr_a,3,3,MAXKLAY_SHELL9,0,"D3");          
amkovr    = am4def("amkovr" ,&amkovr_a,3,3,MAXKLAY_SHELL9,0,"D3");         
amkonr    = am4def("amkonr" ,&amkonr_a,3,3,MAXKLAY_SHELL9,0,"D3");         

akovc     = am4def("akovc"  ,&akovc_a,3,3,MAXKLAY_SHELL9,0,"D3");         
akonc     = am4def("akonc"  ,&akonc_a,3,3,MAXKLAY_SHELL9,0,"D3");         
amkovc    = am4def("amkovc" ,&amkovc_a,3,3,MAXKLAY_SHELL9,0,"D3");        
amkonc    = am4def("amkonc" ,&amkonc_a,3,3,MAXKLAY_SHELL9,0,"D3");        

akovh     = amdef("akovh"  ,&akovh_a,3,3,"DA");         
akonh     = amdef("akonh"  ,&akonh_a,3,3,"DA");         
amkovh    = amdef("amkovh" ,&amkovh_a,3,3,"DA");        
amkonh    = amdef("amkonh" ,&amkonh_a,3,3,"DA");        

gkovr     = amdef("gkovr"  ,&gkovr_a,3,3,"DA");         
gkonr     = amdef("gkonr"  ,&gkonr_a,3,3,"DA");         
gmkovr    = amdef("gmkovr" ,&gmkovr_a,3,3,"DA");        
gmkonr    = amdef("gmkonr" ,&gmkonr_a,3,3,"DA");        

gkovc     = amdef("gkovc"  ,&gkovc_a,3,3,"DA");         
gkonc     = amdef("gkonc"  ,&gkonc_a,3,3,"DA");         
gmkovc    = amdef("gmkovc" ,&gmkovc_a,3,3,"DA");        
gmkonc    = amdef("gmkonc" ,&gmkonc_a,3,3,"DA");        

bop       = amdef("bop"    ,&bop_a ,12,(NUMDOF_SHELL9*MAXNOD_SHELL9),"DA");
C         = amdef("C"      ,&C_a   ,6 ,6                    ,"DA");             
D         = amdef("D"      ,&D_a   ,12,12                   ,"DA");           
work      = amdef("work"   ,&work_a,12,(MAXNOD_SHELL9*NUMDOF_SHELL9),"DA"); 

/* for eas */
P         = amdef("P"      ,&P_a       ,12           ,MAXHYB_SHELL9                ,"DA");         
transP    = amdef("transP" ,&transP_a  ,12           ,MAXHYB_SHELL9                ,"DA"); 
T         = amdef("T"      ,&T_a       ,12           ,12                           ,"DA");
workeas   = amdef("workeas",&workeas_a ,12           ,(MAXNOD_SHELL9*NUMDOF_SHELL9),"DA");
/*workeas   = amdef("workeas", &workeas_a ,MAXHYB_SHELL9        ,(MAXNOD_SHELL9*NUMDOF_SHELL9),"DA");*/
workeas2  = amdef("workeas2",&workeas2_a,(MAXNOD_SHELL9*NUMDOF_SHELL9),MAXHYB_SHELL9        ,"DA");

for (i=0; i<MAXKLAY_SHELL9; i++) alfa[i]       = amdef("alfa"      ,&alfa_a[i]       ,MAXHYB_SHELL9     ,1                ,"DV"); 
for (i=0; i<MAXKLAY_SHELL9; i++) L_kl[i]       = amdef("L_kl"      ,&LKl_a[i]        ,(MAXNOD_SHELL9*NUMDOF_SHELL9),MAXHYB_SHELL9     ,"DA");
for (i=0; i<MAXKLAY_SHELL9; i++) Lt_kl[i]      = amdef("Lt_kl"     ,&LtKl_a[i]       ,MAXHYB_SHELL9     ,(MAXNOD_SHELL9*NUMDOF_SHELL9),"DA");
for (i=0; i<MAXKLAY_SHELL9; i++) Dtild_kl[i]   = amdef("Dtild_kl"  ,&DtildKl_a[i]    ,MAXHYB_SHELL9     ,MAXHYB_SHELL9    ,"DA");
for (i=0; i<MAXKLAY_SHELL9; i++) Dtildinv_kl[i]= amdef("Dtildi_kl" ,&DtildinvKl_a[i] ,MAXHYB_SHELL9     ,MAXHYB_SHELL9    ,"DA");
for (i=0; i<MAXKLAY_SHELL9; i++) Rtild_kl[i]   = amdef("Rtild_kl"  ,&RtildKl_a[i]    ,MAXHYB_SHELL9     ,1                ,"DV");

eashelp   = amdef("eashelp",&eashelp_a ,MAXHYB_SHELL9        ,1                    ,"DV");

for (i=0; i<MAXKLAY_SHELL9; i++) akovr0[i]  = amdef("akovr0"  ,&akovr0_a[i]  ,3,3,"DA");  
for (i=0; i<MAXKLAY_SHELL9; i++) akonr0[i]  = amdef("akonr0"  ,&akonr0_a[i]  ,3,3,"DA");  
for (i=0; i<MAXKLAY_SHELL9; i++) amkovr0[i] = amdef("amkovr0" ,&amkovr0_a[i] ,3,3,"DA");  
for (i=0; i<MAXKLAY_SHELL9; i++) amkonr0[i] = amdef("amkonr0" ,&amkonr0_a[i] ,3,3,"DA");  

/* for ans */
for (i=0; i<nsansmax; i++) funct1q[i]  = amdef("funct1q",&(funct1q_a[i]),MAXNOD_SHELL9,1,"DV");
for (i=0; i<nsansmax; i++) deriv1q[i]  = amdef("deriv1q",&(deriv1q_a[i]),2,MAXNOD_SHELL9,"DA");

for (i=0; i<nsansmax; i++) akovr1q[i]  = am4def("akovr1q" ,&(akovr1q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akonr1q[i]  = am4def("akonr1q" ,&(akonr1q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkovr1q[i] = am4def("amkovr1q",&(amkovr1q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkonr1q[i] = am4def("amkonr1q",&(amkonr1q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) a3kvpr1q[i] = am4def("a3kvpr1q",&(a3kvpr1q_a[i]),3,2,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akovc1q[i]  = am4def("akovc1q" ,&(akovc1q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akonc1q[i]  = am4def("akonc1q" ,&(akonc1q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkovc1q[i] = am4def("amkovc1q",&(amkovc1q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkonc1q[i] = am4def("amkonc1q",&(amkonc1q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) a3kvpc1q[i] = am4def("a3kvpc1q",&(a3kvpc1q_a[i]),3,2,MAXKLAY_SHELL9,0,"D3");

for (i=0; i<nsansmax; i++) funct2q[i]  = amdef("funct2q",&(funct2q_a[i]),MAXNOD_SHELL9,1,"DV");
for (i=0; i<nsansmax; i++) deriv2q[i]  = amdef("deriv2q",&(deriv2q_a[i]),2,MAXNOD_SHELL9,"DA");

for (i=0; i<nsansmax; i++) akovr2q[i]  = am4def("akovr2q" ,&(akovr2q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akonr2q[i]  = am4def("akonr2q" ,&(akonr2q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkovr2q[i] = am4def("amkovr2q",&(amkovr2q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkonr2q[i] = am4def("amkonr2q",&(amkonr2q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) a3kvpr2q[i] = am4def("a3kvpr2q",&(a3kvpr2q_a[i]),3,2,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akovc2q[i]  = am4def("akovc2q" ,&(akovc2q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) akonc2q[i]  = am4def("akonc2q" ,&(akonc2q_a[i]) ,3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkovc2q[i] = am4def("amkovc2q",&(amkovc2q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) amkonc2q[i] = am4def("amkonc2q",&(amkonc2q_a[i]),3,3,MAXKLAY_SHELL9,0,"D3");
for (i=0; i<nsansmax; i++) a3kvpc2q[i] = am4def("a3kvpc2q",&(a3kvpc2q_a[i]),3,2,MAXKLAY_SHELL9,0,"D3");

goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase                                                         */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
amdel(&xrefe_a);
amdel(&xcure_a);
am4del(&a3r_a);
am4del(&a3c_a);

am4del(&a3kvpr_a);
am4del(&a3kvpc_a);

amdel(&funct_a);
amdel(&deriv_a);

am4del(&akovr_a);   
am4del(&akonr_a);   
am4del(&amkovr_a);  
am4del(&amkonr_a);  

am4del(&akovc_a);   
am4del(&akonc_a);   
am4del(&amkovc_a);  
am4del(&amkonc_a);

amdel(&akovh_a);   
amdel(&akonh_a);   
amdel(&amkovh_a);  
amdel(&amkonh_a);  

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

/* for eas */
amdel(&P_a);
amdel(&transP_a);
amdel(&T_a);
amdel(&workeas_a);
amdel(&workeas2_a);

for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&akovr0_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&akonr0_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&amkovr0_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&amkonr0_a[i]);

for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&LKl_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&LtKl_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&DtildKl_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&DtildinvKl_a[i]);
for (i=0; i<MAXKLAY_SHELL9; i++) amdel(&RtildKl_a[i]);

amdel(&eashelp_a);

/* for ans */
for (i=0; i<nsansmax; i++) amdel( &(funct1q_a[i]) );
for (i=0; i<nsansmax; i++) amdel( &(deriv1q_a[i]) );

for (i=0; i<nsansmax; i++) am4del( &(akovr1q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(akonr1q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(amkovr1q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(amkonr1q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(a3kvpr1q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(akovc1q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(akonc1q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(amkovc1q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(amkonc1q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(a3kvpc1q_a[i]));

for (i=0; i<nsansmax; i++) amdel( &(funct2q_a[i]) );
for (i=0; i<nsansmax; i++) amdel( &(deriv2q_a[i]) );

for (i=0; i<nsansmax; i++) am4del( &(akovr2q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(akonr2q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(amkovr2q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(amkonr2q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(a3kvpr2q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(akovc2q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(akonc2q_a[i]) );
for (i=0; i<nsansmax; i++) am4del( &(amkovc2q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(amkonc2q_a[i]));
for (i=0; i<nsansmax; i++) am4del( &(a3kvpc2q_a[i]));

goto end;  
}
/*----------------------------------------------------------------------*/
/* update phase  for material nonlinearity      (init=2)                */
/*----------------------------------------------------------------------*/
else if(init==2)
{
  istore = 1;             /*-- material law is called with this flag----*/
}

/*----------------------------------------------------------------------*/
/* calculation phase                                                    */
/*----------------------------------------------------------------------*/
num_klay = ele->e.s9->num_klay;         /* number of kinematic layers to this element*/
numdf = ele->e.s9->numdf;               /* ndofs per node to this element */
numdof_shell9 = numdf;

/*-------------------------------------------- init the gaussian points */
s9intg(ele,data,0);
/*----------------- some of the fields have to be reinitialized to zero */
/*amzero(&D_a);*/
amzero(estif_global);
estif     = estif_global->a.da;
/*----------------------------------------------- integrationsparameter */
nir     = ele->e.s9->nGP[0];
nis     = ele->e.s9->nGP[1];
nit     = ele->e.s9->nGP[2];
iel     = ele->numnp;
nd      = iel*numdf; 
condfac = ele->e.s9->sdc;
a3ref   = ele->e.s9->a3ref.a.da;
          amzero(&(ele->e.s9->intforce));
intforce= ele->e.s9->intforce.a.dv;
/*------------------------------------------------------- check for eas */
nhyb=ele->e.s9->nhyb;
if (nhyb>0) 
{   
   for (kl=0; kl<num_klay; kl++) 
   { 
     amzero(&LKl_a[kl]);       /* array, that holds the information for ONE kinematic layers*/
     amzero(&LtKl_a[kl]);      /* array, that holds the information for ONE kinematic layers*/
     amzero(&DtildKl_a[kl]);   /* array, that holds the information for ONE kinematic layers*/
     amzero(&RtildKl_a[kl]);   /* array, that holds the information for ONE kinematic layers*/
   }

   mlhgt_eas[0]= 100.0;

   /*----- check if there is a new loadstep -> alfa-values have to be reinitialized to ZERO*/
   if (kstep == 0) ele->e.s9->oldkstep = 0;     /* initialize oldkstep in first kstep*/
   if (kstep != ele->e.s9->oldkstep) 
   {
      for (kl=0; kl<num_klay; kl++) for (i=0; i<nhyb; i++) ele->e.s9->alfa.a.da[kl][i] = 0.0;
      ele->e.s9->oldkstep = kstep;
   }
   
   for (kl=0; kl<num_klay; kl++) 
   { 
      /*--------------------------------------- update of eas strains alfa */
      /*---------------------------- set pointer to actual kinematic layer */
      alfa[kl] = ele->e.s9->alfa.a.da[kl];
      /*-------------------- set pointer to storage of old Dtildinv and Lt */
      oldDtildinv = &ele->e.s9->Dtildinv.a.da[kl*nhyb];
      oldL        = &ele->e.s9->L.a.da[kl*nd];
      oldRtild    = &ele->e.s9->Rtilde.a.dv[kl*nhyb];

      /*------ set number of "displacement - dofs" to this kinematic layer */
      iel = ele->numnp;
      /*------------------- make multiplication eashelp = Lt * disp[kstep] */
      for (i=0; i<nhyb; i++)
      {
         sum=0.0;
         for (j=0; j<iel; j++)
         for (k=0; k<numdf; k++)
         {
            l = j*numdf+k;
            sum += oldL[l][i]*ele->node[j]->sol_residual.a.da[0][k];
         }
         eashelp[i] = sum;
      }

      /*---------------------------------------- add old Rtilde to eashelp */
      for (i=0; i<nhyb; i++) eashelp[i] += oldRtild[i];
      /*----------------- make multiplication alfa -= olDtildinv * eashelp */
      math_matvecdense(alfa[kl],oldDtildinv,eashelp,nhyb,nhyb,1,-1.0);
   }

}
/*----------------------------------------------------- geometry update */
for (kl=0; kl<num_klay; kl++) /*loop over all kinematic layers*/
{  
  klayhgt = ele->e.s9->klayhgt;   /* hgt of kinematic layer on percent of total thickness of shell */
  for (k=0; k<iel; k++)           /*loop over all nodes per layer*/
  {
     hte[k] = ele->e.s9->thick_node.a.dv[k];
     /*if (ele->e.s9->dfield == 0)*/      /*Layerthicknes, norm(a3L) = HL */
     h2 = ele->e.s9->thick_node.a.dv[k] * klayhgt[kl]/100. * condfac;
     /*h2 = 0.5*h2;*/ /*A3_IST_EINHALB halber Direktor*/
     h2 = A3FAC_SHELL9 * h2;
     /*else if (ele->e.s9->dfield == 1)*/ /*half of shell thickness, norm(a3) = H/2*/
     /*  h2 = ele->e.s9->thick_node.a.dv[k]/2. * condfac;*/
 
     a3r[0][k][kl] = a3ref[0][k] * h2;
     a3r[1][k][kl] = a3ref[1][k] * h2;
     a3r[2][k][kl] = a3ref[2][k] * h2;

     xrefe[0][k] = ele->node[k]->x[0];
     xrefe[1][k] = ele->node[k]->x[1];
     xrefe[2][k] = ele->node[k]->x[2];

     xcure[0][k] = xrefe[0][k] + ele->node[k]->sol.a.da[0][0];
     xcure[1][k] = xrefe[1][k] + ele->node[k]->sol.a.da[0][1];
     xcure[2][k] = xrefe[2][k] + ele->node[k]->sol.a.da[0][2];
 
     a3c[0][k][kl] = a3r[0][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+3];
     a3c[1][k][kl] = a3r[1][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+4];
     a3c[2][k][kl] = a3r[2][k][kl]  + ele->node[k]->sol.a.da[0][3*kl+5];
  } /*end loop over nodes*/
} /*end loop over kinematic layers*/
/*============ metric and shape functions at collocation points (ANS=1) */
/*----------------------------------------------------- 4-noded element */
if (ele->e.s9->ans==1)/*------------ querschub_ans */
{
                 ansq=1;
   if (iel==4) nsansq=2;
   if (iel==9) nsansq=6;
   s9_ans_colloqpoints(nsansq,iel,ele->e.s9->ans,ele->distyp,
                       xr1,xs1,xr2,xs2,
                       funct1q,deriv1q,funct2q,deriv2q,
                       xrefe,a3r,xcure,a3c,
                       akovr1q,akonr1q,amkovr1q,amkonr1q,a3kvpr1q,
                       akovc1q,akonc1q,amkovc1q,amkonc1q,a3kvpc1q,
                       akovr2q,akonr2q,amkovr2q,amkonr2q,a3kvpr2q,
                       akovc2q,akonc2q,amkovc2q,amkonc2q,a3kvpc2q,
                       akovh,akonh,amkovh,amkonh,
                       num_klay);
}
else
{
   ansq=0;
   nsansq=0;
}
/*============ metric of element mid point (for eas) -> for each kinematic layer */
if (nhyb>0)
{
   s9_funct_deriv(funct,deriv,0.0,0.0,ele->distyp,1);

   s9_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,akovh,akonh,amkovh,amkonh,
           &detr,funct,deriv,iel,a3kvpr,num_klay);
   /*make hgt at mid point of the element*/
   s9_xint(&hgt,hte,funct,iel);
   
   for (kl=0; kl<num_klay; kl++)
   {
      s9_tmtr(0.0,akovr0[kl],akonr0[kl],amkovr0[kl],amkonr0[kl],&detr0[kl],
                 akovr,a3kvpr,hgt,klayhgt,mlhgt_eas,num_klay,kl,0,condfac);
   }        
}
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/******** equivalent element length -> CARAT! -> s9tvci *******************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*=================================================== integration loops */
ngauss=0;
for (lr=0; lr<nir; lr++)   /* loop in r-direction */
{
   /*================================== gaussian point and weight at it */
   e1   = data->xgpr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)  /* loop in s-direction*/
   {
      /*=============================== gaussian point and weight at it */
      e2   = data->xgps[ls];
      facs = data->wgts[ls];
      /*-------------------- shape functions at gp e1,e2 on mid surface */
      s9_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*----------------------------- shape functions for querschub-ans */
      if (ansq==1) s9_ansq_funct(frq,fsq,e1,e2,iel);
/************/
      /*-------- init mid surface material tensor and stress resultants */
/*      amzero(&D_a);                          */
/*      for (i=0; i<12; i++) stress_r[i]=0.0;  */
/*************/
      /*----- make height at gaussian point -> variable thickness within the element is normally not allowed */
      s9_xint(&hgt,hte,funct,iel);
      /*-- metrics at gaussian point in reference layers of kin. layers */
      s9_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,akovh,akonh,amkovh,amkonh,
              &detr,funct,deriv,iel,a3kvpr,num_klay);
      s9_tvmr(xcure,a3c,akovc,akonc,amkovc,amkonc,akovh,akonh,amkovh,amkonh,
              &detc,funct,deriv,iel,a3kvpc,num_klay);
      /*------------------------- make h as cross product in ref config.
                                       to get area da on shell mid surf */
      h[0] = akovr[1][0][0]*akovr[2][1][0] - akovr[2][0][0]*akovr[1][1][0];
      h[1] = akovr[2][0][0]*akovr[0][1][0] - akovr[0][0][0]*akovr[2][1][0];
      h[2] = akovr[0][0][0]*akovr[1][1][0] - akovr[1][0][0]*akovr[0][1][0];
      /*------------------------------------- make director unit length 
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);
/*---------------------------loop over all kinematic layers (num_klay) */

      actlay = 0;
      for (kl=0; kl<num_klay; kl++)   /*loop over all kinematic layers */
      {
         /*-------- init mid surface material tensor and stress resultants */
         amzero(&D_a);
         for (i=0; i<12; i++) stress_r[i]=0.0;
         /*----------------------------------------------------------------*/
         num_mlay = ele->e.s9->kinlay[kl].num_mlay;
         mlayhgt  = ele->e.s9->kinlay[kl].mlayhgt;   /* hgt of material layer in percent of this kin layer */
         /*--------------------------------------- make eas if switched on */
         if (nhyb>0)
         {
             /*initialize some arrays to zero*/
             amzero(&P_a);
             amzero(&T_a);
             amzero(&transP_a);

             /*------------------- make shape functions for incomp. strains */
             s9_eas(nhyb,e1,e2,iel,ele->e.s9->eas,P);

             /*-------------------- transform basis of Eij to Gausian point */
             s9_transeas(P,transP,T,akovr,akonr0[kl],detr,detr0[kl],nhyb,kl);
             /*------------------------ transform strains to Gaussian point */
             math_matvecdense(epsh,transP,alfa[kl],12,nhyb,0,1.0);
         }
         /*------------------------ make B-operator for compatible strains */
         amzero(&bop_a);  /*initialize bop to ZERO */
         if (kintyp > 0) /*geometric nonlinear -> current configuration */
         {
            s9_tvbo(bop,funct,deriv,iel,numdf,akovc,a3kvpc,num_klay,kl,condfac,nsansq);
         }
         else  /*geometric linear -> reference configuration */
         {
            s9_tvbo(bop,funct,deriv,iel,numdf,akovr,a3kvpr,num_klay,kl,condfac,nsansq);
         }
         /*-------------------------------------- modifications due to ans */
         if (ansq==1) /*Querschub*/
         s9_ans_bbar_q(bop,frq,fsq,funct1q,funct2q,deriv1q,deriv2q,
                       akovc1q,akovc2q,a3kvpc1q,a3kvpc2q,iel,numdf,
                       num_klay,kl,condfac,nsansq);

         for(ml=0; ml<num_mlay; ml++) /*loop over all material layers of aktuel kinematic layer*/
         {
            for (lt=0; lt<nit; lt++)  /*loop in t-direction */
            {
               /*---------------------------- gaussian point and weight at it */
               e3   = data->xgpt[lt];
               fact = data->wgtt[lt];
               /*-------------------- basis vectors and metrics at shell body */ 
               s9_tmtr(e3,gkovr,gkonr,gmkovr,gmkonr,&detr,akovr,a3kvpr,hgt,
                       klayhgt,mlayhgt,num_klay,kl,ml,condfac);

               s9_tmtr(e3,gkovc,gkonc,gmkovc,gmkonc,&detc,akovc,a3kvpc,hgt,
                       klayhgt,mlayhgt,num_klay,kl,ml,condfac);

               /*---------- metric at gp in shell body -> for geo_lin/geo_nl --*/
               if (ansq==0)
               s9_tvhe(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,&detr,&detc,
                       amkovc,amkovr,akovc,akovr,a3kvpc,a3kvpr,e3,kintyp,hgt,
                       klayhgt,mlayhgt,num_klay,kl,ml,condfac);
               /*- modifications to metric of shell body due to querschub-ans */
               else
               s9_ans_tvhe_q(gmkovr,gmkovc,gmkonr,gmkonc,amkovc,amkovr,
                             akovc,akovr,a3kvpc,a3kvpr,&detr,&detc,
                             amkovr1q,amkovc1q,amkovr2q,amkovc2q,
                             frq,fsq,e3,nsansq,hgt,klayhgt,mlayhgt,
                             num_klay,kl,ml,condfac);
               /*----------- calc shell shifter and put it in the weight fact */
               /* xnu = (0.5/condfac)*(detr/da); */               
               /* xnu = (1.0/condfac)*(detr/da); */ /*A3_IST_EINHALB*/
               xnu = (0.5/A3FAC_SHELL9) * (1.0/condfac)*(detr/da);
               fact *= xnu;
               /*----------------------- change to current metrics due to eas */
               if (nhyb>0) s9_vthv(gmkovc,gmkonc,epsh,&detc,e3,hgt,klayhgt,mlayhgt,
                                   num_klay,kl,ml,condfac);
               /*------------------------------------------ call material law */
               actmultimat = &(multimat[ele->e.s9->kinlay[kl].mmatID[ml]-1]);
               rot_axis = mat->m.multi_layer->kinlay[kl].rot[ml];
               phi = mat->m.multi_layer->kinlay[kl].phi[ml];

/*               ip = 2 * ngauss + lt;   if 2 GP per thickness !!*/
               ip = nit * ngauss + lt;
               s9_call_mat(ele,actmultimat,stress,strain,C,gmkovc,gmkovr,gmkonr,
                           gkovr,gkonr,rot_axis,phi,ip,actlay,istore,newval);
               /*---------------- do thickness integration of material tensor */
               s9_tvma(D,C,stress,stress_r,e3,fact,hgt,klayhgt,mlayhgt,
                       num_klay,kl,ml,condfac);
            }/*========================================== end of loop over lt */            

         actlay ++;
         }/*======= end of loop over all material layers of aktual kinematic layer*/

         /*------------ product of all weights and jacobian of mid surface */            
         weight = facr*facs*da;
         /*----------------------------------- elastic stiffness matrix ke */
         s9_BtDB(estif,bop,D,iel,numdf,weight,work);
         
         if (kintyp > 0) /*geometric nonlinear*/
         {
           /*--------------------------------- geometric stiffness matrix kg */
           if (ansq==0)
           s9_tvkg(estif,stress_r,funct,deriv,numdf,iel,weight,kl,num_klay);
           else
           s9_ans_tvkg(estif,stress_r,funct,deriv,numdf,iel,weight,
                       frq,fsq,funct1q,funct2q,deriv1q,deriv2q,
                       ansq,nsansq,kl,num_klay);
         }
         /*-------------------------------- calculation of internal forces */
         if (force) s9_intforce(intforce,stress_r,bop,iel,numdf,12,weight);

         /*----------------------------------- integration of eas matrices */
         if (nhyb>0)
         {
            /* one set of matrizes for each kinematic layer !!! */
            /*=============================================================*/
            /* L(nd,nhyp) = Btrans(nd,12) * D(12,12) * M(12,nhyb)          */
            /*=============================================================*/
            /*----------------------------------------------------DM = D*M */
            math_matmatdense(workeas,D,transP,12,12,nhyb,0,0.0);
            /*--------------------------------------------- L = Bt * D * M */
            math_mattrnmatdense(L_kl[kl],bop,workeas,nd,12,nhyb,1,weight);
            /*=============================================================*/
            /* Ltrans(nhyb,nd) = Mtrans(nhyb,12) * D(12,12) * B_R(12,nd_kl)*/
            /*=============================================================*/
            /*----------------------------------------------------- DB=D*B */
            math_matmatdense(workeas,D,bop,12,12,nd,0,0.0);
            /*----------------------------------- Ltransposed = Mt * D * B */
            math_mattrnmatdense(Lt_kl[kl],transP,workeas,nhyb,12,nd,1,weight);
            /*=============================================================*/
            /* Dtilde(nhyb,nhyb) = Mtrans(nhyb,12) * D(12,12) * M(12,nhyb) */
            /*=============================================================*/
            /*----------------------------------------------------DM = D*M */
            math_matmatdense(workeas,D,transP,12,12,nhyb,0,0.0);
            /*-------------------------------------------- Dtilde = Mt*D*M */
            math_mattrnmatdense(Dtild_kl[kl],transP,workeas,nhyb,12,nhyb,1,weight);
            /*=============================================================*/
            /* Rtilde(nhyb) = Mtrans(nhyb,12) * Forces(12)                 */
            /*=============================================================*/
            /*------------------------- eas part of internal forces Rtilde */
            math_mattrnvecdense(Rtild_kl[kl],transP,stress_r,nhyb,12,1,weight);
          }
      }/*============================== end loop over all kinematic layers */
      ngauss++;
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */


/*----------------- make modifications to stiffness matrices due to eas */
if (nhyb>0)
{
   /*store Dtildinv_kl,Lt_kl,Rtild_kl for each kinematic layer and modify estif*/
   for (kl=0; kl<num_klay; kl++)
   {
     /*------------------------------------ make inverse of matrix Dtilde */
     amcopy(&DtildKl_a[kl],&DtildinvKl_a[kl]);

/**** Test the symmetry of Dtildinv ***/
/*for (i=0; i<nhyb; i++)
for (j=i+1; j<nhyb; j++)
if (FABS(Dtildinv_kl[kl][i][j]-Dtildinv_kl[kl][j][i])>EPS9) 
printf(" Dtild[%d][%d] is not sym with %E\n",   i,j,Dtildinv_kl[kl][i][j]-Dtildinv_kl[kl][j][i]);*/
/**** Test the symmetry of Dtildinv ***/

     math_unsym_inv(Dtildinv_kl[kl],nhyb,nhyb);  /*for unsymmetric D-Matrixes*/
     /*math_sym_inv(Dtildinv_kl[kl],nhyb);*/

     /*------------------------------------------ put Dtildinv_kl to storage */
     for (i=0; i<nhyb; i++)
     for (j=0; j<nhyb; j++) ele->e.s9->Dtildinv.a.da[kl*nhyb + i][j] = Dtildinv_kl[kl][i][j];
     /*------------------------------------------------ put L_kl to storage */
     for (i=0; i<nd; i++)
     for (j=0; j<nhyb; j++) ele->e.s9->L.a.da[kl*nd + i][j] = L_kl[kl][i][j];
     /*-------------------------------------------- put Rtilde_kl to storage */
     for (i=0; i<nhyb; i++) ele->e.s9->Rtilde.a.dv[kl*nhyb + i] = Rtild_kl[kl][i];

     /*===================================================================*/
     /* estif(nd,nd) = estif(nd,nd) - Ltrans(nhyb,nd) * Dtilde^-1(nhyb,nhyb) * L(nd,nhyb) */
     /*===================================================================*/

     /*------------------------------------------- make Ltrans * Dtildinv */
     math_matmatdense(workeas2,L_kl[kl],Dtildinv_kl[kl],nd,nhyb,nhyb,0,0.0);
     /*---------------------------------- make estif -= Lt * Dtildinv * L */
     math_matmatdense(estif,workeas2,Lt_kl[kl],nd,nhyb,nd,1,-1.0);


     /*===================================================================*/
     /* R(12) = R(12) - Ltrans(nhyb,nd) * Dtilde^-1(nhyb,nhyb) * Rtilde(nhyb) */
     /*===================================================================*/
     /*--------------------------- make intforce -= Lt * Dtildinv * Rtild */
     math_matvecdense(intforce,workeas2,Rtild_kl[kl],nd,nhyb,1,-1.0);
   }
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
/**** Test the symmetry of estif ***/
/*for (i=0; i<nd; i++)
for (j=i+1; j<nd; j++)
if (FABS(estif[i][j]-estif[j][i])>EPS12) printf("i %d j %d not sym with %E\n",
                                         i,j,estif[i][j]-estif[j][i]);*/
/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9static_keug */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
