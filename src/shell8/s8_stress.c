#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |routine to calc element stresses                       m.gee 12/01    |
 *----------------------------------------------------------------------*/
void s8_stress(ELEMENT      *ele,
               S8_DATA      *data,
               MATERIAL     *mat,
               int           kstep,
               int           init)
{
int                 i,j,k,l;
int                 nir,nis,nit;
int                 ngauss;

int                 lr,ls,lt;
double              e1,e2,e3;
double              facr,facs,fact;

int                 iel;
int                 nd;

double              condfac;
double              h2;
double            **a3ref;

double              detsmr;
double              detsmc;
double              detsrr;
double              detsrc;
 

double              hte[MAXNOD_SHELL8];
double              hhi;
double              h[3];                                   /* working array */
double              da;                                     /* area on mid surface */

double              stress[6], stress_r[12];                /* stress and stress resultants */
double              strain[6];                              /* strains */

static ARRAY        stress_a;    static double **gp_stress; /* element array for stresses on gaussian points */
                                      double ***ele_stress; /* pointer to array of stress history in element */

static ARRAY        C_a;         static double **C;         /* material tensor */
static ARRAY        D_a;         static double **D;         /* material tensor integrated in thickness direction */

static ARRAY        a3r_a;       static double **a3r;       /* a3 in reference config (lenght h/2) */
static ARRAY        a3c_a;       static double **a3c;       /* a3 in current   config (lenght h/2 + disp) */

static ARRAY        xrefe_a;     static double **xrefe;     /* coords of midsurface in ref config */
static ARRAY        xcure_a;     static double **xcure;     /* coords of midsurface in cur condig */

static ARRAY        funct_a;     static double  *funct;     /* shape functions */
static ARRAY        deriv_a;     static double **deriv;     /* derivatives of shape functions */

static ARRAY        akovr_a;     static double **akovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        akonr_a;     static double **akonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        amkovr_a;    static double **amkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        amkonr_a;    static double **amkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        akovc_a;     static double **akovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        akonc_a;     static double **akonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        amkovc_a;    static double **amkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        amkonc_a;    static double **amkonc;    /* kontravar.--------------"------------ current.config. */

static ARRAY        a3kvpr_a;    static double **a3kvpr;    /* partiel derivatives of normal vector ref.config. */
static ARRAY        a3kvpc_a;    static double **a3kvpc;    /* partiel derivatives of normal vector cur.config. */

/* shell body basis vectors and metric tensors */
static ARRAY        gkovr_a;     static double **gkovr;     /* kovariant basis vectors at Int point ref.config. */
static ARRAY        gkonr_a;     static double **gkonr;     /* kontravar.--------------"----------- ref.config. */
static ARRAY        gmkovr_a;    static double **gmkovr;    /* kovaraiant metric tensor at Int point ref.config. */
static ARRAY        gmkonr_a;    static double **gmkonr;    /* kontravar.--------------"------------ ref.config. */

static ARRAY        gkovc_a;     static double **gkovc;     /* kovariant basis vectors at Int point current.config. */
static ARRAY        gkonc_a;     static double **gkonc;     /* kontravar.--------------"----------- current.config. */
static ARRAY        gmkovc_a;    static double **gmkovc;    /* kovaraiant metric tensor at Int point current.config. */
static ARRAY        gmkonc_a;    static double **gmkonc;    /* kontravar.--------------"------------ current.config. */

#ifdef DEBUG 
dstrc_enter("s8_stress");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
if (init==1)
{
   iel       = MAXNOD_SHELL8; /* maximum number of nodes for this type of shell */
   xrefe     = amdef("xrefe"  ,&xrefe_a,3,MAXNOD_SHELL8,"DA");       
   xcure     = amdef("xcure"  ,&xcure_a,3,MAXNOD_SHELL8,"DA");       
   a3r       = amdef("a3r"    ,&a3r_a  ,3,MAXNOD_SHELL8,"DA");         
   a3c       = amdef("a3c"    ,&a3c_a  ,3,MAXNOD_SHELL8,"DA");         
   funct     = amdef("funct"  ,&funct_a,MAXNOD_SHELL8,1,"DV");       
   deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL8,"DA");       

   akovr     = amdef("akovr"  ,&akovr_a,3,3,"DA");         
   akonr     = amdef("akonr"  ,&akonr_a,3,3,"DA");         
   amkovr    = amdef("amkovr" ,&amkovr_a,3,3,"DA");        
   amkonr    = amdef("amkonr" ,&amkonr_a,3,3,"DA");        

   akovc     = amdef("akovc"  ,&akovc_a,3,3,"DA");         
   akonc     = amdef("akonc"  ,&akonc_a,3,3,"DA");         
   amkovc    = amdef("amkovc" ,&amkovc_a,3,3,"DA");        
   amkonc    = amdef("amkonc" ,&amkonc_a,3,3,"DA");        

   a3kvpr    = amdef("a3kvpr" ,&a3kvpr_a,3,2,"DA");        
   a3kvpc    = amdef("a3kvpc" ,&a3kvpc_a,3,2,"DA");        

   gkovr     = amdef("gkovr"  ,&gkovr_a,3,3,"DA");
   gkonr     = amdef("gkonr"  ,&gkonr_a,3,3,"DA");
   gmkovr    = amdef("gmkovr" ,&gmkovr_a,3,3,"DA");
   gmkonr    = amdef("gmkonr" ,&gmkonr_a,3,3,"DA");
 
   gkovc     = amdef("gkovc"  ,&gkovc_a,3,3,"DA");
   gkonc     = amdef("gkonc"  ,&gkonc_a,3,3,"DA");
   gmkovc    = amdef("gmkovc" ,&gmkovc_a,3,3,"DA");
   gmkonc    = amdef("gmkonc" ,&gmkonc_a,3,3,"DA");
  
   C         = amdef("C"      ,&C_a   ,6 ,6                    ,"DA");             
   D         = amdef("D"      ,&D_a   ,12,12                   ,"DA");   

   gp_stress = amdef("gp_stress",&stress_a,18,MAXGAUSS,"DA");        

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
   amdel(&funct_a);
   amdel(&deriv_a);

   amdel(&akovr_a);   
   amdel(&akonr_a);   
   amdel(&amkovr_a);  
   amdel(&amkonr_a);  

   amdel(&akovc_a);   
   amdel(&akonc_a);   
   amdel(&amkovc_a);  
   amdel(&amkonc_a);

   amdel(&a3kvpr_a);
   amdel(&a3kvpc_a);

   amdel(&gkovr_a);
   amdel(&gkonr_a);
   amdel(&gmkovr_a);
   amdel(&gmkonr_a);

   amdel(&gkovc_a);
   amdel(&gkonc_a);
   amdel(&gmkovc_a);  
   amdel(&gmkonc_a);

   amdel(&C_a);
   amdel(&D_a);

   amdel(&stress_a);
   goto end;
}
/*----------------------------------------------------------------------*/
/* calculation phase                                                    */
/*----------------------------------------------------------------------*/
/*-------------------------------------------- init the gaussian points */
s8intg(ele,data,0);
amzero(&stress_a);
/*----------------------------------------------- integrationsparameter */
nir     = ele->e.s8->nGP[0];
nis     = ele->e.s8->nGP[1];
nit     = ele->e.s8->nGP[2];
iel     = ele->numnp;
nd      = iel*NUMDOF_SHELL8;
condfac = ele->e.s8->sdc;
a3ref   = ele->e.s8->a3ref.a.da;
/*----------------------------------------------------- geometry update */
for (k=0; k<iel; k++)
{
   hte[k] = ele->e.s8->thick_node.a.dv[k];
   h2     = ele->e.s8->thick_node.a.dv[k];
   h2    /= 2.0;
/*   h2 *= condfac; not tested yet*/
   
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
/*=================================================== integration loops */
ngauss=0;
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
      /*--------------------------------- make height at gaussian point */
      s8_xint(&hhi,hte,funct,iel);
      /*---------- built metrics at current and reference configuration */
      s8_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,&detsmr,funct,deriv,iel,a3kvpr,0);
      s8_tvmr(xcure,a3c,akovc,akonc,amkovc,amkonc,&detsmc,funct,deriv,iel,a3kvpc,0);
      /*------------------------- make h as cross product in ref config.
                                       to get area da on shell mid surf */
      h[0] = akovr[1][0]*akovr[2][1] - akovr[2][0]*akovr[1][1];
      h[1] = akovr[2][0]*akovr[0][1] - akovr[0][0]*akovr[2][1];
      h[2] = akovr[0][0]*akovr[1][1] - akovr[1][0]*akovr[0][1];
      /*------------------------------------- make director unit lenght 
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);
      /*------------------------------------------------ clear stresses */
      for (i=0; i<6; i++) stress[i]=0.0;
      /*============================== loop GP in thickness direction t */
      for (lt=0; lt<nit; lt++)
      {
         /*---------------------------- gaussian point and weight at it */
         e3   = data->xgpt[lt];
         fact = data->wgtt[lt];
         /*-------------------- basis vectors and metrics at shell body */ 
         s8_tmtr(xrefe,a3r,e3,gkovr,gkonr,gmkovr,gmkonr,&detsrr,
                    funct,deriv,iel,condfac,0);

         s8_tmtr(xcure,a3c,e3,gkovc,gkonc,gmkovc,gmkonc,&detsrc,
                    funct,deriv,iel,condfac,0);
         /*--------------------------------- metric at gp in shell body */     
         s8_tvhe_lin(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,&detsrr,&detsrc,
                     amkovc,amkovr,akovc,akovr,a3kvpc,a3kvpr,e3);     
         /*------------------------------------------ call material law */
         s8_tmat(ele,mat,stress,strain,C,gmkovc,gmkonc,gmkovr,gmkonr,
                    gkovc,gkonc,gkovr,gkonr,detsrc,detsrr,e3,0);
         /*---------- calculates forces from stresses at gaussian point */
         s8_tfte(gp_stress,ngauss,stress,gkovr,akonr,gmkovr,gmkonr,amkovr,
                 amkonr,hhi,e3,fact,detsmr,detsrr);
      }/*========================================== end of loop over lt */
      /* calculate forces in respect to local/global coordinate systems */
      s8_tforce(gp_stress,ngauss,akovr,akonr,ele);
   /*------------------ set counter for total number of gaussian points */
   ngauss++;
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */
/*------------------------------------------- put stress to the element */
if (ele->e.s8->forces.fdim <= kstep) 
{
   am4redef(&(ele->e.s8->forces),
              ele->e.s8->forces.fdim+3,
              ele->e.s8->forces.sdim,
              ele->e.s8->forces.tdim,
              ele->e.s8->forces.fodim);
}
ele_stress = ele->e.s8->forces.a.d3;
ngauss = nir*nis;
for (i=0; i<ngauss; i++)
for (j=0; j<18; j++)
ele_stress[kstep][j][i] = gp_stress[j][i];
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_stress */





/*----------------------------------------------------------------------*
 | make redundant stresses to the shell8 elements         m.gee 12/01   |
 *----------------------------------------------------------------------*/
void s8_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      INTRA     *actintra,
                      int        kstep)
{
int          i,j,k;
int          imyrank;
int          inprocs;
ELEMENT     *actele;
ARRAY        mpi_buffer;
double     **buffer;

#ifdef DEBUG 
dstrc_enter("s8_stress_reduce");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
buffer = amdef("tmp",&mpi_buffer,18,MAXGAUSS,"DA");
/*----------------------------------------------------------------------*/
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   /*-------------------- there could be other elements in here as well */
   if (actele->eltyp != el_shell8) continue;
   /* check the size of the array to store stresses, should be at least */
   /*                                         of size kstep*18*MAXGAUSS */
   if (actele->e.s8->forces.fdim <= kstep)
   {
      am4redef(&(actele->e.s8->forces),
                 kstep+3,
                 actele->e.s8->forces.sdim,
                 actele->e.s8->forces.tdim,
                 actele->e.s8->forces.fodim);
   }
   /*--- the owner of the element broadcasts the stresses to the others */
   if (actele->proc == imyrank)
   {
      for (k=0; k<18; k++)
      for (j=0; j<MAXGAUSS; j++)
      buffer[k][j] = actele->e.s8->forces.a.d3[kstep][k][j];
   }
#ifdef PARALLEL 
   MPI_Bcast(buffer[0],18*MAXGAUSS,MPI_DOUBLE,actele->proc,actintra->MPI_INTRA_COMM);
/*   MPI_Bcast(actele->e.s8->energy.a.da[0],3*MAXGAUSS,MPI_DOUBLE,actele->proc,actintra->MPI_INTRA_COMM);*/
#endif
   if (actele->proc != imyrank)
   {
      for (k=0; k<18; k++)
      for (j=0; j<MAXGAUSS; j++)
      actele->e.s8->forces.a.d3[kstep][k][j] = buffer[k][j];
   }
} /* end of (i=0; i<actfield->numele; i++) */
/*----------------------------------------------------------------------*/
amdel(&mpi_buffer);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_stress_reduce */


/*----------------------------------------------------------------------*
 | strains nach bischoff seite 37                          m.gee 6/02   |
 *----------------------------------------------------------------------*/
void s8_strains_res(double **akovr, double **akovc, double **a3kvpr, double **a3kvpc,
                   double hh, double *strains)
{
int          i,j,k;
double       a1r[3],a1c[3],a2r[3],a2c[3],a3r[3],a3c[3],a31r[3],a31c[3],a32r[3],a32c[3];
#ifdef DEBUG 
dstrc_enter("s8_strains_res");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<12; i++) strains[i]=0.0;
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++) a1r[i] = akovr[i][0];
for (i=0; i<3; i++) a1c[i] = akovc[i][0];
for (i=0; i<3; i++) a2r[i] = akovr[i][1];
for (i=0; i<3; i++) a2c[i] = akovc[i][1];
for (i=0; i<3; i++) a3r[i] = akovr[i][2];
for (i=0; i<3; i++) a3c[i] = akovc[i][2];
for (i=0; i<3; i++) a31r[i] = a3kvpr[i][0];
for (i=0; i<3; i++) a31c[i] = a3kvpc[i][0];
for (i=0; i<3; i++) a32r[i] = a3kvpr[i][1];
for (i=0; i<3; i++) a32c[i] = a3kvpc[i][1];
/*----------------------------------------------------------------------*/
/* a11 */
for (i=0; i<3; i++) strains[0] += a1c[i]*a1c[i] - a1r[i]*a1r[i];
strains[0] /= 2.0;
/* a12 */
for (i=0; i<3; i++) strains[1] += a1c[i]*a2c[i] - a1r[i]*a2r[i];
strains[1] /= 2.0;
/* a13 */
for (i=0; i<3; i++) strains[2] += a1c[i]*a3c[i] - a1r[i]*a3r[i];
strains[2] /= 2.0;
/* a22 */
for (i=0; i<3; i++) strains[3] += a2c[i]*a2c[i] - a2r[i]*a2r[i];
strains[3] /= 2.0;
/* a23 */
for (i=0; i<3; i++) strains[4] += a2c[i]*a3c[i] - a2r[i]*a3r[i];
strains[4] /= 2.0;
/* a33 */
for (i=0; i<3; i++) strains[5] += a3c[i]*a3c[i] - a3r[i]*a3r[i];
strains[5] /= 2.0;
/* b11 */
for (i=0; i<3; i++) strains[6] += 2.0*a1c[i]*a31c[i] - 2.0*a1r[i]*a31r[i];
strains[6] /= hh;
/* b12 */
for (i=0; i<3; i++) strains[7] += a1c[i]*a32c[i] + a2c[i]*a31c[i] - a1r[i]*a32r[i] - a2r[i]*a31r[i];
strains[7] /= hh;
/* b13 */
for (i=0; i<3; i++) strains[8] += a31c[i]*a3c[i] - a31r[i]*a3r[i];
strains[8] /= hh;
/* b22 */
for (i=0; i<3; i++) strains[9] += 2.0*a2c[i]*a32c[i] - 2.0*a2r[i]*a32r[i];
strains[9] /= hh;
/* b23 */
for (i=0; i<3; i++) strains[10] += a32c[i]*a3c[i] - a32r[i]*a3r[i];
strains[10] /= hh;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_strains_res */
#endif


