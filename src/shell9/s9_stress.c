/*!----------------------------------------------------------------------
\file
\brief calculation of the element stresses
   contains the routines:
   - s9_stress:  calculates the stresses at the GP, extrapolates them to
                 nodal points and stores it in 'ele_stress -ARRAY'
   - s9_stress_reduce: make redundant stresses to the shell9 elements

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 |                                                          sh 10/02    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;



/*!----------------------------------------------------------------------
\brief calculates the element stresses                                      

<pre>                     m.gee 12/01             modified by    sh 01/03
This routine calculates the element stresses at the gauspoints. This is
done with the green-lagrange strains (calculated from metrics in current
and reference configuration) an the constitutive law (-> PK_II stresses).
The PK_II stresses are then transformed to physical stresses, either in
"XYZ", "RST" or "RST_ortho" coordinate system. The stresses at the GPs 
are then extrapolated to the nodal points using the shape functions.
</pre>
\param  ELEMENT   *ele   (i/o) the element structure -> 'ele->e.s9->stresses.a.d3'
\param  S9_DATA   *data   (i)  element integration data
\param  MATERIAL  *mat    (i)  the material structure
\param  INT        kintyp (i)  kintyp=0: geo_lin; =1: upd_lagr; =2: tot_lagr 
\param  INT        kstep  (i)  actual step in nonlinear analysis
\param  INT        init   (i)  init=1 -> init phase / init=0 -> calc. phase / init=-1 -> uninit phase

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9_stress(ELEMENT      *ele,
               S9_DATA      *data,
               MATERIAL     *mat,
               INT           kintyp,   /* typ of kinematic formulation */
               INT           kstep,
               INT           init)
{
INT                 i,j,k,l,kl,ml;

INT                 nir,nis,nit;
INT                 lr,ls,lt;
INT                 iel;
INT                 nd;
INT                 ngauss;

INT                 istore = 0;                             /* controls storing of new stresses to wa */
INT                 newval = 1;                             /* controls evaluation of new stresses    */

INT                 ip;                                     /* actual integration point */
INT                 actlay;                                 /* actual layer */
INT                 nir_x,nis_x,nit_x;                      /*order of extrapolation -> depends on quad4/8/9*/
INT                 gaussperm4[4] = {3,1,0,2};
INT                 gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
DOUBLE              P_u,P_o,m,N_u,N_o;                      /*for interpolation in thickness direction*/
DOUBLE              gp_u[6][9], gp_o[6][9];                 /*interpolated values in plane before interpolation in thickness */
DOUBLE              strK;                                   /*stress at a Node -> extrapolated*/
DOUBLE              strGP[9];                               /*stress at GPs in one plane*/
DIS_TYP             distyp;
INT                 ID_stress,ID_stress_perm;               /*ID, on which the calculated stress has to be written -> gp_stress*/
INT                 sum_lay;                                /*sum of all layers: kinematic and material*/
INT                 num_mlay;                               /* number of material layers to actual kinematic layer */  
INT                 num_klay;                               /* number of kinematic layers to this element*/
INT                 numdf;                                  /* ndofs per node to this element */
DOUBLE             *klayhgt;                                /* hight of kinematic layer in percent of total thicknes of element*/
DOUBLE             *mlayhgt;                                /* hight of material layer in percent of adjacent kinematic layer*/
MULTIMAT           *actmultimat;                            /* material of actual material layer */
INT                 rot_axis;                               /* rotation axis of laminat (1=x-; 2=y-; 3=z-axis) */
DOUBLE              phi;                                    /* angle of rotation about rot_axis */

DOUBLE              e1,e2,e3;
DOUBLE              facr,facs,fact;

DOUBLE              condfac;                                /* sdc conditioning factor */
DOUBLE              h2;                                     /* half nodal height */
DOUBLE              hgt;                                    /* element thickness */
DOUBLE              hte[MAXNOD_SHELL9];                     /* element thickness at nodal points */

DOUBLE              det_dummy;
DOUBLE              detsmr;
DOUBLE              detsmc;
DOUBLE              detsrr;
DOUBLE              detsrc;
 
DOUBLE              h[3];                                   /* working array */
DOUBLE              da;                                     /* area on mid surface */

DOUBLE              stress[6], stress_r[12];                /* stress and stress resultants */
DOUBLE              strain[6];                              /* strains */

static ARRAY        strK_a;      static DOUBLE **gp_strK;   /* element array for stresses on nodal points */

static ARRAY        stress_a;    static DOUBLE **gp_stress; /* element array for stresses on gaussian points */
                                      DOUBLE ***ele_stress; /* pointer to array of stress history in element */
/*static ARRAY        forces_a;    static DOUBLE **gp_forces; /* element array for forces (stress resultants) on gaussian points */
/*                                      DOUBLE ***ele_forces; /* pointer to array of stress history in element */

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

#ifdef DEBUG 
dstrc_enter("s9_stress");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
if (init==1)
{
iel       = MAXNOD_SHELL9; /* maximum number of nodes for this type of shell */

xrefe     = amdef("xrefe"  ,&xrefe_a,3,MAXNOD_SHELL9,"DA");       
xcure     = amdef("xcure"  ,&xcure_a,3,MAXNOD_SHELL9,"DA");       
a3r       = am4def("a3r"    ,&a3r_a,3,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         
a3c       = am4def("a3c"    ,&a3c_a,3,MAXNOD_SHELL9,MAXKLAY_SHELL9,0,"D3");         

a3kvpr    = am4def("a3kvpr" ,&a3kvpr_a,3,2,MAXKLAY_SHELL9,0,"D3");         
a3kvpc    = am4def("a3kvpc" ,&a3kvpc_a,3,2,MAXKLAY_SHELL9,0,"D3");        

funct     = amdef("funct"  ,&funct_a,MAXNOD_SHELL9,1,"DV");       
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_SHELL9,"DA");       

akovr     = am4def("akovr"  ,&akovr_a,3,3,MAXKLAY_SHELL9,0,"D3"); ;         
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
  
C         = amdef("C"      ,&C_a   ,6 ,6                    ,"DA");             
D         = amdef("D"      ,&D_a   ,12,12                   ,"DA");           

gp_strK   = amdef("gp_strK"  ,&strK_a,   6,MAXNODESTRESS_SHELL9,"DA");        
gp_stress = amdef("gp_stress",&stress_a, 6,MAXGAUSS,"DA");        
/*gp_forces = amdef("gp_forces",&forces_a,18,MAXGAUSS,"DA"); */       

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

amdel(&C_a);
amdel(&D_a);

amdel(&strK_a);
amdel(&stress_a);
/*amdel(&forces_a);*/
   
goto end;
}
/*----------------------------------------------------------------------*/
/* calculation phase                                                    */
/*----------------------------------------------------------------------*/
/*------------------------------ count number of layers to this element */
num_klay = ele->e.s9->num_klay;         /* number of kinematic layers to this element*/
numdf = ele->e.s9->numdf;               /* ndofs per node to this element */
sum_lay = 0;
for (kl=0; kl<num_klay; kl++) 
{ 
   num_mlay = ele->e.s9->kinlay[kl].num_mlay;
   sum_lay += num_mlay;
}
/*-------------------------------------------- init the gaussian points */
s9intg(ele,data,0);
amzero(&stress_a);
/*amzero(&forces_a);*/
/*----------------------------------------------- integrationsparameter */
nir     = ele->e.s9->nGP[0];
nis     = ele->e.s9->nGP[1];
nit     = ele->e.s9->nGP[2];
iel     = ele->numnp;
nd      = iel*numdf; 
condfac = ele->e.s9->sdc;
a3ref   = ele->e.s9->a3ref.a.da;
/*----------------------------------------------------- geometry update */
for (kl=0; kl<num_klay; kl++) /*loop over all kinematic layers*/
{  
  klayhgt = ele->e.s9->klayhgt;   /* hgt of kinematic layer on percent of total thickness of shell */
  for (k=0; k<iel; k++)           /*loop over all nodes per layer*/
  {
     hte[k] = ele->e.s9->thick_node.a.dv[k];
     /*if (ele->e.s9->dfield == 0)      /*Layerthicknes, norm(a3L) = HL */
     h2 = ele->e.s9->thick_node.a.dv[k] * klayhgt[kl]/100. * condfac;
     /*h2 = 0.5*h2; /*A3_IST_EINHALB halber Direktor*/
     h2 = A3FAC_SHELL9 * h2;
     /*else if (ele->e.s9->dfield == 1) /*half of shell thickness, norm(a3) = H/2*/
     /*  h2 = ele->e.s9->thick_node.a.dv[k]/2. * condfac;*/
 
     a3r[0][k][kl] = a3ref[0][k] * h2;
     a3r[1][k][kl] = a3ref[1][k] * h2;
     a3r[2][k][kl] = a3ref[2][k] * h2;

     xrefe[0][k] = ele->node[k]->x[0];
     xrefe[1][k] = ele->node[k]->x[1];
     xrefe[2][k] = ele->node[k]->x[2];

     xcure[0][k] = xrefe[0][k] + ele->node[k]->sol.a.da[kstep][0];
     xcure[1][k] = xrefe[1][k] + ele->node[k]->sol.a.da[kstep][1];
     xcure[2][k] = xrefe[2][k] + ele->node[k]->sol.a.da[kstep][2];
 
     a3c[0][k][kl] = a3r[0][k][kl]  + ele->node[k]->sol.a.da[kstep][3*kl+3];
     a3c[1][k][kl] = a3r[1][k][kl]  + ele->node[k]->sol.a.da[kstep][3*kl+4];
     a3c[2][k][kl] = a3r[2][k][kl]  + ele->node[k]->sol.a.da[kstep][3*kl+5];
  } /*end loop over nodes*/
} /*end loop over kinematic layers*/
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
      s9_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      /*----- make height at gaussian point -> variable thickness within the element is normally not allowed */
      s9_xint(&hgt,hte,funct,iel);
      /*------------------------------------- metrics at gaussian point */
      s9_tvmr(xrefe,a3r,akovr,akonr,amkovr,amkonr,akovh,akonh,amkovh,amkonh,
              &det_dummy,funct,deriv,iel,a3kvpr,num_klay);
      s9_tvmr(xcure,a3c,akovc,akonc,amkovc,amkonc,akovh,akonh,amkovh,amkonh,
              &det_dummy,funct,deriv,iel,a3kvpc,num_klay);
      /*------------------------- make h as cross product in ref config.
                                       to get area da on shell mid surf */
      h[0] = akovr[1][0][0]*akovr[2][1][0] - akovr[2][0][0]*akovr[1][1][0];
      h[1] = akovr[2][0][0]*akovr[0][1][0] - akovr[0][0][0]*akovr[2][1][0];
      h[2] = akovr[0][0][0]*akovr[1][1][0] - akovr[1][0][0]*akovr[0][1][0];
      /*------------------------------------- make director unit lenght 
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);
      /*------------------------------------------------ clear stresses */
      for (i=0; i<6; i++) stress[i]=0.0;

/*---------------------------loop over all kinematic layers (num_klay) */
      actlay = 0;
      for (kl=0; kl<num_klay; kl++)   /*loop over all kinematic layers */
      {
         num_mlay = ele->e.s9->kinlay[kl].num_mlay;
         mlayhgt  = ele->e.s9->kinlay[kl].mlayhgt;   /* hgt of material layer in percent of this kin layer */
         for(ml=0; ml<num_mlay; ml++) /*loop over all material layers of aktuel kinematic layer*/
         {
            for (lt=0; lt<nit; lt++)  /*loop in t-direction */
            {
               /*---------------------------- gaussian point and weight at it */
               e3   = data->xgpt[lt];
               fact = data->wgtt[lt];
               /*-------------------- basis vectors and metrics at shell body */ 
               s9_tmtr(xrefe,a3r,e3,gkovr,gkonr,gmkovr,gmkonr,&detsmr,
                          funct,deriv,iel,akovr,a3kvpr,hgt,klayhgt,mlayhgt,
                          num_klay,num_mlay,kl,ml,condfac);

               s9_tmtr(xcure,a3c,e3,gkovc,gkonc,gmkovc,gmkonc,&detsmc,
                          funct,deriv,iel,akovc,a3kvpc,hgt,klayhgt,mlayhgt,
                          num_klay,num_mlay,kl,ml,condfac);
               /*--------------------------------- metric at gp in shell body */
               s9_tvhe(gmkovr,gmkovc,gmkonr,gmkonc,gkovr,gkovc,&detsrr,&detsrc,
                       amkovc,amkovr,akovc,akovr,a3kvpc,a3kvpr,e3,kintyp,hgt,
                       klayhgt,mlayhgt,num_klay,num_mlay,kl,ml,condfac);
               /*------------------------------------------ call material law */
               actmultimat = &(multimat[ele->e.s9->kinlay[kl].mmatID[ml]-1]);
               rot_axis = mat->m.multi_layer->kinlay[kl].rot[ml];
               phi = mat->m.multi_layer->kinlay[kl].phi[ml];
               
               ip = 2* ngauss + lt;
               s9_call_mat(ele,actmultimat,stress,strain,C,gmkovc,gmkonc,gmkovr,gmkonr,
                           gkovc,gkonc,gkovr,gkonr,rot_axis,phi,ip,actlay,istore,newval);
               /*- calculates physical stresses at gaussian point in respect to local/global coordinat system */
               ID_stress = ngauss + (2*actlay + lt) * (nir*nis);  /*write gp's layerwise*/
               s9_tstress(gp_stress,stress,ID_stress,gkovr,ele);
            }/*========================================== end of loop over lt */            
            actlay++;
         }/*======= end of loop over all material layers of aktual kinematic layer*/
      }/*============================== end loop over all kinematic layers */
      /*------------------ set counter for number of gaussian points in plane */
      ngauss++;
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */

/*---------------------- put physical stresses to the element */
if (ele->e.s9->stresses.fdim <= kstep) 
{
   am4redef(&(ele->e.s9->stresses),
              ele->e.s9->stresses.fdim+3,
              ele->e.s9->stresses.sdim,
              ele->e.s9->stresses.tdim,
              ele->e.s9->stresses.fodim);
}
ele_stress = ele->e.s9->stresses.a.d3;

/*extrapolate stresses from GP to nodes*/
if      (ngauss == 4) distyp = quad4;  /*2x2 gp*/
else if (ngauss == 9) distyp = quad9;  /*3x3 gp*/


/*set order for extrapolation*/
if (ele->distyp == quad4)
{ 
   nir_x = 2;
   nis_x = 2;
   nit_x = 2;
   s9intg_str(data,4);     /*integration parameters for stress extrapolation*/
}
else if (ele->distyp == quad8 || ele->distyp == quad9)
{ 
   nir_x = 3;
   nis_x = 3;
   nit_x = 2;
   s9intg_str(data,9);     /*integration parameters for stress extrapolation*/
}


/*--------------------------- loop over all layers (sum_lay) */
for (actlay=0; actlay<sum_lay; actlay++)   /*loop over all layers */
{
   for (lt=0; lt<nit_x; lt++)  /*loop in t-direction */
   {
      ngauss = 0;
      for (lr=0; lr<nir_x; lr++) /*loop in r-direction*/
      {
         /*================================== gaussian point */
         e1   = data->xgpr[lr];
         if (e1 == 0.0) e1 = 0.0;
         else e1   = 1./e1;
         for (ls=0; ls<nis_x; ls++)
         {
            /*=============================== gaussian point */
            e2   = data->xgps[ls];
            if (e2 == 0.0) e2 = 0.0;
            else e2   = 1./e2;
            /*-------------------- shape functions at gp e1,e2 on mid surface */
            s9_funct_deriv(funct,deriv,e1,e2,distyp,0);

            /*interpolate in this plane to nodal points*/
            for (j=0; j<6; j++)    /*6-stress komponents*/
            {
               strK = 0.0;

               for (k=0; k<(nir*nis); k++)
               {
                 ID_stress = k + (2*actlay + lt) * (nir*nis);  /*write gp's layerwise*/
                 strGP[k] = gp_stress[j][ID_stress];
               }

               if      (distyp == quad4) 
                 for (k=0; k<4; k++)    strK += funct[k]*strGP[gaussperm4[k]];
               else if (distyp == quad9) /*write midpoint value for quad8 as well ! */
                 for (k=0; k<9; k++)    strK += funct[k]*strGP[gaussperm9[k]];

               /*store value to finaly interpolate in thickness direction*/
               /*lower part*/
               if (lt == 0)  gp_u[j][ngauss] = strK;
               /*upper part*/
               if (lt == 1)  gp_o[j][ngauss] = strK;
            }
            ngauss++;           
         }/*=============================================== end of loop over ls */ 
      }/*================================================== end of loop over lr */ 
   }/*===================================================== end of loop over lt */
   /* now interpolate the gp_u and gp_o values to the nodes in thickness direction */
   /*interpolate the 6-stress komponents*/
   for (j=0; j<6; j++)
   {
      for (k=0; k<(nir_x*nis_x); k++)
      {
       P_u = gp_u[j][k];
       P_o = gp_o[j][k];
       m   = (P_o - P_u) / (2./sqrt(3.));
       N_u = P_u - m * ( 1. - 1./sqrt(3.));
       N_o = P_o + m * ( 1. - 1./sqrt(3.));    
       /*write to gp_strK array*/
       ID_stress = k + (2*actlay + 0) * (nir_x*nis_x);
       gp_strK[j][ID_stress] = N_u;           
       ID_stress = k + (2*actlay + 1) * (nir_x*nis_x);
       gp_strK[j][ID_stress] = N_o;           
      }
   }
}/*======================================================== end of loop over all layers */

/* now write the stresses at the nodes gp_strK back to ele_stress array */
/*--------------------------- loop over all layers (sum_lay) */
for (actlay=0; actlay<sum_lay; actlay++)   /*loop over all layers */
{
   for (lt=0; lt<nit_x; lt++)  /*loop in t-direction */
   {
      ngauss = 0;
      for (lr=0; lr<nir_x; lr++) /*loop in r-direction*/
      {
         for (ls=0; ls<nis_x; ls++)
         {
            for (j=0; j<6; j++)    /*6-stress komponents*/
            {
             ID_stress = ngauss + (2*actlay + lt) * (nir_x*nis_x);

             if      (ele->distyp == quad4) 
                  ID_stress_perm = gaussperm4[ngauss] + (2*actlay + lt) * (nir_x*nis_x);
             else if (ele->distyp == quad8 ||ele->distyp == quad9) 
                  ID_stress_perm = gaussperm9[ngauss] + (2*actlay + lt) * (nir_x*nis_x);

             ele_stress[kstep][j][ID_stress] = gp_strK[j][ID_stress_perm];
            }
            ngauss++;           
         }/*=============================================== end of loop over ls */ 
      }/*================================================== end of loop over lr */ 
   }/*===================================================== end of loop over lt */
}/*======================================================== end of loop over all layers */
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9_stress */



/*!----------------------------------------------------------------------
\brief make redundant stresses to the shell9 elements                                         

<pre>                     m.gee 12/01             modified by    sh 02/03
This routine make redundant stresses to the shell9 elements, which is only
necessary in parallel.
</pre>
\param *actfield    FIELD       (i)   my field
\param *actpart     PARTITION   (i)   my partition
\param *actintra    INTRA       (i)   my intra-communicator 
\param  kstep       INT         (i)   actual step in nonlinear analysis

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: shell9()   [s9_main.c]

*----------------------------------------------------------------------*/
void s9_stress_reduce(FIELD     *actfield,
                      PARTITION *actpart,
                      INTRA     *actintra,
                      INT        kstep)
{
INT          i,j,k;
INT          imyrank;
INT          inprocs;
ELEMENT     *actele;
ARRAY        mpi_buffer;
DOUBLE     **buffer;
ARRAY        mpi_buffer1;
DOUBLE     **buffer1;

#ifdef DEBUG 
dstrc_enter("s9_stress_reduce");
#endif
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
buffer = amdef("tmp" ,&mpi_buffer ,18,MAXGAUSS,"DA");
buffer1= amdef("tmp1",&mpi_buffer1,18,MAXGAUSS,"DA");
/*----------------------------------------------------------------------*/
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   /*-------------------- there could be other elements in here as well */
   if (actele->eltyp != el_shell9) continue;
   /* check the size of the array to store stresses, should be at least */
   /*                                         of size kstep*18*MAXGAUSS */
   /*---- stress resultants (~ Schnittgroessen) are not implemented yet */
/*   if (actele->e.s9->forces.fdim <= kstep)
   {
      am4redef(&(actele->e.s9->forces),
                 kstep+3,
                 actele->e.s9->forces.sdim,
                 actele->e.s9->forces.tdim,
                 actele->e.s9->forces.fodim);
   }
   /*------------- stresses at nodal points  ---------------------------*/
   if (actele->e.s9->stresses.fdim <= kstep) 
   {
      am4redef(&(actele->e.s9->stresses),
                 kstep+3,
                 actele->e.s9->stresses.sdim,
                 actele->e.s9->stresses.tdim,
                 actele->e.s9->stresses.fodim);
   }
   /*--- the owner of the element broadcasts the stresses to the others */
   if (actele->proc == imyrank)
   {
/*      for (k=0; k<18; k++)
      for (j=0; j<MAXGAUSS; j++)
      buffer[k][j]  = actele->e.s9->forces.a.d3[kstep][k][j];*/
      for (k=0; k<6; k++)
      for (j=0; j<MAXGAUSS; j++)
      buffer1[k][j] = actele->e.s9->stresses.a.d3[kstep][k][j];
   }
#ifdef PARALLEL 
/*   MPI_Bcast(buffer[0] ,18*MAXGAUSS,MPI_DOUBLE,actele->proc,actintra->MPI_INTRA_COMM);*/
   MPI_Bcast(buffer1[0],6*MAXGAUSS,MPI_DOUBLE,actele->proc,actintra->MPI_INTRA_COMM);
#endif
   if (actele->proc != imyrank)
   {
/*      for (k=0; k<18; k++)
      for (j=0; j<MAXGAUSS; j++)
      actele->e.s9->forces.a.d3[kstep][k][j]   = buffer[k][j];*/
      for (k=0; k<6; k++)
      for (j=0; j<MAXGAUSS; j++)
      actele->e.s9->stresses.a.d3[kstep][k][j] = buffer1[k][j];
   }
} /* end of (i=0; i<actfield->numele; i++) */
/*----------------------------------------------------------------------*/
/*amdel(&mpi_buffer);  /* for the forces*/
amdel(&mpi_buffer1);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_stress_reduce */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/


