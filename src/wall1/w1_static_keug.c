#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | integration of nonlinear stiffness keug for wall1 element ah 06/02    |
 *----------------------------------------------------------------------*/
void w1static_keug(ELEMENT   *ele, 
                   W1_DATA   *data, 
                   MATERIAL  *mat,
                   ARRAY     *estif_global, /* global vector for stiffness ke+ku+kg */
                   ARRAY     *emass_global, /* global vector for mass */
                   double    *force,        /* global vector for internal forces (initialized!)*/
                   int        init)
{
int                 i,j,k,a,b;        /* some loopers              */
int                 nir,nis;          /* num GP in r/s/t direction */
int                 lr, ls;           /* loopers over GP           */
int                 iel;              /* numnp to this element     */
int                 nd;               /* dof of this element       */
int                 ip;
int                 istore = 0;          /* controls storing of new stresses to wa */
int                 newval = 0;          /* controls evaluation of new stresses    */
const int           numdf  = 2;          /* number dof per node                    */
const int           numeps = 4;          /* number of strain components            */

double              fac,facm;              /* integration factor                        */
double              e1,e2,e3;         /*GP-coords                                  */
double              facr,facs,fact;   /* weights at GP                             */
double              thick,density;    /* for calculation of mass matrix            */
int                 imass;            /* flag -> mass matrix is/is not to evaluate */

static ARRAY    D_a;      /* material matrix */     
static double **D; 
static ARRAY    stress_a; /* stress matrix (2.PK for total lagr.) */     
static double  **stress; 
        
static ARRAY    xcure_a;  /* coords in current config. */     
static double **xcure; 
static ARRAY    xrefe_a;  /* coords in referenz config. */     
static double **xrefe; 

static ARRAY    funct_a;  /* shape functions */    
static double  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static double **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static double **xjm;         
static ARRAY    boplin_a; /* linear B-operator */   
static double **boplin; 

static ARRAY    F_a;      /* deformation gradient */   
static double  *F; 
static ARRAY    strain_a; /* strain (Green-Lagr for total lagr.) */   
static double  *strain; 
static ARRAY    keu_a;    /* element stiffness matrix keu */
static double **keu;
static ARRAY    kg_a;    /* element stiffness matrix kg */
static double **kg;      


static double **estif;    /* element stiffness matrix keug */
static double **emass;    /* mass matrix */

double fie[18];

double det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1static_keug");
#endif
/*----------------------------------------------------------------------*/
/* init phase        (init=1)                                           */
/*----------------------------------------------------------------------*/
istore = 0;

if (init==1)
{
  xrefe     = amdef("xrefe"  ,&xrefe_a,2,MAXNOD_WALL1,"DA");
  xcure     = amdef("xcure"  ,&xcure_a,2,MAXNOD_WALL1,"DA");       
  funct     = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1,"DV");       
  deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_WALL1,"DA");       
  D         = amdef("D"      ,&D_a   ,6,6             ,"DA");           
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf     ,"DA");           
  boplin    = amdef("boplin" ,&boplin_a ,numeps,(numdf*MAXNOD_WALL1),"DA"); 
  F         = amdef("F"      ,&F_a ,numeps,1,"DV"); 
  strain    = amdef("strain" ,&strain_a ,numeps,1,"DV"); 
  stress    = amdef("stress" ,&stress_a ,numeps,numeps,"DA"); 
  kg        = amdef("kg"     ,&kg_a,2*MAXNOD_WALL1,2*MAXNOD_WALL1,"DA");
  keu       = amdef("keu"    ,&keu_a,2*MAXNOD_WALL1,2*MAXNOD_WALL1,"DA");
          
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
   amdel(&xrefe_a);   
   amdel(&xcure_a);   
   amdel(&funct_a);
   amdel(&deriv_a);
   amdel(&D_a);
   amdel(&xjm_a);
   amdel(&boplin_a);
   amdel(&F_a);
   amdel(&strain_a);
   amdel(&stress_a);
   amdel(&kg_a);
   amdel(&keu_a);
   goto end;  
}

/*----------------------------------------------------------------------*/
/* update phase  for material nonlinearity      (init=2)                */
/*----------------------------------------------------------------------*/
else if(init==2)
{
  istore = 1;             /*-- material law is called with this flag----*/
}
/*--------------------------------------------------------------------- */
/* calculation phase    (init=0)                                        */
/*--------------------------------------------------------------------- */
/*--------------------------------------- get integration parameters ---*/
w1intg(ele,data,1);
/*------------------------------------ check calculation of mass matrix */
if (emass_global) 
{
   imass = 1;
   amzero(emass_global);
   emass = emass_global->a.da;
   w1_getdensity(mat,&density);
   thick = ele->e.w1->thick;
} 
else 
{
   imass   = 0;
   emass   = NULL;
   density = 0.0;
   thick = ele->e.w1->thick;
/*   thick   = NULL;*/
}   
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
iel     = ele->numnp;
nd      = numdf * iel;
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
amzero(&kg_a);
amzero(&keu_a);
if(force)
{
  for(i=0; i<nd; i++)
  {
  force[i]=0.0;
  }
}
estif     = estif_global->a.da;
/*------------------------------ geometry update for total lagrange ----*/

if(ele->e.w1->kintype==total_lagr)
{
  for (k=0; k<iel; k++)
  {
    xrefe[0][k]= ele->node[k]->x[0];                         /* x-direction */
    xrefe[1][k]= ele->node[k]->x[1];                         /* y-direction */
    xcure[0][k]= xrefe[0][k] + ele->node[k]->sol.a.da[0][0]; /* x-direction */
    xcure[1][k]= xrefe[1][k] + ele->node[k]->sol.a.da[0][1]; /* y-direction */
  }
}
else
{
dserror("action unknown");
}  
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
      /*------------------------------------ integration factor  -------*/ 
      fac = facr * facs * det * thick;                        
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
      /*----------------------------------- calculate operator Blin  ---*/
      amzero(&boplin_a);
      w1_boplin(boplin,deriv,xjm,det,iel);
      /*----------------- calculate defgrad F, Green-Lagrange-strain ---*/
      amzero(&F_a);
      amzero(&strain_a);
      w1_defgrad(F,strain,xrefe,xcure,boplin,iel);
      /*------------------------------------------ call material law ---*/
      amzero(&stress_a);
      amzero(&D_a);
      w1_call_matgeononl(ele,mat,ele->e.w1->wtype,boplin,xjm,
                          ip,strain,stress,D,istore,numeps);
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*---------------------- geometric part of stiffness matrix kg ---*/
        w1_kg(kg,boplin,stress,fac,nd,numeps);
      /*------------------ elastic+displacement stiffness matrix keu ---*/
        w1_keu(keu,boplin,D,F,fac,nd,numeps);
      /*------------------------------- stiffness matrix keug=keu+kg ---*/
       for (a=0; a<nd; a++)
       {  
         for (b=0; b<nd; b++)
         {
           estif[a][b] = kg[a][b] + keu[a][b];
         } 
       }
      /*--------------- nodal forces fi from integration of stresses ---*/
        if (force) w1_fint(stress,F,boplin,force,fac,nd);                    
      }
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1static_keug */


#endif
