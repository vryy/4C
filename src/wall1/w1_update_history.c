/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_update_history' which is the main routine
 used to update state variables related with E-M Int. Scheme

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*---------------------------------------------------------------------------*/
#ifdef GEMM
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*--------------------------------------------------------------------------*
 | Update of state variables related with E-M algorithm                     |
 *-------------------------------------------------------------------------*/
void w1_update_history (ELEMENT   *ele, 
                        W1_DATA   *data,
			MATERIAL  *mat)                   
{
INT                 i,j,k,a,b;        /* some loopers              */
INT                 nir,nis;          /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP           */
INT                 iel;              /* numnp to this element     */
INT                 nd;               /* dof of this element       */
INT                 ip;
INT                 istore;
const INT           numdf  = 2;       /* number dof per node          */
const INT           numeps = 4;       /* number of strain components  */

DOUBLE              fac,facm;         /* integration factor             */
DOUBLE              e1,e2,e3;         /*GP-coords                       */
DOUBLE              facr,facs,fact;   /* weights at GP                  */
DOUBLE              thick,density;    /* for calculation of mass matrix */
DOUBLE det;

static ARRAY    D_a;      /* material matrix */     
static DOUBLE **D; 
static ARRAY    stress_a; /* stress matrix (2.PK for total lagr.) */     
static DOUBLE  **stress; 
static ARRAY    xcure_a;  /* coords in current config. */     
static DOUBLE **xcure; 
static ARRAY    xrefe_a;  /* coords in referenz config. */     
static DOUBLE **xrefe; 
static ARRAY    funct_a;  /* shape functions */    
static DOUBLE  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static DOUBLE **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static DOUBLE **xjm;         
static ARRAY    boplin_a; /* linear B-operator */   
static DOUBLE **boplin;
static ARRAY    b_bar_a;  /* B_bar operator*/
static DOUBLE **b_bar; 
static ARRAY    mass_a;   /* Mass matrix*/
static DOUBLE **mass;
static ARRAY    F_a;      /* deformation gradient */   
static DOUBLE  *F; 
static ARRAY    strain_a; /* strain (Green-Lagr for total lagr.) */   
static DOUBLE  *strain; 


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_update_history");
#endif
/*----------------------------------------------------------------------*/
/* init phase                                                           */
/*----------------------------------------------------------------------*/
  istore = 0;
  xrefe     = amdef("xrefe"  ,&xrefe_a,2,MAXNOD_WALL1,"DA");
  xcure     = amdef("xcure"  ,&xcure_a,2,MAXNOD_WALL1,"DA");       
  funct     = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1,"DV");       
  deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_WALL1,"DA");       
  D         = amdef("D"      ,&D_a   ,6,6             ,"DA");           
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf     ,"DA");           
  boplin    = amdef("boplin" ,&boplin_a ,numeps,(numdf*MAXNOD_WALL1),"DA"); 
  b_bar     = amdef("b_bar"  ,&b_bar_a , numeps,(numdf*MAXNOD_WALL1),"DA");
  mass      = amdef("mass"   ,&mass_a ,(numdf*MAXNOD_WALL1),(numdf*MAXNOD_WALL1),"DA" );  
  F         = amdef("F"      ,&F_a ,numeps,1,"DV"); 
  strain    = amdef("strain" ,&strain_a ,numeps,1,"DV"); 
  stress    = amdef("stress" ,&stress_a ,numeps,numeps,"DA"); 
  
/*--------------------------------------- get integration parameters ---*/
w1intg(ele,data,1);
w1_getdensity(mat,&density);
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.w1->nGP[0];
nis     = ele->e.w1->nGP[1];
iel     = ele->numnp;
nd      = numdf * iel;
thick = ele->e.w1->thick;
/*-------------------------------- initialize element energy and momenta*/
ele->e.w1->strain_energy = 0.0;
ele->e.w1->kinetic_energy = 0.0;
ele->e.w1->angular_momentum = 0.0;
ele->e.w1->linmom[0] = 0.0;
ele->e.w1->linmom[1] = 0.0;

for(i=0; i<nd; i++)
  for(j=0; j<nd; j++)
  mass[i][j] = 0.0;
/*------------------------------ geometry update for total lagrange ----*/
  for (k=0; k<iel; k++)
  {
    xrefe[0][k]= ele->node[k]->x[0];                         /* x-direction */
    xrefe[1][k]= ele->node[k]->x[1];                         /* y-direction */
    xcure[0][k]= xrefe[0][k] + ele->node[k]->sol.a.da[0][0]; /* x-direction */
    xcure[1][k]= xrefe[1][k] + ele->node[k]->sol.a.da[0][1]; /* y-direction */
  }

ip = -1;
for (lr=0; lr<nir; lr++)
{
   /*================================ gaussian point and weight at it ===*/
   e1   = data->xgrr[lr];
   facr = data->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      ip++;
      /*============================= gaussian point and weight at it ===*/
      e2   = data->xgss[ls];
      facs = data->wgts[ls];
      
      /*-------------------------- shape functions and their derivatives */
      w1_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
      
      /*-------------------------------------- compute jacobian matrix --*/       
      w1_jaco (funct,deriv,xjm,&det,ele,iel); 
      /*------------------------------------------- integration factor --*/ 
      fac = facr * facs * det * thick;  
      /*------------------------------------------ compute mass matrix---*/
      facm = fac * density;
       for (a=0; a<iel; a++)
       {
        for (b=0; b<iel; b++)
        {
         mass[2*a][2*b]     += facm * funct[a] * funct[b]; /* a,b even */
         mass[2*a+1][2*b+1] += facm * funct[a] * funct[b]; /* a,b odd  */
	}
       }
      /*------------------------------------- calculate operator Blin ---*/
      amzero(&boplin_a);
      w1_boplin(boplin,deriv,xjm,det,iel);
      
      /*------------------ calculate defgrad F, Green-Lagrange-strain ---*/
      amzero(&F_a);
      amzero(&strain_a);
      w1_defgrad(F,strain,xrefe,xcure,boplin,iel);
      
      /*------------------------------------------- call material law ---*/
      amzero(&stress_a);
      amzero(&D_a);
      w1_call_matgeononl(mat,ele->e.w1->wtype,strain,stress,D,numeps);
      
      /*----------------------------------- update b_bar and 2PK stresses*/
#ifdef GEMM
      amzero(&b_bar_a);
      w1_history(ele,b_bar,boplin,stress,F,numeps,nd,ip);      
#endif
      /*----------------------------------------- calculate strain energy*/
      w1_strain_energy(ele,stress,strain,fac);
      /*-----------------------------------------------------------------*/                   
   }/*============================================== end of loop over ls */ 
}/*================================================= end of loop over lr */
      /*---------------------------- calculate kinetic energy and momenta*/
     w1_kinetic_energy(ele,mass);

/*------------------------------------------------------cleaning-up phase*/
   amdel(&xrefe_a);   
   amdel(&xcure_a); 
   amdel(&funct_a);
   amdel(&deriv_a);
   amdel(&D_a);
   amdel(&xjm_a);
   amdel(&boplin_a);
   amdel(&b_bar_a);
   amdel(&mass_a);
   amdel(&F_a);
   amdel(&strain_a);
   amdel(&stress_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_update_history */
#endif /*GEMM*/
/*! @} (documentation module close)*/

