/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1static_keug' which forms the nonlinear
       stiffnes keug for wall1 element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
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
 | integration of nonlinear stiffness keug for wall1 element ah 06/02    |
 *----------------------------------------------------------------------*/
void w1static_keug(ELEMENT   *ele, 
                   W1_DATA   *data, 
                   MATERIAL  *mat,
                   ARRAY     *estif_global, /* global vector for stiffness ke+ku+kg */
                   ARRAY     *emass_global, /* global vector for mass */
                   DOUBLE    *force,        /* global vector for internal forces (initialized!)*/
                   INT        init)
{
INT                 i,k,a,b;          /* some loopers              */
INT                 nir,nis;          /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP           */
INT                 iel;              /* numnp to this element     */
INT                 nd;               /* dof of this element       */
INT                 ip;
INT                 intc;                /* "integration case" for tri-element     */
INT                 istore = 0;          /* controls storing of new stresses to wa */
INT                 newval = 0;          /* controls evaluation of new stresses    */
const INT           numdf  = 2;          /* number dof per node                    */
const INT           numeps = 4;          /* number of strain components            */

DOUBLE              fac,facm;            /* integration factors                    */
DOUBLE              e1,e2;               /*GP-coords                               */
DOUBLE              facr,facs;           /* weights at GP                          */
DOUBLE              thick,density;       /* for calculation of mass matrix         */
INT                 imass;               /* flag -> mass matrix is to evaluate?    */

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
static ARRAY    b_bar_a;     /* B_bar operator */ 
static DOUBLE **b_bar; 
static ARRAY    int_b_bar_a; /* Interpolated B_bar operator   */
static DOUBLE **int_b_bar;
static ARRAY    F_a;      /* deformation gradient */   
static DOUBLE  *F; 
static ARRAY    strain_a; /* strain (Green-Lagr for total lagr.) */   
static DOUBLE  *strain; 
static ARRAY    keu_a;    /* element stiffness matrix keu */
static DOUBLE **keu;
static ARRAY    kg_a;    /* element stiffness matrix kg */
static DOUBLE **kg;      


static DOUBLE **estif;    /* element stiffness matrix keug */
static DOUBLE **emass;    /* mass matrix */

DOUBLE det;

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
  b_bar     = amdef("b_bar"  ,&b_bar_a,numeps,(numdf*MAXNOD_WALL1),"DA");
  int_b_bar = amdef("int_b_bar"  ,&int_b_bar_a,numeps,(numdf*MAXNOD_WALL1),"DA");
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
   amdel(&b_bar_a);
   amdel(&int_b_bar_a);
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
for (k=0; k<iel; k++)
{
 xrefe[0][k]= ele->node[k]->x[0];                         /* x-direction */
 xrefe[1][k]= ele->node[k]->x[1];                         /* y-direction */
 xcure[0][k]= xrefe[0][k] + ele->node[k]->sol.a.da[0][0]; /* x-direction */
 xcure[1][k]= xrefe[1][k] + ele->node[k]->sol.a.da[0][1]; /* y-direction */
}
if(ele->e.w1->kintype==updated_lagr)
{
 dserror("action unknown");
}  
/*------- get integraton data ---------------------------------------- */
switch (ele->distyp)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   nir = ele->e.w1->nGP[0];
   nis = ele->e.w1->nGP[1];
break;
case tri3: /* --> tri - element */  
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
      case tri3:   /* --> tri - element */              
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
      /*----------------------------------------calculate b_bar operator*/
      amzero(&b_bar_a);
      amzero(&int_b_bar_a);
      w1_b_barop(ele,b_bar,int_b_bar,boplin,F,numeps,nd,ip);      
      /*------------------------------------------ call material law ---*/
      amzero(&stress_a);
      amzero(&D_a);
      w1_call_matgeononl(mat,ele->e.w1->wtype,strain,stress,D,numeps);
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*---------------------- geometric part of stiffness matrix kg ---*/
        w1_kg(ele,kg,boplin,stress,fac,nd,numeps,ip);
      /*------------------ elastic+displacement stiffness matrix keu ---*/
        w1_keu(keu,b_bar,int_b_bar,D,fac,nd,numeps);
      /*------------------------------- stiffness matrix keug=keu+kg ---*/
       for (a=0; a<nd; a++)
       {  
         for (b=0; b<nd; b++)
         {
           estif[a][b] = kg[a][b] + keu[a][b];
         } 
       }
      /*--------------- nodal forces fi from integration of stresses ---*/
        if (force) w1_fint(ele,stress,int_b_bar,force,fac,nd,ip);                    
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
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
