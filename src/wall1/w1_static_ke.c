#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | integration of linear stiffness ke for wall1 element      al 9/01    |
 *----------------------------------------------------------------------*/
void w1static_ke(ELEMENT   *ele, 
                 W1_DATA   *data, 
                 MATERIAL  *mat,
                 ARRAY     *estif_global, 
                 double    *force,  /* global vector for internal forces (initialized!) */
                 int        init)
{
int                 i,j,k;            /* some loopers */
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

double              fac;
double              stifac;
double              e1,e2,e3;         /*GP-coords*/
double              facr,facs,fact;   /* weights at GP */
double              xnu;              /* value of shell shifter */
double              weight;

static ARRAY    D_a;      /* material tensor */     
static double **D;         
static ARRAY    funct_a;  /* shape functions */    
static double  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static double **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static double **xjm;         
static ARRAY    bop_a;    /* B-operator */   
static double **bop;       
static double **estif;    /* element stiffness matrix ke */
double F[4];
double fie[18];

double det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1static_ke");
#endif
/*------------------------------------------------- some working arrays */
istore = 0;

if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_WALL1,1 ,"DV");       
deriv     = amdef("deriv"  ,&deriv_a,2,MAXNOD_WALL1 ,"DA");       
D         = amdef("D"      ,&D_a   ,6,6             ,"DA");           
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf     ,"DA");           

bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_WALL1),"DA");           
goto end;
}
else if(init==2)
{
  istore = 1;
}
/*------------------------------------------- integration parameters ---*/
w1intg(ele,data,1);
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
      
      
      fac = facr * facs * det * stifac; 
      /*--------------------------------------- calculate operator B ---*/
      amzero(&bop_a);
      w1_bop(bop,deriv,xjm,det,iel);
      /*------------------------------------------ call material law ---*/
      if(lanz==0)
      {
        w1_call_mat(ele, mat,ele->e.w1->wtype,bop,xjm,ip, F, D,istore,newval);
      }
      else
      {
        w1_mat_rebar(ele,mat,bop,xjm,F,D,ip,lanz,istore); 
      }
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*-------------------------------- element stiffness matrix ke ---*/
        w1_keku(estif,bop,D,fac,nd,numeps);
      /*--------------- nodal forces fi from integration of stresses ---*/
        if (force) w1fi (F,fac,bop,nd,force);                    
      }
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */
}
/*------------------------------- loop concrete reinforcement steel ----*/
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
#endif
