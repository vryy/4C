#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_calc.h"
/*----------------------------------------------------------------------*
 | integration of linear stiffness ke for BRICK1 element     al 9/01    |
 *----------------------------------------------------------------------*/
void b1static_ke(ELEMENT   *ele, 
                    B1_DATA   *data, 
                    MATERIAL  *mat,
                    ARRAY     *estif_global, 
                    int        init)
{
int                 i,j,k;            /* some loopers */
int                 nir,nis,nit;      /* num GP in r/s/t direction */
int                 lr, ls, lt;       /* loopers over GP */
int                 iel;              /* numnp to this element */
int                 nd;
const int           numdf =3;
const int           numeps=6;

double              fac;
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

double det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("b1static_ke");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");       
deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");       
D         = amdef("D"      ,&D_a   ,6,6              ,"DA");           
xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");           

bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");           
goto end;
}
/*------------------------------------------- integration parameters ---*/
b1intg(ele,data,1);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.b1->nGP[0];
nis     = ele->e.b1->nGP[1];
nit     = ele->e.b1->nGP[2];
iel     = ele->numnp;
nd      = numdf * iel;
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgrr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
     e2   = data->xgss[ls];
     facs = data->wgts[ls];
    for (lt=0; lt<nit; lt++)
    {
     /*============================= gaussian point and weight at it ===*/
     e3   = data->xgtt[lt];
     fact = data->wgtt[lt];
     /*-------------------------- shape functions and their derivatives */
     b1_funct_deriv(funct,deriv,e1,e2,e3,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     b1_jaco (deriv,xjm,&det,ele,iel);
     fac = facr * facs *  fact * det;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     b1_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     b1_mat_linel(mat->m.lin_el,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     b1_keku(estif,bop,D,fac,nd,numeps);
  }/*============================================== end of loop over lt */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of b1static_ke */
/*----------------------------------------------------------------------*/
