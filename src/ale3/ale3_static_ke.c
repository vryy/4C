#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"
/*----------------------------------------------------------------------*
 | integration of linear stiffness ke for ALE element     mn 06/02      |
 *----------------------------------------------------------------------*/
void ale_static_ke(ELEMENT   *ele, 
                   ALE3_DATA  *data, 
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
double              vol;              /* element volume */

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
dstrc_enter("ale_static_ke");
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
ale_intg(ele,data,1);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
nir     = ele->e.ale3->nGP[0];
nis     = ele->e.ale3->nGP[1];
nit     = ele->e.ale3->nGP[2];
iel     = ele->numnp;
nd      = numdf * iel;
vol     = 0.0;
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
     e2   = data->xgps[ls];
     facs = data->wgts[ls];
    for (lt=0; lt<nit; lt++)
    {
     /*============================= gaussian point and weight at it ===*/
     e3   = data->xgpt[lt];
     fact = data->wgtt[lt];
     /*-------------------------- shape functions and their derivatives */
     ale_funct_deriv(funct,deriv,e1,e2,e3,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale_jaco (deriv,xjm,&det,ele,iel);
     /*------------------------------------ calculate element volume ---*/
     vol = vol + facr*facs*fact*det;
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale3->jacobi==1)
       fac = facr * facs *  fact * det;
     else
       fac = facr * facs *  fact;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
     if(nir == 1 && nis == 1 && nit ==1)
       ale3_hourglass(ele,estif,vol);
  }/*============================================== end of loop over lt */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of ale_static_ke */
/*----------------------------------------------------------------------*/
#endif
