/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_static_ke' which integrates the 
linear stiffness for the 2d ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element 

<pre>                                                              mn 06/02 
This routine integrates the linear stiffness for the 2d ale element

</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE2_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void                                               
\sa calling: ale2_intg(), ale2_funct_deriv(), ale2_jaco(), ale2_bop(), 
             ale2_mat_linel(), ale_keku(), ale2_hourglass(); 
             called by: ale2()

*----------------------------------------------------------------------*/
void ale2_static_ke(ELEMENT   *ele, 
                   ALE2_DATA  *data, 
                   MATERIAL  *mat,
                   ARRAY     *estif_global, 
                   INT        init)
{
INT                 i,j,k;            /* some loopers */
INT                 nir,nis;      /* num GP in r/s/t direction */
INT                 lr, ls;       /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
const INT           numdf =2;
const INT           numeps=3;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE              xnu;              /* value of shell shifter */
DOUBLE              weight;

static ARRAY    D_a;      /* material tensor */     
static DOUBLE **D;         
static ARRAY    funct_a;  /* shape functions */    
static DOUBLE  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static DOUBLE **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static DOUBLE **xjm;         
static ARRAY    bop_a;    /* B-operator */   
static DOUBLE **bop;       
static DOUBLE **estif;    /* element stiffness matrix ke */

DOUBLE det;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("ale2_static_ke");
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
ale2_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
if(ele->distyp==quad4 || ele->distyp==quad8 || ele->distyp==quad9)
{
nir     = ele->e.ale2->nGP[0];
nis     = ele->e.ale2->nGP[1];
}
else if (ele->distyp==tri3 || ele->distyp==tri6)
{
    nir = ele->e.ale2->nGP[0];
    nis = 1;
}
else
{
   dserror("unknown number of gaussian points in ale2_intg");
}

iel     = ele->numnp;
nd      = numdf * iel;
/*================================================ integration loops ===*/
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgpr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
     /*============================= gaussian point and weight at it ===*/
     if(ele->distyp==quad4 || ele->distyp==quad8 || ele->distyp==quad9)
     {
         e2   = data->xgps[ls];
         facs = data->wgts[ls];
     }
     else if (ele->distyp==tri3 || ele->distyp==tri6)
     {
         e2   = data->xgps[lr];
         facs = ONE;
     }
     else
     {
       dserror("unknown number of gaussian points in ale2_intg");
     }
     /*---------------------- shape functions and their derivatives */
     ale2_funct_deriv(funct,deriv,e1,e2,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale2_jaco (deriv,xjm,&det,ele,iel);
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale2->jacobi==1)
       fac = facr * facs * det;
     else
       fac = facr * facs;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale2_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale2_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
     if(nir == 1 && nis == 1)
       ale2_hourglass(ele,estif);
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of ale2_static_ke */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
