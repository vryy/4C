/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale3_static_ke' which integrates the 
linear stiffness for the 3d ale element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale3.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  integration of linear stiffness ke for ALE element 

<pre>                                                              mn 06/02 
This routine integrates the linear stiffness for the 3d ale element


</pre>
\param *ele           ELEMENT   (i)   my element
\param *data          ALE3_DATA (i)   structure containing gaussian point and weight
\param *mat           MATERIAL  (i)   my material
\param *estif_global  ARRAY     (o)   the stifness matrix
\param init           INT       (i)   flag init == 1 : initialization
                                           init != 1 : integration

\warning There is nothing special to this routine
\return void                                               
\sa calling: ale3_intg(), ale3_funct_deriv(), ale3_jaco(), ale3_bop(),
             ale3_mat_linel(), ale_keku(), ale3_hourglass();
             called by: ale3()

*----------------------------------------------------------------------*/
void ale3_static_ke(ELEMENT   *ele, 
                   ALE3_DATA  *data, 
                   MATERIAL  *mat,
                   ARRAY     *estif_global, 
                   INT        init)
{
INT                 nir,nis,nit;      /* num GP in r/s/t direction */
INT                 lr, ls, lt;       /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
const INT           numdf =3;
const INT           numeps=6;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */
DOUBLE              vol;              /* element volume */

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
dstrc_enter("ale3_static_ke");
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
ale3_intg(ele,data);
/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
/*------------------------------------------- integration parameters ---*/
switch(ele->distyp)
{
    case hex8:
    case hex20:
        nir     = ele->e.ale3->nGP[0];
        nis     = ele->e.ale3->nGP[1];
        nit     = ele->e.ale3->nGP[2];
        break;
    case tet4:
    case tet10:
        nir = ele->e.ale3->nGP[0];
        nis = 1;
        nit = 1;
        break;
    default:
        dserror("unknown number of gaussian points in ale2_intg");
        break;
}
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
    for (lt=0; lt<nit; lt++)
    {
     /*============================= gaussian point and weight at it ===*/
        switch(ele->distyp)
        {
            case hex8:
            case hex20:
                e2   = data->xgps[ls];
                facs = data->wgts[ls];
                e3   = data->xgpt[lt];
                fact = data->wgtt[lt];
                break;
            case tet4:
            case tet10:
                e2   = data->xgps[lr];
                facs = ONE;
                e3   = data->xgpt[lr];
                fact = ONE;
                break;
            default:
                dserror("unknown number of gaussian points in ale2_intg");
                break;
        }
     /*-------------------------- shape functions and their derivatives */
     ale3_funct_deriv(funct,deriv,e1,e2,e3,ele->distyp,1);
     /*------------------------------------- compute jacobian matrix ---*/
     ale3_jaco (deriv,xjm,&det,ele,iel);
     /*------------------------------------ calculate element volume ---*/
     vol = vol + facr*facs*fact*det;
     /*----------------------------- use jacobian determinant or not ---*/
     if(ele->e.ale3->jacobi==1)
       fac = facr * facs *  fact * det;
     else
       fac = facr * facs *  fact;
     /*---------------------------------------- calculate operator B ---*/
     amzero(&bop_a);
     ale3_bop(bop,deriv,xjm,det,iel);
     /*------------------------------------------- call material law ---*/
     ale3_mat_linel(mat->m.stvenant,D);
     /*--------------------------------- elastic stiffness matrix ke ---*/
     ale_keku(estif,bop,D,fac,nd,numeps);
     /*---------------- hourglass stabalization  stiffness matrix ke ---*/
     if(ele->distyp==hex8 && nir == 1 && nis == 1 && nit ==1)
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
} /* end of ale3_static_ke */
/*----------------------------------------------------------------------*/
#endif
/*! @} (documentation module close)*/
