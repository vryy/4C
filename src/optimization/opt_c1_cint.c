/*!----------------------------------------------------------------------
\file
\brief contains the routine 'c1oint' which calclate displacement
       derivatives for a 3D hex element

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/
#ifdef D_BRICK1
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../brick1/brick1.h"
#include "../brick1/brick1_prototypes.h"

/*! 
\addtogroup BRICK1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief integration routine for BRICK1 element

<pre>                                                              al 06/02
This routine performs integration of an 3D-hex-element.

</pre>
\param           *ele ELEMENT  (i)   element data
\param          *data C1_DATA  (i)   hex element data
\param           *mat MATERIAL (i)   material data
\param  *estif_global ARRAY    (o)   element stiffness matrix
\param         *force DOUBLE   (o)   vector for internal forces
\param         *init  INT      (i)   flag for initialization (alloc mem...)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: c1_oint()

*----------------------------------------------------------------------*/
void c1_oint(
             ELEMENT   *ele, 
             C1_DATA   *data, 
             MATERIAL  *mat,
             DOUBLE    *retval,  /* return value */
             INT        init     /* ==2 calc.strain energy */
             )
{
INT                 i,j;              /* some loopers */
INT                 nir,nis,nit;      /* num GP in r/s/t direction */
INT                 lr, ls, lt;       /* loopers over GP */
INT                 ip;
INT                 iel;              /* numnp to this element */
INT                 nd;
INT                 istore = 0;/* controls storing of new stresses to wa */
INT                 newval = 0;/* controls evaluation of new stresses    */
const INT           numdf =3;
const INT           numeps=6;

DOUBLE              fac;
DOUBLE              e1,e2,e3;         /*GP-coords*/
DOUBLE              facr,facs,fact;   /* weights at GP */

DOUBLE disd[9];
DOUBLE F[6];  /* element stress vector   (stress-resultants) */
DOUBLE DF[6]; /* derivative of element stress vector         */
DOUBLE strain[6];
DOUBLE xyze[60];
DOUBLE edis[60];  
DOUBLE  g[6][6]; /* transformation matrix s(glob)= g*s(loc)   */
DOUBLE gi[6][6]; /* inverse of g          s(loc) = gi*s(glob) */
/*-------------------------------------------  for eas elements only ---*/
INT    cc;

static ARRAY    D_a;      /* material tensor */     
static DOUBLE **D;         
static ARRAY    DD_a;      /* material tensor */     
static DOUBLE **DD;         
static ARRAY    funct_a;  /* shape functions */    
static DOUBLE  *funct;     
static ARRAY    deriv_a;  /* derivatives of shape functions */   
static DOUBLE **deriv;     
static ARRAY    xjm_a;    /* jacobian matrix */     
static DOUBLE **xjm;         
static ARRAY    bop_a;    /* B-operator */   
static DOUBLE **bop;       
static ARRAY    bnop_a;   /* BN-operator */   
static DOUBLE **bn;       

DOUBLE det;

INT    iform;             /* index for nonlinear formulation of element */
INT    calstr;            /* flag for stress calculation                */

DOUBLE  ste;  /* strain energy               */
DOUBLE dste;  /* derivative of strain energy */
DOUBLE  stm;  /* mass                        */
DOUBLE  stv;  /* volume                      */

DOUBLE dens;  /* density                     */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("c1_oint");
#endif
/*----------------------------------------------------------------------*/
istore = 0;
calstr = 0;
newval = 0;
/*------------------------------------------------- some working arrays */
if (init==1)
{
  funct     = amdef("funct"  ,&funct_a,MAXNOD_BRICK1,1 ,"DV");       
  deriv     = amdef("deriv"  ,&deriv_a,3,MAXNOD_BRICK1 ,"DA");       
  D         = amdef("D"      ,&D_a   ,6,6              ,"DA");           
  DD        = amdef("DD"     ,&DD_a  ,6,6              ,"DA");           
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");           
  bop       = amdef("bop"  ,&bop_a ,numeps,(numdf*MAXNOD_BRICK1),"DA");           
  bn        = amdef("bnop" ,&bnop_a,3     ,       MAXNOD_BRICK1 ,"DA");           
goto end;
}
/*------------------------------------------- integration parameters ---*/
c1intg(ele,data);
/*------------------------------------------- integration parameters ---*/
  nir     = ele->e.c1->nGP[0];
  nis     = ele->e.c1->nGP[1];
  nit     = ele->e.c1->nGP[2];
  iel     = ele->numnp;
  nd      = numdf * iel;
/*---------------------------------- setup individual element arrays ---*/
  if(iel==8)
  {
    cc=0;
    for (i=0;i<iel;i++) for (j=0;j<3;j++) xyze[cc++] = ele->node[i]->x[j];
    cc=0;
    for (i=0;i<iel;i++) for (j=0;j<3;j++) edis[cc++] = ele->node[i]->sol.a.da[0][j];
  }
  else 
  {/*iel==20*/
   cc=0;
   xyze[cc++] = ele->node[0]->x[0];
   xyze[cc++] = ele->node[0]->x[1];
   xyze[cc++] = ele->node[0]->x[2];
   xyze[cc++] = ele->node[1]->x[0];
   xyze[cc++] = ele->node[1]->x[1];
   xyze[cc++] = ele->node[1]->x[2];
   xyze[cc++] = ele->node[2]->x[0];
   xyze[cc++] = ele->node[2]->x[1];
   xyze[cc++] = ele->node[2]->x[2];
   xyze[cc++] = ele->node[3]->x[0];
   xyze[cc++] = ele->node[3]->x[1];
   xyze[cc++] = ele->node[3]->x[2];
   xyze[cc++] = ele->node[4]->x[0];
   xyze[cc++] = ele->node[4]->x[1];
   xyze[cc++] = ele->node[4]->x[2];
   xyze[cc++] = ele->node[5]->x[0];
   xyze[cc++] = ele->node[5]->x[1];
   xyze[cc++] = ele->node[5]->x[2];
   xyze[cc++] = ele->node[6]->x[0];
   xyze[cc++] = ele->node[6]->x[1];
   xyze[cc++] = ele->node[6]->x[2];
   xyze[cc++] = ele->node[7]->x[0];
   xyze[cc++] = ele->node[7]->x[1];
   xyze[cc++] = ele->node[7]->x[2];
   xyze[cc++] = ele->node[8]->x[0];
   xyze[cc++] = ele->node[8]->x[1];
   xyze[cc++] = ele->node[8]->x[2];
   xyze[cc++] = ele->node[9]->x[0];
   xyze[cc++] = ele->node[9]->x[1];
   xyze[cc++] = ele->node[9]->x[2];
   xyze[cc++] = ele->node[10]->x[0];
   xyze[cc++] = ele->node[10]->x[1];
   xyze[cc++] = ele->node[10]->x[2];
   xyze[cc++] = ele->node[11]->x[0];
   xyze[cc++] = ele->node[11]->x[1];
   xyze[cc++] = ele->node[11]->x[2];
   xyze[cc++] = ele->node[16]->x[0];
   xyze[cc++] = ele->node[16]->x[1];
   xyze[cc++] = ele->node[16]->x[2];
   xyze[cc++] = ele->node[17]->x[0];
   xyze[cc++] = ele->node[17]->x[1];
   xyze[cc++] = ele->node[17]->x[2];
   xyze[cc++] = ele->node[18]->x[0];
   xyze[cc++] = ele->node[18]->x[1];
   xyze[cc++] = ele->node[18]->x[2];
   xyze[cc++] = ele->node[19]->x[0];
   xyze[cc++] = ele->node[19]->x[1];
   xyze[cc++] = ele->node[19]->x[2];
   xyze[cc++] = ele->node[12]->x[0];
   xyze[cc++] = ele->node[12]->x[1];
   xyze[cc++] = ele->node[12]->x[2];
   xyze[cc++] = ele->node[13]->x[0];
   xyze[cc++] = ele->node[13]->x[1];
   xyze[cc++] = ele->node[13]->x[2];
   xyze[cc++] = ele->node[14]->x[0];
   xyze[cc++] = ele->node[14]->x[1];
   xyze[cc++] = ele->node[14]->x[2];
   xyze[cc++] = ele->node[15]->x[0];
   xyze[cc++] = ele->node[15]->x[1];
   xyze[cc++] = ele->node[15]->x[2];
   cc=0;
   edis[cc++] = ele->node[0]->sol.a.da[0][0];
   edis[cc++] = ele->node[0]->sol.a.da[0][1];
   edis[cc++] = ele->node[0]->sol.a.da[0][2];
   edis[cc++] = ele->node[1]->sol.a.da[0][0];
   edis[cc++] = ele->node[1]->sol.a.da[0][1];
   edis[cc++] = ele->node[1]->sol.a.da[0][2];
   edis[cc++] = ele->node[2]->sol.a.da[0][0];
   edis[cc++] = ele->node[2]->sol.a.da[0][1];
   edis[cc++] = ele->node[2]->sol.a.da[0][2];
   edis[cc++] = ele->node[3]->sol.a.da[0][0];
   edis[cc++] = ele->node[3]->sol.a.da[0][1];
   edis[cc++] = ele->node[3]->sol.a.da[0][2];
   edis[cc++] = ele->node[4]->sol.a.da[0][0];
   edis[cc++] = ele->node[4]->sol.a.da[0][1];
   edis[cc++] = ele->node[4]->sol.a.da[0][2];
   edis[cc++] = ele->node[5]->sol.a.da[0][0];
   edis[cc++] = ele->node[5]->sol.a.da[0][1];
   edis[cc++] = ele->node[5]->sol.a.da[0][2];
   edis[cc++] = ele->node[6]->sol.a.da[0][0];
   edis[cc++] = ele->node[6]->sol.a.da[0][1];
   edis[cc++] = ele->node[6]->sol.a.da[0][2];
   edis[cc++] = ele->node[7]->sol.a.da[0][0];
   edis[cc++] = ele->node[7]->sol.a.da[0][1];
   edis[cc++] = ele->node[7]->sol.a.da[0][2];
   edis[cc++] = ele->node[8]->sol.a.da[0][0];
   edis[cc++] = ele->node[8]->sol.a.da[0][1];
   edis[cc++] = ele->node[8]->sol.a.da[0][2];
   edis[cc++] = ele->node[9]->sol.a.da[0][0];
   edis[cc++] = ele->node[9]->sol.a.da[0][1];
   edis[cc++] = ele->node[9]->sol.a.da[0][2];
   edis[cc++] = ele->node[10]->sol.a.da[0][0];
   edis[cc++] = ele->node[10]->sol.a.da[0][1];
   edis[cc++] = ele->node[10]->sol.a.da[0][2];
   edis[cc++] = ele->node[11]->sol.a.da[0][0];
   edis[cc++] = ele->node[11]->sol.a.da[0][1];
   edis[cc++] = ele->node[11]->sol.a.da[0][2];
   edis[cc++] = ele->node[16]->sol.a.da[0][0];
   edis[cc++] = ele->node[16]->sol.a.da[0][1];
   edis[cc++] = ele->node[16]->sol.a.da[0][2];
   edis[cc++] = ele->node[17]->sol.a.da[0][0];
   edis[cc++] = ele->node[17]->sol.a.da[0][1];
   edis[cc++] = ele->node[17]->sol.a.da[0][2];
   edis[cc++] = ele->node[18]->sol.a.da[0][0];
   edis[cc++] = ele->node[18]->sol.a.da[0][1];
   edis[cc++] = ele->node[18]->sol.a.da[0][2];
   edis[cc++] = ele->node[19]->sol.a.da[0][0];
   edis[cc++] = ele->node[19]->sol.a.da[0][1];
   edis[cc++] = ele->node[19]->sol.a.da[0][2];
   edis[cc++] = ele->node[12]->sol.a.da[0][0];
   edis[cc++] = ele->node[12]->sol.a.da[0][1];
   edis[cc++] = ele->node[12]->sol.a.da[0][2];
   edis[cc++] = ele->node[13]->sol.a.da[0][0];
   edis[cc++] = ele->node[13]->sol.a.da[0][1];
   edis[cc++] = ele->node[13]->sol.a.da[0][2];
   edis[cc++] = ele->node[14]->sol.a.da[0][0];
   edis[cc++] = ele->node[14]->sol.a.da[0][1];
   edis[cc++] = ele->node[14]->sol.a.da[0][2];
   edis[cc++] = ele->node[15]->sol.a.da[0][0];
   edis[cc++] = ele->node[15]->sol.a.da[0][1];
   edis[cc++] = ele->node[15]->sol.a.da[0][2];
  }
/*-------------------------------------  type of element formulation ---*/
  iform   = ele->e.c1->form;/*=1:linear:=2 total lagrangian formulation */

  ste  = 0.;
  stv  = 0.;
  dste = 0.;
/*================================================ integration loops ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
  /*================================ gaussian point and weight at it ===*/
  e1   = data->xgrr[lr];
  facr = data->wgtr[lr];
  for (ls=0; ls<nis; ls++)
  {
    /*============================== gaussian point and weight at it ===*/
    e2   = data->xgss[ls];
    facs = data->wgts[ls];
    for (lt=0; lt<nit; lt++)
    {
      ip++;
      /*============================ gaussian point and weight at it ===*/
      e3   = data->xgtt[lt];
      fact = data->wgtt[lt];
      /*------------------------- shape functions and their derivatives */
      c1_funct_deriv(funct,deriv,e1,e2,e3,iel,1);
      /*------------------------------------ compute jacobian matrix ---*/
      c1_jaco (deriv,xjm,&det,xyze,iel);
      fac = facr * facs *  fact * det;
      /*------------------------------------------- calculate volume ---*/
      if (init==5||init==3)
      {
          stv += fac;
          continue;
      }
      amzero(&bop_a);
      /*-- local element coordinate system for anisotropic materials ---*/
      c1tram (xjm,g,gi);
      /*--------------------------------------- calculate operator B ---*/
      c1_bop(bop,bn,deriv,xjm,det,iel);
      /*--------------------------- compute displacement derivatives ---*/        
      c1_disd (bop,edis,disd,iel) ;                  
      /*---------------- include initial displacements to b-operator ---*/        
      if(iform==2 && mat->mattyp!=m_pl_mises_ls)
      {
        c1_bdis (bop,disd,iel) ;                  
      }  
      /*------------------------------- get actual strains -> strain ---*/
      c1_eps (disd,strain,iform);
      /*------------------------------------------ call material law ---*/
      c1_call_mat(ele, mat,ip,F,strain,D,disd,g,gi,istore,newval);
      /*------------------------------------ calculate strain energy ---*/
      if (init==2)
      {
          ste += (F[0]*strain[0] + F[1]*strain[1] + F[2]*strain[2] +
                  F[3]*strain[3] + F[4]*strain[4] + F[5]*strain[5])
                  *fac*0.5;
          continue;
        
      }
      /*----------- calculate derivatives material matrix and stress ---*/
      if(init==4)
      {
        for (i=0;i<6;i++) DF[i] = 0.;
        c1_call_matd(ele, mat,DF,strain,DD,g);
        dste -= (DF[0]*strain[0] + DF[1]*strain[1] + DF[2]*strain[2] +
                 DF[3]*strain[3] + DF[4]*strain[4] + DF[5]*strain[5])
                 *fac*0.5;
      }
  }/*============================================== end of loop over lt */
  }/*============================================== end of loop over ls */
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
  switch(init)
  {
  case 5:
          (*retval) = stv;
  break;
  case 2:
          (*retval) = ste;
  break;
  case 3:
          if(ele->e.c1->elewa->matdata==NULL) c1_getdensity(mat,&dens);
          else
          {
            dens = ele->e.c1->elewa->matdata[0];
          }
          stm       = stv*dens;
          (*retval) = stm;
  break;
  case 4:
          (*retval) = dste;
  break;
  }

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of c1_oint */
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
