/*!----------------------------------------------------------------------
\file
\brief contains the 2D material routines which calclate 
       derivatives for optimization

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0771 - 685-6575
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/
#ifdef D_WALL1
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"

/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_matd(ELEMENT     *ele,
                 MATERIAL     *mat, 
                 WALL_TYPE   wtype,
                 DOUBLE      **bop,
                 DOUBLE    *stress,
                 DOUBLE        **d);
/*----------------------------------------------------------------------*
 | constitutive matrix - linear elastic porous - 2D       al    9/01    |
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_matd_stvpor(MATERIAL  *mat,
                    DOUBLE *matdata,
                    WALL_TYPE wtype, 
                    DOUBLE **d);

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief integration routine for WALL1 element

<pre>                                                              al 06/02
This routine performs integration of an 3D-hex-element.

</pre>
\param           *ele ELEMENT  (i)   element data
\param          *data w1_DATA  (i)   hex element data
\param           *mat MATERIAL (i)   material data
\param        *retval DOUBLE   (o)   return value
\param         *init  INT      (i)   flag for initialization (alloc mem...)

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: w1_oint()

*----------------------------------------------------------------------*/
void w1_oint(
             ELEMENT   *ele, 
             W1_DATA   *data, 
             MATERIAL  *mat,
             DOUBLE    *retval,  /* return value */
             INT        init     /* ==2 calc.strain energy */
             )
{

/*----------------------------------------------------------------------*/
INT                 i;                                /* some loopers */
INT                 nir,nis;          /* num GP in r/s/t direction */
INT                 lr, ls;           /* loopers over GP */
INT                 iel;              /* numnp to this element */
INT                 nd;
INT                 ip;
INT                 istore, newval;
const INT           numdf  = 2;
const INT           numeps = 3;
DOUBLE              fac;              /* integration factors            */
DOUBLE              stifac;           /* thickness                      */
DOUBLE              e1,e2;            /* GP-coords                      */
DOUBLE              facr,facs;        /* weights at GP                  */
/*----------------------------------------------------------------------*/
static ARRAY    D_a;         /* material tensor */     
static DOUBLE **D;         
static ARRAY   DD_a;         /* material tensor */     
static DOUBLE **DD;         
static ARRAY    funct_a;     /* shape functions */    
static DOUBLE  *funct;     
static ARRAY    deriv_a;     /* derivatives of shape functions */   
static DOUBLE **deriv;     
static ARRAY    xjm_a;       /* jacobian matrix */     
static DOUBLE **xjm;         
static ARRAY    bop_a;       /* B-operator */   
static DOUBLE **bop;       

DOUBLE det;             /* Jacobi-det*/
DOUBLE disd[5];
DOUBLE strain[4];
DOUBLE  F[4];                 /* stress */
DOUBLE DF[4];                 /* stress */
/*----------------------------------------------------------------------*/
DOUBLE  ste;  /* strain energy               */
DOUBLE dste;  /* derivative of strain energy */
DOUBLE  stm;  /* mass                        */
DOUBLE  stv;  /* volume                      */

DOUBLE dens;  /* density                     */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_oint");
#endif
/*----------------------------------------------------------------------*/
istore = 0;
newval = 0;
/*------------------------------------------------- some working arrays */
if (init==1)
{
  funct     = amdef("funct"  ,&funct_a ,MAXNOD_WALL1,1 ,"DV");       
  deriv     = amdef("deriv"  ,&deriv_a ,2,MAXNOD_WALL1 ,"DA");       
  D         = amdef("D"      ,&D_a     ,6,6            ,"DA");           
  DD        = amdef("DD"     ,&DD_a  ,6,6              ,"DA");           
  xjm       = amdef("xjm"    ,&xjm_a ,numdf,numdf      ,"DA");           
  bop       = amdef("bop"    ,&bop_a   ,numeps,(numdf*MAXNOD_WALL1),"DA");           
  goto end;
}
/*------------------------------------------- integration parameters ---*/
  w1intg(ele,data,1);
/*------------------------------------------- integration parameters ---*/
  nir     = ele->e.w1->nGP[0];
  nis     = ele->e.w1->nGP[1];
  iel     = ele->numnp;
  nd      = numdf * iel;
/*---------------------------------- setup individual element arrays ---*/
  ste  = 0.;
  stv  = 0.;
  dste = 0.;
/*----------------------------------------------------------------------*/
  stifac = ele->e.w1->thick;
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
      w1_jaco (deriv,xjm,&det,ele,iel);                      
      /*------------------------------------ integration factor  -------*/ 
      fac = facr * facs * det * stifac;
      /*------------------------------------------- calculate volume ---*/
      if (init==5||init==3)
      {
          stv += fac;
          continue;
      }
      amzero(&bop_a);
      /*--------------------------------------- calculate operator B ---*/
      w1_bop(bop,deriv,xjm,det,iel);
      /*--------------------------- compute displacement derivatives ---*/        
      w1_disd (ele,bop,NULL,NULL,ele->e.w1->wtype,disd) ;                  
      /*------------------------------- get actual strains -> strain ---*/
      w1_eps (disd,ele->e.w1->wtype,strain);
      /*------------------------------------------ call material law ---*/
      w1_call_mat(ele,mat,ele->e.w1->wtype,bop,NULL,NULL,ip, F, D,istore,newval);
      /*------------------------------------ calculate strain energy ---*/
      if (init==2)
      {
          ste += (F[0]*strain[0] + F[1]*strain[1] + F[2]*strain[2] +
                  F[3]*strain[3])*fac*0.5;
          continue;
        
      }
      /*----------- calculate derivatives material matrix and stress ---*/
      if(init==4)
      {
        for (i=0;i<4;i++) DF[i] = 0.;
        w1_call_matd(ele,mat,ele->e.w1->wtype,bop,DF, DD);
       
        dste -= (DF[0]*strain[0] + DF[1]*strain[1] + DF[2]*strain[2] +
                 DF[3]*strain[3]) *fac*0.5;
      }
   }/*============================================= end of loop over ls */ 
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
          if(ele->e.w1->elewa->matdata==NULL) w1_getdensity(mat,&dens);
          else
          {
            dens = ele->e.w1->elewa->matdata[0];
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
} /* end of w1_oint */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | select proper material law                               al 01/02    |
 *----------------------------------------------------------------------*/
void w1_call_matd(ELEMENT     *ele,
                 MATERIAL     *mat, 
                 WALL_TYPE   wtype,
                 DOUBLE      **bop,
                 DOUBLE    *stress,
                 DOUBLE        **d)
{
INT i;
INT j;
DOUBLE disd[4];
DOUBLE strain[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_call_matd");
#endif
/*------------------------------------------------ call material law ---*/
  switch(mat->mattyp)
  {
  case m_stvenpor:/*------------------------ porous linear elastic ---*/
    w1_matd_stvpor(mat, ele->e.w1->elewa->matdata, wtype, d);
    w1_disd(ele,bop,NULL,NULL,wtype,disd);
    w1_eps (disd,ele->e.w1->wtype,strain);
    for (i=0; i<4; i++) stress[i] = 0.0;
    for (i=0; i<4; i++) for (j=0; j<4; j++) stress[i] += d[i][j]*strain[j];
  break;
  default:
    dserror(" material law not implemented for optimization");
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of w1_call_matd */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | constitutive matrix - linear elastic porous - 2D       al    9/01    |
 | plane stress, plane strain, rotational symmetry                      |
 *----------------------------------------------------------------------*/
void w1_matd_stvpor(MATERIAL  *mat,
                    DOUBLE *matdata,
                    WALL_TYPE wtype, 
                    DOUBLE **d)
{
DOUBLE e1, e2, e3, a1, b1, c1;
DOUBLE dym, ym, pv, dn, rd, ex;
#ifdef DEBUG 
dstrc_enter("w1_matd_stvpor");
#endif
/*----------------------------------------------------------------------*/
  /* input value: dn = mat->m.stvenpor->density     ;*/
  /* current values */
  dn = matdata[0];


  ym = mat->m.stvenpor->youngs      ;
  pv = mat->m.stvenpor->possionratio;
  rd = mat->m.stvenpor->refdens     ;
  ex = mat->m.stvenpor->exponent    ;

  dym=pow((dn/rd),(ex-1.)); 
  dym=(ym*ex*dym)/rd;
/*----------------------------------------------------- plane stress ---*/
  switch(wtype)
  {
  case plane_stress:
    e1=dym/(1. - pv*pv);
    e2=pv*e1;
    e3=e1*(1. - pv)/2.;

    d[0][0]=e1;
    d[0][1]=e2;
    d[0][2]=0.;
    d[1][0]=e2;
    d[1][1]=e1;
    d[1][2]=0.;
    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=e3;
  break;
  default:
/*-------------------------------- plane strain, rotational symmetry ---*/
    c1=dym/(1.0+pv);
    b1=c1*pv/(1.0-2.0*pv);
    a1=b1+c1;

    d[0][0]=a1;
    d[0][1]=b1;
    d[0][2]=0.;
    d[0][3]=b1;
    
    d[1][0]=b1;
    d[1][1]=a1;
    d[1][2]=0.;
    d[1][3]=b1;

    d[2][0]=0.;
    d[2][1]=0.;
    d[2][2]=c1/2.;
    d[2][3]=0.;
    
    d[3][0]=b1;
    d[3][1]=b1;
    d[3][2]=0.;
    d[3][3]=a1;
  break;
  }
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_matd_stvpor */

/*----------------------------------------------------------------------*/
#endif
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
