/*!----------------------------------------------------------------------
\file
\brief contains the routine 

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
#include "../wallge/wallge.h"
#include "../wallge/wallge_prototypes.h"
#include "s2ml.h"
#include "s2ml_prototypes.h"

/*! 
\addtogroup MLSTRUCT 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  routine for calculation of stiffnesses and internal forces
of a wall-submesh element

\param  *actsmmat         MATERIAL   (I)    actual submesh material
\param  *actmaele         ELEMENT    (I)    actual macro element
\param  *estif_ma_ma      ARRAY      (O)    submesh-element stiffness macro-macro
\param  *estif_ma_mi      ARRAY      (O)    submesh-element stiffness macro-micro
\param  *estif_mi_ma      ARRAY      (O)    submesh-element stiffness micro-macro
\param  *estif_mi_mi      ARRAY      (O)    submesh-element stiffness micro-micro
\param  *intforce_ma      ARRAY      (O)    submesh-element internal force macro
\param  *intforce_i       ARRAY      (O)    submesh-element internal force micro
\param   istore           INT        (I)    update?
\param   init             INT        (I)    allocate arrays? 

\return void                                               

*----------------------------------------------------------------------*/
void s2ml_stiff_wall(MATERIAL  *actsmmat,    /* actual submesh material*/
                     ELEMENT   *actmaele,    /* actual macro element   */
                     ELEMENT   *actsmele,    /* actual submesh element */
                     ARRAY     *estif_ma_ma, /* element stiffness macro-macro  */
                     ARRAY     *estif_ma_mi, /* element stiffness macro-micro  */
                     ARRAY     *estif_mi_ma, /* element stiffness micro-macro  */
                     ARRAY     *estif_mi_mi, /* element stiffness micro-micro  */
                     ARRAY     *intforce_ma, /* element internal force macro   */
                     ARRAY     *intforce_mi, /* element internal force micro   */
                     INT        istore,      /* update?                */
                     INT        init)        /* allocate arrays?       */
{
const INT           numdf  = 2;
const INT           numeps = 3;
INT         lr,ls,i;          /* loopers */
INT         ielmi;            /* numnp to this sm-element */
INT         ndmi;             /* numdof for this sm-element */
INT         ndma;             /* numdof for the actual macroelement */
INT         ip;               /* GP-number */
INT         nir,nis;          /* num GP in r/s direction */
DOUBLE      e1,e2;            /* GP r and s coordinates  */
DOUBLE      facr,facs;        /* weights at GP  */
DOUBLE      detmi;            /* det Jacobian of submesh  */
DOUBLE      facmi;            /* submesh integration factor */
DOUBLE      nue;              /* poisson ratio of actual sm-element */

static ARRAY    functmi_a;     /* submesh shape functions (f prime) */    
static DOUBLE  *functmi;     
static ARRAY    derivmi_a;     /* deriv. of sm-shape functions with respt r,s */   
static DOUBLE **derivmi;     
static ARRAY    functq4_a;     /* submesh quad4 shape functions (f prime) */    
static DOUBLE  *functq4;     
static ARRAY    xjmmi_a;       /* submesh jacobian matrix */     
static DOUBLE **xjmmi;         
static ARRAY    bopmi_a;       /* submesh B-operator */     
static DOUBLE **bopmi;         
static ARRAY    bopma_a;       /* macroelement B-operator */     
static DOUBLE **bopma;         
static ARRAY    D_a;           /* material tensor */     
static DOUBLE **D;         

static DOUBLE **stiff_ma_ma;      /* element stiffness macro-macro  */
static DOUBLE **stiff_ma_mi;      /* element stiffness macro-micro  */
static DOUBLE **stiff_mi_ma;      /* element stiffness micro-macro  */
static DOUBLE **stiff_mi_mi;      /* element stiffness micro-micro  */
static DOUBLE  *fint_ma;          /* element internal force macro   */
static DOUBLE  *fint_mi;          /* element internal force micro   */

W1_DATA      actsmdata;
DOUBLE       strainma[4];       /* macrostrains */
DOUBLE       strainmi[4];       /* small scale strains */
DOUBLE       strain_tot[4];     /* total strain = macrostrain + small scale strain */
DOUBLE       stress[4];         /* stress = stress(total strain) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s2ml_stiff_wall");
#endif
/*----------------------------------------------------------------------*/
if (init==1)
{
functmi     = amdef("functmi"  ,&functmi_a ,MAXNOD_WALL1,1 ,"DV");       
derivmi     = amdef("derivmi"  ,&derivmi_a ,2,MAXNOD_WALL1 ,"DA");       
functq4     = amdef("functq4"  ,&functq4_a ,4,1 ,"DV");       
xjmmi       = amdef("xjmmi"    ,&xjmmi_a   ,numdf,numdf    ,"DA");           
bopmi       = amdef("bopmi"    ,&bopmi_a   ,numeps,(numdf*MAXNOD_WALL1),"DA");           
bopma       = amdef("bopma"    ,&bopma_a   ,numeps,(numdf*MAXNOD_WALL1),"DA");           
D           = amdef("D"        ,&D_a       ,4,4            ,"DA");           

goto end;
}
/*------------------------------------------- integration parameters ---*/
w1intg(actsmele,&actsmdata,1);
nir     = actsmele->e.w1->nGP[0];
nis     = actsmele->e.w1->nGP[1];
ielmi   = actsmele->numnp;
ndmi    = numdf * ielmi;
ndma    = actmaele->numnp * numdf;
/*---- stifness and internal forces have to be reinitialized to zero ---*/

amzero(estif_ma_ma);
amzero(estif_ma_mi);
amzero(estif_mi_ma);
amzero(estif_mi_mi);
amzero(intforce_ma);
amzero(intforce_mi);
/*------------------------------------------------------- short cuts ---*/
stiff_ma_ma     = estif_ma_ma->a.da;
stiff_ma_mi     = estif_ma_mi->a.da;
stiff_mi_ma     = estif_mi_ma->a.da;
stiff_mi_mi     = estif_mi_mi->a.da;
fint_ma         = intforce_ma->a.dv;
fint_mi         = intforce_mi->a.dv;
/*----------------------------------------------------------------------*/
nue     = actsmmat->m.damage->possionratio;
/*================================= submesh element integration loop ===*/
ip = -1;
for (lr=0; lr<nir; lr++)
{
   /*======================= submesh gaussian point and weight at it ===*/
   e1   = (&actsmdata)->xgrr[lr];
   facr = (&actsmdata)->wgtr[lr];
   for (ls=0; ls<nis; ls++)
   {
      ip++;
      /*==================== submesh gaussian point and weight at it ===*/
      e2   = (&actsmdata)->xgss[ls];
      facs = (&actsmdata)->wgts[ls];
      /*-------------- submesh shape functions and their derivatives ---*/
      w1_funct_deriv(functmi,derivmi,e1,e2,actsmele->distyp,1);
      /*--- submesh quad4 shape functions for interpolation of bopma ---*/
      w1_funct_deriv(functq4,NULL,e1,e2,quad4,0);
      /*------------------------------------ submesh jacobian matrix ---*/       
      w1_jaco(derivmi,xjmmi,&detmi,actsmele,ielmi);                         
      /*--------------------------------- submesh integration factor ---*/ 
      facmi = facr * facs * detmi * actsmele->e.w1->thick;
      /*----------------------------------------- submesh B-operator ---*/
      amzero(&bopmi_a);
      w1_bop(bopmi,derivmi,xjmmi,detmi,ielmi);
      /*-- B-operator of macroelement and macrostrains at submesh GP ---*/
      amzero(&bopma_a);
      s2ml_bopstrainma(bopma,strainma,functq4,actsmele);
      /*---------------- submesh (small scale) strains at submesh GP ---*/
      s2ml_strainmi(bopmi,strainmi,actsmele,actmaele,nue);
      /*--------------------------------- total strain at submesh GP ---*/
      for (i=0; i<4; i++)
      {
        strain_tot[i] = strainma[i] + strainmi[i];
      }
      /*-------------------- call material law -> stress and tangent ---*/
      s2ml_callmat(actmaele,actsmele,actsmmat,strain_tot,ip,stress,D,istore);
      /*----------------------- element stiffness matrix macro macro ---*/
      w1_keku(stiff_ma_ma,bopma,D,facmi,ndma,numeps);
      /*----------------------- element stiffness matrix macro micro ---*/
      s2ml_cal_stiff(stiff_ma_mi,bopma,bopmi,D,facmi,ndma,ndmi,numeps);
      /*----------------------- element stiffness matrix micro macro ---*/
      s2ml_cal_stiff(stiff_mi_ma,bopmi,bopma,D,facmi,ndmi,ndma,numeps);
      /*----------------------- element stiffness matrix micro micro ---*/
      w1_keku(stiff_mi_mi,bopmi,D,facmi,ndmi,numeps);
      /*------------------------------ element internal forces macro ---*/
      wge_fintd(stress,facmi,bopma,ndma,numeps,fint_ma);
      /*------------------------------ element internal forces micro ---*/
      wge_fintd(stress,facmi,bopmi,ndmi,numeps,fint_mi);
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of s2ml_stiff_wall */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
