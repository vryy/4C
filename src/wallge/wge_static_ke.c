/*!-----------------------------------------------------------------------
\file
\brief contains the routine 'wge_static_ke' which calculates the 
stiffness matrix and internal forces for a gradient enhanced wall element

*-----------------------------------------------------------------------*/
#ifdef D_WALLGE
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"

/*! 
\addtogroup WALLGE
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates stiffness matrix and internal forcees  

<pre>                                                              mn 05/03
This routine calculates stiffness matrix for small strains formulation.

</pre>
\param **s       DOUBLE    (o)  blablabla 
\param   dl      DOUBLE    (i)  blablabal

\warning There is nothing special to this routine
\return void                                               
\sa calling:  ---; 
    caled by: wallge();

*----------------------------------------------------------------------*/
 /*----------------------------------------------------------------------*
 |          quad8                      quad4                           |
 |                 |                            |                        |
 |         3-------6--------2          3-----------------2               |
 |         |       |        |          |        |        |               |
 |         |  x2   x1   x0  |          |   x2   |    x0  |               |
 |         |       |        |          |        |        |               |
 |         |       |        |          |        |        |               |
 |    r <--7--x5---x4---x3--5--   r<---|--------|--------|--             |
 |         |       |        |          |        |        |               |
 |         |       |        |          |   x3   |    x1  |               |
 |         |  x8   x7   x6  |          |        |        |               |
 |         |       |        |          |        |        |               |
 |         0-------4--------1          0-----------------1               |
 |                 |                            |                        |
 |                 v                            v                        |
 |                 s                            s                        |
 *----------------------------------------------------------------------*/
void wgestatic_ke(ELEMENT       *ele, 
                  WALLGE_DATA   *data, 
                  MATERIAL      *mat,
                  ARRAY         *estif_global,
                  ARRAY         *emass_global,
                  DOUBLE        *force,  /* global int forces (0-initialized in the corrector, not needed for predictor) */
                  INT            init)
{
INT             lr,ls,i,j;       /* looper */
INT             ip;
INT             nir,nis;     /* num GP in r/s direction */
INT             ield;        /* num nodal points for displacements */
INT             iele;        /* num nodal points for nonlocal equiv.strains */
INT             istore = 0;  /* controls storing of new stresses to wa */
INT             newval = 0;  /* controls evaluation of new stresses    */
INT             neps   = 3;  /* number of strain components    */
INT             numdfd;      /* number of displacement DOF of element  */
INT             numdf;       /* total DOF of element  */

DOUBLE          e1,e2;       /* GP-coords                      */
DOUBLE          facr,facs;   /* weights at GP                  */
DOUBLE          detd;        /* Jacobi-det for displacement */
DOUBLE          dete;        /* Jacobi-det for nonlocal strains */
DOUBLE          thick;       /* thickness perpend. to wall plane */
DOUBLE          facd,face;   /* integration factor */
DOUBLE          eps_vl;      /* local equivalent strain at gauss point */
DOUBLE          eps_vnl;     /* nonlocal equivalent strain at gauss point */

static ARRAY    functd_a;     /* shape functions for displacements */    
static DOUBLE  *functd;     
static ARRAY    derivd_a;     /* derivatives of displacement shape functions */   
static DOUBLE **derivd;     
static ARRAY    functe_a;     /* shape functions for nonlocal equiv. strains */    
static DOUBLE  *functe;     
static ARRAY    derive_a;     /* derivatives of nonl. equiv.strain shape funct*/   
static DOUBLE **derive;     
static ARRAY    xjmd_a;       /* jacobian matrix for displacement */     
static DOUBLE **xjmd;
static ARRAY    xjme_a;       /* jacobian matrix for nonlocal strains */     
static DOUBLE **xjme;
static ARRAY    bopd_a;       /* B-operator for displacements */   
static DOUBLE **bopd;       
static ARRAY    bope_a;       /* B-operator for nonlocal equiv. strains */   
static DOUBLE **bope;       
static ARRAY    D_a;          /* 1.material tangente-tensor */     
static DOUBLE **D;         
static ARRAY    F_a;          /* 2.material tangente-tensor */     
static DOUBLE  *F;         
static ARRAY    E_a;          /* 3.material tangente-tensor */     
static DOUBLE  *E;         
static ARRAY    Kdd_a;        /* stiffness displacements - displacements */     
static DOUBLE **Kdd;         
static ARRAY    Kde_a;        /* stiffness displacements - equiv.strain */     
static DOUBLE **Kde;         
static ARRAY    Ked_a;        /* stiffness equiv.strain - displacements */     
static DOUBLE **Ked;         
static ARRAY    Kee_a;        /* stiffness equiv.strain - equiv.strain */     
static DOUBLE **Kee;         
       
static DOUBLE **estif;        /* element stiffness matrix ke */

DOUBLE sig[4];                /* stress */
DOUBLE grad_eps_vnl[2];       /* gradient of nonlocal equivalent strain */
DOUBLE fintd[18];             /* internal forces (displacement dof's)*/
DOUBLE finte[4];              /* internal forces (nonl.equiv.strain dof's) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("wgestatic_ke");
#endif
/*------------------------------------------------- some working arrays */
if (init==1)
{
 functd     = amdef("functd"  ,&functd_a ,MAXNOD_WALL1,1 ,"DV");       
 derivd     = amdef("derivd"  ,&derivd_a ,2,MAXNOD_WALL1 ,"DA");       
 functe     = amdef("functe"  ,&functe_a ,4,1 ,"DV");       
 derive     = amdef("derive"  ,&derive_a ,2,4 ,"DA");       
 xjmd       = amdef("xjmd"     ,&xjmd_a  ,2,2 ,"DA");           
 xjme       = amdef("xjme"     ,&xjme_a  ,2,2 ,"DA");           
 bopd       = amdef("bopd"    ,&bopd_a   ,3,(2*MAXNOD_WALL1),"DA");           
 bope       = amdef("bope"    ,&bope_a   ,2,4,"DA");           
 D          = amdef("D"       ,&D_a      ,4,4,"DA");           
 E          = amdef("E"       ,&E_a      ,4,1,"DV");           
 F          = amdef("F"       ,&F_a      ,4,1,"DV");           
 Kdd       = amdef("Kdd"      ,&Kdd_a    ,(2*MAXNOD_WALL1),(2*MAXNOD_WALL1),"DA");           
 Kde       = amdef("Kde"      ,&Kde_a    ,(2*MAXNOD_WALL1),4,"DA");           
 Kee       = amdef("Kee"      ,&Kee_a    ,4,4,"DA");           
 Ked       = amdef("Ked"      ,&Ked_a    ,4,(2*MAXNOD_WALL1),"DA");           
goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init=-1)                                        */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
  amdel(&functd_a);
  amdel(&derivd_a);
  amdel(&functe_a);
  amdel(&derive_a);
  amdel(&xjmd_a);
  amdel(&xjme_a);
  amdel(&bopd_a);
  amdel(&bope_a);
  amdel(&D_a);
  amdel(&E_a);
  amdel(&F_a);
  amdel(&Kdd_a);
  amdel(&Kde_a);
  amdel(&Kee_a);
  amdel(&Ked_a);
  goto end;  
}
/*--------------------------------------------- action = update_istep ---*/
else if(init==2)
{
  istore = 1;
}
/*-----------------------------------------------------------------------*/
/*--------------------------------------------- integration parameters---*/
wgeintg(ele,data,1);

/*-------------------------------------------- reinitalization to zero---*/
amzero(estif_global);
estif  = estif_global->a.da;
amzero(&Kdd_a);
amzero(&Kde_a);
amzero(&Kee_a);
amzero(&Ked_a);
for (i=0; i<18; i++) fintd[i] = 0.0;
for (i=0; i<4; i++)  finte[i] = 0.0;
/*--------------------------------------------- integration parameters---*/
nir     = ele->e.wallge->nGP[0];
nis     = ele->e.wallge->nGP[1];
ield    = ele->numnp;
numdfd  = ield*2;
iele    = 4;
numdf   = numdfd + iele;
thick   = ele->e.wallge->thick;
/*================================================= integration loops ===*/
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
      w1_funct_deriv(functd,derivd,e1,e2,ele->distyp,1);
      w1_funct_deriv(functe,derive,e1,e2,quad4,1);
      /*------------------------------------ compute jacobian matrix ---*/       
      w1_jaco(derivd,xjmd,&detd,ele,ield);                         
      w1_jaco(derive,xjme,&dete,ele,iele);                         
      /*------------------------------------ integration factor  -------*/ 
      facd = facr * facs * detd * thick;
      face = facr * facs * dete * thick;
      /*--------------------------------------- calculate operator B ---*/
      amzero(&bopd_a);
      w1_bop(bopd,derivd,xjmd,detd,ield);
      amzero(&bope_a);
      /* Verwendung von xjme und dete -> Konvergenzprobleme bei verzerrten Elementen
      wge_bope(bope,derive,xjme,dete);
      */
      wge_bope(bope,derive,xjmd,detd);
      /*------------------------------------------ call material law ---*/
      wge_call_mat(ele,mat,bopd,functe,bope,ip,sig,&eps_vl,&eps_vnl,grad_eps_vnl,D,E,F,istore,newval);
      /*----------------------------------------------------------------*/
      if(istore==0)
      {
      /*---------------------------- element stiffness matrices kele ---*/
         w1_keku(Kdd,bopd,D,facd,numdfd,neps);
         wge_stiff_de(Kde,bopd,E,functe,facd,numdfd,iele,neps);
         wge_stiff_ed(Ked,bopd,F,functe,facd,numdfd,iele,neps);
         /*
         wge_stiff_ee(Kee,bope,mat->m.damage_ge->crad,functe,face,iele);
         */
         wge_stiff_ee(Kee,bope,mat->m.damage_ge->crad,functe,facd,iele);
      /*-------------------------------------------- internal forces ---*/        
       if (force)
       {
         wge_fintd(sig,facd,bopd,numdfd,neps,fintd);
         /*
         wge_finte(eps_vl,eps_vnl,grad_eps_vnl,mat->m.damage_ge->crad,face,bope,functe,iele,finte);
         */
         wge_finte(eps_vl,eps_vnl,grad_eps_vnl,mat->m.damage_ge->crad,facd,bope,functe,iele,finte);
       } 
      } 
   }/*============================================= end of loop over ls */ 
}/*================================================ end of loop over lr */

/*----------------------------------------------------------------------*/
if(istore==0)
{
  wge_permstiff(Kdd,Kde,Ked,Kee,iele,ield,numdf,estif);
  if (force)
  {
    wge_permforce(fintd,finte,iele,ield,numdf,force);
  }
}
/*----------------------------------------------------------------------*/
# if 0
for (i=0;i<20;i++)
{
  for (j=0;j<20;j++)
  {
    printf ("%9.4E",estif[i][j]);
  }
  printf ("\n");
}
# endif

/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of wgestatic_ke */



/*----------------------------------------------------------------------*/
#endif /*D_WALLGE*/
/*! @} (documentation module close)*/
