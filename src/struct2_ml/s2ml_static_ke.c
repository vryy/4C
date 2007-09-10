/*!-----------------------------------------------------------------------
\file
\brief contains the routine 's2static_ke_ml' which calculates the
stiffness matrix and internal forces for a structural multiscale analysis

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*-----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_MLSTRUCT

#include "../headers/standardtypes.h"
#include "../wall1/wall1.h"
#include "../wall1/wall1_prototypes.h"
#include "s2ml_prototypes.h"
#include "s2ml.h"
#include "../solver/solver.h"

/*!
\addtogroup MLSTRUCT
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates stiffness matrix and internal forcees

<pre>                                                            ah 06/04
This routine calculates stiffness matrix for small strains formulation.
</pre>

\param  *actmaele       ELEMENT    (I)  actual macroelement
\param  *data           W1_DATA    (I)  Integration data
\param  *mat            MATERIAL   (I)  material of the actual macroelement
\param  *estif_global   ARRAY      (O)  condensed stiffness matrix for macro-element-DOF's
\param  *emass_global   ARRAY      (O)  condensed mass matrix for macro-element-DOF's if dynamic (not jet!)
\param  *force          DOUBLE     (O)  condensed internal forcees for macro-element-DOF's
\param   init           INT        (I)  init-phase 0:calc_nlnstiff 2: update 3: calc_stress

\warning There is nothing special to this routine
\return void

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/* statvar->praedictor==1 : global predictor */
/* statvar->praedictor==2 : first global correktor iteration */
/* statvar->praedictor==3 : all the following global correktor iterations */
/*----------------------------------------------------------------------*
 |                                                       ah 03/04       |
 | vector of numfld submesh-FIELDs, defined in global_control.c         |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
/*- ist eine globale Variable und dafuer global in global_control.c
    definiert, deshalb hier extern */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | global variable *sm_solv, vector of lenght numfld of structures SOLVAR  |
 | defined in input_submesh.c                                          |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *sm_solv;
/*----------------------------------------------------------------------*/
void s2ml_static_ke(ELEMENT       *actmaele,
                    W1_DATA       *data,
                    MATERIAL      *mat,
                    ARRAY         *estif_global,
                    ARRAY         *emass_global,
                    DOUBLE        *force,  /* global int forces (0-initialized in the corrector, not needed for predictor) */
                    INT            init)/* 1: init-phase 0:calc_nlnstiff 2: update 3: calc_stress */
{
INT                 i,j,lr,ls,a,b,smnodei;    /* some loopers */
INT                 nir,nis;          /* num GP in r/s direction */
INT                 ielma;            /* numnp to this Makro-element */
INT                 ndma;             /* numDOF to this Makro-element */
INT                 multiscale = 0;   /* is multiscale analysis to be done for this element? */
INT                 ip;

INT                 istore = 0;       /* controls storing of new stresses to wa */
INT                 newval = 0;       /* controls evaluation of new stresses    */
const INT           numdf  = 2;
const INT           numeps = 4;

DOUBLE              e1,e2;            /* GP-coords                      */
DOUBLE              detma;            /* Macro-Jacobi-detGP-coords    */
DOUBLE              nue;              /* poission ratio   */
/*DOUBLE              eps_equiv;        equivalent strain    */
DOUBLE              facr,facs;        /* weights at GP                  */
DOUBLE              fac;              /* integration factor            */

static ARRAY    functma_a;     /* Makro-shape functions */
static DOUBLE  *functma;
static ARRAY    derivma_a;     /* derivatives of Makro-shape functions */
static DOUBLE **derivma;
static ARRAY    xjmma_a;       /* Makro-jacobian matrix */
static DOUBLE **xjmma;
static ARRAY    bopma_a;       /*  Makro-B-operator */
static DOUBLE **bopma;
static ARRAY    D_a;           /* material tensor */
static DOUBLE **D;

static DOUBLE **estif;          /* (reduced) element stiffness matrix */

DOUBLE fint[8];               /* (reduced) internal force */
DOUBLE Stress[4];                   /* makro-stress */

/*----------------------------------------------------------------------*/

INT           numeq_sm;        /* number of submesh DOF's (without the dirichlet ones) */
INT           init_sol;        /* init-flag for solver */

static ARRAY         smintforcemi_a;    /* assembled  vector of submesh internal forces micro */
static DOUBLE       *smintforcemi;      /* pointer to smintforcemi_a.da.dv */

static ARRAY         smintforcema_a;    /* "assembled"  vector of submesh internal forces macro */
static DOUBLE       *smintforcema;      /* pointer to smintforcema_a.da.dv */

static ARRAY         smstiffmama_a;    /* "assembled"  vector of submesh stiffness macro */
static DOUBLE      **smstiffmama;      /* pointer to  */

static ARRAY         densesmstiffmimaT_a;   /* dense assebled sm-stiffness_mi_ma*/
static DOUBLE      **densesmstiffmimaT;
static DBCSR        *smstiffmima_csr; /* over submesh assebled stiffness_mi_ma*/
static DBCSR        *smstiffmami_csr; /* over submesh assebled stiffness_ma_mi*/

static ARRAY         workforce_a;    /* stiff_mi_ma*DELTAdbar (needed for calc. of d') */
static DOUBLE       *workforce;      /* pointer to */

SOLVAR        *actsmsolv;       /* the active submeshSOLVAR */
static INTRA  *actsmintra;      /* submesh Pseudo intra-communicator
(static: dann muss ich ihn nur einmal allociern und er wird tortxdem nicht vergessen) */

FIELD        *actsmfield;      /* pointer to active submeshfield */

DOUBLE        DELTA_F[8];      /* DELTA force in static condensation */
DOUBLE        DELTA_K[8][8];   /* DELTA stiffness in static condensation */

DOUBLE        DELTA_Dbar[8];   /* incremental (in 1st korrector)/residual(else) sol of macroDOF's at maelem.nodes */
NODE         *actsmnode;       /* actual submeshnode                    */
INT           ID;              /* ID of actual submeshnode */
INT           dof;             /* degree of freedom number of submesh-dof */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s2ml_static_ke");
#endif
/*----------------------------------------------------------------------*/
/* init phase        (init = 1)                                         */
/*----------------------------------------------------------------------*/
numeq_sm = sm_field[0].dis[0].numeq;
if (init==1)
{
  functma     = amdef("functma"  ,&functma_a ,MAXNOD_WALL1,1 ,"DV");
  derivma     = amdef("derivma"  ,&derivma_a ,2,MAXNOD_WALL1 ,"DA");
  xjmma       = amdef("xjmma"    ,&xjmma_a   ,numdf,numdf    ,"DA");
  bopma       = amdef("bopma"    ,&bopma_a   ,numeps,(numdf*MAXNOD_WALL1),"DA");
  D           = amdef("D"        ,&D_a       ,6,6             ,"DA");
/*----------------------------------------------------------------------*/
  smintforcemi  = amdef("smintforcemi",&smintforcemi_a,numeq_sm,1,"DV");
  smintforcema  = amdef("smintforcema",&smintforcema_a,8,1,"DV");
  smstiffmama   = amdef("smstiffmama", &smstiffmama_a, 8,8,"DA");
  smstiffmima_csr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
  smstiffmami_csr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
  /* Zeilen und Spalten-anzahl absichtlich vertauscht, s.u. -> Transposed */
  densesmstiffmimaT  = amdef("densesmstiffmimaT"  ,&densesmstiffmimaT_a ,8,numeq_sm,"DA");
  /* allocation of 10 rhs and sol vectors for static cond. of sm_dof's    */
  /* has alraedy been done in mask_submesh_matrices()                    */
  workforce  = amdef("workforce",&workforce_a,numeq_sm,1,"DV");
  /* because we are not parallel here, we have to allocate a pseudo-intracommunicator */
  actsmintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  goto end;
}
/*----------------------------------------------------------------------*/
/* uninit phase        (init = -1)                                      */
/*----------------------------------------------------------------------*/
else if (init==-1)
{
   amdel(&functma_a);
   amdel(&derivma_a);
   amdel(&xjmma_a);
   amdel(&bopma_a);
   amdel(&D_a);
}
/*----------------------------------------------------------------------*/
/* action = update_istep  (init = 2)                                    */
/*----------------------------------------------------------------------*/
else if(init==2)
{
  istore = 1;
}
/*----------------------------------------------------------------------*/
/*                  calculation phase (init = 0)                        */
/*----------------------------------------------------------------------*/

/*-------------- some of the fields have to be reinitialized to zero ---*/
amzero(estif_global);
estif     = estif_global->a.da;
for (i=0; i<8; i++) fint[i] = 0.0;

/*------------------------------------------- integration parameters ---*/
w1intg(actmaele,data,1);
nir     = actmaele->e.w1->nGP[0];
nis     = actmaele->e.w1->nGP[1];
ielma   = actmaele->numnp;
if(ielma>4)
{
   dserror("For multiscale analysis, the macroelements must be bilinear");
}
ndma    = numdf * ielma;

nue     = mat->m.damage->possionratio;


# if 0

/*  check if Makroelement has already been in Omega' in last load step  */
if(actmaele->e.w1->isinomegaprime == 0)
{
  /*-- if not: check initial.cond for multiscale -> element in Omega' ---*/
  multiscale = 0;
  for (lr=0; lr<nir; lr++)
  {
    e1   = data->xgrr[lr];
    for (ls=0; ls<nis; ls++)
    {
       e2   = data->xgss[ls];
       /*------------------  Makro shape functions and their derivatives */
       w1_funct_deriv(functma,derivma,e1,e2,actmaele->distyp,1);
       /*------------------------------ compute Makro-jacobian matrix ---*/
       w1_jaco (functma,derivma,xjmma,&detma,actmaele,ielma);
       /*--------------------------------- calculate Makro-operator B ---*/
       amzero(&bopma_a);
       w1_bop(bopma,derivma,xjmma,detma,ielma);
       /*-------------------------- calculate equivalent strain in GP ---*/
       s2ml_aequistrain(actmaele,actmaele->e.w1->wtype,bopma,nue,&eps_equiv);
       if(eps_equiv >= statvar->eps_equiv)  multiscale = 1;
    }
  }
}/*  end: if(actmaele->e.w1->isinomegaprime == 0)                   */
# endif

multiscale = 1;
/*----------------------------------------------------------------------*/
/*     Element is not in Omega'                                         */
/*----------------------------------------------------------------------*/

if(actmaele->e.w1->isinomegaprime == 0 && multiscale ==0)
{
/*     calculate stiffness K-bar-bar and internal force Fint bar        */

  /*============================================== integration loops ===*/
  ip = -1;
  for (lr=0; lr<nir; lr++)
  {
     /*============================= gaussian point and weight at it ===*/
     e1   = data->xgrr[lr];
     facr = data->wgtr[lr];
     for (ls=0; ls<nis; ls++)
     {
        ip++;
        /*========================== gaussian point and weight at it ===*/
        e2   = data->xgss[ls];
        facs = data->wgts[ls];
        /*----------------------- shape functions and their derivatives */
        w1_funct_deriv(functma,derivma,e1,e2,actmaele->distyp,1);
        /*---------------------------------- compute jacobian matrix ---*/
        w1_jaco (derivma,xjmma,&detma,actmaele,ielma);
        /*---------------------------------- integration factor  -------*/
        fac = facr * facs * detma * actmaele->e.w1->thick;
        /*------------------------------------- calculate operator B ---*/
        amzero(&bopma_a);
        w1_bop(bopma,derivma,xjmma,detma,ielma);
        /*---------------------------------------- call material law ---*/
        w1_call_mat(actmaele,mat,actmaele->e.w1->wtype,bopma,NULL,NULL,ip,Stress,D,istore,newval);
        /*--------------------------------------------------------------*/
        if(init==3)
        {
           for (i=0; i<4; i++)
           {
             actmaele->e.w1->stress_GP.a.d3[0][i][ip]= Stress[i];
           }
        }
        else if(istore==0)
        {
         /*----------------------------- element stiffness matrix ke ---*/
         w1_keku(estif,bopma,D,fac,ndma,numeps);
         /*------------ nodal forces fi from integration of stresses ---*/
         if (force)
         {
           w1fi (Stress,fac,bopma,ndma,force);
         }
       }
     }/*=========================================== end of loop over ls */
  }/*============================================== end of loop over lr */
}/*  end: if(actmaele->e.w1->isinomegaprime == 0 && multiscale ==0) */

/*----------------------------------------------------------------------*/
/*     Element is in Omega'                                             */
/*----------------------------------------------------------------------*/
else if(init!=3)/* if action is not calc_stress */
{

  actsmsolv  = &(sm_solv[0]);
  actsmfield = &(sm_field[0]);
/*-------------------------------------------------- !Zeilenanzahl! ----*/
  smstiffmima_csr->numeq       = numeq_sm;
  smstiffmima_csr->numeq_total = numeq_sm;
  smstiffmami_csr->numeq       = ndma;
  smstiffmami_csr->numeq_total = ndma;
  actsmintra->intra_fieldtyp = structure;
  actsmintra->intra_rank     = 0;
  actsmintra->intra_nprocs   = 1;
  /* set RHS and solution vectors for static condens. to zero           */
  for (i=0; i<actsmsolv->nrhs; i++) solserv_zero_vec(&(actsmsolv->rhs[i]));
  for (i=0; i<actsmsolv->nsol; i++) solserv_zero_vec(&(actsmsolv->sol[i]));
/*   allocation of some arrays if element is first time in Omega prime  */
  if(actmaele->e.w1->firstinomegaprime == 1)
  {
    s2ml_init(actmaele,actsmfield->dis[0].numnp,actsmfield->dis[0].numele);

  }
/*----------------------------------------------------------------------*/
/*   calculate Bbar,Nbar and strainbar at submesh-nodes and store there */
   s2ml_bopstraintonode(actsmfield,actmaele,nue,0);
/*----------------------------------------------------------------------*/
/*   update displacements d'-> store at working-array of macro-element  */
/*----------------------------------------------------------------------*/
  if((statvar->praedictor == 2 || statvar->praedictor == 3) && actmaele->e.w1->firstinomegaprime != 1)
  {
    for (i=0; i<ielma; i++)/* this is the incremental solution after the predictor  */
    { /*(i.e. at the beginning of the 1.Korrector), and the residual solution in the following correktor steps !*/
      DELTA_Dbar[2*i]   = actmaele->node[i]->sol_residual.a.da[0][0];
      DELTA_Dbar[2*i+1] = actmaele->node[i]->sol_residual.a.da[0][1];
    }
    /* calculate workforce:= - fint' - K-'bar *  (DELTA d-bar)          */
    for (i=0; i<numeq_sm; i++)
    {
      workforce[i] = 0.0;
      for(j=actmaele->e.w1->stiff_mi_ma_csr->ia.a.iv[i]; j<actmaele->e.w1->stiff_mi_ma_csr->ia.a.iv[i+1]; j++)
      {
        workforce[i] -= actmaele->e.w1->stiff_mi_ma_csr->a.a.dv[j] * DELTA_Dbar[actmaele->e.w1->stiff_mi_ma_csr->ja.a.iv[j]];
      }
      workforce[i] -= actmaele->e.w1->fint_mi[i];
    }
    /* workforce = -(fint'+K'bar*DELTAdbar)->DIST-vector:actsmsolv->rhs[9]*/
    assemble_vec(actsmintra,&(actsmsolv->sysarray_typ[0]),
               &(actsmsolv->sysarray[0]),&(actsmsolv->rhs[9]),workforce,1.0);
    /* solve for DELTA d'= K''^-1 * workforce(do not use an old LU-decomp of K'')  */
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->reuse = 0;
    init_sol = 0;
    solver_control(actsmfield,0,
                   actsmsolv,
                   actsmintra,
                 &(actsmsolv->sysarray_typ[0]),    /* der Typ ist ccf, genau wie bei sysarray */
                 &(actmaele->e.w1->stiff_mi_mi_ccf[0]),
                 &(actsmsolv->sol[9]),
                 &(actsmsolv->rhs[9]),
                   init_sol);
    /* update u' = u'+DELTAu' and store at macroelement                 */
    for (smnodei=0; smnodei<actsmfield->dis[0].numnp; smnodei++)
    {
      actsmnode = &(actsmfield->dis[0].node[smnodei]);
      ID = actsmnode->Id;
      for (j=0; j<numdf; j++)
      {
        dof = actsmnode->dof[j];
        if(dof >= numeq_sm) continue;
        /* incrementelle Verschiebung DELTA d'  */
        actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[j] = actsmsolv->sol[9].vec.a.dv[dof];
        /* Verschiebung d' = d' alt + DELTA d'  */
        actmaele->e.w1->sm_nodaldata[ID].displ_mi[j] += actsmsolv->sol[9].vec.a.dv[dof];
      }
    }
  }/* end: if((statvar->praedictor == 2 || statvar->praedictor == 3) && actmaele->e.w1->firstinomegaprime != 1) */
/*----------------------------------------------------------------------*/

/*--------------------------------- init the dist sparse matrix to zero */
  solserv_zero_mat(actsmintra,
                 &(actsmsolv->sysarray[0]),
                 &(actsmsolv->sysarray_typ[0]));
  amzero(&smintforcemi_a);
  amzero(&smintforcema_a);
  amzero(&smstiffmama_a);
/*----------------------------------------------------------------------*/
/*       loop submesh -> Fint',Fintbar, K'',K'bar, Kbar' Kbarbar        */
  s2ml_cal_smelm(actmaele,actsmfield,actsmsolv,actsmintra,
                 smintforcemi,smstiffmima_csr,smstiffmami_csr,smintforcema,smstiffmama,
                 istore);

/*----------------------------------------------------------------------*/
/* assebled sm-int forces smintforcemi -> DIST-vector: actsmsolv->rhs[0]*/
  assemble_vec(actsmintra,&(actsmsolv->sysarray_typ[0]),
             &(actsmsolv->sysarray[0]),&(actsmsolv->rhs[0]),smintforcemi,1.0);
/*----------------------------------------------------------------------*/
/*   DENSE CSR-Matrix smstiffmima_csr -> dense Matrix densesmstiffmima  */
/*   die Matrix wird mit Vertauschung Zeilen-Spalten dense gespeichert, */
/*   da ich auf die einzelnen SPALTEN als Vektoren zugreifen moechte,   */
/*   aber bei arrays name[i] die ZEILE des arrays name[i][j] ist.       */
/*                    -> T fuer transposed                              */
  for (i=0; i<numeq_sm; i++)
  {
    for(j=smstiffmima_csr->ia.a.iv[i]; j<smstiffmima_csr->ia.a.iv[i+1]; j++)
    {
/* eigentlich:densesmstiffmima[i][smstiffmima_csr->ja.a.iv[j]] = smstiffmima_csr->a.a.dv[j];*/
       densesmstiffmimaT[smstiffmima_csr->ja.a.iv[j]][i] = smstiffmima_csr->a.a.dv[j];
    }
  }
/*----------------------------------------------------------------------*/
/* assebled RHS's sm_stiff_mima(dense) -> DIST-vectors: actsmsolv->rhs[1..8] */
  for (i=0; i<8; i++)
  {
     assemble_vec(actsmintra,&(actsmsolv->sysarray_typ[0]),
                &(actsmsolv->sysarray[0]),&(actsmsolv->rhs[i+1]),densesmstiffmimaT[i],1.0);
  }
/*----------------------------------------------------------------------*/
/*     solve for K'' * sol[0] = Fint'    -> LU- decomposition of K''    */
  init_sol = 0;
  actsmsolv->sysarray[0].ccf->reuse = 0;
  solver_control(actsmfield, 0,
                 actsmsolv,
                 actsmintra,
               &(actsmsolv->sysarray_typ[0]),
               &(actsmsolv->sysarray[0]),
               &(actsmsolv->sol[0]),
               &(actsmsolv->rhs[0]),
                 init_sol);

/*----------------------------------------------------------------------*/
/*solve for K''*sol[i]=K'bar[i.Spalte] LU-dec.of K''not to be done again*/
  actsmsolv->sysarray[0].ccf->reuse = 1;
  for (i=0; i<8; i++)
  {
    solver_control(actsmfield, 0,
                   actsmsolv,
                   actsmintra,
                 &(actsmsolv->sysarray_typ[0]),
                 &(actsmsolv->sysarray[0]),
                 &(actsmsolv->sol[i+1]),
                 &(actsmsolv->rhs[i+1]),
                   init_sol);
  }
/*----------------------------------------------------------------------*/
/*          calculate DELTA_F:= K-bar' *  (K''^-1 * F'int)              */
  for (i=0; i<ndma; i++)
  {
    DELTA_F[i] = 0.0;
    for(j=smstiffmami_csr->ia.a.iv[i]; j<smstiffmami_csr->ia.a.iv[i+1]; j++)
    {
      DELTA_F[i] -= smstiffmami_csr->a.a.dv[j] *
                    actsmsolv->sol[0].vec.a.dv[smstiffmami_csr->ja.a.iv[j]];
    }
  }
/*----------------------------------------------------------------------*/
/*          calculate DELTA_K:= K-bar' *  (K''^-1 * K'bar)              */
  for (b=0; b<ndma; b++)
  {
    for (a=0; a<ndma; a++)
    {
      DELTA_K[a][b] = 0.0;
      for(j=smstiffmami_csr->ia.a.iv[a]; j<smstiffmami_csr->ia.a.iv[a+1]; j++)
      {
        DELTA_K[a][b] += smstiffmami_csr->a.a.dv[j] *
                         actsmsolv->sol[b+1].vec.a.dv[smstiffmami_csr->ja.a.iv[j]];
      }
    }
  }
/*----------------------------------------------------------------------*/
/*             static condensation of submesh-DOF's d' :                */
/*         Kred = K-macro-macro - DELTA K, Fred = F-macro - DELTA F     */
  for (a=0; a<ndma; a++)
  {
    force[a] = smintforcema[a] - DELTA_F[a];
    for (b=0; b<ndma; b++)
    {
      estif[a][b] = smstiffmama[a][b] - DELTA_K[a][b];
    }
  }
/*----------------------------------------------------------------------*/
/*              store stiffnesses k'' and k'bar and fint'               */
/*----------------------------------------------------------------------*/
  if(istore != 1)
  {
    /* store fint' */
    for (i=0; i<numeq_sm; i++)
    {
      actmaele->e.w1->fint_mi[i] = smintforcemi[i];
    }
    if(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ax.Typ == cca_XX)/* first time: allocate and store */
    {
      /* store k'' */
      am_alloc_copy(&(actsmsolv->sysarray->ccf->update),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->update));
      am_alloc_copy(&(actsmsolv->sysarray->ccf->Ap),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ap));
      am_alloc_copy(&(actsmsolv->sysarray->ccf->Ai),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ai));
      am_alloc_copy(&(actsmsolv->sysarray->ccf->Ax),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ax));
      /* store k'bar */
      am_alloc_copy(&(smstiffmima_csr->update),&(actmaele->e.w1->stiff_mi_ma_csr->update));
      am_alloc_copy(&(smstiffmima_csr->a),&(actmaele->e.w1->stiff_mi_ma_csr->a));
      am_alloc_copy(&(smstiffmima_csr->ja),&(actmaele->e.w1->stiff_mi_ma_csr->ja));
      am_alloc_copy(&(smstiffmima_csr->ia),&(actmaele->e.w1->stiff_mi_ma_csr->ia));
    }
    else  /* not the first time: do not allocate again but store */
    {
      /* store k'' */
      amcopy(&(actsmsolv->sysarray->ccf->update),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->update));
      amcopy(&(actsmsolv->sysarray->ccf->Ap),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ap));
      amcopy(&(actsmsolv->sysarray->ccf->Ai),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ai));
      amcopy(&(actsmsolv->sysarray->ccf->Ax),&(actmaele->e.w1->stiff_mi_mi_ccf->ccf->Ax));
      /* store k'bar */
      amcopy(&(smstiffmima_csr->update),&(actmaele->e.w1->stiff_mi_ma_csr->update));
      amcopy(&(smstiffmima_csr->a),&(actmaele->e.w1->stiff_mi_ma_csr->a));
      amcopy(&(smstiffmima_csr->ja),&(actmaele->e.w1->stiff_mi_ma_csr->ja));
      amcopy(&(smstiffmima_csr->ia),&(actmaele->e.w1->stiff_mi_ma_csr->ia));
    }
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->is_init     = actsmsolv->sysarray->ccf->is_init;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->is_factored = actsmsolv->sysarray->ccf->is_factored;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->ncall       = actsmsolv->sysarray->ccf->ncall;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->numeq_total = actsmsolv->sysarray->ccf->numeq_total;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->numeq       = actsmsolv->sysarray->ccf->numeq;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->nnz_total   = actsmsolv->sysarray->ccf->nnz_total;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->nnz         = actsmsolv->sysarray->ccf->nnz;
    actmaele->e.w1->stiff_mi_mi_ccf->ccf->comm        = actsmsolv->sysarray->ccf->comm;
  }
/*----------------------------------------------------------------------*/
/*   update displacements d'-> store at working-array of macro-element  */
/*   if it's the first time multiscale for this macroele,               */
/*   update can only take place after calc. of stif. and int.forces     */
/*----------------------------------------------------------------------*/
  if((statvar->praedictor == 2 || statvar->praedictor == 3) && actmaele->e.w1->firstinomegaprime == 1)
  {
    if(statvar->praedictor == 2)
    {/* -> use incremental solution! */
      for (i=0; i<ielma; i++)
      {
        DELTA_Dbar[2*i]   = actmaele->node[i]->sol_increment.a.da[0][0];
        DELTA_Dbar[2*i+1] = actmaele->node[i]->sol_increment.a.da[0][1];
      }
    }
    else if(statvar->praedictor == 3)
    {
      for (i=0; i<ielma; i++)/* -> use incremental solution! */
      {
        DELTA_Dbar[2*i]   = actmaele->node[i]->sol_residual.a.da[0][0];
        DELTA_Dbar[2*i+1] = actmaele->node[i]->sol_residual.a.da[0][1];
      }
    }
    /* calculate workforce:= - fint' - K-'bar *  (DELTA d-bar)          */
    for (i=0; i<numeq_sm; i++)
    {
      workforce[i] = 0.0;
      for(j=smstiffmima_csr->ia.a.iv[i]; j<smstiffmima_csr->ia.a.iv[i+1]; j++)
      {
        workforce[i] -= smstiffmima_csr->a.a.dv[j] * DELTA_Dbar[smstiffmima_csr->ja.a.iv[j]];
      }
      workforce[i] -= smintforcemi[i];
    }
    /* workforce = -(fint'+K'bar*DELTAdbar)->DIST-vector:actsmsolv->rhs[9]*/
    assemble_vec(actsmintra,&(actsmsolv->sysarray_typ[0]),
               &(actsmsolv->sysarray[0]),&(actsmsolv->rhs[9]),workforce,1.0);
    /* solve for DELTA d'= K''^-1 * workforce(use old LU-decomp of K'')  */
    actsmsolv->sysarray[0].ccf->reuse = 1;
    init_sol = 0;
    solver_control(actsmfield, 0,
                   actsmsolv,
                   actsmintra,
                 &(actsmsolv->sysarray_typ[0]),
                 &(actsmsolv->sysarray[0]),
                 &(actsmsolv->sol[9]),
                 &(actsmsolv->rhs[9]),
                   init_sol);
    /* update u' = u'+DELTAu' and store at macroelement                 */
    for (smnodei=0; smnodei<actsmfield->dis[0].numnp; smnodei++)
    {
      actsmnode = &(actsmfield->dis[0].node[smnodei]);
      ID = actsmnode->Id;
      for (j=0; j<numdf; j++)
      {
        dof = actsmnode->dof[j];
        if(dof >= numeq_sm) continue;
        if(statvar->praedictor == 2)
        {
          /* incrementelle Verschiebung DELTA d'  */
          actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[j] = actsmsolv->sol[9].vec.a.dv[dof];
          /* Verschiebung d' = d' alt + DELTA d'  */
          actmaele->e.w1->sm_nodaldata[ID].displ_mi[j] += actsmsolv->sol[9].vec.a.dv[dof];
        }
        else if(statvar->praedictor == 3)
        {
          /* residuelle Verschiebung DELTA d'  */
          if(istore == 1)
          {
            actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[j] = 0.0;
          }
          else
          {
            actmaele->e.w1->sm_nodaldata[ID].incre_displ_mi[j] += actsmsolv->sol[9].vec.a.dv[dof];
          }
          /* Verschiebung d' = d' alt + DELTA d'  */
          actmaele->e.w1->sm_nodaldata[ID].displ_mi[j] += actsmsolv->sol[9].vec.a.dv[dof];
        }
      }
    }
  }/* end: if(statvar->praedictor==2 &&actmaele->e.w1->firstinomegaprime==1) */
/*----------------------------------------------------------------------*/
  actmaele->e.w1->firstinomegaprime = 0;
/*               store history variables for next load-step             */
  if(istore == 1)
  {
    actmaele->e.w1->isinomegaprime = 1;
  }
}/*  end: else if(init!=3)                                              */
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s2ml_static_ke*/


/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
#endif
