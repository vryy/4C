/*!----------------------------------------------------------------------
\file
\brief contains the routine 'calfrq',
      to control dynamic eigenvalue analysis for optimization

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>


*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../headers/optimization.h"
#include "opt_prototypes.h"
/*!
\addtogroup OPTIMIZATION
*//*! @{ (documentation module open)*/


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                          al 08/02    |
 | pointer to allocate eigensolution variables                          |
 | dedfined in global_control.c                                         |
 | struct _ALLEIG       *alleig;                                        |
 *----------------------------------------------------------------------*/
extern ALLEIG              *alleig;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];

/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;


DOUBLE acttime;

/*----------------------------------------------------------------------*
 |  routine to control dynamic eigenvalue analysis              al 08/02|
 *----------------------------------------------------------------------*/
void calfrq(INT init)
{
INT             i,j;                /* some counters                    */
static INT      numeq;              /* number of equations on this proc */
static INT      numeq_total;        /* total number of equations */
INT             actsysarray;        /* active sparse system matrix in actsolv->sysarray[] */

INT             stiff_array;        /* indice of the active system sparse matrix */
INT             stiffa_copy;        /* indice of the active system sparse matrix */
INT             mass_array;         /* indice of the active system sparse matrix */

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structure cal_action enum */
STRUCT_DYNAMIC *sdyn;               /* pointer to structural dynamic input data */
CONTAINER     container;            /* contains variables defined in container.h */

enum   _SPARSE_TYP     *sysarray_typ;
union  _SPARSE_ARRAY   *sysarray;
struct _SKYMATRIX      *sky;

enum   _SPARSE_TYP     *sysarray_typ_copy;
union  _SPARSE_ARRAY   *sysarray_copy;
struct _SKYMATRIX      *sky_copy;

enum   _SPARSE_TYP     *sysarray_typ_mass;
union  _SPARSE_ARRAY   *sysarray_mass;
struct _SKYMATRIX      *sky_mass;

static DIST_VECTOR    *work;         /* working vectors */

static ARRAY    lumpedmass_a;        /* redundant vector of full length for lumped mass matrix */
static DOUBLE  *lumpedmass;          /* lumped mass matrix */

INT numvec;                          /* number of iteration (eigen)vectors */
static DIST_VECTOR    *itervec;             /* iteration (eigen)vectors */

static ARRAY           eigvalue_a;         /* vector for eigenvalues */
static DOUBLE         *eigvalue;
static ARRAY work1_a; static DOUBLE *work1;           /* working vectors */
static ARRAY work2_a; static DOUBLE *work2;           /* working vectors */
static ARRAY work3_a; static DOUBLE *work3;           /* working vectors */
static ARRAY rtolv_a; static DOUBLE *rtolv;           /* relative tolerances */


static ARRAY  eigfou_a; static DOUBLE   *eigfou;          /*  */
static ARRAY    rfou_a; static DOUBLE   *rfou;            /*  */
static ARRAY   rtfou_a; static DOUBLE   *rtfou;           /*  */

static ARRAY  ar_a; static DOUBLE   *ar;            /* system matrix a reduced problem */
static ARRAY  br_a; static DOUBLE   *br;            /* system matrix b reduced problem */
static ARRAY vec_a; static DOUBLE   *vec;           /* eigenvectors of reduced problem */

static ARRAY ind_a; static INT      *ind;           /* subspace specific */

INT    isiz;
INT    nn  ; /* real problem size (number of unknowns)            */
INT    nnm ; /* length of pointer array maxa                      */
INT    nwk ; /* size of matrix a                                  */
INT    nwm ; /* size of matrix b                                  */
INT    nc  ; /* size of reduced problem (num.of.iterationvectors) */
INT    nsta; /* restart flag                                      */
INT    nnz ; /* number of nonzero elements in system matrix       */
INT    isol;
INT    nroot;
INT    nitem;
INT    nsmax;
INT    ifss ;
INT    isub ;
DOUBLE    toleig;
DOUBLE    toljac;
DOUBLE    shift ;
INT    istldl;
INT    inits ;
INT    iprint;

INT    dof, nvec;
DOUBLE totmass;
#ifdef DEBUG
dstrc_enter("calfrq");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* if (par.myrank==0) out_cgs(1); /* initialize cgs */
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
sdyn        =   alldyn[0].sdyn;
container.isdyn   = 0;            /* static calculation */
container.actndis = 0;            /* only one discretisation */
container.fieldtyp  = actfield->fieldtyp;
#ifdef PARALLEL
actintra    = &(par.intra[0]);
/* if we are not parallel, we have to allocate an alibi intra-communicator structure */
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of structural eigen analysis, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*-------------------- set indice of stiffness and mass sparse matrices */
  stiff_array = 0;
  stiffa_copy = 1;
  mass_array  = 2;
  if(alleig->ilmp) actsolv->nsysarray=2; /*     lumped mass matrix */
  else             actsolv->nsysarray=3; /* consistent mass matrix */
  numvec = alleig->numvec;               /* number of eigenvectors */
/*--------------- stiff_array already exists, so copy the mask of it to */
/*---------------------------------------------------------- mass_array */
/* reallocate the vector of sparse matrices and the vectors*/
/* formerly lenght 1, now lenght 2*/
if(init==1)
{/*init==1*/
  isiz = actsolv->nsysarray*(INT)(sizeof(SPARSE_TYP));
  actsolv->sysarray_typ =
  (SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,isiz);
  if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

  actsolv->sysarray =
  (SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,isiz);
  if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

/*-copy the matrices sparsity mask from stiff_array to a copy of stiff_array */
  solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[stiff_array]),
                            &(actsolv->sysarray[stiff_array]),
                            &(actsolv->sysarray_typ[stiffa_copy]),
                            &(actsolv->sysarray[stiffa_copy]));
/*-copy the matrices sparsity mask from stiff_array to mass_array*/
  if(!alleig->ilmp)
  {
    solserv_alloc_cp_sparsemask(  actintra,
                              &(actsolv->sysarray_typ[stiff_array]),
                              &(actsolv->sysarray[stiff_array]),
                              &(actsolv->sysarray_typ[mass_array]),
                              &(actsolv->sysarray[mass_array]));
  }
/*---------------------------- get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[stiff_array]),
                     actsolv->sysarray_typ[stiff_array],
                     &numeq,
                     &numeq_total);


/*----------------------------------------allocate 1 dist. vectors 'rhs'*/
/*  these hold lumped mass matrix if nescessary */
  if(alleig->ilmp)
  {

    /* used by the element routines to assemble the lumped mass matrix */
    lumpedmass = amdef("lumpedmass",&lumpedmass_a,numeq_total,1,"DV");



    actsolv->nrhs = 1;
    solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
  }
/*--------------------------- there i one solution vectors to hold ...*/
  actsolv->nsol= 1;
  solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
/*--------------------------------------allocate three working vectors */
  solserv_create_vec(&work,3,numeq_total,numeq,"DV");
/*-----------------------------------------------allocate eigenvectors */
  solserv_create_vec(&itervec,1,numeq_total*numvec,numeq_total*numvec,"DV");
/*------------------------------------ create a vector for eigenvalues */
  eigvalue = amdef("eigenvalues",&eigvalue_a,numvec,1,"DV");
/*--------------------------- create some working vectors for subspace */
  work1 = amdef("work1",&work1_a,numeq_total,1,"DV");
  work2 = amdef("work2",&work2_a,numeq_total,1,"DV");
  work3 = amdef("work3",&work3_a,numvec     ,1,"DV");
  rtolv = amdef("rtolv",&rtolv_a,numvec     ,1,"DV");

  ar  = amdef("ar" ,&ar_a ,numeq_total,1,"DV");
  br  = amdef("br" ,&br_a ,numeq_total,1,"DV");
  vec = amdef("vec",&vec_a,numeq_total,1,"DV");

  ind = amdef("ind",&ind_a,numvec*3,1,"IV");

  eigfou = amdef("eigfou" ,&eigfou_a ,numvec,1,"DV");
  rfou   = amdef("rfou"   ,&rfou_a   ,numeq_total*numvec,1,"DV");
  rtfou  = amdef("rtfou"  ,&rtfou_a  ,numvec,1,"DV");
/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init phase
      is cheap in computation (can be costly in memory)
      There will be no solver call on mass or damping array.
*/
/*--------------------------------------------------- initialize solver */
  init=1;
  solver_control(actsolv, actintra,
                 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &work[0],
                 &work[1],
                 init);
  solver_control(actsolv, actintra,
                 &(actsolv->sysarray_typ[stiffa_copy]),
                 &(actsolv->sysarray[stiffa_copy]),
                 &work[0],
                 &work[1],
                 init);
  if(!alleig->ilmp)
  {
    solver_control(actsolv, actintra,
                 &(actsolv->sysarray_typ[mass_array]),
                 &(actsolv->sysarray[mass_array]),
                 &work[0],
                 &work[1],
                 init);
  }
/*----------------- init the assembly for stiffness and for mass matrix */
init_assembly(actpart,actsolv,actintra,actfield,stiff_array,0);
if(!alleig->ilmp)
{
  init_assembly(actpart,actsolv,actintra,actfield,mass_array,0);
}

/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action,&container);

goto end;
}/*init==1*/
/*------------------------------- init the dist sparse matrices to zero */
for (i=0; i<actsolv->nsysarray; i++)
  solserv_zero_mat(
                   actintra,
                   &(actsolv->sysarray[i]),
                   &(actsolv->sysarray_typ[i])
                  );
/*-------------------------------------------------------- init to zero */
  for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));
  for (i=0; i<3; i++) solserv_zero_vec(&(work[i]));
  for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));
  solserv_zero_vec(&(itervec[0]));
  amzero(&eigvalue_a);
  if(alleig->ilmp) amzero(&lumpedmass_a);

  amzero(&work1_a);
  amzero(&work2_a);
  amzero(&work3_a);
  amzero(&rtolv_a);
  amzero(&ar_a);
  amzero(&br_a);
  amzero(&vec_a);
  amzero(&ind_a);
  amzero(&eigfou_a);
  amzero(&rfou_a);
  amzero(&rtfou_a);
/*----------------------- call elements to calculate stiffness and mass */


*action = calc_struct_linstiffmass;
container.dvec         = lumpedmass;
container.dirich       = NULL;
container.global_numeq = actsolv->sol[0].numeq_total;
container.kstep        = 0;
container.actndis      = 0;            /* only one discretisation */
if(alleig->ilmp) *action = calc_struct_linstifflmass;
if(alleig->ilmp)
{
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);

/*
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,
                                  lumpedmass,
                                  NULL,numeq_total,0,action);
*/
}
else
{
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
/*  calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,
                       NULL,
                       NULL,numeq_total,0,action);*/
}
/*----------------------------------------------------------------------*/
sysarray_typ = &actsolv->sysarray_typ[stiff_array];
sysarray     = &actsolv->sysarray[stiff_array];

sysarray_typ_copy = &actsolv->sysarray_typ[stiffa_copy];
sysarray_copy     = &actsolv->sysarray[stiffa_copy];

sysarray_typ_mass = &actsolv->sysarray_typ[mass_array];
sysarray_mass     = &actsolv->sysarray[mass_array];

if(*sysarray_typ==skymatrix)
{
  sky      = sysarray->sky;
  sky_copy = sysarray_copy->sky;
  sky_mass = sysarray_mass->sky;

  if(sky_copy->A.fdim!=sky->A.fdim) dserror("action unknown");
  for (i=0; i<sky_copy->A.fdim; i++)  sky_copy->A.a.dv[i] = sky->A.a.dv[i];
}
/*----------------------------------------------------------------------*/
nn     = numeq_total;     /* real problem size (number of unknowns)            */
nnm    = sky->maxaf.fdim; /* length of pointer array maxa                      */
nwk    = sky->A.fdim;     /* size of matrix a                                  */
if(alleig->ilmp) nwm    = numeq_total;     /* size of matrix b                                  */
else nwm    = sky_mass->A.fdim;            /* size of matrix b                                  */

nc     = numvec;          /* size of reduced problem (num.of.iterationvectors) */
nsta   = 0;               /* restart flag                                      */
nnz    = 0;               /* number of nonzero elements in system matrix       */
isol   = 1;
nroot  = alleig->nroot;
nitem  = alleig->itemax;
nsmax  = 16;             /* maximum number of elemination */
ifss   = alleig->sturm;
isub   = alleig->subtyp;
toleig = alleig->toleig;
toljac = 1E-12; /* tolerance to be used for jacobi iteration (subproblem) */
shift  = alleig->shift;
istldl =0;
inits  =0;
iprint = alleig->ifctr;
/*----------------------------------------------------------------------*/
totmass=0.0;
for (i=0; i<nwm; i++) totmass +=lumpedmass[i];
/*----------------------------------------------------------------------*/
if(alleig->ilmp)
{
  sspace(
         sky->A.a.dv,
         sky_copy->A.a.dv,
         lumpedmass,
         sky->maxaf.a.iv,
         itervec->vec.a.dv,
         eigvalue,
         work1,
         work2,
         ar,
         br,
         vec,
         work3,
         rtolv,
         ind,
         &nn ,
         &nnm,
         &nwk,
         &nwm,
         &nc ,
         &nsta,
         &isol,
         &nnz,
         eigfou,
         rfou,
         &nroot,
         &nitem,
         &nsmax,
         &ifss,
         &isub,
         &toleig,
         &toljac,
         &shift,
         &istldl,
         &inits,
         &iprint
         );
}
else
{
  sspace(
         sky->A.a.dv,
         sky_copy->A.a.dv,
         sky_mass->A.a.dv,
         sky->maxaf.a.iv,
         itervec->vec.a.dv,
         eigvalue,
         work1,
         work2,
         ar,
         br,
         vec,
         work3,
         rtolv,
         ind,
         &nn ,
         &nnm,
         &nwk,
         &nwm,
         &nc ,
         &nsta,
         &isol,
         &nnz,
         eigfou,
         rfou,
         &nroot,
         &nitem,
         &nsmax,
         &ifss,
         &isub,
         &toleig,
         &toljac,
         &shift,
         &istldl,
         &inits,
         &iprint
         );
}




   /*------------------------------------------------ distribute result */
#ifdef PARALLEL
      MPI_Bcast(itervec,numeq_total*numvec,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif
/*--------------------------------------- save data for optimization ---*/
  for (i=0; i<opt->oeig->numeigv; i++)
  {
    opt->oeig->eigv[i] = eigvalue[i];
  }
/*--- scaling factor for eigen frequency optimization ------------------*/
for (i=0;i<opt->oeig->numeigv;i++) opt->oeig->eigs[i] = 0.0;
for (i=0;i<opt->oeig->numeigv;i++)
{
  for (j=0;j<numeq_total;j++)
  {
    opt->oeig->eigs[i] +=
             lumpedmass[j]
           * itervec->vec.a.dv[i*numeq_total+j]
           * itervec->vec.a.dv[i*numeq_total+j];
  }
}
/*------------------------------------- get result for graphical output */
  for (nvec=0; nvec<nroot; nvec++)
  {
    for (i=0; i<actsolv->sol[0].numeq; i++)
    {
        dof      = sky->update.a.iv[i]+nvec*actsolv->sol[0].numeq;
        actsolv->sol[0].vec.a.dv[i] = itervec->vec.a.dv[dof];
    }
    /*----------------------- return vectors of eigenforms to the nodes */
    solserv_result_total(actfield,actintra, &(actsolv->sol[0]),0,
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));

    opt_g_out(gr_disp);
  }
/*---------------------------------------------------------- get result */
  for (nvec=0; nvec<nroot; nvec++)
  {
    for (i=0; i<actsolv->sol[0].numeq; i++)
    {
        dof      = sky->update.a.iv[i]+nvec*actsolv->sol[0].numeq;
        actsolv->sol[0].vec.a.dv[i] = itervec->vec.a.dv[dof];
    }
    /*----------------------- return vectors of eigenforms to the nodes */
    solserv_result_total(actfield,actintra, &(actsolv->sol[0]),nvec,
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));
  }
/*---- ------------------------------------------------------------------*/
end:
/*---- ----------------------------------------------- cleaning up phase */
if(init==2)
{/*init==2*/
  if(alleig->ilmp)
  {
    amdel(&lumpedmass_a);
    solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  }
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);
  solserv_del_vec(&work,3);
}/*init==2*/
/*----------------------------------------------------------------------*/
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of calfrq */
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/





