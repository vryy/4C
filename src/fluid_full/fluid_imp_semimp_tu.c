/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fluid

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
static INT     rans=0;             /* flag, if vel. discr. is used     */
static INT     kapeps=1;           /* flag, if pres. discr. is used    */
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
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

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
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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

static FLUID_DYNAMIC *fdyn;
/*----------------------------------------------------------------------*
 | routine to control implicit and semi-implicit algorithms for fluid   |
 | problems combined with Newton and fixed point iteration schemes.     |
 | IN PROGRESS: ONE-STEP-THETA                                          |
 |              fixed point like iteration                              |
 |              only Euler no ALE                                       |
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief implicit and semi-implicit algorithms for fluid problems

<pre>                                                         genk 03/02

this routine conrols the implicit and semi-implicit algorithms for fluid
problems combined with different nonlinear iteration schemes.

time-discretisation:
fdyn->iop=2: Semi-Implicit-One-Step Method
fdyn->iop=3: Semi-Implicit-Two-Step Method
fdyn->iop=4: One-Step-Theta Scheme
fdyn->iop=5: Fractional-Step-Theta Scheme

see dissertation of W.A. Wall chapter 4.2 'Zeitdiskretisierung'

nonlinear iteration scheme:
fdyn->ite=0: no nonlinear iteration
fdyn->ite=1: fixed-point-like iteration
fdyn->ite=2: Newton iteration
fdyn->ite=3: fixed-point iteration

see dissertation chapter 4.3 'Linearisierung und Iteratonsverfahren'.

</pre>
\param *fdyn	 FLUID_DYNAMIC (i)

\return void
\warning up to now only the One-Step-Theta scheme combined with a
fixed-point-like iteration scheme is tested!

------------------------------------------------------------------------*/
void fluid_isi_tu(void)
{
INT             itnum;              /* counter for nonlinear iteration  */
INT             itnum1;             /* counter for nonlinear iteration  */
INT             itnum2;             /* counter for nonlinear iteration  */
INT             itnumke;            /* counter for nonlinear iteration  */
INT             itnum_n=1;          /* counter for nonlinear iteration  */
INT             itnum_check;        /* counter for nonlinear iteration  */
INT             conv_check_rans;
INT             i,kk;               /* simply a counter                 */
INT             numeq[2];           /* number of equations on this proc */
INT             numeq_total[2];     /* total number of equations        */
INT             numeq_oll;
INT             numeq_total_oll;
INT             init;               /* flag for solver_control call     */
INT             nsysarray=1;        /* two system matrix                */
INT             actsysarray;        /* number of actual sysarray        */
INT             k_array=0;          /* index of K-matrix in solver      */
INT             iststep=0;          /* counter for time integration     */
INT             nfrastep;           /* number of steps for fractional-
                                       step-theta procedure             */
INT             outstep=0;          /* counter for time integration     */
INT             pssstep=0;
INT             actcurve;           /* actual timecurve                 */
INT             converged=0;        /* convergence flag                 */
INT             steady=0;           /* flag for steady state            */
INT             actpos;             /* actual position in sol. history  */
DOUBLE          vrat,prat;          /* convergence ratios               */
DOUBLE          kapepsrat,lenghtrat;/* convergence ratios               */
DOUBLE          lower_limit_kappa;  /* convergence ratio                */
DOUBLE          lower_limit_eps;    /* convergence ratio                */

DOUBLE          t1,ts,te;	    /*					*/
DOUBLE          tes=0.0;            /*					*/
DOUBLE          tss=0.0;            /*					*/
FLUID_STRESS    str;

DIST_VECTOR    *rhs_ke;             /* distr. RHS for solving ke        */
DIST_VECTOR    *sol_ke;             /* distr. ke solution               */

SOLVAR         *actsolv;            /* pointer to active sol. structure */
SOLVAR         *kesolv;             /* solver for kappa-epsilon         */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */

ARRAY           ftimerhs_a;
DOUBLE         *ftimerhs;	    /* time - RHS			*/
ARRAY           ftimerhs_ke_a;
DOUBLE         *ftimerhs_ke;	    /* time - RHS			*/
ARRAY           ftimerhs_kappa_a;
DOUBLE         *ftimerhs_kappa;     /* time - RHS			*/
ARRAY           ftimerhs_epsilon_a;
DOUBLE         *ftimerhs_epsilon;   /* time - RHS			*/
ARRAY           ftimerhs_pro_a;
DOUBLE         *ftimerhs_pro;       /* time - RHS			*/
ARRAY           ftimerhs_pro_kappa_a;
DOUBLE         *ftimerhs_pro_kappa; /* time - RHS			*/
ARRAY           ftimerhs_pro_epsilon_a;
DOUBLE         *ftimerhs_pro_epsilon; /* time - RHS			*/

ARRAY           fiterhs_a;
DOUBLE         *fiterhs;	    /* iteration - RHS  		*/
ARRAY           fiterhs_ke_a;
DOUBLE         *fiterhs_ke;	    /* iteration - RHS  		*/

ARRAY           time_a;             /* stored time                      */

CONTAINER       container;      /* contains variables defined in container.h */

INT             initialisation=0;
INT             kapeps_yeah=0;
#ifdef DEBUG
dstrc_enter("fluid_isi_tu");
#endif


/*======================================================================*
 |                    I N I T I A L I S A T I O N                       |
 *======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
fdyn = alldyn[genprob.numff].fdyn;

actfield           = &(field[0]);
actsolv            = &(solv[0]);
actpart            = &(partition[0]);
action             = &(calc_action[0]);
container.fieldtyp = actfield->fieldtyp;
container.turbu    = fdyn->turbu;
fdyn->washvel    = 1.0;
str                = str_none;
fdyn->acttime    = ZERO;

/*------------------------- this is no relaxation parameter solution ---*/
container.is_relax = 0;

/*--------------------------------------------- set max. iterationsteps */
fdyn->itemax_ke= 5*fdyn->itemax;

/*---------------- if we are not parallel, we have to allocate an alibi *
  ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = fluid;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) goto end;

solv = (SOLVAR*)CCAREALLOC(solv,2*sizeof(SOLVAR));
/*-------------------------------------- set pointer to pressure solver */
actsolv = &(solv[rans]);
kesolv  = &(solv[kapeps]);
kesolv->fieldtyp = fluid;

/*-------------------------------------- now we set some default values */
#ifdef PARALLEL
#ifndef SPOOLES_PACKAGE
dserror("SPOOLES package is not compiled in");
#else
kesolv->solvertyp = SPOOLES_nonsym;
#endif
#else
kesolv->solvertyp = umfpack;
#endif
kesolv->parttyp = cut_elements;
kesolv->matrixtyp = oll_matrix;

/* -------------------------- create solver for pressure discretisation */
kesolv->nsysarray = 1;
kesolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(kesolv->nsysarray,sizeof(SPARSE_TYP));
kesolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(kesolv->nsysarray,sizeof(SPARSE_ARRAY));

kesolv->sysarray_typ[0] = oll;
kesolv->sysarray[0].oll = (OLL*)CCACALLOC(1,sizeof(OLL));


/*------------------------------- loop the matrices and intitialise them */
for(kk=0;kk<nsysarray;kk++)
{
actsysarray=kk;

/*------------------------------- init the dist sparse matrices to zero */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq[kk],
                   &numeq_total[kk]);

} /* end of loop over sys_arrays */

/*-------------------------------------- initialise solver for pressure */
numeq_total_oll = actfield->dis[kapeps].numeq;
oll_numeq(actfield, actpart, actintra, kapeps, &numeq_oll);

oll_open(kesolv->sysarray[0].oll, numeq_oll, numeq_total_oll,
	 actfield, actpart, actintra, kapeps);

numeq[kapeps] = numeq_oll;
numeq_total[kapeps] = numeq_total_oll;

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
printf("PROC  %3d | FIELD FLUID     | number of rans equations        : %10d \n",
        par.myrank,numeq[rans]);
printf("PROC  %3d | FIELD FLUID     | number of kapeps equations      : %10d \n",
        par.myrank,numeq[kapeps]);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD FLUID     | total number of rans equations  : %10d \n",numeq_total[rans]);
printf("          | FIELD FLUID     | total number of kapeps equations: %10d \n",numeq_total[kapeps]);
if (par.myrank==0) printf("\n\n");


/*--------------------------------------- allocate 1 dist. vector 'rhs' */
solserv_create_vec(&(actsolv->rhs),1,numeq_total[rans],numeq[rans],"DV");
solserv_zero_vec(&(actsolv->rhs[0]));

/*------------------------------------------- allocate 1 dist. solution */
solserv_create_vec(&(actsolv->sol),1,numeq_total[rans],numeq[rans],"DV");
solserv_zero_vec(&(actsolv->sol[0]));

/*-------------------------- allocate 1 dist. vector 'rhs' for ke       */
solserv_create_vec(&rhs_ke,1,numeq_total[kapeps],numeq[kapeps],"DV");
solserv_zero_vec(rhs_ke);

/*-------------------------- allocate 1 dist. vector 'sol' for ke       */
solserv_create_vec(&sol_ke,1,numeq_total[kapeps],numeq[kapeps],"DV");
solserv_zero_vec(sol_ke);

/*--------------- allocate one redundant vector ftimerhs of full lenght */
/*        this is used by the element routines to assemble the  Time RHS*/
ftimerhs         = amdef("ftimerhs",   &ftimerhs_a,numeq_total[rans],1,"DV");
ftimerhs_ke      = amdef("ftimerhs_ke",&ftimerhs_ke_a,numeq_total[kapeps],1,"DV");
ftimerhs_kappa   = amdef("ftimerhs_kappa",&ftimerhs_kappa_a,numeq_total[kapeps],1,"DV");
ftimerhs_epsilon = amdef("ftimerhs_epsilon",&ftimerhs_epsilon_a,numeq_total[kapeps],1,"DV");
ftimerhs_pro     = amdef("ftimerhs_pro",&ftimerhs_pro_a,numeq_total[kapeps],1,"DV");
ftimerhs_pro_kappa   = amdef("ftimerhs_pro_kappa",&ftimerhs_pro_kappa_a,numeq_total[kapeps],1,"DV");
ftimerhs_pro_epsilon = amdef("ftimerhs_pro_epsilon",&ftimerhs_pro_epsilon_a,numeq_total[kapeps],1,"DV");

/*---------------  allocate one redundant vector fiterhs of full lenght */
/*   this is used by the element routines to assemble the  Iteration RHS*/
fiterhs    = amdef("fiterhs",&fiterhs_a,numeq_total[rans],1,"DV");
fiterhs_ke = amdef("fiterhs_ke",&fiterhs_ke_a,numeq_total[kapeps],1,"DV");

/*--------------------------- allocate one vector for storing the time */
amdef("time",&time_a,1000,1,"DV");

/*--------------------------------------------- initialise fluid field */
fluid_init(actpart,actintra,actfield, 0,action,&container,4,str);
fluid_init_tu(actfield);
actpos=0;

/*--------------------------------------- init all applied time curves */
for (actcurve=0; actcurve<numcurve; actcurve++)
   dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);

/*-------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield);
/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init
      phase is cheap in computation (can be costly in memory)
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);

solver_control(kesolv, actintra,
               &(kesolv->sysarray_typ[0]),
               &(kesolv->sysarray[0]),
               sol_ke,
               rhs_ke,
               init);

/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,k_array,rans);
init_assembly(actpart,kesolv,actintra,actfield,0,kapeps);

/*---------------------------------- allocate fluid integration data ---*/
alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);
actpos++;

fluid_cons();

/*======================================================================*
 |                         T I M E L O O P                              |
 *======================================================================*/
/* nodal solution history fluid field:                                  *
 * sol[0][j]           ... initial data 				*
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[0][j] ... solution at time (n-1)			*
 * sol_increment[1][j] ... solution at time (n) 			*
 * sol_increment[2][j] ... solution at time (n+g)			*
 * sol_increment[3][j] ... solution at time (n+1)			*
 *======================================================================*/
timeloop:
fdyn->step++;
iststep++;
/*------------------------------------------ check (starting) algorithm */
if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(&nfrastep,0);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons();
if (kapeps_yeah == 1) fluid_tcons_tu();

fdyn->acttime += fdyn->dta;

/*------ flag and counter for convergence check for kappa-epsilon model*/
conv_check_rans = 0;
itnum_check = 1;

nonlniter_check:
container.actndis=0;

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout();

/*--------------------- set dirichlet boundary conditions for  timestep */
fluid_setdirich(actfield,3);

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf("----------------------------------------------------------------\n");
   printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -| \n");
}
itnum=1;
/*======================================================================*
 |           N O N L I N E A R   I T E R A T I O N                      |
 *======================================================================*/
nonlniter:
/*------------------------- calculate constants for nonlinear iteration */
fluid_icons(itnum);

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[k_array]),
                 &(actsolv->sysarray_typ[k_array]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();
container.ftimerhs     = ftimerhs;
container.fiterhs      = fiterhs;
container.global_numeq = numeq_total[rans];
container.nii          = fdyn->nii;
container.nif          = fdyn->nif;
container.nim          = 0;
container.kstep        = 0;
calelm(actfield,actsolv,actpart,actintra,k_array,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[k_array]),
             &(actsolv->sysarray[k_array]),
             &(actsolv->rhs[0]),
             ftimerhs,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[k_array]),
             &(actsolv->sysarray[k_array]),
             &(actsolv->rhs[0]),
             fiterhs,
             1.0
             );

/*-------------------------------------------------------- solve system */
init=0;
t1=ds_cputime();
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);
ts=ds_cputime()-t1;
tss+=ts;

/*-- set flags for stability parameter evaluation and convergence check */
fdyn->ishape=0;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre(actfield, 0,actintra,&(actsolv->sol[k_array]),3,
                     &(actsolv->sysarray[k_array]),
                     &(actsolv->sysarray_typ[k_array]),
		         &vrat,&prat,NULL);

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(vrat,prat,ZERO,itnum,te,ts);

/*--------------------- check if nonlinear iteration has to be finished */
if (converged==0)
{
   itnum++;
   goto nonlniter;
}

/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration   for rans                           |
 *----------------------------------------------------------------------*/

/*======================================================================*
                 K A P P A - E P S I L O N
   if possible: calc. until steady state before kappa, epsilon and
                     eddy-viscosity are calcul.
 *======================================================================*/
if (kapeps_yeah == 1)
{
container.actndis=1;

/*-------------- set dirichlet boundary conditions const. for timestep */
fluid_setdirich_tu(actfield,&lower_limit_kappa,&lower_limit_eps);

/*---------------- initialise the turbulence variables for calculation */
if (initialisation == 0)
{
fdyn->stepke++;
fluid_set_check_tu(actfield,lower_limit_kappa);
/*-------------------------------- amzero timerhs for kappa end epsilon */
amzero(&ftimerhs_kappa_a);
amzero(&ftimerhs_epsilon_a);

initialisation = 1;
}

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout_tu();

/*-------------------- counter for charact. lenght and production terms */
itnum2 = 1;
itnumke= 1;

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_pro_kappa_a);
amzero(&ftimerhs_pro_epsilon_a);

/*======================================================================*
 | do the nonlinear iteration  for turbulence  (kappa)                  |
 *======================================================================*/
kappa:

fdyn->kapeps_flag=0;
itnum1=1;

if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf(" ______________________________________________________________ \n");
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- kappa  error -| \n");
}

if(fdyn->kapeps_flag==1)
{
/*======================================================================*
 | do the nonlinear iteration  for turbulence  (epsilon)                |
 *======================================================================*/
epsilon:

fdyn->kapeps_flag=1;
itnum1=1;
/*----------------------------------------------------------------------*/

if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- epsil  error -| \n");
}

} /* end fdyn->kapeps_flag==1 */

/*---------------------------------------------------------------------*/
nonlniter1:
fluid_icons_tu(itnum1,itnumke,itnum_n);
/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(rhs_ke);
solserv_zero_mat(actintra,&(kesolv->sysarray[0]),
                 &(kesolv->sysarray_typ[0]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_ke_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();

if(fdyn->kapeps_flag==0)
{
 container.ftimerhs     = ftimerhs_kappa;
 container.ftimerhs_pro = ftimerhs_pro_kappa;
}
if(fdyn->kapeps_flag==1)
{
 container.ftimerhs     = ftimerhs_epsilon;
 container.ftimerhs_pro = ftimerhs_pro_epsilon;
}
container.fiterhs      = fiterhs_ke;
container.global_numeq = numeq_total[kapeps];
container.niturbu_pro  = fdyn->niturbu_pro;
container.niturbu_n    = fdyn->niturbu_n;
container.nii          = 0;
container.nim          = 0;
container.nif          = 0;
container.kstep        = 0;
calelm(actfield,kesolv,actpart,actintra,0,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;

if(fdyn->kapeps_flag==0)
{
 ftimerhs_ke = ftimerhs_kappa;
 ftimerhs_pro= ftimerhs_pro_kappa;
}
if(fdyn->kapeps_flag==1)
{
 ftimerhs_ke = ftimerhs_epsilon;
 ftimerhs_pro= ftimerhs_pro_epsilon;
}
/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
assemble_vec(actintra,
             &(kesolv->sysarray_typ[0]),
             &(kesolv->sysarray[0]),
             rhs_ke,
             ftimerhs_ke,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(kesolv->sysarray_typ[0]),
             &(kesolv->sysarray[0]),
             rhs_ke,
             fiterhs_ke,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(kesolv->sysarray_typ[0]),
             &(kesolv->sysarray[0]),
             rhs_ke,
             ftimerhs_pro,
             1.0
             );

/*-------------------------------------------------------- solve system */
init=0;
t1=ds_cputime();
solver_control(kesolv, actintra,
               &(kesolv->sysarray_typ[0]),
               &(kesolv->sysarray[0]),
               sol_ke,
               rhs_ke,
               init);
ts=ds_cputime()-t1;
tss+=ts;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre_tu(actfield,actintra,sol_ke,3,
                     &(kesolv->sysarray[0]),
                     &(kesolv->sysarray_typ[0]),
		     &kapepsrat,fdyn,lower_limit_kappa,lower_limit_eps);

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck_tu(kapepsrat,itnum1,te,ts);

/*---------- check if nonlinear iteration for kapeps has to be finished */
if (converged==0)
{
   itnum1++;
   goto nonlniter1;
}

if (fdyn->kapeps_flag==0)
{
   fluid_eddy_update(actfield,sol_ke);
   goto epsilon;
}

itnumke++;
itnum_n++;

/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration  for turbulence (kappa || epsilon)   |
 *----------------------------------------------------------------------*/
if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf("|____________|___________________|________________| \n");
   printf("\n");
}

/*--------------- iteration convergence check for characteristic lenght */
if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- lenght error -| \n");
}

fluid_eddy_update(actfield,sol_ke);
fluid_lenght_update(actfield,sol_ke,&lenghtrat);

/*--------------------- check if nonlinear iteration has to be finished */
converged = fluid_convcheck_tu(lenghtrat,itnum2,te,ts);

if (converged==0)
{
   itnum2++;
   goto kappa;
}

/*----------------- write eddy (n) to eddy (n+g) due to production-term */
fluid_eddy_pro(actfield);

if (fdyn->itnorm!=fncc_no && par.myrank==0)
{
   printf("|____________|___________________|________________| \n");
   printf(" ______________________________________________________________ \n");
   printf("\n");
}

/*------------------------- check if rans has finished due to kappa-eps */
if (conv_check_rans == 1)
{
 converged = fluid_convcheck_test(actfield,itnum_check);

 if (converged!=0)
 {
/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
  fluid_copysol_tu(actfield,3,1,0);

/*------------------------------------------- counter for timerhs terms */
  itnum_n = 1;
  fdyn->stepke++;

/*-------------------------------- amzero timerhs for kappa end epsilon */
  amzero(&ftimerhs_kappa_a);
  amzero(&ftimerhs_epsilon_a);

/*---------------------------------------- end outer-iteration for RANS */
  goto endrans;
 }

 itnum_check++;
}
/*------------------------------ copy solution from rans for conv-check *
                                  for rans due to kappa-eps             */
fluid_copysol_test(actfield,3,1);

conv_check_rans = 1;

/*--------------------- check if nonlinear iteration for rans converge */
goto nonlniter_check;
} /* endif kapeps_yeah */

/*========================== E N D: ====================================*
 /                  K A P P A - E P S I L O N                          /
*========================== E N D: =====================================*/
endrans:

/*-------------------------------------------------- steady state check */
if (fdyn->stchk==iststep)
{
   iststep=0;
   steady = fluid_steadycheck(actfield,numeq_total[rans]);
}

/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,3,1);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]
           and transform kinematic to real pressure --------------------*/
solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol,3,actpos);
fluid_transpres(actfield,0,0,actpos,fdyn->numdf-1,0);

fluid_copysol_tu(actfield,3,actpos,1);

if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;

/*-------- calculate wall shear velocity and c_f (for TURBULENCE MODEL) */
   fdyn->washvel  = ZERO;
   container.actndis=0;
   *action = calc_fluid_shearvelo;
   container.nii= 0;
   container.nim= 0;
   container.nif= 0;
   calelm(actfield,actsolv,actpart,actintra,k_array,-1,
          &container,action);
#ifdef PARALLEL
fluid_reduceshstr(actintra,actfield);
#endif
if (par.myrank==0) printf("wall shear velocity: %10.3E \n",fdyn->washvel);

/*------------------------------------------------ print out to .out */
 out_sol(actfield,actpart,actintra,fdyn->step,actpos);
/*------------------------------------------------ print out to .tur */
 out_fluidtu(actfield,actintra,fdyn->step,actpos);

#ifdef PARALLEL
fluid_nullshstr(actintra,actpart,actfield);
#endif
}
if (pssstep==fdyn->uppss && ioflags.fluid_vis_file==1 && par.myrank==0)
{
   pssstep=0;
/*----------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim) amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = fdyn->acttime;
   actpos++;
}

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime && steady==0)
{
/*------ check time and number of steps if steady state cannot achieved */
if (kapeps_yeah==0 && (fdyn->step >= 0.5*fdyn->nstep || fdyn->acttime > 0.5*fdyn->maxtime))
{
   kapeps_yeah = 1;
}

   goto timeloop;
}

/*------------------------ if steady state was achieved goto kappa-eps */
if (kapeps_yeah==0 && steady==1)
{
   kapeps_yeah = 1;
   steady = 0;
   goto timeloop;
}
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/

/*======================================================================*
 |                      F I N A L I S I N G                             |
 *======================================================================*/
if (pssstep==0) actpos--;
/*------------------------ print out solution to .out and to .tur file */
if (outstep!=0 && ioflags.fluid_sol_file==1){
 out_sol(actfield,actpart,actintra,fdyn->step,actpos);
 out_fluidtu(actfield,actintra,fdyn->step,actpos);}

/*----------------------------- print out solution to 0.flavia.res file */
if (ioflags.fluid_sol_gid==1 && par.myrank==0)
{
    for(i=0;i<actpos+1;i++)
    {
        out_gid_sol("velocity",actfield,actintra,i,i,ZERO);
        out_gid_sol("pressure",actfield,actintra,i,i,ZERO);
    }
}

/*------------------------------------ print out solution to 0.pss file */
if (ioflags.fluid_vis_file==1 && par.myrank==0)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = fdyn->acttime;
   }
    visual_writepss(actfield,actpos+1,&time_a);
}

/*---------------------------------- print total CPU-time to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
{
if (par.myrank==i)
{
printf("\n");
printf("PROC  %3d | FIELD FLUID     | total time element for calculations: %10.3E \n",
        par.myrank,tes);
printf("PROC  %3d | FIELD FLUID     | total time for solver              : %10.3E \n",
        par.myrank,tss);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&ftimerhs_ke_a);
amdel(&ftimerhs_kappa_a);
amdel(&ftimerhs_epsilon_a);
amdel(&ftimerhs_pro_a);
amdel(&ftimerhs_pro_kappa_a);
amdel(&ftimerhs_pro_epsilon_a);
amdel(&fiterhs_a);

if (par.myrank==0 && ioflags.fluid_vis_file==1)
amdel(&time_a);

solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&rhs_ke,1);
solserv_del_vec(&sol_ke,1);
/*----------------------------------------------------------------------*/
#ifndef PARALLEL
CCAFREE(actintra);
#endif

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_isi_tu */

#endif
