/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fluid

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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

/*----------------------------------------------------------------------*
 | routine to control implicit and semi-implicit algorithms for fluid   |
 | problems combined with Newton and fixed point iteration schemes.     |
 | IN PROGRESS: ONE-STEP-THETA                                          |
 |              fixed point like iteration                              |
 |              only Euler no ALE                                       |
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

static FLUID_DYNAMIC *fdyn;

/*!---------------------------------------------------------------------
\brief implicit and semi-implicit algorithms for fluid problems

<pre>                                                         genk 03/02

this routine conrols the implicit and semi-implicit algorithms for fluid
problems combined with different nonlinear iteration schemes.

Some features have been added:                               chfoe 01/04
Generalised Alpha time integration and BDF2 as well as adaptive time
stepping has been implemented
There is now a choice for two different ways of treating the time rhs.
a) fdyn->time_rhs == 1: -> use 'classic' time rhs term as described in
                           the dissertation of W.A. Wall, time rhs is
                           calculated ONCE per time step
b) fdyn->time_rhs == 0: -> 'time rhs' is calculated by means of the
                           elemental mass matrices and reduces to an
                           assembled matrix vector product. It is
                           calculated newly within every iteration.

time-discretisation:
fdyn->iop=1: Generalised Alpha time integration
fdyn->iop=2: Semi-Implicit-One-Step Method \  not implemented
fdyn->iop=3: Semi-Implicit-Two-Step Method /
fdyn->iop=4: One-Step-Theta Scheme
fdyn->iop=5: Fractional-Step-Theta Scheme
fdyn->iop=7: BDF2 (2nd order backward differencing)

see dissertation of W.A. Wall chapter 4.2 'Zeitdiskretisierung'

nonlinear iteration scheme:
fdyn->ite=0: no nonlinear iteration
fdyn->ite=1: fixed-point-like iteration
fdyn->ite=2: Newton iteration
fdyn->ite=3: fixed-point iteration

see dissertation chapter 4.3 'Linearisierung und Iteratonsverfahren'.

</pre>

\return void
\warning up to now only the One-Step-Theta scheme combined with a
fixed-point-like iteration scheme is tested!

------------------------------------------------------------------------*/
void fluid_isi(void)
{
INT             itnum;              /* counter for nonlinear iteration  */
INT             i;                  /* simply a counter                 */
INT             numeq;              /* number of equations on this proc */
INT             numeq_total;        /* total number of equations        */
INT             init;               /* flag for solver_control call     */
INT             actsysarray=0;      /* number of actual sysarray        */
INT             outstep=0;          /* counter for output control       */
INT             resstep=0;          /* counter for output control       */
INT             pssstep=0;	    /* counter for output control	*/
INT             restartstep=0;
INT             iststep=0;          /* counter for time integration     */
INT             nfrastep;           /* number of steps for fractional-
                                       step-theta procedure             */
INT             actcurve;           /* actual timecurve                 */
INT             converged=0;        /* convergence flag                 */
INT             steady=0;           /* flag for steady state            */
INT             actpos;             /* actual position in sol. history  */
INT             restart;
INT		repeat=0;	    /* flag, repeat time step?		*/
INT		repeated=0;	    /* has time step been repeated?	*/
DOUBLE          vrat,prat;          /* convergence ratios               */
DOUBLE          t1,t2,ts,te,tt;	    /*					*/
DOUBLE          tes=0.0;            /*					*/
DOUBLE          tss=0.0;            /*					*/
DOUBLE          tts=0.0;            /*					*/
DOUBLE          fact1,fact2;
FLUID_STRESS    str;

SOLVAR         *actsolv;            /* pointer to active sol. structure */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */

ARRAY           frhs_a;
DOUBLE         *frhs;	            /* iteration - RHS                  */

CONTAINER       container;          /* contains variables defined in container.h */
FILE           *out = allfiles.out_out;

ARRAY_POSITION *ipos;

#ifdef BINIO
BIN_OUT_FIELD   out_context;
#endif

#ifdef DEBUG
dstrc_enter("fluid_isi");
#endif


/*======================================================================*
 |                    I N I T I A L I S A T I O N                       |
 *======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
fdyn = alldyn[genprob.numff].fdyn;

actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
restart     = genprob.restart;
container.actndis   = 0;
container.turbu     = fdyn->turbu;
container.fieldtyp  = actfield->fieldtyp;

ipos = &(actfield->dis[0].ipos);

/* set flag for stress evaluation */
str         = str_none;


fdyn->acttime=ZERO;

/*---------------- if we are not parallel, we have to allocate an alibi *
  ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
actintra->intra_fieldtyp = fluid;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) goto end;

/*------------------------------------ prepare lift&drag calculation ---*/
if (fdyn->liftdrag==ld_stress)
  str         = str_liftdrag;
if (fdyn->liftdrag==ld_nodeforce)
  fluid_liftdrag(-1,action,&container,actfield,
                 actsolv,actpart,actintra,ipos);

if (ioflags.fluid_stress==1)
  str         = str_all;

/*------------------------------- init the dist sparse matrices to zero */
solserv_zero_mat(
                 actintra,
                &(actsolv->sysarray[actsysarray]),
                &(actsolv->sysarray_typ[actsysarray])
                 );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
if (par.myrank==i)
printf("PROC  %3d | FIELD FLUID     | number of equations      : %10d \n",
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD FLUID     | total number of equations: %10d \n",numeq_total);
if (par.myrank==0) printf("\n\n");




/* write general data to .out */
if (par.myrank==0)
{
  fprintf(out,"max. values:\n");
  fprintf(out,"============\n");


  /* table head */
  fprintf(out," time |            |fluid| fluid error in ");

  switch(fdyn->itnorm)
  {
    case fncc_Linf: /* infinity norm */
      fprintf(out,"inf-norm");
      break;
    case fncc_L1: /* L_1 norm */
      fprintf(out,"L_1-norm");
      break;
    case fncc_L2: /* L_2 norm */
      fprintf(out,"L_2-norm");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
  }  /* switch(fdyn->itnorm) */

  fprintf(out," | steady state in ");

  switch(fdyn->stnorm)
  {
    case fnst_Linf: /* infinity norm */
      fprintf(out,"inf-norm|");
      break;
    case fnst_L1: /* L_1 norm */
      fprintf(out,"L_1-norm|");
      break;
    case fnst_L2: /* L_2 norm */
      fprintf(out,"L_2-norm|");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
  }  /* switch(fdyn->itnorm) */
  fprintf(out,"    total   |\n");

  fprintf(out," step |  sim. time | ite |     vel.   |     pre.   |     vel.   |     pre.   | calc. time |\n");
  fprintf(out,"-------------------------------------------------------------------------------------------\n");



  fprintf(out,"%5d | %10.3E | %3d |        %10.3E       |        %10.3E       |            |\n",
      fdyn->nstep,fdyn->maxtime,fdyn->itemax,fdyn->ittol,fdyn->sttol);
  fprintf(out,"-------------------------------------------------------------------------------------------\n");



  fprintf(out,"\n\ntimeloop:  ");

  switch(fdyn->iop)
  {
    case 1:
      fprintf(out,"Generalised Alpha\n");
      break;
    case 4:
      fprintf(out,"One-Step-Theta\n");
      break;
    case 7:
      fprintf(out,"BDF2\n");
      break;
    default:
      dserror("parameter out of range: IOP\n");
  }  /* switch(fdyn->iop) */

  fprintf(out,"=========\n");


  /* table head */
  fprintf(out," time |            |fluid| fluid error in ");

  switch(fdyn->itnorm)
  {
    case fncc_Linf: /* infinity norm */
      fprintf(out,"inf-norm");
      break;
    case fncc_L1: /* L_1 norm */
      fprintf(out,"L_1-norm");
      break;
    case fncc_L2: /* L_2 norm */
      fprintf(out,"L_2-norm");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
  }  /* switch(fdyn->itnorm) */

  fprintf(out," | steady state in ");

  switch(fdyn->stnorm)
  {
    case fnst_Linf: /* infinity norm */
      fprintf(out,"inf-norm|");
      break;
    case fnst_L1: /* L_1 norm */
      fprintf(out,"L_1-norm|");
      break;
    case fnst_L2: /* L_2 norm */
      fprintf(out,"L_2-norm|");
      break;
    default:
      dserror("Norm for nonlin. convergence check unknown!!\n");
  }  /* switch(fdyn->itnorm) */
  fprintf(out,"    total   |\n");

  fprintf(out," step |  sim. time | ite |     vel.   |     pre.   |     vel.   |     pre.   | calc. time |\n");
  fprintf(out,"-------------------------------------------------------------------------------------------\n");

  fflush(out);


}  /* if (par.myrank==0) */





/*---------------------------------------- allocate dist. vectors 'rhs' */
actsolv->nrhs = 1;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->rhs[0]));

/*---------------------------------- allocate dist. solution vectors ---*/
actsolv->nsol = 1;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->sol[0]));

/*---------------  allocate one redundant vector frhs of full lenght ---*/
/*   this is used by the element routines to assemble the  Iteration RHS*/
/* is this really necessary??? */
frhs = amdef("frhs",&frhs_a,numeq_total,1,"DV");

/*--------------------------------------------- initialise fluid field */
if (restart != 0)
{
  if (fdyn->init>0)
    dserror("Initial field either by restart, or by function or from file ...\n");
  else
  {
    fdyn->resstep=genprob.restart;
    fdyn->init=2;
  }
}

fluid_init_pos_euler(ipos);
fluid_init(actpart,actintra,actfield, 0,action,&container,8,ipos,str);
actpos=0;

/*------------------------------------ initialize multilevel algorithm */
#if defined(FLUID2_ML) || defined(FLUID3_ML)
if (fdyn->mlfem==1) fluid_ml_init(actfield);
#endif

/*--------------------------------------- init all applied time curves -*/
for (actcurve=0; actcurve<numcurve; actcurve++)
  dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);

/*-------------------------------------- init the dirichlet-conditions -*/
fluid_initdirich(actfield, ipos);

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
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);

/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);

/*---------------------------------- allocate fluid integration data ---*/
alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   actfield, actpart, actintra, 0);

#endif

/*--------------------------------------------- calculate nodal normals */
fluid_cal_normal(actfield,1,action);

/*------------------------------------------------- define local co-sys */
fluid_locsys(actfield,fdyn);

/*-------------------------------------- print out initial data to .out */
if (ioflags.output_out==1 && ioflags.fluid_sol==1 && par.myrank==0)
  out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*------------------------------- print out initial data to .flavia.res */
if (ioflags.output_gid==1 && par.myrank==0)
{
  if (ioflags.fluid_sol==1)
  {
    out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
    out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }
  if (ioflags.fluid_stress==1)
  {
    out_gid_sol("stress",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }
}

/*-------------------------------------------- write solution to binary */
#ifdef BINIO
if (ioflags.output_bin==1)
{
  if (ioflags.fluid_sol==1)
  {
    out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY | OUTPUT_PRESSURE);
  }
  if (ioflags.fluid_stress==1)
  {
    out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_STRESS);
  }
}
#endif

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
{
   out_monitor(actfield,genprob.numff,ZERO,1);
   monitoring(actfield,genprob.numff,actpos,fdyn->acttime);
}

/*---------- calculate time independent constants for time algorithm ---*/
fluid_cons();

dsmemreport();

/*======================================================================*
 |                         T I M E L O O P                              |
 *======================================================================*/
/* nodal solution history fluid field:                                  *
 * sol[0][j]           ... initial data 				*
 * sol[0...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[position][j] ... solution at some time level time      *
 * position flags:                                                      *
 *  ipos->velnp ... velocity at time (n+1)                               *
 *  ipos->veln  ... velocity at time (n)                                 *
 *  ipos->velnm ... velocity at time (n-1)                               *
 *  ipos->accn  ... acceleration at time (n)                             *
 *  ipos->accnm ... acceleration at time (n-1)                           *
 *  ipos->hist  ... old solution data depending on time integration      *
 *  ipos->pred  ... predicted solution for new time level                *
 *  ipos->terr  ... local truncation error                               *
 *======================================================================*/
timeloop:
t2=ds_cputime();

fdyn->step++;
iststep++;

if (par.myrank==0)
  fprintf(out,"%5d | %10.3E |",fdyn->step,fdyn->acttime);

/*------------------------------------------ check (starting) algorithm */
if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(&nfrastep,0);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons();

/*-------------------------------------------- set new absolute time ---*/
if (fdyn->iop == 1)/* generalised alpha is solved for time n+alpha_f 	*/
  fdyn->acttime += fdyn->dta * fdyn->alpha_f;
else
  fdyn->acttime += fdyn->dta;

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout();

/*------------------------ predictor step for adaptive time stepping ---*/
if (fdyn->adaptive && fdyn->step > 1)
{
  fluid_predictor(actfield,ipos,fdyn->iop);
}
else
{
/* do explicit predictor step to start iteration from better value */
  fact1 = fdyn->dta*(1.0+fdyn->dta/fdyn->dtp);
  fact2 = DSQR(fdyn->dta/fdyn->dtp);
  solserv_sol_add(actfield,0,node_array_sol_increment,
                             node_array_sol_increment,
                             ipos->accn,ipos->velnp, fact1);
  solserv_sol_add(actfield,0,node_array_sol_increment,
                             node_array_sol_increment,
                             ipos->veln,ipos->velnp,-fact2);
  solserv_sol_add(actfield,0,node_array_sol_increment,
                             node_array_sol_increment,
                             ipos->velnm,ipos->velnp,fact2);
}
/*-------- set dirichlet boundary conditions to sol_increment[velnp] ---*/
fluid_setdirich(actfield,ipos,ipos->velnp);
/*fluid_setdirich_cyl(actfield);*/

/*------------------------------------ prepare time rhs in mass form ---*/
fluid_prep_rhs(actfield, ipos);

/*------------------------------------- start time step on the screen---*/
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
/*fluid_icons(itnum); this is not needed any more */

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));

/*------------------------------------------- re-initialise neumann bcs */
inherit_design_dis_neum(&(actfield->dis[0]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&frhs_a);

/*-------------- form incremental matrices, residual and element forces */
#ifdef PERF
  perf_begin(81);
#endif

*action = calc_fluid;
t1=ds_cputime();
container.dvec         = NULL;
container.frhs         = frhs;
container.global_numeq = numeq_total;
container.nii          = fdyn->nii;
container.kstep        = 0;
container.is_relax     = 0;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;

#ifdef PERF
  perf_end(81);
#endif

/*--------------------------------------------------------- add rhs: ---*/
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[0]),
             frhs,
             1.0
             );
/*-------------------------------------------------------- solve system */
#ifdef PERF
  perf_begin(80);
#endif



#if 0
  fprintf(sol,"\nval:\n====\n");
  for (i=0; i<actsolv->sysarray[actsysarray].msr->val.fdim; i++)
  fprintf(sol,"%12.4e\n",actsolv->sysarray[actsysarray].msr->val.a.dv[i]);

  fprintf(sol,"\nrhs:\n====\n");
  for (i=0; i<actsolv->rhs[0].vec.fdim; i++)
  fprintf(sol,"%12.4e\n",actsolv->rhs[0].vec.a.dv[i]);
#endif


init=0;
t1=ds_cputime();
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);
ts=ds_cputime()-t1;
tss+=ts;


#if 0
  fprintf(sol,"\nsol:\n====\n");
  for (i=0; i<actsolv->sol[0].vec.fdim; i++)
  fprintf(sol,"%12.4e\n",actsolv->sol[0].vec.a.dv[i]);
#endif


#ifdef PERF
  perf_end(80);
#endif

/*-- set flags for stability parameter evaluation and convergence check */
fdyn->ishape=0;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre(actfield, 0,actintra,&(actsolv->sol[0]),
                     ipos->velnp,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
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
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/

/*---------- extrapolate from n+alpha_f to n+1 for generalised alpha ---*/
if (fdyn->iop == 1)
{
  solserv_sol_zero(actfield,0,node_array_sol_increment,ipos->accnm);
  solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->velnp,ipos->accnm,1.0/fdyn->alpha_f);
  solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->veln,ipos->accnm,1.0-1.0/fdyn->alpha_f);
  solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->accnm,ipos->velnp);

  fdyn->acttime += fdyn->dta * (1.0 - fdyn->alpha_f);
}
/*-------------------------- check time step size in adaptive regime ---*/
if(fdyn->adaptive)
{
  if (fdyn->step > 1)
  {
    /*------------------------ evaluate local truncation error (LTE) ---*/
    fluid_lte(actfield,ipos,fdyn->iop);

    /*------------------ evaluate norm of LTE and new time step size ---*/
    fluid_lte_norm(actpart,actintra,ipos,
                   &iststep,&repeat,&repeated,itnum);
    if (repeat)
    {
      repeated++;
      goto timeloop;
    }
  }

  else
  {
    repeated = 0;
    fdyn->dt_prop = fdyn->dta;
  }
}

/*---------------------------------------------- update acceleration ---*/
/*----- copy solution from sol_increment[5][j] to sol_increment[4][j] */
/*--- -> prev. acceleration becomes (n-1)-accel. of next time step ---*/
if (fdyn->step > 1)
   solserv_sol_copy(actfield,0,node_array_sol_increment,
                               node_array_sol_increment,
                               ipos->accn,ipos->accnm);
/*------------------------ evaluate acceleration in this time step ---*
 *-------------------------------- depending on integration method ---*/
if (fdyn->step == 1)
{ /* do just a linear interpolation within the first timestep */
  solserv_sol_zero(actfield,0,node_array_sol_increment,ipos->accnm);
  solserv_sol_add(actfield,0,node_array_sol_increment,
                             node_array_sol_increment,
                             ipos->velnp,ipos->accn, 1.0/fdyn->dta);
  solserv_sol_add(actfield,0,node_array_sol_increment,
                             node_array_sol_increment,
                             ipos->veln,ipos->accn,-1.0/fdyn->dta);
  solserv_sol_copy(actfield,0,node_array_sol_increment,
                              node_array_sol_increment,
                              ipos->accn,ipos->accnm);
}
else
  fluid_acceleration(actfield,ipos,fdyn->iop);

/*------------------------------------------- update time step sizes ---*/
fdyn->dtp = fdyn->dta;
if (fdyn->adaptive)
{
  fdyn->dta = fdyn->dt_prop;
}

/*-------------------------------------------------- steady state check */
if (fdyn->stchk==iststep)
{
  iststep=0;
  steady = fluid_steadycheck(actfield,ipos,numeq_total);
}
else
{
  if (par.myrank==0)
    fprintf(out,"            |            |");
}


/*----------------------------------------------- lift&drag computation */
if (fdyn->liftdrag==ld_stress || fdyn->liftdrag==ld_nodeforce)
{
  container.str = str;
  *action = calc_fluid_liftdrag;
  fluid_liftdrag(1,action,&container,actfield,
                 actsolv,actpart,actintra,ipos);
}


/* stress computation */
if (ioflags.fluid_stress==1)
{
  container.str = str;
  *action = calc_fluid_stress;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
      &container,action);
}


/* error calculation for beltrami and kim-moin */
if (fdyn->init==8 || fdyn->init==9)
{
  container.vel_error  = 0.0;
  container.pre_error  = 0.0;
  container.error_norm = 2;

  *action = calc_fluid_error;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
      &container,action);


  container.vel_error = sqrt(container.vel_error);
  container.pre_error = sqrt(container.pre_error);

  printf("\n  L2_err  :  velocity %f  pressure %f \n\n",
      container.vel_error,container.pre_error);
}





/*------- copy solution from sol_increment[1][j] to sol_increment[0][j] */
/*------- -> prev. solution becomes (n-1)-solution of next time step ---*/
solserv_sol_copy(actfield,0,node_array_sol_increment,
                            node_array_sol_increment,ipos->veln,ipos->velnm);

/* copy solution from sol_increment[velnp][j] to sol_increment[veln][j] */
/*--- -> actual solution becomes previous solution of next time step ---*/
solserv_sol_copy(actfield,0,node_array_sol_increment,
                            node_array_sol_increment,ipos->velnp,ipos->veln);

/* some time could be saved here by swapping flags instead of copying but
   this interferes with restart */                            
 /* leftspace = ipos->velnm;
    ipos->velnm = ipos->veln;*/
 /* ipos->veln = ipos->velnp; */

 /* use remaining space for new solution */
 /* ipos->velnp = leftspace; */

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;
resstep++;
restartstep++;

/*------------------------------------------- write restart to pss file */
if (restartstep==fdyn->uprestart)
{
   restartstep=0;
#ifdef BINIO
   restart_write_bin_fluiddyn(&out_context,fdyn);
#else
   restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);
#endif
}

/*-------- copy solution from sol_increment[ipos->velnp][j] to sol_[actpos][j]
           and transform kinematic to real pressure --------------------*/
solserv_sol_copy(actfield,0,node_array_sol_increment,
                            node_array_sol,ipos->veln,actpos);
fluid_transpres(actfield,0,0,actpos,fdyn->numdf-1,0);

/*-- copy solution on level 2 at (n+1) to place (n) for multi-level FEM */
#if defined(FLUID2_ML) || defined(FLUID3_ML)
if (fdyn->mlfem==1) fluid_smcopy(actpart);
#endif

/*--------------------------------------- write solution to .flavia.res */
if (resstep==fdyn->upres && par.myrank==0 && ioflags.output_gid==1)
{

  if(ioflags.fluid_sol==1)
  {
    out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
    out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }
  if(ioflags.fluid_stress==1)
  {
    out_gid_sol("stress",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }
}

/*--------------------------------------- write solution to binary */
#ifdef BINIO
if (resstep==fdyn->upres && ioflags.output_bin==1)
{

  if (ioflags.fluid_sol==1)
  {
    out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY | OUTPUT_PRESSURE);
  }
  if (ioflags.fluid_stress==1)
  {
    out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_STRESS);
  }
}
#endif

if (resstep==fdyn->upres) {
  resstep=0;
}

/*---------------------------------------------- write solution to .out */
if (outstep==fdyn->upout && ioflags.output_out==1 && ioflags.fluid_sol==1)
{
   outstep=0;
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,genprob.numff,actpos,fdyn->acttime);

tt=ds_cputime()-t2;
tts+=tt;
printf("PROC  %3d | total time for this time step: %10.3e \n",par.myrank,tt);


if (par.myrank==0)
{
  fprintf(out,"            |            |");
  fprintf(out," %10.3E |\n",tt);
  fflush(out);
}

dsmemreport();

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime && steady==0)
{
   goto timeloop;
}
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/

/*======================================================================*
 |                      F I N A L I S I N G                             |
 *======================================================================*/
if (pssstep==0) actpos--;
/*------------------------------------- print out solution to .out file */
if (outstep!=0 && ioflags.output_out==1 && ioflags.fluid_sol==1)
  out_sol(actfield,actpart,actintra,fdyn->step,actpos);

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
printf("PROC  %3d | FIELD FLUID     | total time for time loop           : %10.3E \n",
        par.myrank,tts);
}

/* write total cpu-time to .out */
if (par.myrank==0)
{
  fprintf(out,"\n\n\n");
  fprintf(out," total time element for calculations: %10.3E \n", tes);
  fprintf(out," total time for solver              : %10.3E \n", tss);
  fprintf(out," total time for time loop           : %10.3E \n", tts);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

#ifdef BINIO
destroy_bin_out_field(&out_context);
#endif

/*--------------------------------------------------- cleaning up phase */
amdel(&frhs_a);
if (ioflags.fluid_vis==1 )
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);

/*----------------------------------------------------------------------*/
#ifndef PARALLEL
CCAFREE(actintra);
#endif

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_isi */

#endif
/*! @} (documentation module close)*/
