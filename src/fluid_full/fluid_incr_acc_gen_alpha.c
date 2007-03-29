/*!----------------------------------------------------------------------
\file
\brief implicit, one step time integration algorithm for fluid problems

<pre>



Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/gammi/
            +49-(0)89-289-15235
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

#ifdef D_FLUID2_TDS
#include "../fluid2/fluid2.h"
#include "../fluid2_TDS/fluid2_TDS_prototypes.h"


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
\brief implicit, one step time integration algorithm for fluid problems

<pre>                                                         gammi 12/06

The time integration for the fluid is based on the gen-alpha
scheme proposed by

Kenneth E. Jansen, Christian H. Whiting, Gregory M. Hulbert

"A generalized-\alpha method for integrating the filtered Navier-Stokes
 equations with a stabilized finite element method"

(Elsevier Science)

or the dissertation of Christian H. Whiting.

For this algorithm, the nonlinear equation

 +-du                     dsu                                      -+
R| --(n+a_M) , u(n+a_F) , ---(n+a_M) , su(n+a_F) , p(n+1) , p (n+1) | = 0
 +-dt			   dt                                      -+

is (partially) linearised an solved iteratively.


</pre>

\return void


------------------------------------------------------------------------*/

void fluid_incr_acc_gen_alpha(void)
{
INT             itnum;              /* counter for nonlinear iteration  */
INT             i;                  /* simply a counter                 */
INT             numeq;              /* number of equations on this proc */
INT             numeq_total;        /* total number of equations        */
INT             init;               /* flag for solver_control call     */
INT             actsysarray;        /* number of actual sysarray        */
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
DOUBLE          vrat,prat;          /* convergence ratios               */
DOUBLE          t1,t2,ts,te,tt;	    /*					*/
DOUBLE          tes=0.0;            /*					*/
DOUBLE          tss=0.0;            /*					*/
DOUBLE          tts=0.0;            /*					*/
DOUBLE          fact;
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

INT             disnum_calc;
INT             disnum_io;


#ifdef BINIO
BIN_OUT_FIELD   out_context;
BIN_OUT_FIELD   restart_context;
#endif

#ifdef DEBUG
dstrc_enter("fluid_incr_acc_gen_alpha");
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

disnum_calc = disnum_io = 0;

actsysarray = disnum_calc;

/* Zero everything so we can rely on the values we find. */
memset(&container, 0, sizeof(CONTAINER));
container.disnum    = disnum_calc;
container.turbu     = fdyn->turbu;
container.fieldtyp  = actfield->fieldtyp;

ipos = &(actfield->dis[disnum_calc].ipos);

/* set flag for stress evaluation */
str         = str_none;


fdyn->acttime=ZERO;

if (fdyn->adaptive)
{
    dserror("No adaptive time integration for incr. acc. gen-alpha");
}

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
if (actintra->intra_fieldtyp != fluid)
{
    dserror("there is a non fluid proc for a pure fluid problem!");
}

/*------------------------------------ prepare lift&drag calculation ---*/
if (fdyn->liftdrag!=ld_none)
{
    dserror("lift and drag calculation not included in incr_acc_gen_alpha!");
}

if (ioflags.fluid_stress==1)
{
    dserror("no fluid stresses for incr_acc_gen_alpha!");
}

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
    case 8:
      fprintf(out,"Solving for Incremental Accelerations with Generalised Alpha\n");
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
    /*dserror("Initial field either by restart, or by function or from file ...\n");*/
    printf("Restart: Initial field not used!!!\n\n\n");


  fdyn->resstep=genprob.restart;
  fdyn->init=2;
}


/*------------------------------------------ fixed time step size ---*/
ipos->velnp = 0;  /* most recent velocities and pressure             */
ipos->veln  = 1;  /* previous velocities and pressure (at time n)    */
ipos->velnm = 2;  /* intermediate velocities (at time n-alpha_F)     */
ipos->accnp = 3;  /* most recent acceleration (at time n+1)          */
ipos->accn  = 4;  /* previous acceleration (at time n)               */
ipos->accnm = 5;  /* intermediate acceleration (at time n+alpha_M)   */

fluid_init(actpart, actintra, actfield, disnum_calc, disnum_io,
    action, &container, 6,ipos,str);
actpos=0;

/*--------------------------------------- init all applied time curves -*/
for (actcurve=0; actcurve<numcurve; actcurve++)
  dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);

/*-------------------------------------- init the dirichlet-conditions -*/
fluid_initdirich(actfield, disnum_calc, ipos);

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
solver_control(actfield,disnum_calc,actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);

/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,disnum_calc);

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
                   actfield, actpart, actintra, disnum_io);

if (disnum_io != disnum_calc)
  init_bin_out_field(&restart_context,
      &(actsolv->sysarray_typ[actsysarray]),
      &(actsolv->sysarray[actsysarray]),
      actfield, actpart, actintra, disnum_calc);

#endif
  
/*-------------------------------------- print out initial data to .out */
if (ioflags.output_out==1 && ioflags.fluid_sol==1 && par.myrank==0)
  out_sol(actfield,actpart,disnum_io,actintra,fdyn->step,actpos);


/*------------------------------- print out initial data to .flavia.res */
if (ioflags.output_gid==1 && par.myrank==0)
{
  if (ioflags.fluid_sol==1)
  {
    out_gid_sol("velocity",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
    out_gid_sol("pressure",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
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
   monitoring(actfield,disnum_io,genprob.numff,actpos,fdyn->acttime);
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
 *  ipos->velnp ... velocity at time (n+1,i)                            *
 *  ipos->velnm ... velocity at time (n+alpha_F,i)                      *
 *  ipos->veln  ... velocity at time (n)                                *
 *  ipos->accnp ... acceleration at time (n+1,i)                        *
 *  ipos->accnm ... acceleration at time (n+alpha_M,i)                  *
 *  ipos->accn  ... acceleration at time (n)                            *
 *                                                                      *
 *  ipos->velnp ... pressure at time (n+1,i)                            *
 *  ipos->veln  ... pressure at time (n)                                * 
 *======================================================================*/

while (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime && steady==0)
{
 t2=ds_cputime();

 fdyn->step++;
 iststep++;
 
 if (par.myrank==0)
     fprintf(out,"%5d | %10.3E |",fdyn->step,fdyn->acttime);
 
/*------------------------------------------ check (starting) algorithm */
 if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(&nfrastep,0);
 
/*------------------------------ calculate constants for time algorithm */
 fluid_tcons();
  
 fdyn->acttime += fdyn->dta;
 

/*--------- reset counter and convergence criterion for nonlin. iter. */
 itnum    =1;
 converged=0;

/*------------------------------------------------ output to the screen */
 if (par.myrank==0) fluid_algoout();
 
/*------------------------------------- start time step on the screen---*/
 if (fdyn->itnorm!=fncc_no && par.myrank==0)
 {
     printf("----------------------------------------------------------------\n");
     printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -| \n");
 }

/*-------- set dirichlet boundary conditions to sol_increment[velnp] ---*/
 fluid_setdirich(actfield, disnum_calc,ipos,ipos->velnp);

/*-------------------------- estimate new acceleration (and pressure)
                             --- I assume constant velo's and pressure */

 if(actfield->dis[disnum_calc].element[0].e.f2->stab_type == stab_tds)
 {
     f2_estimate_new_trial_values_for_inc_gen_alpha(
	 actpart,
	 actintra,
	 actfield,
	 ipos,
	 disnum_calc);
 }


/*======================================================================*
 |           N O N L I N E A R   I T E R A T I O N                      |
 *======================================================================*/
 while(converged==0)
 {
     /*----------------------- intitialise global matrix and global rhs */
     solserv_zero_vec(&(actsolv->rhs[0]));
     solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
		      &(actsolv->sysarray_typ[actsysarray]));

     /*-------------------------------------- re-initialise neumann bcs */
     inherit_design_dis_neum(&(actfield->dis[disnum_calc]));

     /*-------------------------------------- initialise iterations-rhs */
     amzero(&frhs_a);

#ifdef QUASI_NEWTON
     fdyn->itnum = itnum;
#endif

     /*--------- form incremental matrices, residual and element forces */
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

     /*---------------------------------------------------- add rhs: ---*/
     assemble_vec(actintra,
		  &(actsolv->sysarray_typ[actsysarray]),
		  &(actsolv->sysarray[actsysarray]),
		  &(actsolv->rhs[0]),
		  frhs,
		  1.0
	 );

     
#if 0						
    for(i=0;i<actsolv->rhs[0].numeq_total;i++)
     {
	 printf("rhs[%3d] %12.5e\n",i,actsolv->rhs[0].vec.a.dv[i]);
     }
#endif
     
    /*---------------------------------------------------- solve system */
#ifdef PERF
     perf_begin(80);
#endif


     init=0;
     t1=ds_cputime();

     solver_control(actfield,disnum_calc,actsolv, actintra,
		    &(actsolv->sysarray_typ[actsysarray]),
		    &(actsolv->sysarray[actsysarray]),
		    &(actsolv->sol[0]),
		    &(actsolv->rhs[0]),
		    init);
     ts=ds_cputime()-t1;
     tss+=ts;


#ifdef PERF
     perf_end(80);
#endif


     /*-- set flags for stability parameter evaluation and convergence
                                                                 check */
     fdyn->ishape=0;
     
     /*-- return solution to the nodes and calculate the convergence
                                                                ratios */
     fluid_result_incre_for_genalpha(
	                actfield, disnum_calc, actintra,
			&(actsolv->sol[0]),&(actsolv->rhs[0]),
			ipos,
			&(actsolv->sysarray[actsysarray]),
			&(actsolv->sysarray_typ[actsysarray]),
			&vrat,&prat,NULL);

     /*---------------- update of estimate for time dependent subscales */
     if(actfield->dis[disnum_calc].element[0].e.f2->stab_type == stab_tds)
     {
	 f2_update_subscale_pres_for_inc_gen_alpha(
	     actpart,
	     actintra,
	     actfield,
	     ipos,
	     disnum_calc);
	 
	 f2_update_subscale_vel_for_inc_gen_alpha(
	     actpart,
	     actintra,
	     actfield,
	     ipos,
	     disnum_calc);
     }

     
     /*----------------------------------- iteration convergence check */
     converged = fluid_convcheck(vrat,prat,ZERO,itnum,te,ts);

     itnum++;

 }
 /*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/

 f2_time_update_subscales_for_incr_gen_alpha (
     actpart,
     actintra,
     actfield,
     ipos,
     disnum_calc);

/*-------------------------------------------------- steady state check */
 if (fdyn->stchk==iststep)
 {
     iststep=0;
     steady = fluid_steadycheck(actfield, disnum_calc, ipos,numeq_total);
 }
 else
 {
     if (par.myrank==0)
	 fprintf(out,"            |            |");
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

/* copy solution from sol_increment[velnp][j] to sol_increment[veln][j] */
/*--- -> actual solution becomes previous solution of next time step ---*/
 solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,
		  node_array_sol_increment,ipos->velnp,ipos->veln);

 /* copy solution from sol_increment[accnp][j] to sol_increment[accn][j] */
/*--- -> actual solution becomes previous solution of next time step ---*/
 solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,
		  node_array_sol_increment,ipos->accnp,ipos->accn);


/* copy solution from sol_increment[ipos->velnp][j] to sol_[actpos][j]
   and transform kinematic to real pressure */
 solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,
                            node_array_sol,ipos->velnp,actpos);
 fluid_transpres(actfield,disnum_calc,0,actpos,fdyn->numdf-1,0);


/* finalise this timestep */
 outstep++;
 pssstep++;
 resstep++;
 restartstep++;
 
 
/* write restart to pss file */
 if (restartstep==fdyn->uprestart)
 {
     restartstep=0;
#ifdef BINIO
     if (disnum_io != disnum_calc)
	 restart_write_bin_fluiddyn(&restart_context,fdyn);
     else
	 restart_write_bin_fluiddyn(&out_context,fdyn);
#else
     restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);
#endif
 }
 
 
/*--------------------------------------- write solution to .flavia.res */
 if (resstep==fdyn->upres && par.myrank==0 && ioflags.output_gid==1)
 {
     
     if(ioflags.fluid_sol==1)
     {
	 out_gid_sol("velocity",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
	 out_gid_sol("pressure",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
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
}
#endif

 if (resstep==fdyn->upres) {
     resstep=0;
 }

/*---------------------------------------------- write solution to .out */
 if (outstep==fdyn->upout && ioflags.output_out==1 && ioflags.fluid_sol==1)
 {
    outstep=0;
    out_sol(actfield,actpart, disnum_io, actintra,fdyn->step,actpos);
 }
 
/*---------------------------------------------------------- monitoring */
 if (ioflags.monitor==1)
     monitoring(actfield, disnum_io, genprob.numff,actpos,fdyn->acttime);
 
 tt=ds_cputime()-t2;
 tts+=tt;
 
 if (par.myrank==0)
 {
     printf("total time for this time step: %10.3e \n",tt);
     
     fprintf(out,"            |            |");
     fprintf(out," %10.3E |\n",tt);
     fflush(out);
 }


 dsmemreport();

}
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/

 dsmemreport();

#ifdef BINIO
destroy_bin_out_field(&out_context);
if (disnum_io != disnum_calc)
  destroy_bin_out_field(&restart_context);
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
} /* end of fluid_incr_acc_gen_alpha */

#endif /*D_FLUID2_TDS*/


#endif /*D_FLUID*/
/*! @} (documentation module close)*/
