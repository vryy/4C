#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
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
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern int            numcurve;
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
 
void fluid_isi(FLUID_DYNAMIC *fdyn)
{
int             itnum;              /* counter for nonlinear iteration          */
int             i;                  /* simply a counter */
int             numeq;              /* number of equations on this proc */
int             numeq_total;        /* total number of equations */
int             init;               /* flag for solver_control call */
int             nsysarray=1;        /* at the moment there's only one system matrix */
int             actsysarray=0;      /* number of actual sysarray */
int             istep=0;            /* counter for time integration */
int             iststep=0;          /* counter for time integration */
int             nfrastep;           /* number of steps for fractional-step-theta procedure */
int             actcurve;           /* actual timecurve */
int             converged=0;        /* convergence flag */
int             steady=0;           /* flag for steady state */
int             actpos;             /* actual position in solution history */
int             diff;
int             max;
double          vrat,prat;          /* convergence ratios */
double          t,t0,t1,ts,te;	    /* variables for time tracing */
double          tfs=0.0;
double          tes=0.0;
double          tas=0.0;
double          tss=0.0;
double          tds=0.0;
double          tcs=0.0;
double          tsts=0.0;

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structure cal_action enum */
FLUID_DYN_CALC *dynvar;             /* pointer to the structure fluid_dyn_calc */

ARRAY           ftimerhs_a;         /* time - RHS */
double         *ftimerhs;
ARRAY           fiterhs_a;          /* iteration - RHS */
double         *fiterhs;
ARRAY           time_a;             /* stored time */

#ifdef DEBUG 
dstrc_enter("fluid_isi");
#endif

/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
dynvar      = &(fdyn->dynvar);

#ifdef PARALLEL 
actintra    = &(par.intra[0]);
/* if we are not parallel, we have to allocate an alibi intra-communicator structure */
#else
actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = fluid;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the fluid */
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all) */
if (actintra->intra_fieldtyp != fluid) goto end;

/*--------------- stiff_array already exists, so copy the mask of it to */
/* reallocate the vector of sparse matrice (length: nsyssarray=1) */
actsolv->sysarray_typ = 
(SPARSE_TYP*)REALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

actsolv->sysarray = 
(SPARSE_ARRAY*)REALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

/*------------------------------- init the dist sparse matrices to zero */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(actsolv->sysarray[actsysarray],
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
printf("number of eqations on PROC %3d : %10d \n", 
        par.myrank,numeq);
if (par.myrank==0)
   printf("total number of equations: %10d \n",numeq_total);	
#else
printf("total number of equations: %10d \n",numeq_total);
#endif

/*---------------------------------------allocate 1 dist. vector 'rhs' */
actsolv->nrhs = 1;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->rhs[0]));

/*------------------------------------------- allocate 1 dist. solution */		       
actsolv->nsol= 1;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->sol[0]));
                                     
/*--------------- allocate one redundant vector ftimerhs of full lenght */
/*        this is used by the element routines to assemble the  Time RHS*/
ftimerhs = amdef("ftimerhs",&ftimerhs_a,numeq_total,1,"DV");

/*---------------  allocate one redundant vector fiterhs of full lenght */
/*   this is used by the element routines to assemble the  Iteration RHS*/
fiterhs = amdef("fiterhs",&fiterhs_a,numeq_total,1,"DV");

/*--------------------------- allocate one vector for storing the time */
max = fdyn->nstep/fdyn->idisp;
amdef("time",&time_a,max,1,"DV");
/*--------------------------------------------- initialise fluid field */
fluid_init(actfield, fdyn);		     
actpos=0;

/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);   

/*-------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield, fdyn);

/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to 
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init phase
      is cheap in computation (can be costly in memory)
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);
	       
/*----------------- init the assembly for stiffness and for mass matrix */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray);
	       	       
/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit_fluid(actfield,actpart,action);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*----------------------------------------------------------------------* 
 | do the time loop                                                     |
 *----------------------------------------------------------------------*/
timeloop:
fdyn->step++;
iststep++;
/*------------------------------------------ check (starting) algorithm */
if (fdyn->step<=(fdyn->numfss+1))
   fluid_startproc(fdyn,&nfrastep);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons(fdyn,dynvar);
fdyn->time += dynvar->dta; 
/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout(fdyn,dynvar);
/*------------------------------ set dirichlet boundary conditions for 
                                                          this timestep */
t1=ds_cputime();
fluid_setdirich(actfield,fdyn);
tds+=ds_cputime()-t1;

/*------------- get maximum velocity for stability parameter definition */
/* fluid_maxvel(fdyn,actfield,&(dynvar->velmax)); 
   not necessary at the moment !!!!!! ----------------------------------*/

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

/*----------------------------------------------------------------------*
 | do the nonlinear iteration                                           |
 *----------------------------------------------------------------------*/
if (fdyn->itchk!=0 && par.myrank==0)
{
   printf(" _____________________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error-| \n");
}
itnum=1;
nonlniter:

/*------------------------- calculate constants for nonlinear iteration */
fluid_icons(fdyn,dynvar,itnum);

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));

/*--------------------------------------- initialise and iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();
calelm_fluid(actfield,actsolv,actpart,actintra,actsysarray,-1,
             ftimerhs,fiterhs,numeq_total,dynvar->nii,
	     dynvar->nif,0,action);
te=ds_cputime()-t1;
tes+=te;	     
/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
t1=ds_cputime();
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[0]),
             ftimerhs,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[0]),
             fiterhs,
             1.0
             );
tas+=ds_cputime()-t1;
/*-------------------------------------------------------- solve system */
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
/*-- set flags for stability parameter evaluation and convergence check */
dynvar->ishape=0;
/*--- return solution to the nodes and calculate the convergence ratios */
t1=ds_cputime();
fluid_result_incre(actfield,actintra,&(actsolv->sol[0]),3,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,fdyn);
tcs+=ds_cputime()-t1;		     
/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(fdyn,vrat,prat,itnum,te,ts);

/*--------------------- check if nonlinear iteration has to be finished */
if (converged==0)
{
   itnum++;
   goto nonlniter;
}
/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/
if (par.myrank==0)
{
   printf("|____________|___________________|______________|_____________| \n");
   printf("\n"); 
}   
/*-------------------------------------------------- steady state check */
t1=ds_cputime();
if (fdyn->stchk==iststep)
{
   iststep=0;
   steady = fluid_steadycheck(fdyn,actfield,numeq_total);
}
tsts+=ds_cputime()-t1;
/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
fluid_copysol(fdyn,actfield,3,1,0);
/*---------------------------------------------- finalise this timestep */
istep++;
t1=ds_cputime();
/*------------------------------------------ store results in node->sol *
 * on position actpos; actpos=0: initial data;                          *
 *                     actpos=1... stored timesteps --------------------*/
if (istep==fdyn->idisp || steady>0)
{
   istep=0;
   actpos++;
   fluid_copysol(fdyn,actfield,3,actpos,1);
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   {
      diff = actpos - time_a.fdim;
      max  = IMAX(diff,5);
      amredef(&(time_a),time_a.fdim+max,1,"DV");
   }
   time_a.a.dv[actpos] = fdyn->time;
   /*------------------------------------------------ print out to .out */
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
} 
tfs+=ds_cputime()-t1;  
/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->time <= fdyn->maxtime && steady==0)
   goto timeloop;

/*----------------------------------------------------------------------*/
#ifdef PARALLEL
MPI_Barrier(MPI_COMM_WORLD);
#endif
for (i=0;i<par.nprocs;i++)
{
if (par.myrank==i)
{
printf("\n");
printf("PROC %3d: TOTAL TIME for setting dirichlet values: %10.3#E \n", par.myrank,tds);
printf("PROC %3d: TOTAL TIME for element calculations: %10.3#E \n", par.myrank,tes);
printf("PROC %3d: TOTAL TIME for assembly of RHS: %10.3#E \n", par.myrank,tas);
printf("PROC %3d: TOTAL TIME for SOLVER: %10.3#E \n", par.myrank,tss);
printf("PROC %3d: TOTAL TIME for calculation of convergence ratios: %10.3#E \n", par.myrank,tcs);
printf("PROC %3d: TOTAL TIME for steady state check: %10.3#E \n", par.myrank,tsts);
printf("PROC %3d: TOTAL TIME finalising of time step: %10.3#E \n", par.myrank,tfs);
}
}
end:
/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&fiterhs_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
} /* end of fluid_isi */ 

