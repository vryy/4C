/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fluid

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
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
void fluid_isi(FLUID_DYNAMIC *fdyn)
{
INT             itnum;              /* counter for nonlinear iteration  */
INT             i;                  /* simply a counter                 */
INT             numeq;              /* number of equations on this proc */
INT             numeq_total;        /* total number of equations        */
INT             init;               /* flag for solver_control call     */
INT             calstress=0;        /* flag for stress calculation      */
INT             nsysarray=1;        /* one system matrix                */
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
INT             step_s;
DOUBLE          vrat,prat;          /* convergence ratios               */
DOUBLE          t1,t2,ts,te,tt;	    /*					*/
DOUBLE          tes=0.0;            /*					*/
DOUBLE          tss=0.0;            /*					*/
DOUBLE          tts=0.0;            /*					*/
FLUID_STRESS    str;           

SOLVAR         *actsolv;            /* pointer to active sol. structure */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */

ARRAY           ftimerhs_a;
DOUBLE         *ftimerhs;	    /* time - RHS			*/
ARRAY           fiterhs_a;
DOUBLE         *fiterhs;	    /* iteration - RHS  		*/
ARRAY           time_a;             /* stored time                      */

CONTAINER       container;          /* contains variables defined in container.h */

#ifdef DEBUG 
dstrc_enter("fluid_isi");
#endif


/*======================================================================* 
 |                    I N I T I A L I S A T I O N                       |
 *======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
dynvar      = &(fdyn->dynvar);
restart     = genprob.restart;
container.actndis  = 0;
container.turbu    = fdyn->turbu;
container.fieldtyp = actfield->fieldtyp;
str         = str_none;
dynvar->acttime=ZERO;

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

/*--------------------------------------- allocate 1 dist. vector 'rhs' */
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
if (par.myrank==0 && ioflags.fluid_vis_file==1 )
amdef("time",&time_a,1000,1,"DV");

/*--------------------------------------------- initialise fluid field */
if (restart>0)
{
   if (fdyn->init>0)
      dserror("Initial field either by restart, or by function or from file ...\n");
   else
   {
      fdyn->resstep=genprob.restart;
      fdyn->init=2;
   }
}
fluid_init(actpart,actintra,actfield,fdyn,action,&container,4,str);		     
actpos=0;

/*--------------------------------------- init all applied time curves */
for (actcurve=0; actcurve<numcurve; actcurve++)
dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);   

/*-------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield, fdyn);

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
	       	       
/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);
if (ioflags.fluid_sol_gid==1 && par.myrank==0) 
{
   out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->time);
   out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->time);
}


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
t2=ds_cputime();

fdyn->step++;
iststep++;

/*------------------------------------------ check (starting) algorithm */
if (restart!=0)
{
   step_s=fdyn->step;
   fdyn->step=1;
   fluid_startproc(fdyn,&nfrastep);
   restart=0; 
   fdyn->step = step_s;
}
if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(fdyn,&nfrastep);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons(fdyn,dynvar);
fdyn->time += dynvar->dta; 
dynvar->acttime=fdyn->time;

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout(fdyn,dynvar);

/*--------------------- set dirichlet boundary conditions for  timestep */
fluid_setdirich(actfield,fdyn,3);

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

if (fdyn->itchk!=0 && par.myrank==0)
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
fluid_icons(fdyn,dynvar,itnum);

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();
container.dvec         = NULL;
container.ftimerhs     = ftimerhs;
container.fiterhs      = fiterhs;
container.global_numeq = numeq_total;
container.nii          = dynvar->nii;
container.nif          = dynvar->nif;
container.kstep        = 0;
container.is_relax     = 0;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;	     

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
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
fluid_result_incre(actfield,actintra,&(actsolv->sol[0]),3,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,NULL,fdyn);    

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(fdyn,vrat,prat,ZERO,itnum,te,ts);

/*--------------------- check if nonlinear iteration has to be finished */
if (converged==0)
{
   itnum++;
   goto nonlniter;
}
/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/

/*-------------------------------------------------- steady state check */
if (fdyn->stchk==iststep)
{
   iststep=0;
   steady = fluid_steadycheck(fdyn,actfield,numeq_total);
}

/*-------- copy solution from sol_increment[3][j] to sol_increment[1[j] */
solserv_sol_copy(actfield,0,1,1,3,1);
/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;
resstep++;
restartstep++;

/*---------------------------------------------- write solution to .pss */
if (pssstep==fdyn->uppss && ioflags.fluid_vis_file==1 && par.myrank==0)
{
   pssstep=0;   
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = fdyn->time;   
   actpos++;
}

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]   
           and transform kinematic to real pressure --------------------*/
fluid_sol_copy(actfield,0,1,0,3,actpos,fdyn->numdf);

/*--------------------------------------- write solution to .flavia.res */
if (resstep==fdyn->upres &&ioflags.fluid_sol_gid==1 && par.myrank==0) 
{
   resstep=0;
   /*out_checkfilesize(1);*/
   out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->time);
   out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->time);
}

/*---------------------------------------------- write solution to .out */
if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
}

/*------------------------------------------- write restart to pss file */
if (restartstep==fdyn->res_write_evry)
{
   restartstep=0;
   restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);   
}

tt=ds_cputime()-t2;
tts+=tt;
printf("PROC  %3d | total time for this time step: %10.3e \n",par.myrank,tt);

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->time <= fdyn->maxtime && steady==0)
   goto timeloop; 
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/
 
/*======================================================================* 
 |                      F I N A L I S I N G                             |
 *======================================================================*/
if (pssstep==0) actpos--;
/*------------------------------------- print out solution to .out file */
if (outstep!=0 && ioflags.fluid_sol_file==1)
out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*------------------------------------ print out solution to 0.pss file */
if (ioflags.fluid_vis_file==1 && par.myrank==0)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = fdyn->time;   
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
printf("PROC  %3d | FIELD FLUID     | total time for time loop           : %10.3E \n", 
        par.myrank,tts);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&fiterhs_a);
if (par.myrank==0 && ioflags.fluid_vis_file==1 )
amdel(&time_a);
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
