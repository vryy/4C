#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "ls_prototypes.h"
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
 | routine to control implicit and semi-implicit algorithms for levelset|
 | problems combined with Newton and fixed point iteration schemes.     | 
 | IN PROGRESS: ONE-STEP-THETA                                          |
 |              fixed point like iteration                              |
 |              only Euler no ALE                                       |
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
extern struct _LS_UPDATE lsupdt;



void ls_isi()
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
INT             pssstep=0;	    		/* counter for output control			*/
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
DOUBLE          lrat;		         /* convergence ratios               */
DOUBLE          t1,ts,te;	         /*												*/
DOUBLE          tes=0.0;            /*												*/
DOUBLE          tss=0.0;            /*												*/

SOLVAR         *actsolv;            /* pointer to active sol. structure */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
FRONTLSFLAG     frontlsflag;
 
ARRAY           ftimerhs_a;
DOUBLE         *ftimerhs;	         /* time - RHS							   */
ARRAY           fiterhs_a;
DOUBLE         *fiterhs;	         /* iteration - RHS  						*/

CONTAINER       container;          /* variables defined in container.h */

LS_DYNAMIC     *lsdyn;              /* pointer to ls_dynamic	  		   */

#ifdef DEBUG 
dstrc_enter("ls_isi");
#endif


/*======================================================================* 
 |                    I N I T I A L I S A T I O N                       |
 *======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actfield    = &(field[genprob.numls]);
actsolv     = &(solv[genprob.numls]);
actpart     = &(partition[genprob.numls]);
action      = &(calc_action[genprob.numls]);
lsdyn       = alldyn[genprob.numls].lsdyn;
restart     = genprob.restart;
container.actndis  = 0;
container.fieldtyp = actfield->fieldtyp;
lsdyn->t    = ZERO;

/*---------------- if we are not parallel, we have to allocate an alibi * 
  ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL 
actintra    = &(par.intra[genprob.numls]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = levelset;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the levelset -----*/
/* intracommunicator (in case of nonlinear levelset dyn., this should be all)*/
if (actintra->intra_fieldtyp != levelset) goto end;

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
printf("PROC  %3d | FIELD LEVELSET     | number of equations      : %10d \n", 
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD LEVELSET     | total number of equations: %10d \n",numeq_total);
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

/*------------------------------------------ initialise level set field */
ls_init(actfield,lsdyn,2);
actpos=0;

/*----------------------------------------- initialize level set update */
frontlsflag = front_ls_set;
ls_update(frontlsflag);

/* --------------- construct initial front and write it to output files */
frontlsflag = front_ls_init;
ls_update(frontlsflag);

/* -------------------------------------- write initial front into file */
frontlsflag = front_ls_write;
ls_update(frontlsflag);

/*---------------------------------------- init all applied time curves */
for (actcurve=0; actcurve<numcurve; actcurve++)
dyn_init_curve(actcurve,lsdyn->nitn,lsdyn->dt,lsdyn->mtim);

/*--------------------------------------- init the dirichlet-conditions */
ls_initdirich(actfield,lsdyn);

/*----------------------------------- initialize solver on all matrices */
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
*action = calc_ls_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,lsdyn->step,actpos);
if (ioflags.lset_sol_gid==1 && par.myrank==0) 
{
   out_gid_sol("levelset"    ,actfield,actintra,lsdyn->step,actpos,lsdyn->t);
   out_gid_sol("zerolevelset",actfield,actintra,lsdyn->step,actpos,lsdyn->t);	
}


/*======================================================================* 
 |                         T I M E L O O P                              |
 *======================================================================*/
/* nodal solution history levelset field:                               *
 * sol[0][j]           ... initial data                     				*
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[0][j] ... solution at time (n) 								*
 * sol_increment[1][j] ... solution at time (n+1)								*
 *======================================================================*/
timeloop:
lsdyn->step++;
iststep++;
lsdyn->t += lsdyn->dt;

/*------------------------------------------------ output to the screen */
if (par.myrank==0) ls_algout(lsdyn);

/*--------------------- set dirichlet boundary conditions for  timestep */
ls_setdirich(actfield,lsdyn,1);

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

if (lsdyn->ichk!=0 && par.myrank==0)
{
   printf("|------------------------------------------------|\n");
   printf("|- step/max -|-  tol     [norm]  -|- lset error -|\n");
}
itnum=1;
/*======================================================================* 
 |           N O N L I N E A R   I T E R A T I O N                      |
 *======================================================================*/
nonlniter:

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));
					  
/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_ls;
t1=ds_cputime();
container.dvec         = NULL;
container.ftimerhs     = ftimerhs;
container.fiterhs      = fiterhs;
container.global_numeq = numeq_total;
container.nii          = lsdyn->nii;
container.nif          = lsdyn->nif;
container.kstep        = 0;
container.is_relax     = 0;
calelm(actfield,actsolv,actpart,actintra,actsysarray,actsysarray,
       &container,action);
lsdyn->nif = 0;
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

/*--- return solution to the nodes and calculate the convergence ratios */
ls_result_incre(actfield,actintra,&(actsolv->sol[0]),1,
                &(actsolv->sysarray[actsysarray]),
                &(actsolv->sysarray_typ[actsysarray]),
	             &lrat,lsdyn);

/*----------------------------------------- iteration convergence check */
converged = ls_convcheck(lsdyn,lrat,itnum,te,ts);


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
if (lsdyn->schk==iststep)
{
   iststep=0;
   steady = ls_steadycheck(lsdyn,actfield,numeq_total);
}

/*----------------------------------------------- update zero level set */
frontlsflag = front_ls_updt;
ls_update(frontlsflag);

/* ---------------------------------------------- write front into file */
frontlsflag = front_ls_write;
ls_update(frontlsflag);

/*------- copy solution from sol_increment[1][j] to sol_increment[0][j] */
solserv_sol_copy(actfield,0,1,1,1,0);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;
resstep++;
restartstep++;

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j] ---*/
solserv_sol_copy(actfield,0,1,0,1,0);

/*--------------------------------------- write solution to .flavia.res */
if (ioflags.lset_sol_gid==1 && par.myrank==0) 
{
   out_gid_sol("levelset",actfield,actintra,lsdyn->step,actpos,lsdyn->t);
   out_gid_sol("zerolevelset",actfield,actintra,lsdyn->step,actpos,lsdyn->t);	
}

/*---------------------------------------------- write solution to .out */
if (ioflags.lset_sol_file==1)
{
   out_sol(actfield,actpart,actintra,lsdyn->step,actpos);
}

/*--------------------- check time and number of steps and steady state */
if (lsdyn->step < lsdyn->nitn && lsdyn->t <= lsdyn->mtim && steady==0)
{
	lsdyn->nif = 1;			 
   goto timeloop; 
}
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/

/*======================================================================* 
 |                      F I N A L I S I N G                             |
 *======================================================================*/

/*---------------------------------- print total CPU-time to the screen */
/* ---------------------------------------------- write front into file */
frontlsflag = front_ls_finalize;
ls_update(frontlsflag);

/*---------------------------------- print total CPU-time to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
{
if (par.myrank==i)
{
printf("\n");
printf("PROC  %3d | FIELD LEVELSET     | total time element for calculations: %10.3E \n", 
        par.myrank,tes);
printf("PROC  %3d | FIELD LEVELSET     | total time for solver              : %10.3E \n", 
        par.myrank,tss);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&fiterhs_a);
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
} /* end of ls_isi */ 

#endif
