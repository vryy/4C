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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
static INT     rans=0;             /* flag, if vel. discr. is used     */
static INT     kapomega=1;           /* flag, if pres. discr. is used    */
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
void fluid_isi_tu_1(FLUID_DYNAMIC *fdyn)
{
INT             itnum;              /* counter for nonlinear iteration  */
INT             itnum1;             /* counter for nonlinear iteration  */
INT             itnum2;             /* counter for nonlinear iteration  */
INT             itnumke;            /* counter for nonlinear iteration  */
INT             itnum_n=1;          /* counter for nonlinear iteration  */
INT             itnum_check;        /* counter for nonlinear iteration  */
INT             kapomegastep;       /* counter for timeloop             */
INT             conv_check_rans;
INT             i,kk;               /* simply a counter                 */
INT             numeq[2];           /* number of equations on this proc */
INT             numeq_total[2];     /* total number of equations        */
INT             numeq_oll;
INT             numeq_total_oll;
INT             init;               /* flag for solver_control call     */
INT             start;               /* flag for solver_control call     */
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
DOUBLE          vrat,prat;            /* convergence ratios               */
DOUBLE          kapomegarat,lenghtrat;/* convergence ratios               */
DOUBLE          lower_limit_kappa;    /* convergence ratio                */
DOUBLE          lower_limit_omega;    /* convergence ratio                */

DOUBLE          t1,ts,te;	     /*					*/
DOUBLE          tes=0.0;             /*					*/
DOUBLE          tss=0.0;             /*					*/
FLUID_STRESS    str;           

DIST_VECTOR    *rhs_ko;              /* distr. RHS for solving ko        */
DIST_VECTOR    *sol_ko;              /* distr. ko solution               */

SOLVAR         *actsolv;            /* pointer to active sol. structure */
SOLVAR         *kosolv;             /* solver for kappa-epsilon         */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */

ARRAY           ftimerhs_a;
DOUBLE         *ftimerhs;	    /* time - RHS			*/
ARRAY           ftimerhs_ko_a;
DOUBLE         *ftimerhs_ko;	    /* time - RHS			*/
ARRAY           ftimerhs_kappa_a;
DOUBLE         *ftimerhs_kappa;    /* time - RHS			*/
ARRAY           ftimerhs_omega_a;
DOUBLE         *ftimerhs_omega;    /* time - RHS			*/
ARRAY           ftimerhs_pro_a;
DOUBLE         *ftimerhs_pro;      /* time - RHS			*/
ARRAY           ftimerhs_pro_kappa_a;
DOUBLE         *ftimerhs_pro_kappa; /* time - RHS			*/
ARRAY           ftimerhs_pro_omega_a;
DOUBLE         *ftimerhs_pro_omega; /* time - RHS			*/

ARRAY           fiterhs_a;
DOUBLE         *fiterhs;	    /* iteration - RHS  		*/
ARRAY           fiterhs_ko_a;
DOUBLE         *fiterhs_ko;	    /* iteration - RHS  		*/

ARRAY           time_a;          /* stored time                      */

CONTAINER       container;      /* contains variables defined in container.h */

INT             initialisation=0;
INT             kapomega_yeah=0;
#ifdef DEBUG 
dstrc_enter("fluid_isi_tu_1");
#endif


/*======================================================================* 
 |                    I N I T I A L I S A T I O N                       |
 *======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actfield           = &(field[0]);
actsolv            = &(solv[0]);
actpart            = &(partition[0]);
action             = &(calc_action[0]);
dynvar             = &(fdyn->dynvar);
container.fieldtyp = actfield->fieldtyp;
container.turbu    = fdyn->turbu;
dynvar->dis_capt   = fdyn->dis_capt;
dynvar->coord_scale[0] = fdyn->coord_scale[0];
dynvar->coord_scale[1] = fdyn->coord_scale[1];
dynvar->washvel    = 1.0;
str                = str_none;
dynvar->acttime    = ZERO;

/*-------- this is no relaxation parameter solution and no gen alpha ---*/
container.is_relax = 0;
dynvar->gen_alpha = 0;

/*--------------------------------------------- set max. iterationsteps */
fdyn->itemax_ke    = 5*fdyn->itemax;

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
kosolv  = &(solv[kapomega]);
kosolv->fieldtyp = fluid;

/*-------------------------------------- now we set some default values */
#ifdef PARALLEL
#ifndef SPOOLES_PACKAGE
dserror("SPOOLES package is not compiled in");
#else
kosolv->solvertyp = SPOOLES_nonsym;
#endif
#else
kosolv->solvertyp = umfpack;
#endif
kosolv->parttyp = cut_elements;
kosolv->matrixtyp = oll_matrix;

/* -------------------------- create solver for pressure discretisation */
kosolv->nsysarray = 1;
kosolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(kosolv->nsysarray,sizeof(SPARSE_TYP));
kosolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(kosolv->nsysarray,sizeof(SPARSE_ARRAY));

kosolv->sysarray_typ[0] = oll;
kosolv->sysarray[0].oll = (OLL*)CCACALLOC(1,sizeof(OLL));


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
solserv_getmatdims(actsolv->sysarray[actsysarray],
                   actsolv->sysarray_typ[actsysarray],
                   &numeq[kk],
                   &numeq_total[kk]);

} /* end of loop over sys_arrays */

/*-------------------------------------- initialise solver for pressure */
numeq_total_oll = actfield->dis[kapomega].numeq;
oll_numeq(actfield, actpart, actintra, kapomega, &numeq_oll);

oll_open(kosolv->sysarray[0].oll, numeq_oll, numeq_total_oll, 
	 actfield, actpart, actintra, kapomega);

numeq[kapomega] = numeq_oll;
numeq_total[kapomega] = numeq_total_oll;

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
printf("PROC  %3d | FIELD FLUID     | number of rans equations        : %10d \n", 
        par.myrank,numeq[rans]);
printf("PROC  %3d | FIELD FLUID     | number of kapome equations      : %10d \n", 
        par.myrank,numeq[kapomega]);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD FLUID     | total number of rans equations  : %10d \n",numeq_total[rans]);
printf("          | FIELD FLUID     | total number of kapome equations: %10d \n",numeq_total[kapomega]);
if (par.myrank==0) printf("\n\n");


/*--------------------------------------- allocate 1 dist. vector 'rhs' */
solserv_create_vec(&(actsolv->rhs),1,numeq_total[rans],numeq[rans],"DV");
solserv_zero_vec(&(actsolv->rhs[0]));

/*------------------------------------------- allocate 1 dist. solution */		       
solserv_create_vec(&(actsolv->sol),1,numeq_total[rans],numeq[rans],"DV");
solserv_zero_vec(&(actsolv->sol[0]));

/*-------------------------- allocate 1 dist. vector 'rhs' for ke       */
solserv_create_vec(&rhs_ko,1,numeq_total[kapomega],numeq[kapomega],"DV");
solserv_zero_vec(rhs_ko);
                                     
/*-------------------------- allocate 1 dist. vector 'sol' for ke       */
solserv_create_vec(&sol_ko,1,numeq_total[kapomega],numeq[kapomega],"DV");
solserv_zero_vec(sol_ko);

/*--------------- allocate one redundant vector ftimerhs of full lenght */
/*        this is used by the element routines to assemble the  Time RHS*/
ftimerhs         = amdef("ftimerhs",   &ftimerhs_a,numeq_total[rans],1,"DV");
ftimerhs_ko      = amdef("ftimerhs_ko",&ftimerhs_ko_a,numeq_total[kapomega],1,"DV");
ftimerhs_kappa   = amdef("ftimerhs_kappa",&ftimerhs_kappa_a,numeq_total[kapomega],1,"DV");
ftimerhs_omega   = amdef("ftimerhs_omega",&ftimerhs_omega_a,numeq_total[kapomega],1,"DV");
ftimerhs_pro     = amdef("ftimerhs_pro",&ftimerhs_pro_a,numeq_total[kapomega],1,"DV");
ftimerhs_pro_kappa   = amdef("ftimerhs_pro_kappa",&ftimerhs_pro_kappa_a,numeq_total[kapomega],1,"DV");
ftimerhs_pro_omega   = amdef("ftimerhs_pro_omega",&ftimerhs_pro_omega_a,numeq_total[kapomega],1,"DV");

/*---------------  allocate one redundant vector fiterhs of full lenght */
/*   this is used by the element routines to assemble the  Iteration RHS*/
fiterhs    = amdef("fiterhs",&fiterhs_a,numeq_total[rans],1,"DV");
fiterhs_ko = amdef("fiterhs_ko",&fiterhs_ko_a,numeq_total[kapomega],1,"DV");

/*--------------------------- allocate one vector for storing the time */
amdef("time",&time_a,1000,1,"DV");

/*--------------------------------------------- initialise fluid field */
fluid_init(actpart,actintra,actfield,fdyn,action,&container,4,str);		     
fluid_init_tu(actfield, fdyn);	
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
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);

solver_control(kosolv, actintra,
               &(kosolv->sysarray_typ[0]),
               &(kosolv->sysarray[0]),
               sol_ko,
               rhs_ko,
               init);
	       
/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,k_array,rans);
init_assembly(actpart,kosolv,actintra,actfield,0,kapomega);
	       	       
/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);
/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);
actpos++;

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
if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(fdyn,&nfrastep,0);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons(fdyn,dynvar);
if (kapomega_yeah == 1) fluid_tcons_tu(fdyn,dynvar); 

fdyn->time += dynvar->dta;
dynvar->acttime=fdyn->time;

/*------ flag and counter for convergence check for kappa-omega model*/
conv_check_rans = 0;
itnum_check = 1;

nonlniter_check:
container.actndis=0;

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
container.nii          = dynvar->nii;
container.nif          = dynvar->nif;
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
dynvar->ishape=0;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre(actfield,actintra,&(actsolv->sol[k_array]),3,
                     &(actsolv->sysarray[k_array]),
                     &(actsolv->sysarray_typ[k_array]),
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
 | -->  end of nonlinear iteration   for rans                           |
 *----------------------------------------------------------------------*/

/*======================================================================* 
                        K A P P A - O M E G A 
   if possible: calc. until steady state before kappa, omega and 
                     eddy-viscosity are calcul.                         
 *======================================================================*/
if (kapomega_yeah == 1)   
{
container.actndis=1;

/*-------------- set dirichlet boundary conditions const. for timestep */
fluid_setdirich_tu_1(actfield,fdyn,&lower_limit_kappa,&lower_limit_omega);

/*---------------- initialise the turbulence variables for calculation */
if (initialisation == 0) 
{
fdyn->stepke++;
fluid_set_check_tu_1(actfield,fdyn,lower_limit_kappa,lower_limit_omega);       
/*-------------------------------- amzero timerhs for kappa end epsilon */
amzero(&ftimerhs_kappa_a);
amzero(&ftimerhs_omega_a);

initialisation = 1;
}

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout_tu(fdyn,dynvar);

/*-------------------- counter for charact. lenght and production terms */
itnum2 = 1;
itnumke= 1;

/*-------------------------------------------- initialise productionrhs */
amzero(&ftimerhs_pro_kappa_a);
amzero(&ftimerhs_pro_omega_a);

/*======================================================================* 
 / do the nonlinear iteration  for turbulence  (kappa)                  |
 *======================================================================*/
kappa:

dynvar->kapomega_flag=0;
itnum1=1;

if (fdyn->itchk!=0 && par.myrank==0)
{
   printf(" ______________________________________________________________ \n");
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- kappa  error -| \n");
}

if(dynvar->kapomega_flag==1)
{
/*======================================================================* 
 | do the nonlinear iteration  for turbulence  (omega)                  |
 *======================================================================*/
omega:

dynvar->kapomega_flag=1;
itnum1=1;
/*----------------------------------------------------------------------*/

if (fdyn->itchk!=0 && par.myrank==0)
{
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- omega  error -| \n");
}

} /* end dynvar->kapeps_flag==1 */
 
/*---------------------------------------------------------------------*/
nonlniter1:
fluid_icons_tu(fdyn,dynvar,itnum1,itnumke,itnum_n);
/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(rhs_ko);
solserv_zero_mat(actintra,&(kosolv->sysarray[0]),
                 &(kosolv->sysarray_typ[0]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_ko_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();

if(dynvar->kapomega_flag==0)
{
 container.ftimerhs     = ftimerhs_kappa;
 container.ftimerhs_pro = ftimerhs_pro_kappa;
}
if(dynvar->kapomega_flag==1)
{
 container.ftimerhs     = ftimerhs_omega;
 container.ftimerhs_pro = ftimerhs_pro_omega;
}
container.fiterhs      = fiterhs_ko;
container.global_numeq = numeq_total[kapomega];
container.niturbu_pro  = dynvar->niturbu_pro; 
container.niturbu_n    = dynvar->niturbu_n; 
container.nii          = 0;
container.nim          = 0;
container.nif          = 0;
container.kstep        = 0;
calelm(actfield,kosolv,actpart,actintra,0,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;	     

if(dynvar->kapomega_flag==0) 
{
 ftimerhs_ko = ftimerhs_kappa;
 ftimerhs_pro= ftimerhs_pro_kappa;
} 
if(dynvar->kapomega_flag==1) 
{
 ftimerhs_ko = ftimerhs_omega;
 ftimerhs_pro= ftimerhs_pro_omega;
} 
/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
assemble_vec(actintra,
             &(kosolv->sysarray_typ[0]),
             &(kosolv->sysarray[0]),
             rhs_ko,
             ftimerhs_ko, 
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(kosolv->sysarray_typ[0]),
             &(kosolv->sysarray[0]),
             rhs_ko,
             fiterhs_ko,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(kosolv->sysarray_typ[0]),
             &(kosolv->sysarray[0]),
             rhs_ko,
             ftimerhs_pro,
             1.0
             );

/*-------------------------------------------------------- solve system */
init=0;
t1=ds_cputime();
solver_control(kosolv, actintra,
               &(kosolv->sysarray_typ[0]),
               &(kosolv->sysarray[0]),
               sol_ko,
               rhs_ko,
               init);
ts=ds_cputime()-t1;
tss+=ts;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre_tu_1(actfield,actintra,sol_ko,3,
                     &(kosolv->sysarray[0]),
                     &(kosolv->sysarray_typ[0]),
		         &kapomegarat,fdyn,lower_limit_kappa,lower_limit_omega);

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck_tu(fdyn,kapomegarat,itnum1,te,ts);

/*---------- check if nonlinear iteration for kapeps has to be finished */
if (converged==0)
{
   itnum1++;
   goto nonlniter1;
}

if (dynvar->kapomega_flag==0)
{
   fluid_eddy_update_1(actfield,sol_ko);
   goto omega;
}

itnumke++;
itnum_n++;

/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration  for turbulence (kappa || epsilon)   |
 *----------------------------------------------------------------------*/
if (fdyn->itchk!=0 && par.myrank==0)
{
   printf("|____________|___________________|________________| \n");
   printf("\n"); 
}   

/*--------------- iteration convergence check for characteristic lenght */
if (fdyn->itchk!=0 && par.myrank==0)
{
   printf(" __________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- lenght error -| \n");
}

fluid_eddy_update_1(actfield,sol_ko);
fluid_lenght_update_1(actfield,sol_ko,&lenghtrat,fdyn);
 
/*--------------------- check if nonlinear iteration has to be finished */
converged = fluid_convcheck_tu(fdyn,lenghtrat,itnum2,te,ts);

if (converged==0)
{
   itnum2++;
   goto kappa;
}
  
/*----------------- write eddy (n) to eddy (n+g) due to production-term */
fluid_eddy_pro(actfield);

if (fdyn->itchk!=0 && par.myrank==0)
{
   printf("|____________|___________________|________________| \n");
   printf(" ______________________________________________________________ \n");
   printf("\n"); 
}   

/*----------------------- check if rans has finished due to kappa-omega */
if (conv_check_rans == 1)
{
 converged = fluid_convcheck_test(fdyn,actfield,itnum_check);

 if (converged!=0) 
 {
/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
  fluid_copysol_tu(fdyn,actfield,3,1,0);

/*------------------------------------------- counter for timerhs terms */
  itnum_n = 1;
  fdyn->stepke++;

/*-------------------------------- amzero timerhs for kappa end epsilon */
  amzero(&ftimerhs_kappa_a);
  amzero(&ftimerhs_omega_a);

/*---------------------------------------- end outer-iteration for RANS */
  goto endrans;
 }

 itnum_check++;
}
/*------------------------------ copy solution from rans for conv-check *
                                 for rans due to kappa-omega            */
fluid_copysol_test(fdyn,actfield,3,1);

conv_check_rans = 1;

/*--------------------- check if nonlinear iteration for rans converge */  
goto nonlniter_check;
} /* endif kapomega_yeah */

/*========================== E N D: ====================================* 
 /                    K A P P A - O M E G A                            /
*========================== E N D: =====================================*/ 
endrans:

/*-------------------------------------------------- steady state check */
if (fdyn->stchk==iststep)
{
   iststep=0;
   steady = fluid_steadycheck(fdyn,actfield,numeq_total[rans]);
}

/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
solserv_sol_copy(actfield,0,1,1,3,1);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]   
           and transform kinematic to real pressure --------------------*/
solserv_sol_copy(actfield,0,1,0,3,actpos);
fluid_transpres(actfield,0,0,actpos,fdyn->numdf-1,0);

fluid_copysol_tu(fdyn,actfield,3,actpos,1);

if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;

/*-------- calculate wall shear velocity and c_f (for TURBULENCE MODEL) */
   dynvar->washvel  = ZERO;
   container.actndis=0;
   *action = calc_fluid_shearvelo;
   container.nii= 0;
   container.nim= 0;
   container.nif= 0;
   calelm(actfield,actsolv,actpart,actintra,k_array,-1,
          &container,action);
#ifdef PARALLEL
fluid_reduceshstr(actintra,actfield,dynvar);
#endif
if (par.myrank==0) printf("wall shear velocity: %10.3E \n",dynvar->washvel);

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
   time_a.a.dv[actpos] = fdyn->time;
   actpos++;
}

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->time <= fdyn->maxtime && steady==0)
{
/*------ check time and number of steps if steady state cannot achieved */
if (kapomega_yeah==0 && (fdyn->step >= 0.5*fdyn->nstep || fdyn->time > 0.5*fdyn->maxtime)) 
{
   kapomega_yeah = 1;
}

   goto timeloop;
}

/*------------------------ if steady state was achieved goto kappa-eps */
if (kapomega_yeah==0 && steady==1) 
{
   kapomega_yeah = 1;
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
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&ftimerhs_ko_a);
amdel(&ftimerhs_kappa_a);
amdel(&ftimerhs_omega_a);
amdel(&ftimerhs_pro_a);
amdel(&ftimerhs_pro_kappa_a);
amdel(&ftimerhs_pro_omega_a);
amdel(&fiterhs_a);

if (par.myrank==0 && ioflags.fluid_vis_file==1)
amdel(&time_a);

solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&rhs_ko,1);
solserv_del_vec(&sol_ko,1);
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_isi_tu_1 */ 

#endif
