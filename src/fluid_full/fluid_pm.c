/*!----------------------------------------------------------------------
\file
\brief projection method algorithm for fluid
 
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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"

static INT       veldis=0;           /* flag, if vel. discr. is used    */
static INT       predis=1;           /* flag, if pres. discr. is used   */
static INT       predof=0;           /* flag for pressure dof           */ 
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
extern INT     numcurve=1;            /*set the number of curves to one */
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];

static FLUID_DYNAMIC *fdyn;

/*!---------------------------------------------------------------------  
\brief projection method algorithm for fluid problems

<pre>                                                        basol 08/02

this routine conrols projection method algorithms for fluid
problems.

see:
GRESHO, Philip M.
"On the theory of Semi-Implicit Projection Methods for viscous 
incompressible flow and its implementation via a finite element
method that also introduces a nearly consistent mass matrix
part 1 & part 2 "
International Journal for Numerical Methods in Fluids, 
VOL 11, pp. 587-659, (1990)

</pre>     

\return void 
\warning this is all in progress 

------------------------------------------------------------------------*/
void fluid_pm(void)
{
INT             actcurve;           /* actual timecurve                 */
INT             i;                  /* simply a counter                 */
INT             kk;                 /* counter for discretisations      */
INT             numeq[3];           /* number of equations on this proc */
INT             numeq_total[2];     /* total number of equations        */
INT             numeq_total_full[2];/* number of dofs on this proc      */
INT             numeq_full[2];      /* total number of dofs             */
INT             numeq_oll;
INT             numeq_total_oll;
INT             numeq_save;
INT             nsysarray=2;        /* two system matrices K and M_hat  */
INT             actsysarray;        /* number of actual sysarray        */
INT             k_array=0;          /* index of K-matrix in solver      */
INT             m_array=1;          /* index of M_hat-matrix in solver  */
INT             init;               /* flag for solver_control call     */
INT             actpos=0;           /* actual position in sol. history  */
INT             outstep=0;          /* counter for output control       */
INT             resstep=0;          /* counter for output control       */
INT             pssstep=0;	    /* counter for output control	*/
ARRAY           fvelrhs1_a;          
DOUBLE         *fvelrhs1;           /* velocity rhs term1=MUn+dt*(fn-N(Un)*Un)  */
ARRAY           fvelrhs2_a;          
DOUBLE         *fvelrhs2;           /* velocity rhs term2= -dt*CPn              */
ARRAY           fdirich_a; 
DOUBLE         *fdirich;            /* dirichlet forces for velocity            */
ARRAY           lmass_a;            /* inverse of reduced lumped mass matrix    */
DOUBLE         *lmass;
ARRAY           resultvel_a;        /* intermediate velocities                  */
DOUBLE         *resultvel;
ARRAY           resultpre_a;        /* lagrange multipliers                     */
DOUBLE         *resultpre;
ARRAY           ctranspufull_a;     /* result of C transpose u_full multipl.    */
DOUBLE         *ctranspufull;
ARRAY           ctimesphi_a;        /* result of C times phi multipl.   */
DOUBLE         *ctimesphi;
ARRAY           fullvel_a;          /* full velocity vector             */
DOUBLE         *fullvel;
ARRAY           time_a;             /* stored time                      */

DIST_VECTOR    *work1;              /* distr. work1 vector for rhs      */
DIST_VECTOR    *work2;              /* distr. work2 vector for rhs      */
DIST_VECTOR    *rhs_v;              /* distr. RHS for solving velocity  */
DIST_VECTOR    *rhs_p;              /* distr. RHS for solving pressure  */
DIST_VECTOR    *sol_v;              /* distr. velocity solution         */
DIST_VECTOR    *sol_p;              /* distr. lagrange multiplier sol.  */
DIST_VECTOR    *sol_pnew;           /* distr. pressure solution         */
SOLVAR         *actsolv;            /* pointer to active sol. structure */
SOLVAR         *presolv;
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
CONTAINER       container;          /* contains variables defined in container.h */      
OLL            *gradmatrix_oll;     /* full gradient matrix             */
OLL            *rgradmatrix_oll;    /* reduced gradient matrix          */
FLUID_STRESS    str;           

#ifdef DEBUG 
dstrc_enter("fluid_pm");
#endif

/*======================================================================*
 |                   I N I T I A L I S A T I O N                        |   
/*======================================================================*/
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
fdyn        = alldyn[genprob.numff].fdyn;
actfield    = &(field[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
str         = str_none;
container.fieldtyp = fluid;
numeq_full[veldis] = actpart->pdis[veldis].numnp*2;
numeq_total_full[veldis] = actfield->dis[veldis].numdf;
numeq_full[predis] = actpart->pdis[predis].numnp;
numeq_total_full[predis] = actfield->dis[predis].numdf;
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

/* it seems to be better if the second discretisation has it's own direct
   solver (colsol or spooles) and the matrix is stored in oll format.
   so we first reallocate the solver: */
if (genprob.numfld!=1)
dserror("Projection method only for single field problem\n");

solv = (SOLVAR*)CCAREALLOC(solv,2*sizeof(SOLVAR));
/*-------------------------------------- set pointer to pressure solver */
actsolv = &(solv[0]);
presolv = &(solv[1]);
presolv->fieldtyp = fluid;


/*-------------------------------------- now we set some default values */
#ifdef PARALLEL
#ifndef SPOOLES_PACKAGE
dserror("SPOOLES package is not compiled in");
#endif
presolv->solvertyp = SPOOLES_nonsym;
#else
presolv->solvertyp = colsol_solver;
presolv->colsolvars = (COLSOLVARS*)CCACALLOC(1,sizeof(COLSOLVARS));
#endif
presolv->parttyp = cut_elements;
presolv->matrixtyp = oll_matrix;

/* -------------------------- create solver for pressure discretisation */
presolv->nsysarray = 1;
presolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(presolv->nsysarray,sizeof(SPARSE_TYP));
presolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(presolv->nsysarray,sizeof(SPARSE_ARRAY));

presolv->sysarray_typ[0] = oll;
presolv->sysarray[0].oll = (OLL*)CCACALLOC(1,sizeof(OLL));


/*-----------------------------------------------------------------------*
  at the moment there exists two sparse matrices with the same mask:
  K-matrix: k_array = 0;
  M_hat-matrix: m_array=1;
 *------------------------------------------------------------------------*/
actsolv->nsysarray=nsysarray;
actsolv->sysarray_typ = 
(SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
actsolv->sysarray = 
(SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));

/*--------------------- copy the mask from the k_array to the mass_array */
solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[k_array]),
                            &(actsolv->sysarray[k_array]),
                            &(actsolv->sysarray_typ[m_array]),
                            &(actsolv->sysarray[m_array]));

/*------------------------------- loop the matrices and intitialise them */
for(kk=0;kk<nsysarray;kk++)
{
actsysarray=kk;
   /*----------------------------- init the dist sparse matrices to zero */
   solserv_zero_mat(
                    actintra,
                    &(actsolv->sysarray[actsysarray]),
                    &(actsolv->sysarray_typ[actsysarray])
                   );
   /*------------------------- get global and local number of equations */
   solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                      actsolv->sysarray_typ[actsysarray],
                      &numeq[kk],
                      &numeq_total[kk]);

} /* end of loop over sys_arrays */

/*-------------------------------------- initialise solver for pressure */
numeq_total_oll = actfield->dis[predis].numeq;
oll_numeq(actfield, actpart, actintra, predis, &numeq_oll);

oll_open(presolv->sysarray[0].oll, numeq_oll, numeq_total_oll, 
	 actfield, actpart, actintra, predis);

numeq[predis] = numeq_oll;
numeq_total[predis] = numeq_total_oll;

/*--------------------------------------------------- open gradmatrices */
gradmatrix_oll  = (OLL*)CCACALLOC(1,sizeof(OLL));
rgradmatrix_oll = (OLL*)CCACALLOC(1,sizeof(OLL));
numeq_save = actfield->dis[veldis].numeq;
actfield->dis[veldis].numeq = numeq_total_full[veldis];
oll_open(gradmatrix_oll,numeq_full[veldis],numeq_total_full[veldis],
         actfield,actpart,actintra,veldis);
actfield->dis[veldis].numeq = numeq_save;	 
oll_open(rgradmatrix_oll,numeq[veldis],numeq_total[veldis],
         actfield,actpart,actintra,veldis);

/*--------------------------------------------- output to the screen */
#ifdef PARALLEL
printf("number of velocity eqations on PROC %3d : %10d \n", 
        par.myrank,numeq[veldis]);
if (par.myrank==0)
   printf("total number of velocity equations: %10d \n",numeq_total[veldis]);	
printf("number of pressure eqations on PROC %3d : %10d \n", 
        par.myrank,numeq[predis]);
if (par.myrank==0)
   printf("total number of pressure equations:       %10d \n",numeq_total[predis]);
#else
printf("total number of velocity equations: %10d \n",numeq_total[veldis]);
printf("total number of pressure equations: %10d \n",numeq_total[predis]);
#endif

/*-------------------------- allocate 1 dist. vector 'rhs' for velocity */
solserv_create_vec(&rhs_v,1,numeq_total[veldis],numeq[veldis],"DV");
solserv_zero_vec(rhs_v);
/*----------------------------- allocate 1 dist. vector 'work1' for rhs */
solserv_create_vec(&work1,1,numeq_total[veldis],numeq[veldis],"DV");
solserv_zero_vec(work1);
/*----------------------------- allocate 1 dist. vector 'work2' for rhs */
solserv_create_vec(&work2,1,numeq_total[veldis],numeq[veldis],"DV");
solserv_zero_vec(work2);
/*-------------------------- allocate 1 dist. vector 'rhs' for pressure */
solserv_create_vec(&rhs_p,1,numeq_total[predis],numeq[predis],"DV");
solserv_zero_vec(rhs_p); 
/*------------------------ allocate 1 dist. vector 'sol_v' for velocity */
solserv_create_vec(&sol_v,1,numeq_total[veldis],numeq[veldis],"DV");
solserv_zero_vec(sol_v);
/*------------------------ allocate 1 dist. vector 'sol_p' for pressure */
solserv_create_vec(&sol_p,1,numeq_total[predis],numeq[predis],"DV");
solserv_zero_vec(sol_p);
/*---------------------------------- allocate 1 dist. vector 'sol_pnew' */
solserv_create_vec(&sol_pnew,1,numeq_total[predis],numeq[predis],"DV");
solserv_zero_vec(sol_pnew);
/*----------------------------------- allocate 1 redundant vector lmass */
lmass = amdef("lmass",&lmass_a,numeq_total[veldis],1,"DV");
amzero(&lmass_a);
/*------------------ allocate redundant vectors fvelrhs1 of full lenght 
/*          these are used by the element routines to assemble the  RHS */
fvelrhs1    = amdef("fvelrhs1",&fvelrhs1_a,numeq_total[veldis],1,"DV");
amzero(&fvelrhs1_a);
/*------------------ allocate redundant vectors fvelrhs2 of full lenght */
fvelrhs2    = amdef("fvelrhs2",&fvelrhs2_a,numeq_total[veldis],1,"DV");
amzero(&fvelrhs2_a);
/*------------------- allocate redundant vectors fdirich of full lenght */
fdirich     = amdef("fdirich",&fdirich_a,numeq_total[veldis],1,"DV");
amzero(&fdirich_a);
/*------------------------------- allocate 1 redundant vector resultvel */
resultvel = amdef("resultvel",&resultvel_a,numeq_total[veldis],1,"DV");
amzero(&resultvel_a);
/*------------------------------- allocate 1 redundant vector ctimesphi */
ctimesphi = amdef("ctimesphi",&ctimesphi_a,numeq_total[veldis],1,"DV");
amzero(&ctimesphi_a);
/*----------------------------- allocate redundant vector full velocity */
fullvel = amdef("fullvel",&fullvel_a,actfield->dis[0].numdf,1,"DV");
amzero(&fullvel_a);
/*------------------------------- allocate 1 redundant vector resultpre */
resultpre = amdef("resultpre",&resultpre_a,numeq_total[predis],1,"DV");
amzero(&resultpre_a);
/*-------------------------------- allocate 1 redundant vector ctranspu */
ctranspufull  = amdef("ctranspufull",&ctranspufull_a,actfield->dis[1].numdf,1,"DV");
amzero(&ctranspufull_a);

/*---------------------------- allocate one vector for storing the time */
if (ioflags.fluid_vis_file==1 )
amdef("time",&time_a,1000,1,"DV");

/*---------------------------------------------- initialise fluid field */
fluid_init(actpart,actintra,actfield,action,&container,4,str);		     
actpos=0;

/*---------------------------------------- init all applied time curves */
for (actcurve=0; actcurve<numcurve; actcurve++)
   dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime); 

/*--------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield);

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
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               sol_v,
               rhs_v,
               init);
solver_control(presolv, actintra,
               &(presolv->sysarray_typ[0]),
               &(presolv->sysarray[0]),
               sol_p,
               rhs_p,
               init);

/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,k_array,veldis);
init_assembly(actpart,presolv,actintra,actfield,0,predis);

/*---------------------------------- allocate fluid integration data ---*/
alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*------------------------------- print out initial data to .flavia.res */
if (ioflags.fluid_sol_gid==1 && par.myrank==0) 
{
   out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
   out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
}

/*----------------------------------------------------------------------*
/*                    END OF THE INITIALISATION       
/*----------------------------------------------------------------------*/

fdyn->theta = ONE;
fdyn->pro_profile = 1;  /* parabolic profile              */
/*fdyn->pro_profile = 2 ; /* constant profile               */
/*fdyn->pro_profile = 3;    /* parabolic profile for cylinder */

/*======================================================================*
 |                         A  -  M A T R I X                            |   
/*======================================================================*/
/*-------------------- calculate the A-matrix : A = CT*ML-1*C-----------*/
/*-------------------- A  corresponds to the laplace operator ----------*/
if (par.myrank==0)
printf("\nCalculate the A-matrix: A = CT*ML-1*C\n");
*action = calc_fluid_amatrix;
fluid_pm_calamatrix( actfield,
		      actintra,
		      actpart,
		      actsolv,
		      &(presolv->sysarray[0]),
		      fdyn,
		      gradmatrix_oll,
		      rgradmatrix_oll,
		      lmass,
		      numeq_total,
		      numeq,
		      numeq_full,
		      numeq_total_full,
		      action,
		     &container		      
		    );

/*----------------- initialize the solver again-------------------------*/
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               sol_v,
               rhs_v,
               init);
solver_control(presolv, actintra,
               &(presolv->sysarray_typ[0]),
               &(presolv->sysarray[0]),
               sol_p,
               rhs_p,
               init);	 

/*----------------------------------------------------------------------*/
fdyn->acttime=ZERO;
fdyn->step=0;
fdyn->dta = fdyn->dt;

/*======================================================================*
 |                         T I M E   L O O P                            |   
/*======================================================================*/
if (par.myrank==0) printf("\nTIMELOOP:\n\n");
timeloop:
fdyn->step++; 
fdyn->acttime += fdyn->dta;
/*----------------------------------------------- output to the screen */
if (par.myrank==0)
printf("TIME: %11.4E/%11.4E  DT = %11.4E  Projection Method  STEP = %4d/%4d \n",
        fdyn->acttime,fdyn->maxtime,fdyn->dta,fdyn->step,fdyn->nstep);

/*--------------------------- intitialise global matrix and global rhs */
solserv_zero_mat(actintra,&(actsolv->sysarray[k_array]),
                 &(actsolv->sysarray_typ[k_array]));
solserv_zero_mat(actintra,&(actsolv->sysarray[m_array]),
                 &(actsolv->sysarray_typ[m_array]));

/*------------------------------------------ initialise working arrays */
solserv_zero_vec(rhs_v);
solserv_zero_vec(sol_v);
solserv_zero_vec(work1);
solserv_zero_vec(work2);
solserv_zero_vec(rhs_p);
solserv_zero_vec(sol_p);

/*--------------------- set dirichlet boundary conditions for  timestep */
switch (fdyn->pro_profile)
{
case 1:
   fluid_setdirich_parabolic(actfield);
break;
case 2:
   fluid_setdirich(actfield,3);
break;
case 3:
   fluid_setdirich_cyl(actfield);
break;
default:
   dserror("unknown velocity profile!\n");
}
/*---------------------------------- form the LHS matrices & RHS vector */
amzero(&fvelrhs1_a);
amzero(&fvelrhs2_a);
amzero(&fdirich_a);

/*----------------------------------------------------------------------*	            
| Stiffness matrix and Mass matrix are also calculated at every         |
| time step because K includes the BTD stabilization term which         |
| is dependent on the velocity and changing at every time step          |
/*----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*
 | instead of two times element calculation for LHS and RHS,             |
 | both actions are combined into one action called                      |
 | calc_fluid_f2pro_rhs_both dirichlet forces are also calculated        |
 | within the time loop.                                                 |
/*-----------------------------------------------------------------------*/
*action = calc_fluid_f2pro_rhs_both;
container.fiterhs      = fvelrhs1;
container.nii          = 1;
container.ftimerhs     = fvelrhs2;
container.nif          = 1;
container.nim          = 0;
container.fidrichrhs   = fdirich;
container.global_numeq = numeq_total[veldis]; 
container.kstep        = 0;
container.actndis      = veldis;
fdyn->pro_calmat     = 1;
fdyn->pro_caldirich  = 1;
fdyn->pro_kvv        = 1;
fdyn->pro_calrhs     = 1;
fdyn->pro_mvv        = 1;
fdyn->pro_calveln    = 1;
fdyn->pro_lum        = 0;
calelm(actfield,actsolv,actpart,actintra,k_array,m_array,&container,action);

/*--------------------------------------------- scale K-matrix K = K*dt */
solserv_scal_mat(&(actsolv->sysarray_typ[k_array]),
                   &(actsolv->sysarray[k_array]),
		   fdyn->dta*fdyn->theta);
/*-------------------------------------------------------- K = M + dt*K */
solserv_add_mat(actintra,
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               &(actsolv->sysarray_typ[m_array]),
               &(actsolv->sysarray[m_array]),
               ONE);      
/*---------------------------------------------------- M_hat = M * ML-1 */
fluid_pm_matlmatmul(&(actsolv->sysarray[m_array]),
                   &(actsolv->sysarray_typ[m_array]),
		   lmass);

/*----------------------------------------------------------------------*
 |           rhs_v = fvelrhs1 + (MMl-1)*fvelrhs2  + fdirich             | 
 *----------------------------------------------------------------------*/
/*------------ assemble fvelrhs1  --------------------------------------*/
assemble_vec(actintra,
             &(actsolv->sysarray_typ[veldis]),
             &(actsolv->sysarray[veldis]),
             rhs_v,
             fdirich,
             1.0);
/*--------------------------- assemble fvelrhs1 = MUn + dt*(fn-N(Un)Un) */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[veldis]),
             &(actsolv->sysarray[veldis]),
             rhs_v,
             fvelrhs1,
             1.0);
/*----------------------------------------- assemble fvelrhs2 = -dt*CPn */	     
assemble_vec(actintra,
             &(actsolv->sysarray_typ[veldis]),
             &(actsolv->sysarray[veldis]),
             work2,
             fvelrhs2,
             1.0);
/*----------------- multiply (MMl-1) by work2(fvelrhs2) in dist. format */
solserv_sparsematvec(actintra,
                     work1,
                     &(actsolv->sysarray[m_array]),
                     &(actsolv->sysarray_typ[m_array]),
                     work2);
/*-------------------------------- rhs term for velocity rhs_v += work1 */
solserv_add_vec(work1,rhs_v,1.0);		     	     	     

/*--------------------------------- solving for intermediate velocities */
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[k_array]),
               &(actsolv->sysarray[k_array]),
               sol_v,
               rhs_v,
               init);

/*----------------------------------------------------------------------*
 | REMARK: We are dealing with the full velocity vector because we use  | 
 | the full gradient matrix in our calculations, also the A matrix is   | 
 | calculated by using the full gradient matrix. One can use the        |
 | reduced gradient matrix but we believe that full A matrix is a       |
 | better approximation to the laplace operator over the whole domain   |                              |
 *----------------------------------------------------------------------*/
solserv_result_incre(actfield,actintra,sol_v,3,
                     &(actsolv->sysarray[k_array]),
                     &(actsolv->sysarray_typ[k_array]),veldis);
		     
/*--------------------------------------- form the full velocity vector */
fluid_pm_fullvel(actfield,veldis,fullvel,3); 

/*------------------------------------- redistribute sol_v to resultvel */
solserv_reddistvec(sol_v,
                   &(actsolv->sysarray[k_array]),
                   &(actsolv->sysarray_typ[k_array]),
                   resultvel,numeq_total[veldis],actintra);


/*======================================================================*
 |              S O L V I N G   F O R    P R E S S U R E                |
 |                       A*(phi)=CT*(U~n+1)-gn+1                        |
 *======================================================================*
 | REMARK: RHS of the PPE equation can be formed in two different ways  |
 | one way is to consider a vector called g				|
 | the other way is to multiply the full divergence operator by the full|
 | velocity vector, both ways end up with the same result and the second|
 | is implemented below. 						|
 *----------------------------------------------------------------------*/

/*----------------------- CT*U~n+1 is formed using full velocity vector */
/*----------------and full grad matrix, then the pressure bc is reduced */
fluid_pm_matvecmul(&ctranspufull_a,&fullvel_a,gradmatrix_oll,
                   actintra,numeq_full[veldis],1);

/*-------------- assemble ctranspu-->CT*(U~n+1) into distributed vector */
assemble_vec(actintra,
             &(presolv->sysarray_typ[0]),
             &(presolv->sysarray[0]),
             rhs_p,
             ctranspufull,
             1.0);

/*----------------------------------- solve for the lagrange parameters */
init=0;
solver_control(presolv, actintra,
               &(presolv->sysarray_typ[0]),
               &(presolv->sysarray[0]),
               sol_p,
               rhs_p,
               init);

/*------------ do the projection and get the divergence free velocities 
               Un+1 = U~n+1 - (ML-1)*C*(phi)                            */
 
/*--------------------------------------- redistribute the sol_p vector */
solserv_reddistvec(sol_p,
                  &(presolv->sysarray[0]),
                  &(presolv->sysarray_typ[0]),
                   resultpre,numeq_total[predis],actintra);

/*------------------------------------------------ C*phi multiplication */
fluid_pm_matvecmul(&ctimesphi_a,&resultpre_a,rgradmatrix_oll,
                    actintra,numeq[veldis],0);		   

/*--------------------------------------- Ml-1*ctimesphi multiplication */
fluid_pm_lmatmulvec(ctimesphi,lmass,numeq_total[veldis]);

/*----------------------------------------------------------------------* 
 | subtract the correction term from intermediate velocities and get    | 
 | the divergence free velocities                                       |
 *----------------------------------------------------------------------*/
amadd(&resultvel_a,&ctimesphi_a,-ONE,0);

/*----- assemble the divergence free velocities into distributed vector */
solserv_zero_vec(sol_v);
assemble_vec(actintra,
             &(actsolv->sysarray_typ[k_array]),
             &(actsolv->sysarray[k_array]),
             sol_v,
             resultvel,
             1.0);

/*---------------------- return divergence free velocities to the nodes */
solserv_result_incre(actfield,actintra,sol_v,3,
                     &(actsolv->sysarray[k_array]),
                     &(actsolv->sysarray_typ[k_array]),veldis);


/*----------------------- update the pressure values-->Pn+1=Pn+2*phi/dt */
solserv_add_vec(sol_p,sol_pnew,TWO/fdyn->dta);

/*--------------------------------- return pressure values to the nodes */
solserv_result_incre(actfield,actintra,sol_pnew,3,
                     &(presolv->sysarray[0]),
                     &(presolv->sysarray_typ[0]),predis);

/*-------- copy solution from sol_increment[3][j] to sol_increment[1[j] */
solserv_sol_copy(actfield,veldis,1,1,3,1);
solserv_sol_copy(actfield,predis,1,1,3,1);


/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]   
           and transform kinematic to real pressure --------------------*/
solserv_sol_copy(actfield,veldis,1,0,3,actpos);
solserv_sol_copy(actfield,predis,1,0,3,actpos);
fluid_transpres(actfield,predis,0,actpos,predof,0);

/*------------------------- copy pressure from pre discr. to vel discr. */
fluid_pm_pretovel(actfield,actpos);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;
resstep++;

/*---------------------------------------------- write solution to .out */
if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
}

/*--------------------------------------- write solution to .flavia.res */
if (resstep==fdyn->upres &&ioflags.fluid_sol_gid==1 && par.myrank==0) 
{
   resstep=0;
   out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
   out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
}

/*---------------------------------------------- write solution to .pss */
if (pssstep==fdyn->uppss && ioflags.fluid_vis_file==1)
{
   pssstep=0;   
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = fdyn->acttime;   
   actpos++;
}

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime)
   goto timeloop; 
/*----------------------------------------------------------------------*
 | -->  end of timeloop                                                 |
 *----------------------------------------------------------------------*/
end:

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
      time_a.a.dv[actpos] = fdyn->acttime;   
   }   
   visual_writepss(actfield,actpos+1,&time_a);
}

/*--------------------------------------------------- cleaning up phase */
solserv_del_vec(&rhs_v,1);
solserv_del_vec(&rhs_p,1);
solserv_del_vec(&sol_v,1);
solserv_del_vec(&sol_p,1);
solserv_del_vec(&work1,1);
solserv_del_vec(&work2,1);
if (par.myrank==0 && ioflags.fluid_vis_file==1 )
amdel(&time_a);
amdel(&fvelrhs1_a);
amdel(&fvelrhs2_a);
amdel(&fdirich_a);
amdel(&lmass_a);
amdel(&resultvel_a);
amdel(&resultpre_a);
amdel(&ctranspufull_a);
amdel(&ctimesphi_a);
amdel(&fullvel_a);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_pm */ 
#endif
/*! @} (documentation module close)*/
