/*!----------------------------------------------------------------------
\file
\brief ale part of fsi-problems

*----------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fsi_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
 extern struct _IO_FLAGS     ioflags;
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
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
 extern struct _FILES  allfiles;
 
 
/*!----------------------------------------------------------------------
\brief  solving for linear mesh displacements 

<pre>                                                          genk 10/02

in this function the mesh of a multifield problem is solved as a 
pseude structure. Displacements are prescribed at the fluid structure
interface and at the free surface as Dirichlet boundary conditions

Routine also includes calculation part for fluid mesh dynamics used to 
determine Relaxation parameter via steepest descent relaxation.

</pre>

\param  *fsidyn     FSI_DYNAMIC    (i)  			  
\param  *adyn       STRUCT_DYNAMIC (i)  			  
\param  *actfield   FIELD          (i)     ale field 			  
\param   mctrl      INT            (i)     control flag		  
\param   numfa      INT            (i)     number of ale field	  
\warning 
\return void  
                                             
\sa   calling: calelm(), monitoring(), ale_setdirich(), ale_rhs()
      called by: fsi_ale()

*----------------------------------------------------------------------*/
void fsi_ale_lin(
                  FSI_DYNAMIC      *fsidyn,
                  ALE_DYNAMIC      *adyn,
                  FIELD            *actfield,
                  INT               mctrl,
                  INT               numfa
	       )
{
INT        	i;		  /* a counter  		                        */
static INT 	numeq;  	  /* number of equations on this proc                   */
static INT 	numeq_total;	  /* total number of equations over all procs           */
INT        	init;		  /* init flag for solver                               */
static INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
static INT      actpos;           /* actual position in nodal solution history          */
static INT      outstep;          /* counter for output to .out                         */
static INT      pssstep;          /* counter for output to .pss                         */
static SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
static PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
static INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
static CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
static ARRAY         dirich_a;    /* red. vector of full length for dirich-part of rhs  */
static DOUBLE       *dirich;
static ARRAY         time_a;      /* stored time                                        */
static CONTAINER     container;   /* contains variables defined in container.h          */
static SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

#ifdef DEBUG 
dstrc_enter("fsi_ale_lin");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1: 
adyn->dt=fsidyn->dt;
adyn->maxtime=fsidyn->maxtime;
adyn->nstep=fsidyn->nstep;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
container.isdyn   = 0;  
container.actndis = 0;    
actpos=1;
outstep=0;
pssstep=0;
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actsolv     = &(solv[numfa]);
actpart     = &(partition[numfa]);
action      = &(calc_action[numfa]);
container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL 
actintra    = &(par.intra[numfa]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = ale;
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the structural */
/*    intracommunicator (in case of linear statics, this should be all) */
if (actintra->intra_fieldtyp != ale) goto end;

/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(actsolv->sysarray[actsysarray],
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
printf("PROC  %3d | FIELD ALE       | number of equations      : %10d \n", 
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD ALE       | total number of equations: %10d \n",numeq_total);	
if (par.myrank==0) printf("\n\n");

/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=1;
actsolv->nsol=1;
solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");

/*------------------------------ init the created dist. vectors to zero */
for (i=0; i<actsolv->nrhs; i++)
   solserv_zero_vec(&(actsolv->rhs[i]));
for (i=0; i<actsolv->nsol; i++)
   solserv_zero_vec(&(actsolv->sol[i]));

/*----------- create a vector of full length for dirichlet part of rhs */
dirich = amdef("intforce",&dirich_a,numeq_total,1,"DV");
/*--------------------------- allocate one vector for storing the time */
if (par.myrank==0) amdef("time",&time_a,1000,1,"DV");

/*--------------------------------------------------- initialize solver */
init=1;
solver_control(  
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[actsysarray]),
                  &(actsolv->rhs[actsysarray]),
                    init
                 );
/*--------------------------------- init the dist sparse matrix to zero */
/*               NOTE: Has to be called after solver_control(init=1) */
solserv_zero_mat(
                    actintra,
                    &(actsolv->sysarray[actsysarray]),
                    &(actsolv->sysarray_typ[actsysarray])
                   );
/*----------------------------- init the assembly for ONE sparse matrix */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);

/*------------------------------- init the element calculating routines */
*action = calc_ale_init;
calinit(actfield,actpart,action,&container);

/*------call element routines to calculate & assemble stiffness matrice */
*action = calc_ale_stiff;
container.dvec         = NULL;
container.dirich       = NULL;
container.global_numeq = 0;
container.kstep        = 0;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action); 

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numfa,0,0,adyn->time);

/*------------------------------------------- print out results to .out */
#ifdef PARALLEL 
/*if (ioflags.ale_disp_gid==1 && par.myrank==0)
out_gid_domains(actfield);*/
#endif

break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        * 
 *======================================================================*/
case 2:
/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != ale) break;

if (par.myrank==0)
{
   printf("Solving ALE (classic linear)...\n");
   printf("\n");
}

/*--------------------------------------- sequential staggered schemes: 
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi<3) solserv_sol_copy(actfield,0,3,3,1,0);

dsassert(fsidyn->ifsi!=3,"ale-solution handling not implemented for algo with DT/2-shift!\n");
   
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[actsysarray]));
solserv_zero_vec(&(actsolv->sol[actsysarray]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);

/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich(actfield,adyn,actpos,0);

/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actfield,actsolv,actpart,actintra,actsysarray,-1,dirich,
        numeq_total,0,&container,action);

/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[actsysarray]),
     dirich,1.0);

/*--------------------------------------------------------- call solver */
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[actsysarray]),
                  &(actsolv->rhs[actsysarray]),
                    init
                 );
/*-------------------------allreduce the result and put it to the nodes */
solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     actpos,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*---------------------allreduce the result and put it to sol_increment */
solserv_result_incre(
                     actfield,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
   
/*----------------- copy from nodal sol_increment[0][j] to sol_mf[1][j] */
solserv_sol_copy(actfield,0,1,3,0,1);

if (fsidyn->ifsi>=4)
break;

/*======================================================================* 
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
/*------------------------------------ for iterative staggared schemes:
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi>=4) solserv_sol_copy(actfield,0,3,3,1,0);

/*------------------------------------------- print out results to .out */
outstep++;
pssstep++;
if (outstep==adyn->updevry_disp && ioflags.ale_disp_file==1)
{ 
    outstep=0;
    out_sol(actfield,actpart,actintra,adyn->step,actpos);
/*    if (par.myrank==0) out_gid_sol("displacement",actfield,actintra,adyn->step,0);*/
}
/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numfa,actpos,adyn->step,adyn->time);

if (pssstep==fsidyn->uppss && ioflags.fluid_vis_file==1 && par.myrank==0)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = adyn->time;   
   actpos++;
} 
/*--------------------------------------------------------------------- */   

break;


/*======================================================================*
 |                A D D I T I O N A L   S O L U T I O N                 |
 |    for determination of relaxation parameter via steepest descent    |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        * 
 * sol_mf[2][i]        ... grid position in relaxation parameter calc.  *
 * sol_increment[0][i] ... displacement used to determine omega_RELAX   *
 *======================================================================*/
case 6:
/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != ale) break;

if (par.myrank==0)
{
   printf("          - Solving ALE ... \n");
}

dsassert(fsidyn->ifsi!=3,"ale-solution handling not implemented for algo with DT/2-shift!\n");
   
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[actsysarray]));
solserv_zero_vec(&(actsolv->sol[actsysarray]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);

/*-------------------------set dirichlet boundary conditions on at time */
/* note: ale Dirichlet boundary conditions different from ZERO are not 
         helpful for fsi coupling problems and would cause trouble here.
	 But there's no test to set all ordinary dbc = 0.0 !!!
	 The required Dirichlet boundary conditions from fsi coupling 
	 are calculated here.						*/
ale_setdirich(actfield,adyn,actpos,6);

/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actfield,actsolv,actpart,actintra,actsysarray,-1,dirich,
        numeq_total,0,&container,action);

/*------ add rhs from fsi coupling (-> prescribed displacements) to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[actsysarray]),
     dirich,1.0);

/*--------------------------------------------------------- call solver */
/*--------- the system matrix has been calculated within the init phase */
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[actsysarray]),
                  &(actsolv->rhs[actsysarray]),
                    init
                 );

/*-------------- allreduce the result and put it to sol_increment[0][i] */
solserv_result_incre(
                     actfield,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );

break;


/*======================================================================* 
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:
/*- there are only procs allowed in here, that belong to the ale -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != ale) break;
if (pssstep==0) actpos--;

/*---------------------------------------------- print out to .mon file */
if (ioflags.monitor==1 && par.myrank==0)
out_monitor(actfield,numfa);


/*------------------------------------------- print out results to .out */
if (outstep!=0 && ioflags.ale_disp_file==1)
out_sol(actfield,actpart,actintra,adyn->step,actpos);

/*------------------------------------------- print out result to 0.pss */
if (ioflags.fluid_vis_file==1 && par.myrank==0)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = adyn->time;   
   }   
   visual_writepss(actfield,actpos+1,&time_a);
}

/*------------------------------------------------------------- tidy up */ 
if (par.myrank==0) amdel(&time_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
break;
default:
   dserror("Parameter out of range: mctrl \n");
} /* end switch (mctrl) */

/*----------------------------------------------------------------------*/
end:

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fsi_ale_lin */

#endif
/*! @} (documentation module close)*/
