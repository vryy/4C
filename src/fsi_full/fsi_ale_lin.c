/*!----------------------------------------------------------------------
\file
\brief ale part of fsi-problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../ale3/ale3.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
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

\param  *actfield   FIELD          (i)     ale field
\param   mctrl      INT            (i)     control flag
\warning
\return void

\sa   calling: calelm(), monitoring(), ale_setdirich(), ale_rhs()
      called by: fsi_ale()

*----------------------------------------------------------------------*/
 void fsi_ale_lin(
     FIELD             *actfield,
     INT                disnum_calc,
     INT                disnum_io,
     INT                mctrl
     )

{

INT        	i;		  /* a counter  		                        */
static INT      numaf;            /* actual number of ale field                         */
static INT 	numeq;  	  /* number of equations on this proc                   */
static INT 	numeq_total;	  /* total number of equations over all procs           */
INT        	init;		  /* init flag for solver                               */
static INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
static INT      actpos;           /* actual position in nodal solution history          */
static INT      outstep;          /* counter for output to .out                         */
static INT      pssstep;          /* counter for output to .pss                         */
static INT      restartstep;      /* counter for output of restart data                 */
static SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
static PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
static INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
static CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
static ARRAY         dirich_a;    /* red. vector of full length for dirich-part of rhs  */
static DOUBLE       *dirich;
static ARRAY         time_a;      /* stored time                                        */
static CONTAINER     container;   /* contains variables defined in container.h          */
static SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

static FSI_DYNAMIC  *fsidyn;
static ALE_DYNAMIC  *adyn;

#ifdef BINIO
static BIN_OUT_FIELD out_context;
static BIN_OUT_FIELD restart_context;
#endif


#ifdef DEBUG
dstrc_enter("fsi_ale_lin");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1:
numaf  = genprob.numaf;
adyn   = alldyn[numaf].adyn;
fsidyn = alldyn[numaf+1].fsidyn;
adyn->dt=fsidyn->dt;
adyn->maxtime=fsidyn->maxtime;
adyn->nstep=fsidyn->nstep;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
container.isdyn   = 0;
container.disnum  = disnum_calc;
actpos=0;
outstep=0;
pssstep=0;
restartstep=0;
/*------------ the distributed system matrix, which is used for solving */
actsysarray = disnum_calc;
/*--------------------------------------------------- set some pointers */
actsolv     = &(solv[numaf]);
actpart     = &(partition[numaf]);
action      = &(calc_action[numaf]);
container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
actintra    = &(par.intra[numaf]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = ale;
actintra->intra_rank     = 0;
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
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));

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
                  &(actsolv->sol[0]),
                  &(actsolv->rhs[0]),
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,disnum_calc);

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

#ifdef D_MORTAR
if(adyn->coupmethod == 0) /* mortar method */
{
  /*------- redefine the size of sol_mf from 3 to 4, the fourth field is */
  /*------------- necessary to store the nodal displacements due to fsi */
  solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, 4);
}
#endif

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   actfield, actpart, actintra, disnum_io);

if(disnum_io != disnum_calc)
  init_bin_out_field(&restart_context,
      &(actsolv->sysarray_typ[actsysarray]),
      &(actsolv->sysarray[actsysarray]),
      actfield, actpart, actintra, disnum_calc);
#endif

/*--------------------------------------------------- check for restart */
if (genprob.restart!=0)
{
#if defined(BINIO)
   restart_read_bin_aledyn(adyn,
                           &(actsolv->sysarray_typ[actsysarray]),
                           &(actsolv->sysarray[actsysarray]),
                           actfield,
                           actpart,
                           disnum_calc,
                           actintra,
                           genprob.restart);
#else
   restart_read_aledyn(genprob.restart,adyn,actfield,actpart,actintra);
#endif
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
{
   out_monitor(actfield,numaf,ZERO,1);
   monitoring(actfield,disnum_calc,numaf,actpos,adyn->time);
}

/*------------------------------------------- print out results to .out */
if (ioflags.ale_disp==1 && ioflags.output_out==1)
{
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);
}


#ifdef PARALLEL
/*if (ioflags.ale_disp==1 && par.myrank==0)
out_gid_domains(actfield,disnum_io);*/
#endif

break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        *
 * sol_mf[2][i]        ... used in fsi_gradient.c                       *
 * sol_mf[3][i]        ... displacements at (n+1) stored in fsi_calc_disp4ale*
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

/*--------------------------------------- sequential staggered schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi<3 && fsidyn->ifsi>0)
solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,1,0);


dsassert(fsidyn->ifsi!=3,"ale-solution handling not implemented for algo with DT/2-shift!\n");

/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);

/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich(actfield,disnum_calc,adyn,0);

/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actsolv,actpart,disnum_calc,actintra,actsysarray,-1,dirich,
        numeq_total,&container,action);

/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);

/*--------------------------------------------------------- call solver */
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[0]),
                  &(actsolv->rhs[0]),
                    init
                 );

/* for nodes with locsys and DBCs the values would become mixed up, since
   the solution is in the xyz* co-sys, but the sol-array is in the XYZ
   co-sys, so transform DBC nodes to xyz*                               */
locsys_trans_sol_dirich(actfield,disnum_calc,1,0,0);

/*---------------------allreduce the result and put it to sol_increment */
solserv_result_incre(
                     actfield,
                     disnum_calc,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                     );

/*--------------------------------- solution has to be in XYZ co-system */
locsys_trans_sol(actfield,disnum_calc,1,0,1);

/*----------------- copy from nodal sol_increment[0][j] to sol_mf[1][j] */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol_mf,0,1);

if (fsidyn->ifsi>=4 || fsidyn->ifsi<0)
break;

/*======================================================================*
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
/*------------------------------------ for iterative staggared schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi>=4 || fsidyn->ifsi==-1)
{
   solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,0,2);
   solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,1,0);
}

/*--------------------- to get the corrected free surface position copy
  --------------------------------- from sol_mf[1][j] to sol[actpos][j] */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol,1,actpos);


#ifdef SUBDIV
      /* transfer the solution to the nodes of the master-dis */
      if (actfield->subdivide > 0)
      {
        solserv_sol_trans(actfield, disnum_calc, node_array_sol, actpos);
      }
#endif



/*---------------------------------------------- increment output flags */
outstep++;
pssstep++;
restartstep++;

if (pssstep==fsidyn->uppss && ioflags.fluid_vis==1 && par.myrank==0)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = adyn->time;
   actpos++;
}

if (outstep==adyn->updevry_disp && ioflags.ale_disp==1 && ioflags.output_out==1)
{
    outstep=0;
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,disnum_calc,numaf,actpos,adyn->time);


/*------------------------------------------------- write restart data */
if (restartstep==fsidyn->uprestart)
{
   restartstep=0;
#ifdef BINIO
   if(disnum_io != disnum_calc)
     restart_write_bin_aledyn(&restart_context, adyn);
   else
     restart_write_bin_aledyn(&out_context, adyn);
#else
   restart_write_aledyn(adyn,actfield,actpart,actintra);
#endif
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
 * sol_mf[3][i]        ... displacements at (n+1) stored in fsi_calc_disp4ale*
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
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);

/*-------------------------set dirichlet boundary conditions on at time */
/* note: ale Dirichlet boundary conditions different from ZERO are not
         helpful for fsi coupling problems and would cause trouble here.
	 But there's no test to set all ordinary dbc = 0.0 !!!
	 The required Dirichlet boundary conditions from fsi coupling
	 are calculated here.						*/
ale_setdirich(actfield,disnum_calc,adyn,6);

/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actsolv,actpart,disnum_calc,actintra,actsysarray,-1,dirich,
        numeq_total,&container,action);

/*------ add rhs from fsi coupling (-> prescribed displacements) to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);

/*--------------------------------------------------------- call solver */
/*--------- the system matrix has been calculated within the init phase */
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[actsysarray]),
                  &(actsolv->sysarray[actsysarray]),
                  &(actsolv->sol[0]),
                  &(actsolv->rhs[0]),
                    init
                 );

/*-------------- allreduce the result and put it to sol_increment[0][i] */
solserv_result_incre(
                     actfield,
                     disnum_calc,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                     );

break;

/*======================================================================*
                            Binary Output
 *======================================================================*/
case 98:
#ifdef BINIO
  if (ioflags.output_bin)
    if (ioflags.ale_disp==1)
      out_results(&out_context, adyn->time, adyn->step, actpos, OUTPUT_DISPLACEMENT);
#endif
  break;

/*======================================================================*
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:
/*- there are only procs allowed in here, that belong to the ale -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != ale) break;
if (pssstep==0) actpos--;

/*------------------------------------------- print out results to .out */
if (outstep!=0 && ioflags.ale_disp==1 && ioflags.output_out==1)
out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);

/*------------------------------------------- print out result to 0.pss */
if (ioflags.fluid_vis==1 && par.myrank==0)
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

#ifdef BINIO
destroy_bin_out_field(&out_context);
if(disnum_io != disnum_calc)
  destroy_bin_out_field(&restart_context);
#endif

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
