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
#include "../ale2/ale2.h"
#include "fsi_prototypes.h"
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
\brief  solving for mesh displacements

<pre>                                                            ck 06/03

in this function the mesh of a multifield problem is solved as a
discontinous pseudo structure. All nodes within an element are connected
by springs while torsional springs are attached at the nodes.
Displacements are prescribed at the fluid structure
interface and at the free surface as Dirichlet boundary conditions

</pre>

\param  *actfield   FIELD          (i)     ale field
\param   mctrl      INT            (i)     control flag
\warning
\return void

\sa   calling: calelm(), monitoring(), ale_setdirich_increment_fsi(),
               plot_ale_quality()
      called by: fsi_ale()

*----------------------------------------------------------------------*/
void fsi_ale_spring(
                     FIELD            *actfield,
                     INT               mctrl
	          )
{
INT        	i;		  /* a counter  		                        */
static INT      numaf;            /* actual number of ale field                         */
static INT 	numeq;  	  /* number of equations on this proc                   */
static INT 	numeq_total;	  /* total number of equations over all procs           */
INT        	init;		  /* init flag for solver                               */
static INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
static INT      constsysarray;    /* second system matrix for steepest descent 		*/
static INT      numsys;           /* number of system matrices				*/
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
#ifdef DEBUG
dstrc_enter("fsi_ale_spring");
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
container.isdyn   = 1;
container.actndis = 0;
container.pos     = 0;
actpos=0;
outstep=0;
pssstep=0;
restartstep=0;
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
numsys = 1;  /* the number of system matrices */
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
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif

/*------------------------------------------------------ check proc ---*/
if (actintra->intra_fieldtyp != ale) goto end;
/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];

/*--- in the case of relaxation parameter calculations we would like to
      have a second system matrix which remains constant over time.     */
if (fsidyn->ifsi == 6)
{
   /*----------- one sysarray already exists, so copy the mask of it to */
   /* reallocate the vector of sparse matrices and the vector of there types
      formerly lenght 1, now lenght 2 */
   numsys++;
   actsolv->nsysarray=2;
   constsysarray = 1;
   actsolv->sysarray_typ =
   (SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
   actsolv->sysarray =
   (SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
   /*-copy the matrices sparsity mask */
   solserv_alloc_cp_sparsemask(  actintra,
                               &(actsolv->sysarray_typ[actsysarray]),
                               &(actsolv->sysarray[actsysarray]),
                               &(actsolv->sysarray_typ[constsysarray]),
                               &(actsolv->sysarray[constsysarray]));
}

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
for (i = 0; i<numsys; i++)
{
   solver_control(
                       actsolv,
                       actintra,
                     &(actsolv->sysarray_typ[i]),
                     &(actsolv->sysarray[i]),
                     &(actsolv->sol[actsysarray]),
                     &(actsolv->rhs[actsysarray]),
                       init
                    );
   /*------------------------------ init the dist sparse matrix to zero */
   /*            NOTE: Has to be called after solver_control(init=1) */
   solserv_zero_mat(
                       actintra,
                       &(actsolv->sysarray[i]),
                       &(actsolv->sysarray_typ[i])
                      );
   /*---------------------------- init the assembly for sparse matrices */
   init_assembly(actpart,actsolv,actintra,actfield,i,0);
}

/*------------------------------- init the element calculating routines */
*action = calc_ale_init_spring;
calinit(actfield,actpart,action,&container);

/*--------------------------------------------------- init ale field ---*/
fsi_init_ale(actfield,2);

/*--------------------------------------------------- check for restart */
if (genprob.restart!=0)
{
   restart_read_aledyn(genprob.restart,adyn,actfield,actpart,actintra);
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
{
   out_monitor(actfield,numaf,ZERO,1);
   monitoring(actfield,numaf,actpos,adyn->time);
}

/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp_gid==1 && par.myrank==0)
out_gid_domains(actfield);
#endif

if (fsidyn->ifsi == 6) /* if we do steepest descent method! */
{
   /* init linear element routines also */
   *action = calc_ale_init;
   calinit(actfield,actpart,action,&container);
   /* evaluate linear element stiffnesses (remain constant) */
   *action = calc_ale_stiff;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.kstep        = 0;
   calelm(actfield,actsolv,actpart,actintra,constsysarray,-1,&container,action);
}
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
/*------------------------------------------------------- check proc ---*/
if (actintra->intra_fieldtyp != ale) break;

if (par.myrank==0)
{
   printf("Solving ALE (springs)...\n");
   printf("\n");
}

/*--------------------------------------- sequential staggered schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi<3 && fsidyn->ifsi>0)
   solserv_sol_copy(actfield,0,3,3,1,0);

dsassert(fsidyn->ifsi!=3,"ale-solution handling not implemented for algo with DT/2-shift!\n");

/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[actsysarray]));
solserv_zero_vec(&(actsolv->sol[actsysarray]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich_increment_fsi(actfield,adyn,actpos);
/*----------------------------------------------------------------------*/
solserv_zero_mat(actintra,
		 &(actsolv->sysarray[actsysarray]),
		 &(actsolv->sysarray_typ[actsysarray])
	         );
/*--- call element routines to calculate & assemble stiffness matrix ---*/
*action = calc_ale_stiff_spring;
container.dvec         = NULL;
container.dirich       = dirich;
container.global_numeq = numeq_total;
switch (adyn->measure_quality)
{
   /*--------------------------------*/
   case no_quality:
      container.quality = 0;
   break;
   /*--------------------------------*/
   case aspect_ratio:
      container.quality = 1;
   break;
   /*--------------------------------*/
   case corner_angle:
      container.quality = 2;
   break;
   /*--------------------------------*/
   case min_detF:
      container.quality = 3;
   break;
   default:
      dserror("element quality unknown");
   break;
}
calelm(actfield,actsolv,actpart,actintra,
       actsysarray,-1,&container,action);
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
/*---------------------allreduce the result and put it to sol_increment */
solserv_result_incre(
                     actfield,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
                     0);

/*---------- add actual solution increment to sol (to serve output): ---*/
/* step 1: */
   /* copy prev. solution form sol_increment[1][j] to sol[actpos][j] ---*/
solserv_sol_copy(actfield,0,1,0,1,actpos);
/* step 2: */
       /*----------- add actual solution increment to sol[actpos][j] ---*/
solserv_sol_add(actfield,0,1,0,0,actpos,1.0);

/* save actual solution to be the previous one in the next time step ---*/
/*--------- copy actual solution from sol[actpos][i] to sol_mf[1][i] ---*/
solserv_sol_copy(actfield,0,0,3,actpos,1);

/*----------------------------------------------------------------------*/
if (fsidyn->ifsi>=4 || fsidyn->ifsi<0)
break;

/*======================================================================*
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
/*------------------------------------ for iterative staggared schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi>=4 || fsidyn->ifsi==-1)
   solserv_sol_copy(actfield,0,3,3,1,0);

/*--------------------- to get the corrected free surface position copy
  --------------------------------- from sol_mf[1][j] to sol[actpos][j] */
solserv_sol_copy(actfield,0,3,0,0,actpos);

/*---------------------------- from sol_mf[1][j] to sol_increment[1][j] */
solserv_sol_copy(actfield,0,3,1,1,1);

/*------------------------------------------- print out results to .out */
outstep++;
pssstep++;
restartstep++;

if (outstep==adyn->updevry_disp && ioflags.ale_disp_file==1)
{
    outstep=0;
    out_sol(actfield,actpart,actintra,adyn->step,actpos);
/*    if (par.myrank==0) out_gid_sol("displacement",actfield,actintra,adyn->step,0);*/
}
/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numaf,actpos,adyn->time);

if (pssstep==fsidyn->uppss && ioflags.fluid_vis_file==1 && par.myrank==0)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = adyn->time;
   actpos++;
}

/*------------------------------------------------- write restart data */
if (restartstep==fsidyn->uprestart)
{
   restartstep=0;
   restart_write_aledyn(adyn,actfield,actpart,actintra);
}

/*--------------------------------------- do mesh quality statistics ---*/
if (container.quality)
  plot_ale_quality(actfield,fsidyn->step,actintra,actpart);
/*--------------------------------------------------------------------- */

break;



/*======================================================================*
 |                A D D I T I O N A L   S O L U T I O N                 |
 |    for determination of relaxation parameter via steepest descent    |
 |                            l i n e a r !                             |
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
ale_setdirich(actfield,adyn,6);

/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actfield,actsolv,actpart,actintra,constsysarray,-1,dirich,
        numeq_total,0,&container,action);

/*------ add rhs from fsi coupling (-> prescribed displacements) to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[constsysarray]),
  &(actsolv->sysarray[constsysarray]),&(actsolv->rhs[actsysarray]),
  dirich,1.0);

/*--------------------------------------------------------- call solver */
/*---- system matrix has been calculated within initialisation phase ---*/
init=0;
solver_control(
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[constsysarray]),
                  &(actsolv->sysarray[constsysarray]),
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
                     &(actsolv->sysarray[constsysarray]),
                     &(actsolv->sysarray_typ[constsysarray]),
                     0);

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
} /* end of fsi_ale_spring */

#endif
/*! @} (documentation module close)*/
