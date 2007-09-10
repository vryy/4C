/*!----------------------------------------------------------------------
\file
\brief contains the routine 'dyn_ale' which controles the calculation
of the problemtype ale

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#include "../headers/standardtypes.h"
#include "ale3.h"
#include "../ale2/ale2.h"
#include "../solver/solver.h"
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

/*!
\addtogroup Ale
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief controls the  execution of pure ale problems

<pre>                                                            ck 05/03
This routine  controls the  execution of pure ale problems. Depending on
the type of the ale problem different dynamic routines are called.

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: dyn_ale_lin(), dyn_ale_nln(), dyn_ale_2step(),
               dyn_ale_spring(), dyn_ale_laplace()
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale()
{
#ifdef D_ALE

#ifdef DEBUG
dstrc_enter("dyn_ale");
#endif
switch (alldyn->adyn->typ)
{
/*---------------------------------------- purely linear calculation ---*/
   case classic_lin:
      dyn_ale_lin();
   break;
/*------------------- incremental calculation stiffened with min J_e ---*/
   case min_Je_stiff:
      dyn_ale_nln();
   break;
/*--------------------------------------------- two step calculation ---*/
/*  calculation in two steps per timestep following Chiandussi et al. in
    'A simple method for automatic update of finite element meshes'
    Commun. Numer. Meth. Engng. 2000; 16: 1-19                          */
   case two_step:
      dyn_ale_2step();
   break;
/*--------------------------------------------------- spring analogy ---*/
/*  calculation following Farhat et al. in 'Torsional springs for
    two-dimensional dynamic unstructured fluid meshes' Comput. Methods
    Appl. Mech. Engrg. 163 (1998) 231-245 */
    case springs:
       dyn_ale_spring();
    break;
/*------------------------------------------------ Laplace smoothing ---*/
/*  calculation following Loehner et al. in 'Improved ALE mesh velocities
    for moving bodies' Commun. num. methd. engng. 12 (1996) 599-608 */
    case laplace:
       dyn_ale_laplace();
    break;
/*---------------------------------------------------------- default ---*/

    case LAS:
       dserror("Aletype LAS not implementd for pure ALE!!");
    break;


   default:
    dserror("unknown ale typ");
   break;
}
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale*/


/*!----------------------------------------------------------------------
\brief controls the  execution of purely linear ale problems

<pre>                                                              mn 06/02
This routine  controls the execution of linear pure ale problems.
Initialization, once the calculation of the stiffness matrix, and multiple
calculation of the rhs and solving.

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: ale_calelm(), ale_setdirich(), ale_rhs();
      called by: dyn_ale()

*----------------------------------------------------------------------*/
#ifdef D_ALE
void dyn_ale_lin()
{
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

INT             actcurve;           /* indice of active time curve */
DOUBLE          t0,t1;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
ALE_DYNAMIC  *adyn;             /* pointer to ale dynamic input data */
CONTAINER     container;        /* contains variables defined in container.h */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

INT           disnum_calc;
INT           disnum_io;
ARRAY_POSITION *ipos;

#ifdef DEBUG
dstrc_enter("dyn_ale_lin");
#endif
/*----------------------------------------------------------------------*/
container.isdyn = 0;            /* static calculation */
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
actsolv             = &(solv[0]);
actpart             = &(partition[0]);
action              = &(calc_action[0]);
adyn                =   alldyn[0].adyn;


#ifdef SUBDIV
if (actfield->subdivide > 0)
{
  disnum_calc = 1;
  disnum_io   = ioflags.output_dis;
}
else
#endif
  disnum_calc = disnum_io = 0;


ipos   = &(actfield->dis[disnum_calc].ipos);

  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;


/*------------ the distributed system matrix, which is used for solving */
actsysarray=disnum_calc;

container.fieldtyp  = actfield->fieldtyp;
container.disnum    = disnum_calc;

#ifdef PARALLEL
actintra    = &(par.intra[0]);
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
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
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
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(     actfield,disnum_calc,
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
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,adyn->nstep,adyn->dt,adyn->maxtime);
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp==1)
{
  if (par.myrank==0)  out_gid_domains(actfield,disnum_io);
}
#endif

if (par.myrank==0) printf("Performing pure ALE problem linear\n\n");
/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
adyn->step++;
/*------------------------------------------------ set new absolue time */
adyn->time += adyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich(actfield,disnum_calc,adyn,0);
/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actsolv,actpart,disnum_calc,actintra,actsysarray,-1,dirich,numeq_total,&container,action);
/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(     actfield,disnum_calc,
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

/* allreduce the result and put it to the nodes */
solserv_result_incre(
                     actfield,
                     disnum_calc,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );

/* solution has to be in XYZ co-system */
locsys_trans_sol(actfield,disnum_calc,1,0,1);

/* copy from nodal sol_increment[0][j] to sol[0][j] */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol,0,0);


#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (actfield->subdivide > 0)
{
  solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
}
#endif




/* print out results */
if (ioflags.ale_disp==1)
{
  if (ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
  if (par.myrank==0 && ioflags.output_gid==1)
  {
    out_gid_sol("displacement",actfield,disnum_io,actintra,adyn->step,0,adyn->time);
  }
}


/* measure time for this step */
t1 = ds_cputime();
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
   printf("TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
}

/*-------------------------------------- check time and number of steps */
if (adyn->step < adyn->nstep-1 && adyn->time <= adyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dyn_ale_lin */
#endif





/*!----------------------------------------------------------------------
\brief controls the  execution of nonlinear ale problems

<pre>                                                            ck 05/03
This routine controls the execution of nonlinear ale problems with
deformation dependent stiffening on elemental basis via minimal Jacobian
determinant.

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: ale_calelm(), ale_setdirich()
               plot_ale_quality();
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale_nln()
{
#ifdef D_ALE
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

INT           actcurve;         /* indice of active time curve */
DOUBLE        t0,t1;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
ALE_DYNAMIC  *adyn;             /* pointer to structural dynamic input data */
CONTAINER     container;        /* contains variables defined in container.h */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

INT           disnum_calc;
INT           disnum_io;
ARRAY_POSITION *ipos;



#ifdef DEBUG
dstrc_enter("dyn_ale_nln");
#endif
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
actsolv             = &(solv[0]);
actpart             = &(partition[0]);
action              = &(calc_action[0]);
adyn                =   alldyn[0].adyn;


#ifdef SUBDIV
if (actfield->subdivide > 0)
{
  disnum_calc = 1;
  disnum_io   = ioflags.output_dis;
}
else
#endif
  disnum_calc = disnum_io = 0;


container.fieldtyp  = actfield->fieldtyp;
container.disnum    = disnum_calc;
container.isdyn     = 1;
container.pos       = 0;

ipos   = &(actfield->dis[disnum_calc].ipos);

  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;

/*------------ the distributed system matrix, which is used for solving */
actsysarray=disnum_calc;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintra    = &(par.intra[0]);
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
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=1;
actsolv->nsol=1;
solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------- create a vector of full length for dirichlet part of rhs ---*/
dirich = amdef("intforce",&dirich_a,numeq_total,1,"DV");
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(     actfield,disnum_calc,
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
*action = calc_ale_init_nln;
calinit(actfield,actpart,action,&container);
/*--------------------------------------- init solution to zero -------*/
/* solution on sol_increment[1][j] */
solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,1);
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,adyn->nstep,adyn->dt,adyn->maxtime);
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp==1)
{
  if (par.myrank==0)  out_gid_domains(actfield,disnum_io);
}
#endif
if (par.myrank==0) printf("Performing pure ALE problem nonlinear and stiffened\n\n");

/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
adyn->step++;
/*------------------------------------------------ set new absolue time */
adyn->time += adyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich_increment(actfield,disnum_calc,adyn);
/*----------------------------------------------------------------------*/
solserv_zero_mat(actintra,
		 &(actsolv->sysarray[actsysarray]),
		 &(actsolv->sysarray_typ[actsysarray])
	         );
/*----call element routines to calculate & assemble stiffness matrix */
if (adyn->step <= adyn->num_initstep)
   *action = calc_ale_stiff_stress;
else
   *action = calc_ale_stiff_nln;
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
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(   actfield,disnum_calc,actsolv,
                  actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[0]),
                 &(actsolv->rhs[0]),
                   init
                 );
/*-------------------------allreduce the result and put it to the nodes */
/*---------------------  write increment to increment vector place 0 ---*/
solserv_result_incre(
		     actfield,
                     disnum_calc,
		     actintra,
		     &(actsolv->sol[0]),
		     0,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]));
/*--------- update nodal solution values by adding actual increments ---*/
   /* sol_increment.a.da[1][j] += sol_increment.a.da[0][j]; */
solserv_sol_add(actfield,disnum_calc,node_array_sol_increment,node_array_sol_increment,0,1,1.0);
   /* sol.a.da[0][j] = sol_increment.a.da[1][j]; */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol,1,0);



#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (actfield->subdivide > 0)
{
  solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
}
#endif



/* print out results */
if (ioflags.ale_disp==1)
{
  if (ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
  if (par.myrank==0 && ioflags.output_gid==1)
  {
    out_gid_sol("displacement",actfield,disnum_io,actintra,adyn->step,0,adyn->time);
  }
}


/*--------------------------------------- do mesh quality statistics ---*/
plot_ale_quality(actfield,disnum_calc,adyn->step,(adyn->time-adyn->dt),actintra,actpart);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
   printf("TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
}
/*-------------------------------------- check time and number of steps */
if (adyn->step < adyn->nstep && adyn->time < adyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale_nln */



/*!----------------------------------------------------------------------
\brief controls the  execution of ale problems in 2 steps per time step

<pre>                                                              ck 06/03
This routine  controls the  execution of ale problems following Chiandussi
et al. in 'A simple method for automatic umpdate of finite element meshes'
Commun. Numer. Engng. 2000; 16: 1-19
The calculation performes two solution steps per timestep. A first one to
get a strain distribution based on spatially constant stiffness within the
increment. In the second step the previous results are used to obtain a
spatially variying stiffness.

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: ale_calelm(), ale_setdirich_increment(),
               plot_ale_quality();
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale_2step()
{
#ifdef D_ALE
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

INT           actcurve;         /* indice of active time curve */
DOUBLE        t0,t1;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
ALE_DYNAMIC  *adyn;             /* pointer to structural dynamic input data */
CONTAINER     container;        /* contains variables defined in container.h */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of sparse system matrix */

INT           disnum_calc;
INT           disnum_io;
ARRAY_POSITION *ipos;



#ifdef DEBUG
dstrc_enter("dyn_ale_2step");
#endif

/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
actsolv             = &(solv[0]);
actpart             = &(partition[0]);
action              = &(calc_action[0]);
adyn                =   alldyn[0].adyn;

/* to follow Chiandussi et al. calculation is performed in two steps,
   first step linear, second step with modified stiffness
   reference calculation is performed incremental                       */

#ifdef SUBDIV
if (actfield->subdivide > 0)
{
  disnum_calc = 1;
  disnum_io   = ioflags.output_dis;
}
else
#endif
  disnum_calc = disnum_io = 0;


container.fieldtyp  = actfield->fieldtyp;
container.disnum    = disnum_calc;
container.isdyn     = 1;
container.pos       = 0;

ipos   = &(actfield->dis[disnum_calc].ipos);

  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;

/*------------ the distributed system matrix, which is used for solving */
actsysarray=disnum_calc;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintra    = &(par.intra[0]);
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
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=1;
actsolv->nsol=1;
solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------- create a vector of full length for dirichlet part of rhs ---*/
dirich = amdef("intforce",&dirich_a,numeq_total,1,"DV");
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(     actfield,disnum_calc,
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
*action = calc_ale_init_step2;
calinit(actfield,actpart,action,&container);
/*--------------------------------------- init solution to zero -------*/
/* solution on sol_increment[1][j] */
solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,1);
/*---------------- put actual min and max stiffening factor in place ---*/
container.min = 0.10;
container.max = 500.0;
container.min_stiff = container.min;
container.max_stiff = container.max;
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,adyn->nstep,adyn->dt,adyn->maxtime);
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp==1)
{
  if (par.myrank==0)  out_gid_domains(actfield,disnum_io);
}
#endif
if (par.myrank==0) printf("Performing pure ALE problem in two steps\n\n");

/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
adyn->step++;
/*------------------------------------------------ set new absolue time */
adyn->time += adyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
/* write dirichlet cond. from dt on sol_increment */
ale_setdirich_increment(actfield,disnum_calc,adyn);
/*----------------------------------------------------------------------*/
solserv_zero_mat(
		    actintra,
		    &(actsolv->sysarray[actsysarray]),
		    &(actsolv->sysarray_typ[actsysarray])
		   );
/*----call element routines to calculate & assemble stiffness matrix */
*action = calc_ale_stiff;   /* first (linear) reference step */
container.dvec         = NULL;
container.dirich       = dirich;
container.global_numeq = numeq_total;
calelm(actfield,actsolv,actpart,actintra,
       actsysarray,-1,&container,action);
/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(   actfield,disnum_calc,actsolv,
                  actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[0]),
                 &(actsolv->rhs[0]),
                   init
                 );
/*-------------------- allreduce the result and put it to sol_increment */
solserv_result_incre(
		     actfield,
                     disnum_calc,
		     actintra,
		     &(actsolv->sol[0]),
		     0,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]));
/*-------------------------------------- reset system matrix to zero ---*/
solserv_zero_mat(
     		 actintra,
     		 &(actsolv->sysarray[actsysarray]),
     		 &(actsolv->sysarray_typ[actsysarray])
     		);
/*------------------------------------------ rezero dirichlet vector ---*/
amzero(&dirich_a);
/*---- put actual min and max stiffening factor in place for scaling ---*/
container.min = container.min_stiff; /* min stiffness of prev time step */
container.max = container.max_stiff; /* max stiffness of prev time step */
container.max_stiff = 0.0;
container.min_stiff = 1.E10;
/*-call element routines to calculate & assemble stiffness matrix */
*action = calc_ale_stiff_step2;
container.kstep        = adyn->step;
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
/*--------------------------- init the created dist. vectors to zero ---*/
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------- add rhs from prescribed displacements to rhs ---*/
assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
   &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
   dirich,1.0);
/*------------------------------------------------------ call solver ---*/
init=0;
solver_control(actfield,disnum_calc,actsolv,
     	       actintra,
     	      &(actsolv->sysarray_typ[actsysarray]),
     	      &(actsolv->sysarray[actsysarray]),
     	      &(actsolv->sol[0]),
     	      &(actsolv->rhs[0]),
     		init
     	      );
/*----------- allreduce the result and write it to the sol_increment ---*/
solserv_result_incre(actfield,
                     disnum_calc,
     		     actintra,
     		    &(actsolv->sol[0]),
     		     0,
     		    &(actsolv->sysarray[actsysarray]),
     		    &(actsolv->sysarray_typ[actsysarray]));
/*--------- update nodal solution values by adding actual increments ---*/
   /* sol_increment.a.da[1][j] += sol_increment.a.da[0][j]; */
solserv_sol_add(actfield,disnum_calc,node_array_sol_increment,node_array_sol_increment,0,1,1.0);
   /* sol.a.da[0][j] = sol_increment.a.da[1][j]; */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol,1,0);



#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (actfield->subdivide > 0)
{
  solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
}
#endif



/* print out results */
if (ioflags.ale_disp==1)
{
  if (ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
  if (par.myrank==0 && ioflags.output_gid==1)
  {
    out_gid_sol("displacement",actfield,disnum_io,actintra,adyn->step,0,adyn->time);
  }
}


/*--------------------------------------- do mesh quality statistics ---*/
plot_ale_quality(actfield,disnum_calc,adyn->step,adyn->time,actintra,actpart);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
   printf("TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
}
/*-------------------------------------- check time and number of steps */
if (adyn->step < adyn->nstep && adyn->time < adyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale_2step */



/*!----------------------------------------------------------------------
\brief controls the  execution of pure ale problems

<pre>                                                            ck 05/03
This routine controls the execution of nonlinear pure ale problems that
are solved by means of spring analogy with lineal and torsional springs.
The calculation follows Farhat et al. in 'Torsional springs for two-
dimensional dynamic unstructured fluid meshes' in Comput. Methods Appl. Mech.
Engrg. 163 (1998) 231-245

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: ale_calelm(), ale_setdirich(),
               plot_ale_quality();
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale_spring()
{
#ifdef D_ALE
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

INT           actcurve;         /* indice of active time curve */
DOUBLE        t0,t1;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
ALE_DYNAMIC  *adyn;             /* pointer to structural dynamic input data */
CONTAINER     container;        /* contains variables defined in container.h */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

INT           disnum_calc;
INT           disnum_io;
ARRAY_POSITION *ipos;



#ifdef DEBUG
dstrc_enter("dyn_ale_spring");
#endif
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
actsolv             = &(solv[0]);
actpart             = &(partition[0]);
action              = &(calc_action[0]);
adyn                =   alldyn[0].adyn;

#ifdef SUBDIV
if (actfield->subdivide > 0)
{
  disnum_calc = 1;
  disnum_io   = ioflags.output_dis;
}
else
#endif
  disnum_calc = disnum_io = 0;


container.fieldtyp  = actfield->fieldtyp;
container.disnum    = disnum_calc;
container.isdyn     = 1;
container.pos       = 0;


ipos   = &(actfield->dis[disnum_calc].ipos);

  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;


/*------------ the distributed system matrix, which is used for solving */
actsysarray=disnum_calc;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintra    = &(par.intra[0]);
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
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=1;
actsolv->nsol=1;
solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------- create a vector of full length for dirichlet part of rhs ---*/
dirich = amdef("intforce",&dirich_a,numeq_total,1,"DV");
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(     actfield,disnum_calc,
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
 *action = calc_ale_init_spring;
calinit(actfield,actpart,action,&container);
/*--------------------------------------- init solution to zero -------*/
/* solution on sol_increment[1][j] */
solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,1);
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,adyn->nstep,adyn->dt,adyn->maxtime);
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp==1)
{
  if (par.myrank==0)  out_gid_domains(actfield,disnum_io);
}
#endif

if (par.myrank==0) printf("Performing pure ALE problem with springs\n\n");

/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
adyn->step++;
/*------------------------------------------------ set new absolue time */
adyn->time += adyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich_increment(actfield,disnum_calc,adyn);
/*----------------------------------------------------------------------*/
solserv_zero_mat(actintra,
		 &(actsolv->sysarray[actsysarray]),
		 &(actsolv->sysarray_typ[actsysarray])
	         );
/*----call element routines to calculate & assemble stiffness matrix */
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
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(   actfield,disnum_calc,actsolv,
                  actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[0]),
                 &(actsolv->rhs[0]),
                   init
                 );
/*-------------------------allreduce the result and put it to the nodes */
/*---------------------  write increment to increment vector place 0 ---*/
solserv_result_incre(
		     actfield,
                     disnum_calc,
		     actintra,
		     &(actsolv->sol[0]),
		     0,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]));
/*--------- update nodal solution values by adding actual increments ---*/
   /* sol_increment.a.da[1][j] += sol_increment.a.da[0][j]; */
solserv_sol_add(actfield,disnum_calc,node_array_sol_increment,node_array_sol_increment,0,1,1.0);
   /* sol.a.da[0][j] = sol_increment.a.da[1][j]; */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol,1,0);



#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (actfield->subdivide > 0)
{
  solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
}
#endif



/* print out results */
if (ioflags.ale_disp==1)
{
  if (ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
  if (par.myrank==0 && ioflags.output_gid==1)
  {
    out_gid_sol("displacement",actfield,disnum_io,actintra,adyn->step,0,adyn->time);
  }
}


/*--------------------------------------- do mesh quality statistics ---*/
plot_ale_quality(actfield,disnum_calc,adyn->step,adyn->time,actintra,actpart);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
   printf("TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
}
/*-------------------------------------- check time and number of steps */
if (adyn->step < adyn->nstep && adyn->time < adyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale_spring */




/*!----------------------------------------------------------------------
\brief controls the  execution of pure ale problems

<pre>                                                            ck 06/03
This routine controls the execution of nonlinear pure ale problems
with Laplace smoothing.
The calculation follows Loehner et al. in 'Improved ale mesh velocities
for moving bodies' commun. numer. Methods engineering Vol. 12, 599-608
(1996)
spatially varying diffusivity, however, is not implemented as described
there.

</pre>

\warning There is nothing special to this routine
\return void
\sa   calling: ale_calelm(), ale_setdirich(),
               plot_ale_quality();
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale_laplace()
{
#ifdef D_ALE
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

INT           actcurve;         /* indice of active time curve */
DOUBLE        t0,t1;

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */
ALE_DYNAMIC  *adyn;             /* pointer to structural dynamic input data */
CONTAINER     container;        /* contains variables defined in container.h */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

INT           disnum_calc;
INT           disnum_io;
ARRAY_POSITION *ipos;



#ifdef DEBUG
dstrc_enter("dyn_ale_laplace");
#endif
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
actsolv             = &(solv[0]);
actpart             = &(partition[0]);
action              = &(calc_action[0]);
adyn                =   alldyn[0].adyn;

#ifdef SUBDIV
if (actfield->subdivide > 0)
{
  disnum_calc = 1;
  disnum_io   = ioflags.output_dis;
}
else
#endif
  disnum_calc = disnum_io = 0;



container.fieldtyp  = actfield->fieldtyp;
container.disnum    = disnum_calc;
container.isdyn     = 1;
container.pos       = 0;

ipos   = &(actfield->dis[disnum_calc].ipos);

  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;

/*------------ the distributed system matrix, which is used for solving */
actsysarray=disnum_calc;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintra    = &(par.intra[0]);
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
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=1;
actsolv->nsol=1;
solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------- create a vector of full length for dirichlet part of rhs ---*/
dirich = amdef("intforce",&dirich_a,numeq_total,1,"DV");
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(     actfield,disnum_calc,
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
*action = calc_ale_init_laplace;
calinit(actfield,actpart,action,&container);
/*--------------------------------------- init solution to zero -------*/
/* solution on sol_increment[1][j] */
solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,1);
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,adyn->nstep,adyn->dt,adyn->maxtime);
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
if (ioflags.ale_disp==1)
{
  if (par.myrank==0)  out_gid_domains(actfield,disnum_io);
}
#endif
if (par.myrank==0) printf("Performing pure ALE problem with Laplace smoothing\n\n");

/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
adyn->step++;
/*------------------------------------------------ set new absolue time */
adyn->time += adyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich_increment(actfield,disnum_calc,adyn);
/*----------------------------------------------------------------------*/
solserv_zero_mat(actintra,
		 &(actsolv->sysarray[actsysarray]),
		 &(actsolv->sysarray_typ[actsysarray])
	         );
/*----call element routines to calculate & assemble stiffness matrix */
*action = calc_ale_stiff_laplace;
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
     &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
     dirich,1.0);
/*--------------------------------------------------------- call solver */
init=0;
solver_control(   actfield,disnum_calc,actsolv,
                  actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[0]),
                 &(actsolv->rhs[0]),
                   init
                 );
/*-------------------------allreduce the result and put it to the nodes */
/*---------------------  write increment to increment vector place 0 ---*/
solserv_result_incre(
		     actfield,
                     disnum_calc,
		     actintra,
		     &(actsolv->sol[0]),
		     0,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]));
/*--------- update nodal solution values by adding actual increments ---*/
   /* sol_increment.a.da[1][j] += sol_increment.a.da[0][j]; */
solserv_sol_add(actfield,disnum_calc,node_array_sol_increment,node_array_sol_increment,0,1,1.0);
   /* sol.a.da[0][j] = sol_increment.a.da[1][j]; */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,node_array_sol,1,0);



#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (actfield->subdivide > 0)
{
  solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
}
#endif



/* print out results */
if (ioflags.ale_disp==1)
{
  if (ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
  if (par.myrank==0 && ioflags.output_gid==1)
  {
    out_gid_sol("displacement",actfield,disnum_io,actintra,adyn->step,0,adyn->time);
  }
}


/*--------------------------------------- do mesh quality statistics ---*/
plot_ale_quality(actfield,disnum_calc,adyn->step,adyn->time,actintra,actpart);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
if (par.myrank==0)
{
   fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
   printf("TIME for ALE step %d is %f sec\n",adyn->step,t1-t0);
}
/*-------------------------------------- check time and number of steps */
if (adyn->step < adyn->nstep && adyn->time < adyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale_laplace */

/*! @} (documentation module close)*/
#endif
