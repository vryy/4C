/*!----------------------------------------------------------------------
\file
\brief contains the routine 'dyn_ale' which controles the calculation
of the problemtype ale

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "ale3.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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

<pre>                                                              mn 06/02 
This routine  controls the  execution of pure ale problems.
Initialization, once the calculation of the stiffness matrix, and multiple
calculation of the rhs and solving.

</pre>

\warning There is nothing special to this routine
\return void                                               
\sa   calling: ale_calelm(), ale_setdirich(), ale_rhs(); 
      called by: caldyn()

*----------------------------------------------------------------------*/
void dyn_ale()
{
#ifdef D_ALE
INT           i;                /* a counter */
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
STRUCT_DYNAMIC *sdyn;               /* pointer to structural dynamic input data */

ARRAY         dirich_a;         /* redundant vector of full length for dirichlet-part of rhs*/
DOUBLE       *dirich;

SPARSE_TYP    array_typ;        /* type of psarse system matrix */
#ifdef DEBUG 
dstrc_enter("dyn_ale");
#endif
/*----------------------------------------------------------------------*/
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
sdyn        =   alldyn[0].sdyn;
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray);
/*------------------------------- init the element calculating routines */
*action = calc_ale_init;
calinit(actfield,actpart,action);
/*------call element routines to calculate & assemble stiffness matrice */
*action = calc_ale_stiff;
ale_calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,action);
/*--------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);   
/*------------------------------------------- print out results to .out */
#ifdef PARALLEL 
if (ioflags.struct_disp_file==1)
{
  if (par.myrank==0)  out_gid_domains(actfield);
}
#endif
/*===================================================================== */
/*                      T I M E L O O P                                 */
/*===================================================================== */
timeloop:
t0 = ds_cputime();
/*--------------------------------------------- increment step and time */
sdyn->step++;
/*------------------------------------------------ set new absolue time */
sdyn->time += sdyn->dt;
/*------------------------------ init the created dist. vectors to zero */
solserv_zero_vec(&(actsolv->rhs[actsysarray]));
solserv_zero_vec(&(actsolv->sol[actsysarray]));
/*--------------------------------------------------------------------- */
amzero(&dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
ale_setdirich(actfield,sdyn);
/*------------------------------- call element-routines to assemble rhs */
*action = calc_ale_rhs;
ale_rhs(actfield,actsolv,actpart,actintra,actsysarray,-1,dirich,numeq_total,0,action);
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
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*------------------------------------------- print out results to .out */
if (ioflags.struct_disp_file==1)
{
    out_sol(actfield,actpart,actintra,sdyn->step,0);
/*    out_gid_sol_init(); */
/*    out_gid_msh();*/
    if (par.myrank==0) out_gid_sol("displacement",actfield,actintra,sdyn->step,0);
}
/*--------------------------------------------------------------------- */
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
fprintf(allfiles.out_err,"TIME for ALE step %d is %f sec\n",sdyn->step,t1-t0);
printf("TIME for ALE step %d is %f sec\n",sdyn->step,t1-t0);
/*-------------------------------------- check time and number of steps */
if (sdyn->step < sdyn->nstep-1 && sdyn->time <= sdyn->maxtime)
goto timeloop;

/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
#endif
return;
} /* end of dyn_ale */
/*! @} (documentation module close)*/
