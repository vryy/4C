/*!----------------------------------------------------------------------
\file
\brief contains the routines 'opt_calsta',to control static execution
       and 'opt_stalin' to control linear static structural analysis
       for optimization

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef D_OPTIM                   /* include optimization code to ccarat */
/*----------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../headers/optimization.h"
#include "opt_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
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
/*!----------------------------------------------------------------------
\brief the optimization main structure
<pre>                                                            al 06/01   
defined in opt_cal_main.c
</pre>
*----------------------------------------------------------------------*/
 struct _OPTI *opt;


/*----------------------------------------------------------------------*
 |  routine to control static execution                  m.gee 6/01     |
 *----------------------------------------------------------------------*/
void opt_calsta(CALSTA_EXEC stalact)
{
#ifdef DEBUG 
dstrc_enter("opt_calsta");
#endif
/*----------------------------------------------------------------------*/

if (statvar->linear==1 && statvar->nonlinear==1)
   dserror("linear and nonlinear static analysis on");
   
if (statvar->linear==1) 
{
   opt_stalin(stalact);
}
if (statvar->nonlinear==1) 
{
   opt_stanln(stalact);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of opt_calsta */




/*----------------------------------------------------------------------*
 |  routine to control linear static structural analysis    m.gee 6/01  |
 *----------------------------------------------------------------------*/
void opt_stalin(CALSTA_EXEC stalact) 
{
INT           i;                /* a counter */
INT           numeq;            /* number of equations on this proc */
INT           numeq_total;      /* total number of equations over all procs */
INT           init;             /* init flag for solver */
INT           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */

CONTAINER     container;        /* contains variables defined in container.h */

SPARSE_TYP    array_typ;        /* type of psarse system matrix */

container.isdyn   = 0;            /* static calculation */
container.actndis = 0;            /* only one discretisation */

#ifdef DEBUG 
dstrc_enter("opt_stalin");
#endif
/*----------------------------------------------------------------------*/
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield            = &(field[0]);
container.fieldtyp  = actfield->fieldtyp;
/*----------------------------------------------------------------------*/
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/*    intracommunicator (in case of linear statics, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
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
if(stalact==calsta_init || stalact==calsta_init_solve)
{
  actsolv->nrhs=2;
  actsolv->nsol=2;
  solserv_create_vec(&(actsolv->rhs),2,numeq_total,numeq,"DV");
  solserv_create_vec(&(actsolv->sol),2,numeq_total,numeq,"DV");
}
/*------------------------------ init the created dist. vectors to zero */
for (i=0; i<actsolv->nrhs; i++)
   solserv_zero_vec(&(actsolv->rhs[i]));
for (i=0; i<actsolv->nsol; i++)
   solserv_zero_vec(&(actsolv->sol[i]));
/*--------------------------------------------------- initialize solver */
if(stalact==calsta_init || stalact==calsta_init_solve)
{
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
}
/*--------------------------------- init the dist sparse matrix to zero */
/*               NOTE: Has to be called after solver_control(init=1) */
solserv_zero_mat(
                    actintra,
                    &(actsolv->sysarray[actsysarray]),
                    &(actsolv->sysarray_typ[actsysarray])
                   );
/*----------------------------- init the assembly for ONE sparse matrix */
if(stalact==calsta_init || stalact==calsta_init_solve)
{
  init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
}
/*------------------------------- init the element calculating routines */
if(stalact==calsta_init || stalact==calsta_init_solve)
{
  *action = calc_struct_init;
  calinit(actfield,actpart,action,&container);
}
/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0&&opt->optstep==0)
if (ioflags.struct_disp_gid||ioflags.struct_stress_gid) 
   out_gid_msh();
/*------call element routines to calculate & assemble stiffness matrice */
if(stalact==calsta_init_solve || stalact==calsta_solve)
{
  *action = calc_struct_linstiff;
  container.dvec         = NULL;
  container.dirich       = NULL;
  container.actndis      = 0;
  container.global_numeq = 0;
  container.kstep        = 0;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
}
/*----------------------------------- call rhs-routines to assemble rhs */
/*-------------------------- the approbiate action is set inside calrhs */
if(stalact==calsta_init_solve || stalact==calsta_solve)
{
  *action = calc_struct_eleload; 
  container.kstep   = 0;
  container.inherit = 1;
  calrhs(actfield,actsolv,actpart,actintra,actsysarray,
       &(actsolv->rhs[actsysarray]),action,&container);
}
/*--------------------------------------------------------- call solver */
if(stalact==calsta_init_solve || stalact==calsta_solve)
{
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
}
/*-------------------------allreduce the result and put it to the nodes */
if(stalact==calsta_init_solve || stalact==calsta_solve)
{
  solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[actsysarray]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
}
/*--------------------------------------------- printout results to gid */
if (ioflags.struct_disp_gid==1 && par.myrank==0 && stalact!=calsta_init)
{
   out_gid_sol("displacement",actfield,actintra,opt->optstep,0);
   /* out_gid_domains(actfield); */
}
/*------------------------------------------ perform stress calculation */
if (ioflags.struct_stress_file==1 || ioflags.struct_stress_gid==1 && stalact!=calsta_init)
{
   *action = calc_struct_stress;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.kstep        = 0;
   calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);   
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   container.kstep = 0;
   calreduce(actfield,actpart,actintra,action,&container);
   out_sol(actfield,actpart,actintra,0,0);
   if (par.myrank==0 && stalact!=calsta_init) out_gid_sol("stress",actfield,actintra,opt->optstep,0);
}
/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of opt_stalin */
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/
