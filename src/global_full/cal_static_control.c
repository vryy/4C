#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
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
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];


/*----------------------------------------------------------------------*
 |  routine to control static execution                  m.gee 6/01     |
 *----------------------------------------------------------------------*/
void calsta()
{
#ifdef DEBUG 
dstrc_enter("calsta");
#endif
/*----------------------------------------------------------------------*/

if (statvar->geolinear==1 && statvar->geononlinear==1)
   dserror("linear and nonlinear static analysis on");
   
if (statvar->geolinear==1) 
{
   stalin();
}
if (statvar->geononlinear==1) 
{
   stanln(); 
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calsta */




/*----------------------------------------------------------------------*
 |  routine to control linear static structural analysis    m.gee 6/01  |
 *----------------------------------------------------------------------*/
void stalin() 
{
int           i;                /* a counter */
int           numeq;            /* number of equations on this proc */
int           numeq_total;      /* total number of equations over all procs */
int           numdf;            /* number of dofs over all procs */
int           init;             /* init flag for solver */
int           actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

SOLVAR       *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION    *actpart;          /* pointer to the fields PARTITION structure */
FIELD        *actfield;         /* pointer to the structural FIELD */
INTRA        *actintra;         /* pointer to the fields intra-communicator structure */
CALC_ACTION  *action;           /* pointer to the structures cal_action enum */

SPARSE_TYP    array_typ;     /* type of psarse system matrix */
#ifdef DEBUG 
dstrc_enter("stalin");
#endif
/*----------------------------------------------------------------------*/
/*------------ the distributed system matrix, which is used for solving */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)calloc(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/*    intracommunicator (in case of linear statics, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
numdf       = actfield->numdf;
/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];
switch(array_typ)
{
case mds:/*--------------------------------- system array is mds matrix */
   numeq       = actsolv->sysarray[actsysarray].mds->numeq;
   numeq_total = numeq;
break;
case msr:/*--------------------------------- system array is msr matrix */
   numeq       = actsolv->sysarray[actsysarray].msr->numeq;
   numeq_total = actsolv->sysarray[actsysarray].msr->numeq_total;
break;
case parcsr:/*--------------------------- system array is parcsr matrix */
   numeq       = actsolv->sysarray[actsysarray].parcsr->numeq;
   numeq_total = actsolv->sysarray[actsysarray].parcsr->numeq_total;
break;
case ucchb:/*----------------------------- system array is ucchb matrix */
   numeq       = actsolv->sysarray[actsysarray].ucchb->numeq;
   numeq_total = actsolv->sysarray[actsysarray].ucchb->numeq_total;
break;
case dense:/*----------------------------- system array is dense matrix */
   numeq       = actsolv->sysarray[actsysarray].dense->numeq;
   numeq_total = actsolv->sysarray[actsysarray].dense->numeq_total;
break;
case rc_ptr:/*----------------------- system array is row/column matrix */
   numeq       = actsolv->sysarray[actsysarray].rc_ptr->numeq;
   numeq_total = actsolv->sysarray[actsysarray].rc_ptr->numeq_total;
break;
default:
   dserror("unknown type of global matrix");
break;
}
/*---------------------------------- number of rhs and solution vectors */
actsolv->nrhs=2;
actsolv->nsol=2;
solserv_create_vec(&(actsolv->rhs),2,numeq_total,numeq,"DV");
solserv_create_vec(&(actsolv->sol),2,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
for (i=0; i<actsolv->nrhs; i++)
   solserv_zero_vec(&(actsolv->rhs[i]));
for (i=0; i<actsolv->nsol; i++)
   solserv_zero_vec(&(actsolv->sol[i]));
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
*action = calc_struct_init;
calinit(actfield,actpart,action);
/*------call element routines to calculate & assemble stiffness matrice */
*action = calc_struct_linstiff;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,NULL,0,0,action);
/*----------------------------------- call rhs-routines to assemble rhs */
*action = calc_struct_eleload;
calrhs(
          actfield,
          actsolv,
          actpart,
          actintra,
          actsysarray,
          &(actsolv->rhs[actsysarray]),
          &(actsolv->rhs[actsysarray+1]),
          0,
          action
         );
/*--------------------------------------------- add the two rhs vectors */
solserv_add_vec(&(actsolv->rhs[actsysarray+1]),&(actsolv->rhs[actsysarray]));
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
/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL 
free(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of stalin */
