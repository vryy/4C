#include "../headers/standardtypes.h"
#include "../headers/solution.h"
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
int        i;                /* a counter */
int        numeq;            /* number of equations on this proc */
int        numeq_total;      /* total number of equations over all procs */
int        numdf;            /* number of dofs over all procs */
int        init;             /* init flag for solver */
int        actsysarray;      /* active sparse system matrix in actsolv->sysarray[] */

SOLVAR    *actsolv;          /* pointer to the fields SOLVAR structure */
PARTITION *actpart;          /* pointer to the fields PARTITION structure */
FIELD     *actfield;         /* pointer to the structural FIELD */
INTRA     *actintra;         /* pointer to the fields intra-communicator structure */

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
calinit(actfield,actpart);
/*------call element routines to calculate & assemble stiffness matrice */
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,NULL,0,0,1);
/*----------------------------------- call rhs-routines to assemble rhs */
calrhs(
          actfield,
          actsolv,
          actpart,
          actintra,
          actsysarray,
          &(actsolv->rhs[actsysarray]),
          &(actsolv->rhs[actsysarray+1]),
          0,
          6
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
