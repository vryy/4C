/*!----------------------------------------------------------------------
\file
\brief contains the routines to control nonlinear static execution
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
/*! 
\addtogroup OPTIMIZATION 
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;
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
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
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

/*************** prototypes ***********************/
/*----------------------------------------------------------------------*
 |  PERFORM LINEAR PREDICTOR STEP                               al 11/01|
 *----------------------------------------------------------------------*/
void nlninc(
            FIELD         *actfield,     /* the actual physical field */
            SOLVAR        *actsolv,      /* the field-corresponding solver */
            PARTITION     *actpart,      /* the partition of the proc */
            INTRA         *actintra,     /* the intra-communicator of this field */
            CALC_ACTION   *action,       /* calculation flag */
            INT            kstep,        /* the load or time step we are in */
            INT            actsysarray,  /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,          /* dist. vector of incremental residual forces used for iteration in conequ */
            STANLN        *nln_data,     /* data of the Newton-Raphson method */
            INT           *iconv,        /* type of control algorithm */
            CONTAINER     *container 	 /* contains variables defined in container.h */
          );
/*----------------------------------------------------------------------*
 |  PERFORM EQUILLIBRIUM ITERATION                           m.gee 11/01|
 |  within Newton Raphson                                               |
 *----------------------------------------------------------------------*/
void nlnequ(
            FIELD         *actfield,      /* the actual physical field */
            SOLVAR        *actsolv,       /* the field-corresponding solver */
            PARTITION     *actpart,       /* the partition of the proc */
            INTRA         *actintra,      /* the intra-communicator of this field */
            CALC_ACTION   *action,        /* calculation flag */
            INT            kstep,         /* the load or time step we are in */
            INT           *itnum,         /* number of corrector steps taken by this routine */
            INT            actsysarray,   /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,           /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,         /* dist. vector of incremental displacements */
            DIST_VECTOR   *re,            /* re[0..2] 3 vectors for residual displacements */
            INT            cdof,          /* number of dof to be controlled */
            STANLN        *nln_data,      /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,    /* type of control algorithm */
            INT            oflag,         /* type of control algorithm */
            CONTAINER     *container      /* contains variables defined in container.h */
          );
















/*----------------------------------------------------------------------*
 |  routine to control nonlinear static execution       m.gee 11/01     |
 *----------------------------------------------------------------------*/
void opt_stanln(CALSTA_EXEC stalact) 
{
INT           i;                  /* simply a counter */
INT           numeq;              /* number of equations on this proc */
INT           numeq_total;        /* total number of equations */
INT           init;               /* flag for solver_control call */
INT           actsysarray;        /* number of the active system sparse matrix */
INT           cdof;               /* the Id of the controlled dof */

INT           kstep;              /* active step in nonlinear analysis */
INT           nstep;              /* total number of steps */
INT           itnum;              /* number of iterations taken by corrector */
INT           maxiter;            /* max. number of iterations allowed */
DOUBLE        stepsi;             /* load step size for temp. storage*/

INT           mod_displ;          /* write or not dicplacements of this step in output */
INT           mod_stress;         /* write or not stresses of this step in output */

INT           mod_restart;        /* write or not valuesof this step needed for restart in output */
INT           restart;            /* if this is a restart calc -> restart =!0 elseif ==0*/
INT           restartevery;       /* how often to write restart (temp. storage) */
INT           iconv, oflag;

static DOUBLE act_loadfac;        /* current load factor, if load controlled */
DOUBLE disval, dnorm, dinorm;
DOUBLE step_fac;

NR_CONTROLTYP controltyp;         /* type of path following technic */

static DIST_VECTOR  *re;          /* vector of out of balance loads */
static DIST_VECTOR  *rsd;         /* vector for incremental displacements */
static DIST_VECTOR  *dispi;       /* incrementel displacements */              

SOLVAR       *actsolv;            /* pointer to active solution structure */
PARTITION    *actpart;            /* pointer to active partition */
FIELD        *actfield;           /* pointer to active field */
INTRA        *actintra;           /* pointer to active intra-communicator */

SPARSE_TYP    array_typ;          /* type of sparse matrix */

static STANLN        nln_data;    /* structure to store and pass data for opt_stanln */

CALC_ACTION  *action;             /* pointer to the structures cal_action enum */

CONTAINER     container;          /* contains variables defined in container.h */
container.isdyn = 0;              /* static computation */

#ifdef DEBUG 
dstrc_enter("opt_stanln");
#endif
/*--------------- check wether this calculation is a restart or not ---*/
restart = genprob.restart;
/*----------------------------------------------------------------------*/
oflag=0;
/*------------ the distributed system matrix, which is used for solving */
/* 
NOTE: This routine only uses 1 global sparse matrix, which was created
      in global_mask_matrices. If You need more, it is necessary to allocate
      more matrices in actsolv->sysarray.
      Then one has to either copy the sparisty mask from 
      actsolv->actsysarray[0] to the others if using the same matrix format,
      or one has to call global mask matrices again for another type of 
      sparsity mask.
      Normally a control routine should not mix up different storage formats
      for system matrices
*/         
/*---- indize in the vector actsolv->sysarray[] of the stiffness matrix */
/*   the sparse matrix actsolv->sysarray[actsysarray] is used as system */
/*                                                               matrix */
/*              the sparsity type is actsolv->sysarray_typ[actsysarray] */
actsysarray=0;
/*--------------------------------------------------- set some pointers */
actfield    = &(field[0]);                      /* the structural field */
actsolv     = &(solv[0]);         /* the corresponednt SOLVAR structure */
actpart     = &(partition[0]);    /*     the partitioning of this field */
action      = &(calc_action[0]);  /*        the calculation action enum */
container.fieldtyp = actfield->fieldtyp;
#ifdef PARALLEL 
actintra    = &(par.intra[0]);    /*     the field's intra-communicator */
/*---------------- if we are not parallel, we have to allocate a pseudo */
/*                                         intra-communicator structure */
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of nonlinear statics, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*================================ start with preparing the calculation */
/*--------------------------------------------------------- free memory */
if(stalact == calsta_free) goto end;
/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(actsolv->sysarray[actsysarray],
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);
/*--------------------------------------------------- find control node */
/*---- cdof is the dof number of the degree of freedom to be controlled */
calstatserv_findcontroldof(actfield,
                           statvar->control_node_global,
                           statvar->control_dof,
                           &(statvar->controlnode),
                           &cdof); 
/*------------------------------------------------- get type of control */
/*              type of control algorithm (displacement, arclenght ...) */
controltyp = statvar->nr_controltyp;
/*------------------------------------------------- get number of steps */
nstep      = statvar->nstep; 
/*------------------------------------ get maximum number of iterations */
maxiter    = statvar->maxiter;
/*----------------------------------------------- number of rhs vectors */
/*------- the iteration uses rhs[0] for calculations and rhs[1] to hold */
/*----------------------------------------allocate 2 dist. load vectors */
actsolv->nrhs = 2;
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
}
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));
/*------------------------------------------ number of solution vectors */
/*-------------------- there is one solution vector to hold total displ.*/
actsolv->nsol= 1;
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
  for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));
}
/*------------------------ create vector re[0] for out of balance loads */
/*------------------ re[0] is used to hold residual forces in iteration */
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  solserv_create_vec(&(re),1,numeq_total,numeq,"DV");
}
solserv_zero_vec(&(re[0]));
/*--------------------- create 3 vectors rsd for residual displacements */
/*------------ the iteration uses rsd[0] to hold actual residual displ. */
/*--------------------- rsd[1] and rsd[2] are additional working arrays */
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  solserv_create_vec(&(rsd),3,numeq_total,numeq,"DV");
}
for (i=0; i<3; i++) solserv_zero_vec(&(rsd[i]));
/*------------------- create vector dispi for incremental displacements */
/*------------------ dispi[0] holds converged incremental displacements */
/*------------------ dispi[1] holds converged displacements last step   */
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  solserv_create_vec(&(dispi),2,numeq_total,numeq,"DV");
}
solserv_zero_vec(&(dispi[0]));
solserv_zero_vec(&(dispi[1]));
/*--------------------------------------------------- initialize solver */
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  init=1;
  solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[actsysarray]),
               &(actsolv->rhs[actsysarray]),
                init);
/*----------------------------- init the assembly for ONE sparse matrix */
  init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
/*------------------------------- init the element calculating routines */
  *action = calc_struct_init;
  calinit(actfield,actpart,action,&container);
}
/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0)
if (ioflags.struct_disp_gid||ioflags.struct_stress_gid) 
   out_gid_msh();
/*-------------------------------------- create the original rhs vector */
/*-------------------------- the approbiate action is set inside calrhs */
container.kstep = 0;
container.inherit = 1;
container.isdyn   = 0;              /* static computation */
container.actndis = 0;              /* only one discretisation */
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,actsysarray,
        &(actsolv->rhs[actsysarray]),action,&container);
/*--------------------------------------------------copy the rhs vector */
solserv_copy_vec(&(actsolv->rhs[actsysarray]),&(actsolv->rhs[actsysarray+1]));
/*----------------------------------------------------------------------*/
/*          The original rhs vector is now on actsolv->rhs[actsysarray] */
/*          AND in actsolv->rhs[actsysarray+1]                          */
/* the original load vector is held on [actsysarray+1]  all the time    */
/*----------------------------------------------------------------------*/
/*------------------ calculate euclidian vector norm of external forces */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[actsysarray+1]),&(nln_data.rinorm));
/*----------------------------------------------------------------------*/
/*--------------------------------------- set working variables to zero */
nln_data.sp1   = 0.0;
nln_data.csp   = 0.0;
nln_data.rlold = 0.0;
nln_data.rlnew = 0.0;
nln_data.rlpre = 0.0;
if(stalact == calsta_init || stalact==calsta_init_solve)
{
  amdef("arcfac",&(nln_data.arcfac),nstep,1,"DV");
}
amzero(&(nln_data.arcfac));
/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   /* colors the partitions in the postprocessing, if this is sequentiell*/
   /*              the picture of the partitions can be quite boring... */
   out_gid_domains(actfield);
}
/*------------------------------------------ initialialization finished */
if(stalact == calsta_init || stalact==calsta_init_solve) goto end;
/*======================================================================*/
/*                     START LOOP OVER ALL STEPS                        */
/*======================================================================*/
/*===================== predictor variables:============================*/
/* actsolv->rhs[actsysarray]      = Fext (will be changed)              */
/* actsolv->rhs[actsysarray+1]    = Fext (all the time)                 */
/* actsolv->sysarray[actsysarray] = Keug (nonl. Stiffness Matrix)       */
/* statvar->stepsize              = Delta l (control presc.deformation) */
/* rldiff                         = Delta mu                            */
/* nln_data->rlnew                = mu (actuel step)                    */
/* nln_data->rlold                = mu (last step)                      */
/* rsd[0]                         = first:Delta df / later:Delta d      */
/* dispi[0]                       = Delta d                             */
/* cdof                           = controled DOF                       */
/*===================== corrector variables:============================*/
/* actsolv->rhs[actsysarray]      = first:Fext / later:mu*Fext          */
/* actsolv->rhs[actsysarray+1]    = Fext (all the time)                 */
/* re[0]                          = Fext / mu*Fext / mu*Fext - Fint     */
/* actsolv->sysarray[actsysarray] = Keug (nonl. Stiffness Matrix)       */
/* intforce                       = Fint                                */
/* nln_data->rlnew                = mu (actuel step)                    */
/* nln_data->rlold                = mu (last step)                      */
/* rli                            = delta mu                            */
/* rsd[0]                         = delta d =delta d0+delta mu*delta df */
/* rsd[1]                         = delta df                            */
/* rsd[2]                         = delta d0                            */
/* dispi[0]                       = Delta d =Delta d + sum (delta d)    */
/* actsolv->sol[0]                = d act. (is as well known at nodes)  */
/* cdof                           = controled DOF*/
/*======================================================================*/
/* d=displacement, mu=load factor, Delta=incremental, delta= residual   */
/*======================================================================*/
/*--- load controlled ---*/
step_fac = curve[0].value.a.da[0][0]; /* first val. of curve for incr. */
act_loadfac = 0.;
/*
if(oflag==1)
{
   act_loadfac = (fin_loadfac-1)*step_fac; /* final value ! */ 
/*}
/*--- load controlled ---*/

for (kstep=0; kstep<nstep; kstep++)
{
   /*---------------------------------------------- write memory report */
   if (par.myrank==0) dsmemreport();
   /*---------- write report about all ARRAYs and ARRAY4Ds to .err file */
   /*dstrace_to_err();*/
   /*---- if this is a restart calc. get necessary info about last step */
   if (restart)
   {
   /*- tempor storage of new input (will be overwritten in rest._read_.)*/
   /*- (nstep and maxiter are already stored)                           */
     stepsi=statvar->stepsize;
     restartevery = statvar->resevery_restart;
   /*------------------------------------ read restart from pss- file --*/
     restart_read_nlnstructstat(restart,
                           statvar,
                           &nln_data,
                           actfield,
                           actpart,
                           actintra,
                           action,
                           actsolv->nrhs, actsolv->rhs,
                           actsolv->nsol, actsolv->sol,
                           1            , dispi,
                           &container);  /* contains variables defined in container.h */ 
   /*--------------- give the restartinfo to the associated variables --*/
     kstep = restart; 
   /*--------------------------------------------- switch restart off --*/
     restart = 0; 
   /*----------------------- give back temporary storage of new input --*/
     statvar->stepsize = stepsi;
     statvar->nstep = nstep; 
     statvar->maxiter = maxiter;
     statvar->resevery_restart = restartevery;
   }
                                            
   /*----------------------------------------------- LOAD-CONTROLLED ---*/
   if(controltyp==control_load)
   {
     /* copy original load vector from [actsysarray+1] to [actsysarray] */
     solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(actsolv->rhs[actsysarray]));
     /* get first load factor */
     act_loadfac += step_fac;
     
     nln_data.rlnew = act_loadfac;
     /*---------------------- if load controlled; multiply with factor  */
     solserv_scalarprod_vec(&(actsolv->rhs[actsysarray]),act_loadfac);
     /*------------- calculate euclidian vector norm of external forces */
     solserv_vecnorm_euclid(actintra,&(actsolv->rhs[actsysarray]),
                                                      &(nln_data.rinorm));
     nln_data.rrnorm = nln_data.rinorm;
     iconv=0;
     nlninc(
           actfield,
           actsolv,
           actpart,
           actintra,
           action,
           kstep,
           actsysarray,
           rsd,
           &nln_data,
           &iconv,
           &container
          );
     if(iconv==1)
     {
     /*--------------------------- update total displacements on sol[0] */
       solserv_copy_vec(&(rsd[0]),&(dispi[0]));
       /*-------------------------------- update total displacements on sol[0] */
       solserv_add_vec(&(dispi[0]),&(actsolv->sol[0]),1.0);
       /*---------- put actsolv->sol[kstep] to the nodes (total displacements) */ 
       solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
       /*----------- printf out iteration heading to err and to shell */
       if (actintra->intra_rank==0)
                      conequ_printhead(kstep,controltyp,cdof,nln_data.csp);
       solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sol[0]),
                   cdof,
                   &disval);
       if (actintra->intra_rank==0)
       /*---------------------------- norm of incremental displacements */
       solserv_vecnorm_euclid(actintra,&(dispi[0]),&dnorm);
       /*------------------------------- norm of residual displacements */
       solserv_vecnorm_euclid(actintra,&(rsd[0]),&dinorm);

       conequ_printiter(kstep,
                        disval,
                        nln_data.rlnew,
                        dinorm,
                        nln_data.renorm,
                        nln_data.renergy,
                        dnorm,
                        nln_data.rrnorm);
     }
     else
     {
   /*-------------------------------------- make equillibrium iteration */  
     nlnequ(
           actfield,
           actsolv,
           actpart,
           actintra,
           action,
           kstep,
          &itnum,
           actsysarray,
           rsd,
           dispi,
           re,
           cdof,
          &nln_data,
           controltyp,
           oflag,
           &container
         );    
     }
   }
   /*--------------------------------------- DISPLACEMENT-CONTROLLED ---*/
   if(controltyp==control_disp)
   {
   /*--------------------------------------------------- make predictor */
   conpre(actfield,
          actsolv,
          actpart,
          actintra,
          action,
          kstep,
          actsysarray,
          rsd,
          dispi,
          cdof,
         &nln_data,
          controltyp,
         &container);
   /*-------------------------------------- make equillibrium iteration */  
   conequ(actfield,
          actsolv,
          actpart,
          actintra,
          action,
          kstep,
          &itnum,
          actsysarray,
          rsd,
          dispi,
          re,
          cdof,
          &nln_data,
          controltyp,
          &container);    
   }
   /*-- update for nonlinear material models - new stress/strain values */  
   *action = calc_struct_update_istep;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.kstep        = 0;
   container.inherit = 1;
   container.isdyn   = 0;              /* static computation */
   container.actndis = 0;              /* only one discretisation */
   calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
   /*-------------- check, if output is to be written is this loadstep */
   mod_displ =    kstep  % statvar->resevry_disp;        
   mod_stress =   kstep  % statvar->resevry_stress;
   mod_restart =  kstep  % statvar->resevery_restart;
   /*-------------------------------------- perform stress calculation */
   if (mod_stress==0 || mod_displ==0)
   {
      if (ioflags.struct_stress_file==1 || ioflags.struct_stress_gid==1)
      {
         *action = calc_struct_stress;
         container.dvec         = NULL;
         container.dirich       = NULL;
         container.global_numeq = 0;
         container.kstep        = 0;
         container.inherit = 1;
         container.isdyn   = 0;              /* static computation */
         container.actndis = 0;              /* only one discretisation */
         calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,&container,action);
         /*---------------------- reduce stresses, so they can be written */
         /*- this makes the result of the stress calculation redundant on */
         /*                                                         procs */
         *action = calc_struct_stressreduce;
         container.kstep = 0;
         calreduce(actfield,actpart,actintra,action,&container);
      }
   }
   /*--------------------------------------- print out results to .out */
   if (mod_stress==0 || mod_displ==0)
   {
      if (ioflags.struct_stress_file==1 && ioflags.struct_disp_file==1)
      {
        out_sol(actfield,actpart,actintra,kstep,0);
      }
   }
   /*----------------------------------------- printout results to gid */
   if (par.myrank==0) 
   {
       if (ioflags.struct_disp_gid==1 && mod_displ==0)
       out_gid_sol("displacement",actfield,actintra,kstep,0);
       if (ioflags.struct_stress_gid==1 && mod_stress==0)
       out_gid_sol("stress"      ,actfield,actintra,kstep,0);
   }
   /*----------------------------------- printout results to pss file */
   if(mod_restart==0)
  { 
   switch(controltyp)
    {
    case control_disp:                        /* displacement control */
      restart_write_nlnstructstat(statvar,
                           &nln_data,
                           actfield,
                           actpart,
                           actintra,
                           action,
                           kstep,
                           actsolv->nrhs, actsolv->rhs,
                           actsolv->nsol, actsolv->sol,
                           1            , dispi,
                           &container);  /* contains variables defined in container.h */ 
    break;
    case control_arc:                            /* arclenght control */
      dserror("restart for arclengh not yet impl.");
    break;
    case control_load:

    break;
    case control_none:
      dserror("Unknown typ of path following technique");
    break;
    default:
      dserror("Unknown typ of path following technique");
    break;
    }
  }
} /* end of (kstep=0; kstep<nstep; kstep++) */
/*----------------------------------------------------------------------*/
opt->nlndat = nln_data.rlnew;    /* store load factor for optimization */ 
/*----------------------------------------------------------------------*/
end:
/*---------------------------------------------------- make cleaning up */
if(stalact == calsta_free)
{
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);
  solserv_del_vec(&(re),1);
  solserv_del_vec(&(rsd),3);
  solserv_del_vec(&(dispi),2);
}
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
/*----------------------------------------------------------------------*/
return;
} /* end of opt_stanln */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  PERFORM LINEAR PREDICTOR STEP                               al 11/01|
 *----------------------------------------------------------------------*/
void nlninc(
            FIELD         *actfield,     /* the actual physical field */
            SOLVAR        *actsolv,      /* the field-corresponding solver */
            PARTITION     *actpart,      /* the partition of the proc */
            INTRA         *actintra,     /* the intra-communicator of this field */
            CALC_ACTION   *action,       /* calculation flag */
            INT            kstep,        /* the load or time step we are in */
            INT            actsysarray,  /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,          /* dist. vector of incremental residual forces used for iteration in conequ */
            STANLN        *nln_data,     /* data of the Newton-Raphson method */
            INT           *iconv,        /* type of control algorithm */
            CONTAINER     *container 	 /* contains variables defined in container.h */
          )
{
INT                  init;             /* solver flag */
DOUBLE               dnorm;
DOUBLE               dinorm;

ARRAY                intforce_a;             /* global redundant vector of internal forces */
DOUBLE              *intforce;               /* pointer to intforce_a.a.dv */

#ifdef DEBUG 
dstrc_enter("nlninc");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------- create an array for the internal forces */
intforce = amdef("intforc1",&intforce_a,actsolv->sol[0].numeq_total,1,"DV");
/*------------------------------- initialize vector for internal forces */
amzero(&intforce_a);
/*--------------------------------- init the dist sparse matrix to zero */
/*                  NOTE: Has to be called after solver_control(init=1) */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );
/*----------------------- calculate tangential stiffness in actsysarray */
*action = calc_struct_nlnstiff;
container->inherit = 1;
container->isdyn   = 0;              /* static computation */
container->actndis = 0;              /* only one discretisation */
container->dvec         = intforce;
container->dirich       = NULL;
container->global_numeq = actsolv->sol[0].numeq_total;
container->kstep        = kstep;
calelm(actfield,        /* active field                          */
       actsolv,         /* active solver typ                     */
       actpart,         /* my partition of this field            */
       actintra,        /* my intra-comunicators                 */
       actsysarray,     /* system-stiffness matrix               */
       -1,              /*system-mass matrix (there is no)       */
       container,       /* contains variables defined in container.h */
       action);         /* what to do                            */

/* add internal forces to scaled external forces to get residual forces */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[actsysarray]),
             intforce,
             -1.0
             );
/*------------- calculate euclidian vector norm of external forces */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[actsysarray]),
                                                 &(nln_data->rinorm));
nln_data->rrnorm = nln_data->rinorm;
/*------------------------------------------ solve for incremental load */
init=0;
solver_control(  
                 actsolv,
                 actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(rsd[0]),
               &(actsolv->rhs[actsysarray]),
                 init
              );

/*----------------------- topic: check for concergence in displacements */
/*-------------------------------------- norm of residual displacements */
solserv_vecnorm_euclid(actintra,&(rsd[0]),&dinorm);
/*----------------------------------- norm of incremental displacements */
dnorm = dinorm;
*iconv = 0;
if (dinorm < statvar->toldisp) *iconv = 1;
/*----------------------------------------------------------------------*/
amdel(&intforce_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of nlninc */
/*----------------------------------------------------------------------*
 |  PERFORM EQUILLIBRIUM ITERATION                           m.gee 11/01|
 |  within Newton Raphson                                               |
 *----------------------------------------------------------------------*/
void nlnequ(
            FIELD         *actfield,      /* the actual physical field */
            SOLVAR        *actsolv,       /* the field-corresponding solver */
            PARTITION     *actpart,       /* the partition of the proc */
            INTRA         *actintra,      /* the intra-communicator of this field */
            CALC_ACTION   *action,        /* calculation flag */
            INT            kstep,         /* the load or time step we are in */
            INT           *itnum,         /* number of corrector steps taken by this routine */
            INT            actsysarray,   /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,           /* dist. vector of incremental residual forces used for iteration in conequ */
            DIST_VECTOR   *dispi,         /* dist. vector of incremental displacements */
            DIST_VECTOR   *re,            /* re[0..2] 3 vectors for residual displacements */
            INT            cdof,          /* number of dof to be controlled */
            STANLN        *nln_data,      /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,    /* type of control algorithm */
            INT            oflag,         /* type of control algorithm */
            CONTAINER     *container      /* contains variables defined in container.h */
          )
{
INT                  init;                   /* init flag for solver */
INT                  itemax;                 /* max. number of corrector steps allowed */
DOUBLE               stepsize;               /* stepsize */
INT                  convergence=0;          /* test flag for convergence */
DOUBLE               rl0;                    /* some norms and working variables */                   
DOUBLE               rlold;
DOUBLE               rlnew;
DOUBLE               rlpre;
DOUBLE               renorm=0.0;
DOUBLE               energy=0.0;
DOUBLE               dnorm;
DOUBLE               dinorm;
DOUBLE               told;
DOUBLE               disval;

static DOUBLE dsx,dsi,alfa0;
static DOUBLE rni, rn1, sp, sp1; 
DOUBLE spi;

ARRAY                intforce_a;             /* global redundant vector of internal forces */
DOUBLE              *intforce;               /* pointer to intforce_a.a.dv */

#ifdef DEBUG 
dstrc_enter("nlnequ");
#endif
/*----------------------------------------------------------------------*/
*itnum=0;
if (kstep == 0)
{
  dsx=0.;
  if(oflag!=1) rn1=0.;
}
/*----------------------------------- get load factor of last increment */
if (kstep != 0) rl0 = nln_data->arcfac.a.dv[kstep-1];
else            rl0 = 0.0;
nln_data->rlold     = rl0;
rlold               = nln_data->rlold;
rlnew               = nln_data->rlnew;
rlpre               = nln_data->rlpre;
stepsize            = statvar->stepsize;
itemax              = statvar->maxiter;
/*-------------------------------- update total displacements on sol[0] */
solserv_copy_vec(&(rsd[0]),&(dispi[0]));
/*----------------------------------------------- store old load factor */
rlold               = rlnew;
nln_data->rlold     = rlold;
/*-------------------------------- update total displacements on sol[0] */
solserv_add_vec(&(dispi[0]),&(actsolv->sol[0]),1.0);
/*---------- put actsolv->sol[kstep] to the nodes (total displacements) */ 
solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*----------- printf out iteration heading to err and to shell */
if (actintra->intra_rank==0) conequ_printhead(kstep,controltyp,cdof,nln_data->csp);
/*--------------------- make iteration parameter before first iteration */
/*----------------------------------- norm of incremental displacements */
solserv_vecnorm_euclid(actintra,&(dispi[0]),&dnorm);
/*-------------------------------------- norm of residual displacements */
solserv_vecnorm_euclid(actintra,&(rsd[0]),&dinorm);
/*------------------------------------- total disp value of control dof */ 
solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(dispi[0]),
                   cdof,
                   &disval);
/*---  evaluate incremental arc-length "dsi" */
/*---  evaluate first scaling factor "alfa0 = dsi/ds0" */
/*---  itiarc */
  dsi=disval;
  if (dsx==0.) dsx=dsi;
  alfa0 = dsx/dsi;

/*---------------------------------- create current stiffness parameter */
/*--- itcsp ---*/
solserv_dot_vec(actintra,&(actsolv->rhs[actsysarray]),&(rsd[0]),&spi);
/*????????*/
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[actsysarray]),&(rni));


/*------------------------------------- initialize values if zero */
if (ABS(rn1) <= EPS14) rn1=rni;
if (ABS(nln_data->sp1) <= EPS14) nln_data->sp1 = spi;

/*------------------------------- evaluate current stiffness parameterr */
if (statvar->iarc==0) 
{
  sp=nln_data->sp1/spi*(rni/rn1);
  sp*=sp;
}
/*------------------------------------ save current stiffness parameter */
nln_data->csp = sp;


/*----------------------------- create an array for the internal forces */
intforce = amdef("intforce",&intforce_a,actsolv->sol[0].numeq_total,1,"DV");
/*======================================================================*/
/*                        start of the iteration loop of this increment */
/*======================================================================*/
while (*itnum < itemax)
{
/*---------------------------- copy the initial rhs to rhs[actsysarray] */
solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(actsolv->rhs[actsysarray]));
/*------------------------------- scale the rhs[actsysarray] with rlnew */
solserv_scalarprod_vec(&(actsolv->rhs[actsysarray]),rlnew);
/*--------------------------------------- copy the initial rhs to re[0] */
solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(re[0]));
/*------------------------------------------ scale the re[0] with rlnew */
solserv_scalarprod_vec(&(re[0]),rlnew);
/* put residual displacements to the nodes (needed for material, eas ..)*/              
solserv_result_resid(
                     actfield,
                     actintra,
                     &(rsd[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*-- put incremental displacements to the nodes (needed for material...)*/              
solserv_result_incre(
                     actfield,
                     actintra,
                     &(dispi[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*---------------------------------------------initialize system matrix */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );
/*------------------------------- initialize vector for internal forces */
amzero(&intforce_a);
/*----------------------- calculate tangential stiffness in actsysarray */
/*----- calculate the rhs resulting from dirichlet conditions in dirich */
*action = calc_struct_nlnstiff;
container->inherit = 1;
container->isdyn   = 0;              /* static computation */
container->actndis = 0;              /* only one discretisation */
container->dvec         = intforce;
container->dirich       = NULL;
container->global_numeq = actsolv->sol[0].numeq_total;
container->kstep        = kstep;
calelm(actfield,        /* active field                          */
       actsolv,         /* active solver typ                     */
       actpart,         /* my partition of this field            */
       actintra,        /* my intra-comunicators                 */
       actsysarray,     /* system-stiffness matrix               */
       -1,              /*system-mass matrix (there is no)       */
       container,       /* contains variables defined in container.h */
       action);         /* what to do                            */
/* add internal forces to scaled external forces to get residual forces */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(re[0]),
             intforce,
             -1.0
             );
/*-------------- solve for out-of balance loads, put solution to rsd[2] */
/*                                                         K * du2 = -R */
/*                                                initial guess is zero */
init=0;
solserv_zero_vec(&(rsd[2]));
solver_control(  
                 actsolv,
                 actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(rsd[2]),
               &(re[0]),
                 init
              );
/*------------------------------------------- calculate residual energy */
solserv_dot_vec(actintra,&(rsd[2]),&(re[0]),&energy);
nln_data->renergy = energy;
/*------------------------------ calculate norm of out-of-balance-loads */
solserv_vecnorm_euclid(actintra,&(re[0]),&renorm);
nln_data->renorm = renorm;
/*------------------------------------ update of load and displacements */
solserv_add_vec(&(rsd[2]),&(actsolv->sol[0]),1.0);
/*----------------------------- put actual total displacements to nodes */
solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*------------------------------------- calculate norm of displacements */
solserv_vecnorm_euclid(actintra,&(dispi[0]),&dnorm);
/*----------------- calculate norm of updated incremental displacements */
solserv_vecnorm_euclid(actintra,&(rsd[2]),&dinorm);
/*----------------------------------------------- check for convergence */
told=0.0;
if (dinorm > statvar->toldisp) told = dinorm/dnorm;
if (told   > statvar->toldisp) convergence=0;
else                        convergence=1;
/*------------------------------------------------------- make printout */
solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sol[0]),
                   cdof,
                   &disval);
if (actintra->intra_rank==0)
conequ_printiter(*itnum,disval,rlnew,dinorm,renorm,energy,
                                                  dnorm,nln_data->rrnorm);
/*----------------------------- decide for or against another iteration */
if (convergence)/*------------------------------------------convergence */
{
   nln_data->arcfac.a.dv[kstep] = rlnew;
   nln_data->rlnew              = rlnew;
   break;
}
else/*---------------------------------------------------no convergence */
{
   if (*itnum < itemax) *itnum = *itnum + 1;
   else 
   {
      if (actintra->intra_rank==0)
      {
         printf("WARNING: No convergence in global NR in maxiter steps! continue....\n");
         fprintf(allfiles.out_err,"WARNING: No convergence in global NR in maxiter steps! continue....\n");
      }
      break;
   }
}
} /*--------------------------------------------- end of iteration loop */
/*---------------------------------- update data after incremental step */           
/*--------------------- copy converged solution from sol[0] to dispi[1] */
solserv_copy_vec(&(actsolv->sol[0]),&(dispi[1]));
/*----------------- put converged solution to the elements in next step */
solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*----------------------------------------------- update of load factor */                                  
nln_data->rlold = nln_data->rlnew;
/*--------------- update of parameters for material law and eas etc.... */
/* nothing needed yet..... */
/*----------------------------------------------------------------------*/
amdel(&intforce_a);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of nlnequ */
/*----------------------------------------------------------------------*/
#endif /* stop including optimization code to ccarat :*/
/*----------------------------------------------------------------------*/

/*! @} (documentation module close)*/

