/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../io/io.h"
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

/*----------------------------------------------------------------------*
 |  routine to control nonlinear static execution       m.gee 11/01     |
 *----------------------------------------------------------------------*/
void stanln()
{
INT           i;                  /* simply a counter */
INT           numeq;              /* number of equations on this proc */
INT           numeq_total;        /* total number of equations */
INT           init;               /* flag for solver_control call */
INT           actsysarray;        /* number of the active system sparse matrix */
INT           cdof;               /* the Id of the controlled dof */

INT           reldof[6];          /* the Id's of the dofs for output */

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

NR_CONTROLTYP controltyp;         /* type of path following technic */

DIST_VECTOR  *re;                 /* vector of out of balance loads */
DIST_VECTOR  *rsd;                /* vector for incremental displacements */
DIST_VECTOR  *dispi;              /* incrementel displacements */

SOLVAR       *actsolv;            /* pointer to active solution structure */
PARTITION    *actpart;            /* pointer to active partition */
FIELD        *actfield;           /* pointer to active field */
INTRA        *actintra;           /* pointer to active intra-communicator */

SPARSE_TYP    array_typ;          /* type of sparse matrix */

STANLN        nln_data;           /* structure to store and pass data for stanln */

CALC_ACTION  *action;             /* pointer to the structures cal_action enum */

#ifdef BINIO
BIN_OUT_FIELD out_context;
#endif

CONTAINER     container;          /* contains variables defined in container.h */
container.isdyn   = 0;              /* static computation */
container.actndis = 0;              /* only one discretisation */

#ifdef DEBUG
dstrc_enter("stanln");
#endif
/*--------------- check wether this calculation is a restart or not ---*/
restart = genprob.restart;
/*----------------------------------------------------------------------*/
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
/*------------------------------------ put kintyp to container sh 03/03 */
switch(statvar->kintyp)
{
case geo_linear:
   container.kintyp = 0;
break;
case upd_lagr:
   container.kintyp = 1;
break;
case tot_lagr:
   container.kintyp = 2;
break;
default:
   dserror("Unknown typ of kinematic");
break;
}

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
/*------------------------------------------------ typ of global matrix */
array_typ   = actsolv->sysarray_typ[actsysarray];
/*---------------------------- get global and local number of equations */
/*   numeq equations are on this proc, the total number of equations is */
/*                                                          numeq_total */
solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
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
/*--------------- reldof are dof numbers wanted for output  to err file */
if(ioflags.relative_displ>0)
{
  calstatserv_findreldofs(actfield,
                          statvar->reldisnode_ID,
                          statvar->reldis_dof,
                          ioflags.relative_displ,
                          reldof);
}
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
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));
/*------------------------------------------ number of solution vectors */
/*-------------------- there is one solution vector to hold total displ.*/
actsolv->nsol= 1;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));
/*------------------------ create vector re[0] for out of balance loads */
/*------------------ re[0] is used to hold residual forces in iteration */
/*----- re[1] holds up-to-date residual forces for check of convergence */
solserv_create_vec(&(re),2,numeq_total,numeq,"DV");

solserv_zero_vec(&(re[0]));
solserv_zero_vec(&(re[1]));
/*--------------------- create 3 vectors rsd for residual displacements */
/*------------ the iteration uses rsd[0] to hold actual residual displ. */
/*--------------------- rsd[1] and rsd[2] are additional working arrays */
solserv_create_vec(&(rsd),3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(rsd[i]));
/*------------------- create vector dispi for incremental displacements */
/*------------------ dispi[0] holds converged incremental displacements */
/*------------------ dispi[1] holds converged displacements last step   */
solserv_create_vec(&(dispi),2,numeq_total,numeq,"DV");
solserv_zero_vec(&(dispi[0]));
solserv_zero_vec(&(dispi[1]));
/*--------------------------------------------------- initialize solver */
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

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   &(actsolv->sysarray_typ[actsysarray]), &(actsolv->sysarray[actsysarray]),
                   actfield, actpart, actintra, 0);
#endif

/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0)
if (ioflags.struct_disp_gid||ioflags.struct_stress_gid)
   out_gid_msh();
/*-------------------------------------- write output of submesh to gid */
#ifdef D_MLSTRUCT
if (genprob.multisc_struct == 1 && (ioflags.struct_sm_disp_gid ||ioflags.struct_sm_stress_gid ))
{
   out_gid_smsol_init();
   out_gid_submesh();
}
#endif
/*-------------------------------------- create the original rhs vector */
/*-------------------------- set action before of calrhs --- genk 10/02 */
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
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
amdef("arcfac",&(nln_data.arcfac),nstep,1,"DV");
amzero(&(nln_data.arcfac));
/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0)
{
   /* colors the partitions in the postprocessing, if this is sequentiell*/
   /*              the picture of the partitions can be quite boring... */
   out_gid_domains(actfield);
}
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

for (kstep=0; kstep<nstep; kstep++)
{
   /*---------------------------------------------- write memory report */
   if (par.myrank==0) dsmemreport();
   /*---------------------- if stepsize is variable get actual stepsize */
   if (statvar->isrelstepsize==1)
   {
     get_stepsize(kstep,statvar);
   }
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
#if defined(BINIO) && defined(NEW_RESTART_READ)
     restart_read_bin_nlnstructstat(statvar,
                                    &nln_data,
                                    &array_typ,
                                    &(actsolv->sysarray[actsysarray]),
                                    actfield,
                                    actpart,
                                    0,
                                    actintra,
                                    actsolv->nrhs, actsolv->rhs,
                                    actsolv->nsol, actsolv->sol,
                                    1            , dispi,
                                    restart);
#else
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
#endif
   /*--------------- give the restartinfo to the associated variables --*/
   /* kstep = restart; --*/
     kstep = restart+1;
   /*--------------------------------------------- switch restart off --*/
     restart = 0;
   /*----------------------- give back temporary storage of new input --*/
     statvar->stepsize = stepsi;
     statvar->nstep = nstep;
     statvar->maxiter = maxiter;
     statvar->resevery_restart = restartevery;
   }
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
          reldof,
          &nln_data,
          controltyp,
          &container);
   /*-- update for nonlinear material models - new stress/strain values */
   *action = calc_struct_update_istep;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.kstep        = 0;
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
       out_gid_sol("displacement",actfield,actintra,kstep,0,ZERO);
       if (ioflags.struct_stress_gid==1 && mod_stress==0)
       out_gid_sol("stress"      ,actfield,actintra,kstep,0,ZERO);
   }
   /*------------------------------ printout multiscale results to gid */
#ifdef D_MLSTRUCT
   if (genprob.multisc_struct && ioflags.struct_sm_disp_gid && mod_displ==0)
     out_gid_smdisp("displacement",kstep);
   if (genprob.multisc_struct && ioflags.struct_sm_stress_gid && mod_stress==0)
     out_gid_smstress("stress",kstep); 
#endif

#ifdef BINIO
   if (ioflags.struct_disp_gid==1 && mod_displ==0)
     out_results(&out_context, 0, kstep, 0, OUTPUT_DISPLACEMENT);

   if (ioflags.struct_stress_gid==1 && mod_stress==0)
     out_results(&out_context, 0, kstep, 0, OUTPUT_STRESS);
#endif

   /*----------------------------------- printout results to pss file */
   if(mod_restart==0 && genprob.restart != 0)
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
#ifdef BINIO
      restart_write_bin_nlnstructstat(&out_context,
                                      statvar,
                                      &nln_data,
                                      kstep,
                                      actsolv->nrhs, actsolv->rhs,
                                      actsolv->nsol, actsolv->sol,
                                      1            , dispi);
#endif
    break;
    case control_arc:                            /* arclenght control */
      dserror("restart for arclengh not yet impl.");
    break;
    case control_load:
      dserror("load control not yet impl.");
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
/*---------------------------------------------------- make cleaning up */
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&(re),1);
solserv_del_vec(&(rsd),3);
solserv_del_vec(&(dispi),2);
/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL
CCAFREE(actintra);
#endif
amdel(&(nln_data.arcfac));
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of stanln */








/*----------------------------------------------------------------------*
 |  PERFORM LINEAR PREDICTOR STEP                            m.gee 11/01|
 *----------------------------------------------------------------------*/
void conpre(  /*!< rueckgabe INT da kein prototyp definiert wurde! */
            FIELD         *actfield,     /* the actual physical field */
            SOLVAR        *actsolv,      /* the field-corresponding solver */
            PARTITION     *actpart,      /* the partition of the proc */
            INTRA         *actintra,     /* the intra-communicator of this field */
            CALC_ACTION   *action,       /* calculation flag */
            INT            kstep,        /* the load or time step we are in */
            INT            actsysarray,  /* number of the system matrix in actsolv->sysarray[actsysarray] to be used */
            DIST_VECTOR   *rsd,          /* dist. vector of incremental residual displ. used for iteration in conequ */
            DIST_VECTOR   *dispi,        /* dist. vector of incremental displacements */
            INT            cdof,         /* number of the dof to be controlled */
            STANLN        *nln_data,     /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,   /* type of control algorithm */
            CONTAINER     *container 	 /* contains variables defined in container.h */
           )
{
INT                  init;             /* solver flag */
DOUBLE               prenorm;          /* norm of predictor step */
DOUBLE               controldisp;      /* displacment value at controled dof */
DOUBLE               rldiff;           /* forgot.... */
DOUBLE               spi;              /* forgot.... */


#ifdef DEBUG
dstrc_enter("conpre");
#endif
/*----------------------------------------------------------------------*/
statvar->praedictor = 1;
/*--------------------------------- init the dist sparse matrix to zero */
/*                  NOTE: Has to be called after solver_control(init=1) */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );
/*----------------------- calculate tangential stiffness in actsysarray */
/*----- calculate the rhs resulting from dirichlet conditions in dirich */
*action = calc_struct_nlnstiff;
container->dvec         = NULL;
container->dirich       = NULL;
container->global_numeq = 0;
container->kstep        = kstep;
calelm(actfield,        /* active field                          */
       actsolv,         /* active solver typ                     */
       actpart,         /* my partition of this field            */
       actintra,        /* my intra-comunicators                 */
       actsysarray,     /* system-stiffness matrix               */
       -1,              /*system-mass matrix (there is no)       */
       container,       /* contains variables defined in container.h */
       action);         /* what to do                            */
/*----- copy original load vector from [actsysarray+1] to [actsysarray] */
solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(actsolv->rhs[actsysarray]));
/*---------- LGS: K*u=F -------------------- solve for incremental load */
init=0;
/*
oll_print(actsolv->sysarray[actsysarray].oll,200);
*/
solver_control(
                 actsolv,                             /* active solver typ        */
                 actintra,                            /* my intra-comunicators    */
               &(actsolv->sysarray_typ[actsysarray]), /* Systemmatrixtyp          */
               &(actsolv->sysarray[actsysarray]),     /* Systemstema K of LGS     */
               &(rsd[0]),                             /* solution vector u of LGS */
               &(actsolv->rhs[actsysarray]),          /* RHS F of LGS             */
                 init                                 /* option                   */
              );
/*--------------------------- calculate euclidian norm of displacements */
solserv_vecnorm_euclid(actintra,&(rsd[0]),&prenorm);
/*--------------------------------------- do scaling for load parameter */
switch(controltyp)
{
case control_disp:
/*---------------------------------- get control dof from distr. vector */
   solserv_getele_vec( actintra,                         /* my intra-comunicators */
                      &(actsolv->sysarray_typ[actsysarray]),    /* sparsity typ   */
                      &(actsolv->sysarray[actsysarray]),        /* sparse matrix  */
                      &(rsd[0]),         /* vector the value shall be taken from  */
                       cdof,             /* field-local (unsupported) dof number  */
                      &controldisp);     /* value in vector at the given dof      */
   rldiff = (statvar->stepsize)/controldisp;
break;
case control_arc:
   if (statvar->iarc==0) statvar->arcscl=0.0;
   else if (kstep==0 && statvar->iarc==1)
   {
      statvar->arcscl = prenorm / sqrt((DOUBLE)(actsolv->sol[0].numeq_total));
   }
   rldiff = (statvar->stepsize) / (sqrt(prenorm*prenorm + statvar->arcscl*statvar->arcscl));
break;
case control_none:
   rldiff = 0.0;
   dserror("Unknown typ of path following control");
break;
default:
   rldiff = 0.0;
   dserror("Unknown typ of path following control");
break;
}
/*---------------------------------- create current stiffness parameter */
solserv_dot_vec(actintra,                    /* my intra-comunicators           */
                &(actsolv->rhs[actsysarray]),/* first vector a to be multiplied */
                &(rsd[0]),                   /* secnd vector b to be multiplied */
                &spi);                       /* result a*b                      */
if (ABS(nln_data->sp1) <= EPS14) nln_data->sp1 = spi;
if (ABS(spi)           <= EPS14) nln_data->csp = 1.0;
else                             nln_data->csp = nln_data->sp1 / spi;
/*----------- use current stiffness parameter to check for limit points */
/*
sp1 = gradient of controlled load path in first increment
spi = gradient of controlled load path in current increment

csp =  sp1/spi is a relation of the current gradient of the load curve to
the initial gradient in the first step, so
if csp turns to zero, a limit or saddle point is reached, and if
csp takes negative values we are on a branch which is going downwards
(neg. gradient)
Dependent on the type of constraint that it used (e.g. displacement/arclenght control)
there is no sign to the predictor, so the predictor is always positive,
that means pointing upwards. To get a better predictor on a descending path
rldiff takes the sign of csp if flag signchcsp is set:
*/
if (statvar->signchcsp==1)
{
  rldiff *= nln_data->csp/ABS(nln_data->csp);
}
/*--------------------------------------------------------- save values */
nln_data->rlnew = nln_data->rlold + rldiff;
/*---------- create the correct displacmements after predictor solution */
/* displacements of increment are stored in dispi[0]                    */
/* displacements of iteration are in rsd[0]                             */
solserv_scalarprod_vec(&(rsd[0]),rldiff);
solserv_copy_vec(&(rsd[0]),&(dispi[0]));
/*---------------------------------------- make new norm of load vector */
nln_data->rrnorm = nln_data->rinorm * nln_data->rlnew;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of conpre */


/*----------------------------------------------------------------------*
 |  PERFORM EQUILLIBRIUM ITERATION                           m.gee 11/01|
 |  within Newton Raphson                                               |
 *----------------------------------------------------------------------*/
void conequ( /*!< rueckgabe INT da kein prototyp definiert wurde! */
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
            INT           *reldof,        /* numbers of dofs for output */
            STANLN        *nln_data,      /* data of the Newton-Raphson method */
            NR_CONTROLTYP  controltyp,    /* type of control algorithm */
            CONTAINER     *container	  /* contains variables defined in container.h */
           )
{
INT                  init;                   /* init flag for solver */
INT                  itemax;                 /* max. number of corrector steps allowed */
DOUBLE               stepsize;               /* stepsize */
INT                  convergence=0;          /* test flag for convergence */
INT                  j;
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
DOUBLE               rli;

INT                  rel_dof;
DOUBLE               reldisval[6];
DOUBLE               reldis;
DOUBLE               resnorm_uptodate=0.0;


ARRAY                intforce_a;             /* global redundant vector of internal forces */
DOUBLE              *intforce;               /* pointer to intforce_a.a.dv */
#ifdef PARALLEL
static ARRAY intforcerecv_a;
static DOUBLE *intforcerecv;
#endif


#ifdef DEBUG
dstrc_enter("conequ");
#endif
/*----------------------------------------------------------------------*/
*itnum=0;
/*----------------------------------------------------------------------*/
statvar->praedictor = 2;
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
                   &(actsolv->sol[0]),
                   cdof,
                   &disval);
/*------------------------------------------- print out first iteration */
if (actintra->intra_rank==0)
conequ_printiter(*itnum,disval,rlnew,dinorm,renorm,energy,dnorm,nln_data->rrnorm);
/*----------------------------------------------------------------------*/
/*----------------------------- create an array for the internal forces */
intforce = amdef("intforce",&intforce_a,actsolv->sol[0].numeq_total,1,"DV");
#ifdef PARALLEL
if (intforcerecv_a.Typ != cca_DV)
{
   intforcerecv = amdef("tmpintf",&intforcerecv_a,actsolv->sol[0].numeq_total,1,"DV");
}
if (intforcerecv_a.fdim < actsolv->sol[0].numeq_total)
{
   amdel(&intforcerecv_a);
   intforcerecv = amdef("tmpintforce",&intforcerecv_a,actsolv->sol[0].numeq_total,1,"DV");
}
amzero(&intforcerecv_a);
#endif
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
                     &(actsolv->sysarray_typ[actsysarray]),
		     actsysarray
		     );
/*---------------------------------------------initialize system matrix */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );
/*------------------------------- initialize vector for internal forces */
amzero(&intforce_a);
/*------------------------- calculate new stiffness and internal forces */
*action = calc_struct_nlnstiff;
container->dvec         = intforce;
container->dirich       = NULL;
container->global_numeq = actsolv->sol[0].numeq_total;
container->kstep        = kstep;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,container,action);
/*----------------------------------------------------------------------*/
statvar->praedictor = 3;
/*
oll_print(actsolv->sysarray[actsysarray].oll,7);
*/
/*--------------------------------------- allreduce the vector intforce */
#ifdef PARALLEL
amzero(&intforcerecv_a);
MPI_Allreduce(intforce,intforcerecv,actsolv->sol[0].numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);

/* add internal forces to scaled external forces to get residual forces */
assemble_vec(actintra, &(actsolv->sysarray_typ[actsysarray]), &(actsolv->sysarray[actsysarray]), &(re[0]), intforcerecv, -1.0 );
#else
/* add internal forces to scaled external forces to get residual forces */
assemble_vec(actintra, &(actsolv->sysarray_typ[actsysarray]), &(actsolv->sysarray[actsysarray]), &(re[0]), intforce, -1.0 );
#endif
/*-------------- solve for out-of balance loads, put solution to rsd[2] */
/*                                                         K * du2 = -R */
/*                                                initial guess is zero */
init=0;
solserv_zero_vec(&(rsd[2]));
solver_control(  actsolv,
                 actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(rsd[2]),
               &(re[0]),
                 init);
/*-------------------------------------- solve for original load vector */
/*                                                          K * du1 = P */
/*                                 initial guess is value of last solve */
init=0;
solver_control(  actsolv,
                 actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(rsd[1]),
               &(actsolv->rhs[actsysarray+1]),
                 init);
/*===============================make increment of load and displacment */
switch(controltyp)
{
case control_disp:/*===============================displacement control */
   increment_controldisp(actintra,actsolv,actsysarray,cdof,rsd,dispi,&rli);
break;
case control_arc:/*=================================== arlenght control */
   increment_controlarc(actintra,actsolv,actsysarray,rsd,dispi,&(actsolv->sol[0]),
                        rlnew,rlold,stepsize,&rli);
break;
case control_load:
   dserror("load control not yet impl.");
break;
case control_none:
   dserror("Unknown typ of path following technique");
break;
default:
   dserror("Unknown typ of path following technique");
break;
}
/*----------------------------------------- update of load factor rlnew */
rlnew += rli;
nln_data->rlnew = rlnew;
/*------------------------------------ update of load and displacements */
solserv_add_vec(&(rsd[0]),&(actsolv->sol[0]),1.0);
/*----------------------------- put actual total displacements to nodes */
solserv_result_total(
                     actfield,
                     actintra,
                     &(actsolv->sol[0]),
                     0,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray])
                    );
/*------------------------------------------- calculate residual energy */
solserv_dot_vec(actintra,&(rsd[0]),&(re[0]),&energy);
nln_data->renergy = energy;
/*------------------------------ calculate norm of out-of-balance-loads */
solserv_vecnorm_euclid(actintra,&(re[0]),&renorm);
nln_data->renorm = renorm;
/*--------------- calculate the up-to-date norm of out-of-balance-loads */
solserv_copy_vec(&(actsolv->rhs[actsysarray+1]),&(re[1]));
solserv_scalarprod_vec(&(re[1]),rlnew);  /* up-to-date int - ext forces */
assemble_vec(actintra, &(actsolv->sysarray_typ[actsysarray]), &(actsolv->sysarray[actsysarray]), &(re[1]), intforce, -1.0 );
solserv_vecnorm_euclid(actintra,&(re[1]),&resnorm_uptodate);
/*-------------------------------------------------------- update norms */
nln_data->rrnorm = nln_data->rinorm * rlnew;
/*------------------------------------- calculate norm of displacements */
solserv_vecnorm_euclid(actintra,&(dispi[0]),&dnorm);
/*----------------------- calculate norm of updated incr. displacements */
solserv_vecnorm_euclid(actintra,&(rsd[0]),&dinorm);
/*----------------------------------------------- check for convergence */
told=0.0;
if (dinorm > statvar->toldisp) told = dinorm/dnorm;
if (told >statvar->toldisp || resnorm_uptodate > statvar->tolresid) convergence=0;
else                        convergence=1;

/*------------------------------------------------------- make printout */
solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sol[0]),
                   cdof,
                   &disval);
if (actintra->intra_rank==0)
conequ_printiter(*itnum,disval,rlnew,dinorm,resnorm_uptodate,energy,dnorm,nln_data->rrnorm);
if (convergence)
{
  if (ioflags.relative_displ>0)
  {
    for (j=0; j<ioflags.relative_displ; j++)
    {
       rel_dof = reldof[j];
       solserv_getele_vec(actintra,
                          &(actsolv->sysarray_typ[actsysarray]),
                          &(actsolv->sysarray[actsysarray]),
                          &(actsolv->sol[0]),
                          rel_dof,
                          &reldis);
       reldisval[j] = reldis;
    }
    fprintf(allfiles.out_err,"final state: disval %20.10f reldisA %20.10f reldisB %20.10f reldisC %20.10f reldisD %20.10f rlnew %20.10f\n",disval,reldisval[0],reldisval[1],reldisval[2],reldisval[3],rlnew);
  }
  else
  {
    fprintf(allfiles.out_err,"final state: disval %20.10f rlnew %20.10f\n",disval,rlnew);
  }
}
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
/*------------------------------- output for load - displacement curve */
if (par.myrank==0)
{
   fprintf(allfiles.out_cur,"%f %f\n",disval, rlnew);
   fflush(allfiles.out_cur);
}
/*----------------------------------------------------------------------*/
amdel(&intforce_a);
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of conequ */



/*----------------------------------------------------------------------*
 |  make load factor for arclenght control                    m.gee 1/02|
 |                                                                      |
 *----------------------------------------------------------------------*/
/*
C-----------------------------------------------------------------------
C        TOPIC: NEWTON-RAPHSON ITERATION (CRISFIELD'S METHOD)
C-----------------------------------------------------------------------
C        DISP   .. ACTUAL DISPLACEMENTS
C        DISP0  .. DISPLACEMENT OF LAST CONVERGED SOLUTION
C        DISPI  .. DISPLACEMENTS OF INCREMENT
C        RSD    .. DISPLACEMENTS OF ITERATION
C        RSD1   .. SOLUTION WITH LOAD AS RHS-VECTOR
C        RSD2   .. SOLUTION WITH OUT-OF-BALANCE-LOAD P AS RHS-VECTOR
C        RLNEW  .. ACTUAL LOAD PARAMETER
C        RL0    .. LOAD PARAMETER OF LAST CONVERGED SOLUTION
C        RLI    .. LOAD PARAMETER INCREMENT
C        DS0    .. ARC LENGTH
C        NEQ    .. NUMBER OF EQUATIONS
C-----------------------------------------------------------------------
*/
void increment_controlarc(INTRA         *actintra,
                          SOLVAR        *actsolv,
                          INT            actsysarray,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,
                          DIST_VECTOR   *disp,
                          DOUBLE         rlnew,
                          DOUBLE         rlold,
                          DOUBLE         stepsize,
                          DOUBLE        *rli)
{
INT                  i;
DOUBLE               rsd1,rsd2;
DOUBLE               val1,val2,val3,val4,val5,val6;
DOUBLE               a1,a2,a3,a4,a5;
DOUBLE               ddisp;
DOUBLE               actdisp;
DOUBLE               convdisp;
DOUBLE               drl;
DOUBLE               valsqr;
DOUBLE               rli1,rli2;
DOUBLE               valcos1,valcos2;

#ifdef DEBUG
dstrc_enter("increment_controlarc");
#endif
/*----------------------------------------------------------------------*/
val1 = val2 = val3 = val4 = val5 = val6 = 0.0;
   for (i=0; i<dispi[0].numeq_total; i++)
   {
      solserv_getele_vec(actintra,
                         &(actsolv->sysarray_typ[actsysarray]),
                         &(actsolv->sysarray[actsysarray]),
                         disp,
                         i,
                         &actdisp);
      solserv_getele_vec(actintra,
                         &(actsolv->sysarray_typ[actsysarray]),
                         &(actsolv->sysarray[actsysarray]),
                         &(dispi[1]),
                         i,
                         &convdisp);
      solserv_getele_vec(actintra,
                         &(actsolv->sysarray_typ[actsysarray]),
                         &(actsolv->sysarray[actsysarray]),
                         &(rsd[1]),
                         i,
                         &rsd1);
      solserv_getele_vec(actintra,
                         &(actsolv->sysarray_typ[actsysarray]),
                         &(actsolv->sysarray[actsysarray]),
                         &(rsd[2]),
                         i,
                         &rsd2);

      ddisp  = actdisp - convdisp;
      val1  += rsd1 * rsd1;
      val2  += rsd1 * (ddisp+rsd2);
      val3  += (ddisp+rsd2) * (ddisp+rsd2);
      val4  += ddisp * ddisp;
      val5  += ddisp * rsd2;
      val6  += ddisp * rsd1;
   }/* end of (i=0; i<dispi[0].numeq_total; i++) */
/*
C-----------------------------------------------------------------------
C            PARAMETER FOR QUADRATIC FUNCTION
C               A1*RLI*RLI + A2*RLI +A3 = 0
C               A3,A4  TO CHOOSE ONE OF BOTH SOLUTIONS FOR RLI
C            ARCSCL --> SCALING FACTOR
C               = 1  SPERICAL ARC-LENGTH-METHOD
C               = 0  CYLINDRICAL ARC-LENGTH-METHOD
C-----------------------------------------------------------------------
*/
   drl = rlnew - rlold;
   a1 =     val1 + statvar->arcscl*statvar->arcscl;
   a2 = 2.0*val2 + 2.0*drl*statvar->arcscl*statvar->arcscl;
   a3 =     val3 + drl*drl*statvar->arcscl*statvar->arcscl - stepsize*stepsize;
   a4 =     val4 + val5;
   a5 =     val6;
   valsqr = a2*a2 - 4.0*a1*a3;
   if (valsqr < 0.0) dserror("No solution for qudratic function in arclenght-method, try new stepsize");
   else if (valsqr<EPS8) *rli = -a2 / (2.0*a1);
   else
   {
      rli1 = (-a2 + sqrt(valsqr)) / (2.0*a1);
      rli2 = (-a2 - sqrt(valsqr)) / (2.0*a1);
      valcos1 = a4+a5*rli1;
      valcos2 = a4+a5*rli2;
      if (valcos1>valcos2) *rli = rli1;
      else                 *rli = rli2;
   }
/*----------------------------------- make rsd[0] = rsd[1]*rli + rsd[2] */
/*---------------------- make dispi[0] = dispi[0] + rsd[1]*rli + rsd[2] */
solserv_copy_vec(&(rsd[1]),&(rsd[0]));
solserv_scalarprod_vec(&(rsd[0]),*rli);
solserv_add_vec(&(rsd[0]),&(dispi[0]),1.0);
solserv_add_vec(&(rsd[2]),&(rsd[0]),1.0);
solserv_add_vec(&(rsd[2]),&(dispi[0]),1.0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of increment_controlarc */




/*----------------------------------------------------------------------*
 |  make load factor for displacement control                 m.gee 1/02|
 |                                                                      |
 *----------------------------------------------------------------------*/
void increment_controldisp(INTRA *actintra,
                          SOLVAR        *actsolv,
                          INT            actsysarray,
                          INT            cdof,
                          DIST_VECTOR   *rsd,
                          DIST_VECTOR   *dispi,
                          DOUBLE        *rli)
{
DOUBLE rsd1,rsd2;

#ifdef DEBUG
dstrc_enter("increment_controldisp");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------ get values of control node */
solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(rsd[2]),
                   cdof,
                   &rsd2);
solserv_getele_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(rsd[1]),
                   cdof,
                   &rsd1);
*rli = -rsd2/rsd1;
/*----------------------------------- make rsd[0] = rsd[1]*rli + rsd[2] */
/*---------------------- make dispi[0] = dispi[0] + rsd[1]*rli + rsd[2] */
solserv_copy_vec(&(rsd[1]),&(rsd[0]));
solserv_scalarprod_vec(&(rsd[0]),*rli);
solserv_add_vec(&(rsd[0]),&(dispi[0]),1.0);
solserv_add_vec(&(rsd[2]),&(rsd[0]),1.0);
solserv_add_vec(&(rsd[2]),&(dispi[0]),1.0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of increment_controldisp */




/*----------------------------------------------------------------------*
 |  print out iteration head                                 m.gee 11/01|
 |  kstep                           load or time step we are in         |
 |  controltyp                      type of control algorithm           |
 |  cdof                            number of dof that is controlled    |
 *----------------------------------------------------------------------*/
void conequ_printhead(INT kstep, NR_CONTROLTYP  controltyp, INT cdof, DOUBLE csp)
{
#ifdef DEBUG
dstrc_enter("conequ_printhead");
#endif
/*----------------------------------------------------------------------*/
   printf("----------------------------------------------------------------------------------------------\n");
   printf("Incemental load step No. %d\n",kstep);
   printf("----------------------------------------------------------------------------------------------\n");
   switch(controltyp)
   {
   case control_disp:
   printf("Displacement control\n");
   break;
   case control_arc:
   printf("Arclenght control\n");
   break;
   case control_load:
   printf("Load control\n");
   break;
   case control_none:
      dserror("Unknown typ of path following technique");
   break;
   default:
      dserror("Unknown typ of path following technique");
   break;
   }
   printf("----------------------------------------------------------------------------------------------\n");
   printf("DISVAL....Total Displacement at equation %d\n",cdof);
   printf("----------------------------------------------------------------------------------------------\n");
   printf("RLNEW.....Actual Load factor\n");
   printf("DINORM....Norm of Residual Displacements\n");
   printf("RENORM....Norm of Out-of-Balance Loads\n");
   printf("ENERGY....Product of Out-of-Balance Loads and Residual Displacements\n");
   printf("DNORM.....Norm of incemental Displacements\n");
   printf("RNORM.....Norm of total Load\n");
   printf("CSP=%f (Current Stiffness to Reference Stiffness of first step)\n",csp);
   printf("----------------------------------------------------------------------------------------------\n");
   printf("ITE  DISVAL       RLNEW        DINORM       RENORM       ENERGY       DNORM        RNORM\n");
/*----------------------------------------------------------------------*/
   fprintf(allfiles.out_err,"----------------------------------------------------------------------------------------------\n");
   fprintf(allfiles.out_err,"Incemental load step No. %d\n",kstep);
   fprintf(allfiles.out_err,"----------------------------------------------------------------------------------------------\n");
   switch(controltyp)
   {
   case control_disp:
   fprintf(allfiles.out_err,"Displacement control\n");
   break;
   case control_arc:
   fprintf(allfiles.out_err,"Arclenght control\n");
   break;
   case control_load:
   fprintf(allfiles.out_err,"Load control\n");
   break;
   case control_none:
      dserror("Unknown typ of path following technique");
   break;
   default:
      dserror("Unknown typ of path following technique");
   break;
   }
   fprintf(allfiles.out_err,"----------------------------------------------------------------------------------------------\n");
   fprintf(allfiles.out_err,"DISVAL....Total Displacement at equation %d\n",cdof);
   fprintf(allfiles.out_err,"----------------------------------------------------------------------------------------------\n");
   fprintf(allfiles.out_err,"RLNEW.....Actual Load factor\n");
   fprintf(allfiles.out_err,"DINORM....Norm of Residual Displacements\n");
   fprintf(allfiles.out_err,"RENORM....Norm of Out-of-Balance Loads\n");
   fprintf(allfiles.out_err,"ENERGY....Product of Out-of-Balance Loads and Residual Displacements\n");
   fprintf(allfiles.out_err,"DNORM.....Norm of incemental Displacements\n");
   fprintf(allfiles.out_err,"RNORM.....Norm of total Load\n");
   fprintf(allfiles.out_err,"CSP=%f (Current Stiffness to Reference Stiffness of first step)\n",csp);
   fprintf(allfiles.out_err,"----------------------------------------------------------------------------------------------\n");
   fprintf(allfiles.out_err,"ITE  DISVAL       RLNEW        DINORM       RENORM       ENERGY       DNORM        RNORM\n");
   fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of conequ_printhead */
/*----------------------------------------------------------------------*
 |  print out iteration info                                 m.gee 11/01|
 |                                                                      |
 | itnum                                   number of actual iteration   |
 | disval               displcament value of controled dof in this step |
 | rlnew                                       actual total load factor |
 | dinorm                    norm of residual incremental displacements |
 | renorm                                  Norm of out-of-balance-loads |
 | energy   norm of product of out-of-balance-loads and residual displ. |
 | dnorm                              norm of incremental displacements |
 | rrnorm                                            norm of total load |
 *----------------------------------------------------------------------*/
void conequ_printiter(INT itnum, DOUBLE disval, DOUBLE rlnew, DOUBLE dinorm,
                     DOUBLE renorm, DOUBLE energy, DOUBLE dnorm, DOUBLE rrnorm)
{
#ifdef DEBUG
dstrc_enter("conequ_printiter");
#endif
/*----------------------------------------------------------------------*/
    printf("%4d %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E\n",
   itnum,disval,rlnew,dinorm,renorm,energy,dnorm,rrnorm);
/*----------------------------------------------------- printout to err */
   fprintf(allfiles.out_err,"%4d %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E %-12.5E\n",
   itnum,disval,rlnew,dinorm,renorm,energy,dnorm,rrnorm);
   fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of conequ_printiter */
