/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../io/io.h"
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
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
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
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];


extern DOUBLE acttime;

/*---------------------------------------------------------------------------*
 |  routine to control nonlinear dyn explicit structural analysis m.gee 05/02|
 *---------------------------------------------------------------------------*/
void dyn_nln_stru_expl()
{
INT             i;                  /* simply a counter */
INT             numeq;              /* number of equations on this proc */
INT             numeq_total;        /* total number of equations */
INT             init;               /* flag for solver_control call */
INT             mod_disp,mod_stress;
INT             mod_res_write;
INT             restart;
DOUBLE          t0_res,t1_res;

DOUBLE          dt;
DOUBLE          dthalf;
INT             nstep;

DOUBLE          t0,t1;

INT             stiff_array;        /* indice of the active system sparse matrix */
INT             mass_array;         /* indice of the active system sparse matrix */
INT             damp_array;         /* indice of the active system sparse matrix */
INT             actcurve;           /* indice of active time curve */

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structure cal_action enum */
STRUCT_DYNAMIC *sdyn;               /* pointer to structural dynamic input data */

DIST_VECTOR    *vel;                /* total velocities */
DIST_VECTOR    *acc;                /* total accelerations */
DIST_VECTOR    *fie;                /* internal forces and working array */
DIST_VECTOR    *dispi;              /* distributed vector to hold incremental displacments */
DIST_VECTOR    *work;               /* working vectors */

ARRAY           intforce_a;         /* redundant vector of full length for internal forces */
DOUBLE         *intforce;
ARRAY           dirich_a;           /* redundant vector of full length for dirichlet-part of rhs */
DOUBLE         *dirich;

STRUCT_DYN_CALC dynvar;             /* variables to perform dynamic structural simulation */

#ifdef BINIO
BIN_OUT_FIELD   out_context;
#endif

INT             disnum = 0;

CONTAINER       container;          /* contains variables defined in container.h */
container.isdyn = 1;                /* dynamic calculation */
#ifdef DEBUG
dstrc_enter("dyn_nln_stru_expl");
#endif
/*----------------------------------------------------------------------*/
restart = genprob.restart;
/*--------------------------------------------------- set some pointers */
actfield           = &(field[0]);
actsolv            = &(solv[0]);
actpart            = &(partition[0]);
action             = &(calc_action[0]);
sdyn               =   alldyn[0].sdyn;
container.fieldtyp = actfield->fieldtyp;
container.disnum   = disnum;
/*----------------------------------------------------------------------*/
#ifdef PARALLEL
actintra    = &(par.intra[0]);
/* if we are not parallel, we have to allocate an alibi intra-communicator structure */
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = structure;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*-------------------------------- init the variables in dynvar to zero */
dynvar.rldfac = 0.0;
dynvar.rnorm  = 0.0;
dynvar.epot   = 0.0;
dynvar.eout   = 0.0;
dynvar.etot   = 0.0;
dynvar.ekin   = 0.0;
dynvar.dinorm = 0.0;
dynvar.dnorm  = 0.0;
for (i=0; i<20; i++) dynvar.constants[i] = 0.0;
acttime=0.0;
/*------------------------------------ check presence of damping matrix */
/*                 and set indice of stiffness and mass sparse matrices */
   stiff_array = 0;
   mass_array  = 1;

if (sdyn->damp==1)
{
   damp_array  = 2;
   actsolv->nsysarray=3;
}
else
{
   damp_array  =-1;
   actsolv->nsysarray=2;
}

/*--------------- stiff_array already exists, so copy the mask of it to */
/*------------------------------- mass_array (and damp_array if needed) */
/* reallocate the vector of sparse matrices and the vector of there types */
/* formerly lenght 1, now lenght 2 or 3 dependent on presence of damp_array */
actsolv->sysarray_typ =
(SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

actsolv->sysarray =
(SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
if (!actsolv->sysarray) dserror("Allocation of memory failed");

/*-copy the matrices sparsity mask from stiff_array to mass_array (and to damp_array) */
solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[stiff_array]),
                            &(actsolv->sysarray[stiff_array]),
                            &(actsolv->sysarray_typ[mass_array]),
                            &(actsolv->sysarray[mass_array]));
if (damp_array>0)
{
solserv_alloc_cp_sparsemask(  actintra,
                            &(actsolv->sysarray_typ[stiff_array]),
                            &(actsolv->sysarray[stiff_array]),
                            &(actsolv->sysarray_typ[damp_array]),
                            &(actsolv->sysarray[damp_array]));
}

/*------------------------------- init the dist sparse matrices to zero */
for (i=0; i<actsolv->nsysarray; i++)
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[i]),
                 &(actsolv->sysarray_typ[i])
                );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(&(actsolv->sysarray[stiff_array]),
                   actsolv->sysarray_typ[stiff_array],
                   &numeq,
                   &numeq_total);


/*---------------------------------------allocate 4 dist. vectors 'rhs' */
/*  these hold original load vector, load vector at time t and t-dt and */
/*                                             interpolated load vector */
actsolv->nrhs = 4;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));

/*-------------------- there are 2 solution vector to hold total displ.*/
/*                                  one at time t and one at time t-dt */
actsolv->nsol= 2;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));

/*-------------- there is one vector to hold incremental displacements */
solserv_create_vec(&dispi,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(dispi[i]));

/*-------------------------------------------- allocate one vector vel */
solserv_create_vec(&vel,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(vel[i]));

/*-------------------------------------------- allocate one vector acc */
solserv_create_vec(&acc,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(acc[i]));

/*-------------- allocate one redundant vector intforce of full lenght */
/* this is used by the element routines to assemble the internal forces*/
intforce = amdef("intforce",&intforce_a,numeq_total,1,"DV");
/*----------- create a vector of full length for dirichlet part of rhs */
dirich = amdef("dirich",&dirich_a,numeq_total,1,"DV");
/*----------------------------------------- allocate 3 DIST_VECTOR fie */
/*                    to hold internal forces at t, t-dt and inbetween */
solserv_create_vec(&fie,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(fie[i]));

/*--------------------------------------allocate three working vectors */
/*   By optimizing this routine one could live with one or two working */
/*    vectors, I needed three to make things straight-forward and easy */
solserv_create_vec(&work,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(work[i]));

/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init phase
      is cheap in computation (can be costly in memory)
      There will be no solver call on mass or damping array.
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[0]),
               init);
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[mass_array]),
               &(actsolv->sysarray[mass_array]),
               &work[0],
               &work[1],
               init);
if (damp_array>0)
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[damp_array]),
               &(actsolv->sysarray[damp_array]),
               &work[0],
               &work[1],
               init);
/*----------------- init the assembly for stiffness and for mass matrix */
/*                                           (damping is not assembled) */
init_assembly(actpart,actsolv,actintra,actfield,stiff_array,0);
init_assembly(actpart,actsolv,actintra,actfield,mass_array,0);

/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action,&container);

/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0 && ioflags.output_gid==1)
   out_gid_msh();

/*----------------------- call elements to calculate stiffness and mass */
*action = calc_struct_nlnstiffmass;
container.dvec         = NULL;
container.dirich       = NULL;
container.global_numeq = 0;
container.dirichfacs   = NULL;
container.kstep        = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);


/*------------------------------- make eigenvalue analysis if requested */
if (sdyn->eigen==1)
{
   if (actintra->intra_nprocs > 1) dserror("Eigenanalysis only on 1 processor");
   dyn_eigen(actfield,actpart,actsolv,actintra,stiff_array,mass_array);
}
/*-------------------------------------------- calculate damping matrix */
if (damp_array>0)
{
   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   sdyn->k_damp);

   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[mass_array]),
                   &(actsolv->sysarray[mass_array]),
                   sdyn->m_damp);
}


/*-------------------------------------- create the original rhs vector */
/*-------------------------- the approbiate action is set inside calrhs */
/*---------------------- this vector holds loads due to external forces */
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[1]),action,&container);
/*------------------------------------------------- copy the rhs vector */
solserv_copy_vec(&(actsolv->rhs[1]),&(actsolv->rhs[3]));
/*----------------------- init the time curve applied to the loads here */
/*-------------- this control routine at the moment always uses curve 0 */
/*-------------------------------------------------- init the timecurve */
actcurve = 0;
dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);

/*---------------------------------- get factor at a certain time t=0.0 */
dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));

/*-------------------------------------- multiply load vector by rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);

/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[1]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);

sdyn->step = -1;
sdyn->time = 0.0;

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   &(actsolv->sysarray_typ[stiff_array]), &(actsolv->sysarray[stiff_array]),
                   actfield, actpart, actintra, 0);
#endif

/*----------------------------------------- output to GID postprozessor */
if (par.myrank==0 && ioflags.output_gid==1)
{
   out_gid_domains(actfield, disnum);
}
/*------------------------------------------------------- printout head */
if (par.myrank==0) dyn_nlnstru_outhead_expl();
/*-------------------------------------------------- set some constants */
dyn_setconstants_expl(&dynvar,sdyn,sdyn->dt);
/*--------------------------------------- form effective left hand side */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
dyn_keff_expl(actintra,actsolv->sysarray_typ,actsolv->sysarray,
              stiff_array, mass_array, damp_array,
              &dynvar, sdyn);

/*-------------------------------- make triangulation of left hand side */
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[2]),
               init);
/*---------------------- set incremental displacements dispi[0] to zero */
solserv_zero_vec(&dispi[0]);
/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,disnum,actintra, &dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,disnum,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                     START LOOP OVER ALL STEPS                        */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
timeloop:
t0 = ds_cputime();
/*------------------------------------------------- write memory report */
if (par.myrank==0) dsmemreport();
/*--------------------------------------------------- check for restart */
if (restart)
{
   t0_res = ds_cputime();
   /*-------------- save the stepsize as it will be overwritten in sdyn */
   dt    = sdyn->dt;
   /*------ save the number of steps, as it will be overwritten in sdyn */
   nstep = sdyn->nstep;
   /*------------- save the restart interval, as it will be overwritten */
   mod_res_write = sdyn->res_write_evry;
   /*----------------------------------- the step to read in is restart */
#if defined(BINIO)
   restart_read_bin_nlnstructdyn(sdyn, &dynvar,
                                 &(actsolv->sysarray_typ[stiff_array]),
                                 &(actsolv->sysarray[stiff_array]),
                                 actfield, actpart, 0, actintra,
                                 actsolv->nrhs, actsolv->rhs,
                                 actsolv->nsol, actsolv->sol,
                                 1            , dispi       ,
                                 1            , vel         ,
                                 1            , acc         ,
                                 3            , fie         ,
                                 3            , work        ,
                                 restart);
#else
   restart_read_nlnstructdyn(restart,
                             sdyn,
                             &dynvar,
                             actfield,
                             actpart,
                             actintra,
                             action,
                             actsolv->nrhs, actsolv->rhs,
                             actsolv->nsol, actsolv->sol,
                             1            , dispi       ,
                             1            , vel         ,
                             1            , acc         ,
                             3            , fie         ,
                             3            , work        ,
                             &intforce_a,
                             &dirich_a,
                             &container /* contains variables defined in container.h */
                             );
#endif
   /*-------------------------------------- put the dt to the structure */
   sdyn->dt = dt;
   /*--------------------------------------- put nstep to the structure */
   sdyn->nstep = nstep;
   /*-------------------------------- put restart interval to structure */
   sdyn->res_write_evry = mod_res_write;
   /*------------------------------------------- switch the restart off */
   restart=0;
   /*----------------------------------------------------- measure time */
   t1_res = ds_cputime();
   fprintf(allfiles.out_err,"TIME for restart reading of step %d is %f sec\n",sdyn->step,t1_res-t0_res);
}
/*--------------------------------------------- increment step and time */
sdyn->step++;
/*------------------------------------------------ set new absolue time */
sdyn->time += sdyn->dt;
/*--- put time to global variable for time-dependent load distributions */
acttime = sdyn->time;
/*------------------------------------- do starting procedure in step 0 */
if (sdyn->step==0)
{
   /*------------------- set incremental displacements dispi[0] to zero */
   solserv_zero_vec(&dispi[0]);
   /*---------------------- set residual displacements in nodes to zero */
   solserv_result_resid(actfield,disnum,actintra,&dispi[0],0,
                        &(actsolv->sysarray[stiff_array]),
                        &(actsolv->sysarray_typ[stiff_array]));
   /*---------------------------- calculate internal forces at time 0.0 */
   amzero(&intforce_a);
   *action = calc_struct_internalforce;
   container.dvec         = intforce;
   container.dirich       = NULL;
   container.global_numeq = numeq_total;
   container.dirichfacs   = NULL;
   container.kstep        = 0;
   calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,
              &container,action);
   /*------------------------- store positive internal forces on fie[1] */
   solserv_zero_vec(&fie[1]);
   assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(fie[1]),intforce,1.0);
}
else
{
   solserv_copy_vec(&(fie[2]),&(fie[1]));
}
/*--------------- make sol[1] = sol[1] + dt * vel[0] + dt*dt/2 * acc[0] */
/*          dt = dynvar.constants[4]        dt*dt/2=dynvar.constants[2] */
solserv_add_vec(&(vel[0]),&(actsolv->sol[1]),dynvar.constants[4]);
solserv_add_vec(&(acc[0]),&(actsolv->sol[1]),dynvar.constants[2]);
/*---------------------------------- put the displacements to the nodes */
solserv_result_total(actfield,disnum,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------- copy rhs in rhs[1] at time t to rhs[2] at time t-dt */
solserv_copy_vec(&(actsolv->rhs[1]),&(actsolv->rhs[2]));
/*--------------------------------------- make load at time t in rhs[1] */
solserv_zero_vec(&(actsolv->rhs[1]));
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[1]),action,&container);
dyn_facfromcurve(actcurve,sdyn->time,&(dynvar.rldfac));
solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);
/*-------------------------------------- make internal forces at time t */
amzero(&intforce_a);
*action = calc_struct_internalforce;
container.dvec         = intforce;
container.dirich       = NULL;
container.global_numeq = numeq_total;
container.dirichfacs   = NULL;
container.kstep        = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,
           &container,action);
/*---------------------------- store positive internal forces on fie[2] */
solserv_zero_vec(&fie[2]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),&(fie[2]),intforce,1.0);
/*--------------------------------------- make rhs[0] = rhs[1] - fie[2] */
solserv_copy_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]));
solserv_add_vec(&(fie[2]),&(actsolv->rhs[0]),-1.0);
/*------------------------------------------ make effective load vector */
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
              mass_array,damp_array);
/*---------------------------------------------------- solve for system */
solserv_zero_vec(&work[0]);
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(work[0]),
               &(actsolv->rhs[0]),
               init);
/*-------------------------- update accelerations and velocities etc... */
/* sol[1] -> sol[0] */
/* t      -> t-dt   */
solserv_copy_vec(&(actsolv->sol[1]),&(actsolv->sol[0]));
/* vel[0] += dt/2 * (acc[0] + work[0]) */
dthalf = sdyn->dt / 2.0;
solserv_add_vec(&(acc[0]),&(vel[0]),dthalf);
solserv_add_vec(&(work[0]),&(vel[0]),dthalf);
/* work[0] -> acc[0] */
solserv_copy_vec(&(work[0]),&(acc[0]));
/* put sol[1] to nodes */
solserv_result_total(actfield,disnum,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*------------------------------------------ make all types of energies */
dynnle(&dynvar,sdyn,actintra,actsolv,&dispi[0],&fie[1],&fie[2],
       &(actsolv->rhs[1]),&(actsolv->rhs[2]),&work[0]);
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);
dynvar.etot = dynvar.epot + dynvar.ekin;
/*-------------------------------------- return velocities to the nodes */
solserv_result_total(actfield,disnum,actintra, &vel[0],1,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*------------------------------------------ return accel. to the nodes */
solserv_result_total(actfield,disnum,actintra, &acc[0],2,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*------------------------------- check whether to write results or not */
mod_disp      = sdyn->step % sdyn->updevry_disp;
mod_stress    = sdyn->step % sdyn->updevry_stress;
/*------------------------------- check whether to write restart or not */
mod_res_write = sdyn->step % sdyn->res_write_evry;
/*------------------------------------------ perform stress calculation */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress==1)
{
   *action = calc_struct_stress;
   container.dvec         = NULL;
   container.dirich       = NULL;
   container.global_numeq = 0;
   container.dirichfacs   = NULL;
   container.kstep        = 0;
   calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   container.kstep = 0;
   calreduce(actfield,actpart,disnum,actintra,action,&container);
}
/*-------------------------------------------- print out results to out */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress==1 && ioflags.struct_disp==1 && ioflags.output_out==1)
{
  out_sol(actfield,actpart,disnum,actintra,sdyn->step,0);
}
/*--------------------------------------------- printout results to gid */
#ifdef BINIO
if (ioflags.output_bin==1)
{
if (mod_disp==0)
  if (ioflags.struct_disp==1) {
    out_results(&out_context, sdyn->time, sdyn->step, 0, OUTPUT_DISPLACEMENT);
    out_results(&out_context, sdyn->time, sdyn->step, 1, OUTPUT_VELOCITY);
    out_results(&out_context, sdyn->time, sdyn->step, 2, OUTPUT_ACCELERATION);
  }
if (mod_stress==0)
  if (ioflags.struct_stress==1)
    out_results(&out_context, sdyn->time, sdyn->step, 0, OUTPUT_STRESS);
}
#endif

if (par.myrank==0 && ioflags.output_gid==1)
{
   if (mod_disp==0)
   if (ioflags.struct_disp==1)
   {
      out_gid_sol("displacement",actfield,disnum,actintra,sdyn->step,0,ZERO);
      out_gid_sol("velocities",actfield,disnum,actintra,sdyn->step,1,ZERO);
      out_gid_sol("accelerations",actfield,disnum,actintra,sdyn->step,2,ZERO);
   }
   if (mod_stress==0)
   if (ioflags.struct_stress==1)
   out_gid_sol("stress"      ,actfield,disnum,actintra,sdyn->step,0,ZERO);
}
/*-------------------------------------- write restart data to pss file */
if (mod_res_write==0) {
#ifdef BINIO
restart_write_bin_nlnstructdyn(&out_context,
                               sdyn, &dynvar,
                               actsolv->nrhs, actsolv->rhs,
                               actsolv->nsol, actsolv->sol,
                               1            , dispi       ,
                               1            , vel         ,
                               1            , acc         ,
                               3            , fie         ,
                               3            , work);
#else
restart_write_nlnstructdyn(sdyn,
                           &dynvar,
                           actfield,
                           actpart,
                           actintra,
                           action,
                           actsolv->nrhs, actsolv->rhs,
                           actsolv->nsol, actsolv->sol,
                           1            , dispi       ,
                           1            , vel         ,
                           1            , acc         ,
                           3            , fie         ,
                           3            , work        ,
                           &intforce_a,
                           &dirich_a,
                           &container   /* contains variables defined in container.h */
                           );
#endif
}
/*----------------------------------------------------- print time step */
if (par.myrank==0) dyn_nlnstruct_outstep(&dynvar,sdyn,0,sdyn->dt);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
fprintf(allfiles.out_err,"TIME for step %d is %f sec\n",sdyn->step,t1-t0);
/*-------------------------------------- check time and number of steps */
if (sdyn->step < sdyn->nstep-1 && sdyn->time <= sdyn->maxtime)
goto timeloop;
/*----------------------------------------------------------------------*/
end:
/*--------------------------------------------------- cleaning up phase */
amdel(&intforce_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&dispi,1);
solserv_del_vec(&vel,1);
solserv_del_vec(&acc,1);
solserv_del_vec(&fie,3);
solserv_del_vec(&work,3);

#ifdef BINIO
destroy_bin_out_field(&out_context);
#endif

/*----------------------------------------------------------------------*/
#ifndef PARALLEL
CCAFREE(actintra);
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dyn_nln_stru_expl */
