#include "../headers/standardtypes.h"
#include "../headers/solution.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | global variable *partition, vector of lenght numfld of structures    |
 | PARTITION is defined in global_control.c                             |
 *----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | ranks and communicators                                              |
 | This structure struct _PAR par; is defined in main_ccarat.c
 *----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | struct _DYNAMIC      *dyn;                                           |
 *----------------------------------------------------------------------*/
#ifdef SUSE73
extern DYNAMIC *dyn;   
#else
extern struct _DYNAMIC *dyn;   
#endif
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _DYNAMIC      *curve;                                         |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];




/*----------------------------------------------------------------------*
 |  routine to control nonlinear dynamic structural analysis m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_nln_structural() 
{
int             i;                  /* simply a counter */
int             numeq;              /* number of equations on this proc */
int             numeq_total;        /* total number of equations */
int             init;               /* flag for solver_control call */
int             itnum;              /* counter for NR-Iterations */
int             convergence;        /* convergence flag */
double          dmax;               /* infinity norm of residual displacements */

int             stiff_array;        /* indice of the active system sparse matrix */
int             mass_array;         /* indice of the active system sparse matrix */
int             damp_array;         /* indice of the active system sparse matrix */
int             num_array;          /* indice of global stiffness matrices */
int             actcurve;           /* indice of active time curve */

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structure cal_action enum */
STRUCT_DYNAMIC *sdyn;               /* pointer to structural dynamic input data */

DIST_VECTOR    *vel;                /* total velocities */              
DIST_VECTOR    *acc;                /* total accelerations */
DIST_VECTOR    *fie;                /* internal forces and working array */
ARRAY           intforce_a;
double         *intforce;
DIST_VECTOR    *dispi;              /* distributed vector to hold incremental displacments */ 
DIST_VECTOR    *work;               /* working vectors */
 
STRUCT_DYN_CALC dynvar;             /* variables to perform dynamic structural simulation */              

#ifdef DEBUG 
dstrc_enter("dyn_nln_structural");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------- set some pointers */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
sdyn        =   dyn[0].sdyn;
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
/* if we are not parallel, we have to allocate an alibi intra-communicator structure */
#else
actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
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
(SPARSE_TYP*)REALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

actsolv->sysarray = 
(SPARSE_ARRAY*)REALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

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
solserv_getmatdims(actsolv->sysarray[stiff_array],
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
init_assembly(actpart,actsolv,actintra,actfield,stiff_array);
init_assembly(actpart,actsolv,actintra,actfield,mass_array);

/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action);

/*----------------------- call elements to calculate stiffness and mass */
*action = calc_struct_nlnstiffmass;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,NULL,0,0,action);

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
                   sdyn->k_damp);
}
/*-------------------------------------- create the original rhs vector */
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[2]),&(actsolv->rhs[3]),0,action);

/*--------------------------------------------- add the two rhs vectors */
solserv_add_vec(&(actsolv->rhs[3]),&(actsolv->rhs[2]),1.0);
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[3]));

/*----------------------- init the time curve applied to the loads here */
/*-------------- this control routine at the moment always uses curve 0 */
/*-------------------------------------------------- init the timecurve */
actcurve = 0;
dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);

/*---------------------------------- get factor at a certain time t=0.0 */
dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));

/*-------------------------------------- multiply load vector by rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[2]),dynvar.rldfac);

/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[2]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);

sdyn->step = -1;
sdyn->time = 0.0;

/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains(actfield);
}
/*------------------------------------------------------- printout head */
if (par.myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);
/*----------------------------------------------------------------------*/
/*                     START LOOP OVER ALL STEPS                        */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                     Predictor                                        */
/*----------------------------------------------------------------------*/
/*
   rhs[3]    original load vector
   rhs[2]             load vector at time t-dt
   rhs[1]             load vector at time t
   rhs[0]    interpolated load vector and working array

   fie[2]    internal forces at step t
   fie[1]    internal forces at step t-dt
   fie[0]    interpolated internal forces and working array

   dispi[0]  displacement increment from t-dt to t

   sol[0]    total displacements at time t-dt
   sol[1]    total displacements at time t
   
   vel[0]    velocities    at t-dt
   acc[0]    accelerations at t-dt

   work[2]   working vector for sums and matrix-vector products 
   work[1]   working vector for sums and matrix-vector products 
   work[0]   working vector for sums and matrix-vector products 
   work[0]   is used to hold residual displacements in corrector 
             iteration
             
   in the nodes, displacements are kept in node[].sol[0][0..numdf-1]
                 velocities    are kept in node[].sol[1][0..numdf-1]
                 accelerations are kept in node[].sol[2][0..numdf-1]

*/
timeloop:
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);

/*---------------------- set incremental displacements dispi[0] to zero */
solserv_zero_vec(&dispi[0]);

/*------------------------- set residual displacements in nodes to zero */
solserv_result_resid(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*--------------------------------------------- increment step and time */
sdyn->step++;
sdyn->time += sdyn->dt;

/*----------------------------------------------------------------------*/
/*                     PREDICTOR                                        */
/*----------------------------------------------------------------------*/
/*---------------------- copy initial load vector from rhs[3] to rhs[1] */
solserv_copy_vec(&(actsolv->rhs[3]),&(actsolv->rhs[1]));

/*------------------------------------------------ get factor at time t */
dyn_facfromcurve(actcurve,sdyn->time,&(dynvar.rldfac));

/*------------------------ multiply rhs[1] by actual load factor rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);

/*----- calculate tangential stiffness and internal forces at time t-dt */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
amzero(&intforce_a);
*action = calc_struct_nlnstiff;
calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,intforce,numeq_total,0,action);

/*---------------------------- store positive internal forces on fie[1] */
solserv_zero_vec(&fie[1]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(fie[1]),intforce,1.0);

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/*---------- subtract internal forces from interpolated external forces */
solserv_add_vec(&(fie[1]),&(actsolv->rhs[0]),-1.0);

/*--------------------- create effective load vector (rhs[0]-fie[2])eff */
/*
  Peff = rhs[0] - fie[0] 
         + M*(-a1*dispi[0]+a2*vel[0]+a3*acc[0]) 
         + D*(-a4*dispi[0]+a5*vel[0]+a6*acc[0]) (if present)
    
    a1 =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
    a2 = ((1.0-alpham) * (1.0/beta)/(DSQR(dt)))*dt
    a3 =  (1.0-alpham) / (2.0*beta) - 1.0
    a4 =  (1.0-alphaf) * ((gamma/beta)/dt)
    a5 = ((1.0-alphaf) * ((gamma/beta)/dt))*dt - 1.0
    a6 =  (gamma/beta)/2.0 - 1.0) * dt * (1.0-alphaf)
*/
/*----------------------------------------------------------------------*/
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
              mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
/*
  keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
*/  
/*----------------------------------------------------------------------*/
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
              damp_array);
/*------------- call for solution of system dispi[0] = Keff^-1 * rhs[0] */
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[0]),
               init);
/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);
/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------------------------------------------------------*/
/*                     PERFORM EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/
itnum=0;
iterloop:
/*------------ zero the stiffness matrix and vector for internal forces */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
amzero(&intforce_a);

/* call element routines for calculation of tangential stiffness and intforce */
*action = calc_struct_nlnstiff;
calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,intforce,numeq_total,0,action);

/*---------------------------- store positive internal forces on fie[2] */
solserv_zero_vec(&fie[2]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(fie[2]),intforce,1.0);

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/* interpolate internal forces fie[0] = (1-alfaf)fie[2] + alphaf*fie[1] */
solserv_copy_vec(&fie[2],&fie[0]);
solserv_scalarprod_vec(&fie[0],(1.0-sdyn->alpha_f));
solserv_add_vec(&fie[1],&fie[0],sdyn->alpha_f);

/*-- subtract interpolated internal forces from interp. external forces */
solserv_add_vec(&fie[0],&(actsolv->rhs[0]),-1.0);

/*--------------------- create effective load vector (rhs[0]-fie[0])eff */
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
              mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
              damp_array);

/*---------------------------------------- solve keff * rsd[0] = rhs[0] */
/* solve for residual displacements to correct incremental displacements*/
init=0;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(work[0]),
               &(actsolv->rhs[0]),
               init);

/*-------------------------- return residual displacements to the nodes */
solserv_result_resid(actfield,actintra,&work[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*-- update the incremental displacements by the residual displacements */
solserv_add_vec(&work[0],&dispi[0],1.0);

/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);
/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*----------------------------------------------- check for convergence */
convergence = 0;
dmax        = 0.0;
solserv_vecnorm_euclid(actintra,&(work[0]),&(dynvar.dinorm));
solserv_vecnorm_euclid(actintra,&(dispi[0]),&(dynvar.dnorm));
solserv_vecnorm_Linf(actintra,&(work[0]),&dmax);
if (dynvar.dinorm < sdyn->toldisp ||
    dynvar.dnorm  < EPS14 ||
    dynvar.dinorm < EPS14 && dmax < EPS12)
{
   convergence = 1;
}    
else
{
   itnum++;
   if (itnum==sdyn->maxiter) dserror("No convergence in maxiter steps");
   goto iterloop;
}
/*----------------------------------------------------------------------*/
/*                      END OF EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/

/*----------- make temporary copy of actsolv->rhs[2] to actsolv->rhs[0] */
/*                                   (load at t-dt)                     */
/* because in  dyn_nlnstructupd actsolv->rhs[2] is overwritten but is   */
/* still needed to compute energies in dynnle                           */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
/*------------------ update displacements, velocities and accelerations */
dyn_nlnstructupd(&dynvar,sdyn,actsolv,
                 &(actsolv->sol[0]),   /* total displacements at time t-dt */
                 &(actsolv->sol[1]),   /* total displacements at time t    */
                 &(actsolv->rhs[1]),   /* load vector         at time t    */
                 &(actsolv->rhs[2]),   /* load vector         at time t-dt */
                 &vel[0],              /* velocities          at time t    */
                 &acc[0],              /* accelerations       at time t    */
                 &work[0],             /* working arrays                   */
                 &work[1],             /* working arrays                   */
                 &work[2]);            /* working arrays                   */

/*-------------------------------------- return velocities to the nodes */
solserv_result_total(actfield,actintra, &vel[0],1,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*------------------------------------------ return accel. to the nodes */
solserv_result_total(actfield,actintra, &acc[0],2,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*------------------------------------------ make all types of energies */
dynnle(&dynvar,sdyn,actintra,actsolv,&dispi[0],&fie[1],&fie[2],
       &(actsolv->rhs[1]),&(actsolv->rhs[0]),&work[0]);
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);
dynvar.etot = dynvar.epot + dynvar.ekin;

/*------------------------------------------ perform stress calculation */
if (ioflags.struct_stress_file==1 || ioflags.struct_stress_gid==1)
{
   *action = calc_struct_stress;
   calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,NULL,0,0,action);
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   calreduce(actfield,actpart,actintra,action,0);
}
/*-------------------------------------------- print out results to out */
if (ioflags.struct_stress_file==1 && ioflags.struct_disp_file==1)
{
  out_sol(actfield,actpart,actintra,sdyn->step,0);
}
/*--------------------------------------------- printout results to gid */
if (par.myrank==0) 
{
   if (ioflags.struct_disp_gid==1)
   {
      out_gid_sol("displacement",actfield,actintra,sdyn->step,0);
      out_gid_sol("velocities",actfield,actintra,sdyn->step,1);
      out_gid_sol("accelerations",actfield,actintra,sdyn->step,2);
   }
   if (ioflags.struct_stress_gid==1)
   out_gid_sol("stress"      ,actfield,actintra,sdyn->step,0);
}

/*----------------------------------------------------- print time step */
if (par.myrank==0) dyn_nlnstruct_outstep(&dynvar,sdyn,itnum);

/*-------------------------------------- check time and number of steps */
if (sdyn->step < sdyn->nstep && sdyn->time <= sdyn->maxtime)
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
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */





