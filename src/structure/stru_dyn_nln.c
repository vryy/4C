#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
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
 |  routine to control nonlinear dynamic structural analysis m.gee 11/01|
 *----------------------------------------------------------------------*/
void dyn_nln_structural() 
{
int             i;                  /* simply a counter */
int             numeq;              /* number of equations on this proc */
int             numeq_total;        /* total number of equations */
int             init;               /* flag for solver_control call */
int             itnum;              /* counter for NR-Iterations */
int             convergence;
double          dmax;

int             stiff_array;        /* number of the active system sparse matrix */
int             mass_array;         /* number of the active system sparse matrix */
int             damp_array;         /* number of the active system sparse matrix */
int             num_array;          /* number of global stiffness matrices */
int             actcurve;           /* number of active time curve */

SOLVAR         *actsolv;            /* pointer to active solution structure */
PARTITION      *actpart;            /* pointer to active partition */
FIELD          *actfield;           /* pointer to active field */
INTRA          *actintra;           /* pointer to active intra-communicator */
CALC_ACTION    *action;             /* pointer to the structures cal_action enum */
STRUCT_DYNAMIC *sdyn;

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
/*------------ the distributed system matrix, which is used for solving */
/* 
   NOTE: This routine uses 1 global sparse matrix, which was created
         in global_mask_matrices. The others are made by copying the
         sparsity pattern from the first one
         (Normally a control routine should not mix up different storage formats
         for system matrices)
*/         
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
/* intracommunicator (in case of nonlinear statics, this should be all) */
if (actintra->intra_fieldtyp != structure) goto end;
/*------------------------------------ check presence of damping matrix */
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

/* stiff_array already exists, so copy it to mass_array (and damp_array if needed) */
actsolv->sysarray_typ = 
(SPARSE_TYP*)REALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

actsolv->sysarray = 
(SPARSE_ARRAY*)REALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

/*-copy the matrices sparsity mask from stiff_array (and to damp_array) */
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


/*----------------------------------------allocate 4 dist. load vectors */
actsolv->nrhs = 4;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));

/*-------------------- there is one solution vector to hold total displ.*/
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
/*----------------------------------------- allocate 3 DIST_VECTOR fie */ 
intforce = amdef("intforce",&intforce_a,numeq_total,1,"DV");
solserv_create_vec(&fie,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(fie[i]));

/*--------------------------------------allocate three working vectors */
solserv_create_vec(&work,3,numeq_total,numeq,"DV");
for (i=0; i<3; i++) solserv_zero_vec(&(work[i]));


/*------------------------ initialize solver on the matrix stiff_array */
/*
NOTE: solver init phase has to be called with each matrix one wants to solve with
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(dispi[0]),
               &(actsolv->rhs[0]),
               init);

/*----------------- init the assembly for stiffness and for mass matrix */
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
   if (ABS(sdyn->k_damp) > EPS12)
   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   sdyn->k_damp);
   if (ABS(sdyn->m_damp) > EPS12)
   solserv_add_mat(actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &(actsolv->sysarray_typ[mass_array]),
                   &(actsolv->sysarray[mass_array]),
                   sdyn->k_damp);
}

/*-------------------------------------- create the original rhs vector */
*action = calc_struct_eleload;
calrhs(
          actfield,
          actsolv,
          actpart,
          actintra,
          stiff_array,
          &(actsolv->rhs[2]),
          &(actsolv->rhs[3]),
          0,
          action
      );

/*--------------------------------------------- add the two rhs vectors */
solserv_add_vec(&(actsolv->rhs[3]),&(actsolv->rhs[2]),1.0);
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[3]));

/*----------------------- init the time curve applied to the loads here */
/*-------------- this control routine at the moment always uses curve 0 */
/*-------------------------------------------------- init the timecurve */
actcurve = 0;
dyn_init_curve(actcurve,
               sdyn->nstep,
               sdyn->dt,
               sdyn->maxtime);

/*-------------------------------------- get factor at a certain time t */
dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));

/*-------------------------------------- multiply load vector by rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[2]),dynvar.rldfac);

/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[2]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actfield,actsolv,stiff_array,mass_array,&vel[0],
     &(actsolv->sol[0]),&(actsolv->rhs[2]),0,sdyn->nstep);

sdyn->step = -1;
sdyn->time = 0.0;

/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains(actfield);
}

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
   rhs[0]    interpolated load vector and is working array

   fie[2]    internal forces at step t
   fie[1]    internal forces at step t-dt
   fie[0]    interpolated internal forces and is working array

   dispi[0]  displacement increment from t-dt to t

   sol[0]    total displacements at time t-dt
   sol[1]    total displacements at time t

   work[2]   working vector for sums and matrix-vector products 
   work[1]   working vector for sums and matrix-vector products 
   work[0]   working vector for sums and matrix-vector products 
   work[0]   is used to hold residual displacements in corrector 
             iteration

*/
timeloop:
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);
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
if (dynvar.dinorm < (dynvar.dnorm*sdyn->toldisp) ||
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

/*!!!!!!!!!!!!!!!!!!!!!!! make stress calculation and write output here */

/*------------------ update displacements, velocities and accelerations */













/*----------------------------------------------------------------------*/
end:
/*--------------------------------------------------- cleaning up phase */
amdel(&intforce_a);
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */





/*----------------------------------------------------------------------*
 |  make effective stiffness matrix                          m.gee 02/02|
 *----------------------------------------------------------------------*/
/*
  keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
*/  
int kefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *work,
                  int              stiff_array,
                  int              mass_array,
                  int              damp_array)
{

double a0,a1,a4;

#ifdef DEBUG 
dstrc_enter("kefnln_struct");
#endif
/*----------------------------------------------------------------------*/
a0 = dynvar->constants[6];
a1 = dynvar->constants[0];
a4 = dynvar->constants[3];
/*---------------------------------------------- make stiff_array *= a0 */
solserv_scal_mat(&(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 a0);
/*--------------------------------- make stiff_array += mass_array * a1 */
solserv_add_mat(actintra,
                &(actsolv->sysarray_typ[stiff_array]),
                &(actsolv->sysarray[stiff_array]),
                &(actsolv->sysarray_typ[mass_array]),
                &(actsolv->sysarray[mass_array]),
                a1);
/*--------------------------------- make stiff_array += damp_array * a4 */
if (damp_array>-1)
solserv_add_mat(actintra,
                &(actsolv->sysarray_typ[stiff_array]),
                &(actsolv->sysarray[stiff_array]),
                &(actsolv->sysarray_typ[damp_array]),
                &(actsolv->sysarray[damp_array]),
                a4);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of kefnln_struct */




/*----------------------------------------------------------------------*
 |  make effective load vector                               m.gee 02/02|
C     *                     EFFECTIVE LOAD VECTOR                      *
C     *                                                                *
C     *     RHSEFF=RHS + MASS(-DISPI,VEL,ACC) + DAMP(-DISPI,VEL,ACC)   *
 *----------------------------------------------------------------------*/
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
int pefnln_struct(STRUCT_DYN_CALC *dynvar, 
                  STRUCT_DYNAMIC  *sdyn,
                  FIELD           *actfield,
                  SOLVAR          *actsolv,
                  INTRA           *actintra,
                  DIST_VECTOR     *dispi,
                  DIST_VECTOR     *vel,
                  DIST_VECTOR     *acc,
                  DIST_VECTOR     *work,
                  int              mass_array,
                  int              damp_array)
{

double a1,a2,a3,a4,a5,a6;

#ifdef DEBUG 
dstrc_enter("pefnln_struct");
#endif
/*----------------------------------------------------------------------*/
a1 = -dynvar->constants[0];
a2 =  dynvar->constants[1];
a3 =  dynvar->constants[2];
a4 = -dynvar->constants[3];
a5 =  dynvar->constants[4];
a6 =  dynvar->constants[5];
/*----------------------------------------------------------------------*/
/*                   make work[0] = dispi[0]*a1 + vel[0]*a2 + acc[0]*a3 */
/*----------------------------------------------------------------------*/
/* work[0] = dispi[0]*a1 */
solserv_copy_vec(&dispi[0],&work[0]);
solserv_scalarprod_vec(&work[0],a1);
/* work[0] += vel[0]*a2 */
solserv_add_vec(&vel[0],&work[0],a2);
/* work[0] += acc[0]*a3 */
solserv_add_vec(&acc[0],&work[0],a3);

/*----------------------------------------------------------------------*/
/* make work[1] = dispi[0]*a4 + vel[0]*a5 + acc[0]*a6 if damping present */
/*----------------------------------------------------------------------*/
if (sdyn->damp==1)
{
   /* work[1] = dispi[0]*a4 */
   solserv_copy_vec(&dispi[0],&work[1]);
   solserv_scalarprod_vec(&work[1],a4);
   /* work[1] += vel[0]*a5 */
   solserv_add_vec(&vel[0],&work[1],a5);
   /* work[1] += acc[0]*a6 */
   solserv_add_vec(&acc[0],&work[1],a6);
}
/*----------------------------------------------------------------------*/
/* make rhs[0] = rhs[0] + mass_array * work[0] + damp_array * work[1]   */
/*----------------------------------------------------------------------*/
/* work[2] = mass_array * work[0] */
solserv_sparsematvec(actintra,
                     &work[2],
                     &(actsolv->sysarray[mass_array]),
                     &(actsolv->sysarray_typ[mass_array]),
                     &work[0]);
/* rhs[0] += work[2] */
solserv_add_vec(&work[2],&(actsolv->rhs[0]),1.0);
/* if damping matrix is present */
if (sdyn->damp==1)
{
   /* work[2] = damp_array * work[1] */
   solserv_sparsematvec(actintra,
                        &work[2],
                        &(actsolv->sysarray[damp_array]),
                        &(actsolv->sysarray_typ[damp_array]),
                        &work[1]);
   /* rhs[0] += work[2] */
   solserv_add_vec(&work[2],&(actsolv->rhs[0]),1.0);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of pefnln_struct */

