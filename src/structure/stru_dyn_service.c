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
extern struct _DYNAMIC *dyn;   
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
 |  make effective stiffness matrix                          m.gee 02/02|

  keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
 *----------------------------------------------------------------------*/
  
void kefnln_struct(STRUCT_DYN_CALC *dynvar, 
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
     *     RHSEFF=RHS + MASS(-DISPI,VEL,ACC) + DAMP(-DISPI,VEL,ACC)   *

  Peff = rhs[0] - fie[0] 
         + M*(-a1*dispi[0]+a2*vel[0]+a3*acc[0]) 
         + D*(-a4*dispi[0]+a5*vel[0]+a6*acc[0]) (if present)
    
    a1 =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
    a2 = ((1.0-alpham) * (1.0/beta)/(DSQR(dt)))*dt
    a3 =  (1.0-alpham) / (2.0*beta) - 1.0
    a4 =  (1.0-alphaf) * ((gamma/beta)/dt)
    a5 = ((1.0-alphaf) * ((gamma/beta)/dt))*dt - 1.0
    a6 =  (gamma/beta)/2.0 - 1.0) * dt * (1.0-alphaf)

 *----------------------------------------------------------------------*/
void pefnln_struct(STRUCT_DYN_CALC *dynvar, 
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






/*----------------------------------------------------------------------*
 |  compute system energies                                  m.gee 02/02|
     ******************************************************************
     *                                                                *
     *                       +---------------+                        *
     *                       I    dynnle     I                        *
     *                       +---------------+                        *
     *                                                                *
     *  INCREMENT POTENTIAL AND EXTERNAL ENERGY IN NONLINEAR DYNAMIC  *
     *                     NO SUPPORTED DOF'S!!!                      *
     *                                                                *
     *      **  DELTA EPOT =  1/2 * (Fi(T)+Fi(T+DT)) * DISPI **       *
     *      **  DELTA EOUT =  1/2 * (R(T) +R(T+DT) ) * DISPI **       *
     ******************************************************************
(ported from fortran)
 *----------------------------------------------------------------------*/
void dynnle(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INTRA *actintra, SOLVAR *actsolv, 
           DIST_VECTOR *dispi, /* converged incremental displacements */
           DIST_VECTOR *fie1,  /* internal forces at time t-dt */
           DIST_VECTOR *fie2,  /* internal forces at time t */
           DIST_VECTOR *rhs1,  /* load at time t                      */ 
           DIST_VECTOR *rhs2,  /* load at time t-dt                   */ 
           DIST_VECTOR *work0)
           
{
double deltae;
double deltaa;
#ifdef DEBUG 
dstrc_enter("dyne");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- work0 = fie1 + fie2 */
solserv_zero_vec(work0);
solserv_add_vec(fie1,work0,1.0);
solserv_add_vec(fie2,work0,1.0);
/*-------------------------------------- deltae = (fie1 + fie2) * dispi */
solserv_dot_vec(actintra,work0,dispi,&deltae);
deltae /= 2.0;
/*------------------------------------------------- work0 = rhs1 + rhs2 */
solserv_zero_vec(work0);
solserv_add_vec(rhs1,work0,1.0);
solserv_add_vec(rhs2,work0,1.0);
/*------------------------------------------ deltaa = (rhs1+rhs2)*dispi */
solserv_dot_vec(actintra,work0,dispi,&deltaa);
deltaa /= 2.0;
/*----------------------------- increment potential and external energy */
dynvar->epot += deltae;
dynvar->eout += deltaa;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyne */




/*----------------------------------------------------------------------*
 |  compute system energy                                    m.gee 02/02|
     ******************************************************************
     *                                                                *
     *                       +---------------+                        *
     *                       I    D Y N E    I                        *
     *                       +---------------+                        *
     *                                                                *
     *                         kinetic ENERGY                         *
     *     **  EKIN = 1/2*VEL^T*M*VEL                         **      *
     ******************************************************************
 *----------------------------------------------------------------------*/
void dyne(STRUCT_DYN_CALC *dynvar,
         INTRA           *actintra,
         SOLVAR          *actsolv,
         int              mass_array,
         DIST_VECTOR     *vel,
         DIST_VECTOR     *work)
{
#ifdef DEBUG 
dstrc_enter("dyne");
#endif
/*------------------------------------------------- make kinetic energy */
/*------------------------------------------- perform help = mass * vel */
solserv_sparsematvec(actintra,
                     work,
                     &(actsolv->sysarray[mass_array]),
                     &(actsolv->sysarray_typ[mass_array]),
                     vel);
/*------------------------------------------------ perform vel^T * help */
solserv_dot_vec(actintra,vel,work,&(dynvar->ekin));
dynvar->ekin /= 2.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyne */




/*----------------------------------------------------------------------*
 |  set time integration constants                           m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_setconstants(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, double dt)
{
double *constants;
double  beta;
double  gamma;
double  alpham;
double  alphaf;
#ifdef DEBUG 
dstrc_enter("dyn_setconstants");
#endif
/*----------------------------------------------------------------------*/
constants = dynvar->constants;
beta      = sdyn->beta;
gamma     = sdyn->gamma;
alpham    = sdyn->alpha_m;
alphaf    = sdyn->alpha_f;
/*----------------------------------------------------------------------*/
constants[6]  = 1.0-alphaf;
constants[7]  = 1.0-alpham;
constants[8]  = (1.0/beta)/(DSQR(dt));
constants[9]  = -constants[8] * dt;
constants[10] = 1.0 - 0.5/beta;
constants[11] = (gamma/beta)/dt;
constants[12] = 1.0 - gamma/beta;
constants[13] = (1.0-(gamma/beta)/2.0)*dt;
constants[0]  = constants[7]*constants[8];
constants[1]  = constants[0]*dt;
constants[2]  = (constants[7]/2.0)/beta - 1.0;
constants[3]  = constants[6]*constants[11];
constants[4]  = constants[3]*dt - 1.0;
constants[5]  = ((gamma/beta)/2.0 - 1.0)*dt*constants[6];
constants[14] = dt;
constants[15] = beta;
constants[16] = gamma;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_setconstants */




/*----------------------------------------------------------------------*
 |  update the vectors                                       m.gee 02/02|
 | displacements sol_new -> sol_old                                     |
 |               rhs_new -> rhs_old                                     |
 | new accelerations at time t                                          |
 | new velocities    at time t                                          |
 *----------------------------------------------------------------------*/
void dyn_nlnstructupd(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, SOLVAR *actsolv,
                     DIST_VECTOR *sol_old, DIST_VECTOR *sol_new,
                     DIST_VECTOR *rhs_new, DIST_VECTOR *rhs_old,
                     DIST_VECTOR *vel,     DIST_VECTOR *acc,
                     DIST_VECTOR *work0,   DIST_VECTOR *work1,
                     DIST_VECTOR *work2)
{

double a1,a2,a3,a4,a5,a6;

#ifdef DEBUG 
dstrc_enter("dyn_nlnstructupd");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------- get some constants */
a1 = dynvar->constants[8];
a2 = dynvar->constants[9];
a3 = dynvar->constants[10];
a4 = dynvar->constants[11];
a5 = dynvar->constants[12];
a6 = dynvar->constants[13];
/*----------------------------------------------------------------------*/
/*----------------- work0 = sol_new - sol_old == displacement increment */
solserv_copy_vec(sol_new,work0);
solserv_add_vec(sol_old,work0,-1.0);
/*-------------------------- copy displacements from sol_new to sol_old */
solserv_copy_vec(sol_new,sol_old);
/*---------------------------- copy load vector from rhs_new to rhs_old */
solserv_copy_vec(rhs_new,rhs_old);
/*--------------------------------- make temporary copy of acc to work1 */
solserv_copy_vec(acc,work1);
/*--------------------------------- make temporary copy of vel to work2 */
solserv_copy_vec(vel,work2);
/*--------------- make new acc = a1*(sol_new-sol_old)+a2*vel+a3*acc_old */
/* this is equiv. to       acc = a1*work0            +a2*vel+a3*work1   */
solserv_copy_vec(work0,acc);
solserv_scalarprod_vec(acc,a1);
solserv_add_vec(vel,acc,a2);
solserv_add_vec(work1,acc,a3);
/*--------------- make new vel = a4*(sol_new-sol_old)+a5*vel+a6*acc_old */
/* this is equiv. to       vel = a4*work0            +a5*work2+a6*work1 */
solserv_copy_vec(work0,vel);
solserv_scalarprod_vec(vel,a4);
solserv_add_vec(work2,vel,a5);
solserv_add_vec(work1,vel,a6);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nlnstructupd */



/*----------------------------------------------------------------------*
 |  print head of nln structural dynamics                    m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_nlnstruct_outhead(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn)
{
double  beta;
double  gamma;
double  alpham;
double  alphaf;
#ifdef DEBUG 
dstrc_enter("dyn_nlnstruct_outhead");
#endif
/*----------------------------------------------------------------------*/
beta      = sdyn->beta;
gamma     = sdyn->gamma;
alpham    = sdyn->alpha_m;
alphaf    = sdyn->alpha_f;
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
fprintf(allfiles.out_err,"       Nonlinear Structural Dynamics\n");
fprintf(allfiles.out_err,"\n");
fprintf(allfiles.out_err,"       Generalized Alpha Method\n");
fprintf(allfiles.out_err,"  beta=%f gamma=%f alpha_m=%f alpha_f=%f\n",beta,gamma,alpham,alphaf);
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
printf("--------------------------------------------------------------------------\n");
printf("       Nonlinear Structural Dynamics\n");
printf("\n");
printf("       Generalized Alpha Method\n");
printf("  beta=%f gamma=%f alpha_m=%f alpha_f=%f\n",beta,gamma,alpham,alphaf);    
printf("--------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nlnstruct_outhead */



/*----------------------------------------------------------------------*
 |  print head of nln structural dynamics                    m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_nlnstruct_outstep(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, int numiter)
{
#ifdef DEBUG 
dstrc_enter("dyn_nlnstruct_outstep");
#endif
/*----------------------------------------------------------------------*/
printf("STEP=%6d | NSTEP=%6d | TIME=%-14.8E | NUMITER=%3d | ETOT=%-14.8E | \n",
       sdyn->step,sdyn->nstep,sdyn->time,numiter+1,dynvar->etot);
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"STEP=%6d | NSTEP=%6d | TIME=%-14.8E | NUMITER=%3d | ETOT=%-14.8E | EPOT=%-14.8E | EKIN=%-14.8E | EOUT=%-14.8E\n",
                         sdyn->step,sdyn->nstep,sdyn->time,numiter+1,dynvar->etot,dynvar->epot,dynvar->ekin,dynvar->eout);
/*----------------------------------------------------------------------*/
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nlnstruct_outstep */
