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
 |  update the vectors                                       m.gee 03/02|
 | displacements sol_new -> sol_old                                     |
 |               rhs_new -> rhs_old                                     |
 | new accelerations at time t                                          |
 | new velocities    at time t                                          |
 |                                                                      |
 | in node->sol[place][0..numdf-1] is stored:                           |
 |                                                                      |
 | in the nodes the results are stored the following way:               |
 | place 0 holds total displacements of free dofs at time t             |
 | place 1 holds velocities at time t                                   |
 | place 2 holds accels at time t                                       |
 | place 3 holds prescribed displacements at time t-dt                  |
 | place 4 holds prescribed displacements at time t                     |
 | place 5 holds place 4 - place 3                                      |
 | place 6 holds the  velocities of prescribed dofs                     |
 | place 7 holds the  accels of prescribed dofs                         |
 | place 8 is working space                                             |
 |                                                                      |
 *----------------------------------------------------------------------*/
void dyn_nlnstructupd(FIELD *actfield,      STRUCT_DYN_CALC *dynvar, 
                      STRUCT_DYNAMIC *sdyn, SOLVAR *actsolv,
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
/*----------------------------------------------------------------------*/
solserv_copy_vec(work0,acc);
solserv_scalarprod_vec(acc,a1);
solserv_add_vec(vel,acc,a2);
solserv_add_vec(work1,acc,a3);

/*--------------------------------------- for prescribed displacements: */
/*  sol[7] = a1*(sol[4]-sol[3])+a2*sol[6]+a3*sol[7] */
/*----------------------------------------------------------------------*/
solserv_cpdirich(actfield,0,7,8);               /* make a copy of sol[7] in sol[8] */              
solserv_zerodirich(actfield,0,7);               /* init sol[7] to zero             */              
solserv_assdirich_fac(actfield,0,4,3,7,a1,-a1); /* sol[7]+=a1*sol[4]-a1*sol[3]     */
solserv_assdirich_fac(actfield,0,6,0,7,a2,0.0); /* sol[7]+=a2*sol[6]               */
solserv_assdirich_fac(actfield,0,8,0,7,a3,0.0); /* sol[7]+=a3*sol[8]               */

/*--------------- make new vel = a4*(sol_new-sol_old)+a5*vel+a6*acc_old */
/* this is equiv. to       vel = a4*work0            +a5*work2+a6*work1 */
/*----------------------------------------------------------------------*/
solserv_copy_vec(work0,vel);
solserv_scalarprod_vec(vel,a4);
solserv_add_vec(work2,vel,a5);
solserv_add_vec(work1,vel,a6);

/*--------------------------------------- for prescribed displacements: */
/* sol[6] = a4*(sol[4]-sol[3])+a5*sol[6]+a6*sol[7]_old */
/*----------------------------------------------------------------------*/
solserv_adddirich(actfield,0,6,0,6,a5,0.0);     /* sol[6] = sol[6] * a5                 */ 
solserv_assdirich_fac(actfield,0,4,3,6,a4,-a4); /* sol[6] += a4 * sol[4] + -a4 * sol[3] */
solserv_assdirich_fac(actfield,0,8,0,6,a6,0.0); /* sol[6] += a6 * sol[8] (=sol[7]_old)  */

/*--------------------------------------- zero the working place sol[8] */
solserv_zerodirich(actfield,0,8);  

/*--------------------------------- update the prescribed displacements */
solserv_adddirich(actfield,0,4,0,3,1.0,0.0);  
solserv_zerodirich(actfield,0,4);

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


/*---------------------------------------------------------------------------*
 |  dirichlet conditions to an element vector elevec_a       m.gee 3/02      |
 |  and then assembles this element vector of cond. dirich.conditions to the |
 |  global vector fullvec                                                    |
 |  this assembly is a special design for structural dynamics with           |
 |  with generalized alfa method                                             |
 |                                                                           |
 |  facs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))                            |
 |  facs[1] =  (1.0-alpham)*(1.0/beta)/dt                                    |
 |  facs[2] =  (1.0-alpham)/(2*beta) - 1                                     |
 |  facs[3] = -(1.0-alphaf)*(gamma/beta)/dt                                  |
 |  facs[4] =  (1.0-alphaf)*gamma/beta - 1                                   |
 |  facs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)                               |
 |  facs[6] = -(1.0-alphaf) or 0                                             |
 |  facs[7] =  raleigh damping factor for mass                               |
 |  facs[8] =  raleigh damping factor for stiffness                          |
 |  facs[9] =  dt                                                            |
 |                                                                           |
 | in the nodes the results are stored the following way:                    |
 | place 0 holds total displacements of free dofs at time t                  |
 | place 1 holds velocities at time t                                        |
 | place 2 holds accels at time t                                            |
 | place 3 holds prescribed displacements at time t-dt                       |
 | place 4 holds prescribed displacements at time t                          |
 | place 5 holds place 4 - place 3                                           |
 | place 6 holds the  velocities of prescribed dofs                          |
 | place 7 holds the  accels of prescribed dofs                              |
 | place 8 is working space                                                  |
 |                                                                           |
 |  fullvec contains the rhs mass, damping and stiffness parts due to        |
 |  prescribed displacements (dirichlet conditions !=0)                      |
 |                                                                           |
 | see PhD theses Mok page 165                                               |
 *---------------------------------------------------------------------------*/
void assemble_dirich_dyn(ELEMENT *actele, double *fullvec, int dim,
                         ARRAY *estif_global, ARRAY *emass_global, double *facs)
{
int                   i,j;
int                   counter;
int                   dof;
int                   numdf;
int                   iel;
int                   nd=0;
int                   idamp=0;
double                mdamp,kdamp;
double              **estif;
double              **emass;
double                dirich[MAXDOFPERELE];
double                dforces[MAXDOFPERELE];
int                   dirich_onoff[MAXDOFPERELE];
int                   lm[MAXDOFPERELE];
GNODE                *actgnode;
#ifdef DEBUG 
dstrc_enter("assemble_dirich_dyn");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------- check for presence of damping */
if (ABS(facs[7]) > EPS13 || ABS(facs[8]) > EPS13) 
{
   idamp=1;
   mdamp = facs[7];
   kdamp = facs[8];
}
/*----------------------------------------------------------------------*/
estif  = estif_global->a.da;
emass  = emass_global->a.da;
/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;
/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dforces[i] = 0.0;
   dirich_onoff[i] = 0;
}
/*-------------------------------- fill vectors dirich and dirich_onoff */
for (i=0; i<actele->numnp; i++)
{
   numdf    = actele->node[i]->numdf;
   actgnode = actele->node[i]->gnode;
   for (j=0; j<numdf; j++)
   {
      lm[i*numdf+j] = actele->node[i]->dof[j];
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
   }
}
/*----------------------------------------------------------------------*/
/*------------------------------ make the entries -K * (u(t) - u(t-dt)) */
/*--------------------------------                 K * facs[6]          */
/*                                       facs[6] is 0 in corrector call */
/*----------------------------------------------------------------------*/
/* 
prescribed displacement increment (u(t) - u(t-dt)) are in node->sol[5][0..numdf-1]
*/
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[5][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += estif[i][j] * dirich[j] * facs[6];
   }/* loop j over columns */
}/* loop i over rows */




/*----------------------------------------------------------------------*/
/*--------------------------- make the entries M  * h(-du,udot,udotdot) */
/* 
this consists of three parts:
   M * facs[0] * sol[5][0..numdf-1]       sol[5] holds (u(t) - u(t-dt))
   M * facs[1] * sol[6][0..numdf-1]       sol[6] holds ddot(t)
   M * facs[2] * sol[7][0..numdf-1]       sol[7] holds ddotdot(t)
*/   
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   M * facs[0] * sol[5][0..numdf-1]       sol[5] holds (u(t) - u(t-dt))
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[5][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += emass[i][j] * dirich[j] * facs[0];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   M * facs[1] * sol[6][0..numdf-1]       sol[6] holds ddot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[6][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += emass[i][j] * dirich[j] * facs[1];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   M * facs[2] * sol[7][0..numdf-1]       sol[7] holds ddotdot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[7][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += emass[i][j] * dirich[j] * facs[2];
   }/* loop j over columns */
}/* loop i over rows */



if (idamp)
{
/*----------------------------------------------------------------------*/
/*--------------------------- make the entries C  * e(-du,udot,udotdot) */
/* 
this consists of three parts:
   C * facs[3] * sol[5][0..numdf-1]       sol[5] holds (u(t) - u(t-dt))
   C * facs[4] * sol[6][0..numdf-1]       sol[6] holds ddot(t)
   C * facs[5] * sol[7][0..numdf-1]       sol[7] holds ddotdot(t)
*/   
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[3] * sol[5][0..numdf-1]       sol[5] holds (u(t) - u(t-dt))
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[5][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * facs[3];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[4] * sol[6][0..numdf-1]       sol[6] holds ddot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[6][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * facs[4];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
   C * facs[5] * sol[7][0..numdf-1]       sol[7] holds ddotdot(t)
*/
for (i=0; i<nd; i++) dirich[i]=0.0;
counter=0;
for (i=0; i<actele->numnp; i++)
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      if (dirich_onoff[counter]==0)
      {
         counter++;
         continue;
      }
      dirich[counter] = actele->node[i]->sol.a.da[7][j];
      counter++;
   }
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * facs[5];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
} /* end of if (idamp) */




/*-------- now assemble the vector dforces to the global vector fullvec */
for (i=0; i<nd; i++)
{
   if (lm[i] >= dim) continue;
   fullvec[lm[i]] += dforces[i];
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_dirich_dyn */
