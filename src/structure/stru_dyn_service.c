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
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#ifdef GEMM
#include "../wall1/wall1.h"
#endif
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
 | struct _DYNAMIC      *dyn;                                           |
 *----------------------------------------------------------------------*/
extern struct _DYNAMIC *dyn;   
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
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/

#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif

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
                  INT              stiff_array,
                  INT              mass_array,
                  INT              damp_array)
{
DOUBLE a0,a1,a4;
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
                  INT              mass_array,
                  INT              damp_array)
{

DOUBLE a1,a2,a3,a4,a5,a6,a7;
#ifdef DEBUG 
dstrc_enter("pefnln_struct");
#endif
/*------------------------------------------------------Chung Hulbert --*/
if (sdyn->Typ == gen_alfa || sdyn->Typ == Gen_EMM)
{
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
}
else if (sdyn->Typ == centr_diff)/*----------- central difference method */
{
/*----------------------------------------------------------------------*/
a5 =  dynvar->constants[4];
a7 =  -1.0;
a5 =  -a5/2.0;
/*----------------------------------------------------------------------*/
if (sdyn->damp==1)
{
   /* copy vel[0] to work[0] */
   solserv_copy_vec(&vel[0],&work[0]);
   /* work[0] = vel[0] * a7 */
   solserv_scalarprod_vec(&work[0],a7);
   /* work[0] = vel[0] * a7 + acc[0] * a5 */
   solserv_add_vec(&acc[0],&work[0],a5);
   /* make rhs[0] += damp_array * work[0] */
   solserv_sparsematvec(actintra,
                        &work[1],
                        &(actsolv->sysarray[damp_array]),
                        &(actsolv->sysarray_typ[damp_array]),
                        &work[0]);
   solserv_add_vec(&work[1],&(actsolv->rhs[0]),1.0);
}
/*----------------------------------------------------------------------*/
} /* end of else if (sdyn->Typ = centr_diff) */
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
DOUBLE deltae;
DOUBLE deltaa;
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
         INT              mass_array,
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
 |                                                           m.gee 04/03|
 |  make incremental potential energy                                   |
 *----------------------------------------------------------------------*/
void dyn_epot(FIELD *actfield, INT disnum, INTRA *actintra, STRUCT_DYN_CALC *dynvar, DOUBLE *deltaepot)
{
INT               i,j;
INT               myrank;
ARRAY            *array;
NODE             *actnode;
DISCRET          *actdis;
DOUBLE            fint[MAXDOFPERNODE];
DOUBLE            send;
#ifdef DEBUG 
dstrc_enter("dyn_epot");
#endif
/*----------------------------------------------------------------------*/
actdis     = &(actfield->dis[disnum]);
myrank     = actintra->intra_rank;
*deltaepot = 0.0;
send       = 0.0;
/*----------------------------------------------------------------------*/
for (i=0; i<actdis->numnp; i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->proc != myrank) continue;
   array   = &(actnode->sol_increment);
   for (j=0; j<actnode->numdf; j++)
   fint[j] = 0.5*(array->a.da[1][j] + array->a.da[2][j]);
   for (j=0; j<actnode->numdf; j++)
      send += fint[j] * array->a.da[0][j];
}
#ifdef PARALLEL
MPI_Allreduce(&send,deltaepot,1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
*deltaepot = send;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_epot */

/*----------------------------------------------------------------------*
 |                                                           m.gee 04/03|
 |  make kinetic energy                                                 |
 *----------------------------------------------------------------------*/
void dyn_ekin(FIELD *actfield, SOLVAR *actsolv, PARTITION *actpart, INTRA *actintra, CALC_ACTION *action,
             CONTAINER *container, INT stiff_array, INT mass_array)
{
DOUBLE send;
#ifdef DEBUG 
dstrc_enter("dyn_ekin");
#endif
/*----------------------------------------------------------------------*/
*action = calc_struct_update_istep;
container->dvec          = NULL;
container->dirich        = NULL;
container->global_numeq  = 0;
container->dirichfacs    = NULL;
container->kstep         = 0;
container->ekin          = 0.0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,container,action);
#ifdef PARALLEL
send = container->ekin;
MPI_Allreduce(&send,&(container->ekin),1,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_ekin */


/*----------------------------------------------------------------------*
 |                                                           m.gee 04/03|
 |  make incremental kinetic energy at element level energy             |
 | all velocities are in node->sol.a.da[1][j] (including prescribed dofs)|
 | ekin = 0.5 * v^T * M * v                                             |
 *----------------------------------------------------------------------*/
void dyn_ekin_local(ELEMENT *actele,ARRAY *emass, CONTAINER  *container)
{
INT               i,j,nd;
DOUBLE            deltaekin = 0.0;
DOUBLE          **mass;
DOUBLE            sum;
DOUBLE            vel[MAXDOFPERELE];
DOUBLE            work[MAXDOFPERELE];
INT               counter=0;
#ifdef DEBUG 
dstrc_enter("dyn_ekin_local");
#endif
/*----------------------------------------------------------------------*/
mass = emass->a.da;
nd   = actele->numnp * actele->node[0]->numdf;
/*-------------------------------------------- get velocity to a vector */
for (i=0; i<actele->numnp; i++)
for (j=0; j<actele->node[i]->numdf; j++)
   vel[counter++] = actele->node[i]->sol.a.da[1][j];
/*---------------------------------------------------------- work = M*v */   
for (i=0; i<nd; i++)
{
   sum = 0.0;
   for (j=0; j<nd; j++) 
      sum += mass[i][j] * vel[j];
   work[i] = sum;
}
/*----------------------------------------- make deltaekin = vel * work */
for (i=0; i<nd; i++) deltaekin += vel[i]*work[i];
/*-------------------------------- add deltaekin to the container->ekin */
deltaekin *= 0.5;
container->ekin += deltaekin;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_ekin_local */
/*----------------------------------------------------------------------*
 |                                                           m.gee 04/03|
 |  make incremental external energy  eout = 0.5*(rhs(t)+rhs(t-dt))*dispi|
 *----------------------------------------------------------------------*/
void dyn_eout(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INTRA *actintra, SOLVAR *actsolv, 
              DIST_VECTOR *dispi, /* converged incremental displacements */
              DIST_VECTOR *rhs1,  /* load at time t                      */ 
              DIST_VECTOR *rhs2,  /* load at time t-dt                   */ 
              DIST_VECTOR *work0)
           
{
DOUBLE deltaa;
#ifdef DEBUG 
dstrc_enter("dyn_eout");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- work0 = rhs1 + rhs2 */
solserv_zero_vec(work0);
solserv_add_vec(rhs1,work0,1.0);
solserv_add_vec(rhs2,work0,1.0);
/*------------------------------------------ deltaa = (rhs1+rhs2)*dispi */
solserv_dot_vec(actintra,work0,dispi,&deltaa);
deltaa /= 2.0;
/*------------------------------------------- increment external energy */
dynvar->eout += deltaa;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_eout */



/*----------------------------------------------------------------------*
 |  set time integration constants for Gen alfa method       m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_setconstants(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, DOUBLE dt)
{
DOUBLE *constants;
DOUBLE  beta;
DOUBLE  gamma;
DOUBLE  alpham;
DOUBLE  alphaf;
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
 |  set time integration constants for centr. diff method    m.gee 05/02|
 *----------------------------------------------------------------------*/
void dyn_setconstants_expl(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, DOUBLE dt)
{
DOUBLE *constants;
#ifdef DEBUG 
dstrc_enter("dyn_setconstants_expl");
#endif
/*----------------------------------------------------------------------*/
constants = dynvar->constants;
/*----------------------------------------------------------------------*/
constants[0]  = 1.0/(dt*dt);
constants[1]  = 2.0 * constants[0];
constants[2]  = 1.0 / constants[1];
constants[3]  = 1.0 / (2.0*dt);
constants[4]  = dt;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_setconstants_expl */




/*----------------------------------------------------------------------*
 |  update the vectors                                       m.gee 03/02|
 | displacements sol_new -> sol_old                                     |
 |               rhs_new -> rhs_old                                     |
 | new accelerations at time t                                          |
 | new velocities    at time t                                          |
 |                                                                      |
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

DOUBLE a1,a2,a3,a4,a5,a6;

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
solserv_cpdirich(actfield,0,0,7,8);               /* make a copy of sol[7] in sol[8] */              
solserv_zerodirich(actfield,0,0,7);               /* init sol[7] to zero             */              
solserv_assdirich_fac(actfield,0,0,4,3,7,a1,-a1); /* sol[7]+=a1*sol[4]-a1*sol[3]     */
solserv_assdirich_fac(actfield,0,0,6,0,7,a2,0.0); /* sol[7]+=a2*sol[6]               */
solserv_assdirich_fac(actfield,0,0,8,0,7,a3,0.0); /* sol[7]+=a3*sol[8]               */

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
solserv_adddirich(actfield,0,0,6,0,6,a5,0.0);     /* sol[6] = sol[6] * a5                 */ 
solserv_assdirich_fac(actfield,0,0,4,3,6,a4,-a4); /* sol[6] += a4 * sol[4] + -a4 * sol[3] */
solserv_assdirich_fac(actfield,0,0,8,0,6,a6,0.0); /* sol[6] += a6 * sol[8] (=sol[7]_old)  */

/*--------------------------------------- zero the working place sol[8] */
solserv_zerodirich(actfield,0,0,8);  

/*--------------------------------- update the prescribed displacements */
/*------------------------------------------------ copy u(t) -> u(t-dt) */
/*                                                copy sol[4] to sol[3] */
solserv_adddirich(actfield,0,0,4,0,3,1.0,0.0);  
/*----------------------------------------------------------- zero u(t) */
/*                                                          zero sol[4] */
solserv_zerodirich(actfield,0,0,4);
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
DOUBLE  beta;
DOUBLE  gamma;
DOUBLE  alpham;
DOUBLE  alphaf;
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
 |  print head of nln structural dynamics (GEMM)              m.gee 02/02|
 *----------------------------------------------------------------------*/
#ifdef GEMM
void dyn_nlngemm_outhead(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn)
{
DOUBLE  beta;
DOUBLE  gamma;
DOUBLE  alpham;
DOUBLE  alphaf;
DOUBLE  xsi;
#ifdef DEBUG 
dstrc_enter("dyn_nlngemm_outhead");
#endif
/*----------------------------------------------------------------------*/
beta      = sdyn->beta;
gamma     = sdyn->gamma;
alpham    = sdyn->alpha_m;
alphaf    = sdyn->alpha_f;
xsi       = sdyn->xsi;
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
fprintf(allfiles.out_err,"       Nonlinear Structural Dynamics\n");
fprintf(allfiles.out_err,"\n");
fprintf(allfiles.out_err,"       Generalized Energy-Momentum Method\n");
fprintf(allfiles.out_err,"  beta=%f gamma=%f alpha_m=%f alpha_f=%f xsi=%f\n",beta,gamma,alpham,alphaf,xsi);
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
printf("--------------------------------------------------------------------------\n");
printf("       Nonlinear Structural Dynamics\n");
printf("\n");
printf("       Generalized Energy-Momentum Method\n");
printf("  beta=%f gamma=%f alpha_m=%f alpha_f=%f xsi=%f\n",beta,gamma,alpham,alphaf,xsi);    
printf("--------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nlngemm_outhead */
#endif

/*----------------------------------------------------------------------*
 |  print head of nln structural explicit dynamics           m.gee 05/02|
 *----------------------------------------------------------------------*/
void dyn_nlnstru_outhead_expl()
{
#ifdef DEBUG 
dstrc_enter("dyn_nlnstru_outhead_expl");
#endif
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
fprintf(allfiles.out_err,"       Nonlinear Structural Dynamics\n");
fprintf(allfiles.out_err,"\n");
fprintf(allfiles.out_err,"       Central Difference Method\n");
fprintf(allfiles.out_err,"----------------------------------------------------------------------\n");
printf("--------------------------------------------------------------------------\n");
printf("       Nonlinear Structural Dynamics\n");
printf("\n");
printf("       Central Difference Method\n");
printf("--------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nlnstru_outhead_expl */



/*----------------------------------------------------------------------*
 |  print head of nln structural dynamics                    m.gee 02/02|
 *----------------------------------------------------------------------*/
void dyn_nlnstruct_outstep(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, INT numiter, DOUBLE dt)
{
#ifdef DEBUG 
dstrc_enter("dyn_nlnstruct_outstep");
#endif
/*----------------------------------------------------------------------*/
printf("STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | NUMITER=%3d | ETOT=%-14.8E | \n",
       sdyn->step,sdyn->nstep,sdyn->time,dt,numiter,dynvar->etot);
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | NUMITER=%3d | ETOT=%-14.8E | EPOT=%-14.8E | EKIN=%-14.8E | EOUT=%-14.8E\n",
                         sdyn->step,sdyn->nstep,sdyn->time,dt,numiter,dynvar->etot,dynvar->epot,dynvar->ekin,dynvar->eout);
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
void assemble_dirich_dyn(ELEMENT *actele, ARRAY *estif_global, 
                         ARRAY *emass_global, CONTAINER *container)
{
INT                   i,j;
INT                   counter,hasdirich;
INT                   numdf;
INT                   nd=0;
INT                   idamp=0;
DOUBLE                mdamp,kdamp;
DOUBLE              **estif;
DOUBLE              **emass;
DOUBLE                dirich[MAXDOFPERELE];
DOUBLE                dforces[MAXDOFPERELE];
INT                   dirich_onoff[MAXDOFPERELE];
INT                   lm[MAXDOFPERELE];
GNODE                *actgnode;
#ifdef DEBUG 
dstrc_enter("assemble_dirich_dyn");
#endif
/*----------------------------------------------------------------------*/
/*---------- check presence of any dirichlet conditions to this element */
hasdirich=0;
for (i=0; i<actele->numnp; i++)
{
   if (actele->node[i]->gnode->dirich != NULL) 
   {
      hasdirich=1;
      break;
   }
}
/*--------------------- there are no dirichlet conditions here so leave */
if (hasdirich==0) goto end;
/*--------------------------------------- check for presence of damping */
if (ABS(container->dirichfacs[7]) > EPS13 || ABS(container->dirichfacs[8]) > EPS13) 
{
   idamp=1;
   mdamp = container->dirichfacs[7];
   kdamp = container->dirichfacs[8];
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
      dforces[i] += estif[i][j] * dirich[j] * container->dirichfacs[6];
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
      dforces[i] += emass[i][j] * dirich[j] * container->dirichfacs[0];
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
      dforces[i] += emass[i][j] * dirich[j] * container->dirichfacs[1];
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
      dforces[i] += emass[i][j] * dirich[j] * container->dirichfacs[2];
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
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * container->dirichfacs[3];
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
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * container->dirichfacs[4];
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
      dforces[i] += (mdamp*emass[i][j]+kdamp*estif[i][j]) * dirich[j] * container->dirichfacs[5];
   }/* loop j over columns */
}/* loop i over rows */
/*----------------------------------------------------------------------*/
} /* end of if (idamp) */




/*-------- now assemble the vector dforces to the global vector fullvec */
for (i=0; i<nd; i++)
{
   /*if (lm[i] >= dim) continue;
     fullvec[lm[i]] += dforces[i];*/
   if (lm[i] >= container->global_numeq) continue;
   container->dirich[lm[i]] += dforces[i];
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_dirich_dyn */



/*----------------------------------------------------------------------*
 |  make eff. left hand side for central diff method         m.gee 05/02|
 *----------------------------------------------------------------------*/
void dyn_keff_expl(INTRA *actintra,
                   SPARSE_TYP *sysarray_typ, SPARSE_ARRAY *sysarray,
                   INT stiff_array, INT mass_array, INT damp_array,
                   STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn)
{
DOUBLE a1;

#ifdef DEBUG 
dstrc_enter("dyn_keff_expl");
#endif
/*----------------------------------------------------------------------*/
a1 = dynvar->constants[4]/2.0;
/*----------------------------------------------------------------------*/
if (damp_array>0)
{
   solserv_add_mat(actintra,&sysarray_typ[stiff_array],&sysarray[stiff_array],
                            &sysarray_typ[damp_array],&sysarray[damp_array],
                            a1);
}
   solserv_add_mat(actintra,&sysarray_typ[stiff_array],&sysarray[stiff_array],
                            &sysarray_typ[mass_array],&sysarray[mass_array],
                            1.0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_keff_expl */


/*----------------------------------------------------------------------*
 |  make eigenanalysis of dynamic system                     m.gee 06/02|
 *----------------------------------------------------------------------*/
void dyn_eigen(FIELD *actfield, PARTITION *actpart, SOLVAR *actsolv,
               INTRA *actintra, INT stiff, INT mass)
{
INT       i,j,k;
SPOOLMAT *spostiff;
SPOOLMAT *spomass;
INT       numeq;
INT       itype=1;
char      jobz[1];
char      uplo[1];
ARRAY     A_a,B_a;
ARRAY     EW_a;
DOUBLE  **A,**B, *EW;
INT       lwork;
ARRAY     WORK_a;
DOUBLE   *WORK;
INT       info=1;
INT      *irn,*jcn, nnz;
#ifdef DEBUG 
dstrc_enter("dyn_eigen");
#endif
/*----------------------------------------------------------------------*/
if (actsolv->sysarray_typ[stiff] != spoolmatrix)
   dserror("Eigenanalysis only with spooles");
/*-------------------------------------------------------some variables */
spostiff = actsolv->sysarray[stiff].spo;
spomass  = actsolv->sysarray[mass].spo;
numeq    = spostiff->numeq_total;
lwork    = (INT)(3.2*numeq);
jobz[0]  = 'N';
uplo[0]  = 'L';
A    = amdef("A",&A_a,numeq,numeq,"DA");
B    = amdef("B",&B_a,numeq,numeq,"DA");
EW   = amdef("EW",&EW_a,numeq,1,"DV");
WORK = amdef("WORK",&WORK_a,lwork,1,"DV");
amzero(&A_a);
amzero(&B_a);
irn  = spostiff->irn_loc.a.iv;
jcn  = spostiff->jcn_loc.a.iv;
nnz  = spostiff->irn_loc.fdim;
/*--- fill the dense arrays A and B with the sparse stiffness and mass */
for (k=0; k<nnz; k++)
{
   i = irn[k];
   j = jcn[k];
   A[i][j] = spostiff->A_loc.a.dv[k];
   B[i][j] = spomass->A_loc.a.dv[k];
}
/*--------------------------- call lapack to calculate the eigenvalues */
printf("Reached eigensolver\n");
fflush(stdout);
dsygv(&itype,jobz,uplo,&numeq,A[0],&numeq,B[0],&numeq,EW,WORK,&lwork,&info);
printf("info=%d\n",info);
printf("Largest Eigenvalue is %30.15E\n",EW[numeq-1]);
printf("Time step should then be smaller then %30.15E\n",2.0/(sqrt(EW[numeq-1])));
fflush(stdout);

/*----------------------------------------------------------------------*/
amdel(&A_a);
amdel(&B_a);
amdel(&EW_a);
amdel(&WORK_a);
dserror("Eigensolution finished "); 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_eigen */



/*----------------------------------------------------------------------*
 |  make eigenanalysis of system matrice                     m.gee 06/02|
 *----------------------------------------------------------------------*/
void solserv_eigen(FIELD *actfield, PARTITION *actpart, SOLVAR *actsolv,
               INTRA *actintra, INT stiff)
{
INT       i,j,k;
SPOOLMAT *spostiff;
INT       numeq;
char      jobz[1];
char      uplo[1];
ARRAY     A_a;
ARRAY     EW_a;
DOUBLE  **A, *EW;
INT       lwork;
ARRAY     WORK_a, IWORK_a;
DOUBLE   *WORK;
INT      *IWORK,liwork;
INT       info=1;
INT      *irn,*jcn, nnz;
DIST_VECTOR    *distvecs;
#ifdef DEBUG 
dstrc_enter("solserv_eigen");
#endif
/*----------------------------------------------------------------------*/
if (actsolv->sysarray_typ[stiff] != spoolmatrix)
   dserror("Eigenanalysis only with spooles");
/*-------------------------------------------------------some variables */
spostiff = actsolv->sysarray[stiff].spo;
numeq    = spostiff->numeq_total;
solserv_create_vec(&distvecs,numeq,numeq,numeq,"DV");
for (k=0; k<numeq; k++) solserv_zero_vec(&(distvecs[k]));
lwork    = 1 + 6*numeq + 2*numeq*numeq;
liwork   = 3 + 5*numeq;
jobz[0]  = 'V';
uplo[0]  = 'L';
A     = amdef("A",&A_a,numeq,numeq,"DA");
EW    = amdef("EW",&EW_a,numeq,1,"DV");
WORK  = amdef("WORK",&WORK_a,lwork,1,"DV");
IWORK = amdef("IWORK",&IWORK_a,liwork,1,"IV");
amzero(&A_a);
irn  = spostiff->irn_loc.a.iv;
jcn  = spostiff->jcn_loc.a.iv;
nnz  = spostiff->irn_loc.fdim;
/*--- fill the dense arrays A and B with the sparse stiffness and mass */
for (k=0; k<nnz; k++)
{
   i = irn[k];
   j = jcn[k];
   A[i][j] = spostiff->A_loc.a.dv[k];
} 
/*--------------------------- call lapack to calculate the eigenvalues */
printf("Reached eigensolver\n");
fflush(stdout);
dsyevd(jobz,uplo,&numeq,A[0],&numeq,EW,WORK,&lwork,IWORK,&liwork,&info);
printf("info=%d\n",info);
printf("Largest Eigenvalue is %30.15E\n",EW[numeq-1]);
printf("Condition number is %30.15E\n",EW[numeq-1]/EW[0]);
printf("Time step should then be smaller then %30.15E\n",2.0/(sqrt(EW[numeq-1])));
fflush(stdout);
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   for (k=0; k<numeq; k++)
   distvecs[i].vec.a.dv[k] = A[i][k];
}
for (i=0; i<numeq; i++)
solserv_result_total(actfield,actintra, &(distvecs[i]),i,
                     &(actsolv->sysarray[0]),
                     &(actsolv->sysarray_typ[0]));
out_gid_sol("eigenmodes",actfield,actintra,0,numeq,ZERO);
/*----------------------------------------------------------------------*/
fprintf(allfiles.out_err,"------------------Eigenanalysis of SYMMETRIC (K-lambda*I)*phi=0-----\n");
fprintf(allfiles.out_err,"Eigenvalues in ascending order:\n");
for (k=0; k<numeq; k++)
fprintf(allfiles.out_err," %d %40.20f \n",k+1,EW[k]);
fprintf(allfiles.out_err,"------------------End Eigenanalysis---------------------------------\n");
fflush(allfiles.out_err);
/*----------------------------------------------------------------------*/
amdel(&A_a);
amdel(&EW_a);
amdel(&WORK_a);
amdel(&IWORK_a);
dserror("Eigensolution finished "); 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of solserv_eigen */


/*----------------------------------------------------------------------*
 |  Calculate total energy(for GEMM)                         m.gee 06/02|
 *----------------------------------------------------------------------*/
#ifdef GEMM
void total_energy(PARTITION *actpart,INTRA *actintra, STRUCT_DYN_CALC *dynvar)
{
INT          i;
ELEMENT     *actele;

#ifdef PARALLEL
DOUBLE       send[5],recv[5];
#endif

#ifdef DEBUG 
dstrc_enter("total_energy");
#endif

/*----------------------------------------------------------------------*/
dynvar->total_strain_energy    = 0.0;
dynvar->total_kinetic_energy   = 0.0;
dynvar->total_angular_momentum = 0.0;
dynvar->total_linmom[0] = 0.0;
dynvar->total_linmom[1] = 0.0;
dynvar->local_strain_energy    = 0.0;
dynvar->local_kinetic_energy   = 0.0;
dynvar->local_angular_momentum = 0.0;
dynvar->local_linmom[0] = 0.0;
dynvar->local_linmom[1] = 0.0;

/*----------------------------------------------------------------------*/

for (i=0; i<actpart->pdis[0].numele; i++)
{
  actele = actpart->pdis[0].element[i];
  
  if (actele->eltyp != el_wall1) continue;
  if (actele->proc != actintra->intra_rank) continue;
  dynvar->local_strain_energy    = dynvar->local_strain_energy  + actele->e.w1->strain_energy;
  dynvar->local_kinetic_energy   = dynvar->local_kinetic_energy + actele->e.w1->kinetic_energy;
  dynvar->local_linmom[0]        = dynvar->local_linmom[0] +  actele->e.w1->linmom[0];
  dynvar->local_linmom[1]        = dynvar->local_linmom[1] +  actele->e.w1->linmom[1];
  dynvar->local_angular_momentum = dynvar->local_angular_momentum + actele->e.w1->angular_momentum;  
}
#ifdef PARALLEL
send[0] = dynvar->local_strain_energy;
send[1] = dynvar->local_kinetic_energy;
send[2] = dynvar->local_linmom[0];
send[3] = dynvar->local_linmom[1];
send[4] = dynvar->local_angular_momentum;
MPI_Allreduce (send,recv,5,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
recv[0] = send[0];
recv[1] = send[1];
recv[2] = send[2];
recv[3] = send[3];
recv[4] = send[4];
#endif
  
#ifdef DEBUG 
dstrc_exit();
#endif
return;
}/* end of total_energy */
#endif/*GEMM*/
