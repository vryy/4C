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
 |  compute system energy                                    m.gee 02/02|
C     ******************************************************************
C     *                                                                *
C     *                       +---------------+                        *
C     *                       I    D Y N E    I                        *
C     *                       +---------------+                        *
C     *                                                                *
C     *                         SYSTEM ENERGY                          *
C     *     **  EKIN = 1/2*VEL^T*M*VEL                         **      *
C     *     **  ETOT = 1/2*VEL^T*M*VEL + 1/2*DISP^T*K*DISP     **      *
C     ******************************************************************
   in here, only kinetic energy is computed
 *----------------------------------------------------------------------*/
int dyne(STRUCT_DYN_CALC *dynvar,
         INTRA           *actintra,
         FIELD           *actfield,
         SOLVAR          *actsolv,
         int              stiff_array,
         int              mass_array,
         DIST_VECTOR     *vel,
         DIST_VECTOR     *disp,
         DIST_VECTOR     *rhs,
         int              actstep,
         int              nstep
        )
{
DIST_VECTOR *help;
int          numeq,numeq_total;

#ifdef DEBUG 
dstrc_enter("dyne");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------- create a working distvec */
numeq       = vel->numeq;
numeq_total = vel->numeq_total;
solserv_create_vec(&help,1,numeq_total,numeq,"DV");
/*------------------------------------------------- make kinetic energy */
/*------------------------------------------- perform help = mass * vel */
solserv_sparsematvec(actintra,
                     &help[0],
                     &(actsolv->sysarray[mass_array]),
                     &(actsolv->sysarray_typ[mass_array]),
                     vel);
/*------------------------------------------------ perform vel^T * help */
solserv_dot_vec(actintra,vel,help,&(dynvar->ekin));
dynvar->ekin /= 2.0;
/*-------------------------------------------- delete the working array */
solserv_del_vec(&help,1);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyne */




/*----------------------------------------------------------------------*
 |  set time integration constants                           m.gee 02/02|
 *----------------------------------------------------------------------*/
int dyn_setconstants(STRUCT_DYN_CALC *dynvar, STRUCT_DYNAMIC *sdyn, double dt)
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
