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
/*----------------------------------------allocate 2 dist. load vectors */
actsolv->nrhs = 2;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));
/*------------------------------------------ number of solution vectors */
/*-------------------- there is one solution vector to hold total displ.*/
actsolv->nsol= 1;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nsol; i++) solserv_zero_vec(&(actsolv->sol[i]));
/*-------------------------------------------- allocate one vector vel */
solserv_create_vec(&vel,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(vel[i]));
/*-------------------------------------------- allocate one vector acc */
solserv_create_vec(&acc,1,numeq_total,numeq,"DV");
for (i=0; i<1; i++) solserv_zero_vec(&(acc[i]));
/*------------------------ initialize solver on the matrix stiff_array */
/*
NOTE: solver init phase has to be called with each matrix one wants to solve with
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(  
                    actsolv,
                    actintra,
                  &(actsolv->sysarray_typ[stiff_array]),
                  &(actsolv->sysarray[stiff_array]),
                  &(actsolv->sol[0]),
                  &(actsolv->rhs[0]),
                    init
              );
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
          &(actsolv->rhs[0]),
          &(actsolv->rhs[1]),
          0,
          action
      );
/*--------------------------------------------- add the two rhs vectors */
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]));
solserv_copy_vec(&(actsolv->rhs[0]),&(actsolv->rhs[1]));
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
solserv_scalarprod_vec(&(actsolv->rhs[0]),dynvar.rldfac);
/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[0]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actfield,actsolv,stiff_array,mass_array,&vel[0],
     &(actsolv->sol[0]),&(actsolv->rhs[0]),0,sdyn->nstep);
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);
/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains(actfield);
}


























/*----------------------------------------------------------------------*/
end:
#ifndef PARALLEL 
FREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */




