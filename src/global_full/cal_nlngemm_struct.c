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
#ifdef GEMM
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

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

/*!------------------------------------------------------------------------
\brief main structure for 2-D contact (bilinear discretization) m.gee 10/02

<pre>
defined in wall_contact_detection.c                             
</pre>

*-------------------------------------------------------------------------*/
#ifdef WALLCONTACT
extern struct _WALL_CONTACT contact;
#endif

extern enum _CALC_ACTION calc_action[MAXFIELD];


DOUBLE acttime;

/*----------------------------------------------------------------------------*
 |  routine to control nonlinear dynamic structural analysis(GEMM) m.gee 02/02|
 *---------------------------------------------------------------------------*/
void dyn_nln_gemm() 
{
INT             i;                  /* simply a counter */
INT             numeq;              /* number of equations on this proc */
INT             numeq_total;        /* total number of equations */
INT             init;               /* flag for solver_control call */
INT             itnum,itstore;      /* counter for NR-Iterations */
INT             convergence;        /* convergence flag */
INT             mod_disp,mod_stress;
INT             mod_res_write;
INT             updevry_disp;
INT             restart;
DOUBLE          maxtime;
DOUBLE          t0_res,t1_res;

DOUBLE          dt;
INT             nstep;
DOUBLE          deltaepot=0.0;

DOUBLE          t0,t1;

DOUBLE          dmax;               /* infinity norm of residual displacements */

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
DIST_VECTOR    *midcon;             /* distributed vector for generalized mid-configuration*/
INT ttt;                            /* counter for the output file *.gemm */   
ARRAY           intforce_a;         /* redundant vector of full length for internal forces */
DOUBLE         *intforce;
ARRAY           dirich_a;           /* redundant vector of full length for dirichlet-part of rhs */
DOUBLE         *dirich;             
DOUBLE          dirichfacs[10];     /* factors needed for dirichlet-part of rhs */
 
STRUCT_DYN_CALC dynvar;             /* variables to perform dynamic structural simulation */              

INT             contactflag = 0;    /* flag for contact onoff */
ARRAY           contactforce_a;     /* redundant vector of full length for contact forces */
DOUBLE         *cforce;
DIST_VECTOR    *con;                /*  contact forces */              
INT             augon=0;
INT             actaug=0;
INT             aug_number=20;       /*  # of augmentation loops */
DOUBLE          contactdt;

INT             timeadapt;          /* flag to switch time adaption on/off */
INT             itwant;
DOUBLE          maxdt;
DOUBLE          resultdt;
DOUBLE          newdt;
DOUBLE          olddt;
DOUBLE          eta;
INT             repeatcount;
DOUBLE          lowtime,uptime,writetime;
DOUBLE          remain;
DOUBLE          low,up;
INT             ilow,iup;
DOUBLE          tau,tau2,tau3,fac;

CONTAINER       container;          /* contains variables defined in container.h */
container.isdyn = 1;                /* dynamic calculation */
container.actndis = 0;              /* only one discretisation */

#ifdef DEBUG 
dstrc_enter("dyn_nln_gemm");
#endif
/*----------------------------------------------------------------------*/
restart = genprob.restart;
/*--------------------------------------------------- set some pointers */
actfield           = &(field[0]);
actsolv            = &(solv[0]);
actpart            = &(partition[0]);
action             = &(calc_action[0]);
sdyn               = alldyn[0].sdyn;
contactflag        = sdyn->contact;
container.fieldtyp = actfield->fieldtyp;
timeadapt          = sdyn->timeadapt;
itwant             = sdyn->itwant;
maxdt              = sdyn->maxdt;
resultdt           = sdyn->resultdt;
sdyn->writecounter = 0;
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
dynvar.rldfac       = 0.0;
dynvar.rnorm        = 0.0;
dynvar.epot         = 0.0;
dynvar.eout         = 0.0;
dynvar.etot         = 0.0;
dynvar.ekin         = 0.0;
dynvar.dinorm       = 0.0;
dynvar.dnorm        = 0.0;

/*--------------------------------------------Initialization of strain energy*/
dynvar.total_strain_energy = 0.0;

/*-----------------------------------Write out the initial energy and momenta*/
   ttt = 0;
   fprintf(allfiles.out_gemm, " step       Total En.         Strain           Kinetic            Lx               Ly             AngMom \n");
   total_energy(actpart,actintra,&dynvar);
   fprintf(allfiles.out_gemm,"%5d %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f\n", ttt, 
           dynvar.total_strain_energy+dynvar.total_kinetic_energy,
           dynvar.total_strain_energy, dynvar.total_kinetic_energy, 
	   dynvar.total_linmom[0], dynvar.total_linmom[1], dynvar.total_angular_momentum);
   ttt++;
   fflush(allfiles.out_gemm);   
/*--------------------------------------------------------------------------*/ 

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
solserv_zero_mat(actintra,&(actsolv->sysarray[i]),&(actsolv->sysarray_typ[i]));

/*---------------------------- get global and local number of equations */
solserv_getmatdims(&(actsolv->sysarray[stiff_array]),actsolv->sysarray_typ[stiff_array],
                   &numeq,&numeq_total);

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
/*------------ there are two vectors for generalized mid-configurations*/
solserv_create_vec(&midcon,2,numeq_total,numeq,"DV");
for (i=0; i<2; i++) solserv_zero_vec(&(midcon[i]));

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
/*------------------ create a vector of full length for contact forces */
if (contactflag)
{
   cforce = amdef("contact",&contactforce_a,numeq_total,1,"DV");
   solserv_create_vec(&con,1,numeq_total,numeq,"DV");
   solserv_zero_vec(&(con[0]));
}
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
solver_control(actsolv, actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),
               &(dispi[0]),&(actsolv->rhs[0]),init);

solver_control(actsolv, actintra,&(actsolv->sysarray_typ[mass_array]),&(actsolv->sysarray[mass_array]),
               &work[0],&work[1],init);

if (damp_array>0)
solver_control(actsolv, actintra, &(actsolv->sysarray_typ[damp_array]),&(actsolv->sysarray[damp_array]),
               &work[0],&work[1],init);
/*----------------- init the assembly for stiffness and for mass matrix */
/*                                           (damping is not assembled) */
init_assembly(actpart,actsolv,actintra,actfield,stiff_array,0);
init_assembly(actpart,actsolv,actintra,actfield,mass_array,0);

/*------------------------------- init the element calculating routines */
*action = calc_struct_init;
calinit(actfield,actpart,action,&container);

/*----------------------------------------- write output of mesh to gid */
if (par.myrank==0)
if (ioflags.struct_disp_gid||ioflags.struct_stress_gid) 
   out_gid_msh();

/*----------- init the contact algorithms for contact with wall element */
#ifdef WALLCONTACT
if (contactflag)
wallcontact_init(actfield);
#endif
/*----------------------- call elements to calculate stiffness and mass */
*action = calc_struct_nlnstiffmass;
container.dvec          = NULL;
container.dirich        = NULL;
container.global_numeq  = 0;
container.dirichfacs    = NULL;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);

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
/*----------------------- init the time curve applied to the loads here */
/*-------------- this control routine at the moment always uses curve 0 */
/*-------------------------------------------------- init the timecurve */
actcurve = 0;
dyn_init_curve(actcurve,sdyn->nstep,sdyn->dt,sdyn->maxtime);

/* put a zero to the place 11 in node->sol to init the velocities and accels */
/* of prescribed displacements */
solserv_putdirich_to_dof(actfield,0,0,0.0,12);

/*-------------------- put a zero to the place 1 and 2 in sol_increment */
/*------------------ later this will hold internal forces at t and t-dt */
solserv_sol_zero(actfield,0,1,2);
solserv_sol_zero(actfield,0,1,1);

/*------------------------------------------- set initial step and time */
sdyn->step = -1;
sdyn->time = 0.0;

/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains(actfield);
}
/*------------------------------------------------------- printout head */
if (par.myrank==0) dyn_nlngemm_outhead(&dynvar,sdyn);

/*----------------------------------------------------------------------*/
/*                     START LOOP OVER ALL STEPS                        */
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

Values of the different vectors from above in one loop:
  /......no change in this step
  =,+=...evaluation in this step

   vector	Predictor - Start     Precictor - End	  Corrector - Start	Corrector - End			Update - End

   rhs[3]  	/{=orig. load vect.} 	/	 		/		/				/
   rhs[2]       /{=rhs(t-dt)}     	/			/		/				=rhs[1]{=rhs(t)}
   rhs[1]       =rhs(t)      		/			/		/				/
   rhs[0]    	/{=rhs(t-2dt)}		=feff_p			/		=feff_c				=rhs[2]{=rhs(t-dt)}

   fie[2]    	/			/			=fint(t)	/				/
   fie[1]    	=fint(t-dt)		/			/		/				/
   fie[0]    	/			/			/		=(1-alpha_f)*fie[2]+alpha_f*fie[1]	/

   dispi[0]  	=0			=Keff^-1*feff-p		/		+=work[0]			/

   sol[0]    	/{=d(t-dt)}		/			/		/				=sol[1]{=d(t)}
   sol[1]    	{=d(t-dt)}		=sol[0]+dispi[0]{=d(t)}	/		=sol[0]+dispi[0]		/
   
   vel[0]    	/{=v(t-dt)}		/			/		/				=v(t)
   acc[0]    	/{=a(t-dt)}		/			/		/				=a(t)	

   work[2]    	/{=v(t-2dt)}		/			/		/				=v(t-dt)
   work[1]    	/{=a(t-2dt)}		/			/		/				=a(t-dt)
   work[0]    	/ 			/			/		=Keff^-1*feff-c			=M*vel[0]

*/
timeloop:
t0 = ds_cputime();
/*------------------------------------------------- write memory report */
/*if (par.myrank==0) dsmemreport();*/
/*--------------------------------------------------- check for restart */
if (restart)
{
   t0_res = ds_cputime();
   /*-------------- save the stepsize as it will be overwritten in sdyn */
   dt    = sdyn->dt;
   /*------ save the number of steps, as it will be overwritten in sdyn */
   nstep = sdyn->nstep;
   maxtime = sdyn->maxtime;
   /*------------- save the restart interval, as it will be overwritten */
   mod_res_write = sdyn->res_write_evry;
   updevry_disp  = sdyn->updevry_disp;
   /*----------------------------------- the step to read in is restart */
   restart_read_nlnstructdyn(restart,sdyn,&dynvar,actfield,actpart,actintra,action,
                             actsolv->nrhs, actsolv->rhs,
                             actsolv->nsol, actsolv->sol,
                             1            , dispi       ,
                             1            , vel         ,
                             1            , acc         ,
                             3            , fie         ,
                             3            , work        ,
                             &intforce_a,
                             &dirich_a,
                             &container);     /* contains variables defined in container.h */
   /*-------------------------------------- put the dt to the structure */
   sdyn->dt = dt;
   /*--------------------------------------- put nstep to the structure */
   sdyn->nstep = nstep;
   sdyn->maxtime = maxtime;
   /*-------------------------------- put restart interval to structure */
   sdyn->res_write_evry = mod_res_write;
   sdyn->updevry_disp   = updevry_disp;
   /*------------------------------------------- switch the restart off */
   restart=0;
   /*----------------------------------------------------- measure time */
   t1_res = ds_cputime();
   fprintf(allfiles.out_err,"TIME for restart reading is %f sec\n",t1_res-t0_res);
}
/*--------------------------------------------- increment step and time */
sdyn->step++;
repeatcount = 0;
/*------------------- modifications to time steps size can be done here */
adapt_top:
/*------------------------------------------------ set new absolue time */
sdyn->time += sdyn->dt;
/*--- put time to global variable for time-dependent load distributions */
acttime = sdyn->time;
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);

/*store displacement of step n at midcon[1] and ARRAY.sol.a.da[13][....] */
solserv_zero_vec(&(midcon[1]));
solserv_copy_vec(&(actsolv->sol[1]),&midcon[1]);
solserv_result_total(actfield,actintra, &midcon[1],13,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*----------------------------------------------------------------------*/
/*---------------------- set incremental displacements dispi[0] to zero */
solserv_zero_vec(&dispi[0]);

/*------------------------- set residual displacements in nodes to zero */
solserv_result_resid(actfield,actintra,&dispi[0],0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));

/*----------------------------------------------------------------------*/
/*                     PREDICTOR                                        */
/*----------------------------------------------------------------------*/
/*---------------------- this vector holds loads due to external forces */
solserv_zero_vec(&(actsolv->rhs[1]));
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[1]),action,&container);
/*------------------------------------------------ get factor at time t */
dyn_facfromcurve(actcurve,sdyn->time,&(dynvar.rldfac));

/*------------------------ multiply rhs[1] by actual load factor rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,0);

/* put the prescribed scaled displacements to the nodes in field sol at */
/*                                  place 4 separate from the free dofs */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,4);

/*-------- put presdisplacements(t) - presdisplacements(t-dt) in place 5 */
solserv_adddirich(actfield,0,0,3,4,5,-1.0,1.0);

/*----- set factors needed for prescribed displacement terms on rhs eff */
/*
dirichfacs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))         
dirichfacs[1] =  (1.0-alpham)*(1.0/beta)/dt                 
dirichfacs[2] =  (1.0-alpham)/(2*beta) - 1                  
dirichfacs[3] = -(1.0-alphaf)*(gamma/beta)/dt               
dirichfacs[4] =  (1.0-alphaf)*gamma/beta - 1                
dirichfacs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)            
dirichfacs[6] = -(1.0-alphaf) or 0                          
dirichfacs[7] =  raleigh damping factor for mass            
dirichfacs[8] =  raleigh damping factor for stiffness       
dirichfacs[9] =  dt     
see phd theses Mok page 165: generalized alfa time integration with prescribed displ.                                    
*/
dirichfacs[0] = -dynvar.constants[0];
dirichfacs[1] =  dynvar.constants[1];
dirichfacs[2] =  dynvar.constants[2];
dirichfacs[3] = -dynvar.constants[3];
dirichfacs[4] =  dynvar.constants[4];
dirichfacs[5] =  dynvar.constants[5];
dirichfacs[6] = -dynvar.constants[6]; 
dirichfacs[9] =  sdyn->dt; 
if (damp_array>0) {
   dirichfacs[7] =  sdyn->m_damp;
   dirichfacs[8] =  sdyn->k_damp;}
else {
   dirichfacs[7] =  0.0;
   dirichfacs[8] =  0.0;}

/*- calculate tangential stiffness/mass and internal forces at time t-dt */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),&(actsolv->sysarray_typ[mass_array]));
amzero(&dirich_a);
amzero(&intforce_a);
/*----------------------------------------------------- contact detection */
#ifdef WALLCONTACT
if (contactflag)
{
   amzero(&contactforce_a);
   wall_contact_mid(actfield,actintra);
   wall_contact_detection_em(actfield,actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]),cforce);
}
#endif

/*-------------------------------------------------------- call elements */
*action = calc_struct_nlnstiffmass;
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);

/*---------------------------- store positive internal forces on fie[1] */
solserv_zero_vec(&fie[1]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(fie[1]),intforce,1.0);

/*-------------------------------- put contact forces to internal forces */
#ifdef WALLCONTACT
if (contactflag)
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(fie[1]),cforce,1.0);
#endif

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/*---------- subtract internal forces from interpolated external forces */
solserv_add_vec(&(fie[1]),&(actsolv->rhs[0]),-1.0);

/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);

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
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
/*
  keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
*/  
/*----------------------------------------------------------------------*/
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,damp_array);
/*------------- call for solution of system dispi[0] = Keff^-1 * rhs[0] */
init=0;
solver_control(actsolv, actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),
               &(dispi[0]),&(actsolv->rhs[0]),init);

/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);

/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));

/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]),0);

/* here put incremental prescribed displacements from sol[5] to sol_increment[0] ? */

/*--------- copy predicted displacements and generate generalized mid-configuration*/
solserv_zero_vec(&(midcon[0]));
solserv_copy_vec(&dispi[0],&midcon[0]);
solserv_scalarprod_vec(&midcon[0],sdyn->alpha_f);
solserv_result_total(actfield,actintra, &midcon[0],8,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*---------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/*                     AUGMENTATION FOR CONTACT                         */
/*----------------------------------------------------------------------*/
#ifdef WALLCONTACT
if (contactflag)
augstart:;
#endif

/*----------------------------------------------------------------------*/
/*                     PERFORM EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/
itnum=0;
iterloop:
/*----- set factors needed for prescribed displacement terms on rhs eff */
dirichfacs[0] = -dynvar.constants[0];
dirichfacs[1] =  dynvar.constants[1];
dirichfacs[2] =  dynvar.constants[2];
dirichfacs[3] = -dynvar.constants[3];
dirichfacs[4] =  dynvar.constants[4];
dirichfacs[5] =  dynvar.constants[5];
dirichfacs[6] =  0.0; 
dirichfacs[9] =  sdyn->dt; 
if (damp_array>0) {
   dirichfacs[7] =  sdyn->m_damp;
   dirichfacs[8] =  sdyn->k_damp;}
else {
   dirichfacs[7] =  0.0;
   dirichfacs[8] =  0.0;}
/* zero the stiffness matrix and vector for internal forces and dirichlet forces */
solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),&(actsolv->sysarray_typ[mass_array]));
amzero(&intforce_a);
amzero(&dirich_a);
/*------------------------------------------------------ detect contact */
#ifdef WALLCONTACT
if (contactflag)
{
   amzero(&contactforce_a);
   wall_contact_mid(actfield,actintra);   
   wall_contact_detection_em(actfield,actintra,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]),cforce);
}
#endif

/* call element routines for calculation of tangential stiffness and intforce */
*action = calc_struct_nlnstiffmass;
solserv_sol_zero(actfield,0,1,2);
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);

/*---------------------------- store positive internal forces on fie[2] */
solserv_zero_vec(&fie[2]);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(fie[2]),intforce,1.0);

/*-------------------------------- put contact forces to internal forces */
#ifdef WALLCONTACT
if (contactflag)
{
solserv_zero_vec(&(con[0]));
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(con[0]),cforce,1.0);
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(fie[2]),cforce,1.0);
}
#endif

/* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));

/* ---------------------------------------------- create internal forces */
solserv_copy_vec(&fie[2],&fie[0]);

/*---------------- subtract internal forces from interp. external forces */
solserv_add_vec(&fie[0],&(actsolv->rhs[0]),-1.0);

/*------------------ add dirichlet forces from prescribed displacements */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);

/*--------------------- create effective load vector (rhs[0]-fie[0])eff */
pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,mass_array,damp_array);

/*----------------------------------- create effective stiffness matrix */
kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,damp_array);

/*---------------------------------------- solve keff * rsd[0] = rhs[0] */
/* solve for residual displacements to correct incremental displacements*/
init=0;
solver_control(actsolv, actintra,&(actsolv->sysarray_typ[stiff_array]),&(actsolv->sysarray[stiff_array]),
               &(work[0]),&(actsolv->rhs[0]),init);

/*-------------------------- return residual displacements to the nodes */
solserv_result_resid(actfield,actintra,&work[0],0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));

/*-- update the incremental displacements by the residual displacements */
solserv_add_vec(&work[0],&dispi[0],1.0);

/*------------------------------------------------ update displacements */
/*------------------------------------------ sol[1] = sol[0] + dispi[0] */
solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);

/*----------------------------- return total displacements to the nodes */
solserv_result_total(actfield,actintra, &(actsolv->sol[1]),0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));

/*----------------------- return incremental displacements to the nodes */
solserv_result_incre(actfield,actintra,&dispi[0],0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]),0);

/*--------------------------------- update generalized mid-configuration*/
solserv_zero_vec(&(midcon[0]));
solserv_copy_vec(&dispi[0],&midcon[0]);
solserv_scalarprod_vec(&midcon[0],sdyn->alpha_f);
solserv_result_total(actfield,actintra, &midcon[0],8,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*----------------------------------------------- check for convergence */
convergence = 0;
dmax        = 0.0;
solserv_vecnorm_euclid(actintra,&(work[0]),&(dynvar.dinorm));
solserv_vecnorm_euclid(actintra,&(dispi[0]),&(dynvar.dnorm));
solserv_vecnorm_Linf(actintra,&(work[0]),&dmax);
if (par.myrank==0) printf("                                                   Residual %10.5E\n",dynvar.dinorm);fflush(stdout);
if (dynvar.dinorm < sdyn->toldisp ||
    dynvar.dnorm  < EPS14 ||
   (dynvar.dinorm < EPS14 && dmax < EPS12) )
{
   itnum++;
   convergence = 1;
#ifdef WALLCONTACT
   wall_contact_mid(actfield,actintra);
   wall_contact_update_em_a(actfield,actintra);
#endif
}    
else
{
   itnum++;
   if (itnum==sdyn->maxiter && !timeadapt)
      dserror("No convergence in maxiter steps");   
   goto iterloop;
}
/*----------------------------------------------------------------------*/
/*                      END OF EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                     AUGMENTATION FOR CONTACT                         */
/*----------------------------------------------------------------------*/
#ifdef WALLCONTACT
if (contactflag)
{
  if(actaug >= aug_number) goto augend;
 
  wall_contact_flag(actintra);
  if(aug_number == 1) goto augend;
  if(contact.contactflag == contact_off && actaug == 0) goto augend;
  
  wall_contact_set();
  wall_contact_augmentation_em(actintra);
  actaug++;
  
  contact.contact_set = NULL;
  contact.set_size = 0;   
  goto augstart;
}
#endif

/*------------------------------ make time adaption iteration indicator */
#if 0
adapt_bot:
if (timeadapt)
{
   eta    = ((DOUBLE)itwant)/((DOUBLE)itnum); 
   olddt  = sdyn->dt;
   newdt  = sdyn->dt * pow(eta,0.3333333);
   if (newdt>1.5*(sdyn->dt)) 
      newdt = 1.5*(sdyn->dt);
   if (newdt>maxdt)
      newdt = maxdt;
   /* repeat step, if no convergence was reached */
   if (convergence==0)
   {
      repeatcount++;
      if (repeatcount >= 10) dserror("No convergence after repeating step 10 times");
      if (par.myrank==0) printf("No convergence in maxiter steps - repeat step\n");
      /* reduce the step size drastically */
      eta         = ((DOUBLE)itwant)/40.0; 
      newdt       = sdyn->dt * pow(eta,0.3333333);
      sdyn->time -= sdyn->dt;
      sdyn->dt    = newdt;
      /* zero the incremental displacements */
      solserv_zero_vec(&dispi[0]);
      /* zero the residual displacements */
      solserv_result_incre(actfield,actintra,&dispi[0],0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
      /* set old total displacements */
      solserv_result_total(actfield,actintra,&(actsolv->sol[0]),0,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
      /* set element internal variables bacl to the last step */
      *action = calc_struct_update_stepback;
      container.dvec          = NULL;
      container.dirich        = NULL;
      container.global_numeq  = 0;
      container.dirichfacs    = NULL;
      container.kstep         = 0;
      container.ekin          = 0.0;
      calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);
      goto adapt_top;
   }
#ifdef PARALLEL
   MPI_Bcast(&newdt,1,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);
#endif
   sdyn->dt = newdt;
}
#endif
/*----------------------------------------------------------------------*/
/*                     AUGMENTATION FOR CONTACT                         */
/*----------------------------------------------------------------------*/
#ifdef WALLCONTACT
/*------------------------------------------------- make contact history*/
augend:;
if (contactflag){
actaug = 0;
/*------------------------ write contact forces to the nodes in place 9 */
solserv_scalarprod_vec(&(con[0]),(-1.0));
solserv_result_total(actfield,actintra, &(con[0]),9,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
/*----------------------------------------------------------------------*/
}
#endif

#ifdef WALLCONTACT
/*------------------------------------------------Update contact history*/
wall_contact_update_em(actfield,actintra);
wall_contact_his_update_em(actintra);
CCAFREE(contact.contact_set);
#endif

/*----------- make temporary copy of actsolv->rhs[2] to actsolv->rhs[0] */
/*                                   (load at t-dt)                     */
/* because in  dyn_nlnstructupd actsolv->rhs[2] is overwritten but is   */
/* still needed to compute energies                                     */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
/*------------------------------ copy disp from sol place 0 to place 10 */
solserv_sol_copy(actfield,0,0,0,0,10);
/*------------------------------ copy vels from sol place 1 to place 10 */
solserv_sol_copy(actfield,0,0,0,1,11);
/*------------------------------ copy accs from sol place 2 to place 11 */
solserv_sol_copy(actfield,0,0,0,2,12);
/*------------------ update displacements, velocities and accelerations */
dyn_nlnstructupd(actfield,
                 &dynvar,sdyn,actsolv,
                 &(actsolv->sol[0]),/* total displacements at time t-dt */
                 &(actsolv->sol[1]),/* total displacements at time t    */
                 &(actsolv->rhs[1]),/* load vector         at time t    */
                 &(actsolv->rhs[2]),/* load vector         at time t-dt */
                 &vel[0],           /* velocities          at time t    */
                 &acc[0],           /* accelerations       at time t    */ 
                 &work[0],          /* working arrays                   */
                 &work[1],          /* working arrays                   */
                 &work[2]);         /* working arrays                   */
/*-------------------------------------- return velocities to the nodes */
solserv_result_total(actfield,actintra, &vel[0],1,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
/*-------------------------velocities for prescribed dofs to velocities */
solserv_adddirich(actfield,0,0,6,0,1,1.0,0.0);
/*------------------------------------------ return accel. to the nodes */
solserv_result_total(actfield,actintra, &acc[0],2,&(actsolv->sysarray[stiff_array]),&(actsolv->sysarray_typ[stiff_array]));
/*-------------------------------------------accel. for prescribed dofs */
solserv_adddirich(actfield,0,0,7,0,2,1.0,0.0);
/* 
It is a bit messed up, but anyway:
in the nodes the results are stored the following way: 

in ARRAY sol.a.da[place][0..numdf-1]:
place 0  holds total displacements  time t      (free/prescr) 
place 1  holds velocities           time t      (free/prescr) 
place 2  holds accels               time t      (free/prescr) 
place 3  holds displacements        time t-dt   (prescr only) 
place 4  holds displacements        time t      (prescr only) 
place 5  holds place 4 - place 3
place 6  holds velocities           time t      (prescr only)
place 7  holds accels               time t      (prescr only)
place 8  holds generalized mid-configuration disp.
place 9  holds contact forces       time t      (free only) 
place 10 holds total displacements  time t-dt   (free/prescr) 
place 11 holds velocities           time t-dt   (free/prescr) 
place 12 holds accels               time t-dt   (free/prescr) 

in ARRAY sol_increment.a.da[place][0..numdf-1]
place 0 holds converged incremental displacements (without prescribed dofs) 
place 1 holds converged internal forces at time t-dt 
place 2 holds converged internal forces at time t

in ARRAY sol_residual
place 0 holds residual displacements during iteration (without prescribed dofs)
*/
/*-----------------------update state variables for EM conserving scheme*/
*action = calc_struct_update_istep;
container.dvec          = NULL;
container.dirich        = NULL;
container.global_numeq  = 0;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);
/* Calculate total energy and momenta(for wall element only for GEMM) and write to the output file*/
total_energy(actpart,actintra,&dynvar);
   fprintf(allfiles.out_gemm,"%5d %16.9f %16.9f %16.9f %16.9f %16.9f %16.9f\n", ttt, 
           dynvar.total_strain_energy+dynvar.total_kinetic_energy,
           dynvar.total_strain_energy, dynvar.total_kinetic_energy, 
           dynvar.total_linmom[0], dynvar.total_linmom[1], dynvar.total_angular_momentum);  
           ttt++;
           fflush(allfiles.out_gemm);
/*------------------------------------------------------------------------------------------------*/
/*---------------------- make incremental potential energy at the nodes */
dyn_epot(actfield,0,actintra,&dynvar,&deltaepot);
dynvar.epot += deltaepot;
/*-------------------------------- make kinetic energy at element level */
dyn_ekin(actfield,actsolv,actpart,actintra,action,&container,stiff_array,mass_array);
dynvar.ekin = container.ekin;
/*------------------------------------------------ make external energy */
dyn_eout(&dynvar,sdyn,actintra,actsolv,&dispi[0],&(actsolv->rhs[1]),&(actsolv->rhs[0]),&work[0]);
/*--------------------------------------------------- make total energy */
dynvar.etot = dynvar.epot + dynvar.ekin;
/*------------------------- update the internal forces in sol_increment */
/*       copy from sol_increment.a.da[2][i] to sol_increment.a.da[1][i] */
solserv_sol_copy(actfield,0,1,1,2,1);
/*------------------------------- check whether to write results or not */
mod_disp      = sdyn->step % sdyn->updevry_disp;
mod_stress    = sdyn->step % sdyn->updevry_stress;
/*------------------------------- check whether to write restart or not */
mod_res_write = sdyn->step % sdyn->res_write_evry;
/*------------------------------------------ perform stress calculation */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress_file==1 || ioflags.struct_stress_gid==1)
{
   *action = calc_struct_stress;
   container.dvec          = NULL;
   container.dirich        = NULL;
   container.global_numeq  = 0;
   container.dirichfacs    = NULL;
   container.kstep         = 0;
   calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);
   /*-------------------------- reduce stresses, so they can be written */
   *action = calc_struct_stressreduce;
   container.kstep = 0;
   calreduce(actfield,actpart,actintra,action,&container);
}
/*-------------------------------------------- print out results to out */
if (mod_stress==0 || mod_disp==0)
if (ioflags.struct_stress_file==1 && ioflags.struct_disp_file==1)
{
  out_sol(actfield,actpart,actintra,sdyn->step,0);
}
/*-------------------------- printout results to gid no time adaptivity */
if (timeadapt==0)
if (par.myrank==0) 
{
   if (mod_disp==0)
   if (ioflags.struct_disp_gid==1)
   {
      out_gid_soldyn("displacement",actfield,actintra,sdyn->step,0,sdyn->time);
      /*out_gid_soldyn("velocity",actfield,actintra,sdyn->step,1,sdyn->time);*/
      /*out_gid_soldyn("accelerations",actfield,actintra,sdyn->step,2,sdyn->time);*/
#ifdef WALLCONTACT
      if (contactflag)
      out_gid_soldyn("contact",actfield,actintra,sdyn->step,9,sdyn->time);
#endif

   }
   if (mod_stress==0)
   if (ioflags.struct_stress_gid==1)
   {
      out_gid_soldyn("stress"      ,actfield,actintra,sdyn->step,0,sdyn->time);
   }
}
/*----------------------------- printout results to gid time adaptivity */
#if 0
if (timeadapt)
{
   lowtime = sdyn->time - olddt;
   uptime  = sdyn->time;
   low     = lowtime / resultdt;
   up      =  uptime / resultdt;
   ilow    = (INT)floor(low);
   iup     = (INT)ceil(up);
   remain  = fmod(uptime,iup*resultdt);
   if (remain<EPS12) iup++;
   for (i=ilow+1; i<iup; i++) 
   {
       writetime = i*resultdt; 
       tau       = writetime - lowtime;
       tau2      = tau*tau;
       tau3      = tau2*tau;
       /* u_write = u_old */
       solserv_sol_copy(actfield,0,0,0,10,8);
       /* u_write += udot_old * tau */
       solserv_sol_add (actfield,0,0,0,11,8,tau);
       /* u_write += udotdot_old * fac */
       fac = 0.5*tau2 - (sdyn->beta/olddt)*tau3;
       solserv_sol_add (actfield,0,0,0,12,8,fac);
       /* u_write += udotdot_new * fac */
       fac = (sdyn->beta/olddt)*tau3;
       solserv_sol_add (actfield,0,0,0,2,8,fac);
       /* write displacements u_write */
       if (par.myrank==0)
       out_gid_soldyn("displacement",actfield,actintra,sdyn->writecounter,8,writetime);
       /* write contact forces */
       if (contactflag && par.myrank==0)
       out_gid_soldyn("contact",actfield,actintra,sdyn->writecounter,9,writetime);
       sdyn->writecounter++;
   }
}
#endif
/*-------------------------------------- write restart data to pss file */
if (mod_res_write==0)
{
restart_write_nlnstructdyn(sdyn,&dynvar,actfield,actpart,actintra,action,
                           actsolv->nrhs, actsolv->rhs,
                           actsolv->nsol, actsolv->sol,
                           1            , dispi       ,
                           1            , vel         ,
                           1            , acc         ,
                           3            , fie         ,
                           3            , work        ,
                           &intforce_a,
                           &dirich_a,
                           &container);     /* contains variables defined in container.h */
}                           
/*----------------------------------------------------- print time step */
#if 0
if (par.myrank==0 &&  timeadapt)  dyn_nlnstruct_outstep(&dynvar,sdyn,itnum,olddt);
#endif
if (par.myrank==0 && !timeadapt) dyn_nlnstruct_outstep(&dynvar,sdyn,itnum,sdyn->dt);
/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
fprintf(allfiles.out_err,"TIME for step %d is %f sec\n",sdyn->step,t1-t0);
/*-------------------------------------- check time and number of steps */
if (sdyn->step < sdyn->nstep-1 && sdyn->time <= sdyn->maxtime)
goto timeloop;
/*----------------------------------------------------------------------*/
end:
/*--------------------------------------------------- cleaning up phase */
if (contactflag) amdel(&contactforce_a);
amdel(&intforce_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);
solserv_del_vec(&dispi,1);
solserv_del_vec(&midcon,2);
solserv_del_vec(&vel,1);
solserv_del_vec(&acc,1);
solserv_del_vec(&fie,3);
solserv_del_vec(&work,3);
/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
CCAFREE(actintra);
#endif
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */
#endif /*GEMM*/
/*! @} (documentation module close)*/
