/*!----------------------------------------------------------------------
\file
\brief generalsed alfa time integration algorithm for ssi

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_SSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../output/output_prototypes.h"
#include "ssi_prototypes.h"
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*-----------------------------------------------------------------------*/
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
extern enum _CALC_ACTION calc_action[MAXFIELD];


DOUBLE acttime;

/* due to some problems with the transfer of some pointers to the function
   ssi_struct_decide() the following variables are defined globally    */
   
INT numeq;                           /* number of equations on this proc */
INT numeq_total;                     /* total number of equations */
INT itnum;                           /* counter for NR-Iterations */
INT mod_disp, mod_stress;                                                                                   
INT mod_res_write; 
INT restart; 
INT nstep; 
INT outstep;                         /* counter for output control                          */
INT restartstep;                     /* counter for restart control                         */ 
DOUBLE maxtime; 
DOUBLE t0_res; 
DOUBLE t1_res; 
DOUBLE dt; 
DOUBLE t0;                                        
DOUBLE t1; 
DOUBLE dmax;                         /*infinity norm of residual displacements*/
INT stiff_array;                     /* indice of the active system sparse matrix           */
INT mass_array;                      /* indice of the active system sparse matrix           */
INT damp_array;                      /* indice of the active system sparse matrix           */
INT num_array;                       /* indice of global stiffness matrices                 */
INT actcurve;                        /* indice of active time curve                         */
SOLVAR *actsolv;                     /* pointer to active solution structure                */
PARTITION *actpart;                  /* pointer to active partition                         */
INTRA *actintra;                     /* pointer to active intra-communicator                */
CALC_ACTION *action;                 /* pointer to the structure cal_action enum            */
DIST_VECTOR *vel;                    /* total velocities                                    */
DIST_VECTOR *acc;                    /* total accelerations                                 */
DIST_VECTOR *fie;                    /* internal forces and working array                   */
DIST_VECTOR *dispi;                  /* distributed vector to hold incremental displacments */
DIST_VECTOR *work;                   /* working vectors                                     */
DIST_VECTOR *coup_force;
ARRAY intforce_a;                    /* redundant vector of full length for internal forces */
DOUBLE *intforce;                                                                             
ARRAY dirich_a;                      /* red. vector of full length for dirich-part of rhs   */
DOUBLE *dirich;                      
INT length_of_coup_force;
STRUCT_DYN_CALC dynvar;              /* variables to perform dynamic structural simulation  */
CONTAINER container;                 /* contains variables defined in container.h           */
INT  init;	                     /* flag for solver_control call           */
INT  convergence;                    /* convergence flag                  */
INT  numaf;                          /* number of actual field            */
                                   



/*!---------------------------------------------------------------------                                         
\brief structural control algorithm for ssi problems

<pre>                                                         genk 09/02

This function solves for the structural displacements within an 
multifield problem. The loads are transfared from the fluidfield
as Neumann boundary conditions
			     
</pre>   

\param *ssidyn   SSI_DYNAMIC	(i)				
\param *sdyn	 STRUCT_DYNAMIC	(i)				
\param *actfield FIELD          (i)     actual field		
\param  mctrl    INT            (i)     evaluation flag		
\param  numfs    INT            (i)     number of actual field	
\param  ssiitnum INT            (i)     counter for Iterations over fields
\return void 

------------------------------------------------------------------------*/
void ssi_struct(   
                   SSI_DYNAMIC       *ssidyn,
                   STRUCT_DYNAMIC    *sdyn, 
		   FIELD             *actfield, 
		   INT                mctrl, 
		   INT                ssiitnum,
		   enum _SSI_COUPTYP  ssi_couptyp
	       )
{
INT                  i;		         /* simply a counter		                        */

DOUBLE        dirichfacs[10];	 /* factors needed for dirichlet-part of rhs            */
/*STRUCT_DYN_CALC dynvar;        */   /* variables to perform dynamic structural simulation  */              
/*CONTAINER       container;     */   /* contains variables defined in container.h           */

#ifdef DEBUG 
dstrc_enter("ssi_struct");
#endif


switch (ssi_couptyp)
{
case ssi_master: numaf=0; break;
case ssi_slave:  numaf=1; break;
default: dserror("SSI COUPTYP out of range!\n");
}

if (mctrl > 1) 
{
  ssi_struct_decide(mctrl, numaf);
}

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1: 
sdyn->dt=ssidyn->dt;
sdyn->maxtime=ssidyn->maxtime;
sdyn->nstep=ssidyn->nstep;
container.isdyn   = 1;    
container.actndis = 0;    
container.coupl_typ = ssidyn->conformmesh;
outstep=0;
restartstep=0;

/*----------------------------------------------------------------------*/
/*   ALLOCATE A DIRICH CONDITION FOR EACH SLAVE GNODE AT THE INTERFACE  */
/*----------------------------------------------------------------------*/
/* check if actfield == slave field */
if(numaf == 1)
{
  ssi_alloc_dirich(actfield);
}
/*----------------------------------------------------------------------*/
restart = genprob.restart;
/*--------------------------------------------------- set some pointers */
actsolv            = &(solv[numaf]);
actpart            = &(partition[numaf]);
action             = &(calc_action[numaf]);
container.fieldtyp = actfield->fieldtyp;
/*----------------------------------------- check for explicit dynamics */
dsassert(sdyn->Typ==gen_alfa,"Structural DYNAMICTYP not possible for SSI\n");

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
if (actintra->intra_fieldtyp != structure) break;

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
solserv_getmatdims(&(actsolv->sysarray[stiff_array]),
                   actsolv->sysarray_typ[stiff_array],
                   &numeq,
                   &numeq_total);

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
if (par.myrank==i)
printf("PROC  %3d | FIELD STRUCTURE | number of equations      : %10d \n", 
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD STRUCTURE | total number of equations: %10d \n",numeq_total);
if (par.myrank==0) printf("\n\n");

/*---------------------------------------allocate 5 dist. vectors 'rhs' */
/*  these hold ssi-load vector original load vector, load vector 
    at time t and t-dt and interpolated load vector                     */
actsolv->nrhs = 5;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
for (i=0; i<actsolv->nrhs; i++) solserv_zero_vec(&(actsolv->rhs[i]));

/*--initialize sol.a.da[0..7][numdf-1] at every node, set value to zero */
/*---------------------------------------only necessary for slave field */
#if 0
if (numaf == 1)
  init_sol_a_tozero(actfield);
/*-initialize sol_mf.a.da[0..6][numdf] at every node, set value to zero */
/*--------------------------------------only necessary for master field */
if (numaf == 0)
  init_sol_mf_a_tozero(actfield);
#endif
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

/*------------------------------------ allocate two vectors coup_force */
if (numaf == 1)
{
  length_of_coup_force = actfield->dis->numdf - numeq_total;  
  solserv_create_vec(&coup_force,2,numeq_total,length_of_coup_force,"DV");
  for (i=0; i<2; i++) solserv_zero_vec(&(coup_force[i]));
}

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
/*-------------------------------------- create the original rhs vector */
/*-------------------------- the approbiate action is set inside calrhs */
/*---------------------- this vector holds loads due to external forces */
container.kstep = 0;
container.inherit = 1;
container.point_neum = 1;
*action = calc_struct_eleload;
calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[2]),action,&container);

/*------------------------------------------------- copy the rhs vector */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[3]));

/* put a zero the the place 7 in node->sol to init the velocities and accels */
/* of prescribed displacements */
solserv_putdirich_to_dof(actfield,0,0,0.0,8);

/*---------------------------------- get factor at a certain time t=0.0 */
dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));

/*-------------------------------------- multiply load vector by rldfac */
solserv_scalarprod_vec(&(actsolv->rhs[2]),dynvar.rldfac);

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,0);


/*----- also put prescribed displacements to the nodes in field sol at  */
/*                                  place 3 separate from the free dofs */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,3);

/*-------------------------------------------- make norm of initial rhs */
solserv_vecnorm_euclid(actintra,&(actsolv->rhs[2]),&(dynvar.rnorm));

/*---------------------------------------------- compute initial energy */
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);

/*----------------------------------------- output to GID postprozessor */
if (ioflags.struct_disp_gid==1 || ioflags.struct_stress_gid==1)
if (par.myrank==0) 
{
   out_gid_domains_ssi(actfield, numaf);
}

/*--------------------------------------------------- check for restart */
if (restart)
{
/*   dserror("restart not implemented yet for fsi-problems!\n"); */
   t0_res = ds_cputime();
   /*-------------- save the stepsize as it will be overwritten in sdyn */
   dt    = sdyn->dt;
   /*------ save the number of steps, as it will be overwritten in sdyn */
   nstep = sdyn->nstep;
   maxtime = sdyn->maxtime;
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
   /*------------------------------------------- switch the restart off */
   restart=0;
   /*----------------------------------------------------- measure time */
   t1_res = ds_cputime();
   fprintf(allfiles.out_err,"TIME for restart reading is %f sec\n",t1_res-t0_res);
}
/*------------------------------------------------------- printout head */
/* if (par.myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);*/

/* ssi_struct_decide is defined in ssi_service.c */

  ssi_struct_decide(mctrl, numaf);


break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history structural field:                             * 
 * sol[0][j]           ... total displacements at time (t)              * 
 * sol[1][j]           ... velocities at time (t)		        *
 * sol[2][j]           ... accels at time (t)         		        * 
 * sol[3][j]           ... prescribed displacements at time (t-dt)      * 
 * sol[4][j]           ... prescribed displacements at time (t)	        *
 * sol[5][j]           ... place 4 - place 3                	        *
 * sol[6][j]           ... the  velocities of prescribed dofs  	        *
 * sol[7][j]           ... the  accels of prescribed dofs  	        *
 * sol[8][j]           ... working space   	                        * 
 * sol[9][j]           ... total displacements at time (t-dt)  	        *
 * sol[10][j]          ... velocities at time (t-dt)        	        * 
 *======================================================================*/

/*
   rhs[4]    load vector due to ssi loads at time t
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

   rhs[4]       / ssi-loads             / ssi-loads             / ssi-loads     / ssi-loads
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
case 2:
/*------------------------------------------------- write memory report */
/* if (par.myrank==0) dsmemreport();                                    */
/*------------------------------------------------ output to the screen */
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
if (actintra->intra_fieldtyp != structure) break;
if (par.myrank==0) 
{   
  if (numaf == 0)
  {
    printf("Solving MASTERFIELD ...\n"); 
    printf("---------------------------------------------------------------- \n");
  }
  if (numaf == 1)
  {
    printf("Solving SLAVEFIELD ...\n"); 
    printf("---------------------------------------------------------------- \n");
  }
}
/* copy actual relaxation parameter omega into the container */
container.relax_param = ssidyn->relax;
/*------------------------ copy solution from sol[9][j] to sol[0][j] ---*/
if (ssidyn->ifsi>=4 && ssiitnum>0) solserv_sol_copy(actfield,0,0,0,9,0);

/*----------------------------- copy from nodal sol[1][j] to sol[10][j] */
if (ssidyn->ifsi>=4 && ssiitnum==0) solserv_sol_copy(actfield,0,0,0,1,10);

/*--- put time to global variable for time-dependent load distributions */
acttime = sdyn->time;
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);

/*---------------------- set incremental displacements dispi[0] to zero */
solserv_zero_vec(&dispi[0]);

/*------------------------- set residual displacements in nodes to zero */
solserv_result_resid(actfield,actintra,&dispi[0],0,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
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


/*-------------- assemble external forces due to ssi coupling to rhs[4] */ 
/*if (ssiitnum == 0)*/
solserv_zero_vec(&(actsolv->rhs[4]));
if (numaf == 0)
{
  /* copy old coupling forces sol_mf.a.da[5][..] to sol_mf.a.da[4][..] */
  /* only for non-conforming discretization */

/*  if((ssidyn->conformmesh == 1)&&(sdyn->step > 1))*/
    solserv_sol_add(actfield,0,3,3,5,4,1.0);

  /*--- set factors needed for prescribed displacement terms on rhs eff */
  /*
  dirichfacs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))         
  dirichfacs[1] =  (1.0-alpham)*(1.0/beta)/dt                 
  dirichfacs[2] =  (1.0-alpham)/(2*beta) - 1                  
  dirichfacs[3] = -(1.0-alphaf)*(gamma/beta)/dt               
  dirichfacs[4] =  (1.0-alphaf)*gamma/beta - 1                
  dirichfacs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)*dt            
  dirichfacs[6] =  (1.0-alphaf) or 0                          
  dirichfacs[7] =  raleigh damping factor for mass            
  dirichfacs[8] =  raleigh damping factor for stiffness       
  dirichfacs[9] =  dt     
  see phd theses Mok page 165: generalized alfa time integration with 
  prescribed displ.                                    
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

/* ----------------- add old coupling force to the new one, scaled with */
/* --------------------------------------------- -alpha_f / (1-alpha_f) */
  container.dirichfacs    = dirichfacs;
 /* ssiserv_update_coupforc(actfield, &container); */
  container.inherit = 0;
  container.point_neum = 0;
  *action = calc_struct_ssiload;
  calrhs(actfield,actsolv,actpart,actintra,stiff_array,
       &(actsolv->rhs[4]),action,&container);
}

/*------------------------ add up the two parts of the external forces 
                           and store them in rhs[1]                     */

solserv_add_vec(&(actsolv->rhs[4]),&(actsolv->rhs[1]),ONE);			   

/*------- put conv. displ. from master node to corresponding slave node */
/*--- only necessary for slave field (numaf = 1), used for ssi problems */
/* WARNING! With this approach Dirichlet b.c. on coupling nodes are over-
            written, this is not physically consistent !!!              */
/*
dirichfacs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))         
dirichfacs[1] =  (1.0-alpham)*(1.0/beta)/dt                 
dirichfacs[2] =  (1.0-alpham)/(2*beta) - 1                  
dirichfacs[3] = -(1.0-alphaf)*(gamma/beta)/dt               
dirichfacs[4] =  (1.0-alphaf)*gamma/beta - 1                
dirichfacs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)*dt            
dirichfacs[6] = -(1.0-alphaf) or 0                          
dirichfacs[7] =  raleigh damping factor for mass            
dirichfacs[8] =  raleigh damping factor for stiffness       
dirichfacs[9] =  dt     
see phd theses Mok page 165: generalized alfa time integration with prescribed displ.                                    
*/
if (numaf == 1)
{
  dirichfacs[0] = -dynvar.constants[0];
  dirichfacs[1] =  dynvar.constants[1];
  dirichfacs[2] =  dynvar.constants[2];
  dirichfacs[3] = -dynvar.constants[3];
  dirichfacs[4] =  dynvar.constants[4];
  dirichfacs[5] =  dynvar.constants[5];
  dirichfacs[6] = -dynvar.constants[6]; 
  dirichfacs[9] =  sdyn->dt; 
  if (damp_array>0) 
  {
    dirichfacs[7] =  sdyn->m_damp;
    dirichfacs[8] =  sdyn->k_damp;
  }
  else 
  {
    dirichfacs[7] =  0.0;
    dirichfacs[8] =  0.0;
  }
  container.dirichfacs    = dirichfacs;

  ssiserv_put_disp2slave(actfield, &container, &dynvar.rldfac, ssidyn);
}

/*---------------- put the scaled prescribed displacements to the nodes */
/*             in field sol at place 0 together with free displacements */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,0);

/* put the prescribed scaled displacements to the nodes in field sol at */
/*                                  place 4 separate from the free dofs */
solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,4);

/*-------- put presdisplacements(t) - presdisplacements(t-dt) in place 5 */
solserv_adddirich(actfield,0,0,3,4,5,-1.0,1.0);

/*----- set factors needed for prescribed displacement terms on rhs eff */
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
*action = calc_struct_nlnstiffmass;
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);
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

/*------------------------ add rhs from prescribed displacements to rhs */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);

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
                     &(actsolv->sysarray_typ[stiff_array]),0);
/* here put incremental prescribed displacements from sol[5] to sol_increment[0] ? */

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
/* call element routines for calculation of tangential stiffness and intforce */
*action = calc_struct_nlnstiffmass;
container.dvec          = intforce;
container.dirich        = dirich;
container.global_numeq  = numeq_total;
container.dirichfacs    = dirichfacs;
container.kstep         = 0;
calelm(actfield,actsolv,actpart,actintra,stiff_array,mass_array,&container,action);

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

/*------------------ add dirichlet forces from prescribed displacements */
assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
             &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);
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
                     &(actsolv->sysarray_typ[stiff_array]),0);

/*----------------------------------------------- check for convergence */
convergence = 0;
dmax        = 0.0;
solserv_vecnorm_euclid(actintra,&(work[0]),&(dynvar.dinorm));
solserv_vecnorm_euclid(actintra,&(dispi[0]),&(dynvar.dnorm));
solserv_vecnorm_Linf(actintra,&(work[0]),&dmax);
if (dynvar.dinorm < sdyn->toldisp ||
    dynvar.dnorm  < EPS14 ||
    (dynvar.dinorm < EPS14 && dmax < EPS12) )
{
   convergence = 1;
}    
else
{
   itnum++;
   if (itnum==sdyn->maxiter) dserror("No convergence in maxiter steps");
   goto iterloop;
}

/* non-conforming discretization, slave domain */
if ((ssidyn->conformmesh == 1)&&(numaf==1)) 
{
  /* --------- clean up vector sol_mf[4], this vector contains internal */
  /* --------------------------------------------forces for slave nodes */
  solserv_sol_zero(actfield, 0, 3, 4);
}

if (numaf == 1)  /* only for slave domain */
{
  if (ssidyn->conformmesh == 0)
  {
    /* clean up vector of coupling forces, before it is assembled again */
    ssiserv_erase_coupforc(actfield, &container);
  }
  /* ---------------------to decide in calelm and in wall1 what to do */
  *action = calc_struct_ssi_coup_force;
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

   container.dvec          = intforce;
   container.dirich        = dirich;
   container.global_numeq  = numeq_total;
   container.dirichfacs    = dirichfacs;
   container.kstep         = 0;
/* ---------- compute new coupling forces for nodes of coupled elements */
   calelm(actfield, actsolv, actpart, actintra, stiff_array, mass_array, 
          &container, action);

}
/*----------------------------------------------------------------------*/
/*                      END OF EQUILLIBRIUM ITERATION                   */
/*----------------------------------------------------------------------*/

/*----- for iterative staggered schemes save solution of last iteration */
/*- copy old total displacments from nodal sol_mf[0][j] to sol_mf[1][j] */
if (ssidyn->ifsi>=4) solserv_sol_copy(actfield,0,3,3,0,1);

/*------- copy total displacments from nodal sol[0][j] to sol_mf[0][j]  */
solserv_sol_copy(actfield,0,0,3,0,0);

/*------- for sequential staggered schemes sol_mf[0][j] == sol_mf[1][j] */
if (ssidyn->ifsi<4) solserv_sol_copy(actfield,0,3,3,0,1);

/*----------------------------------------------------- print time step */
if (par.myrank==0) 
{
   printf("| NUMITER = %3d                                                | \n",
             itnum+1);
   printf("---------------------------------------------------------------- \n");
   printf("\n"); 
}
if (ssidyn->ifsi>=4)
break;

/*======================================================================* 
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
/*-------------------------------------------------- set some constants */
dyn_setconstants(&dynvar,sdyn,sdyn->dt);
/*----------- make temporary copy of actsolv->rhs[2] to actsolv->rhs[0] */
/*                                   (load at t-dt)                     */
/* because in  dyn_nlnstructupd actsolv->rhs[2] is overwritten but is   */
/* still needed to compute energies in dynnle                           */
solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
/*------------------ update displacements, velocities and accelerations */
dyn_nlnstructupd(actfield,
                 &dynvar,sdyn,actsolv,
                 &(actsolv->sol[0]),   /* total displacements at time t-dt */
                 &(actsolv->sol[1]),   /* total displacements at time t    */
                 &(actsolv->rhs[1]),   /* load vector         at time t    */
                 &(actsolv->rhs[2]),   /* load vector         at time t-dt */
                 &vel[0],              /* velocities          at time t    */
                 &acc[0],              /* accelerations       at time t    */
                 &work[0],             /* working arrays                   */
                 &work[1],             /* working arrays                   */
                 &work[2]);            /* working arrays                   */

/*-------------------------------------- return velocities to sol[1][j] */
solserv_result_total(actfield,actintra, &vel[0],1,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
/*---------------------------------------velocities for prescribed dofs */
solserv_adddirich(actfield,0,0,6,0,1,1.0,0.0);
/*------------------------------------------ return accel. to sol[2][j] */
solserv_result_total(actfield,actintra, &acc[0],2,
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));

/*-------------------------------------------accel. for prescribed dofs */
solserv_adddirich(actfield,0,0,7,0,2,1.0,0.0);
/*------------------------------------------ make all types of energies */
dynnle(&dynvar,sdyn,actintra,actsolv,&dispi[0],&fie[1],&fie[2],
       &(actsolv->rhs[1]),&(actsolv->rhs[0]),&work[0]);
dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);
dynvar.etot = dynvar.epot + dynvar.ekin;

if (numaf == 0) /*masterfield*/
{
/* ------------- copy actual coupling forces sol_mf[4] to vector of old */
/* ------------------------------------------ coupling forces sol_mf[5] */
   /*solserv_sol_add(actfield,0,3,3,4,5,1.0);*/
   solserv_sol_copy(actfield,0,3,3,4,5);
}
/*--- copy the displacement of the master nodes at the coupling surface */
/*--- to the corresponding slave nodes, here we copy the values at time */
/*------- t, up to this routine there are the values resulting from the */
/*--- generalized alpha time integration scheme at the nodes, these are */
/*------------------------------------------------- interpolated values */ 

ssiserv_put_true_disp2slave(actfield, &container);

/*-------------------------------- save actual solution as old solution */
/*------------------------------ copy from nodal sol[0][j] to sol[9][j] */
solserv_sol_copy(actfield,0,0,0,0,9);

/*---------------- save actual relaxed solution as old relaxed solution */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[2][j] */
solserv_sol_copy(actfield,0,3,3,1,2);
/*----------- perform stress calculation  and print out results to .out */
outstep++;
if (outstep==sdyn->updevry_disp)
{   
   outstep=0;
   if (ioflags.struct_stress_file==1)
   {
      *action = calc_struct_stress;
      container.dvec          = NULL;
      container.dirich        = NULL;
      container.global_numeq  = 0;
      container.dirichfacs    = NULL;
      container.kstep         = 0;
      calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);
      /*----------------------- reduce stresses, so they can be written */
      *action = calc_struct_stressreduce;
      container.kstep = 0;
      calreduce(actfield,actpart,actintra,action,&container);
   }
   if (ioflags.struct_stress_file==1 || ioflags.struct_disp_file==1)
      out_sol(actfield,actpart,actintra,sdyn->step,0);
}
/*-------------------------------------- write restart data to pss file */
restartstep++;
if (restartstep==ssidyn->res_write_evry)
{
   restartstep=0;
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
                              &container);  /* contains variables defined in container.h */
}

/*------------------------------------------ measure time for this step */
t1 = ds_cputime();
fprintf(allfiles.out_err,"TIME for step %d is %f sec\n",sdyn->step,t1-t0);
break;

/*======================================================================* 
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:
/*- there are only procs allowed in here, that belong to the structural */
/* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
if (actintra->intra_fieldtyp != structure) break;

/*---------------------------------------- print out the final solution */
if (outstep!=0)
{   
   outstep=0;
   if (ioflags.struct_stress_file==1)
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
   if (ioflags.struct_stress_file==1 || ioflags.struct_disp_file==1)
      out_sol(actfield,actpart,actintra,sdyn->step,0);
}

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
CCAFREE(actintra);
#endif
break;
default:
dserror("Parameter mctrl out of range\n");
}
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dyn_nln_structural */


/*!----------------------------------------------------------------------*
\brief function to set and read variables

<pre>
                                                   firl / genk 10/03    
   
function to set or read variables either for the master field or for
slave field, there exists one special set of variables for each
field, (.)_m and (.)_s corresponds to the master field and to the
slave field respectively,
routine is called in ssi_struct


  The routine is placed in ssi_service.c

\param  mctrl  INT      (i)
        numaf  INT      (i)

</pre>                                                                       
                                                                       
------------------------------------------------------------------------*/
void ssi_struct_decide(INT mctrl, INT numaf) 
{

static INT     numeq_m; 	         /* number of equations on this proc	*/
static INT     numeq_total_m;	         /* total number of equations		*/
static INT     itnum_m; 		 /* counter for NR-Iteration n	*/
static INT     mod_disp_m; 
static INT     mod_stress_m; 
static INT     mod_res_write_m;
static INT     restart_m;
static INT     nstep_m;
static INT     outstep_m;	         /* counter for output control		*/
static INT     restartstep_m;	         /* counter for restart control		*/ 
static DOUBLE  maxtime_m;
static DOUBLE  t0_res_m;
static DOUBLE  t1_res_m;
static DOUBLE  dt_m;
static DOUBLE  t0_m;  
static DOUBLE  t1_m;  
static DOUBLE  dmax_m;  	         /* infinity norm of residual displacements	*/
static INT     stiff_array_m;	         /* indice of the active system sparse matrix */
static INT     mass_array_m;	         /* indice of the active system sparse matrix	*/
static INT     damp_array_m;	         /* indice of the active system sparse matrix	*/
static INT     num_array_m;	         /* indice of global stiffness matrices 	*/
static INT     actcurve_m;	         /* indice of active time curve 		*/

static SOLVAR      *actsolv_m;	         /* pointer to active solution structure	*/
static PARTITION   *actpart_m;	         /* pointer to active partition 		*/
static INTRA       *actintra_m;	         /* pointer to active intra-communicator	*/
static CALC_ACTION *action_m;    	 /* pointer to the structure cal_action enum	*/

static DIST_VECTOR  *vel_m;		 /* total velocities				*/	   
static DIST_VECTOR  *acc_m;		 /* total accelerations 			*/
static DIST_VECTOR  *fie_m;		 /* internal forces and working array		*/
static DIST_VECTOR  *dispi_m;		 /* distributed vector to hold incremental displacments */ 
static DIST_VECTOR  *work_m;		 /* working vectors				*/
static ARRAY        intforce_a_m;	 /* redundant vector of full length for internal forces */
static DOUBLE       *intforce_m;
static ARRAY        dirich_a_m;	         /* red. vector of full length for dirich-part of rhs */
static DOUBLE       *dirich_m;
static STRUCT_DYN_CALC dynvar_m;	 /* variables to perform dynamic structural simulation */	 
static CONTAINER       container_m;	 /* contains variables defined in container.h	*/

static INT     numeq_s; 	         /* number of equations on this proc	*/
static INT     numeq_total_s;	         /* total number of equations		*/
static INT     itnum_s; 		 /* counter for NR-Iteration n	*/
static INT     mod_disp_s; 
static INT     mod_stress_s; 
static INT     mod_res_write_s;
static INT     restart_s;
static INT     nstep_s;
static INT     outstep_s;	         /* counter for output control		*/
static INT     restartstep_s;	         /* counter for restart control		*/ 
static DOUBLE  maxtime_s;
static DOUBLE  t0_res_s;
static DOUBLE  t1_res_s;
static DOUBLE  dt_s;
static DOUBLE  t0_s;  
static DOUBLE  t1_s;  
static DOUBLE  dmax_s;  	         /* infinity norm of residual displacements	*/
static INT     stiff_array_s;	         /* indice of the active system sparse matrix */
static INT     mass_array_s;	         /* indice of the active system sparse matrix	*/
static INT     damp_array_s;	         /* indice of the active system sparse matrix	*/
static INT     num_array_s;	         /* indice of global stiffness matrices 	*/
static INT     actcurve_s;	         /* indice of active time curve 		*/

static SOLVAR       *actsolv_s;	         /* pointer to active solution structure	*/
static PARTITION    *actpart_s;	         /* pointer to active partition 		*/
static INTRA        *actintra_s;          /* pointer to active intra-communicator	*/
static CALC_ACTION  *action_s;  	 /* pointer to the structure cal_action enum	*/
static DIST_VECTOR  *vel_s;		 /* total velocities				*/	   
static DIST_VECTOR  *acc_s;		 /* total accelerations 			*/
static DIST_VECTOR  *fie_s;		 /* internal forces and working array		*/
static DIST_VECTOR  *dispi_s;		 /* distributed vector to hold incremental displacments */ 
static DIST_VECTOR  *work_s;		 /* working vectors				*/
static ARRAY        intforce_a_s;        /* redundant vector of full length for internal forces */
static DOUBLE       *intforce_s;
static ARRAY        dirich_a_s;	         /* red. vector of full length for dirich-part of rhs */
static DOUBLE       *dirich_s;
static STRUCT_DYN_CALC dynvar_s;	 /* variables to perform dynamic structural simulation */	 
static CONTAINER       container_s;	 /* contains variables defined in container.h	*/

#ifdef DEBUG 
dstrc_enter("ssi_struct_decide");
#endif

if (mctrl==1) 
{
  /* for mctrl = 1 : definition of the variables */
  if (numaf == 0) 
  {
    /* numaf == 0 --> master field */
    /* copy data from variables defined in ssi_struct to the static variables in ssi_struct_decide */ 
    numeq_m		     =    numeq; 	       
    numeq_total_m	     =    numeq_total;	       
    itnum_m		     =    itnum; 	       
    mod_disp_m               =    mod_disp;  
    mod_stress_m             =    mod_stress;  
    mod_res_write_m	     =    mod_res_write;        
    restart_m  	             =    restart;	       
    nstep_m		     =    nstep; 	       
    outstep_m  	             =    outstep;	       
    restartstep_m	     =    restartstep;	       
    maxtime_m  	             =    maxtime;	       
    t0_res_m                 =    t0_res;
    t1_res_m	             =    t1_res;        
    dt_m		     =    dt;		       
    t0_m		     =    t0; 	       
    t1_m		     =    t1; 	       
    dmax_m		     =    dmax;  	       
    stiff_array_m	     =    stiff_array;	       
    mass_array_m	     =    mass_array;	       
    damp_array_m	     =    damp_array;	       
    num_array_m	             =    num_array;	       
    actcurve_m 	             =    actcurve;	       

    actsolv_m 	             =    actsolv;	       
    actpart_m 	             =    actpart;	       
    actintra_m	             =    actintra;	       
    action_m  	             =    action;	       
    intforce_m	             =    intforce;	       
    dirich_m  	             =    dirich;	       
    vel_m		     =    vel; 	       
    acc_m		     =    acc;  	       
    fie_m		     =    fie;  	       
    dispi_m		     =    dispi;	       
    work_m		     =    work; 	       

    intforce_a_m	     =    intforce_a;	       
    dirich_a_m 	             =    dirich_a;	       
    dynvar_m                 =    dynvar;      
    container_m	             =    container;	       
  }

  if (numaf == 1)
  {   
    /* numaf == 1 --> slave field */ 
    /* copy data from variables defined in ssi_struct to the static variables in ssi_struct_decide */ 
    numeq_s		     =    numeq;	 	
    numeq_total_s	     =    numeq_total;  	
    itnum_s		     =    itnum;	 	
    mod_disp_s               =    mod_disp;  
    mod_stress_s             =    mod_stress;  
    mod_res_write_s	     =    mod_res_write;	 
    restart_s 	             =    restart;	 	
    nstep_s		     =    nstep;	 	
    outstep_s  	             =    outstep;	 	
    restartstep_s	     =    restartstep;  	
    maxtime_s  	             =    maxtime;	 	
    t0_res_s	             =    t0_res;	  
    t1_res_s	             =    t1_res;	  
    dt_s		     =    dt;  	 	
    t0_s                     =    t0;         
    t1_s		     =    t1;         
    dmax_s		     =    dmax;	 	
    stiff_array_s	     =    stiff_array;  	
    mass_array_s	     =    mass_array;   	
    damp_array_s	     =    damp_array;   	
    num_array_s	             =    num_array;	 	
    actcurve_s 	             =    actcurve;	 	

    actsolv_s 	             =    actsolv;	 	
    actpart_s 	             =    actpart;	 	
    actintra_s	             =    actintra;	 	
    action_s  	             =    action;	 	
    vel_s		     =    vel;	 	
    acc_s		     =    acc;	 	
    fie_s		     =    fie;           	
    dispi_s		     =    dispi;  	 	
    work_s		     =    work;  	 	
    intforce_s	             =    intforce;	 	
    dirich_s  	             =    dirich; 	 	

    intforce_a_s	     =    intforce_a;	 	
    dirich_a_s 	             =    dirich_a;	 	
    dynvar_s                 =    dynvar;      
    container_s 	     =    container;	 	

  }
}  

if (mctrl > 1) 
{
  /* for mctrl > 1 desicion between master field and slave field */
  
  if (numaf == 0) 
  {
    /* numaf == 0 --> master field */ 
    /* for ssi_couptyp = ssi_master set the master variables (static) to the nonstatic variables of function ssi_struct.c */
    numeq    		     =    numeq_m;	    
    numeq_total     	     =    numeq_total_m;     
    itnum		     =    itnum_m;	      
    mod_disp                 =    mod_disp_m; 
    mod_stress               =    mod_stress_m; 
    mod_res_write	     =    mod_res_write_m;
    restart		     =    restart_m;
    nstep		     =    nstep_m;
    outstep	      	     =    outstep_m;	    
    restartstep     	     =    restartstep_m;     
    maxtime		     =    maxtime_m;
    t0_res	             =    t0_res_m;
    t1_res	             =    t1_res_m;
    dt 		             =    dt_m;
    t0  		     =    t0_m;  
    t1  		     =    t1_m;  
    dmax        	     =    dmax_m;      
    stiff_array      	     =    stiff_array_m;      
    mass_array  	     =    mass_array_m;        
    damp_array       	     =    damp_array_m;     
    num_array   	     =    num_array_m;  	
    actcurve		     =    actcurve_m;		

    actsolv		     =	  actsolv_m;	      
    actpart		     =	  actpart_m;
    actintra		     =	  actintra_m;	      
    action		     =	  action_m; 
    vel 		     =	  vel_m;     
    acc 		     =	  acc_m;     
    fie 		     =	  fie_m;     
    dispi		     =	  dispi_m;  
    work		     =	  work_m;   
    intforce		     =	  intforce_m;
    dirich	             =	  dirich_m;

    intforce_a       	     =    intforce_a_m;     
    dirich_a		     =    dirich_a_m;		
    dynvar                   =    dynvar_m;   
    container      	     =    container_m;    

  }
  
  if (numaf == 1) 
  {
    /* numaf == 1 --> slave field */ 
    /* for ssi_couptyp = ssi_slave set the slave variables (static) to the nonstatic variables of function ssi_struct.c */
    numeq	  	     =    numeq_s;	    
    numeq_total  	     =    numeq_total_s;     
    itnum	  	     =    itnum_s;	      
    mod_disp	  	     =    mod_disp_s; 
    mod_stress   	     =    mod_stress_s; 
    mod_res_write	     =    mod_res_write_s;
    restart	  	     =    restart_s;
    nstep	  	     =    nstep_s;
    outstep	  	     =    outstep_s;	    
    restartstep  	     =    restartstep_s;     
    maxtime	  	     =    maxtime_s;
    t0_res	  	     =    t0_res_s;
    t1_res	  	     =    t1_res_s;
    dt 	  	             =    dt_s;
    t0 	  	             =    t0_s;  
    t1 	  	             =    t1_s;  
    dmax	  	     =    dmax_s;      
    stiff_array  	     =    stiff_array_s;      
    mass_array   	     =    mass_array_s;        
    damp_array   	     =    damp_array_s;     
    num_array    	     =    num_array_s;  	
    actcurve	  	     =    actcurve_s;		

    actsolv		     =	  actsolv_s;	      
    actpart		     =	  actpart_s;
    actintra		     =	  actintra_s;	      
    action		     =	  action_s; 
    vel 		     =	  vel_s;     
    acc 		     =	  acc_s;     
    fie 		     =	  fie_s;     
    dispi		     =	  dispi_s;  
    work		     =	  work_s;   
    intforce		     =	  intforce_s;
    dirich		     =	  dirich_s;

    intforce_a       	     =    intforce_a_s;     
    dirich_a		     =    dirich_a_s;		
    dynvar                   =    dynvar_s;   
    container      	     =    container_s;    

  }
}
#ifdef DEBUG 
dstrc_exit();
#endif
}

#endif

/*! @} (documentation module close)*/
