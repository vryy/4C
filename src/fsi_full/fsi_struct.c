/*!----------------------------------------------------------------------
\file
\brief generalsed alfa time integration algorithm for fsi

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


#ifdef D_FSI


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
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



/*!---------------------------------------------------------------------
\brief structural control algorithm for fsi problems

<pre>                                                         genk 09/02

This function solves for the structural displacements within an
multifield problem. The loads are transfared from the fluidfield
as Neumann boundary conditions

</pre>

\param *fsidyn   FSI_DYNAMIC    (i)
\param *sdyn     STRUCT_DYNAMIC (i)
\param *actfield FIELD          (i)     actual field
\param  mctrl    INT            (i)     evaluation flag
\param  fsiitnum INT            (i)     counter for Iterations over fields
\return void

------------------------------------------------------------------------*/
void fsi_struct(
    FIELD              *actfield,
    INT                 disnum_calc,
    INT                 disnum_io,
    INT                 mctrl,
    INT                 fsiitnum
    )

{

  INT                  i;              /* simply a counter                  */
  static INT           numsf;          /* actual number of struct field     */
  static INT           numeq;          /* number of equations on this proc  */
  static INT           numeq_total;    /* total number of equations         */
  INT                  init;           /* flag for solver_control call      */
  static INT           itnum;          /* counter for NR-Iterations         */
  INT                  convergence;    /* convergence flag                  */
  static INT           restart;
  static INT           nstep;
  static INT           outstep;        /* counter for output control        */
  static INT           restartstep;    /* counter for restart control       */
  static DOUBLE        maxtime;
  static DOUBLE        t0_res,t1_res;
  static DOUBLE        dt;
  static DOUBLE        t0,t1;
  static DOUBLE        dmax;           /* infinity norm of residual displacements    */

  static INT           stiff_array;    /* indice of the active system sparse matrix  */
  static INT           mass_array;     /* indice of the active system sparse matrix  */
  static INT           damp_array;     /* indice of the active system sparse matrix  */
  static INT           actcurve=0;     /* indice of active time curve                */

  static SOLVAR       *actsolv;        /* pointer to active solution structure       */
  static PARTITION    *actpart;        /* pointer to active partition                */
  static INTRA        *actintra;       /* pointer to active intra-communicator       */
  static CALC_ACTION  *action;         /* pointer to the structure cal_action enum   */

  static DIST_VECTOR  *vel;            /* total velocities                           */
  static DIST_VECTOR  *acc;            /* total accelerations                        */
  static DIST_VECTOR  *fie;            /* internal forces and working array          */
  static DIST_VECTOR  *dispi;          /* distributed vector to hold incr. disp.     */
  static DIST_VECTOR  *work;           /* working vectors                            */

  static ARRAY         intforce_a;     /* red. vector of full length for int forces  */
  static DOUBLE       *intforce;
  static ARRAY         fsiforce_a;     /* red. vector of full length for int forces  */
  static DOUBLE       *fsiforce;
  static ARRAY         dirich_a;       /* red. vector of full length for dirich rhs  */
  static DOUBLE       *dirich;
  static DOUBLE        dirichfacs[10]; /* factors needed for dirichlet-part of rhs   */
  static STRUCT_DYN_CALC dynvar;       /* variables to perform dynamic struct sim    */
  static CONTAINER       container;    /* contains variables defined in container.h  */

  static FSI_DYNAMIC       *fsidyn;
  static STRUCT_DYNAMIC    *sdyn;
  FILE           *out = allfiles.out_out;

#ifdef BINIO
  static BIN_OUT_FIELD out_context;
  static BIN_OUT_FIELD restart_context;
#endif




#ifdef DEBUG
  dstrc_enter("fsi_struct");
#endif




  switch (mctrl)
  {
    /*======================================================================*
      |                      I N I T I A L I S A T I O N                     |
     *======================================================================*/
    case 1:
      numsf             = genprob.numsf;
      fsidyn            = alldyn[genprob.numaf+1].fsidyn;
      sdyn              = alldyn[numsf].sdyn;

      sdyn->dt          = fsidyn->dt;
      sdyn->maxtime     = fsidyn->maxtime;
      sdyn->nstep       = fsidyn->nstep;

      container.isdyn   = 1;
      container.disnum  = disnum_calc;

      outstep           = 0;
      restartstep       = 0;
      restart           = genprob.restart;


      /* set some pointers */
      actsolv            = &(solv[numsf]);
      actpart            = &(partition[numsf]);
      action             = &(calc_action[numsf]);
      container.fieldtyp = actfield->fieldtyp;


      /* check for explicit dynamics */
      dsassert(sdyn->Typ==gen_alfa,"Structural DYNAMICTYP not possible for FSI\n");



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



      /* init the variables in dynvar to zero */
      dynvar.rldfac = 0.0;
      dynvar.rnorm  = 0.0;
      dynvar.epot   = 0.0;
      dynvar.eout   = 0.0;
      dynvar.etot   = 0.0;
      dynvar.ekin   = 0.0;
      dynvar.dinorm = 0.0;
      dynvar.dnorm  = 0.0;
      for (i=0; i<20; i++) dynvar.constants[i] = 0.0;
      acttime       = 0.0;



      /* check presence of damping matrix *
       *   and set indice of stiffness and mass sparse matrices */

#ifdef SUBDIV
      if (actfield->subdivide > 0)
      {

        stiff_array = 1;
        mass_array  = 2;

        if (sdyn->damp==1)
        {
          damp_array  = 3;
          actsolv->nsysarray=4;
        }
        else
        {
          damp_array  =-1;
          actsolv->nsysarray=3;
        }

      }
      else
#endif
      {

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

      }





      /* stiff_array already exists, so copy the mask of it to
       * mass_array (and damp_array if needed)
       * reallocate the vector of sparse matrices and the vector of their types
       * formerly lenght 1, now lenght 2 or 3 dependent on presence of damp_array */
      actsolv->sysarray_typ =
        (SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
      if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

      actsolv->sysarray =
        (SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
      if (!actsolv->sysarray_typ) dserror("Allocation of memory failed");

      /* copy the matrices sparsity mask from stiff_array to mass_array (and to damp_array) */
      solserv_alloc_cp_sparsemask(
          actintra,
          &(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[mass_array]),
          &(actsolv->sysarray[mass_array]));

      if (damp_array>0)
        solserv_alloc_cp_sparsemask(
            actintra,
            &(actsolv->sysarray_typ[stiff_array]),
            &(actsolv->sysarray[stiff_array]),
            &(actsolv->sysarray_typ[damp_array]),
            &(actsolv->sysarray[damp_array]));



      /* init the dist sparse matrices to zero */
      for (i=0; i<actsolv->nsysarray; i++)
        solserv_zero_mat(
            actintra,
            &(actsolv->sysarray[i]),
            &(actsolv->sysarray_typ[i])
            );


      /* get global and local number of equations */
      solserv_getmatdims(&(actsolv->sysarray[stiff_array]),
          actsolv->sysarray_typ[stiff_array],
          &numeq,
          &numeq_total);


      /* output to the screen */
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



      /* allocate 5 dist. vectors 'rhs':
       *   these hold fsi-load vector original load vector, load vector
       *   at time t and t-dt and interpolated load vector */
      actsolv->nrhs = 5;
      solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
      for (i=0; i<actsolv->nrhs; i++)
        solserv_zero_vec(&(actsolv->rhs[i]));


      /* there are 2 solution vector to hold total displ.
       * one at time t and one at time t-dt */
      actsolv->nsol= 2;
      solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
      for (i=0; i<actsolv->nsol; i++)
        solserv_zero_vec(&(actsolv->sol[i]));


      /* there is one vector to hold incremental displacements */
      solserv_create_vec(&dispi,1,numeq_total,numeq,"DV");
      solserv_zero_vec(&(dispi[0]));


      /* allocate one vector vel */
      solserv_create_vec(&vel,1,numeq_total,numeq,"DV");
      solserv_zero_vec(&(vel[0]));


      /* allocate one vector acc */
      solserv_create_vec(&acc,1,numeq_total,numeq,"DV");
      solserv_zero_vec(&(acc[0]));


      /* allocate one redundant vector intforce of full lenght
       * this is used by the element routines to assemble the internal forces*/
      intforce = amdef("intforce",&intforce_a,numeq_total,1,"DV");


      /* create a vector of full length for dirichlet part of rhs */
      dirich = amdef("dirich",&dirich_a,numeq_total,1,"DV");


      /* this is used to assemble the fsi-coupling forces */
      fsiforce = amdef("fsiforce",&fsiforce_a,numeq_total,1,"DV");


      /* allocate 3 DIST_VECTOR fie
       * to hold internal forces at t, t-dt and inbetween */
      solserv_create_vec(&fie,3,numeq_total,numeq,"DV");
      for (i=0; i<3; i++)
        solserv_zero_vec(&(fie[i]));


      /* allocate three working vectors
       *   By optimizing this routine one could live with one or two working
       *   vectors, I needed three to make things straight-forward and easy */
      solserv_create_vec(&work,3,numeq_total,numeq,"DV");
      for (i=0; i<3; i++)
        solserv_zero_vec(&(work[i]));




      /* initialize solver on all matrices:
       *-----------------------------------
       *  NOTE: solver init phase has to be called with each matrix one wants to
       *  solve with. Solver init phase has to be called with all matrices
       *  one wants to do matrix-vector products and matrix scalar products.
       *  This is not needed by all solver libraries, but the solver-init phase
       *  is cheap in computation (can be costly in memory)
       *  There will be no solver call on mass or damping array.
       */

      /* initialize solver */
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


      /* init the assembly for stiffness and for mass matrix
       * (damping is not assembled) */
      init_assembly(actpart,actsolv,actintra,actfield,stiff_array,disnum_calc);
      init_assembly(actpart,actsolv,actintra,actfield,mass_array,disnum_calc);


      /* init the element calculating routines */
      *action = calc_struct_init;
      calinit(actfield,actpart,action,&container);


      /* call elements to calculate stiffness and mass */
      *action = calc_struct_nlnstiffmass;
      container.dvec          = NULL;
      container.dirich        = NULL;
      container.global_numeq  = 0;
      container.dirichfacs    = NULL;
      container.kstep         = 0;
      calelm(actfield, actsolv, actpart, actintra, stiff_array, mass_array,
          &container,action);



      /* calculate damping matrix */
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



      /* create the original rhs vector:
       *   the approbiate action is set inside calrhs
       *   this vector holds loads due to external forces */
      container.kstep = 0;
      container.inherit = 1;
      container.point_neum = 1;
      *action = calc_struct_eleload;
      calrhs(actfield,actsolv,actpart,actintra,stiff_array,
          &(actsolv->rhs[2]),action,&container);



      /* copy the rhs vector */
      solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[3]));


      /* put a zero the the place 7 in node->sol to init the velocities and accels
       * of prescribed displacements */
      solserv_sol_zero(actfield, disnum_calc, node_array_sol, 8);


      /* get factor at a certain time t=0.0 */
      dyn_facfromcurve(actcurve,0.0,&(dynvar.rldfac));


      /* multiply load vector by rldfac */
      solserv_scalarprod_vec(&(actsolv->rhs[2]),dynvar.rldfac);


      /* put the scaled prescribed displacements to the nodes
       * in field sol at place 0 together with free displacements */
      solserv_putdirich_to_dof(actfield,disnum_calc,0,0,sdyn->time);


      /* also put prescribed displacements to the nodes in field sol at
       * place 3 separate from the free dofs */

      /*solserv_putdirich_to_dof(actfield,0,0,dynvar.rldfac,3);*/
      solserv_putdirich_to_dof(actfield,disnum_calc,0,3,sdyn->time);


      /* make norm of initial rhs */
      solserv_vecnorm_euclid(actintra,&(actsolv->rhs[2]),&(dynvar.rnorm));


      /* compute initial energy */
      dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);


      /* initialise FSI-predictor */
      if (fsidyn->ifsi == 2  ||  fsidyn->ifsi >= 4)
        fsi_structpredictor(actfield,disnum_calc,1);


      /* initialise energy check */
      if (fsidyn->ichecke > 0)
        fsi_dyneint(actfield,disnum_calc,1);


      if (fsidyn->ifsi >= 4) {
        /* entry 10 is used below. So make sure it exists now. */
        solserv_sol_zero(actfield, disnum_calc, node_array_sol, 10);
      }
      else {
        /* Make sure the entry for our old solution is there. */
        solserv_sol_zero(actfield, disnum_calc, node_array_sol, 9);
      }
      if (fsidyn->ifsi == 6) {
        solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, 6);
      }

#ifdef BINIO

      /* initialize binary output
       * It's important to do this only after all the node arrays are set
       * up because their sizes are used to allocate internal memory. */
      init_bin_out_field(&out_context,
          &(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),
          actfield, actpart, actintra, disnum_io);

      if(disnum_io != disnum_calc)
        init_bin_out_field(&restart_context,
            &(actsolv->sysarray_typ[stiff_array]),
            &(actsolv->sysarray[stiff_array]),
            actfield, actpart, actintra, disnum_calc);
#endif



#ifdef PARALLEL
      /* output to GID postprozessor */
      if (ioflags.output_gid==1 && par.myrank==0)
      {
        out_gid_domains(actfield, disnum_io);
      }
#endif



      /* check for restart */
      if (restart)
      {
        t0_res = ds_cputime();

        /* save the stepsize as it will be overwritten in sdyn */
        dt    = sdyn->dt;

        /* save the number of steps, as it will be overwritten in sdyn */
        nstep = sdyn->nstep;
        maxtime = sdyn->maxtime;

        /* the step to read in is restart */
#if defined(BINIO)
        restart_read_bin_nlnstructdyn(sdyn,
            &dynvar,
            &(actsolv->sysarray_typ[stiff_array]),
            &(actsolv->sysarray[stiff_array]),
            actfield,
            actpart,
            disnum_calc,
            actintra,
            actsolv->nrhs, actsolv->rhs,
            actsolv->nsol, actsolv->sol,
            1            , dispi       ,
            1            , vel         ,
            1            , acc         ,
            3            , fie         ,
            3            , work        ,
            restart);
#else
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
#endif


        /* put the dt to the structure */
        sdyn->dt = dt;

        /* put nstep to the structure */
        sdyn->nstep = nstep;
        sdyn->maxtime = maxtime;

        /* switch the restart off */
        restart=0;

        /* measure time */
        t1_res = ds_cputime();
        fprintf(allfiles.out_err,"TIME for restart reading is %f sec\n",t1_res-t0_res);

        /* initialise predictor */
        if (fsidyn->ifsi == 2  ||  fsidyn->ifsi >= 4)
          fsi_structpredictor(actfield,disnum_calc,2);

      }  /* if (restart) */



      /* printout head */
      /* if (par.myrank==0) dyn_nlnstruct_outhead(&dynvar,sdyn);*/

      break;












      /*======================================================================*
       *                      S O L U T I O N    P H A S E                    *
       *======================================================================*
       * nodal solution history structural field:                             *
       * sol[0][j]           ... total displacements at time (t)              *
       * sol[1][j]           ... velocities at time (t)                       *
       * sol[2][j]           ... accels at time (t)                           *
       * sol[3][j]           ... prescribed displacements at time (t-dt)      *
       * sol[4][j]           ... prescribed displacements at time (t)         *
       * sol[5][j]           ... place 4 - place 3                            *
       * sol[6][j]           ... the  velocities of prescribed dofs           *
       * sol[7][j]           ... the  accels of prescribed dofs               *
       * sol[8][j]           ... working space                                *
       * sol[9][j]           ... total displacements at time (t-dt)           *
       * sol[10][j]          ... velocities at time (t-dt)                    *
       * sol_mf[0][j]        ... latest struct-displacements                  *
       * sol_mf[1][j]        ... (relaxed) displ. of the last iteration step  *
       * sol_mf[2][j]        ... converged relaxed displ. at time (t-dt)      *
       * sol_mf[3][j]        ... actual dispi                                 *
       * sol_mf[4][j]        ... FSI coupl.-forces at the end of the timestep *
       * sol_mf[5][j]        ... FSI coupl.-forces at beginning of the timest.*
       * sol_mf[6][j]        ... used in fsi_gradient.c                       *
       *======================================================================*/

      /*
         rhs[4]    load vector due to fsi loads at time t
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

         rhs[4]       / fsi-loads             / fsi-loads             / fsi-loads     / fsi-loads
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


      t0 = ds_cputime();

      /* write memory report */
      if (par.myrank == 0) dsmemreport();


      /* there are only procs allowed in here, that belong to the structural */
      /* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
      if (actintra->intra_fieldtyp != structure) break;


      /* output to the screen */
      if (par.myrank==0)
      {
        printf("Solving STRUCTURE ...\n");
        printf("---------------------------------------------------------------- \n");
      }


      /* copy solution from sol[9][j] to sol[0][j] */
      if (fsidyn->ifsi >= 4  &&  fsiitnum > 0)
        solserv_sol_copy(actfield,disnum_calc,node_array_sol,node_array_sol,9,0);


      /* copy from nodal sol[1][j] to sol[10][j] */
      if (fsidyn->ifsi >= 4  &&  fsiitnum == 0)
        solserv_sol_copy(actfield,disnum_calc,node_array_sol,node_array_sol,1,10);


      /* increment step and time */
      /*sdyn->step++;*/


      /* modifications to time steps size can be done here */
      /* set new absolue time */
      /*sdyn->time += sdyn->dt;*/


      /* put time to global variable for time-dependent load distributions */
      acttime = sdyn->time;


      /* set some constants */
      dyn_setconstants(&dynvar,sdyn,sdyn->dt);


      /* set incremental displacements dispi[0] to zero */
      solserv_zero_vec(&dispi[0]);


      /* set residual displacements in nodes to zero */
      solserv_result_resid(actfield,disnum_calc,actintra,&dispi[0],0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));





      /*----------------------------------------------------------------------*/
      /*                     PREDICTOR                                        */
      /*----------------------------------------------------------------------*/

      /* this vector holds loads due to external forces */
      solserv_zero_vec(&(actsolv->rhs[1]));
      container.kstep = 0;
      container.inherit = 1;
      container.point_neum = 1;
      *action = calc_struct_eleload;
      calrhs(actfield,actsolv,actpart,actintra,stiff_array,
          &(actsolv->rhs[1]),action,&container);


      /* get factor at time t */
      dyn_facfromcurve(actcurve,sdyn->time,&(dynvar.rldfac));


      /* multiply rhs[1] by actual load factor rldfac */
      solserv_scalarprod_vec(&(actsolv->rhs[1]),dynvar.rldfac);


      /* calculate external forces due to fsi coupling */
      /* -> conforming discretization, take values from coincodent nodes */
      if (fsidyn->coupmethod == 1)
      {
        solserv_zero_vec(&(actsolv->rhs[4]));

        switch (fsidyn->coupforce)
        {
          case cf_nodeforce:
          /* determine coupling forces from consistent nodal fluid forces */

            /* initialise full vector to store consistent nodal fluid forces */
            amzero(&fsiforce_a);

            /* get consistent nodal forces from respective fluid node */
            fsi_load(actpart,disnum_calc,fsiforce,numeq_total);

            /* assemble the forces into the rhs[4] */
            assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
                &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[4]),fsiforce,1.0);

            break;


          case cf_stress:
            /* determine coupling forces from stresses */

            container.inherit    = 0;
            container.point_neum = 0;

            *action = calc_struct_fsiload;
            calrhs(actfield,actsolv,actpart,actintra,stiff_array,
                &(actsolv->rhs[4]),action,&container);

            break;


          default: dserror("FSI coupling force type unknown");
        }

      }  /* if (fsidyn->coupmethod == 1) */


      /* mortar method (mtr) */
      else if(fsidyn->coupmethod == 0)
      {
        solserv_zero_vec(&(actsolv->rhs[4]));

        container.inherit = 0;
        container.point_neum = 0;

        *action = calc_struct_fsiload_mtr;
        calrhs(actfield,actsolv,actpart,actintra,stiff_array,
            &(actsolv->rhs[4]),action,&container);

      }
      else
        dserror("fsidyn->coupmethod unknown in fsi_struct! ");



      /* add up the two parts of the external forces
       * and store them in rhs[1]                     */
      solserv_add_vec(&(actsolv->rhs[4]),&(actsolv->rhs[1]),ONE);


      /* set factors needed for prescribed displacement terms on rhs eff */
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



      /* put the prescribed scaled displacements to the nodes in field sol at
       * place 4 separate from the free dofs
       * these are used to calculate the rhs due to dirichlet conditions */
      solserv_putdirich_to_dof(actfield,disnum_calc,0,4,sdyn->time);


      /* put presdisplacements(t) - presdisplacements(t-dt) in place 5 */
      solserv_adddirich(actfield,disnum_calc,0,3,4,5,-1.0,1.0);


      /* calculate tangential stiffness/mass and internal forces at time t-dt */
      solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));
      solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),
          &(actsolv->sysarray_typ[mass_array]));
      amzero(&dirich_a);
      amzero(&intforce_a);

      *action = calc_struct_nlnstiffmass;
      container.dvec          = intforce;
      container.dirich        = dirich;
      container.global_numeq  = numeq_total;
      container.dirichfacs    = dirichfacs;
      container.kstep         = 0;
      calelm(actfield, actsolv, actpart, actintra, stiff_array, mass_array,
          &container, action);


      /* store positive internal forces on fie[1] */
      solserv_zero_vec(&fie[1]);
      assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),&(fie[1]),intforce,1.0);


      /* interpolate external forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
      solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));
      solserv_scalarprod_vec(&(actsolv->rhs[0]),sdyn->alpha_f);
      solserv_add_vec(&(actsolv->rhs[1]),&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));


      /* subtract internal forces from interpolated external forces */
      solserv_add_vec(&(fie[1]),&(actsolv->rhs[0]),-1.0);


      /* add rhs from prescribed displacements to rhs */
      assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);


      /* create effective load vector (rhs[0]-fie[2])eff */
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
      pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
          mass_array,damp_array);



      /* create effective stiffness matrix */
      /*
         keff = constants[6] * K + constants[0] * M + constants[3] * D
         constants[6] =  (1.0-alphaf)
         constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
         constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt)
         */
      kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
          damp_array);



      /* call for solution of system dispi[0] = Keff^-1 * rhs[0] */
      init=0;
      solver_control(actsolv, actintra,
          &(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),
          &(dispi[0]),
          &(actsolv->rhs[0]),
          init);



      /* update displacements
       *   sol[1] = sol[0] + dispi[0] */
      solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
      solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);


      /* put the scaled prescribed displacements to the nodes
       * in field sol at place 0 together with free displacements
       * these are used to calculate the stiffness matrix */
      solserv_putdirich_to_dof(actfield,disnum_calc,0,0,sdyn->time);


      /* return total displacements to the nodes */
      solserv_result_total(actfield,disnum_calc,actintra, &(actsolv->sol[1]),0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* return incremental displacements to the nodes */
      solserv_result_incre(actfield,
          disnum_calc,actintra,&dispi[0],0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));



      /* TODO: */
      /* here put incremental prescribed displacements from sol[5] to sol_increment[0] ? */



      /*----------------------------------------------------------------------*/
      /*                     PERFORM EQUILLIBRIUM ITERATION                   */
      /*----------------------------------------------------------------------*/

      itnum=0;


iterloop:


      /* set factors needed for prescribed displacement terms on rhs eff */
      dirichfacs[0] = -dynvar.constants[0];
      dirichfacs[1] =  dynvar.constants[1];
      dirichfacs[2] =  dynvar.constants[2];
      dirichfacs[3] = -dynvar.constants[3];
      dirichfacs[4] =  dynvar.constants[4];
      dirichfacs[5] =  dynvar.constants[5];
      dirichfacs[6] =  0.0;
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


      /* zero the stiffness matrix and vector for internal forces and dirichlet forces */
      solserv_zero_mat(actintra,&(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));
      solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),
          &(actsolv->sysarray_typ[mass_array]));
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


      /* store positive internal forces on fie[2] */
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


      /* subtract interpolated internal forces from interp. external forces */
      solserv_add_vec(&fie[0],&(actsolv->rhs[0]),-1.0);


      /* add dirichlet forces from prescribed displacements */
      assemble_vec(actintra,&(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),&(actsolv->rhs[0]),dirich,1.0);


      /* create effective load vector (rhs[0]-fie[0])eff */
      pefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,dispi,vel,acc,work,
          mass_array,damp_array);


      /* create effective stiffness matrix */
      kefnln_struct(&dynvar,sdyn,actfield,actsolv,actintra,work,stiff_array,mass_array,
          damp_array);


      /* solve keff * rsd[0] = rhs[0] */
      /* solve for residual displacements to correct incremental displacements*/
      init=0;
      solver_control(actsolv, actintra,
          &(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),
          &(work[0]),
          &(actsolv->rhs[0]),
          init);


      /* return residual displacements to the nodes */
      solserv_result_resid(actfield,disnum_calc,actintra,&work[0],0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* update the incremental displacements by the residual displacements */
      solserv_add_vec(&work[0],&dispi[0],1.0);


      /* update displacements
       *   sol[1] = sol[0] + dispi[0] */
      solserv_copy_vec(&(actsolv->sol[0]),&(actsolv->sol[1]));
      solserv_add_vec(&dispi[0],&(actsolv->sol[1]),1.0);


      /* return total displacements to the nodes */
      solserv_result_total(actfield,disnum_calc,actintra, &(actsolv->sol[1]),0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* return incremental displacements to the nodes */
      solserv_result_incre(actfield,disnum_calc,actintra,&dispi[0],0,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* check for convergence */
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



      /*----------------------------------------------------------------------*
       *                      END OF EQUILLIBRIUM ITERATION                   *
       *----------------------------------------------------------------------*/



      /* for iterative staggered schemes save solution of last iteration
       * copy old total displacments from nodal sol_mf[0][j] to sol_mf[1][j] */
      if (fsidyn->ifsi >= 4)
        solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,0,1);


      /* copy total displacments from nodal sol[0][j] to sol_mf[0][j]  */
      solserv_sol_copy(actfield,disnum_calc,node_array_sol,node_array_sol_mf,0,0);



#ifdef SUBDIV
      /* transfer the solution to the nodes of the master-dis */
      if (actfield->subdivide > 0)
      {
        solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
      }
#endif



      /* for sequential staggered schemes sol_mf[0][j] == sol_mf[1][j] */
      if (fsidyn->ifsi < 4)
        solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,0,1);


      /* print time step */
      if (par.myrank==0)
      {
        printf("| NUMITER = %3d                                                | \n",
            itnum+1);
        printf("---------------------------------------------------------------- \n");
        printf("\n");

        fprintf(out," %3d |",itnum+1);


      }
      if (fsidyn->ifsi>=4)
        break;




      /*======================================================================*
       *                       F I N A L I S I N G                            *
       *======================================================================*/

    case 3:


      /* make temporary copy of actsolv->rhs[2] to actsolv->rhs[0]
       *                        (load at t-dt)
       * because in  dyn_nlnstructupd actsolv->rhs[2] is overwritten but is
       * still needed to compute energies in dynnle */
      solserv_copy_vec(&(actsolv->rhs[2]),&(actsolv->rhs[0]));


      /* update displacements, velocities and accelerations */
      dyn_nlnstructupd(
          actfield,
          disnum_calc,
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


      /* return velocities to sol[1][j] */
      solserv_result_total(actfield,disnum_calc,actintra, &vel[0],1,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* velocities for prescribed dofs */
      solserv_adddirich(actfield,disnum_calc,0,6,0,1,1.0,0.0);


      /* return accel. to sol[2][j] */
      solserv_result_total(actfield,disnum_calc,actintra, &acc[0],2,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      /* accel. for prescribed dofs */
      solserv_adddirich(actfield,disnum_calc,0,7,0,2,1.0,0.0);


      /* make all types of energies */
      dynnle(&dynvar,sdyn,actintra,actsolv,&dispi[0],&fie[1],&fie[2],
          &(actsolv->rhs[1]),&(actsolv->rhs[0]),&work[0]);
      dyne(&dynvar,actintra,actsolv,mass_array,&vel[0],&work[0]);
      dynvar.etot = dynvar.epot + dynvar.ekin;


      if (fsidyn->ichecke>0)
      {

        /* write fsi-forces to nodal sol_mf[4][j] */
        solserv_result_mf(actfield,disnum_calc,actintra, &(actsolv->rhs[4]),4,
            &(actsolv->sysarray[stiff_array]),
            &(actsolv->sysarray_typ[stiff_array]));


        /* write dispi to nodal sol_mf[3][j] */
        solserv_result_mf(actfield,disnum_calc,actintra,dispi,3,
            &(actsolv->sysarray[stiff_array]),
            &(actsolv->sysarray_typ[stiff_array]));


        /* energy transported over the interface */
        fsi_dyneint(actfield,disnum_calc,0);


        /* copy old fsi-forces from nodal sol_mf[4][j] to sol_mf[5][j] */
        solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,4,5);
      }



      /* save actual solution as old solution
       * copy from nodal sol[0][j] to sol[9][j] */
      solserv_sol_copy(actfield,disnum_calc,node_array_sol,node_array_sol,0,9);


      /* save actual relaxed solution as old relaxed solution
       * copy from nodal sol_mf[1][j] to sol_mf[2][j] */
      solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,1,2);


      /* perform stress calculation  and print out results to .out */
      outstep++;
      if (outstep == sdyn->updevry_disp)
      {
        outstep=0;
        if (ioflags.struct_stress==1)
        {
          *action = calc_struct_stress;
          container.dvec          = NULL;
          container.dirich        = NULL;
          container.global_numeq  = 0;
          container.dirichfacs    = NULL;
          container.kstep         = 0;
          calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);

          /* reduce stresses, so they can be written */
          *action = calc_struct_stressreduce;
          container.kstep = 0;
          calreduce(actfield,actpart,disnum_calc,actintra,action,&container);
        }

        if ( (ioflags.struct_stress == 1  ||  ioflags.struct_disp == 1) &&
            ioflags.output_out == 1)
          out_sol(actfield,actpart,disnum_calc,actintra,sdyn->step,0);

      }


      /* write restart data to pss file */
      restartstep++;
      if (restartstep==fsidyn->uprestart)
      {
        restartstep=0;

#ifdef BINIO
        if(disnum_io != disnum_calc)
          restart_write_bin_nlnstructdyn(&restart_context,
            sdyn,
            &dynvar,
            actsolv->nrhs, actsolv->rhs,
            actsolv->nsol, actsolv->sol,
            1            , dispi       ,
            1            , vel         ,
            1            , acc         ,
            3            , fie         ,
            3            , work);
        else
          restart_write_bin_nlnstructdyn(&out_context,
            sdyn,
            &dynvar,
            actsolv->nrhs, actsolv->rhs,
            actsolv->nsol, actsolv->sol,
            1            , dispi       ,
            1            , vel         ,
            1            , acc         ,
            3            , fie         ,
            3            , work);
#else
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
#endif
      }


      /* print time step */
      /*if (par.myrank==0)
        {
        dyn_nlnstruct_outstep(&dynvar,sdyn,itnum);
        printf("--------------------------------------------------------------- \n");
        printf("\n");
        } */


      /* measure time for this step */
      t1 = ds_cputime();
      fprintf(allfiles.out_err,"TIME for step %d is %f sec\n",sdyn->step,t1-t0);


      break;





      /*======================================================================*
       *  Binary Output                                                       *
       *======================================================================*/

    case 98:


#ifdef BINIO
      if (ioflags.output_bin)
      {

        if (ioflags.struct_disp==1)
        {
          out_results(&out_context, sdyn->time, sdyn->step, 0, OUTPUT_DISPLACEMENT);

          /* This was not implemented before. I doubt it works. */
          /*
          out_results(&out_context, sdyn->time, sdyn->step, 1, OUTPUT_VELOCITY);
          out_results(&out_context, sdyn->time, sdyn->step, 2, OUTPUT_ACCELERATION);
          */
        }

        if (ioflags.struct_stress==1)
        {
          out_results(&out_context, sdyn->time, sdyn->step, 0, OUTPUT_STRESS);
        }
      }
#endif

      break;



      /*======================================================================*
       *                C L E A N I N G   U P   P H A S E                     *
       *======================================================================*/

    case 99:


      /* there are only procs allowed in here, that belong to the structural */
      /* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
      if (actintra->intra_fieldtyp != structure) break;


      /* print out the final solution */
      if (outstep!=0)
      {

        outstep = 0;

        if (ioflags.struct_stress==1)
        {
          *action = calc_struct_stress;
          container.dvec          = NULL;
          container.dirich        = NULL;
          container.global_numeq  = 0;
          container.dirichfacs    = NULL;
          container.kstep         = 0;
          calelm(actfield,actsolv,actpart,actintra,stiff_array,-1,&container,action);


          /* reduce stresses, so they can be written */
          *action = calc_struct_stressreduce;
          container.kstep = 0;
          calreduce(actfield,actpart,disnum_calc,actintra,action,&container);

        }

        if ( (ioflags.struct_stress==1 || ioflags.struct_disp==1) && ioflags.output_out==1)
          out_sol(actfield,actpart,disnum_io,actintra,sdyn->step,0);
      }


      amdel(&intforce_a);
      amdel(&fsiforce_a);
      solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
      solserv_del_vec(&(actsolv->sol),actsolv->nsol);
      solserv_del_vec(&dispi,1);
      solserv_del_vec(&vel,1);
      solserv_del_vec(&acc,1);
      solserv_del_vec(&fie,3);
      solserv_del_vec(&work,3);


#ifdef BINIO
      destroy_bin_out_field(&out_context);
      if(disnum_io != disnum_calc)
        destroy_bin_out_field(&restart_context);
#endif


#ifndef PARALLEL
      CCAFREE(actintra);
#endif

      break;



      /*======================================================================*
       *               F S I - P R E D I C T O R   P H A S E                  *
       *======================================================================*/

    case 4:

      fsi_structpredictor(actfield,disnum_calc,0);

      break;





      /*======================================================================*
        |   S O L U T I O N   F O R   R E L A X A T I O N   P A R A M E T E R  |
        |                      using steepest descent method                   |
       *======================================================================*
       * nodal solution history structural field:                             *
       * sol[0][j]           ... total displacements at time (t)              *
       * sol[1][j]           ... velocities at time (t)                       *
       * sol[2][j]           ... accels at time (t)                           *
       * sol[3][j]           ... prescribed displacements at time (t-dt)      *
       * sol[4][j]           ... prescribed displacements at time (t)         *
       * sol[5][j]           ... place 4 - place 3                            *
       * sol[6][j]           ... the  velocities of prescribed dofs           *
       * sol[7][j]           ... the  accels of prescribed dofs               *
       * sol[8][j]           ... working space                                *
       * sol[9][j]           ... total displacements at time (t-dt)           *
       * sol[10][j]          ... velocities at time (t-dt)                    *
       * sol_mf[0][j]        ... latest struct-displacements                  *
       * sol_mf[1][j]        ... (relaxed) displ. of the last iteration step  *
       * sol_mf[2][j]        ... converged relaxed displ. at time (t-dt)      *
       * sol_mf[3][j]        ... actual dispi                                 *
       * sol_mf[4][j]        ... FSI coupl.-forces at the end of the timestep *
       * sol_mf[5][j]        ... FSI coupl.-forces at beginning of the timest.*
       *======================================================================*/

      /*

         rhs[4]    load vector due to fsi loads at time t
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

         This calculation in performed in one single step. There are no external
         forces and no time dependencies.

*/

    case 6:


      if (fsidyn->ifsi != 6)
        dserror("No auxiliary structure solution within this coupling scheme");


      /* there are only procs allowed in here, that belong to the structural */
      /* intracommunicator (in case of nonlinear struct. dyn., this should be all) */
      if (actintra->intra_fieldtyp != structure) break;


      /* output to the screen */
      if (par.myrank==0)
      {
        printf("          - Solving STRUCTURE ...\n");
      }


      /* set incremental displacements dispi[0] to zero */
      solserv_zero_vec(&dispi[0]);


      /* calculate rhs from external forces due to fsi coupling */
      solserv_zero_vec(&(actsolv->rhs[0]));

      container.inherit = 1;
      container.point_neum = 1;

      *action = calc_struct_fsiload;
      calrhs(actfield,actsolv,actpart,actintra,stiff_array,
          &(actsolv->rhs[0]),action,&container);


      /* note: this calculation is performed on the initial configuration
               changed by the increment step vector g_i only. There are no
               Dirichlet boundary conditions different from zero */

      /* multiply rhs by (1.0-alpha_f) */
      solserv_scalarprod_vec(&(actsolv->rhs[0]),(1.0-sdyn->alpha_f));


      /* solve with Keff from previous system soluiton */

      /* call for solution of system dispi[0] = Keff^-1 * rhs[0] */
      init=0;
      solver_control(actsolv, actintra,
          &(actsolv->sysarray_typ[stiff_array]),
          &(actsolv->sysarray[stiff_array]),
          &(dispi[0]),
          &(actsolv->rhs[0]),
          init);


      /* write solution to the nodes (sol[8][i]) */
      solserv_result_total(actfield,disnum_calc,actintra, &(dispi[0]),8,
          &(actsolv->sysarray[stiff_array]),
          &(actsolv->sysarray_typ[stiff_array]));


      break;



    default:
      dserror("Parameter out of range: mctrl \n");
  } /* end switch (mctrl) */



#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of dyn_nln_structural */


#endif   /* ofdef D_FSI */



/*! @} (documentation module close)*/





