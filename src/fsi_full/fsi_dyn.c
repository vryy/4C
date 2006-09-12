/*!----------------------------------------------------------------------
\file
\brief control function for FSI

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/


#ifdef D_FSI


/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/



#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#include "../io/io.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;


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


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*!---------------------------------------------------------------------
\brief routine to control fsi dynamic analyis

<pre>                                                         genk 09/02

Nonlinear dynamical algorithms for FSI-problems:
Implemented Algorithms:
 - basic sequentiel staggered scheme
 - sequential staggered scheme with predictor
 - iterative staggered scheme with fixed relaxation parameter
 - iterative staggered scheme with relaxation parameter via AITKEN
   iteration

</pre>

\param mctrl    INT   (i)     evaluation flag

\return void

------------------------------------------------------------------------*/
void dyn_fsi(
    INT                 mctrl
    )

{

  static INT             numfld;       /* number of fiels */
  static INT             numsf;
  static INT             numff;
  static INT             numaf;        /* actual number of fields */
  static INT             resstep=0;    /* counter for output control  */
  static INT             restartstep=0;/* counter for restart control */
  INT                    actcurve;     /* actual curve */
  INT                    itnum=0;      /* iteration counter */
  INT                    converged;    /* convergence flag */

  static FIELD          *fluidfield;
  static FIELD          *structfield;
  static FIELD          *alefield;
  static FLUID_DYNAMIC  *fdyn;
  static FSI_DYNAMIC    *fsidyn;
  static STRUCT_DYNAMIC *sdyn;
  static ALE_DYNAMIC    *adyn;


  static FSI_STRUCT_WORK struct_work;
  static FSI_FLUID_WORK  fluid_work;
  static FSI_ALE_WORK    ale_work;

  FILE                  *out = allfiles.out_out;

  DOUBLE                 t2,tt;


  INT                    f_disnum_calc = 0;
  INT                    f_disnum_io   = 0;

  INT                    s_disnum_calc = 0;
  INT                    s_disnum_io   = 0;

  INT                    a_disnum_calc = 0;
  INT                    a_disnum_io   = 0;



#ifdef D_MORTAR
  INTERFACES            *int_faces; /* interface information for mortar   */
#endif



#ifdef DEBUG
  dstrc_enter("dyn_fsi");
#endif



#ifdef PERF
  perf_begin(40);
#endif



  if (mctrl==99) goto cleaningup;


  /* check input -------------------------------------------------------- *
   * convention used at the moment:                                       *
   *   FIELD 0: structure                                                 *
   *   FIELD 1: fluid                                                     *
   *   FIELD 2: mesh / ale                                                *
   * -------------------------------------------------------------------- */
  numfld = genprob.numfld;
  dsassert(numfld==3,"Three fields needed for FSI-problem!\n");

  structfield = &(field[0]);
  numsf       = 0;
  fluidfield  = &(field[1]);
  numff       = 1;
  alefield    = &(field[2]);
  numaf       = 2;

  /* plausibility check */
  dsassert(structfield->fieldtyp==structure,"FIELD 0 has to be structure\n");
  dsassert(fluidfield->fieldtyp==fluid,"FIELD 1 has to be fluid\n");
  dsassert(alefield->fieldtyp==ale,"FIELD 2 has to be ale\n");





  /*======================================================================*
    I N I T I A L I S A T I O N
   *======================================================================*/

#ifdef PERF
  perf_begin(41);
#endif



  /* check for subdivision in fluid field */
#ifdef SUBDIV
  if (fluidfield->subdivide > 0)
  {
    f_disnum_calc = 1;
    f_disnum_io   = ioflags.output_dis;
  }
  else
#endif
    f_disnum_calc = f_disnum_io = 0;


  /* check for subdivision in structure field */
#ifdef SUBDIV
  if (structfield->subdivide > 0)
  {
    s_disnum_calc = 1;
    s_disnum_io   = ioflags.output_dis;
  }
  else
#endif
    s_disnum_calc = s_disnum_io = 0;


  /* check for subdivision in ale field */
#ifdef SUBDIV
  if (alefield->subdivide > 0)
  {
    a_disnum_calc = 1;
    a_disnum_io   = ioflags.output_dis;
  }
  else
#endif
    a_disnum_calc = a_disnum_io = 0;



  mctrl=1;
  sdyn= alldyn[0].sdyn;
  fdyn= alldyn[1].fdyn;
  adyn= alldyn[2].adyn;
  fsidyn= alldyn[3].fsidyn;
  fsidyn->time=ZERO;
  fsidyn->step=0;

  adyn->coupmethod = fsidyn->coupmethod;

  if (fdyn->freesurf==1)
    dserror("No explicit free surface combined with FSI!");


  /* initialise fsi coupling conditions */
#ifdef SUBDIV
  if (alefield->subdivide > 0)
  {
    fsi_initcoupling(
        structfield, s_disnum_io,
        fluidfield,  f_disnum_io,
        alefield,    a_disnum_io);

  }
#endif

  fsi_initcoupling(
      structfield, s_disnum_calc,
      fluidfield,  f_disnum_calc,
      alefield,    a_disnum_calc);



#ifdef D_MORTAR
  int_faces = (INTERFACES*)CCACALLOC(1,sizeof(INTERFACES));


  if (fsidyn->coupmethod == 0) /* mortar method */
  {
    /*-------------------------------------------- intitialize interfaces */
    fsi_initcoupling_intfaces(structfield, s_disnum_calc,
        fluidfield, f_disnum_calc, int_faces);
    /* ------------------------allocate memory for the several interfaces */
    int_faces->interface = (INTERFACE*)CCACALLOC(int_faces->numint,
        sizeof(INTERFACE));
  }
#endif


  /* determine structural interface dofs */
  fsi_struct_intdofs(structfield, s_disnum_calc);


  /* init all applied time curves */
  for (actcurve = 0;actcurve<numcurve;actcurve++)
    dyn_init_curve(actcurve,fsidyn->nstep,fsidyn->dt,fsidyn->maxtime);


  /* initialise ale */
  fsi_ale_setup(&ale_work,alefield,a_disnum_calc,a_disnum_io);


  /* initialise fluid */
  fsi_fluid_setup(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);


  /* initialise structure */
  fsi_struct_setup(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);




  if (genprob.restart!=0)
  {

#if defined(BINIO)
    restart_read_bin_fsidyn(fsidyn, genprob.restart);
#else
    restart_read_fsidyn(genprob.restart,fsidyn);
#endif

    /* plausibility check */
    if (fsidyn->time != adyn->time ||
        fsidyn->time != fdyn->acttime ||
        fsidyn->time != sdyn->time   )
      dserror("Restart problem: Time not identical in fields!\n");

    if (fsidyn->step != fdyn->step ||
        fsidyn->step != adyn->step ||
        fsidyn->step != sdyn->step   )
      dserror("Restart problem: Step not identical in fields!\n");
  }



  /* initialise AITKEN iteration */
  if (fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param)
  {
    fsi_aitken(structfield,s_disnum_calc,itnum,0);
  }
  else if (fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force)
  {
    fsi_aitken_force(structfield,s_disnum_calc,fluidfield,f_disnum_calc,itnum,numff,0);
  }


#ifdef SUBDIV
/* transfer the solution to the nodes of the master-dis */
if (fluidfield->subdivide > 0)
  solserv_sol_trans(fluidfield,  f_disnum_calc, node_array_sol, 0);
if (structfield->subdivide > 0)
  solserv_sol_trans(structfield, s_disnum_calc, node_array_sol, 0);
if (alefield->subdivide > 0)
  solserv_sol_trans(alefield,    a_disnum_calc, node_array_sol, 0);
#endif



  /* print the mesh to GiD */
  if (ioflags.output_gid==1 && par.myrank==0)
    out_gid_msh();



  /* write initial solution to gid */
  /*   print out solution to 0.flavia.res file */
  if (ioflags.output_gid==1 && par.myrank==0)
    out_gid_sol_fsi(fluidfield,structfield,f_disnum_io,s_disnum_io);



  /*
   * Binary output has to be done by the algorithms because the
   * contexts are there. */
  if (ioflags.output_bin==1)
  {
    mctrl = 98;
    fsi_ale_output(&ale_work,alefield,a_disnum_calc,a_disnum_io);
    fsi_fluid_output(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
    fsi_struct_output(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
  }



#ifdef D_MORTAR
  if (fsidyn->coupmethod == 0) /* mortar method */
  {
    /*----------------------------- fill structure interfaces with values */
    fsi_init_interfaces(structfield,s_disnum_calc,fluidfield,f_disnum_calc,int_faces);
  }
#endif



  /* write general data to .out */
  if (par.myrank==0)
  {
    fprintf(out,"max. values:\n");
    fprintf(out,"============\n");


    /* table head */
    fprintf(out," time |            |field|fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
      case fncc_Linf: /* infinity norm */
        fprintf(out,"inf-norm");
        break;
      case fncc_L1: /* L_1 norm */
        fprintf(out,"L_1-norm");
        break;
      case fncc_L2: /* L_2 norm */
        fprintf(out,"L_2-norm");
        break;
      default:
        dserror("Norm for nonlin. convergence check unknown!!\n");
    }  /* switch(fdyn->itnorm) */

    fprintf(out," |struc| convergence| relaxation |   total    |\n");

    fprintf(out," step |  sim. time | ite | ite |     vel.   |     pre.   | ite | over fields|  parameter | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");



    /* max values */
    fprintf(out,"%5d | %10.3f | %3d | %3d |        %10.3E       | %3d | %10.3E |            |            |\n",
        fdyn->nstep,fdyn->maxtime,fsidyn->itemax,fdyn->itemax,fdyn->ittol,sdyn->maxiter,fsidyn->convtol);
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");



    fprintf(out,"\n\ntimeloop:  ");

    switch (fsidyn->ifsi)
    {
      case fsi_basic_sequ_stagg:
        fprintf(out,"Basic Sequential Staggered Scheme\n");
        break;
      case fsi_sequ_stagg_pred:
        fprintf(out,"Sequential Staggered Scheme with Predictor\n");
        break;
      case fsi_iter_stagg_fixed_rel_param:
        fprintf(out,"Iterative Staggered Scheme with Fixed Relaxation Parameter\n");
        break;
      case fsi_iter_stagg_AITKEN_rel_param:
        fprintf(out,"Iterative Staggered Scheme with Relaxation Parameter via Aitken Iteration\n L: %d \n", fsidyn->ifsi);
        break;
      case fsi_iter_stagg_steep_desc:
        fprintf(out,"Iterative Staggered Scheme with Relaxation Parameter via Steepest Descent Method\n");
        break;
      case fsi_iter_stagg_AITKEN_rel_force:
        fprintf(out,"Iterative Staggered Scheme with relaxation of interface forces via Aitken parameter\n");
        break;
      case fsi_iter_stagg_steep_desc_force:
        fprintf(out,"Iterative Staggered Scheme with relaxation of interface forces via Steepest Descent Method\n");
        break;
#ifdef FSI_NEWTONCOUPLING
      case fsi_iter_stagg_Newton_FD:
        fprintf(out,"Iterative Staggered Scheme with Newton-Method - Approximation by finite Differenc\n");
        break;
      case fsi_iter_stagg_Newton_I:
        fprintf(out,"Iterative Staggered Scheme with Newton-Method - Approximation by Identity Matrix\n");
        break;
#endif
      default:
        dserror("algoout not implemented yet\n");
    }
    fprintf(out,"=========\n");



    /* table head */
    fprintf(out," time |            |field|fluid| fluid error in ");

    switch(fdyn->itnorm)
    {
      case fncc_Linf: /* infinity norm */
        fprintf(out,"inf-norm");
        break;
      case fncc_L1: /* L_1 norm */
        fprintf(out,"L_1-norm");
        break;
      case fncc_L2: /* L_2 norm */
        fprintf(out,"L_2-norm");
        break;
      default:
        dserror("Norm for nonlin. convergence check unknown!!\n");
    }  /* switch(fdyn->itnorm) */

    fprintf(out," |struc| convergence| relaxation |   total    |\n");

    fprintf(out," step |  sim. time | ite | ite |     vel.   |     pre.   | ite | over fields|  parameter | calc. time |\n");
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");


  }  /* if (par.myrank==0) */

  fflush(out);


#ifdef PERF
  perf_end(41);
#endif





  /*======================================================================*
    T I M E L O O P
   *======================================================================*/
timeloop:

  t2=ds_cputime();


  mctrl=2;

  fsidyn->step++;
  fsidyn->time += fsidyn->dt;

  fdyn->step=fsidyn->step;
  sdyn->step=fsidyn->step;
  adyn->step=fsidyn->step;

  fdyn->acttime = fsidyn->time;
  sdyn->time = fsidyn->time;
  adyn->time = fsidyn->time;



  /*======================================================================*
    do the iteration over the fields within one timestep
   *======================================================================*/
  itnum=0;
fielditer:


  /* output to the screen */
  if (par.myrank==0) fsi_algoout(itnum);

  if (par.myrank==0)
    fprintf(out,"%5d | %10.3f | %3d |",fsidyn->step,fsidyn->time,itnum);



  switch (fsidyn->ifsi)
  {

    /* basic sequential staggered scheme */
    /*===================================*/

  case fsi_basic_sequ_stagg:
  case fsi_sequ_stagg_shift:
  {
    dsassert(fsidyn->ifsi!=fsi_sequ_stagg_shift,"Scheme with DT/2-shift not implemented yet!\n");

    /*------------------------------- CFD -------------------------------*/
    perf_begin(42);
    fsi_fluid_calc(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
    fsi_fluid_final(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
    perf_end(42);

    /*------------------------------- CSD -------------------------------*/
    perf_begin(43);
    fsi_struct_calc(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum,fluidfield,f_disnum_calc);
    fsi_struct_final(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
    perf_end(43);

    /*------------------------------- CMD -------------------------------*/
    perf_begin(44);
    fsi_ale_calc(&ale_work,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);
    fsi_ale_final(&ale_work,alefield,a_disnum_calc,a_disnum_io);
    perf_end(44);

    break;
  }


  /* schemes with predictor */
  /*========================*/

  case fsi_sequ_stagg_pred:
  case fsi_iter_stagg_fixed_rel_param:
  case fsi_iter_stagg_AITKEN_rel_param:
  case fsi_iter_stagg_steep_desc:
  case fsi_iter_stagg_CHEB_rel_param:
#ifdef FSI_NEWTONCOUPLING
  case fsi_iter_stagg_Newton_FD:
  case fsi_iter_stagg_Newton_I:
#endif
  {
#ifdef D_MORTAR
    if (fsidyn->coupmethod == 0) /* mortar method */
    {
      /*-- computation of mortar approximation, only for interfaces with */
      /* --------------nonconforming discretizations */
      fsi_mortar_coeff(fsidyn, int_faces);
    }
#endif

    /*----------------- CSD - predictor for itnum==0 --------------------*/
    if (itnum == 0)
    {
      perf_begin(43);
      fsi_structpredictor(structfield,s_disnum_calc,0);
      perf_end(43);
    }

#ifdef D_MORTAR
    if (fsidyn->coupmethod == 0) /* mortar method */
    {
      /*-- computation of interface displacements for ale nodes, only for*/
      /* ----------------- interfaces with nonconforming discretizations */
      fsi_calc_disp4ale(fsidyn, int_faces);
    }
#endif

#ifdef FSI_NEWTONCOUPLING
    /* backup for Newton via finite differences */
    {
      ARRAY_POSITION *ipos;
      ipos = &(structfield->dis[s_disnum_calc].ipos);
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, ipos->mf_dispnp, 8);
      solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, ipos->mf_reldisp, 9);
      if (itnum > 0)
        solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,9,11);
      if (itnum == 0)
        solserv_sol_copy(structfield,s_disnum_calc,node_array_sol,node_array_sol,1,12);
    }
#endif

    /*------------------------------- CMD -------------------------------*/
    perf_begin(44);
    fsi_ale_calc(&ale_work,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);
    perf_begin(44);

    /*------------------------------- CFD -------------------------------*/
    perf_begin(42);
    fsi_fluid_calc(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
    perf_end(42);

#ifdef D_MORTAR
    if (fsidyn->coupmethod == 0) /* mortar method */
    {
      fsi_calc_intforces(int_faces);
      /* --------put coupling forces from fluid nodes to structure nodes */
      fsi_put_coupforc2struct(structfield, s_disnum_calc, int_faces);
    }
#endif

    /*------------------------------- CSD -------------------------------*/
    perf_begin(43);
    fsi_struct_calc(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum,fluidfield,f_disnum_calc);
    perf_end(43);

    break;
  }

  case fsi_iter_stagg_AITKEN_rel_force:
  case fsi_iter_stagg_steep_desc_force:
  {
    /*------------------------------- CSD -------------------------------*/
    perf_begin(43);
    fsi_struct_calc(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum,fluidfield,f_disnum_calc);
    perf_end(43);

    /*------------------------------- CMD -------------------------------*/
    perf_begin(44);
    fsi_ale_calc(&ale_work,alefield,a_disnum_calc,a_disnum_io,structfield,s_disnum_calc);
    perf_begin(44);

    /*------------------------------- CFD -------------------------------*/
    perf_begin(42);
    fsi_fluid_calc(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io,alefield,a_disnum_calc);
    perf_end(42);

    break;
  }

  default:
    dserror("coupling method %d not supported", fsidyn->ifsi);
  }

#ifdef FSI_NEWTONCOUPLING

  switch (fsidyn->ifsi)
  {

    /* Newton methods */
    /*================*/

  case fsi_iter_stagg_Newton_I:
  case fsi_iter_stagg_Newton_FD:
  {
/*     perf_begin(61); */
    solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, ipos->mf_reldisp, 7);

    /* Calculate right Side*/
    double *rS;/* Vector for right Side */
    rS = fsi_newton_rightSide(structfield, s_disnum_calc);

    /* Solving the linear Syste*/
    double *deltaD;
    deltaD = fsi_newton_SolveSystem(rS, structfield, fluidfield, alefield, &struct_work, &fluid_work, &ale_work, s_disnum_calc, s_disnum_io, f_disnum_calc, f_disnum_io, a_disnum_calc, a_disnum_io, itnum);

    /* Update of the interface displacement */
    fsi_newton_final(deltaD, structfield,  s_disnum_calc);

    solserv_sol_copy(structfield, s_disnum_calc, node_array_sol_mf, node_array_sol_mf, 7, ipos->mf_reldisp);

/*     perf_end(61); */
    break;
  }
  default:
  {
  }
  }

#endif

  switch (fsidyn->ifsi)
  {

    /* strong coupling schemes */
    /*=========================*/

  case fsi_iter_stagg_fixed_rel_param:
  case fsi_iter_stagg_AITKEN_rel_param:
  case fsi_iter_stagg_steep_desc:
  case fsi_iter_stagg_CHEB_rel_param:
#ifdef FSI_NEWTONCOUPLING
  case fsi_iter_stagg_Newton_FD:
  case fsi_iter_stagg_Newton_I:
#endif
  {
    DOUBLE resnorm;

    /* iteration convergence check */
    converged = fsi_convcheck(structfield, s_disnum_calc, itnum, &resnorm);


    if (converged == 0) /* no convergence */
    {
      /* compute optimal relaxation parameter */

      perf_begin(45);

      switch (fsidyn->ifsi)
      {
      case fsi_iter_stagg_fixed_rel_param:
      case fsi_iter_stagg_Newton_FD:
      case fsi_iter_stagg_Newton_I:
        /* Nothing to do. */
        break;
      case fsi_iter_stagg_AITKEN_rel_param:
        fsi_aitken(structfield,s_disnum_calc,itnum,1);
        break;
      case fsi_iter_stagg_steep_desc:
        fsi_gradient(&struct_work,&fluid_work,&ale_work,
		     alefield,structfield,fluidfield,
                     a_disnum_io,a_disnum_calc,
                     s_disnum_io,s_disnum_calc,
                     f_disnum_io,f_disnum_calc,
                     numaf,numff,numsf);
        break;
      case fsi_iter_stagg_CHEB_rel_param:
        dserror("RELAX via CHEBYCHEV not implemented yet!\n");
        break;
      default:
        dserror("Ups");
      }

      if (par.myrank==0)
        fprintf(out,"   %7.5f  |",fsidyn->relax);

      /* relaxation of structural interface displacements */
      fsi_relax_intdisp(structfield,s_disnum_calc);
      perf_end(45);


      itnum++;

      if (par.myrank==0)
        fprintf(out,"            |\n");
      fflush(out);

      goto fielditer;
    }
    else /* convergence */
    {

      if (par.myrank==0)
        fprintf(out,"            |");

      mctrl=3;

      /*--------------------- update MESH data -------------------------*/
      perf_begin(44);
      fsi_ale_final(&ale_work,alefield,a_disnum_calc,a_disnum_io);
      perf_end(44);

      /*-------------------- update FLUID data -------------------------*/
      perf_begin(42);
      fsi_fluid_final(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
      perf_end(42);

      /*------------------ update STRUCTURE data -----------------------*/
      perf_begin(43);
      fsi_struct_final(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
      perf_end(43);


    }
    break;
  }

    /* strong coupling schemes with force relaxation */
    /*===============================================*/

  case fsi_iter_stagg_AITKEN_rel_force:
  case fsi_iter_stagg_steep_desc_force:
  {
    /* iteration convergence check */
    converged = fsi_convcheck_force(structfield, s_disnum_calc, fluidfield, f_disnum_calc, itnum, numff);

    if (converged == 0) /* no convergence */
    {
      /* compute optimal relaxation parameter */
      perf_begin(45);

      switch (fsidyn->ifsi)
      {
      case fsi_iter_stagg_AITKEN_rel_force:
        fsi_aitken_force(structfield,s_disnum_calc,fluidfield,f_disnum_calc,itnum,numff,1);
        break;
      case fsi_iter_stagg_steep_desc_force:
        fsi_gradient_force(&struct_work,&fluid_work,&ale_work,
                           alefield,structfield,fluidfield,
                           a_disnum_io,a_disnum_calc,
                           s_disnum_io,s_disnum_calc,
                           f_disnum_io,f_disnum_calc,
                           numaf,numff,numsf);
        break;
      default:
        dserror("Ups");
      }

      if (par.myrank==0)
        fprintf(out,"   %7.5f  |",fsidyn->relax);

      fsi_relax_intdisp_force(structfield,s_disnum_calc,fluidfield,f_disnum_calc,numff);

      perf_end(45);

      itnum++;
      if (par.myrank==0)
        fprintf(out,"            |\n");
      fflush(out);

      goto fielditer;
    }
    else /* convergence of iteration with relaxed interface forces */
    {
      if (par.myrank==0)
        fprintf(out,"            |");
      mctrl=33;

      /*--------------------- update MESH data -------------------------*/
      perf_begin(44);
      fsi_ale_final(&ale_work,alefield,a_disnum_calc,a_disnum_io);
      perf_end(44);

      /*-------------------- update FLUID data -------------------------*/
      perf_begin(42);
      fsi_fluid_final(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
      perf_end(42);

      /*------------------ update STRUCTURE data -----------------------*/
      perf_begin(43);
      fsi_struct_final(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
      perf_end(43);

      if (fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force)
      {
        fsi_aitken_force(structfield,s_disnum_calc,fluidfield,f_disnum_calc,itnum,numff,2);
      }
    }
    break;
  }
  default:
  {
    if (par.myrank==0)
      fprintf(out,"            |            |");
  }
  }


  tt=ds_cputime()-t2;
  if (par.myrank==0)
  {
    fprintf(out," %10.3f |\n",tt);
    fprintf(out,"-------------------------------------------------------------------------------------------------------\n");
  }
  fflush(out);



  /* write current solution */
  resstep++;
  restartstep++;

  if (resstep == fsidyn->upres)
  {
    resstep=0;

    /* print out solution to GiD */
    if (ioflags.output_gid==1 && par.myrank==0)
    {
      out_checkfilesize(1);
      out_gid_sol_fsi(fluidfield,structfield,f_disnum_io,s_disnum_io);
    }


#ifdef BINIO
    /*
     * Binary output has to be done by the algorithms because the
     * contexts are there. */
    if (ioflags.output_bin==1)
    {
      mctrl = 98;
      fsi_ale_output(&ale_work,alefield,a_disnum_calc,a_disnum_io);
      fsi_fluid_output(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
      fsi_struct_output(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
    }
#endif

  }  /* if (resstep == fsidyn->upres) */



  /* write fsi-restart data */
  if (restartstep==fsidyn->uprestart)
  {
    restartstep=0;
#ifdef BINIO
    restart_write_bin_fsidyn(fsidyn);
#else
    restart_write_fsidyn(fsidyn);
#endif
  }



  /* energy check */
  if (fsidyn->ichecke>0) fsi_energycheck();



  /* finalising this time step */
  if (fsidyn->step < fsidyn->nstep && fsidyn->time <= fsidyn->maxtime)
    goto timeloop;



  /*======================================================================*
    C L E A N I N G   U P   P H A S E
   *======================================================================*/
cleaningup:

  mctrl=99;
  fsi_fluid_cleanup(&fluid_work,fluidfield,f_disnum_calc,f_disnum_io);
  fsi_struct_cleanup(&struct_work,structfield,s_disnum_calc,s_disnum_io,itnum);
  fsi_ale_cleanup(&ale_work,alefield,a_disnum_calc,a_disnum_io);




#ifdef PERF
  perf_end(40);
#endif



#ifdef DEBUG
  dstrc_exit();
#endif

  return;

} /* end of dyn_fsi */


#endif  /* ifdef D_FSI */


/*! @} (documentation module close)*/
