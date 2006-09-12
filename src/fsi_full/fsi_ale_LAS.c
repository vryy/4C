/*!----------------------------------------------------------------------
\file
\brief ale part of fsi-problems

*----------------------------------------------------------------------*/
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
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
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
 extern struct _FILES  allfiles;


void fsi_ale_LAS_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT              i,j;
  INT              mone=-1;
  DOUBLE     H = THREE*TEN;
  DOUBLE           x1,y1,x2,y2;
  INT       numaf;
  INT       numnp_total;
  INT       numff;
  INT      *index;
  DISCRET  *actdis;
  NODE            *actnode;
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  FSI_DYNAMIC  *fsidyn;
  ALE_DYNAMIC  *adyn;
  ARRAY_POSITION *ipos;

#ifndef PARALLEL
  work->dummy_intra.intra_fieldtyp = ale;
  work->dummy_intra.intra_rank     = 0;
  work->dummy_intra.intra_nprocs   = 1;
#endif

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  actpart     = &(partition[numaf]);
  ipos   = &(actfield->dis[disnum_calc].ipos);

  /*====================================================================*
   * nodal solution history ale field:                                  *
   * sol[1...0][j]  ... solution for visualisation (real pressure)	*
   * sol_mf[0][i]        ... displacements at (n)			*
   * sol_mf[1][i]        ... displacements at (n+1) 		        *
   * sol_mf[2][i]        ... used in fsi_gradient.c                     *
   * sol_mf[3][i]        ... displacements at (n+1) stored in fsi_calc_disp4ale*
   *====================================================================*/
  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  numff        = genprob.numff;
  numnp_total = actfield->dis[disnum_calc].numnp;
  actdis       = &(actfield->dis[disnum_calc]);
  work->outstep=0;
  work->restartstep=0;

/*-------------------------------------------------- define index array */
  index = amdef("index",&work->index_a,numnp_total,1,"IV");
  aminit(&work->index_a,&mone);

  for (i=0;i<numnp_total;i++)
  {
    actnode = &(actdis->node[i]);
    x1 = actnode->x[0];
    y1 = actnode->x[1];
    if (FABS(y1-H)<EPS8)
      index[i]=i;
    else
    {
      for (j=0;j<numnp_total;j++)
      {
        actnode = &(actdis->node[j]);
        x2 = actnode->x[0];
        y2 = actnode->x[1];
        if (FABS(y2-H)<EPS8 && FABS(x1-x2)<EPS8)
          index[i]=j;
      }
    }
  }

/*-------------------------------------------------------------- check */
  for (i=0;i<numnp_total;i++)
    if (index[i]==mone)
      dserror("something went wrong!");

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
  init_bin_out_field(&work->out_context,
                     NULL,
                     NULL,
                     actfield, actpart, actintra, disnum_io);
  if (disnum_io != disnum_calc)
    init_bin_out_field(&work->restart_context,
                       NULL,
                       NULL,
                       actfield, actpart, actintra, disnum_calc);
#endif

/*---------------------------------------------------------- monitoring */
  if (ioflags.monitor==1)
  {
    out_monitor(actfield,numaf,ZERO,1);
    monitoring(actfield,disnum_calc,numaf,0,adyn->time);
  }

/*------------------------------------------- print out results to .out */
  if (ioflags.ale_disp==1 && ioflags.output_out==1)
  {
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }
}


void fsi_ale_LAS_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT              i,j;
  DOUBLE     H = THREE*TEN;
  DOUBLE           x1,y1, dx,dy;
  INT       numaf;
  INT       numnp_total;
  INT       numff;
  INT      *index;
  DISCRET  *actdis;
  NODE            *actnode, *actnode2, *actfnode;
  GNODE           *actgnode;
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  FSI_DYNAMIC  *fsidyn;
  ALE_DYNAMIC  *adyn;
  ARRAY_POSITION *ipos;

  index = work->index_a.a.iv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  actpart     = &(partition[numaf]);

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  numff        = genprob.numff;
  numnp_total = actfield->dis[disnum_calc].numnp;
  actdis       = &(actfield->dis[disnum_calc]);
  ipos   = &(actfield->dis[disnum_calc].ipos);

  if (par.myrank==0)
    printf("Solving ALE (LAS)...\n\n");

  for (i=0;i<numnp_total;i++)
  {
    actnode = &(actdis->node[i]);
    x1 = actnode->x[0];
    y1 = actnode->x[1];
    j = index[i];
    actnode2 =  &(actdis->node[j]);
    actgnode = actnode2->gnode;
    actfnode = actgnode->mfcpnode[numff];
    dx = actfnode->xfs[0]-actfnode->x[0];
    dy = actfnode->xfs[1]-actfnode->x[1];
    actnode->sol_mf.a.da[ipos->mf_dispnp][0] = (y1/H)*dx;
    actnode->sol_mf.a.da[ipos->mf_dispnp][1] = (y1/H)*dy;
  }
}


void fsi_ale_LAS_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT       numaf;
  INT       numnp_total;
  INT       numff;
  INT      *index;
  DISCRET  *actdis;
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  FSI_DYNAMIC  *fsidyn;
  ALE_DYNAMIC  *adyn;
  ARRAY_POSITION *ipos;

  index = work->index_a.a.iv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  actpart     = &(partition[numaf]);
  ipos   = &(actfield->dis[disnum_calc].ipos);

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  numff        = genprob.numff;
  numnp_total = actfield->dis[disnum_calc].numnp;
  actdis       = &(actfield->dis[disnum_calc]);

/*------------------------------------ for iterative staggared schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
  if (fsidyn->ifsi==fsi_iter_stagg_fixed_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc ||
      fsidyn->ifsi==fsi_iter_stagg_CHEB_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc_force ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_FD ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_I ||
      fsidyn->ifsi==fsi_coupling_freesurface)
  {
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_mf,
		     node_array_sol_mf,
		     ipos->mf_dispn,
		     ipos->mf_dispnm);
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_mf,
		     node_array_sol_mf,
		     ipos->mf_dispnp,
		     ipos->mf_dispn);
  }

/*--------------------- to get the corrected free surface position copy
  --------------------------------- from sol_mf[1][j] to sol[0][j] */
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_mf,
		   node_array_sol,
		   ipos->mf_dispn,
		   0);

/*---------------------------------------------- increment output flags */
  work->outstep++;
  work->restartstep++;

  if (work->outstep==adyn->updevry_disp && ioflags.ale_disp==1 && ioflags.output_out==1)
  {
    work->outstep=0;
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }

/*---------------------------------------------------------- monitoring */
  if (ioflags.monitor==1)
    monitoring(actfield,disnum_calc,numaf,0,adyn->time);


/*------------------------------------------------- write restart data */
  if (work->restartstep==fsidyn->uprestart)
  {
    work->restartstep=0;
#ifdef BINIO
    if(disnum_io != disnum_calc)
      restart_write_bin_aledyn(&work->restart_context, adyn);
    else
      restart_write_bin_aledyn(&work->out_context, adyn);
#else
    restart_write_aledyn(adyn,actfield,actpart,actintra);
#endif
  }
}


void fsi_ale_LAS_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  )
{
  dserror("not supported");
}


void fsi_ale_LAS_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
#ifdef BINIO
  if (ioflags.output_bin)
    if (ioflags.ale_disp==1)
    {
      ALE_DYNAMIC  *adyn;
      adyn   = alldyn[genprob.numaf].adyn;
      out_results(&work->out_context, adyn->time, adyn->step, 0, OUTPUT_DISPLACEMENT);
    }
#endif
}


void fsi_ale_LAS_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT       numaf;
  INT       numnp_total;
  INT       numff;
  INT      *index;
  DISCRET  *actdis;
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  FSI_DYNAMIC  *fsidyn;
  ALE_DYNAMIC  *adyn;

  index = work->index_a.a.iv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  actpart     = &(partition[numaf]);

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  numff        = genprob.numff;
  numnp_total = actfield->dis[disnum_calc].numnp;
  actdis       = &(actfield->dis[disnum_calc]);

  /*------------------------------------------- print out results to .out */
  if (work->outstep!=0 && ioflags.ale_disp==1 && ioflags.output_out==1)
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);

  /*------------------------------------------------------------- tidy up */
  amdel(&work->index_a);

#ifdef BINIO
  destroy_bin_out_field(&work->out_context);
  if(disnum_io != disnum_calc)
    destroy_bin_out_field(&work->restart_context);
#endif
}


#if 0
/*!----------------------------------------------------------------------
\brief   interpolating mesh displacements for LAS

<pre>                                                          genk 02/03


</pre>

\param  *fsidyn     FSI_DYNAMIC    (i)
\param  *adyn       STRUCT_DYNAMIC (i)
\param  *actfield   FIELD          (i)     ale field
\param   mctrl      INT            (i)     control flag
\warning
\return void

\sa   calling: calelm(), monitoring(), ale_setdirich(),
      called by: fsi_ale()

*----------------------------------------------------------------------*/
void fsi_ale_LAS(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  INT                mctrl
  )
{
#ifdef DEBUG
dstrc_enter("fsi_ale_LAS");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1:
  fsi_ale_LAS_setup(work,actfield,disnum_calc,disnum_io);
  break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...0][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        *
 *======================================================================*/
case 2:
  fsi_ale_LAS_calc(work,actfield,disnum_calc,disnum_io);
  if (alldyn[genprob.numaf+1].fsidyn->ifsi>=fsi_iter_stagg_fixed_rel_param ||
      alldyn[genprob.numaf+1].fsidyn->ifsi<fsi_coupling_undefined)
    break;

/*======================================================================*
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
  fsi_ale_LAS_final(work,actfield,disnum_calc,disnum_io);
  break;


/*======================================================================*
                            Binary Output
 *======================================================================*/
case 98:
  fsi_ale_LAS_output(work,actfield,disnum_calc,disnum_io);
  break;


/*======================================================================*
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:
  fsi_ale_LAS_cleanup(work,actfield,disnum_calc,disnum_io);
  break;
default:
   dserror("Parameter out of range: mctrl \n");
} /* end switch (mctrl) */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fsi_ale_lin */
#endif

#endif
/*! @} (documentation module close)*/
