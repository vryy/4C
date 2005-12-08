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
     FIELD             *actfield,
     INT                disnum_calc,
     INT                disnum_io,
     INT                mctrl
     )
{
INT              i,j;
INT              mone=-1;
const DOUBLE     H = THREE*TEN;
DOUBLE           x1,y1,x2,y2, dx,dy;
static INT       numaf;
static INT       numnp_total;
static INT       numff;
static INT      actpos;           /* actual position in nodal solution history          */
static INT      outstep;          /* counter for output to .out                         */
static INT      pssstep;          /* counter for output to .pss                         */
static INT      restartstep;      /* counter for output of restart data                 */
static INT      *index;
static ARRAY     index_a;
static DISCRET  *actdis;
NODE            *actnode, *actnode2, *actfnode;
GNODE           *actgnode;
static PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
static INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
static ARRAY         time_a;      /* stored time                                        */
static FSI_DYNAMIC  *fsidyn;
static ALE_DYNAMIC  *adyn;

#ifdef BINIO
static BIN_OUT_FIELD out_context;
#endif


#ifdef DEBUG
dstrc_enter("fsi_ale_LAS");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1:
numaf  = genprob.numaf;
adyn   = alldyn[numaf].adyn;
fsidyn = alldyn[numaf+1].fsidyn;
actpart     = &(partition[numaf]);
#ifdef PARALLEL
actintra    = &(par.intra[numaf]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = ale;
actintra->intra_rank   = 0;
actintra->intra_nprocs   = 1;
#endif

numff        = genprob.numff;
numnp_total = actfield->dis[disnum_calc].numnp;
actdis       = &(actfield->dis[disnum_calc]);
actpos=0;
outstep=0;
pssstep=0;
restartstep=0;

/*--------------------------- allocate one vector for storing the time */
if (par.myrank==0) amdef("time",&time_a,1000,1,"DV");

/*-------------------------------------------------- define index array */
index = amdef("index",&index_a,numnp_total,1,"IV");
aminit(&index_a,&mone);

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
      dserror("something went wrong!\n");

#ifdef BINIO

/* initialize binary output
 * It's important to do this only after all the node arrays are set
 * up because their sizes are used to allocate internal memory. */
init_bin_out_field(&out_context,
                   NULL,
                   NULL,
                   actfield, actpart, actintra, disnum_io);
#endif

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
{
   out_monitor(actfield,numaf,ZERO,1);
   monitoring(actfield,disnum_calc,numaf,actpos,adyn->time);
}

/*------------------------------------------- print out results to .out */
if (ioflags.ale_disp==1 && ioflags.output_out==1)
{
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);
}


break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        *
 *======================================================================*/
case 2:

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
   actnode->sol_mf.a.da[1][0] = (y1/H)*dx;
   actnode->sol_mf.a.da[1][1] = (y1/H)*dy;
}

if (fsidyn->ifsi>=4 || fsidyn->ifsi<0)
break;

/*======================================================================*
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
/*------------------------------------ for iterative staggared schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
if (fsidyn->ifsi>=4 || fsidyn->ifsi==-1)
{
   solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,0,2);
   solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol_mf,1,0);
}

/*--------------------- to get the corrected free surface position copy
  --------------------------------- from sol_mf[1][j] to sol[actpos][j] */
solserv_sol_copy(actfield,disnum_calc,node_array_sol_mf,node_array_sol,0,actpos);

/*---------------------------------------------- increment output flags */
outstep++;
pssstep++;
restartstep++;

if (pssstep==fsidyn->uppss && ioflags.fluid_vis==1 && par.myrank==0)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   amredef(&(time_a),time_a.fdim+1000,1,"DV");
   time_a.a.dv[actpos] = adyn->time;
   actpos++;
}

if (outstep==adyn->updevry_disp && ioflags.ale_disp==1 && ioflags.output_out==1)
{
    outstep=0;
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,disnum_calc,numaf,actpos,adyn->time);


/*------------------------------------------------- write restart data */
if (restartstep==fsidyn->uprestart)
{
   restartstep=0;
#ifdef BINIO
   restart_write_bin_aledyn(&out_context, adyn);
#else
   restart_write_aledyn(adyn,actfield,actpart,actintra);
#endif
}

/*--------------------------------------------------------------------- */

break;


/*======================================================================*
                            Binary Output
 *======================================================================*/
case 98:
#ifdef BINIO
  if (ioflags.output_bin)
    if (ioflags.ale_disp==1)
      out_results(&out_context, adyn->time, adyn->step, actpos, OUTPUT_DISPLACEMENT);
#endif
  break;


/*======================================================================*
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:

if (pssstep==0) actpos--;

/*------------------------------------------- print out results to .out */
if (outstep!=0 && ioflags.ale_disp==1 && ioflags.output_out==1)
  out_sol(actfield,actpart,disnum_io,actintra,adyn->step,actpos);

/*------------------------------------------- print out result to 0.pss */
if (ioflags.fluid_vis==1 && par.myrank==0)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = adyn->time;
   }
   visual_writepss(actfield,actpos+1,&time_a);
}

/*------------------------------------------------------------- tidy up */
amdel(&index_a);

#ifdef BINIO
destroy_bin_out_field(&out_context);
#endif

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
/*! @} (documentation module close)*/
