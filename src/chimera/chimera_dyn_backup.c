#ifdef D_CHIMERA
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "chimera_prototypes.h"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA                      *alldyn;
/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01   |
  | general problem data                                                |
  | global variable GENPROB genprob is defined in global_control.c      |
  *---------------------------------------------------------------------*/
extern struct _GENPROB               genprob;
/*----------------------------------------------------------------------*
  |                                                       m.gee 02/02   |
  | number of load curves numcurve                                      |
  | vector of structures of curves                                      |
  | defined in input_curves.c                                           |
  | INT                   numcurve;                                     |
  | struct _CURVE      *curve;                                          |
  *---------------------------------------------------------------------*/
extern INT                           numcurve;
extern struct _CURVE                *curve;
/*-----------------------------------------------------------------------*
  | enum _CALC_ACTION                                      m.gee 1/02    |
  | command passed from control routine to the element level             |
  | to tell element routines what to do                                  |
  | defined globally in global_calelm.c                                  |
  *----------------------------------------------------------------------*/
extern enum _CALC_ACTION             calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD                *field;
/*----------------------------------------------------------------------*
  | global variable *solv, vector of lenght numfld of structures SOLVAR  |
  | defined in solver_control.c                                          |
  |                                                                      |
  |                                                       m.gee 11/00    |
  *----------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR                   par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS              ioflags;
extern struct _PARTITION            *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR               *solv;
/*extern struct _CHIMERA_DATA         *chm_data;*/





void chimera_dyn()
{
  INT                i;
  INT                dummy = 0;
  FIELD             *fluidfield;
  FLUID_DYNAMIC     *fdyn;
  CALC_ACTION       *action;

#ifdef DEBUG
  dstrc_enter("chimera_dyn");
#endif
/*----------------------------------------------------------------------*/

  /* plausibility check */
  dsassert(genprob.numfld==1,"**ERROR** Number of fields is to be ONE!\n");
  dsassert(genprob.numff==0,"**ERROR** Field number is to be ZERO!\n");

  /* access to the fluid field */
  fluidfield = &(field[genprob.numff]);
  fdyn = alldyn[genprob.numff].fdyn;
  action  = &(calc_action[genprob.numff]);

  /* init all applied time curves */
  for (i=0; i<numcurve; i++)
    dyn_init_curve(i,fdyn->nstep,fdyn->dt,fdyn->maxtime);

  /* initialize the fluid problem on overlapping domains */
  chimera_initialize();
#if 0
  if (ioflags.chimera_to_matlab==1)
  {
    chimera_to_matlab();
  }
#endif

/************************************************************************/
/*----------------------------------------------- start of time process */
/************************************************************************/
  timeloop:
    fdyn->step++;
    fdyn->acttime += fdyn->dt;

    /* solve fluid problem on each discretization */
    for (i=0; i<fluidfield->ndis; i++)               /* => i == numdis  */
    {
      /* solve time step till convergence */
      chimera_solve(i);
      /* finalize time step */
      chimera_finalize(i);
    }

    /* check time step */
    if (fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime)
      goto timeloop;
/************************************************************************/
/*------------------------------------------------- end of time process */
/************************************************************************/

    /* clean up */
    chimera_clean();

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of chimera_dyn */
#endif /* D_CHIMERA */
