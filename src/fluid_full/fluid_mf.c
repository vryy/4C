/*!----------------------------------------------------------------------
\file
\brief fluid multifield (freesurface) algorithm

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "../fsi_full/fsi_prototypes.h"
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
/*!---------------------------------------------------------------------
\brief routine to control fluid with freesurface

<pre>                                                         genk 09/02

this function controls the algorithm for multifield fluid problems,
especially fluid problems with free surfaces

</pre>

\param  mctrl   INT  (i)   control flag
\return void

------------------------------------------------------------------------*/
void fluid_mf(INT mctrl)
{

#ifdef D_FSI

static INT      numfld;       /* number of fiels                        */
static INT      numff;
static INT      numaf;        /* actual number of fields                */
INT             resstep=0;    /* counter to control output              */
INT             actcurve;     /* actual curve                           */
static FIELD          *fluidfield;
static FIELD          *alefield;
static FLUID_DYNAMIC  *fdyn;
static ALE_DYNAMIC    *adyn;
static FSI_DYNAMIC    *fsidyn;

static FSI_FLUID_WORK  fluid_work;
static FSI_ALE_WORK    ale_work;

INT                    disnumf = 0;
INT                    disnums = 0;
INT                    disnuma = 0;

#ifdef DEBUG
dstrc_enter("fluid_mf");
#endif

/*----------------------------------------------------------------------*/
if (mctrl==99) goto cleaningup;

/*--------------------------------------------------------- check input *
 | convention used at the moment:                                       |
 | FIELD 0: fluid                                                       |
 | FIELD 1: mesh / ale                                                  |
 *----------------------------------------------------------------------*/
numfld=genprob.numfld;
numff =genprob.numff;
numaf =genprob.numaf;

dsassert(numfld==2,"TWO fields needed for FSI-problem!\n");
fluidfield  = &(field[numff]);
alefield    = &(field[numaf]);

/*-------------------------------------------------- plausibility check */
dsassert(fluidfield->fieldtyp==fluid,"FIELD 0 has to be fluid\n");
dsassert(alefield->fieldtyp==ale,"FIELD 1 has to be ale\n");

/*======================================================================*
                    I N I T I A L I S A T I O N
 *======================================================================*/
mctrl=1;
fdyn= alldyn[numff].fdyn;
adyn= alldyn[numaf].adyn;
fsidyn=alldyn[numaf+1].fsidyn;
fsidyn->time=ZERO;
fsidyn->step=0;

fsidyn->ichecke=0;
fsidyn->ifsi=fsi_coupling_freesurface;

fdyn->dt=fsidyn->dt;
adyn->dt=fsidyn->dt;

/*--------------------- initialise fluid multifield coupling conditions */
fluid_initmfcoupling(fluidfield,alefield);

/*---------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,fsidyn->nstep,fsidyn->dt,fsidyn->maxtime);

/*------------------------------------------------------ initialise ale */
fsi_ale_setup(&ale_work,alefield,disnuma,disnuma);
/*---------------------------------------------------- initialise fluid */
fsi_fluid_setup(&fluid_work,fluidfield,disnumf,disnumf);

if (genprob.restart>0)
{
   if (fdyn->acttime != adyn->time)
   dserror("Restart problem: Time not identical in fields!\n");
   if (fdyn->step != adyn->step)
   dserror("Restart problem: Step not identical in fields!\n");
   fsidyn->time = fdyn->acttime;
   fsidyn->step = fdyn->step;
}

/*----------------------------------------------------------------------*/
if (par.myrank==0) out_gid_msh();
/*--------------------------------------- write initial solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
if (par.myrank==0)
  out_gid_sol_fsi(fluidfield,NULL,disnumf, disnums);

/*======================================================================*
                              T I M E L O O P
 *======================================================================*/
timeloop:
mctrl=2;
fsidyn->step++;
fsidyn->time += fsidyn->dt;
fdyn->step=fsidyn->step;
adyn->step=fsidyn->step;
fdyn->acttime = fsidyn->time;
adyn->time = fsidyn->time;
/*------------------------------------------------ output to the screen */
if (par.myrank==0)
{
printf("TIME: %11.4E/%11.4E   DT = %11.4E   STEP = %4d/%4d \n",
          fsidyn->time,fsidyn->maxtime,fsidyn->dt,fsidyn->step,fsidyn->nstep);
printf("\n");
}

/*------------------------------- CMD ----------------------------------*/
fsi_ale_calc(&ale_work,alefield,disnuma,disnuma,NULL,0);
/*------------------------------- CFD ----------------------------------*/
fsi_fluid_calc(&fluid_work,fluidfield,disnumf,disnumf,alefield,disnuma);
fsi_fluid_final(&fluid_work,fluidfield,disnumf,disnumf);
/*------------------------------- CMD ----------------------------------*/
mctrl=3; /*------------------------------------ finalise ALE - timestep */
fsi_ale_final(&ale_work,alefield,disnuma,disnuma);

/*--------------------------------------- write current solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
resstep++;
if (resstep==fsidyn->upres && par.myrank==0)
{
   resstep=0;
   out_checkfilesize(1);
   out_gid_sol_fsi(fluidfield,NULL,disnumf, disnums);
}

/*------------------------------------------- finalising this time step */
if (fsidyn->step < fsidyn->nstep && fsidyn->time <= fsidyn->maxtime)
   goto timeloop;

/*======================================================================*
                   C L E A N I N G   U P   P H A S E
 *======================================================================*/
cleaningup:
mctrl=99;
fsi_fluid_cleanup(&fluid_work,fluidfield,disnumf,disnumf);
fsi_ale_cleanup(&ale_work,alefield,disnuma,disnuma);

#else
dserror("FSI-functions not compiled in!\n");
#endif

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fluid_mf */


#endif
/*! @} (documentation module close)*/
