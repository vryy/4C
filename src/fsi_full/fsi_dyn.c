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
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
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
void dyn_fsi(INT mctrl)
{

static INT      numfld;       /* number of fiels                        */
static INT      numsf;
static INT      numff;
static INT      numaf;        /* actual number of fields                */
static INT      resstep=0;    /* counter for output control             */
static INT      restartstep=0;/* counter for restart control            */
const  INT      mctrlpre=4;   /* control flag                           */
INT             actcurve;     /* actual curve                           */
INT             itnum=0;      /* iteration counter                      */
INT             converged;    /* convergence flag                       */
static FIELD          *fluidfield;
static FIELD          *structfield;
static FIELD          *alefield;
static FLUID_DYNAMIC  *fdyn;
static FSI_DYNAMIC    *fsidyn;
static STRUCT_DYNAMIC *sdyn;
static ALE_DYNAMIC    *adyn;

#ifdef D_MORTAR
INTERFACES            *int_faces; /* interface information for mortar   */
#endif

#ifdef DEBUG
dstrc_enter("dyn_fsi");
#endif

/*----------------------------------------------------------------------*/
#ifdef D_FSI
if (mctrl==99) goto cleaningup;
/*--------------------------------------------------------- check input *
 | convention used at the moment:                                       |
 | FIELD 0: structure                                                   |
 | FIELD 1: fluid                                                       |
 | FIELD 2: mesh / ale                                                  |
 *----------------------------------------------------------------------*/
numfld=genprob.numfld;
dsassert(numfld==3,"Three fields needed for FSI-problem!\n");
structfield = &(field[0]);
numsf       = 0;
fluidfield  = &(field[1]);
numff       = 1;
alefield    = &(field[2]);
numaf       = 2;

#ifdef D_MORTAR
int_faces = (INTERFACES*)CCACALLOC(1,sizeof(INTERFACES));
#endif

/*-------------------------------------------------- plausibility check */
dsassert(structfield->fieldtyp==structure,"FIELD 0 has to be structure\n");
dsassert(fluidfield->fieldtyp==fluid,"FIELD 1 has to be fluid\n");
dsassert(alefield->fieldtyp==ale,"FIELD 2 has to be ale\n");
/*======================================================================*
                    I N I T I A L I S A T I O N
 *======================================================================*/
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

/*---------------------------------- initialise fsi coupling conditions */
fsi_initcoupling(structfield,fluidfield,alefield);

#ifdef D_MORTAR
if (fsidyn->coupmethod == 0) /* mortar method */
{
  /*-------------------------------------------- intitialize interfaces */
  fsi_initcoupling_intfaces(structfield,fluidfield, int_faces);
  /* ------------------------allocate memory for the several interfaces */
  int_faces->interface = (INTERFACE*)CCACALLOC(int_faces->numint,
                         sizeof(INTERFACE));
}
#endif

/*--------------------------------- determine structural interface dofs */
fsi_struct_intdofs(structfield);
/*---------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,fsidyn->nstep,fsidyn->dt,fsidyn->maxtime);

/*------------------------------------------------------ initialise ale */
fsi_ale(alefield,mctrl);
/*---------------------------------------------------- initialise fluid */
fsi_fluid(fluidfield,mctrl);
/*------------------------------------------------ initialise structure */
fsi_struct(structfield,mctrl,itnum);

if (genprob.restart!=0)
{
#if defined(BINIO) && defined(NEW_RESTART_READ)
   restart_read_bin_fsidyn(fsidyn, genprob.restart);
#else
   restart_read_fsidyn(genprob.restart,fsidyn);
#endif
   /*----------------------------------------------- plausibility check */
   if (fsidyn->time != adyn->time ||
       fsidyn->time != fdyn->acttime ||
       fsidyn->time != sdyn->time   )
       dserror("Restart problem: Time not identical in fields!\n");
   if (fsidyn->step != fdyn->step ||
       fsidyn->step != adyn->step ||
       fsidyn->step != sdyn->step   )
       dserror("Restart problem: Step not identical in fields!\n");
}

/*----------------------------------------- initialise AITKEN iteration */
if (fsidyn->ifsi==5)
{
   fsi_aitken(structfield,itnum,0);
}

/*----------------------------------------------------------------------*/
if (par.myrank==0) out_gid_msh();
/*--------------------------------------- write initial solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
if (par.myrank==0) out_gid_sol_fsi(fluidfield,structfield);

#ifdef D_MORTAR
if (fsidyn->coupmethod == 0) /* mortar method */
{
  /*----------------------------- fill structure interfaces with values */
  fsi_init_interfaces(structfield,fluidfield,int_faces);
}
#endif

/*======================================================================*
                              T I M E L O O P
 *======================================================================*/
timeloop:
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
/*------------------------------------------------ output to the screen */
if (par.myrank==0) fsi_algoout(itnum);
/*----------------------------------- basic sequential staggered scheme */
if (fsidyn->ifsi==1 || fsidyn->ifsi==3)
{

   dsassert(fsidyn->ifsi!=3,"Scheme with DT/2-shift not implemented yet!\n");

   /*------------------------------- CFD -------------------------------*/
   fsi_fluid(fluidfield,mctrl);
   /*------------------------------- CSD -------------------------------*/
   fsi_struct(structfield,mctrl,itnum);
   /*------------------------------- CMD -------------------------------*/
   fsi_ale(alefield,mctrl);
} /* endif (fsidyn->ifsi==1 || fsidyn->ifsi==3) */

/*---------------------------------------------- schemes with predictor */
else if (fsidyn->ifsi==2 || fsidyn->ifsi>=4)
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
   if (itnum==0)
   {
      fsi_struct(structfield,mctrlpre,itnum);
   }
   #ifdef D_MORTAR
   if (fsidyn->coupmethod == 0) /* mortar method */
   {
     /*-- computation of interface displacements for ale nodes, only for*/
     /* ----------------- interfaces with nonconforming discretizations */
     fsi_calc_disp4ale(fsidyn, int_faces);
   }
   #endif

   /*------------------------------- CMD -------------------------------*/
   fsi_ale(alefield,mctrl);
   /*------------------------------- CFD -------------------------------*/
   fsi_fluid(fluidfield,mctrl);

   #ifdef D_MORTAR
   if (fsidyn->coupmethod == 0) /* mortar method */
   {
     fsi_calc_intforces(int_faces);
     /* --------put coupling forces from fluid nodes to structure nodes */
     fsi_put_coupforc2struct(structfield, int_faces);
   }
   #endif

   /*------------------------------- CSD -------------------------------*/
   fsi_struct(structfield,mctrl,itnum);
}
/*--------------------------------------------- strong coupling schemes */
if (fsidyn->ifsi>=4)
{
   /*-------------------------------------- iteration convergence check */
   converged=fsi_convcheck(structfield, itnum);

   if (converged==0) /*--------------------------------- no convergence */
   {
   /*----------------------------- compute optimal relaxation parameter */
      if (fsidyn->ifsi==5)
      {
         fsi_aitken(structfield,itnum,1);
      }
      else if (fsidyn->ifsi==6)
      {
         fsi_gradient(alefield,structfield,fluidfield,numaf,numff,numsf);
      }
      else if (fsidyn->ifsi==7)
      {
         dserror("RELAX via CHEBYCHEV not implemented yet!\n");
      }
      /*-------------- relaxation of structural interface displacements */
      fsi_relax_intdisp(structfield);
      itnum++;
      goto fielditer;
   }
   else /*------------------------------------------------- convergence */
   {
      mctrl=3;
      /*--------------------- update MESH data -------------------------*/
      fsi_ale(alefield,mctrl);
      /*-------------------- update FLUID data -------------------------*/
      fsi_fluid(fluidfield,mctrl);
      /*------------------ update STRUCTURE data -----------------------*/
      fsi_struct(structfield,mctrl,itnum);
   }
}
/*--------------------------------------- write current solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
resstep++;
restartstep++;

if (resstep==fsidyn->upres)
{
   resstep=0;
   if (par.myrank==0) {
     out_checkfilesize(1);
     out_gid_sol_fsi(fluidfield,structfield);
   }

   /*
    * Binary output has to be done by the algorithms because the
    * contexts are there. */
   mctrl = 98;
   fsi_ale(alefield,mctrl);
   fsi_fluid(fluidfield,mctrl);
   fsi_struct(structfield,mctrl,itnum);
}

/*---------------------------------------------- write fsi-restart data */
if (restartstep==fsidyn->uprestart)
{
   restartstep=0;
   restart_write_fsidyn(fsidyn);
#ifdef BINIO
   restart_write_bin_fsidyn(fsidyn);
#endif
}

/*-------------------------------------------------------- energy check */
if (fsidyn->ichecke>0) fsi_energycheck();

/*------------------------------------------- finalising this time step */
if (fsidyn->step < fsidyn->nstep && fsidyn->time <= fsidyn->maxtime)
   goto timeloop;

/*======================================================================*
                   C L E A N I N G   U P   P H A S E
 *======================================================================*/
cleaningup:
mctrl=99;
fsi_fluid(fluidfield,mctrl);
fsi_struct(structfield,mctrl,itnum);
fsi_ale(alefield,mctrl);

/*----------------------------------------------------------------------*/
#else
dserror("FSI routines are not compiled in!\n");
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dyn_fsi */


/*! @} (documentation module close)*/
