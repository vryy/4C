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
FILE           *out = allfiles.out_out;

DOUBLE          t2,tt;

#ifdef D_MORTAR
INTERFACES            *int_faces; /* interface information for mortar   */
#endif

#ifdef DEBUG
dstrc_enter("dyn_fsi");
#endif

#ifdef PERF
perf_begin(40);
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

#ifdef PERF
perf_begin(41);
#endif

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
#if defined(BINIO)
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

/*
 * Binary output has to be done by the algorithms because the
 * contexts are there. */
mctrl = 98;
fsi_ale(alefield,mctrl);
fsi_fluid(fluidfield,mctrl);
fsi_struct(structfield,mctrl,itnum);

#ifdef D_MORTAR
if (fsidyn->coupmethod == 0) /* mortar method */
{
  /*----------------------------- fill structure interfaces with values */
  fsi_init_interfaces(structfield,fluidfield,int_faces);
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
    case 1:
      fprintf(out,"Basic Sequential Staggered Scheme\n");
      break;
    case 2:
      fprintf(out,"Sequential Staggered Scheme with Predictor\n");
      break;
    case 4:
      fprintf(out,"Iterative Staggered Scheme with Fixed Relaxation Parameter\n");
      break;
    case 5:
      fprintf(out,"Iterative Staggered Scheme with Relaxation Parameter via Aitken Iteration\n");
      break;
    case 6:
      fprintf(out,"Iterative Staggered Scheme with Relaxation Parameter via Steepest Descent Method\n");
      break;
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
/*------------------------------------------------ output to the screen */
if (par.myrank==0) fsi_algoout(itnum);

if (par.myrank==0)
  fprintf(out,"%5d | %10.3f | %3d |",fsidyn->step,fsidyn->time,itnum);


/*----------------------------------- basic sequential staggered scheme */
if (fsidyn->ifsi==1 || fsidyn->ifsi==3)
{

   dsassert(fsidyn->ifsi!=3,"Scheme with DT/2-shift not implemented yet!\n");



   /*------------------------------- CFD -------------------------------*/
#ifdef PERF
perf_begin(42);
#endif

   fsi_fluid(fluidfield,mctrl);

#ifdef PERF
perf_end(42);
#endif



   /*------------------------------- CSD -------------------------------*/
#ifdef PERF
perf_begin(43);
#endif

   fsi_struct(structfield,mctrl,itnum);

#ifdef PERF
perf_end(43);
#endif




   /*------------------------------- CMD -------------------------------*/
#ifdef PERF
perf_begin(44);
#endif

   fsi_ale(alefield,mctrl);

#ifdef PERF
perf_end(44);
#endif

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

#ifdef PERF
perf_begin(43);
#endif

      fsi_struct(structfield,mctrlpre,itnum);

#ifdef PERF
perf_end(43);
#endif

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
#ifdef PERF
perf_begin(44);
#endif

   fsi_ale(alefield,mctrl);

#ifdef PERF
perf_begin(44);
#endif



   /*------------------------------- CFD -------------------------------*/
#ifdef PERF
perf_begin(42);
#endif

   fsi_fluid(fluidfield,mctrl);

#ifdef PERF
perf_end(42);
#endif


   #ifdef D_MORTAR
   if (fsidyn->coupmethod == 0) /* mortar method */
   {
     fsi_calc_intforces(int_faces);
     /* --------put coupling forces from fluid nodes to structure nodes */
     fsi_put_coupforc2struct(structfield, int_faces);
   }
   #endif




   /*------------------------------- CSD -------------------------------*/
#ifdef PERF
perf_begin(43);
#endif

   fsi_struct(structfield,mctrl,itnum);

#ifdef PERF
perf_end(43);
#endif

}
/*--------------------------------------------- strong coupling schemes */
if (fsidyn->ifsi>=4)
{
   /*-------------------------------------- iteration convergence check */
   converged=fsi_convcheck(structfield, itnum);

   if (converged==0) /*--------------------------------- no convergence */
   {


   /*----------------------------- compute optimal relaxation parameter */

#ifdef PERF
perf_begin(45);
#endif

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


      if (par.myrank==0)
        fprintf(out,"   %7.5f  |",fsidyn->relax);


      /*-------------- relaxation of structural interface displacements */
      fsi_relax_intdisp(structfield);

#ifdef PERF
perf_end(45);
#endif

      itnum++;
      if (par.myrank==0)
        fprintf(out,"            |\n");
      fflush(out);
      goto fielditer;
   }
   else /*------------------------------------------------- convergence */
   {
     if (par.myrank==0)
       fprintf(out,"            |");


      mctrl=3;


      /*--------------------- update MESH data -------------------------*/
#ifdef PERF
perf_begin(44);
#endif

      fsi_ale(alefield,mctrl);

#ifdef PERF
perf_end(44);
#endif




      /*-------------------- update FLUID data -------------------------*/
#ifdef PERF
perf_begin(42);
#endif

      fsi_fluid(fluidfield,mctrl);

#ifdef PERF
perf_end(42);
#endif




      /*------------------ update STRUCTURE data -----------------------*/
#ifdef PERF
perf_begin(43);
#endif

      fsi_struct(structfield,mctrl,itnum);

#ifdef PERF
perf_end(43);
#endif

   }
}
else
{
  if (par.myrank==0)
    fprintf(out,"            |            |");
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

if (resstep==fsidyn->upres)
{
   resstep=0;

   /* print out solution to GiD */
   if (ioflags.output_gid==1 && par.myrank==0) {
     out_checkfilesize(1);
     out_gid_sol_fsi(fluidfield,structfield);
   }


#ifdef BINIO
   /*
    * Binary output has to be done by the algorithms because the
    * contexts are there. */
   if (ioflags.output_bin==1)
   {
     mctrl = 98;
     fsi_ale(alefield,mctrl);
     fsi_fluid(fluidfield,mctrl);
     fsi_struct(structfield,mctrl,itnum);
   }
#endif

}

/*---------------------------------------------- write fsi-restart data */
if (restartstep==fsidyn->uprestart)
{
   restartstep=0;
#ifdef BINIO
   restart_write_bin_fsidyn(fsidyn);
#else
   restart_write_fsidyn(fsidyn);
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




#ifdef PERF
perf_end(40);
#endif



#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dyn_fsi */


/*! @} (documentation module close)*/
