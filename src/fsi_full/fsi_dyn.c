/*!----------------------------------------------------------------------
\file
\brief control function for FSI

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fsi_prototypes.h"    
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

\param mctrl    int   (i)     evaluation flag

\return void                                                                             

------------------------------------------------------------------------*/
void dyn_fsi(int mctrl)
{

static int      numfld;       /* number of fiels                        */
static int      numsf;        
static int      numff;
static int      numaf;        /* actual number of fields                */
const  int      mctrlpre=4;   /* control flag                           */
int             actcurve;     /* actual curve                           */
int             itnum=0;      /* iteration counter                      */
int             converged;    /* convergence flag                       */
static FIELD          *fluidfield;   
static FIELD          *structfield;
static FIELD          *alefield;
static FLUID_DYNAMIC  *fdyn;
static FSI_DYNAMIC    *fsidyn;
static STRUCT_DYNAMIC *sdyn;
static STRUCT_DYNAMIC *adyn;

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
adyn= alldyn[2].sdyn;
fsidyn= alldyn[3].fsidyn;
fsidyn->time=ZERO;
fsidyn->step=0;
/*---------------------------------- initialise fsi coupling conditions */
fsi_initcoupling(structfield,fluidfield,alefield);
/*--------------------------------- determine structural interface dofs */
fsi_struct_intdofs(structfield,fsidyn);
/*---------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,fsidyn->nstep,fsidyn->dt,fsidyn->maxtime);

/*---------------------------------------------------- initialise fluid */
fsi_fluid(fsidyn,fdyn,fluidfield,mctrl,numff);
/*------------------------------------------------ initialise structure */
fsi_struct(fsidyn,sdyn,structfield,mctrl,numsf,itnum);
/*------------------------------------------------------ initialise ale */
fsi_ale(fsidyn,adyn,alefield,mctrl,numaf);
 
/*----------------------------------------------------------------------*/
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
fdyn->time = fsidyn->time;
sdyn->time = fsidyn->time;
adyn->time = fsidyn->time;

/*======================================================================*
   do the iteration over the fields within one timestep
 *======================================================================*/
itnum=0;
fielditer:
/*------------------------------------------------ output to the screen */
if (par.myrank==0) fsi_algoout(fsidyn,itnum);
/*----------------------------------- basic sequential staggered scheme */
if (fsidyn->ifsi==1 || fsidyn->ifsi==3)
{

   dsassert(fsidyn->ifsi!=3,"Scheme with DT/2-shift not implemented yet!\n");

   /*------------------------------- CFD -------------------------------*/
   fsi_fluid(fsidyn,fdyn,fluidfield,mctrl,numff);
   /*------------------------------- CSD -------------------------------*/
   fsi_struct(fsidyn,sdyn,structfield,mctrl,numsf,itnum);
   /*------------------------------- CMD -------------------------------*/
   fsi_ale(fsidyn,adyn,alefield,mctrl,numaf);
} /* endif (fsidyn->ifsi==1 || fsidyn->ifsi==3) */

/*---------------------------------------------- schemes with predictor */
else if (fsidyn->ifsi==2 || fsidyn->ifsi>=4)
{
   /*------------------------------------------------- nested iteration */
   if (fsidyn->inest==1 && fsidyn->ifsi>=4)
      dserror("COUALGO not implemented yet!\n");
   else if (fsidyn->inest==2 && fsidyn->ifsi>=4)
     dserror("COUALGO not implemented yet!\n");
   /*----------------- CSD - predictor for itnum==0 --------------------*/
   if (itnum==0)
   {
      fsi_struct(fsidyn,sdyn,structfield,mctrlpre,numsf,itnum);
   }
   /*------------------------------- CMD -------------------------------*/
   fsi_ale(fsidyn,adyn,alefield,mctrl,numaf);
   /*------------------------------- CFD -------------------------------*/
   fsi_fluid(fsidyn,fdyn,fluidfield,mctrl,numff);
   /*------------------------------- CSD -------------------------------*/
   fsi_struct(fsidyn,sdyn,structfield,mctrl,numsf,itnum);
}
/*--------------------------------------------- strong coupling schemes */
if (fsidyn->ifsi>=4)
{
   /*-------------------------------------- iteration convergence check */
   converged=fsi_convcheck(structfield, fsidyn, itnum);   

   if (converged==0) /*--------------------------------- no convergence */
   {
   /*----------------------------- compute optimal relaxation parameter */
      if (fsidyn->ifsi==5)
      {
         fsi_aitken(structfield,fsidyn,itnum);   
      }
      else if (fsidyn->ifsi==6)
      {
         dserror("RELAX via steepest descent not implemented yet!\n");
      }
      else if (fsidyn->ifsi==7)
      {
         dserror("RELAX via CHEBYCHEV not implemented yet!\n");
      }   
      /*-------------- relaxation of structural interface displacements */
      fsi_relax_intdisp(structfield,fsidyn);
      itnum++;
      goto fielditer;
   }
   else /*------------------------------------------------- convergence */
   {
      mctrl=3;
      /*--------------------- update MESH data -------------------------*/
      fsi_ale(fsidyn,adyn,alefield,mctrl,numaf);
      /*-------------------- update FLUID data -------------------------*/
      fsi_fluid(fsidyn,fdyn,fluidfield,mctrl,numff);      
      /*------------------ update STRUCTURE data -----------------------*/
      fsi_struct(fsidyn,sdyn,structfield,mctrl,numsf,itnum);
   }
}

/*-------------------------------------------------------- energy check */
if (fsidyn->ichecke>0) fsi_energycheck(fsidyn);

/*------------------------------------------- finalising this time step */
if (fsidyn->step < fsidyn->nstep && fsidyn->time <= fsidyn->maxtime)
   goto timeloop;

/*======================================================================*
                   C L E A N I N G   U P   P H A S E
 *======================================================================*/
cleaningup:
mctrl=99;
fsi_fluid(fsidyn,fdyn,fluidfield,mctrl,numff);
fsi_struct(fsidyn,sdyn,structfield,mctrl,numsf,itnum);
fsi_ale(fsidyn,adyn,alefield,mctrl,numaf);   

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
