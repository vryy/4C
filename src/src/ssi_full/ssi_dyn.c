/*!----------------------------------------------------------------------
\file
\brief control function for SSI

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup SSI
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../output/output_prototypes.h"
#include "ssi_prototypes.h"
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
\brief routine to control ssi dynamic analyis

<pre>                                                         genk 10/03

Nonlinear dynamical algorithms for ssi-problems:
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
void dyn_ssi()
{

static INT             numfld;       /* number of fiels                 */
static INT             resstep=0;    /* counter for output control      */
static INT             restartstep=0;/* counter for restart control     */
static INT             mctrl;
INT                    actcurve;     /* actual curve                    */
INT                    itnum=0;      /* iteration counter               */
INT                    converged;    /* convergence flag                */

static FIELD          *masterfield;
static FIELD          *slavefield;
static SSI_DYNAMIC    *ssidyn;
static STRUCT_DYNAMIC *slave_sdyn;
static STRUCT_DYNAMIC *master_sdyn;

ARRAY                 slv_nods_onint; /* vector of slv nods on int */
ARRAY                 mst_nods_onint; /* vector of mst nods on int */
INT                   **slv_nods_on_int;
INT                   **mst_nods_on_int;
DENSE                 *lhs_dens;     /* pointer to dense structure */
DENSE                 *lhs_dens1;    /* pointer to dense structure */
INTERFACES            *int_faces;

INT                    disnumm = 0;
INT                    disnums = 0;


#ifdef DEBUG
dstrc_enter("dyn_ssi");
#endif

/*----------------------------------------------------------------------*/
#ifdef D_SSI
numfld=genprob.numfld;
dsassert(numfld==2,"Only two fields allowd for SSI-problem!\n");
masterfield = &(field[0]);
slavefield  = &(field[1]);

slv_nods_on_int = amdef("slave_nods_int",&slv_nods_onint, 1, 1, "IV");
mst_nods_on_int = amdef("master_nods_int",&mst_nods_onint, 1, 1, "IV");
lhs_dens = (DENSE*)CCACALLOC(1,sizeof(DENSE));
lhs_dens1 = (DENSE*)CCACALLOC(1,sizeof(DENSE));
int_faces = (INTERFACES*)CCACALLOC(1,sizeof(INTERFACES));

/*======================================================================*
                    I N I T I A L I S A T I O N
 *======================================================================*/
mctrl=1;
slave_sdyn    = alldyn[0].sdyn;
master_sdyn   = alldyn[1].sdyn;
ssidyn        = alldyn[2].ssidyn;
ssidyn->time  = ZERO;
ssidyn->step  = 0;

/*---------------------------------- initialise ssi coupling conditions */
ssi_initcoupling(masterfield,slavefield, int_faces);

/* --------------------------allocate memory for the several interfaces */
int_faces->interface = (INTERFACE*)CCACALLOC(int_faces->numint,
                       sizeof(INTERFACE));
/*--------------------------------- determine structural interface dofs */
ssi_master_intdofs(masterfield,ssidyn);

/*---------------------------------------- init all applied time curves */
for (actcurve = 0;actcurve<numcurve;actcurve++)
   dyn_init_curve(actcurve,ssidyn->nstep,ssidyn->dt,ssidyn->maxtime);

/*--------------------------------------------- initialise master field */
ssi_struct(ssidyn,master_sdyn,masterfield,disnumm,mctrl,0,ssi_master);
/*---------------------------------------------- initialise slave field */
ssi_struct(ssidyn,slave_sdyn,slavefield,disnums,mctrl,0,ssi_slave);

/*----------------------------------------- initialise AITKEN iteration */
if (ssidyn->ifsi==5)
{
   ssi_aitken(masterfield,ssidyn,itnum,0);
}

/*----------------------------------------------------------------------*/
if (par.myrank==0) out_gid_msh();
/*--------------------------------------- write initial solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
if (par.myrank==0) out_gid_sol_ssi(slavefield,masterfield,disnums,disnumm);

/* computation of the number of nodes on the slave and master interface */
/* ------------------------------only for nonconforming discretizations */

if (ssidyn->conformmesh == 1)
{
  ssi_detect_nmb_nods_int(masterfield, slavefield, int_faces);
}


/*======================================================================*
                              T I M E L O O P
 *======================================================================*/
/*----------------------------------------------------------------------*/
/* structure of nodal solutions in array sol_mf.a.da[0..6][0..1]        */
/* for non-conforming discretization of ssi-problems                    */
/*             mrt_coeff     frc2mst   mst_field   disp4slv   slv_field */
/*  [0][..]                                                             */
/*  mst                                                                 */
/*    free                              d(n+1)                 d(n+1)   */
/*    Dbc                              d(n+0.5)               d(n+0.5)  */
/*  slv                                                                 */
/*    free                              d(n+1)                 d(n+1)   */
/*    Dbc                              d(n+0.5)               d(n+0.5)  */
/*  [1][..]                                                             */
/*  mst                                                                 */
/*    free                               d(n)                   d(n)    */
/*    Dbc                              d(n-0.5)               d(n-0.5)  */
/*  slv                                                                 */
/*    free                               d(n)                   d(n)    */
/*    Dbc                              d(n-0.5)               d(n-0.5)  */
/*  [2][..]                                                             */
/*  [3][..]                                                             */
/*  [4][..]                                                             */
/*  mst                    f_ext(n+1)                                   */
/*  slv                                                      f_int(n+1) */
/*  [5][..]                                                             */
/*  [6][..]                                                             */
/*  mst       d_relax                                                   */
/*  slv                                            d_relax              */
/*----------------------------------------------------------------------*/

timeloop:
mctrl=2;
ssidyn->step++;
ssidyn->time += ssidyn->dt;
master_sdyn->step=ssidyn->step;
slave_sdyn->step=ssidyn->step;
master_sdyn->time = ssidyn->time;
slave_sdyn->time = ssidyn->time;

/*======================================================================*
   do the iteration over the fields within one timestep
 *======================================================================*/
itnum=0;
fielditer:
/*------------------------------------------------ output to the screen */
#if 0
if (par.myrank==0) ssi_algoout(ssidyn,itnum);
#endif

dsassert(ssidyn->ifsi!=3,"Scheme with DT/2-shift not implemented yet!\n");

/* computation of mortar approximation, only for interfaces with noncon-*/
/* ---------------------------------------------forming discretizations */

if (ssidyn->conformmesh == 1) /* nonconforming discretization */
{
  if (ssidyn->coupmethod == 0) /* mortar method */
  {
    ssi_mortar_coeff(masterfield, slavefield, ssidyn->step,
                     ssidyn, int_faces);
    /* ------------put coupling forces from slave nodes to master nodes */
    ssi_put_coupforc2mst(masterfield, slavefield, int_faces);
  }
  if (ssidyn->coupmethod == 1) /* interpolation method */
  {
    /* ------------put coupling forces from slave nodes to master nodes */
    ssi_interpolation_frc(masterfield, slavefield, 4, ssidyn);
  }
}

/*######################################################################*/
/*------------------------- MASTER FIELD -------------------------------*/
/*######################################################################*/
ssi_struct(ssidyn,master_sdyn,masterfield,disnumm,mctrl,itnum,ssi_master);

if ((ssidyn->ifsi==5)&&(itnum > 0))
{
  /* compute relaxation parameter with AITKENS DELTA METHOD */
  ssi_aitken(masterfield,ssidyn,itnum,1);
}

/* ----compute coupling displacements for non-conforming discretization */

if (ssidyn->conformmesh == 1)
{
  if(ssidyn->coupmethod == 0) /* mortar method */
  {
    /* -put displacements from the master nodes at the interface to the */
    /* -----------slave nodes at the interface, relaxation of interface */
    /* --------------------------------------displacements is done here */
    ssi_calc_disp4slv(masterfield, slavefield, ssidyn, int_faces);
  }
  if (ssidyn->coupmethod == 1) /* interpolation method */
  {
    /* -----put coupling displacements from master nodes to slave nodes */
    /* -------------relaxation of interface  displacements is done here */
    ssi_interpolation_disp(slavefield, masterfield, 6, ssidyn);
  }
}
/*######################################################################*/
/*------------------------- SLAVE FIELD --------------------------------*/
/*######################################################################*/
ssi_struct(ssidyn,master_sdyn,slavefield,disnums,mctrl,itnum,ssi_slave);

if (ssidyn->ifsi<4)
{
   mctrl=3;
   /*------------------------ MASTER FIELD -----------------------------*/
   ssi_struct(ssidyn,master_sdyn,masterfield,disnumm,mctrl,itnum,ssi_master);
   /*------------------------- SLAVE FIELD -----------------------------*/
   ssi_struct(ssidyn,master_sdyn,slavefield,disnums,mctrl,itnum,ssi_slave);
}

/*--------------------------------------------- strong coupling schemes */
if (ssidyn->ifsi>=4)
{
   /*-------------------------------------- iteration convergence check */
   converged=ssi_convcheck(masterfield, ssidyn, itnum);

   if (converged==0) /*--------------------------------- no convergence */
   {
   /*----------------------------- compute optimal relaxation parameter */
      if (ssidyn->ifsi==5)
      {
         /*ssi_aitken(masterfield,ssidyn,itnum,1);*/
      }
      else if (ssidyn->ifsi==7)
      {
         dserror("RELAX via CHEBYCHEV not implemented yet!\n");
      }
      /*-------------- relaxation of structural interface displacements */
      /*ssi_relax_intdisp(slavefield,ssidyn);*/
      itnum++;
      goto fielditer;
   }
   else /*------------------------------------------------- convergence */
   {
      mctrl=3;
      /*------------------------- MASTER FIELD -------------------------*/
      ssi_struct(ssidyn,master_sdyn,masterfield,disnumm,mctrl,itnum,ssi_master);
      /*------------------------- SLAVE FIELD --------------------------*/
      ssi_struct(ssidyn,master_sdyn,slavefield,disnums,mctrl,itnum,ssi_slave);
   }
}

/*--------------------------------------- write current solution to gid */
/*----------------------------- print out solution to 0.flavia.res file */
resstep++;
restartstep++;


if (resstep==ssidyn->upres && par.myrank==0)
{
   resstep=0;
   out_checkfilesize(1);
   out_gid_sol_ssi(slavefield,masterfield,disnums,disnumm);
}


/*------------------------------------------- finalising this time step */
if (ssidyn->step < ssidyn->nstep && ssidyn->time <= ssidyn->maxtime)
   goto timeloop;


/*======================================================================*
                  C L E A N I N G   U P   P H A S E
 *======================================================================*/
mctrl=99;
/*---------------------------- MASTER FIELD ----------------------------*/
ssi_struct(ssidyn,master_sdyn,masterfield,disnumm,mctrl,itnum,ssi_master);
/*---------------------------- SLAVE FIELD -----------------------------*/
ssi_struct(ssidyn,master_sdyn,slavefield,disnums,mctrl,itnum,ssi_slave);



/*----------------------------------------------------------------------*/
#else
dserror("SSI routines are not compiled in!\n");
#endif
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of dyn_ssi */


/*! @} (documentation module close)*/
#endif
