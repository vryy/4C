/*!----------------------------------------------------------------------
\file
\brief calling time algorithms (stationary/pm/isi) for fluid

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
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"
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

/*!---------------------------------------------------------------------
\brief routine to control fluid dynamic analyis

<pre>                                                         genk 03/02

In this routine the different control programs for fluid-problems are
called. This depends on the input file paremeter TIMEINTEGR,
which is stored in fdyn->iop:
iop=0: Stationary Solution
iop=1: Projection Method
iop=2: Semi-Implicit-One-Step Method
iop=3: Semi-Implicit-Two-Step Method
iop=4: One-Step-Theta Scheme
iop=5: Fractional-Step-Theta Scheme

see dissertation of W.A. WALL, chapter 4.2 'Zeitdiskretisierung'
</pre>


\return void

------------------------------------------------------------------------*/
void dyn_fluid()
{
INT dyntyp;
INT iop   ;                         /* flag for time algorithm          */
INT freesurf;                       /* flag for fluid problem w/ freesurface */
FLUID_DYNAMIC *fdyn;                /* pointer to fluid dyn. inp.data   */

#ifdef DEBUG
dstrc_enter("dyn_fluid");
#endif

/*----------------------------------------------------------------------*/
#ifdef D_FLUID
/*----------------------------------------------------------------------*/

/*--------------------------------------------------- set some pointers */
fdyn = alldyn[genprob.numff].fdyn;
iop = fdyn->iop;
dyntyp = fdyn->dyntyp;
freesurf = fdyn->freesurf;

/*------------------------------------------------------ initialisation */
if (fdyn->init!=1) fdyn->acttime=0.0;
fdyn->step=0;

/*----------------------------------------------------------------------*
|  call algorithms                                                      |
 *----------------------------------------------------------------------*/
if (dyntyp==1) fluid_pm();
else if (dyntyp==0)
{
   switch (iop)
   {
   case 0:		/* stationary solution algorithm		*/
      fluid_stat();
   break;

   case 1:		/* Generalised alpha time integration		*/
      fdyn->time_rhs = 0;
      fluid_isi();
   break;

   case 4:		/* One step Theta 				*/
      if(freesurf==0 && fdyn->adaptive==0)
      {
	/* implicit and semi-implicit algorithms 			*/
         if(fdyn->turbu == 0 || fdyn->turbu ==1) fluid_isi();
	/* implicit and semi-implicit algorithms with turbulence-model	*/
         if(fdyn->turbu == 2)                    fluid_isi_tu();
	/* implicit and semi-implicit algorithms with turbulence-model	*/
         if(fdyn->turbu == 3)                    fluid_isi_tu_1();
      }
      else
      {
	/* fluid multiefield algorithm      				*/
         if(freesurf && fdyn->adaptive==0)		fluid_mf(0);
   	/* for adaptive time stepping 					*/
         else if(freesurf==0 && fdyn->adaptive)	fluid_isi();
	/* adaptive time stepping fuer multifield 			*/
         else if(freesurf && fdyn->adaptive)
            dserror("free surface and adaptive time stepping not yet combined");
      }
   break;

   case 7:		/* 2nd order backward differencing (BDF2)	*/
      fluid_isi();
   break;

   default:
      dserror("Unknown time integration scheme");
   }		/* end switch						*/
}
else     dserror("Unknown dynamic type");

/*----------------------------------------------------------------------*/
#else
dserror("FLUID routines are not compiled in!\n");
#endif
/*----------------------------------------------------------------------*/

#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of dyn_fluid */
/*! @} (documentation module close)*/


