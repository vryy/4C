/*!----------------------------------------------------------------------
\file
\brief calling time algorithms (stationary/pm/isi) for fluid

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"    
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
iop=1: Predictor Multicorrector scheme
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
int iop   ;                         /* flag for time algorithm          */
int freesurf;                       /* flag for fluid problem w/ freesurface */ 
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
freesurf = fdyn->freesurf;

/*------------------------------------------------------ initialisation */
fdyn->time=0.0;
fdyn->step=0;

/*----------------------------------------------------------------------*
|  call algorithms                                                      |
 *----------------------------------------------------------------------*/
if      (iop ==0) fluid_stat();     /* stationary solution algorithm         */
else if (iop ==1) fluid_pm();       /* predictor-multicorrector algorithm    */
else if (iop >=2 && freesurf==0) 
                  fluid_isi(fdyn);  /* implicit and semi-implicit algorithms */
else if (iop == 4 && freesurf>0)
                  fluid_mf(0);      /* fluid multiefield algorithm      */
else     dserror("unknown time algorithm"); 
/*----------------------------------------------------------------------*/

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


