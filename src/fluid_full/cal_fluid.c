#include "../headers/standardtypes.h"
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
 |  routine to control fluid dynamic analyis                 genk  03/02|
 *----------------------------------------------------------------------*/
void dyn_fluid()
{
int iopfsi;                         /* flag for time algorithm */ 
FLUID_DYNAMIC *fdyn;               /* pointer to fluid dynamic input data */
#ifdef DEBUG 
dstrc_enter("dyn_fluid");
#endif
/*----------------------------------------------------------------------*/
/*--------------------------------------------------- set some pointers */
fdyn = alldyn[0].fdyn;

iopfsi = fdyn->iopfsi;
/*------------------------------------------------------ initialisation */
fdyn->time=0.0;
fdyn->step=0;
/*----------------------------------------------- get split dof numbers */
/* it seems that one does not need this routine any more!!!! 
/* fluid_splitdof(fdyn);
/*----------------------------------------------------------------------*
|  call algorithms                                                      |
 *----------------------------------------------------------------------*/
if      (iopfsi==0) fluid_stat();      /* stationary solution algorithm         */
else if (iopfsi==1) fluid_pm();        /* predictor-multicorrector algorithm    */
else if (iopfsi>=2) fluid_isi(fdyn);   /* implicit and semi-implicit algorithms */
else     dserror("unknown time algorithm"); 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of dynfluid */



