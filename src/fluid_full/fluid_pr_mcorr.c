/*!----------------------------------------------------------------------
\file
\brief predictor-multicorrector time integration algorithm for fluid

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"

/*----------------------------------------------------------------------*
 | routine to control predictor-multicorrector algorithm for fluid      |
 | problems                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_pm()
{
#ifdef DEBUG 
dstrc_enter("fluid_pm");
#endif

dserror("PM algorithm for fluid problems not implemented yet!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_pm */ 
#endif
