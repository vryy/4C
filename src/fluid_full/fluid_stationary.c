/*!----------------------------------------------------------------------
\file
\brief stationary solution algorithm for fluid

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
/*----------------------------------------------------------------------*
 | routine to control stationary algorithm for fluid for fluid problems |
 | combined with Newton, fixed point iteration and fixed point like     |
 | schemes.                                                 genk  03/02 |
 *----------------------------------------------------------------------*/
void fluid_stat()
{

#ifdef DEBUG 
dstrc_enter("fluid_stat");
#endif

dserror("algorithm for stationary fluid problems not implemented yet!\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_stat */ 



#endif
/*! @} (documentation module close)*/
