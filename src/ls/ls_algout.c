#ifdef D_LS
#include "../headers/standardtypes.h" 
#include "ls_prototypes.h" 





/************************************************************************
 ----------------------------------------- last checked by Irhan 26.04.04
 ************************************************************************/
void ls_algout(       
  LS_DYNAMIC	*lsdyn
  )
{
#ifdef DEBUG 
  dstrc_enter("ls_algoout");
#endif
/*----------------------------------------------------------------------*/  
  
  switch(lsdyn->iop)
  {
      case 0:
        printf("\nTIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
               lsdyn->time,lsdyn->maxtime,lsdyn->dt,lsdyn->step,lsdyn->nstep);
        break;
      default:
        dserror("unknown integration scheme: iop\n");
  }
  
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
  dstrc_exit();
#endif
  
  return;
} /* end of ls_algoout*/
#endif
