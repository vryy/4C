/*!----------------------------------------------------------------------
\file
\brief ls_algout.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/


extern struct _LS_DYNAMIC         *lsdyn;        /* pointer to dynamic control struct   */


/*!----------------------------------------------------------------------
\brief printing subroutine for the time integration scheme used

<pre>                                                            irhan 05/0
printing subroutine for the time integration scheme used
</pre>

*----------------------------------------------------------------------*/
void ls_algout()
{
#ifdef DEBUG
  dstrc_enter("ls_algoout");
#endif
/*----------------------------------------------------------------------*/

  switch(lsdyn->iop)
  {
      case 0:
        printf("\nTIME: %11.4E/%11.4E  DT = %11.4E  One-Step-Theta  STEP = %4d/%4d \n",
               lsdyn->acttime,lsdyn->maxtime,lsdyn->dt,lsdyn->step,lsdyn->nstep);
        break;
      default:
        dserror("unknown integration scheme: iop\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_algout*/
/*! @} (documentation module close)*/
#endif
