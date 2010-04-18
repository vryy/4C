/*======================================================================*/
/*!
\file
\brief Element local system
<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 12/06
*/
#ifndef CCADISCRET
#ifdef D_WALL1

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*======================================================================*/
/*!
\brief Check if a WALL1 element has an local system

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 06/07
*/
#ifndef LOCALSYSTEMS_ST
void w1_locsys_checkoff(PARTITION* actpart)
{
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT* actele;  /* pointer to current element */
  WALL1* actw1;  /* pointer to current WALL1 element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_locsys_checkoff");
#endif

  /*--------------------------------------------------------------------*/
  /* Check if an WALL1 element has an local system */
  /* loop over all discretisations of partition structural field */
  for (jdis=0; jdis<actpart->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<actpart->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = actpart->pdis[jdis].element[iele];
      if (actele->eltyp == el_wall1)
      {
        /* set pointer to WALL1 */
        actw1 = actele->e.w1;
        /* verify */
        dsassert(actele->locsys==locsys_no,
            "locsys not compiled for WALL1 element, define LOCALSYSTEMS_ST\n");
      }  /* end if */
    }  /* end for */
  }  /* end for */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}
#endif


#endif
#endif
