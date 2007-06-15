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
#ifdef D_SOLID3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "solid3.h"

/*======================================================================*/
/*!
\brief Check if an SOLID3 element has an local system

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 03/06
*/
#ifndef LOCALSYSTEMS_ST
void so3_locsys_check(PARTITION* actpart)
{
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT* actele;  /* pointer to current element */
  SOLID3* actso3;  /* pointer to current SOLID3 element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("so3_locsys_check");
#endif
  
  /*--------------------------------------------------------------------*/
  /* Check if an SOLID3 element has an local system */
  /* loop over all discretisations of partition structural field */
  for (jdis=0; jdis<actpart->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<actpart->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = actpart->pdis[jdis].element[iele];
      if (actele->eltyp == el_solid3)
      {
        /* set pointer to SOLID3 */
        actso3 = actele->e.so3;
        dsassert(actele->locsys==locsys_no,
                 "locsys not compiled for SOLID3 element, define LOCALSYSTEMS_ST\n");
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
