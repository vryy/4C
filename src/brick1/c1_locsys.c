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
#ifdef D_BRICK1

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "brick1.h"
#include "brick1_prototypes.h"

/*======================================================================*/
/*!
\brief Check if a BRICK1 element has an local system

\param  *actpart  PARTITION   (i)   pointer to current partition
\return void

\author bborn
\date 06/07
*/
#ifndef LOCALSYSTEMS_ST
void c1_locsys_checkoff(PARTITION* actpart)
{
  INT jdis;  /* discretisation loop jndex */
  INT iele;  /* element loop index */
  ELEMENT* actele;  /* pointer to current element */
  BRICK1* actc1;  /* pointer to current BRICK1 element */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("c1_locsys_checkoff");
#endif
  
  /*--------------------------------------------------------------------*/
  /* Check if an BRICK1 element has an local system */
  /* loop over all discretisations of partition structural field */
  for (jdis=0; jdis<actpart->ndis; jdis++)
  {
    /* loop over all elements of current discretisation */
    for (iele=0; iele<actpart->pdis[jdis].numele; iele++)
    {
      /* set current element */
      actele = actpart->pdis[jdis].element[iele];
      if (actele->eltyp == el_brick1)
      {
        /* set pointer to BRICK1 */
        actc1 = actele->e.c1;
        /* verify */
        dsassert(actele->locsys==locsys_no,
                 "locsys not compiled for BRICK1 element, define LOCALSYSTEMS_ST\n");
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
