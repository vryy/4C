/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_call_matgeononl' which selects the proper
       material law for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | select proper material law                               ah 06/02    |
 *----------------------------------------------------------------------*/
void w1_call_matgeononl(MATERIAL  *mat,
                        WALL_TYPE  wtype,
                        DOUBLE    *strain,
                        DOUBLE   **stress,
                        DOUBLE   **d,
                        INT        numeps)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_call_matgeononl");
#endif
/*--------------------------- call material law -> get tangent modulus--*/
  switch(mat->mattyp)
  {
  case m_stvenant:/*--------------------------------- linear elastic ---*/
    w1_mat_linelgeonon(mat->m.stvenant->youngs,
                       mat->m.stvenant->possionratio,
                       wtype,
                       strain,
                       d,
                       stress,
                       numeps);
  break;
  default:
    dserror(" unknown type of material law");
  break;
  }
/*---------------------------------------------- evaluate stress (2.PK)--*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_call_matgeononl */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/

