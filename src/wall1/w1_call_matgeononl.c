#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*----------------------------------------------------------------------*
 | select proper material law                               ah 06/02    |
 *----------------------------------------------------------------------*/
void w1_call_matgeononl(ELEMENT   *ele,
                        MATERIAL  *mat, 
                        WALL_TYPE  wtype,
                        double   **boplin,
                        double   **xjm,
                        int        ip,       
                        double    *strain,
                        double   **stress,
                        double   **d,
                        int        istore,
                        int        numeps)
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

