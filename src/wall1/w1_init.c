#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_calc.h"
/*----------------------------------------------------------------------*
 | initialize the element                                    al 6/01    |
 *----------------------------------------------------------------------*/
void w1init(PARTITION *actpart)
{
int          i,j,k;
ELEMENT     *actele;
NODE        *actnode;
W1_DATA      data;

#ifdef DEBUG 
dstrc_enter("w1init");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actpart->numele; i++)
{
   actele = actpart->element[i];
   /*------------------------------------------ init integration points */
   w1intg(actele,&data,0);
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1init */
 
