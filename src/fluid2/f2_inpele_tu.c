/*!----------------------------------------------------------------------
\file
\brief read fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU 
#include "../headers/standardtypes.h" 
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
#include "fluid2.h"
/*!--------------------------------------------------------------------- 
\brief read fluid2 element from input-file

<pre>                                                        he  11/02
</pre>
\param  *ele0	   ELEMENT	   (i)	 
\param  *ele1	   ELEMENT	   (o)
\return void                                                                       
------------------------------------------------------------------------*/
void f2tu_dis(
    ELEMENT *ele0,
    ELEMENT *ele1, 
    INT numele,
    INT nnode)
{
  INT        i;             /* simply a counter                           */

#ifdef DEBUG 
  dstrc_enter("f2tu_dis");
#endif

  /*------------------------------------------------ allocate the element */      
  ele1->e.f2_tu = (FLUID2_TU*)CCACALLOC(1,sizeof(FLUID2_TU));                         
  if (ele1->e.f2_tu==NULL) dserror("Allocation of element FLUID2 failed\n");
  /*---------------------------------------------- read the element nodes */
  ele1->Id = ele0->Id + numele;
  ele1->mat=ele0->mat;
  ele1->e.f2_tu->is_ale=ele0->e.f2->is_ale;
  ele1->numnp=ele0->numnp;
  ele1->distyp=ele0->distyp;

  ele1->lm = (INT*)CCACALLOC(ele1->numnp,sizeof(INT));
  if (ele1->lm==NULL) dserror("Allocation of lm in ELEMENT failed\n");

  for (i=0; i<ele1->numnp; i++)
  {
    ele1->lm[i] = ele0->lm[i] + nnode;
  } /* end of loop over all nodes */	

  ele1->e.f2_tu->nGP[0]=ele0->e.f2->nGP[0];
  ele1->e.f2_tu->nGP[1]=ele0->e.f2->nGP[1]; 

#ifdef DEBUG 
  dstrc_exit();
#endif

  return;
} /* end of f2tu_dis */

#endif
