/*!----------------------------------------------------------------------
\file
\brief contains the routine 'wge_setdof' to define dof of the gradient
enhanced wall element

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "wallge.h"
#include "wallge_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/

/*! 
\addtogroup WALLGE 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief  define dofs of the gradient enhanced wall element     mn 05/03

<pre>                                                              
This routine .

</pre>
\param *actpart      PARTITION   (i)   my partition

\return void                                               
\sa calling:   ;

*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*/
void wge_setdof() 
{
INT i,j;
FIELD *structfield;
ELEMENT *actele;
NODE *actnode;

#ifdef DEBUG 
dstrc_enter("wge_setdof");
#endif

#ifdef D_WALLGE
/*----------------------------------------------------------------------*/
structfield = &(field[genprob.numsf]);
for (i=0;i<structfield->dis[0].numele;i++)
{
  actele = &(structfield->dis[0].element[i]);   
  for (j=0;j<4;j++)
  {
     actnode=actele->node[j];
     actnode->numdf=3;
  }
  for (j=4;j<actele->numnp;j++)
  {
     actnode=actele->node[j];
     actnode->numdf=2;
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#endif /*D_WALLGE*/
return; 
} /* end of wge_setdof */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/
