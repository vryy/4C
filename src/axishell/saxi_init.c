/*!----------------------------------------------------------------------
\file
\brief contains the routine 'saxiinit' which initializes the element

*----------------------------------------------------------------------*/
#ifdef D_AXISHELL
#include "../headers/standardtypes.h"
#include "axishell.h"
#include "axishell_prototypes.h"

/*! 
\addtogroup AXISHELL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  initialization routine for the axisymmetric shell element

<pre>                                                              mn 05/03 
This routine acts according to the action and either initializes the element
or computes the linear stiffness matrix, the stresses or the right hand side 
vector.

</pre>
\param *actpart      PARTITION   (i)   my partition

\warning There is nothing special to this routine.
\return void                                               
\sa calling:   saxiintg;
    called by: axishell();

*----------------------------------------------------------------------*/
void saxiinit(PARTITION *actpart)
{
INT          i;
ELEMENT     *actele;
SAXI_DATA      data;

for (i=0; i<actpart->pdis[0].numele; i++)
{
  actele = actpart->pdis[0].element[i];
  if (actele->eltyp != el_axishell) continue;
  /*---------------------------------------- init integration points ---*/
  saxiintg(&data);
  /*-------------------------------- allocate the space for stresses ---*/
  am4def("stress_GP",&(actele->e.saxi->stress_GP),1,5,1,0,"D3");
  am4def("stress_ND",&(actele->e.saxi->stress_ND),1,5,1,0,"D3");
}
#ifdef DEBUG 
dstrc_enter("saxiinit");
#endif
return;
} /* end of saxiinit */
/*----------------------------------------------------------------------*/
#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
