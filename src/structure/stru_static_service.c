/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |  routine to find the node and dof to control in field  m.gee 11/01   |
 |                                                                      |
 | actfield                             the physical field to search in |
 | control_node_global                               global node number |
 | control_dof                                number of dof to look for |
 | **node           address of pointer to hold controlled node (output) |
 | *cdof                      adress of INT to hold dof number (output) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(FIELD     *actfield,
                                INT        control_node_global,
                                INT        control_dof,
                                NODE     **node,
                                INT       *cdof) 
{
INT        i;
#ifdef DEBUG 
dstrc_enter("calstatserv_findcontroldof");
#endif
/*----------------------------------------------------------------------*/
*node = NULL;
for (i=0; i<actfield->dis[0].numnp; i++)
{
   if (actfield->dis[0].node[i].Id == control_node_global)
   {
      *node = &(actfield->dis[0].node[i]);
      *cdof = actfield->dis[0].node[i].dof[control_dof-1];
      break;
   }
}
if (!(*node)) dserror("Cannot find control node");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calstatserv_findcontroldof */
