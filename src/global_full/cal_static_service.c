#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |  routine to find the node and dof to control in field  m.gee 11/01   |
 *----------------------------------------------------------------------*/
void calstatserv_findcontroldof(FIELD     *actfield,
                                int        control_node_global,
                                int        control_dof,
                                NODE     **node,
                                int       *cdof) 
{
int        i;
#ifdef DEBUG 
dstrc_enter("calstatserv_findcontroldof");
#endif
/*----------------------------------------------------------------------*/
*node = NULL;
for (i=0; i<actfield->numnp; i++)
{
   if (actfield->node[i].Id == control_node_global)
   {
      *node = &(actfield->node[i]);
      *cdof = actfield->node[i].dof[control_dof-1];
      break;
   }
}
if (!(*node)) dserror("Cannot find control node");
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of calstatserv_findcontroldof */
