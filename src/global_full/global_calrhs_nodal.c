#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  routine to assemble nodal neumann conditions         m.gee 10/01    |
 *----------------------------------------------------------------------*/
void assemble_nn(
                    PARTITION    *actpart,   /* my partition of the active field */
                    DIST_VECTOR  *rhs,       /* dist.vector */
                    double       *drhs       /* global redundant vector of size numeq_total */
                   )
{
int                   i,j;
int                   dof;
NODE                 *actnode;
COND_NODE            *actcond;

#ifdef DEBUG 
dstrc_enter("assemble_nn");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------- check type of field */
switch (actpart->fieldtyp)
{
/*----------------------------------------------------------------------*/
/*                      assembly of structural nodal neumann conditions */
/*----------------------------------------------------------------------*/
case structure:
   /***************************************************************************/
   /* Attention!!! this only works correct with partitions cut "Cut_Elements" */
   /*                                                             m.gee 01/02 */
   /***************************************************************************/
   /*--------------------------------------------------- start assembly */
   for (i=0; i<actpart->pdis[0].numnp; i++)
   {
      /*---------------------------------- check presence of conditions */
      if (!actpart->pdis[0].node[i]->c) continue;
      /*-------------------------- check presence of neumann conditions */
      if (!actpart->pdis[0].node[i]->c->isneum) continue;
      actnode = actpart->pdis[0].node[i];
      actcond = actnode->c;
      /*--------------------------------------- loop neumann conditions */
      for (j=0; j<actcond->neum_onoff.fdim; j++)
      {
         if (!(actcond->neum_onoff.a.iv[j])) continue;
         dof = actnode->dof[j];
         if (dof>=rhs->numeq_total) continue;
         /*----------------------------------------- assemble the value */
         drhs[dof] += actcond->neum_val.a.dv[j];
      }/* end loop j over dofs */
   }/* end loop i over nodes */
break;
/*----------------------------------------------------------------------*/
/*                           assembly of fluid nodal neumann conditions */
/*----------------------------------------------------------------------*/
case fluid:
   dserror("Fluidal nodal neumann conditions not yet implemented ");
break;
/*----------------------------------------------------------------------*/
/*                                                 unknown typ of field */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown typ of field");
/*----------------------------------------------------------------------*/
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of assemble_nn */
