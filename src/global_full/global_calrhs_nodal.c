#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  routine to assemble nodal neumann conditions         m.gee 10/01    |
 *----------------------------------------------------------------------*/
void assemble_nn(
                    PARTITION    *actpart,
                    DIST_VECTOR  *rhs,
                    double       *drhs
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
   /*--------------------------------------------------- start assembly */
   for (i=0; i<actpart->numnp; i++)
   {
      /*---------------------------------- check presence of conditions */
      if (!actpart->node[i]->c) continue;
      /*-------------------------- check presence of neumann conditions */
      if (!actpart->node[i]->c->isneum) continue;
      actnode = actpart->node[i];
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
