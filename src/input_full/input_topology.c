#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | create the node-element topology for this field        m.gee 4/01    |
 *----------------------------------------------------------------------*/
void inp_topology(FIELD *field)
{
int  i,j,k;
int  node_id;
ELEMENT *actele;
NODE    *actnode;
#ifdef DEBUG 
dstrc_enter("inp_topology");
#endif
/*------------------------------- create pointer from elements to nodes */
for (i=0; i<field->dis[0].numele; i++)
{
   actele = &(field->dis[0].element[i]);
/*--------------------------------- allocate the ELEMENTs node-pointers */
   actele->node = (NODE**)CALLOC(actele->numnp,sizeof(NODE*));
   if (actele->node==NULL) dserror("Allocation of node pointers failed");
   
   for (j=0; j<actele->numnp; j++)
   {
      node_id = actele->lm[j];
      
      for (k=0; k<field->dis[0].numnp; k++)
      {
         if (field->dis[0].node[k].Id == node_id)
         {
            actele->node[j] = &(field->dis[0].node[k]);
            break;
         }
      } /* end of loop over all nodes */
   } /* end of loop over elements nodes */
}/* end of loop over elements */

/*------------------------------ create pointers from nodes to elements */
for (i=0; i<field->dis[0].numnp; i++) field->dis[0].node[i].numele=0;
/*---------------------------- count the number of elements to one node */
for (i=0; i<field->dis[0].numele; i++)
{
   actele = &(field->dis[0].element[i]);
   for (j=0; j<actele->numnp; j++)
   {
      (actele->node[j]->numele)++;
   }
}
/*------------------------- allocate space for element pointers in NODE */
for (i=0; i<field->dis[0].numnp; i++)
{
   actnode = &(field->dis[0].node[i]);
   actnode->element = (ELEMENT**)CALLOC(actnode->numele,sizeof(ELEMENT*));
   if (actnode->element==NULL) dserror("Allocation of element pointers failed");
   for (j=0; j<actnode->numele; j++) actnode->element[j]=NULL;
}
/*---------------- loop elements and point from their nodes to themself */
for (i=0; i<field->dis[0].numele; i++)
{
   actele = &(field->dis[0].element[i]);
   for (j=0; j<actele->numnp; j++)
   {
      actnode = actele->node[j];
      for (k=0; k<actnode->numele; k++)
      {
         if (actnode->element[k]==NULL) break;
      }
      actnode->element[k]=actele;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inp_topology */
