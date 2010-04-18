/*!----------------------------------------------------------------------
\file
\brief contains init phase of the shell contact routines

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef S8CONTACT
#include "../headers/standardtypes.h"
#include "s8contact.h"
#include "shell8.h"

/*!
\addtogroup CONTACTS8
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
struct _SHELLCONTACT shellcontact;
/*!---------------------------------------------------------------------
\brief initialization of shell8 contact

<pre>                                                        m.gee 2/03
</pre>
\param actfield   FIELD*       (i)   the discretization
\param actpart    PARTITION*   (i)   my partition of the discretization
\param actintra   INTRA*       (i)   the intra-communicator of this field
\return void

------------------------------------------------------------------------*/
void s8contact_init(FIELD *actfield, PARTITION* actpart, INTRA *actintra)
{
INT           i,j,k,l;
INT           myrank,nproc;
INT           numnp;
SHELLNODE    *actcnode;
NODE         *actnode;
ELEMENT      *actele;
DOUBLE        halfthick;
INT           foundit;
DOUBLE        diag;
DOUBLE        dx,dy,dz;
#ifdef DEBUG
dstrc_enter("s8contact_init");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*-------------------------------- allocate vector of all contact nodes */
shellcontact.numnp  = actfield->dis[0].numnp;
shellcontact.cnode  = (SHELLNODE*)CCACALLOC(shellcontact.numnp,sizeof(SHELLNODE));
/*--------------- loop nodes and init coordindates and ptrs to fe-nodes */
numnp  = shellcontact.numnp;
for (i=0; i<numnp; i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   actcnode = &(shellcontact.cnode[i]);
   /* set ptr to node */
   actcnode->node = actnode;
   /* find the director of this node */
   actele = actnode->element[0];
   for (j=0; j<actele->numnp; j++)
      if (actele->node[j] == actnode) break;
   if (j==actele->numnp) dserror("Cannot find director to node ");
   halfthick = (actele->e.s8->thick_node.a.dv[j])/2.0;
   actcnode->xr[0] = actnode->x[0];
   actcnode->xr[1] = actnode->x[1];
   actcnode->xr[2] = actnode->x[2];
   actcnode->xr[3] = actele->e.s8->a3ref.a.da[0][j]*halfthick;
   actcnode->xr[4] = actele->e.s8->a3ref.a.da[1][j]*halfthick;
   actcnode->xr[5] = actele->e.s8->a3ref.a.da[2][j]*halfthick;
}
/*------------------- loop nodes again and init the neighbourhood stack */
for (i=0; i<numnp; i++)
{
   actcnode = &(shellcontact.cnode[i]);
   if (actcnode->node->proc != myrank) continue;
   actcnode->nneigh = 0;
   for (j=0; j<actcnode->node->numele; j++)
   for (k=0; k<actcnode->node->element[j]->numnp; k++)
   {
      actnode = actcnode->node->element[j]->node[k];
      foundit = 0;
      for (l=0; l<actcnode->nneigh; l++)
      {
         if (actnode == actcnode->neighbours[l])
         {
            foundit = 1;
            break;
         }
      }
      if (!foundit)
      {
         actcnode->neighbours[actcnode->nneigh] = actnode;
         actcnode->nneigh++;
         if (actcnode->nneigh == 30) dserror("Neighbourhood stack full");
      }
   }
}
/*------------------------ calculate the maximum diagonal of an element */
shellcontact.maxdiag = 0.0;
for (i=0; i<actfield->dis[0].numele; i++)
{
   actele = &(actfield->dis[0].element[i]);
   dx     = actele->node[2]->x[0] - actele->node[0]->x[0];
   dy     = actele->node[2]->x[1] - actele->node[0]->x[1];
   dz     = actele->node[2]->x[2] - actele->node[0]->x[2];
   diag   = dx*dx+dy*dy+dz*dz;
   diag   = sqrt(diag);
   if (shellcontact.maxdiag < diag)
   shellcontact.maxdiag = diag;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s8contact_init */











/*! @} (documentation module close)*/

#endif
#endif
