#include "../headers/standardtypes.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
 | create the connectivity of the design                           1/02 |
 *----------------------------------------------------------------------*/
void inpdesign_topology_design()
{
int       i,j,k,l;
DNODE    *actdnode;
DLINE    *actdline;
DSURF    *actdsurf;
DVOL     *actdvol;
int       nodeid;
int       lineid;
int       surfid;

#ifdef DEBUG 
dstrc_enter("inpdesign_topology_design");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------- init the topology of the design to zero */
for (i=0; i<design->ndnode; i++) 
{
   design->dnode[i].ndline = 0;
   design->dnode[i].dline  = NULL;
}
for (i=0; i<design->ndline; i++) 
{
   design->dline[i].ndsurf = 0;
   design->dline[i].dsurf  = NULL;
   design->dline[i].dnode  = NULL;
}
for (i=0; i<design->ndsurf; i++) 
{
   design->dsurf[i].ndvol  = 0;
   design->dsurf[i].dvol   = NULL;
   design->dsurf[i].dline  = NULL;
}
for (i=0; i<design->ndvol; i++)
{
   design->dvol[i].dsurf   = NULL;
}
/*----------------------------- make connectivity from dlines to dnodes */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   actdline->dnode = (DNODE**)CALLOC(2,sizeof(DNODE*));
   if (!actdline->dnode) dserror("Allocation of memory failed");
   for (k=0; k<2; k++)
   {
      nodeid=actdline->my_dnodeId[k];
      for (j=0; j<design->ndnode; j++)
      {
         if (design->dnode[j].Id==nodeid) break;
      }
      actdline->dnode[k] = &(design->dnode[j]);
      actdline->dnode[k]->ndline++;
   }
}
/*----------------------------- make connectivity from dnodes to dlines */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   for (j=0; j<2; j++)
   {
      actdnode = actdline->dnode[j];
      if (actdnode->dline==NULL)
      {
         actdnode->dline = (DLINE**)CALLOC(actdnode->ndline,sizeof(DLINE*));
         if (!actdnode->dline) dserror("Allocation of memory failed");
         actdnode->dline[0] = actdline;
      }
      else
      {
         k=0;
         while (k<actdnode->ndline && actdnode->dline[k]!=NULL) k++;
         if (k==actdnode->ndline-1 && actdnode->dline[k]!=NULL)
         dserror("Cannot make dnode to dline topology");
         actdnode->dline[k] = actdline;
      }
   }
}
/*-------------------------- make connectivity from dsurfaces to dlines */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   actdsurf->dline = (DLINE**)CALLOC(actdsurf->my_dlineId.fdim,sizeof(DLINE*));
   if (!actdsurf->dline) dserror("Allocation of memory failed");
   for (k=0; k<actdsurf->my_dlineId.fdim; k++)
   {
      lineid = actdsurf->my_dlineId.a.ia[k][0];
      for (j=0; j<design->ndline; j++)
      {
         if (design->dline[j].Id==lineid) break;
      }
      actdsurf->dline[k] = &(design->dline[j]);
      actdsurf->dline[k]->ndsurf++;
   }
}
/*-------------------------- make connectivity from dlines to dsurfaces */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   for (j=0; j<actdsurf->my_dlineId.fdim; j++)
   {
      actdline = actdsurf->dline[j];
      if (actdline->dsurf==NULL)
      {
         actdline->dsurf = (DSURF**)CALLOC(actdline->ndsurf,sizeof(DSURF*));
         if (!actdline->dsurf) dserror("Allocation of memory failed");
         actdline->dsurf[0] = actdsurf;
      }
      else
      {
         k=0;
         while (k<actdline->ndsurf && actdline->dsurf[k]!=NULL) k++;
         if (k==actdline->ndsurf-1 && actdline->dsurf[k]!=NULL)
         dserror("Cannot make dline to dsurf topology");
         actdline->dsurf[k] = actdsurf;
      }
   }
}
/*----------------------  make connectivity from dvolumes to dsurfaces */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   actdvol->dsurf = (DSURF**)CALLOC(actdvol->my_dsurfId.fdim,sizeof(DSURF**));
   if (!actdvol->dsurf) dserror("Allocation of memory failed");
   for (k=0; k<actdvol->my_dsurfId.fdim; k++)
   {
      surfid = actdvol->my_dsurfId.a.ia[k][0];
      for (j=0; j<design->ndsurf; j++)
      {
         if (design->dsurf[j].Id==surfid) break;
      }
      actdvol->dsurf[k] = &(design->dsurf[j]);
      actdvol->dsurf[k]->ndvol++;
   }
}
/*----------------------- make connectivity from dsurfaces to dvolumes */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   for (j=0; j<actdvol->my_dsurfId.fdim; j++)
   {
      actdsurf = actdvol->dsurf[j];
      if (actdsurf->dvol==NULL)
      {
         actdsurf->dvol = (DVOL**)CALLOC(actdsurf->ndvol,sizeof(DVOL*));
         if (!actdsurf->dvol) dserror("Allocation of memory failed");
         actdsurf->dvol[0] = actdvol;
      }
      else
      {
         k=0;
         while (k<actdsurf->ndvol && actdsurf->dvol[k]!=NULL) k++;
         if (k==actdsurf->ndvol-1 && actdsurf->dvol[k]==NULL)
         dserror("Cannot make dsurf to dvol topology");
         actdsurf->dvol[k] = actdvol;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_topology_design */



/*----------------------------------------------------------------------*
 | create the connectivity between the design and the fields m.gee 4/01 |
 *----------------------------------------------------------------------*/
void inpdesign_topology_fe()
{
int i,j,k,l;
DNODE *actdnode;
DLINE *actdline;
DSURF *actdsurf;
DVOL  *actdvol;
FIELD *actfield;
NODE  *actnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_topology_fe");
#endif
/*----------------------------------------------- loop the design nodes */
for (i=0; i<design->ndnode; i++)
{
   actdnode = &(design->dnode[i]);
/*------------------------------- find the field this design is part off */   
   for (j=0; j<genprob.numfld; j++)
   {
      actfield = &(field[j]);
/*----------------------------------------- loop the nodes in this field */      
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdnode->mynode)
         {
            actdnode->node = &(actfield->node[k]);
            actdnode->field = actfield;
            goto exit1;
         }
      }
   }
exit1:;
}
/*----------------------------------------------- loop the design lines */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   actdline->node = (NODE**)CALLOC(actdline->mynode.fdim,sizeof(NODE*));
   if (actdline->node==NULL) dserror("Allocation of DLINE ptr to NODE failed");
/*------------------------------- find the field this design is part off */   
   for (j=0; j<genprob.numfld; j++)
   {
      actfield = &(field[j]);
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdline->mynode.a.iv[0])
         {
            actdline->field = actfield;
            goto exit2;
         }
      }
   }
   exit2:
/*--------------------- found the right field, so now find all the nodes */
   for (l=0; l<actdline->mynode.fdim; l++)
   {
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdline->mynode.a.iv[l])
         {
            actdline->node[l] = &(actfield->node[k]);
            break;
         }
      }
   }
}
/*--------------------------------------------- loop the design surfaces */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   actdsurf->node = (NODE**)CALLOC(actdsurf->mynode.fdim,sizeof(NODE*));
   if (actdsurf->node==NULL) dserror("Allocation of DSURF ptr to NODE failed");
/*------------------------------- find the field this design is part off */   
   for (j=0; j<genprob.numfld; j++)
   {
      actfield = &(field[j]);
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdsurf->mynode.a.iv[0])
         {
            actdsurf->field = actfield;
            goto exit3;
         }
      }
   }
   exit3:
/*--------------------- found the right field, so now find all the nodes */
   for (l=0; l<actdsurf->mynode.fdim; l++)
   {
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdsurf->mynode.a.iv[l])
         {
            actdsurf->node[l] = &(actfield->node[k]);
            break;
         }
      }
   }
}
/*--------------------------------------------- loop the design volumes */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   actdvol->node = (NODE**)CALLOC(actdvol->mynode.fdim,sizeof(NODE*));
   if (actdvol->node==NULL) dserror("Allocation of DVOL ptr to NODE failed");
/*------------------------------- find the field this design is part off */   
   for (j=0; j<genprob.numfld; j++)
   {
      actfield = &(field[j]);
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdvol->mynode.a.iv[0])
         {
            actdvol->field = actfield;
            goto exit4;
         }
      }
   }
   exit4:
/*--------------------- found the right field, so now find all the nodes */
   for (l=0; l<actdvol->mynode.fdim; l++)
   {
      for (k=0; k<actfield->numnp; k++)
      {
         if (actfield->node[k].Id == actdvol->mynode.a.iv[l])
         {
            actdvol->node[l] = &(actfield->node[k]);
            break;
         }
      }
   }
}
/*-- we made all pointers design -> FE-nodes, now make FE-nodes->design */
/*--------------------------------------- start with the design volumes */
for (i=0; i<design->ndvol; i++)
{
   actdvol = &(design->dvol[i]);
   for (j=0; j<actdvol->mynode.fdim; j++)
   {
      actnode            = actdvol->node[j];
      actnode->downertyp = dvol_owned;
      actnode->d.dvol   = actdvol;
   }
}
/*------------------------------------- make FE-nodes owned by surfaces */
for (i=0; i<design->ndsurf; i++)
{
   actdsurf = &(design->dsurf[i]);
   for (j=0; j<actdsurf->mynode.fdim; j++)
   {
      actnode            = actdsurf->node[j];
      actnode->downertyp = dsurf_owned;
      actnode->d.dsurf   = actdsurf;
   }
}
/*--------------------------------------- make FE-nodes owned by dlines */
for (i=0; i<design->ndline; i++)
{
   actdline = &(design->dline[i]);
   for (j=0; j<actdline->mynode.fdim; j++)
   {
      actnode            = actdline->node[j];
      actnode->downertyp = dline_owned;
      actnode->d.dline   = actdline;
   }
}
/*--------------------------------------- make FE-nodes owned by dnodes */
for (i=0; i<design->ndnode; i++)
{
   actdnode                  = &(design->dnode[i]);
   actdnode->node->downertyp = dnode_owned;
   actdnode->node->d.dnode   = actdnode;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_topology_fe */
