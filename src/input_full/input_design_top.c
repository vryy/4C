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
 | create the connectivity between the design and the fields m.gee 4/01 |
 *----------------------------------------------------------------------*/
void inpdesign_topology()
{
int i,j,k,l;
DNODE *actdnode;
DLINE *actdline;
DSURF *actdsurf;
DVOL  *actdvol;
FIELD *actfield;
NODE  *actnode;
#ifdef DEBUG 
dstrc_enter("inpdesign_topology");
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
   actdline->node = (NODE**)calloc(actdline->mynode.fdim,sizeof(NODE*));
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
   actdsurf->node = (NODE**)calloc(actdsurf->mynode.fdim,sizeof(NODE*));
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
   actdvol->node = (NODE**)calloc(actdvol->mynode.fdim,sizeof(NODE*));
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
} /* end of inpdesign_topology */
