#include "../headers/standardtypes.h"

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
}
exit1:
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
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of inpdesign_topology */
