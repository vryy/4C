#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 |  routine to write all fields to gid - meshes          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh()
{
int           i,j,k,l,n;
FILE         *out = allfiles.gidmsh;
FIELD        *actfield;
ELEMENT      *actele;
NODE         *actnode;
INTRA        *actintra;
FIELDTYP     fieldtyp;

int           nummesh;
int          *discrets;
int           ndiscrets=10;
DIS_TYP       distyp;

char          sign='"';

#ifdef DEBUG 
dstrc_enter("out_gid_msh");
#endif
/*----------------------------------------------------------------------*/
discrets = (int*)malloc(ndiscrets*sizeof(int));
if (!discrets) dserror("Allocation of memory failed");
/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
/*-------------------------------------------------------- print a head */
fprintf(out,"#--------------------------------------------------------------------------------\n");
fprintf(out,"# P_CARAT postprocessing output to GID\n");
fprintf(out,"#--------------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
   actfield = &(field[i]);
   fieldtyp = actfield->fieldtyp;
/*------------ check how many types of discrets there are in this field */
   nummesh=0;
   for (j=0; j<ndiscrets; j++) discrets[j]=0;
   for (j=0; j<actfield->numele; j++)
   {
      actele = &(actfield->element[j]);
      switch(actele->distyp)
      {
      case quad4:discrets[0]=1;break;
      case quad8:discrets[1]=1;break;
      case quad9:discrets[2]=1;break; 
      case tri3 :discrets[3]=1;break;
      case tri6 :discrets[4]=1;break;
      case hex8 :discrets[5]=1;break;
      case hex20:discrets[6]=1;break;
      case hex27:discrets[7]=1;break;
      case tet4 :discrets[8]=1;break;
      case tet10:discrets[9]=1;break;
      default:
         dserror("Unknown type of element");
      break;
      }
   }
   for (j=0; j<ndiscrets; j++) 
   if (discrets[j]) nummesh++;
/*--------------------------------------- loop j as meshes in one field */
   for (j=0; j<nummesh; j++)
   {
/*------------------------------------------ print the different fields */
         fprintf(out,"#--------------------------------------------------------------------------------\n");
      switch(fieldtyp)
      {
      case structure:
         fprintf(out,"MESH %cstructure%d%c DIMENSION 3 ",sign,j,sign);
      break;
      case fluid:
         fprintf(out,"MESH %cfluid%d%c DIMENSION 3 ",sign,j,sign);
      break;
      case ale:
         fprintf(out,"MESH %cale%d%c DIMENSION 3 ",sign,j,sign);
      break;
      default:
         dserror("Unknown type of field");
      break;
      }
/*----------------------------------------------------- check for quad4 */ 
      if (discrets[0])
      {
         distyp=quad4;
         fprintf(out,"ELEMTYPE Quadrilateral NNODE 4\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[0]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for quad8 */ 
      if (discrets[1])
      {
         distyp=quad8;
         fprintf(out,"ELEMTYPE Quadrilateral NNODE 8\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[1]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for quad9 */ 
      if (discrets[2])
      {
         distyp=quad9;
         fprintf(out,"ELEMTYPE Quadrilateral NNODE 9\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[2]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for tri3  */ 
      if (discrets[3])
      {
         distyp=tri3;
         fprintf(out,"ELEMTYPE Triangle      NNODE 3\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[3]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for tri6  */ 
      if (discrets[4])
      {
         distyp=tri6;
         fprintf(out,"ELEMTYPE Triangle      NNODE 6\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         discrets[4]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for hex8  */ 
      if (discrets[5])
      {
         distyp=hex8;
         fprintf(out,"ELEMTYPE Hexahedra     NNODE 8\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[5]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for hex20 */ 
      if (discrets[6])
      {
         distyp=hex20;
         fprintf(out,"ELEMTYPE Hexahedra     NNODE 20\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[6]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for hex27 */ 
      if (discrets[7])
      {
         distyp=hex27;
         fprintf(out,"ELEMTYPE Hexahedra     NNODE 27\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[7]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for tet4  */ 
      if (discrets[8])
      {
         distyp=tet4;
         fprintf(out,"ELEMTYPE Tetrahedra    NNODE 4\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         discrets[8]=0;
         goto nextmesh;
      }
/*----------------------------------------------------- check for tet10 */ 
      if (discrets[9])
      {
         distyp=tet10;
         fprintf(out,"ELEMTYPE Tetrahedra    NNODE 10\n");
         fprintf(out,"# Nodal coordinates of all meshes, if this is first mesh\n");
         fprintf(out,"COORDINATES\n");
         /*---- check whether this is the very first mesh to be written */
         /*                          if so, write ALL nodal coordinates */
         if (i==0 && j==0) out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
         /*------------------- write element connectivity for this mesh */
         fprintf(out,"# Element connectivity of this mesh\n");
         fprintf(out,"ELEMENTS\n");
         out_gid_elements(out,actfield,distyp);
         fprintf(out,"END ELEMENTS\n");
         /*------------------ set flag, that this mesh has been printed */
         discrets[9]=0;
         goto nextmesh;
      }
/*----------------------------------------------------------------------*/
nextmesh:;
/*----------------------------------------------------------------------*/









/*----------------------------------------------------------------------*/
   } /* end of (j=0; j<nummesh; j++) */
}/* end of (i=0; i<genprob.numfld; i++)
/*----------------------------------------------------------------------*/
fflush(out);
free(discrets);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_msh */



/*----------------------------------------------------------------------*
 |  routine to write element connectivity of a certain   m.gee 12/01    |
 |  discretization type and a certain field to GID-file                 |
 *----------------------------------------------------------------------*/
int out_gid_elements(FILE *out, FIELD *actfield, DIS_TYP distyp)
{
int           i,j;
ELEMENT      *actele;
NODE         *actnode;
#ifdef DEBUG 
dstrc_enter("out_gid_elements");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actfield->numele; i++)
{
   if (actfield->element[i].distyp != distyp) continue;
   actele = &(actfield->element[i]);
   fprintf(out,"%6d ",actele->Id+1);
   for (j=0; j<actele->numnp; j++)
   fprintf(out,"%6d ",actele->node[j]->Id+1);
   fprintf(out,"\n");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_elements */



/*----------------------------------------------------------------------*
 |  routine to write all nodal coordinates of all fields m.gee 12/01    |
 *----------------------------------------------------------------------*/
int out_gid_allcoords(FILE *out)
{
int           i,j;
FIELD        *actfield;
NODE         *actnode;
#ifdef DEBUG 
dstrc_enter("out_gid_allcoords");
#endif
/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   for (j=0; j<actfield->numnp; j++)
   {
      actnode = &(actfield->node[j]);
      fprintf(out,"%6d %-18.5#f %-18.5#f %-18.5#f\n",
                                                     actnode->Id+1,
                                                     actnode->x[0],
                                                     actnode->x[1],
                                                     actnode->x[2]);
   }
}/* end of (i=0; i<genprob.numfld; i++)
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_allcoords */
