#include "../headers/standardtypes.h"
#include "gid.h"
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
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |  routine to write all fields to gid - meshes          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh()
{
int           i,j,k;
FILE         *out = allfiles.gidmsh;
FIELD        *actfield;
GIDSET       *actgid;
ELEMENT      *actele;
ELEMENT      *firstele;
NODE         *actnode;

int           is_firstmesh;
char          sign='"';

#ifdef DEBUG 
dstrc_enter("out_gid_msh");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------- print a head */
fprintf(out,"#--------------------------------------------------------------------------------\n");
fprintf(out,"# P_CARAT postprocessing output to GiD - MESH file\n");
fprintf(out,"#--------------------------------------------------------------------------------\n");
/*---------------------------------------------------- loop over fields */
is_firstmesh=1;
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   actgid   = &(gid[i]);
   /*----------------------------------- print the meshes of this field */

   if (actgid->is_shell8_22)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 2x2 GP\n",actgid->shell8_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 4\n",actgid->shell8_22_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_shell8_33)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 3x3 GP\n",actgid->shell8_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 9\n",actgid->shell8_33_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell8 || actele->numnp !=9) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_brick1_222)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BRICK1 2x2x2 GP\n",actgid->brick1_222_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->brick1_222_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_brick1 || actele->numnp !=8) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_brick1_333)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BRICK1 3x3x3 GP\n",actgid->brick1_333_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 27\n",actgid->brick1_333_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_brick1 || actele->numnp !=27) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_fluid3_222)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s FLUID3 2x2x2 GP\n",actgid->fluid3_222_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->fluid3_222_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_fluid3 || actele->numnp !=8) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_fluid3_333)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s FLUID3 3x3x3 GP\n",actgid->fluid3_333_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 27\n",actgid->fluid3_333_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_fluid3 || actele->numnp !=27) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_ale_22)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s ALE 2x2 GP\n",actgid->ale_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE 4\n",actgid->ale_22_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_ale2 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_ale_222)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s ALE 2x2x2 GP\n",actgid->ale_222_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_222_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_ale3 || actele->numnp !=8) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }

/*---------------------------------------------------------fh 06/02------*/     
   if (actgid->is_wall1_22)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALL1 2x2 GP\n",actgid->wall1_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 4\n",actgid->wall1_22_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_wall1 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_wall1_33)
   {  
      /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALL1 3x3 GP\n",actgid->wall1_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE  %d \n",actgid->wall1_33_name,firstele->numnp);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         out_gid_allcoords(out);
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_wall1 || actele->numnp < 8) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }
/*----------------------------------------------------------------------*/       
   
   
} /* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_msh */




/*----------------------------------------------------------------------*
 |  routine to write all nodal coordinates of all fields m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_allcoords(FILE *out)
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
   for (j=0; j<actfield->dis[0].numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);
      switch (genprob.ndim)
      {
	case 3:
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                                     actnode->Id+1,
                                                     actnode->x[0],
                                                     actnode->x[1],
                                                     actnode->x[2]);
        break;
	case 2:
        fprintf(out,"%6d %-18.5f %-18.5f \n",
                                                     actnode->Id+1,
                                                     actnode->x[0],
                                                     actnode->x[1]);
        break;
	default:
        dserror("Unknown number of dimensions");
	break;
      }
   }
}/* end of (i=0; i<genprob.numfld; i++)
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of out_gid_allcoords */

