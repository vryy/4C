/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

------------------------------------------------------------------------*/

#include "../headers/standardtypes.h"
#include "gid.h"
#ifdef D_SHELL8
  #include "../shell8/shell8.h"
#endif /*D_SHELL8*/
#ifdef D_SHELL9
  #include "../shell9/shell9.h"
#endif /*D_SHELL9*/
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
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
INT           i,j,k;
FILE         *out = allfiles.gidmsh;
FIELD        *actfield;
GIDSET       *actgid;
ELEMENT      *actele;
ELEMENT      *firstele;
NODE         *actnode;

DOUBLE        a1,a2,a3,thick,scal;
INT           tot_numnp;

INT           is_firstmesh;
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
#ifdef D_SHELL8
#if 0
   /* this is the quadrilateral version */
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
#endif
   /* this is the hexahedra version */
   if (actgid->is_shell8_22)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 2x2x2 GP\n",actgid->shell8_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->shell8_22_name);
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         tot_numnp = genprob.nnode;
         for (i=0; i<actfield->dis[0].numnp; i++)
         {
            actnode = &(actfield->dis[0].node[i]);
            /* find the director */
            actele  = actnode->element[0];
            for (j=0; j<actele->numnp; j++)
               if (actele->node[j] == actnode) break;
            if (j==4) dserror("Cannot find matching node in element");
            thick = actele->e.s8->thick;
            scal  = 1.0;
            a1 = actele->e.s8->a3ref.a.da[0][j]*thick*scal/2.0;
            a2 = actele->e.s8->a3ref.a.da[1][j]*thick*scal/2.0;
            a3 = actele->e.s8->a3ref.a.da[2][j]*thick*scal/2.0;
            /* the lower surface coordinate*/
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",actnode->Id+1,
                                                        actnode->x[0]-a1,
                                                        actnode->x[1]-a2,
                                                        actnode->x[2]-a3);
            /* the upper surface coordinate */
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",actnode->Id+1+tot_numnp,
                                                        actnode->x[0]+a1,
                                                        actnode->x[1]+a2,
                                                        actnode->x[2]+a3);
         }
         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell8 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         /* lower surface nodes */
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         /* upper surface nodes */
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1+tot_numnp);
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
#endif /*D_SHELL8*/

#ifdef D_SHELL9
   /* 4-noded shell9 element -> Hex8 */
   if (actgid->is_shell9_4_22 || actgid->is_shell9_4_33)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      if (actgid->is_shell9_4_22) /*2x2 gp*/
      {   
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->shell9_4_22_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->shell9_4_22_name);
      }
      else if (actgid->is_shell9_4_33) /*3x3 gp*/
      {   
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->shell9_4_33_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->shell9_4_33_name);
      }
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");
         
         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_allcoords_unsmo(out,4);    /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_allcoords(out,4);          /*  smoothed stresses to gid*/

         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell9 || actele->numnp !=4) continue;

         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_eletop_unsmo(out,4,j);    /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_eletop(out,4,actele);     /*  smoothed stresses to gid*/

      }
      fprintf(out,"END ELEMENTS\n");
   }

   /* 8-noded shell9 element -> Hex20 */
   if (actgid->is_shell9_8_22 || actgid->is_shell9_8_33)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      if (actgid->is_shell9_8_22) /*2x2 gp*/
      {
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->shell9_8_22_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->shell9_8_22_name);
      }
      else if (actgid->is_shell9_8_33) /*3x3 gp*/
      {
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->shell9_8_33_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->shell9_8_33_name);
      }
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");

         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_allcoords_unsmo(out,8);   /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_allcoords(out,8);         /*  smoothed stresses to gid*/

         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell9 || actele->numnp !=8) continue;

         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_eletop_unsmo(out,8,j);    /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_eletop(out,8,actele);     /*  smoothed stresses to gid*/

      }
      fprintf(out,"END ELEMENTS\n");
   }

   /* 9-noded shell9 element -> Hex27 */
   if (actgid->is_shell9_9_22 || actgid->is_shell9_9_33)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      if (actgid->is_shell9_9_22) /*2x2 gp*/
      {
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 2x2x2 GP\n",actgid->shell9_9_22_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 27\n",actgid->shell9_9_22_name);
      }
      if (actgid->is_shell9_9_33) /*3x3 gp*/
      {
       fprintf(out,"# MESH %s FOR FIELD %s SHELL9 3x3x3 GP\n",actgid->shell9_9_33_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 27\n",actgid->shell9_9_33_name);
      }
      /*-------------- if this is first mesh, print coodinates of nodes */
      if (is_firstmesh)
      {
         is_firstmesh=0;
         fprintf(out,"# printout ALL nodal coordinates of ALL fields in first mesh only\n");
         fprintf(out,"COORDINATES\n");

         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_allcoords_unsmo(out,9); /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_allcoords(out,9);       /*  smoothed stresses to gid*/

         fprintf(out,"END COORDINATES\n");
      }
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_shell9 || actele->numnp !=9) continue;

         if (ioflags.struct_stress_gid_smo == 0 ) s9_out_gid_eletop_unsmo(out,9,j);    /*unsmoothed stresses to gid*/
         if (ioflags.struct_stress_gid_smo == 1 ) s9_out_gid_eletop(out,9,actele);     /*  smoothed stresses to gid*/

      }
      fprintf(out,"END ELEMENTS\n");
   }
#endif /*D_SHELL9*/


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
      fprintf(out,"# MESH %s FOR FIELD %s BRICK1 3x3x3 GP as HEX8!\n",actgid->brick1_333_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      /*fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->brick1_333_name);*/
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->brick1_333_name);
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
      /* gid shows Hex20 as Hex with eight nodes only, 
                                               so the output is reduced */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
         actele = &(actfield->dis[0].element[j]);
         if (actele->eltyp != el_brick1 || actele->numnp !=20) continue;
         fprintf(out," %6d ",actele->Id+1);
         /*for (k=0; k<actele->numnp; k++)*/
         for (k=0; k<8; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_fluid2_22)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s FLUID2 2x2 GP\n",actgid->fluid2_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE 4\n",actgid->fluid2_22_name);
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
         if (actele->eltyp != el_fluid2 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_fluid2_33)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s FLUID2 3x3 GP\n",actgid->fluid2_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE 9\n",actgid->fluid2_33_name);
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
         if (actele->eltyp != el_fluid2 || actele->numnp !=9) continue;
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

   if (genprob.probtyp != prb_fsi)
   {
     if (actgid->is_ale_11)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 1x1 GP\n",actgid->ale_11_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE 4\n",actgid->ale_11_name);
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


     if (actgid->is_ale_tri_1)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 1 GP\n",actgid->ale_tri_1_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Triangle NNODE 3\n",actgid->ale_tri_1_name);
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
         if (actele->eltyp != el_ale2 || actele->numnp !=3) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
           fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
       }
       fprintf(out,"END ELEMENTS\n");
     }


     if (actgid->is_ale_tri_3)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 3 GP\n",actgid->ale_tri_3_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Triangle NNODE 3\n",actgid->ale_tri_3_name);
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
         if (actele->eltyp != el_ale2 || actele->numnp !=3) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
           fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
       }
       fprintf(out,"END ELEMENTS\n");
     }


     if (actgid->is_ale_111)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 1x1x1 GP\n",actgid->ale_111_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_111_name);
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


     if (actgid->is_ale_tet_1)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 1 GP\n",actgid->ale_tet_1_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Tetrahedra NNODE 4\n",actgid->ale_tet_1_name);
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
         if (actele->eltyp != el_ale3 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
           fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
       }
       fprintf(out,"END ELEMENTS\n");
     }


     if (actgid->is_ale_tet_4)
     {
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"# MESH %s FOR FIELD %s ALE 4 GP\n",actgid->ale_tet_4_name,actgid->fieldname);
       fprintf(out,"#-------------------------------------------------------------------------------\n");
       fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Tetrahedra NNODE 4\n",actgid->ale_tet_4_name);
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
         if (actele->eltyp != el_ale3 || actele->numnp !=4) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
           fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
       }
       fprintf(out,"END ELEMENTS\n");
     }
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

/*---------------------------------------------------------fh 05/03------*/     
   if (actgid->is_beam3_21)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BEAM3 1 GP\n",actgid->beam3_21_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Linear NNODE 2\n",actgid->beam3_21_name);
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
         if (actele->eltyp != el_beam3 || actele->numnp !=2) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_beam3_22)
   {  
      /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BEAM3 2 GP\n",actgid->beam3_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Linear NNODE  %d \n",actgid->beam3_22_name,firstele->numnp);
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
         if (actele->eltyp != el_beam3 || actele->numnp !=2) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_beam3_32)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BEAM3 2 GP\n",actgid->beam3_32_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Linear NNODE 3\n",actgid->beam3_32_name);
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
         if (actele->eltyp != el_beam3 || actele->numnp !=3) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }


   if (actgid->is_beam3_33)
   {  
      /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s BEAM3 3 GP\n",actgid->beam3_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Linear NNODE  %d \n",actgid->beam3_33_name,firstele->numnp);
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
         if (actele->eltyp != el_beam3 || actele->numnp !=3) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }

/*----------------------------------------------------------------------*/       
   if (actgid->is_axishell)
   {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s AXISHELL 5 GP\n",actgid->axishell_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Linear NNODE 2\n",actgid->ale_tet_4_name);
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
         if (actele->eltyp != el_axishell || actele->numnp !=2) continue;
         fprintf(out," %6d ",actele->Id+1);
         for (k=0; k<actele->numnp; k++)
         fprintf(out,"%6d ",actele->node[k]->Id+1);
         fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
   }

   
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
INT           i,j;
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

