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
 |                                                       ah 03/04       |
 | vector of numfld submesh-FIELDs, defined in global_control.c         |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *sm_field;
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

/*------------------------ global variable needed by gid postprocessing */
GIDSET *gid;

#ifdef D_MLSTRUCT
/*----------------------------- the same for submesh gid postprocessing */
GIDSET *sm_gid;
#endif /* D_MLSTRUCT */

/*----------------------------------------------------------------------*
 |  routine to write all fields to gid - meshes          m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh()
{

#ifndef NO_TEXT_OUTPUT


  INT           i,j,k,kk;
  FILE         *out = allfiles.gidmsh;
  FIELD        *actfield;
  GIDSET       *actgid;
  ELEMENT      *actele;
  ELEMENT      *firstele;

#ifdef D_SHELL8
  NODE         *actnode = NULL;
  DOUBLE        a1,a2,a3,thick,scal;
  INT           tot_numnp;
#endif

  INT           is_firstmesh;

#ifdef DEBUG
  dstrc_enter("out_gid_msh");
#endif


  /* print a head */
  fprintf(out,"#--------------------------------------------------------------------------------\n");
  fprintf(out,"# CCARAT postprocessing output to GiD - MESH file\n");
  fprintf(out,"#--------------------------------------------------------------------------------\n");



  /* write out al nodes first */
  /*==========================*/

  /* loop over fields */
  is_firstmesh=1;
  for (kk=0; kk<genprob.numfld; kk++)
  {
    actfield = &(field[kk]);
    actgid   = &(gid[kk]);


    /* write the openeing for the coordinates */
    if (is_firstmesh)
    {
      fprintf(out,"MESH node_mesh DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n");
      fprintf(out,"# printout ALL coordinates of ALL fields in first mesh only\n");
      fprintf(out,"COORDINATES\n");
    }


#ifdef D_SHELL9
    /* deal with shell9 elements */
    if (is_firstmesh)
    {
      if (actgid->is_shell9_4_22 || actgid->is_shell9_4_33)
      {
        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_allcoords_unsmo(out,4);    /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_allcoords(out,4);          /*  smoothed stresses to gid*/
        goto elements;
      }
      else if (actgid->is_shell9_8_22 || actgid->is_shell9_8_33)
      {
        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_allcoords_unsmo(out,8);   /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_allcoords(out,8);         /*  smoothed stresses to gid*/
        goto elements;
      }
      else if (actgid->is_shell9_9_22 || actgid->is_shell9_9_33)
      {
        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_allcoords_unsmo(out,9); /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_allcoords(out,9);       /*  smoothed stresses to gid*/
        goto elements;
      }
    }
    else
    {
      if (actgid->is_shell9_4_22 || actgid->is_shell9_4_33 ||
          actgid->is_shell9_8_22 || actgid->is_shell9_8_33 ||
          actgid->is_shell9_9_22 || actgid->is_shell9_9_33)
        dserror("shell9 elements must be on the first field for output to gid!!!");
    }
#endif


    /* write all "normal" coordinates */
    if (is_firstmesh)
    {
      is_firstmesh=0;
      out_gid_allcoords(out);
    }



#ifdef SPLIT_HEX20
    /* if the mesh has hex20 elements create additional nodes */
    if (actgid->is_ale_20_111 ||
        actgid->is_ale_20_222 ||
        actgid->is_ale_20_333 ||
        actgid->is_f3f_20_222 ||
        actgid->is_f3f_20_333 )
    {

      if ( genprob.nele > ELESHIFT_HEX20 )
        dserror("The value of ELESHIFT_HEX20 is too small. It must be larger than %d!!\n",
            genprob.nele);

      if ( genprob.nnode > NODESHIFT_HEX20 )
        dserror("The value of NODESHIFT_HEX20 is too small. It must be larger than %d!!\n",
            genprob.nnode);



      for (j=0; j<actfield->dis[0].numele; j++)
      {

        INT l;
        INT eleid,nodebase;
        DOUBLE x[3];
        actele = &(actfield->dis[0].element[j]);
        if ( !(actele->eltyp == el_ale3 || actele->eltyp == el_fluid3_fast)
            || actele->numnp !=20) continue;

        eleid = actele->Id+1;
        nodebase = eleid * NODESHIFT_HEX20;

        /* node 01 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[ 8]->x[l] + actele->node[ 9]->x[l]
              + actele->node[10]->x[l] + actele->node[11]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 1, x[0], x[1], x[2]);

        /* node 02 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[16]->x[l] + actele->node[17]->x[l]
              + actele->node[18]->x[l] + actele->node[19]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 2, x[0], x[1], x[2]);

        /* node 03 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[ 8]->x[l] + actele->node[13]->x[l]
              + actele->node[16]->x[l] + actele->node[12]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 3, x[0], x[1], x[2]);

        /* node 04 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[ 9]->x[l] + actele->node[14]->x[l]
              + actele->node[17]->x[l] + actele->node[13]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 4, x[0], x[1], x[2]);

        /* node 05 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[10]->x[l] + actele->node[14]->x[l]
              + actele->node[18]->x[l] + actele->node[15]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 5, x[0], x[1], x[2]);

        /* node 06 */
        for (l=0;l<3;l++)
          x[l] = 0.25*(actele->node[11]->x[l] + actele->node[15]->x[l]
              + actele->node[19]->x[l] + actele->node[12]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 6, x[0], x[1], x[2]);

        /* node 07 */
        for (l=0;l<3;l++)
          x[l] = 0.0833333333
            *(actele->node[ 8]->x[l] + actele->node[ 9]->x[l]
                + actele->node[10]->x[l] + actele->node[11]->x[l]
                + actele->node[12]->x[l] + actele->node[13]->x[l]
                + actele->node[14]->x[l] + actele->node[15]->x[l]
                + actele->node[16]->x[l] + actele->node[17]->x[l]
                + actele->node[18]->x[l] + actele->node[19]->x[l] );
        fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
            nodebase + 7, x[0], x[1], x[2]);

      }  /* for (j=0; j<actfield->dis[0].numele; j++) */

    }  /*     actgid->is_f3f_20_333 ) */
#endif


  }  /* for (kk=0; kk<genprob.numfld; kk++) */

  fprintf(out,"END COORDINATES\n");



elements:
  /* now print the elements of all fields */



  /* loop over fields */
  for (kk=0; kk<genprob.numfld; kk++)
  {
    actfield = &(field[kk]);
    actgid   = &(gid[kk]);


#ifdef D_SHELL8
    if (actgid->is_shell8_4_22)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 2x2 GP\n",actgid->shell8_4_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 4\n",actgid->shell8_4_22_name);

      /* print elements */
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

    if (actgid->is_shell8_4_33)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 3x3 GP\n",actgid->shell8_4_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 4\n",actgid->shell8_4_33_name);

      /* print elements */
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

    if (actgid->is_shell8_8_22)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 2x2 GP\n",actgid->shell8_8_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 8\n",actgid->shell8_8_22_name);

      /* print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_shell8 || actele->numnp !=8) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }

    if (actgid->is_shell8_8_33)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 3x3 GP\n",actgid->shell8_8_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 8\n",actgid->shell8_8_33_name);

      /* print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_shell8 || actele->numnp !=8) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }

    if (actgid->is_shell8_9_22)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 2x2 GP\n",actgid->shell8_9_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 8\n",actgid->shell8_9_22_name);

      /* print elements */
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

    if (actgid->is_shell8_9_33)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s SHELL8 3x3 GP\n",actgid->shell8_9_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE 8\n",actgid->shell8_9_33_name);

      /* print elements */
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
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_shell9 || actele->numnp !=4) continue;

        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_eletop_unsmo(out,4,j);    /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_eletop(out,4,actele);     /*  smoothed stresses to gid*/

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
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_shell9 || actele->numnp !=8) continue;

        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_eletop_unsmo(out,8,j);    /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_eletop(out,8,actele);     /*  smoothed stresses to gid*/

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
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_shell9 || actele->numnp !=9) continue;

        if (ioflags.struct_stress_smo == 0 ) s9_out_gid_eletop_unsmo(out,9,j);    /*unsmoothed stresses to gid*/
        if (ioflags.struct_stress_smo == 1 ) s9_out_gid_eletop(out,9,actele);     /*  smoothed stresses to gid*/

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



#ifdef D_FLUID3_F
    if (actgid->is_f3f_8_222 ||
        actgid->is_f3f_8_333)
    {

      fprintf(out,"#-------------------------------------------------------------------------------\n");

      if (actgid->is_f3f_8_222)
      {
        fprintf(out,"# MESH %s FOR FIELD %s FLUID3_FAST 2x2x2 GP\n",actgid->f3f_8_222_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->f3f_8_222_name);
      }
      if (actgid->is_f3f_8_333)
      {
        fprintf(out,"# MESH %s FOR FIELD %s FLUID3_FAST 3x3x3 GP\n",actgid->f3f_8_333_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->f3f_8_333_name);
      }

      /* print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_fluid3_fast || actele->numnp !=8) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }




    if (actgid->is_f3f_20_222 ||
        actgid->is_f3f_20_333)
    {

      fprintf(out,"#-------------------------------------------------------------------------------\n");

      if (actgid->is_f3f_20_222)
      {
        fprintf(out,"# MESH %s FOR FIELD %s FLUID3_FAST 2x2x2 GP\n",actgid->f3f_20_222_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
#ifndef SPLIT_HEX20
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->f3f_20_222_name);
#else
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->f3f_20_222_name);
#endif
      }
      if (actgid->is_f3f_20_333)
      {
        fprintf(out,"# MESH %s FOR FIELD %s FLUID3_FAST 3x3x3 GP\n",actgid->f3f_20_333_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
#ifndef SPLIT_HEX20
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->f3f_20_333_name);
#else
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->f3f_20_333_name);
#endif
      }

      /* print elements */
      fprintf(out,"ELEMENTS\n");

#ifndef SPLIT_HEX20
      /* write hex20 als hex20 elements */
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_fluid3_fast || actele->numnp !=20) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }

#else

      /* write hex20 elements as 8 hex8 elements */
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        INT l;
        INT eleid,nodebase,elebase;
        DOUBLE x[3];
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_fluid3_fast || actele->numnp !=20) continue;
        eleid = actele->Id+1;
        elebase  = eleid * ELESHIFT_HEX20;
        nodebase = eleid * NODESHIFT_HEX20;

        /* ele 01 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 1,
            actele->node[0]->Id+1,
            actele->node[8]->Id+1,
            nodebase + 1,
            actele->node[11]->Id+1,
            actele->node[12]->Id+1,
            nodebase + 3,
            nodebase + 7,
            nodebase + 6);

        /* ele 02 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 2,
            actele->node[8]->Id+1,
            actele->node[1]->Id+1,
            actele->node[9]->Id+1,
            nodebase + 1,
            nodebase + 3,
            actele->node[13]->Id+1,
            nodebase + 4,
            nodebase + 7);

        /* ele 03 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 3,
            nodebase + 1,
            actele->node[9]->Id+1,
            actele->node[2]->Id+1,
            actele->node[10]->Id+1,
            nodebase + 7,
            nodebase + 4,
            actele->node[14]->Id+1,
            nodebase + 5);

        /* ele 04 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 4,
            actele->node[11]->Id+1,
            nodebase + 1,
            actele->node[10]->Id+1,
            actele->node[3]->Id+1,
            nodebase + 6,
            nodebase + 7,
            nodebase + 5,
            actele->node[15]->Id+1);

        /* ele 05 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 5,
            actele->node[12]->Id+1,
            nodebase + 3,
            nodebase + 7,
            nodebase + 6,
            actele->node[4]->Id+1,
            actele->node[16]->Id+1,
            nodebase + 2,
            actele->node[19]->Id+1);

        /* ele 06 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 6,
            nodebase + 3,
            actele->node[13]->Id+1,
            nodebase + 4,
            nodebase + 7,
            actele->node[16]->Id+1,
            actele->node[5]->Id+1,
            actele->node[17]->Id+1,
            nodebase + 2);

        /* ele 07 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 7,
            nodebase + 7,
            nodebase + 4,
            actele->node[14]->Id+1,
            nodebase + 5,
            nodebase + 2,
            actele->node[17]->Id+1,
            actele->node[6]->Id+1,
            actele->node[18]->Id+1);

        /* ele 08 */
        fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
            elebase + 8,
            nodebase + 6,
            nodebase + 7,
            nodebase + 5,
            actele->node[15]->Id+1,
            actele->node[19]->Id+1,
            nodebase + 2,
            actele->node[18]->Id+1,
            actele->node[7]->Id+1);

      }  /* for (j=0; j<actfield->dis[0].numele; j++) */

#endif

      fprintf(out,"END ELEMENTS\n");

    }  /*     actgid->is_f3f_20_333) */


#endif


    if (genprob.probtyp != prb_fsi)
    {
      if (actgid->is_ale_11)
      {
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# MESH %s FOR FIELD %s ALE 1x1 GP\n",actgid->ale_11_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE 4\n",actgid->ale_11_name);
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


      if (actgid->is_ale_8_111 ||
          actgid->is_ale_8_222 ||
          actgid->is_ale_8_333)
      {
        fprintf(out,"#-------------------------------------------------------------------------------\n");

        if (actgid->is_ale_8_111)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 1x1x1 GP\n",actgid->ale_8_111_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_8_111_name);
        }
        if (actgid->is_ale_8_222)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 2x2x2 GP\n",actgid->ale_8_222_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_8_222_name);
        }
        if (actgid->is_ale_8_333)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 3x3x3 GP\n",actgid->ale_8_333_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_8_333_name);
        }

        /* print elements */
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




      if (actgid->is_ale_20_111 ||
          actgid->is_ale_20_222 ||
          actgid->is_ale_20_333)
      {
        fprintf(out,"#-------------------------------------------------------------------------------\n");

        if (actgid->is_ale_20_111)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 1x1x1 GP\n",actgid->ale_20_111_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
#ifndef SPLIT_HEX20
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->ale_20_111_name);
#else
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_20_111_name);
#endif
        }
        if (actgid->is_ale_20_222)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 2x2x2 GP\n",actgid->ale_20_222_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
#ifndef SPLIT_HEX20
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->ale_20_222_name);
#else
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_20_222_name);
#endif
        }
        if (actgid->is_ale_20_333)
        {
          fprintf(out,"# MESH %s FOR FIELD %s ALE 3x3x3 GP\n",actgid->ale_20_333_name,actgid->fieldname);
          fprintf(out,"#-------------------------------------------------------------------------------\n");
#ifndef SPLIT_HEX20
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 20\n",actgid->ale_20_333_name);
#else
          fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Hexahedra NNODE 8\n",actgid->ale_20_333_name);
#endif
        }



        /* print elements */
        fprintf(out,"ELEMENTS\n");

#ifndef SPLIT_HEX20
        /* write hex20 als hex20 elements */
        for (j=0; j<actfield->dis[0].numele; j++)
        {
          actele = &(actfield->dis[0].element[j]);
          if (actele->eltyp != el_ale3 || actele->numnp !=20) continue;
          fprintf(out," %6d ",actele->Id+1);
          for (k=0; k<actele->numnp; k++)
            fprintf(out,"%6d ",actele->node[k]->Id+1);
          fprintf(out,"\n");
        }

#else

        /* write hex20 elements as 8 hex8 elements */
        for (j=0; j<actfield->dis[0].numele; j++)
        {
          INT l;
          INT eleid,nodebase,elebase;
          DOUBLE x[3];
          actele = &(actfield->dis[0].element[j]);
          if (actele->eltyp != el_ale3 || actele->numnp !=20) continue;
          eleid = actele->Id+1;
          elebase  = eleid * ELESHIFT_HEX20;
          nodebase = eleid * NODESHIFT_HEX20;

          /* ele 01 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 1,
              actele->node[0]->Id+1,
              actele->node[8]->Id+1,
              nodebase + 1,
              actele->node[11]->Id+1,
              actele->node[12]->Id+1,
              nodebase + 3,
              nodebase + 7,
              nodebase + 6);

          /* ele 02 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 2,
              actele->node[8]->Id+1,
              actele->node[1]->Id+1,
              actele->node[9]->Id+1,
              nodebase + 1,
              nodebase + 3,
              actele->node[13]->Id+1,
              nodebase + 4,
              nodebase + 7);

          /* ele 03 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 3,
              nodebase + 1,
              actele->node[9]->Id+1,
              actele->node[2]->Id+1,
              actele->node[10]->Id+1,
              nodebase + 7,
              nodebase + 4,
              actele->node[14]->Id+1,
              nodebase + 5);

          /* ele 04 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 4,
              actele->node[11]->Id+1,
              nodebase + 1,
              actele->node[10]->Id+1,
              actele->node[3]->Id+1,
              nodebase + 6,
              nodebase + 7,
              nodebase + 5,
              actele->node[15]->Id+1);

          /* ele 05 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 5,
              actele->node[12]->Id+1,
              nodebase + 3,
              nodebase + 7,
              nodebase + 6,
              actele->node[4]->Id+1,
              actele->node[16]->Id+1,
              nodebase + 2,
              actele->node[19]->Id+1);

          /* ele 06 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 6,
              nodebase + 3,
              actele->node[13]->Id+1,
              nodebase + 4,
              nodebase + 7,
              actele->node[16]->Id+1,
              actele->node[5]->Id+1,
              actele->node[17]->Id+1,
              nodebase + 2);

          /* ele 07 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 7,
              nodebase + 7,
              nodebase + 4,
              actele->node[14]->Id+1,
              nodebase + 5,
              nodebase + 2,
              actele->node[17]->Id+1,
              actele->node[6]->Id+1,
              actele->node[18]->Id+1);

          /* ele 08 */
          fprintf(out," %6d  %6d %6d %6d %6d %6d %6d %6d %6d\n",
              elebase + 8,
              nodebase + 6,
              nodebase + 7,
              nodebase + 5,
              actele->node[15]->Id+1,
              actele->node[19]->Id+1,
              nodebase + 2,
              actele->node[18]->Id+1,
              actele->node[7]->Id+1);

        }
#endif

        fprintf(out,"END ELEMENTS\n");

      }


      if (actgid->is_ale_tet_1)
      {
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"# MESH %s FOR FIELD %s ALE 1 GP\n",actgid->ale_tet_1_name,actgid->fieldname);
        fprintf(out,"#-------------------------------------------------------------------------------\n");
        fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Tetrahedra NNODE 4\n",actgid->ale_tet_1_name);
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
    if (actgid->is_wall1_11)
    {
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALL1 1x1 GP\n",actgid->wall1_11_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Triangle NNODE 3\n",actgid->wall1_11_name);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_wall1 || actele->numnp !=3) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }


    if (actgid->is_wall1_22)
    {
      /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALL1 2x2 GP\n",actgid->wall1_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE %d \n",actgid->wall1_22_name,firstele->numnp);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_wall1 ) continue;
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
    /*----------------------------------------------------------------------*/
    if (actgid->is_interf_22)
    {
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s INERFACE 2(+2) GP\n",actgid->interf_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actgid->interf_22_name,firstele->numnp);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_interf) continue;
        /*----------------------------------------- (|| actele->numnp !=4) */
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }
    if (actgid->is_interf_33)
    {
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s INERFACE 3(+2x3) GP\n",actgid->interf_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actgid->interf_33_name,firstele->numnp);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_interf) continue;
        /*----------------------------------------- (|| actele->numnp !=8) */
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }
    /*----------------------------------------------------------------------*/
    if (actgid->is_wallge_22)
    {
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALLGE 2x2 GP\n",actgid->wallge_22_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actgid->wallge_22_name,firstele->numnp);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_wallge) continue;
        fprintf(out," %6d ",actele->Id+1);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id+1);
        fprintf(out,"\n");
      }
      fprintf(out,"END ELEMENTS\n");
    }
    /*----------------------------------------------------------------------*/
    if (actgid->is_wallge_33)
    {
      firstele = &(actfield->dis[0].element[0]);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"# MESH %s FOR FIELD %s WALLGE 3x3 GP\n",actgid->wallge_33_name,actgid->fieldname);
      fprintf(out,"#-------------------------------------------------------------------------------\n");
      fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actgid->wallge_33_name,firstele->numnp);
      /*------------------------------------------------ print elements */
      fprintf(out,"ELEMENTS\n");
      for (j=0; j<actfield->dis[0].numele; j++)
      {
        actele = &(actfield->dis[0].element[j]);
        if (actele->eltyp != el_wallge) continue;
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
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_msh */









#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  routine to write the submesh to gid                      ah 8/04    |
 *----------------------------------------------------------------------*/
void out_gid_submesh()
{
#ifndef NO_TEXT_OUTPUT
INT           i,j,k;
FILE         *out = allfiles.gidsubmesh;
FIELD        *actsmfield;     /* the submesh-field */
FIELD        *actmafield;     /* the "macro"-field */
GIDSET       *actsmgid;
ELEMENT      *actsmele;       /* the actual submesh-element */
ELEMENT      *firstsmele;     /* the first submesh-element */
INT           GlobalID;       /* pseudo global element number of submesh elements */

INT           is_firstmesh;

#ifdef DEBUG
dstrc_enter("out_gid_submesh");
#endif
/*----------------------------------------------------------------------*/
/*-------------------------------------------------------- print a head */
fprintf(out,"#--------------------------------------------------------------------------------\n");
fprintf(out,"# CCARAT postprocessing output to GiD - SUBMESH file\n");
fprintf(out,"# For postprocessing with gid, this file has to be copied to 'newprojectname.flavia.msh'\n");
fprintf(out,"#--------------------------------------------------------------------------------\n");
/*----------------------------------------------------------------------*/
actsmfield = &(sm_field[0]);
actmafield = &(field[0]);
actsmgid   = &(sm_gid[0]);
/*---------------------------------------- outpunt of nodal coordinates */
is_firstmesh=1;
/*----------------------------------------------------------------------*/
if (actsmgid->is_wall1_22)
{
   /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
   firstsmele = &(actsmfield->dis[0].element[0]);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# MESH %s FOR FIELD %s WALL1 2x2 GP\n",actsmgid->wall1_22_name,actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE %d \n",actsmgid->wall1_22_name,firstsmele->numnp);
   /*-------------- if this is first mesh, print coodinates of nodes */
   if (is_firstmesh)
   {
      is_firstmesh=0;
      fprintf(out,"# printout ALL nodal coordinates of ALL submeshes in first element type only\n");
      fprintf(out,"COORDINATES\n");
      out_gid_allsmcoords(out);
      fprintf(out,"END COORDINATES\n");
   }
   /*------------------------------------------------ print elements */
   fprintf(out,"ELEMENTS\n");
   /*--------------------------------------- loop over macroelements */
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     for (j=0; j<actsmfield->dis[0].numele; j++)
     {
       actsmele = &(actsmfield->dis[0].element[j]);
       if (actsmele->eltyp != el_wall1 ) continue;
       GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
       fprintf(out," %6d ",GlobalID);
       for (k=0; k<actsmele->numnp; k++)
          fprintf(out,"%6d ",i * actsmfield->dis[0].numnp + (actsmele->node[k]->Id+1));
       fprintf(out,"\n");
     }
   }
   fprintf(out,"END ELEMENTS\n");
}
/*----------------------------------------------------------------------*/
if (actsmgid->is_wall1_33)
{
   /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
   firstsmele = &(actsmfield->dis[0].element[0]);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# MESH %s FOR FIELD %s WALL1 3x3 GP\n",actsmgid->wall1_33_name,actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"MESH %s DIMENSION 3 ELEMTYPE Quadrilateral NNODE  %d \n",actsmgid->wall1_33_name,firstsmele->numnp);
   /*-------------- if this is first mesh, print coodinates of nodes */
   if (is_firstmesh)
   {
      is_firstmesh=0;
      fprintf(out,"# printout ALL nodal coordinates of ALL submeshes in first element type only\n");
      fprintf(out,"COORDINATES\n");
      out_gid_allsmcoords(out);
      fprintf(out,"END COORDINATES\n");
   }
   /*------------------------------------------------ print elements */
   fprintf(out,"ELEMENTS\n");
   /*--------------------------------------- loop over macroelements */
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     for (j=0; j<actsmfield->dis[0].numele; j++)
     {
       actsmele = &(actsmfield->dis[0].element[j]);
       if (actsmele->eltyp != el_wall1 ) continue;
       GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
       fprintf(out," %6d ",GlobalID);
       for (k=0; k<actsmele->numnp; k++)
          fprintf(out,"%6d ",i * actsmfield->dis[0].numnp + (actsmele->node[k]->Id+1));
       fprintf(out,"\n");
     }
   }
   fprintf(out,"END ELEMENTS\n");
}
/*----------------------------------------------------------------------*/
if (actsmgid->is_interf_22)
{
   firstsmele = &(actsmfield->dis[0].element[0]);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# MESH %s FOR FIELD %s INERFACE 2(+2) GP\n",actsmgid->interf_22_name,actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actsmgid->interf_22_name,firstsmele->numnp);
   /*-------------- if this is first mesh, print coodinates of nodes */
   if (is_firstmesh)
   {
      is_firstmesh=0;
      fprintf(out,"# printout ALL nodal coordinates of ALL submeshes in first element type only\n");
      fprintf(out,"COORDINATES\n");
      out_gid_allsmcoords(out);
      fprintf(out,"END COORDINATES\n");
   }
   /*------------------------------------------------ print elements */
   fprintf(out,"ELEMENTS\n");
   /*--------------------------------------- loop over macroelements */
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     for (j=0; j<actsmfield->dis[0].numele; j++)
     {
       actsmele = &(actsmfield->dis[0].element[j]);
       if (actsmele->eltyp != el_interf) continue;
       GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
       fprintf(out," %6d ",GlobalID);
       for (k=0; k<actsmele->numnp; k++)
          fprintf(out,"%6d ",i * actsmfield->dis[0].numnp + (actsmele->node[k]->Id+1));
       fprintf(out,"\n");
     }
   }
   fprintf(out,"END ELEMENTS\n");
}
/*----------------------------------------------------------------------*/
if (actsmgid->is_interf_33)
{
   firstsmele = &(actsmfield->dis[0].element[0]);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"# MESH %s FOR FIELD %s INERFACE 3(+2x3) GP\n",actsmgid->interf_33_name,actsmgid->fieldname);
   fprintf(out,"#-------------------------------------------------------------------------------\n");
   fprintf(out,"MESH %s DIMENSION 2 ELEMTYPE Quadrilateral NNODE  %d \n",actsmgid->interf_33_name,firstsmele->numnp);
   /*-------------- if this is first mesh, print coodinates of nodes */
   if (is_firstmesh)
   {
      is_firstmesh=0;
      fprintf(out,"# printout ALL nodal coordinates of ALL submeshes in first element type only\n");
      fprintf(out,"COORDINATES\n");
      out_gid_allsmcoords(out);
      fprintf(out,"END COORDINATES\n");
   }
   /*------------------------------------------------ print elements */
   fprintf(out,"ELEMENTS\n");
   /*--------------------------------------- loop over macroelements */
   for (i=0; i<actmafield->dis[0].numele; i++)
   {
     for (j=0; j<actsmfield->dis[0].numele; j++)
     {
       actsmele = &(actsmfield->dis[0].element[j]);
       if (actsmele->eltyp != el_interf) continue;
       GlobalID = i * actsmfield->dis[0].numele + (actsmele->Id + 1);
       fprintf(out," %6d ",GlobalID);
       for (k=0; k<actsmele->numnp; k++)
          fprintf(out,"%6d ",i * actsmfield->dis[0].numnp + (actsmele->node[k]->Id+1));
       fprintf(out,"\n");
     }
   }
   fprintf(out,"END ELEMENTS\n");
}
/*----------------------------------------------------------------------*/
fflush(out);
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_submesh */
#endif /* D_MLSTRUCT */

/*----------------------------------------------------------------------*
 |  routine to write all nodal coordinates of all fields m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_allcoords(FILE *out)
{
#ifndef NO_TEXT_OUTPUT
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
}/* end of (i=0; i<genprob.numfld; i++) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_allcoords */


#ifdef D_MLSTRUCT
/*----------------------------------------------------------------------*
 |  routine to write all nodal coordinates of all submeshes  ah 08/04   |
 *----------------------------------------------------------------------*/
void out_gid_allsmcoords(FILE *out)
{
#ifndef NO_TEXT_OUTPUT
INT           i,j,k;
FIELD        *actmafield;
FIELD        *actsmfield;
NODE         *actsmnode;
INT           ielema;          /* number of nodes per macroelement */
INT           GlobalID;        /* pseudo global node number of submesh nodes */
DOUBLE        translation[2];/* x/y Translations of submesh-prototyp
                              (prototyp has always x=0/y=0 at left down corner) */

#ifdef DEBUG
dstrc_enter("out_gid_allsmcoords");
#endif
/*----------------------------------------------------------------- */
actsmfield = &(sm_field[0]);
actmafield = &(field[0]);

ielema = actmafield->dis[0].element[0].numnp;
/*----------------------------------------- loop over macroelements */
for (i=0; i<actmafield->dis[0].numele; i++)
{
  translation[0] = actmafield->dis[0].element[i].node[0]->x[0];
  translation[1] = actmafield->dis[0].element[i].node[0]->x[1];
  /*----------------------- find coord. of nodes with lowest values */
  for (k=0; k<ielema; k++)
  {
    if (actmafield->dis[0].element[i].node[k]->x[0]<=translation[0] &&
        actmafield->dis[0].element[i].node[k]->x[1]<=translation[1])
    {
      translation[0] = actmafield->dis[0].element[i].node[k]->x[0];
      translation[1] = actmafield->dis[0].element[i].node[k]->x[1];
    }
  }
  /*---------------------- write NodeID and Nodecoordinates to file */
  for (j=0; j<actsmfield->dis[0].numnp; j++)
  {
     actsmnode = &(actsmfield->dis[0].node[j]);
     GlobalID  = i * actsmfield->dis[0].numnp + (actsmnode->Id + 1);
     fprintf(out,"%6d %-18.5f %-18.5f \n",
                                                  GlobalID,
                                                  (actsmnode->x[0]+translation[0]),
                                                  (actsmnode->x[1]+translation[1]));
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
#endif /* NO_TEXT_OUTPUT */
} /* end of out_gid_allsmcoords */

#endif /*D_MLSTRUCT*/

