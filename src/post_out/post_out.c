/*!
\file
\brief Postprocessing utility that takes ccarat output and produces
plain text.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

The general idea is that we cannot load the whole result data into
memory at once.

\author u.kue
\date 10/04

*/

#include "post_out.h"

static CHAR* UNDERLINE = "________________________________________________________________________________\n\n";
static CHAR* DBLLINE   = "================================================================================\n";


/*----------------------------------------------------------------------*/
/*!
  \brief Output general information.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_out_general(FILE* out, PROBLEM_DATA* problem)
{
  INT        i,j,k;

  INT numnp;
  INT numele;

#ifdef DEBUG
  dstrc_enter("out_general");
#endif

  /*-------------------------------------------------------- print header */
  fprintf(out, UNDERLINE);
  fprintf(out,"CCARAT outputfile\n");
  fprintf(out, UNDERLINE);

  /*------------------------------------------ print general problem data */
  fprintf(out, UNDERLINE);
  numnp = 0;
  numele = 0;
  for (i=0; i<problem->num_discr; ++i) {
    FIELD_DATA* field = &(problem->discr[i]);
    numnp += field->numnp;
    numele += field->numele;
  }
  fprintf(out,"Total number of Discretizations : %d\n", problem->num_discr);
  fprintf(out,"Total number of Elements        : %d\n", numele);
  fprintf(out,"Total number of Nodes           : %d\n", numnp);
  /*fprintf(out,"Total number of Materials       : %d\n", genprob.nmat);*/

  switch (problem->type) {
  case prb_fsi:
    fprintf(out,"Type of Problem                 : Fluid-Structure-Interaction\n");
    break;
  case prb_ssi:
    fprintf(out,"Type of Problem                 : Structure-Structure-Interaction\n");
    break;
  case prb_structure:
    fprintf(out,"Type of Problem                 : Structural\n");
    break;
  case prb_fluid:
    fprintf(out,"Type of Problem                 : Fluid\n");
    break;
  case prb_opt:
    fprintf(out,"Type of Problem                 : Optimization\n");
    break;
  case prb_ale:
    fprintf(out,"Type of Problem                 : Ale\n");
    break;
  case prb_twophase:
    fprintf(out,"Type of Problem                 : Two-Phase-Fluid-Flow\n");
    break;
  case prb_levelset:
    fprintf(out,"Type of Problem                 : Levelset \n");
    break;
  default:
    dserror("Cannot print problem type");
    break;
  }

  fprintf(out,UNDERLINE);

  for (i=0; i<problem->num_discr; ++i) {
    FIELD_DATA* field = &(problem->discr[i]);
    fprintf(out,DBLLINE);

    switch (field->type) {
    case fluid:
      fprintf(out,"FIELD: fluid\n");
      break;
    case ale:
      fprintf(out,"FIELD: ale\n");
      break;
    case structure:
      fprintf(out,"FIELD: structure\n");
      break;
    case levelset:
      fprintf(out,"FIELD: levelset\n");
      break;
    default:
      dserror("Cannot print fieldtype");
      break;
    }
    fprintf(out,DBLLINE);

    fprintf(out,UNDERLINE);
    fprintf(out,"Number of Elements  in this field : %d\n", field->numele);
    fprintf(out,"Number of Nodes     in this field : %d\n", field->numnp);
    fprintf(out,"Number of Dofs      in this field : %d\n", field->numdf);
    /*fprintf(out,"Number of Equations in this field : %d\n", field->numeq);*/
    fprintf(out,UNDERLINE);

#ifdef DEBUG

    {
      MAP* mesh_group;
      INT size_entry_length;
      INT size_offset;
      INT* mesh;

      MAP* coords_group;
      INT value_entry_length;
      INT value_offset;
      DOUBLE* coords;
      INT* node_ids;

      mesh_group = map_read_map(field->table, "mesh");
      size_entry_length = map_read_int(mesh_group, "size_entry_length");
      size_offset = map_read_int(mesh_group, "size_offset");

      /*
       * Have the whole mesh in memory at once: We are not going to
       * use this filter with very huge (gigantic) meshes anyway. */
      mesh = (INT*)CCACALLOC(size_entry_length*field->numele, sizeof(INT));

      fseek(field->size_file, size_offset, SEEK_SET);
      if (fread(mesh, sizeof(INT),
                size_entry_length*field->numele,
                field->size_file) != size_entry_length*field->numele) {
        dserror("reading mesh of discretization %s failed", field->name);
      }

      fprintf(out,"Element connectivity in global Ids:\n");
      for (j=0; j<field->numele; j++) {
        INT numnp;
        INT major;
        INT minor;
        INT* ele = &(mesh[j*size_entry_length]);

        major = ele[1];
        minor = ele[2];

        numnp = element_info[major].variant[minor].node_number;
        fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ", ele[0], numnp);
        for (k=0; k<numnp; k++)
          fprintf(out,"%6d ",ele[3+k]);
        fprintf(out,"\n");
      }

#if 0
      /*
       * To output these we'd have to search the coordinate array. Not
       * worth the trouble. */
      fprintf(out,"Element connectivity in field-local Ids:\n");
      for (j=0; j<actfield->dis[0].numele; j++) {
        actele = &(actfield->dis[0].element[j]);
        fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",actele->Id_loc,actele->numnp);
        for (k=0; k<actele->numnp; k++)
          fprintf(out,"%6d ",actele->node[k]->Id_loc);
        fprintf(out,"\n");
      }
#endif

      fprintf(out,UNDERLINE);
      fprintf(out,"Element types:\n");
      for (j=0; j<field->numele; j++) {
        INT major;
        INT minor;
        INT* ele = &(mesh[j*size_entry_length]);

        major = ele[1];
        minor = ele[2];

        switch (major) {
        case el_shell8:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL8\n",ele[0],j);
          break;
        case el_shell9:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL9\n",ele[0],j);
          break;
        case el_brick1:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d BRICK1\n",ele[0],j);
          break;
        case el_wall1:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d WALL1\n",ele[0],j);
          break;
        case el_fluid3:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3\n",ele[0],j);
          break;
        case el_fluid2:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2\n",ele[0],j);
          break;
        case el_fluid2_pro:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2_PRO\n",ele[0],j);
          break;
        case el_ale3:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE3\n",ele[0],j);
          break;
        case el_ale2:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE2\n",ele[0],j);
          break;
        case el_beam3:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d BEAM3\n",ele[0],j);
          break;
        case el_axishell:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d AXISHELL\n",ele[0],j);
          break;
        case el_interf:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d INTERFACE\n",ele[0],j);
          break;
        case el_wallge:
          fprintf(out,"ELE glob_Id %6d loc_Id %6d WALLGE\n",ele[0],j);
          break;
        default:
          dserror("Cannot print elementtype");
          break;
        }
      }

      fprintf(out,UNDERLINE);

      coords_group = map_read_map(field->table, "coords");

      value_entry_length = map_read_int(coords_group, "value_entry_length");
      value_offset = map_read_int(coords_group, "value_offset");

      size_entry_length = map_read_int(mesh_group, "size_entry_length");
      size_offset = map_read_int(mesh_group, "size_offset");

      /* Again, read everything at once. */
      coords = (DOUBLE*)CCACALLOC(value_entry_length*field->numnp, sizeof(DOUBLE));

      fseek(field->value_file, value_offset, SEEK_SET);
      if (fread(coords, sizeof(DOUBLE),
                value_entry_length*field->numnp,
                field->value_file) != value_entry_length*field->numnp) {
        dserror("reading node coordinates of discretization %s failed", field->name);
      }

      node_ids = (INT*)CCACALLOC(size_entry_length*field->numnp, sizeof(INT));

      fseek(field->size_file, size_offset, SEEK_SET);
      if (fread(node_ids, sizeof(INT),
                size_entry_length*field->numnp,
                field->size_file) != size_entry_length*field->numnp) {
        dserror("reading node ids of discretization %s failed", field->name);
      }

      fprintf(out,"Nodal Coordinates:\n");
      for (j=0; j<field->numnp; j++) {
        DOUBLE x[3];
        x[0] = coords[j*value_entry_length];
        x[1] = coords[j*value_entry_length+1];
        if (value_entry_length==3) {
          x[2] = coords[j*value_entry_length+2];
        }
        else {
          x[2] = 0;
        }
        fprintf(out,"NODE glob_Id %6d loc_Id %6d    % 18.5f % 18.5f % 18.5f \n",
                node_ids[j],j,x[0],x[1],x[2]);
      }

      fprintf(out,UNDERLINE);

#if 0
      /* No way to know them here. Right? */
      fprintf(out,"Degrees of Freedom:\n");
#endif

      CCAFREE(coords);
      CCAFREE(node_ids);
      CCAFREE(mesh);
    }
    fprintf(out,UNDERLINE);

#endif /*ifdef DEBUG */
  }

  fflush(out);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief The filter's main function.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  CHAR basename[100];
  CHAR filename[100];
  MAP control_table;
  PROBLEM_DATA problem;
  FILE* f;

  if (argc != 2) {
    printf("usage: %s control-file\n", argv[0]);
    return 1;
  }

  setup_filter(argv[1], &control_table, basename);

  dsassert(map_has_string(&control_table, "version", "0.1"),
           "expect version 0.1 control file");

  /* Debug output */
  /*map_print(stdout, &control_table, 0);*/

  init_problem_data(&problem, &control_table);

  sprintf(filename, "%s.out", basename);
  f = fopen(filename, "w");

  post_out_general(f, &problem);

  fclose(f);
  printf("Done.\n");
  return 0;
}
