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

Filters are independent programs, thus they have their own main
function. This filter's main function is the last function in this
file (to keep the number of needed function prototypes as small as
possible). You might want to start reading the file from there.

\author u.kue
\date 10/04

*/

#include "post_out.h"

#ifdef D_SHELL8
#include "../shell8/shell8.h"
#endif /*D_SHELL8*/
#ifdef D_SHELL9
#include "../shell9/shell9.h"
#endif /*D_SHELL9*/
#ifdef D_WALL1
#include "../wall1/wall1.h"
#endif /*D_WALL1*/
#ifdef D_BEAM3
#include "../beam3/beam3.h"
#endif /*D_BEAM3*/
#ifdef D_BRICK1
#include "../brick1/brick1.h"
#endif /*D_BRICK1*/
#ifdef D_INTERF
#include "../interf/interf.h"
#endif /*D_INTERF*/
#ifdef D_WALLGE
#include "../wallge/wallge.h"
#endif /*D_WALLGE*/


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
      INT* mesh;
      DOUBLE* coords;
      INT* node_ids;

      CHUNK_DATA mesh_chunk;
      CHUNK_DATA coords_chunk;

      if (!read_chunk_group(&mesh_chunk, field->table, "mesh")) {
        dserror("no mesh chunk found");
      }

      /*
       * Have the whole mesh in memory at once: We are not going to
       * use this filter with very huge (gigantic) meshes anyway. */
      mesh = (INT*)CCACALLOC(mesh_chunk.size_entry_length*field->numele, sizeof(INT));

      fseek(field->size_file, mesh_chunk.size_offset, SEEK_SET);
      if (fread(mesh, sizeof(INT),
                mesh_chunk.size_entry_length*field->numele,
                field->size_file) != mesh_chunk.size_entry_length*field->numele) {
        dserror("reading mesh of discretization %s failed", field->name);
      }

      fprintf(out,"Element connectivity in global Ids:\n");
      for (j=0; j<field->numele; j++) {
        INT numnp;
        INT major;
        INT minor;
        INT* ele = &(mesh[j*mesh_chunk.size_entry_length]);

        /* take the internal type numbers */
        major = field->element_type[j].major;
        minor = field->element_type[j].minor;

        numnp = element_info[major].variant[minor].node_number;
        fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ", ele[0], numnp);
        for (k=0; k<numnp; k++)
          fprintf(out,"%6d ",field->node_ids[ele[3+k]]);
        fprintf(out,"\n");
      }

      fprintf(out,UNDERLINE);
      fprintf(out,"Element connectivity in field-local Ids:\n");
      for (j=0; j<field->numele; j++) {
        INT numnp;
        INT major;
        INT minor;
        INT* ele = &(mesh[j*mesh_chunk.size_entry_length]);

        /* take the internal type numbers */
        major = field->element_type[j].major;
        minor = field->element_type[j].minor;

        numnp = element_info[major].variant[minor].node_number;
        fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",j,numnp);
        for (k=0; k<numnp; k++)
          fprintf(out,"%6d ",ele[3+k]);
        fprintf(out,"\n");
      }

      fprintf(out,UNDERLINE);
      fprintf(out,"Element types:\n");
      for (j=0; j<field->numele; j++) {
        INT major;
        INT minor;
        INT* ele = &(mesh[j*mesh_chunk.size_entry_length]);

        /* take the internal type numbers */
        major = field->element_type[j].major;
        minor = field->element_type[j].minor;

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

      if (!read_chunk_group(&coords_chunk, field->table, "coords")) {
        dserror("no coords chunk found");
      }

      /* Again, read everything at once. */
      coords = (DOUBLE*)CCACALLOC(coords_chunk.value_entry_length*field->numnp, sizeof(DOUBLE));

      fseek(field->value_file, coords_chunk.value_offset, SEEK_SET);
      if (fread(coords, sizeof(DOUBLE),
                coords_chunk.value_entry_length*field->numnp,
                field->value_file) != coords_chunk.value_entry_length*field->numnp) {
        dserror("reading node coordinates of discretization %s failed", field->name);
      }

      node_ids = (INT*)CCACALLOC(coords_chunk.size_entry_length*field->numnp, sizeof(INT));

      fseek(field->size_file, coords_chunk.size_offset, SEEK_SET);
      if (fread(node_ids, sizeof(INT),
                coords_chunk.size_entry_length*field->numnp,
                field->size_file) != coords_chunk.size_entry_length*field->numnp) {
        dserror("reading node ids of discretization %s failed", field->name);
      }

      fprintf(out,"Nodal Coordinates:\n");
      for (j=0; j<field->numnp; j++) {
        DOUBLE x[3];
        x[0] = coords[j*coords_chunk.value_entry_length];
        x[1] = coords[j*coords_chunk.value_entry_length+1];
        if (coords_chunk.value_entry_length==3) {
          x[2] = coords[j*coords_chunk.value_entry_length+2];
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
  \brief Print out solution of a certain step.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_out_sol(FILE* out,
                         PROBLEM_DATA* problem,
                         FIELD_DATA* field,
                         MAP* result_group)
{
  INT j;

  /* print header */
  fprintf(out, DBLLINE);
  switch (field->type){
  case fluid:
    fprintf(out,"FIELD: fluid\n");
    break;
  case ale:
    fprintf(out,"FIELD: ale\n");
    break;
  case structure:
    fprintf(out,"FIELD: structure\n");
    break;
  default:
    dserror("Cannot print fieldtype");
    break;
  }

  /* print nodal values */
  switch (field->type) {
  case structure: {
    CHUNK_DATA chunk;
    if (read_chunk_group(&chunk, result_group, "displacement")) {
      DOUBLE disp[3];

      if ((chunk.value_entry_length != 2) && (chunk.value_entry_length != 3)) {
        dserror("illegal displacement entry length %d", chunk.value_entry_length);
      }

      fprintf(out,DBLLINE);
      fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
              field->disnum, map_read_int(result_group, "step"));
      fprintf(out,DBLLINE);

      fseek(field->value_file, chunk.value_offset, SEEK_SET);
      for (j=0; j<field->numnp; j++) {
        INT k;

        if (fread(disp, sizeof(DOUBLE),
                  chunk.value_entry_length,
                  field->value_file) != chunk.value_entry_length) {
          dserror("reading node displacements of discretization %s failed", field->name);
        }

        fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",field->node_ids[j],j);
        for (k=0; k<chunk.value_entry_length; k++) {
          fprintf(out,"%20.7E ", disp[k]);
        }
        fprintf(out,"\n");
      }
      fprintf(out,UNDERLINE);
    }
    break;
  }
  case fluid: {
    CHUNK_DATA vel_chunk;
    CHUNK_DATA press_chunk;

    /* This needs to be changed if we want to read velocities witout
     * pressure. I put both together here because that's how it's done
     * in the old output routines. */
    if (read_chunk_group(&vel_chunk, result_group, "velocity") &&
        read_chunk_group(&press_chunk, result_group, "pressure")) {
      DOUBLE vel[3];
      DOUBLE pressure;

      if ((vel_chunk.value_entry_length != 2) && (vel_chunk.value_entry_length != 3)) {
        dserror("illegal velocity entry length %d", vel_chunk.value_entry_length);
      }

      if (press_chunk.value_entry_length != 1) {
        dserror("illegal pressure entry length %d", press_chunk.value_entry_length);
      }

      fprintf(out,DBLLINE);
      fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
              field->disnum, map_read_int(result_group, "step"));
      fprintf(out,DBLLINE);

      for (j=0; j<field->numnp; j++) {
        INT k;

        /* read velocity */
        /*
         * We need to fseek all the time because we are reading at two
         * places within one loop. */
        fseek(field->value_file,
              vel_chunk.value_offset + j*vel_chunk.value_entry_length*sizeof(DOUBLE),
              SEEK_SET);
        if (fread(vel, sizeof(DOUBLE),
                  vel_chunk.value_entry_length,
                  field->value_file) != vel_chunk.value_entry_length) {
          dserror("reading node velocity of discretization %s failed", field->name);
        }

        /* read pressure */
        fseek(field->value_file,
              press_chunk.value_offset + j*press_chunk.value_entry_length*sizeof(DOUBLE),
              SEEK_SET);
        if (fread(&pressure, sizeof(DOUBLE),
                  1,
                  field->value_file) != 1) {
          dserror("reading node pressure of discretization %s failed", field->name);
        }

        fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",field->node_ids[j],j);
        for (k=0; k<vel_chunk.value_entry_length; k++) {
          fprintf(out,"%20.7E ", vel[k]);
        }
        fprintf(out,"%20.7E ", pressure);
        fprintf(out,"\n");
      }
    }
    fprintf(out,UNDERLINE);
    break;
  }

    /* This is just like the structure case. We could merge both. */
  case ale: {
    CHUNK_DATA chunk;
    if (read_chunk_group(&chunk, result_group, "displacement")) {
      DOUBLE disp[3];

      if ((chunk.value_entry_length != 2) && (chunk.value_entry_length != 3)) {
        dserror("illegal displacement entry length %d", chunk.value_entry_length);
      }

      fprintf(out,DBLLINE);
      fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
              field->disnum, map_read_int(result_group, "step"));
      fprintf(out,DBLLINE);

      fseek(field->value_file, chunk.value_offset, SEEK_SET);
      for (j=0; j<field->numnp; j++) {
        INT k;

        if (fread(disp, sizeof(DOUBLE),
                  chunk.value_entry_length,
                  field->value_file) != chunk.value_entry_length) {
          dserror("reading node displacements of discretization %s failed", field->name);
        }

        fprintf(out,"NODE glob_Id %6d loc_Id %6d    ",field->node_ids[j],j);
        for (k=0; k<chunk.value_entry_length; k++) {
          fprintf(out,"%20.7E ", disp[k]);
        }
        fprintf(out,"\n");
      }
      fprintf(out,UNDERLINE);
    }
    break;
  }
  default:
    dserror("Cannot print fieldtype");
  }

  /* see if this discretization has a special type */

#ifdef D_SHELL8
  /* If there are shell8 elements it must be a discretization with
   * shell8 elements only. */
  if (map_symbol_count(field->table, "shell8_minor") > 0) {
    CHUNK_DATA chunk;

    /* This is interessting for the stress output. */
    if (read_chunk_group(&chunk, result_group, "stress")) {
      CHAR* forcetype[] = S8_FORCETYPE;
      CHAR* ft;
      INT i;
      INT force;
      INT minor;
      INT ngauss;
      DOUBLE* stress;

      minor = map_read_int(field->table, "shell8_minor");
      ngauss = element_info[el_shell8].variant[minor].gauss_number;

      /* find the force type */
      ft = map_read_string(field->table, "shell8_forcetype");
      for (force=0; forcetype[force]!=NULL; ++force) {
        if (strcmp(ft, forcetype[force])==0) {
          break;
        }
      }
      if (forcetype[force]==NULL) {
        dserror("unknown force type '%s'", ft);
      }

      /* the memory for one chunk entry */
      stress = (DOUBLE*)CCACALLOC(chunk.value_entry_length, sizeof(DOUBLE));

      /* loop all elements */
      fseek(field->value_file, chunk.value_offset, SEEK_SET);
      for (j=0; j<field->numele; j++) {

        if (fread(stress, sizeof(DOUBLE),
                  chunk.value_entry_length,
                  field->value_file) != chunk.value_entry_length) {
          dserror("reading shell8 stesses of discretization %s failed", field->name);
        }

        fprintf(out,UNDERLINE);
        fprintf(out,"Element glob_Id %d loc_Id %d                SHELL8\n",field->element_type[j].Id,j);
        switch (force) {
        case s8_xyz:
          fprintf(out,"Gaussian     Force-xx     Force-xy     Force-yx      Force-yy     Force-xz     Force-zx      Force-yz    Force-zy     Force-zz\n");
          break;
        case s8_rst:
          fprintf(out,"Gaussian     Force-rr     Force-rs     Force-sr      Force-ss     Force-rt     Force-tr      Force-st    Force-ts     Force-tt\n");
          break;
        case s8_rst_ortho:
          fprintf(out,"Gaussian     Force-rr     Force-rs     Force-sr      Force-ss     Force-rt     Force-tr      Force-st    Force-ts     Force-tt\n");
          break;
        default:
          dserror("Unknown type of element stresses");
        }

        for (i=0; i<ngauss; i++) {
          /* There are 18 values to one gauss point */
          fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                  i,
                  stress[ 0+18*i],
                  stress[ 2+18*i],
                  stress[ 8+18*i],
                  stress[ 1+18*i],
                  stress[ 3+18*i],
                  stress[16+18*i],
                  stress[ 4+18*i],
                  stress[17+18*i],
                  stress[ 9+18*i]
            );
        }

        switch (force) {
        case s8_xyz:
          fprintf(out,"Gaussian     Moment-xx    Moment-xy    Moment-yx     Moment-yy    Moment-xz    Moment-zx     Moment-yz    Moment-zy    Moment-zz\n");
          break;
        case s8_rst:
          fprintf(out,"Gaussian     Moment-rr    Moment-rs    Moment-sr     Moment-ss    Moment-rt    Moment-tr     Moment-st    Moment-ts    Moment-tt\n");
          break;
        case s8_rst_ortho:
          fprintf(out,"Gaussian     Moment-rr    Moment-rs    Moment-sr     Moment-ss    Moment-rt    Moment-tr     Moment-st    Moment-ts    Moment-tt\n");
          break;
        default:
          dserror("Unknown type of element stresses");
        }

        for (i=0; i<ngauss; i++) {
          /* There are 18 values to one gauss point */
          fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                  i,
                  stress[ 5+18*i],
                  stress[ 7+18*i],
                  stress[14+18*i],
                  stress[ 6+18*i],
                  stress[10+18*i],
                  stress[12+18*i],
                  stress[11+18*i],
                  stress[13+18*i],
                  stress[15+18*i]
            );
        }
      }

      CCAFREE(stress);
    }

    goto end;
  }
#endif

#ifdef D_SHELL9
  /* If there are shell9 elements it must be a discretization with
   * shell9 elements only. */
  if (map_symbol_count(field->table, "shell9_minor") > 0) {

/*NOTE: It does not seem to be very interesting to write the forces of each kinematic layer; it is more
  usefull to look at the stresses. These are writen to the flavia.res-File, so that they could be
  visualized with gid.   sh 02/03 */

    goto end;
  }
#endif

  /* print element values */
  {
    CHUNK_DATA chunk;

    if (read_chunk_group(&chunk, result_group, "stress")) {
      DOUBLE* stress;

      /* the memory for one chunk entry */
      stress = (DOUBLE*)CCACALLOC(chunk.value_entry_length, sizeof(DOUBLE));

      /* loop all elements */
      fseek(field->value_file, chunk.value_offset, SEEK_SET);
      for (j=0; j<field->numele; j++) {
        INT major;
        INT minor;
        INT ngauss;

        major = field->element_type[j].major;
        minor = field->element_type[j].minor;

        /*
         * Read the current stress entry. Independent of the element type. */
        if (fread(stress, sizeof(DOUBLE),
                  chunk.value_entry_length,
                  field->value_file) != chunk.value_entry_length) {
          dserror("reading stesses of discretization %s failed", field->name);
        }

        switch (major) {
#ifdef D_WALL1
        case el_wall1: {
          INT i;

          ngauss = element_info[el_wall1].variant[minor].gauss_number;
          dsassert(ngauss == 4, "only four gauss point versions supported currently");

          fprintf(out,DBLLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n",field->element_type[j].Id,j);
          fprintf(out,"\n");

          /* check whether stresses at Gauss Points are presented in
           * global xy- or local rs-coordinate system and write stress type */
          if ((map_symbol_count(chunk.group, "wall1_stresstype")>0) &&
              (strcmp(map_read_string(chunk.group, "wall1_stresstype"), "w1_rs")==0)) {
            fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
          }
          else {
            fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
          }

          for (i=0; i<ngauss; i++) {
            /*
             * This is a special case. Four stress values, three other
             * values to be output here and two values to be output
             * using GiD. The wall element needs to be cleaned up! */
            fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                    i,
                    stress[0*9+i],
                    stress[1*9+i],
                    stress[2*9+i],
                    stress[3*9+i],
                    stress[4*9+i],
                    stress[5*9+i],
                    stress[6*9+i]
              );
          }
          break;
        }
#endif /*D_WALL1*/

#ifdef D_BRICK1
        case el_brick1: {
          INT i;
          CHAR* stresstype[] = BRICK1_STRESSTYPE;
          CHAR* st;
          INT stress_type;

          ngauss = element_info[el_brick1].variant[minor].gauss_number;

          /* Find the stress type's number. */
          st = map_read_string(chunk.group, "brick1_stresstype");
          for (stress_type=0; stresstype[stress_type]!=NULL; ++stress_type) {
            if (strcmp(stresstype[stress_type], st)==0) {
              break;
            }
          }

          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                BRICK1\n",field->element_type[j].Id,j);
          fprintf(out,"\n");

          switch (stress_type) {
          case c1_gpxyz:
            fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-xx    stress-yy    stress-zz    stress-xy    stress-xz    stress-yz\n");
            for (i=0; i<ngauss; i++) {
              /* There are 27 values to one gauss point */
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                      i,
                      stress[24+27*i],
                      stress[25+27*i],
                      stress[26+27*i],
                      stress[ 6+27*i],
                      stress[ 7+27*i],
                      stress[ 8+27*i],
                      stress[ 9+27*i],
                      stress[10+27*i],
                      stress[11+27*i]
                );
            }
            break;
          case c1_gprst:
            fprintf(out,"r,s,t    ---> local system on element level \n");
            fprintf(out,"rr,ss,tt ---> normal-stresses               \n");
            fprintf(out,"rs,st,tr ---> shear -stresses               \n\n");
            fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
            for (i=0; i<ngauss; i++) {
              /* There are 27 values to one gauss point */
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                      i,
                      stress[24+27*i],
                      stress[25+27*i],
                      stress[26+27*i],
                      stress[ 0+27*i],
                      stress[ 1+27*i],
                      stress[ 2+27*i],
                      stress[ 3+27*i],
                      stress[ 4+27*i],
                      stress[ 5+27*i]
                );
            }
            break;
          case c1_gp123:
            fprintf(out,"11,22,33 ---> principal-stresses                       \n");
            fprintf(out,"r1,s1,t1 ---> angles to the first  principal direction \n");
            fprintf(out,"r2,s2,t2 ---> angles to the second principal direction \n");
            fprintf(out,"r3,s3,t3 ---> angles to the third  principal direction \n\n");
            fprintf(out,"INT.point   stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
            for (i=0; i<ngauss; i++) {
              /* There are 27 values to one gauss point */
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                      i,
                      stress[12+27*i],
                      stress[13+27*i],
                      stress[14+27*i],
                      stress[15+27*i],
                      stress[16+27*i],
                      stress[17+27*i],
                      stress[18+27*i],
                      stress[19+27*i],
                      stress[20+27*i],
                      stress[21+27*i],
                      stress[22+27*i],
                      stress[23+27*i]
                );
            }
            break;
#if 0
            /* Nodal stresses are currently not available in the binary output. */
          case c1_nprst:
            fprintf(out,"elenode     stress-rr    stress-ss    stress-tt    stress-rs    stress-st    stress-tr\n");
            for (i=0; i<actele->numnp; i++)
            {
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                      i,
                      actele->e.c1->stress_ND.a.d3[place][0][i],
                      actele->e.c1->stress_ND.a.d3[place][1][i],
                      actele->e.c1->stress_ND.a.d3[place][2][i] ,
                      actele->e.c1->stress_ND.a.d3[place][3][i],
                      actele->e.c1->stress_ND.a.d3[place][4][i],
                      actele->e.c1->stress_ND.a.d3[place][5][i]
                );
            }
            break;
          case c1_np123:
            fprintf(out,"elenode     stress-11    stress-22    stress-33  ang-r1  ang-s1   ang-t1    ang-r2   ang-s2   ang-t2   ang-r3   ang-s3   ang-t3\n");
            for (i=0; i<actele->numnp; i++)
            {
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f    %5.2f \n",
                      i,
                      actele->e.c1->stress_ND.a.d3[place][12][i],
                      actele->e.c1->stress_ND.a.d3[place][13][i],
                      actele->e.c1->stress_ND.a.d3[place][14][i] ,
                      actele->e.c1->stress_ND.a.d3[place][15][i],
                      actele->e.c1->stress_ND.a.d3[place][16][i],
                      actele->e.c1->stress_ND.a.d3[place][17][i],
                      actele->e.c1->stress_ND.a.d3[place][18][i],
                      actele->e.c1->stress_ND.a.d3[place][19][i],
                      actele->e.c1->stress_ND.a.d3[place][20][i],
                      actele->e.c1->stress_ND.a.d3[place][21][i],
                      actele->e.c1->stress_ND.a.d3[place][22][i],
                      actele->e.c1->stress_ND.a.d3[place][23][i]
                );
            }
            break;
          case c1_npxyz:
            fprintf(out,"elenode     stress-xx    stress-yy    stress-zz    stress-xy    stress-yz    stress-xz\n");
            for (i=0; i<actele->numnp; i++)
            {
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                      i,
                      actele->e.c1->stress_ND.a.d3[place][ 6][i],
                      actele->e.c1->stress_ND.a.d3[place][ 7][i],
                      actele->e.c1->stress_ND.a.d3[place][ 8][i] ,
                      actele->e.c1->stress_ND.a.d3[place][ 9][i],
                      actele->e.c1->stress_ND.a.d3[place][10][i],
                      actele->e.c1->stress_ND.a.d3[place][11][i]
                );
            }
            break;
#endif
          default:
            fprintf(out,"no stresses available\n");
          }
          break;
        }
#endif /*D_BRICK1*/

#ifdef D_BEAM3
        case el_beam3: {
          INT i;
          ngauss = element_info[el_beam3].variant[minor].gauss_number;
          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                BEAM3\n",field->element_type[j].Id,j);
          fprintf(out,"\n");
          fprintf(out,"Gaussian         Nx           Vy           Vz           Mx           My           Mz\n");
          for (i=0; i<ngauss; i++) {
            /* There are 6 values to one gauss point */
            fprintf(out,"Gauss %d       %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                    i,
                    stress[0+6*i],
                    stress[1+6*i],
                    stress[2+6*i],
                    stress[3+6*i],
                    stress[4+6*i],
                    stress[5+6*i]
              );
          }
          break;
        }
#endif /*D_BEAM3*/

#ifdef D_INTERF
        case el_interf: {
          INT i;
          ngauss = element_info[el_interf].variant[minor].gauss_number;
          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                INTERF\n",field->element_type[j].Id,j);

          if ((map_symbol_count(chunk.group, "interf_orient")>0) &&
              (strcmp(map_read_string(chunk.group, "interf_orient"), "local")==0)) {
            fprintf(out,"Gaussian     stresses-tangential     stresses-normal\n");
          }
          else {
            fprintf(out,"Gaussian     stresses-xx     stresses-yy    stresses-xy\n");
          }

          for (i=0; i<ngauss; i++) {
            /* There are 5 values to one gauss point */
            fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E\n",
                    i,
                    stress[0+5*i],
                    stress[1+5*i],
                    stress[2+5*i] );
          }
          break;
        }
#endif /*D_INTERF*/

#ifdef D_WALLGE
        case el_wallge: {
          INT i;
          ngauss = element_info[el_wallge].variant[minor].gauss_number;
          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n",field->element_type[j].Id,j);
          fprintf(out,"\n");

          /* check wether stresses at Gauss Points are presented in
           * global xy- or local rs-coordinate system and write stress type */
          if ((map_symbol_count(chunk.group, "wallge_stresstype")>0) &&
              (strcmp(map_read_string(chunk.group, "wallge_stresstype"), "wge_rs")==0)) {
            fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
          }
          else {
            fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
          }

          for (i=0; i<ngauss; i++) {
#if 0
            /* These stresses are never exported! Do we need them? */
            fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                    i,
                    actele->e.wallge->stress_GP.a.d3[place][0][i],
                    actele->e.wallge->stress_GP.a.d3[place][1][i],
                    actele->e.wallge->stress_GP.a.d3[place][2][i],
                    actele->e.wallge->stress_GP.a.d3[place][3][i],
                    actele->e.wallge->stress_GP.a.d3[place][4][i],
                    actele->e.wallge->stress_GP.a.d3[place][5][i],
                    actele->e.wallge->stress_GP.a.d3[place][6][i]
              );
#endif
          }
          break;
        }
#endif /*D_WALLGE*/

        default:
          dserror("unknown type of element");
          break;
        }

      }
    }
  }

end:

  fprintf(out,"\n");
  fprintf(out,"\n");
}


#ifdef D_FSI

/*----------------------------------------------------------------------*/
/*!
  \brief Print fsi coupling informations.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void post_out_fsi(FILE* out, PROBLEM_DATA* problem)
{
  INT        i;
  FIELD_DATA* struct_field = NULL;
  FIELD_DATA* fluid_field;
  FIELD_DATA* ale_field;

  INT* fluid_struct_connect;
  INT* fluid_ale_connect;

  fprintf(out,DBLLINE);
  fprintf(out,"FSI node connectivity global Ids:\n");
  fprintf(out,DBLLINE);
  fprintf(out,"\n");
  fprintf(out,"FLUID          ALE          STRUCTURE\n");

  /* Find the corresponding discretizations. We don't rely on any order. */
  for (i=0; i<problem->num_discr; ++i) {
    if (problem->discr[i].type == structure) {
      struct_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == fluid) {
      fluid_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == ale) {
      ale_field = &(problem->discr[i]);
    }
    else {
      dserror("unknown field type %d", problem->discr[i].type);
    }
  }

  post_find_fsi_coupling(problem,
                         struct_field, fluid_field, ale_field,
                         &fluid_struct_connect, &fluid_ale_connect);

  for (i=0; i<fluid_field->numnp; i++) {
    if ((fluid_ale_connect[i]!=-1) && (fluid_struct_connect[i]==-1)) {
      fprintf(out,"%-6d         %-6d       ------\n",
              fluid_field->node_ids[i],
              ale_field->node_ids[fluid_ale_connect[i]]);
    }
    else if ((fluid_ale_connect[i]==-1) && (fluid_struct_connect[i]!=-1)) {
      fprintf(out,"%-6d         ------       %-6d\n",
              fluid_field->node_ids[i],
              struct_field->node_ids[fluid_struct_connect[i]]);
    }
    else if ((fluid_ale_connect[i]!=-1) && (fluid_struct_connect[i]!=-1)) {
      fprintf(out,"%-6d         %-6d       %-6d\n",
              fluid_field->node_ids[i],
              ale_field->node_ids[fluid_ale_connect[i]],
              struct_field->node_ids[fluid_struct_connect[i]]);
    }
    else
      fprintf(out,"%-6d         ------       ------\n",
              fluid_field->node_ids[i]);

  }
  fprintf(out,UNDERLINE);

#if 0
  /* No worth it (?) */
  fprintf(out,DBLLINE);
  fprintf(out,"FSI node connectivity local Ids:\n");
  fprintf(out,DBLLINE);
  fprintf(out,"\n");
  fprintf(out,"FLUID          ALE          STRUCTURE\n");
  for (i=0;i<fluidfield->dis[0].numnp;i++)
  {
    actfnode  = &(fluidfield->dis[0].node[i]);
    actfgnode = actfnode->gnode;
    actsnode  = actfgnode->mfcpnode[numsf];
    actanode  = actfgnode->mfcpnode[numaf];
    if (actsnode==NULL && actanode!=NULL)
      fprintf(out,"%-6d         %-6d       ------\n",actfnode->Id_loc,actanode->Id_loc);
    else if (actanode==NULL && actsnode!=NULL)
      fprintf(out,"%-6d         ------       %-6d\n",actfnode->Id_loc,actsnode->Id_loc);
    else if (actanode!=NULL && actsnode!=NULL)
      fprintf(out,"%-6d         %-6d       %-6d\n",actfnode->Id_loc,actanode->Id,actsnode->Id_loc);
    else
      fprintf(out,"%-6d         ------       ------\n",actfnode->Id_loc);

  }
  fprintf(out,UNDERLINE);
#endif
}

#endif


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
  INT i;

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
  f = fopen(filename, "wb");

  post_out_general(f, &problem);

#ifdef D_FSI
  /* Output the coupling information if it's a fsi problem. It's a
   * little harder to discover a fluid-ale problem, but the function
   * is supposed to work the same way. Maybe we could find a new
   * problemtype for fluid-ale? */
  if (strcmp(map_read_string(&control_table, "problem_type"), "fsi")==0) {
    post_out_fsi(f, &problem);
  }
#endif

  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i) {
    INT res;
    FIELD_DATA* field;
    field = &(problem.discr[i]);

    /* Iterate all results. */
    for (res=0; res<problem.num_results; ++res) {
      MAP* result_group;
      result_group = problem.result_group[res];

      /* We iterate the list of all results. Here we are interested in
       * the results of this discretization. */
      if (match_field_result(field, result_group)) {
        post_out_sol(f, &problem, field, result_group);
      }
    }
  }

  fclose(f);
  fprintf(allfiles.out_err, "Done.\n");
  return 0;
}
