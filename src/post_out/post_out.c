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
  for (i=0; i<problem->num_discr; ++i)
  {
    FIELD_DATA* field = &(problem->discr[i]);
    numnp += field->numnp;
    numele += field->numele;
  }
  fprintf(out,"Total number of Discretizations : %d\n", problem->num_discr);
  fprintf(out,"Total number of Elements        : %d\n", numele);
  fprintf(out,"Total number of Nodes           : %d\n", numnp);
  /*fprintf(out,"Total number of Materials       : %d\n", genprob.nmat);*/

  switch (problem->type)
  {
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

  /* Iterate all discretizations */

  for (i=0; i<problem->num_discr; ++i)
  {
    FIELD_DATA* field = &(problem->discr[i]);
    fprintf(out,DBLLINE);

    switch (field->type)
    {
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

    fprintf(out,"Element connectivity in global Ids:\n");
    for (j=0; j<field->numele; j++)
    {
      INT numnp;
      INT Id;

      /* read the element's data */
      chunk_read_size_entry(&(field->ele_param), j);
      chunk_read_size_entry(&(field->mesh), j);

      Id    = field->ele_param.size_buf[element_variables.ep_size_Id];
      numnp = field->ele_param.size_buf[element_variables.ep_size_numnp];

      fprintf(out,"glob_Id %6d Nnodes %2d Nodes: ", Id, numnp);
      for (k=0; k<numnp; k++)
      {
        /* read the global node ids */
        chunk_read_size_entry(&(field->coords), field->mesh.size_buf[k]);
        fprintf(out,"%6d ",field->coords.size_buf[node_variables.coords_size_Id]);
      }
      fprintf(out,"\n");
    }

    fprintf(out,UNDERLINE);
    fprintf(out,"Element connectivity in field-local Ids:\n");
    for (j=0; j<field->numele; j++)
    {
      INT numnp;

      /* read the element's data */
      chunk_read_size_entry(&(field->ele_param), j);
      chunk_read_size_entry(&(field->mesh), j);

      numnp = field->ele_param.size_buf[element_variables.ep_size_numnp];

      fprintf(out,"loc_Id %6d Nnodes %2d Nodes: ",j,numnp);
      for (k=0; k<numnp; k++)
        fprintf(out,"%6d ", field->mesh.size_buf[k]);
      fprintf(out,"\n");
    }

    fprintf(out,UNDERLINE);
    fprintf(out,"Element types:\n");
    for (j=0; j<field->numele; j++)
    {
      INT el_type;
      INT Id;

      /* read the element's data */
      chunk_read_size_entry(&(field->ele_param), j);

      Id      = field->ele_param.size_buf[element_variables.ep_size_Id];
      el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];

      /* translate to internal value */
      el_type = problem->element_type.table[el_type];

      switch (el_type)
      {
      case el_shell8:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL8\n",Id,j);
        break;
      case el_shell9:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d SHELL9\n",Id,j);
        break;
      case el_brick1:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d BRICK1\n",Id,j);
        break;
      case el_wall1:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d WALL1\n",Id,j);
        break;
      case el_fluid3:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID3\n",Id,j);
        break;
      case el_fluid2:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2\n",Id,j);
        break;
      case el_fluid2_pro:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d FLUID2_PRO\n",Id,j);
        break;
      case el_ale3:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE3\n",Id,j);
        break;
      case el_ale2:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d ALE2\n",Id,j);
        break;
      case el_beam3:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d BEAM3\n",Id,j);
        break;
      case el_axishell:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d AXISHELL\n",Id,j);
        break;
      case el_interf:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d INTERFACE\n",Id,j);
        break;
      case el_wallge:
        fprintf(out,"ELE glob_Id %6d loc_Id %6d WALLGE\n",Id,j);
        break;
      default:
        dserror("Cannot print elementtype");
        break;
      }
    }

    fprintf(out,UNDERLINE);

    fprintf(out,"Nodal Coordinates:\n");
    for (j=0; j<field->numnp; j++)
    {
      DOUBLE x[3];

      chunk_read_size_entry(&(field->coords), j);
      chunk_read_value_entry(&(field->coords), j);

      x[0] = field->coords.value_buf[0];
      x[1] = field->coords.value_buf[1];
      if (field->coords.value_entry_length==3)
      {
        x[2] = field->coords.value_buf[2];
      }
      else
      {
        x[2] = 0;
      }
      fprintf(out,"NODE glob_Id %6d loc_Id %6d    % 18.5f % 18.5f % 18.5f \n",
              field->coords.size_buf[node_variables.coords_size_Id],
              j,x[0],x[1],x[2]);
    }

    fprintf(out,UNDERLINE);

#if 0
    /* No way to know them here. Right? */
    fprintf(out,"Degrees of Freedom:\n");
#endif

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
static void post_out_sol(FILE* out, RESULT_DATA* result)
{
  INT j;
  FIELD_DATA* field;

#ifdef DEBUG
  dstrc_enter("post_out_sol");
#endif

  field = result->field;

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
  switch (field->type)
  {
  case structure:
  {
    CHUNK_DATA chunk;

    init_chunk_data(result, &chunk, "displacement");

    if ((chunk.value_entry_length != 2) && (chunk.value_entry_length != 3))
    {
      dserror("illegal displacement entry length %d", chunk.value_entry_length);
    }

    fprintf(out,DBLLINE);
    fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
            result->field->disnum, map_read_int(result->group, "step"));

    fprintf(out,DBLLINE);

    for (j=0; j<field->numnp; j++)
    {
      INT k;
      INT Id;

      /* again the global node id is needed */
      chunk_read_size_entry(&(field->coords), j);
      Id = field->coords.size_buf[node_variables.coords_size_Id];

      /* read the displacements to this node */
      chunk_read_value_entry(&chunk, j);

      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ", Id, j);
      for (k=0; k<chunk.value_entry_length; k++)
      {
        fprintf(out,"%20.7E ", chunk.value_buf[k]);
      }
      fprintf(out,"\n");
    }
    fprintf(out,UNDERLINE);

    destroy_chunk_data(&chunk);
    break;
  }
  case fluid:
  {
    CHUNK_DATA vel_chunk;
    CHUNK_DATA press_chunk;

    init_chunk_data(result, &vel_chunk, "velocity");
    init_chunk_data(result, &press_chunk, "pressure");

    if ((vel_chunk.value_entry_length != 2) && (vel_chunk.value_entry_length != 3))
    {
      dserror("illegal velocity entry length %d", vel_chunk.value_entry_length);
    }

    if (press_chunk.value_entry_length != 1)
    {
      dserror("illegal pressure entry length %d", press_chunk.value_entry_length);
    }

    fprintf(out,DBLLINE);
    fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
            field->disnum, map_read_int(result->group, "step"));
    fprintf(out,DBLLINE);

    for (j=0; j<field->numnp; j++)
    {
      INT k;
      INT Id;

      /* again the global node id is needed */
      chunk_read_size_entry(&(field->coords), j);
      Id = field->coords.size_buf[node_variables.coords_size_Id];

      /* read velocity */
      chunk_read_value_entry(&vel_chunk, j);

      /* read pressure */
      chunk_read_value_entry(&press_chunk, j);

      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ", Id, j);
      for (k=0; k<vel_chunk.value_entry_length; k++)
      {
        fprintf(out,"%20.7E ", vel_chunk.value_buf[k]);
      }
      fprintf(out,"%20.7E ", press_chunk.value_buf[0]);
      fprintf(out,"\n");
    }
    fprintf(out,UNDERLINE);

    destroy_chunk_data(&vel_chunk);
    destroy_chunk_data(&press_chunk);
    break;
  }

    /* This is just like the structure case. We could merge both. */
  case ale: {
    CHUNK_DATA chunk;

    init_chunk_data(result, &chunk, "displacement");

    if ((chunk.value_entry_length != 2) && (chunk.value_entry_length != 3))
    {
      dserror("illegal displacement entry length %d", chunk.value_entry_length);
    }

    fprintf(out,DBLLINE);
    fprintf(out,"Converged Solution of Discretisation %d in step %d \n",
            field->disnum, map_read_int(result->group, "step"));
    fprintf(out,DBLLINE);

    for (j=0; j<field->numnp; j++) {
      INT k;
      INT Id;

      /* again the global node id is needed */
      chunk_read_size_entry(&(field->coords), j);
      Id = field->coords.size_buf[node_variables.coords_size_Id];

      /* read the displacements to this node */
      chunk_read_value_entry(&chunk, j);

      fprintf(out,"NODE glob_Id %6d loc_Id %6d    ", Id, j);
      for (k=0; k<chunk.value_entry_length; k++) {
        fprintf(out,"%20.7E ", chunk.value_buf[k]);
      }
      fprintf(out,"\n");
    }
    fprintf(out,UNDERLINE);

    destroy_chunk_data(&chunk);
    break;
  }
  default:
    dserror("Cannot print fieldtype");
  }

  /* see if this discretization has a special type */

#ifdef D_SHELL8
  /* If there are shell8 elements it must be a discretization with
   * shell8 elements only. */
  if (map_has_string(field->group, "shell8_problem", "yes"))
  {
    CHUNK_DATA chunk;

    /* This is interessting for the stress output. */
    if (map_has_map(result->group, "stress"))
    {
      CHAR* forcetypenames[] = S8_FORCETYPE;
      INT i;
      TRANSLATION_TABLE forcetype;

      /* We want to interpret the forces in the right way. */
      init_translation_table(&forcetype,
                             map_read_map(field->group, "shell8_forcetype"),
                             forcetypenames);

      /* The chunk we are going to read. */
      init_chunk_data(result, &chunk, "stress");

      /* loop all elements */
      for (j=0; j<field->numele; j++)
      {
        INT Id;
        INT nGP[3];
        INT ngauss;
        INT force;

        /* we need the element parameters here */
        chunk_read_size_entry(&(field->ele_param), j);

        Id     = field->ele_param.size_buf[element_variables.ep_size_Id];

        nGP[0] = field->ele_param.size_buf[shell8_variables.ep_size_nGP0];
        nGP[1] = field->ele_param.size_buf[shell8_variables.ep_size_nGP1];
        nGP[2] = field->ele_param.size_buf[shell8_variables.ep_size_nGP2];
        ngauss = nGP[0]*nGP[1]*nGP[2];

        force  = field->ele_param.size_buf[shell8_variables.ep_size_forcetyp];

        /* translate to internal */
        force = forcetype.table[force];

        /* finally read the stresses */
        chunk_read_value_entry(&chunk, j);

        fprintf(out,UNDERLINE);
        fprintf(out,"Element glob_Id %d loc_Id %d                SHELL8\n", Id, j);
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
                  chunk.value_buf[ 0+18*i],
                  chunk.value_buf[ 2+18*i],
                  chunk.value_buf[ 8+18*i],
                  chunk.value_buf[ 1+18*i],
                  chunk.value_buf[ 3+18*i],
                  chunk.value_buf[16+18*i],
                  chunk.value_buf[ 4+18*i],
                  chunk.value_buf[17+18*i],
                  chunk.value_buf[ 9+18*i]
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
                  chunk.value_buf[ 5+18*i],
                  chunk.value_buf[ 7+18*i],
                  chunk.value_buf[14+18*i],
                  chunk.value_buf[ 6+18*i],
                  chunk.value_buf[10+18*i],
                  chunk.value_buf[12+18*i],
                  chunk.value_buf[11+18*i],
                  chunk.value_buf[13+18*i],
                  chunk.value_buf[15+18*i]
            );
        }
      }

      destroy_chunk_data(&chunk);
      destroy_translation_table(&forcetype);
    }

    goto end;
  }
#endif

#ifdef D_SHELL9
  /* If there are shell9 elements it must be a discretization with
   * shell9 elements only. */
  if (map_has_string(field->group, "shell9_problem", "yes"))
  {

/*NOTE: It does not seem to be very interesting to write the forces of each kinematic layer; it is more
  usefull to look at the stresses. These are writen to the flavia.res-File, so that they could be
  visualized with gid.   sh 02/03 */

    goto end;
  }
#endif

  /* print element values */
  {
    CHUNK_DATA chunk;

    if (map_has_map(result->group, "stress"))
    {
#ifdef D_WALL1
      CHAR* w1_stresstypenames[] = WALL1_STRESSTYPE;
      TRANSLATION_TABLE w1_stresstype;
#endif
#ifdef D_BRICK1
      CHAR* c1_stresstypenames[] = BRICK1_STRESSTYPE;
      TRANSLATION_TABLE c1_stresstype;
#endif
#ifdef D_INTERF
      CHAR* interf_stresstypenames[] = INTERF_STRESSTYPE;
      TRANSLATION_TABLE interf_stresstype;
#endif
#ifdef D_WALLGE
      CHAR* wallge_stresstypenames[] = WALLGE_STRESSTYPE;
      TRANSLATION_TABLE wallge_stresstype;
#endif


#ifdef D_WALL1
      init_translation_table(&w1_stresstype,
                             map_read_map(field->group, "w1_stresstypes"),
                             w1_stresstypenames);
#endif
#ifdef D_BRICK1
      init_translation_table(&c1_stresstype,
                             map_read_map(field->group, "c1_stresstypes"),
                             c1_stresstypenames);
#endif
#ifdef D_INTERF
      init_translation_table(&interf_stresstype,
                             map_read_map(field->group, "interf_stresstypes"),
                             interf_stresstypenames);
#endif
#ifdef D_WALLGE
      init_translation_table(&wallge_stresstype,
                             map_read_map(field->group, "wallge_stresstypes"),
                             wallge_stresstypenames);
#endif


      /* The chunk we are going to read. */
      init_chunk_data(result, &chunk, "stress");

      /* loop all elements */
      for (j=0; j<field->numele; j++)
      {
        INT Id;
        INT el_type;

        chunk_read_value_entry(&chunk, j);

        /* we need the element parameters here */
        chunk_read_size_entry(&(field->ele_param), j);

        Id      = field->ele_param.size_buf[element_variables.ep_size_Id];
        el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];

        /* translate to internal */
        el_type = field->problem->element_type.table[el_type];

        switch (el_type)
        {
#ifdef D_WALL1
        case el_wall1:
        {
          INT i;

          INT nGP[4];           /* why that large? */
          INT stresstype;
          INT ngauss;

          nGP[0] = field->ele_param.size_buf[wall1_variables.ep_size_nGP0];
          nGP[1] = field->ele_param.size_buf[wall1_variables.ep_size_nGP1];
          ngauss = nGP[0]*nGP[1];

          dsassert(ngauss == 4, "only four gauss point versions supported currently");

          /*
           * get the element's stress type and translate in to
           * internal */
          stresstype = field->ele_param.size_buf[wall1_variables.ep_size_stresstyp];
          stresstype = w1_stresstype.table[stresstype];

          fprintf(out,DBLLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n", Id, j);
          fprintf(out,"\n");

          switch (stresstype)
          {
          case w1_xy:
            fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
            break;
          case w1_rs:
            fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
            break;
          default:
            dserror("stress type %d unknown", stresstype);
          }

          for (i=0; i<ngauss; i++) {
            /*
             * This is a special case. Four stress values, three other
             * values to be output here and two values to be output
             * using GiD. The wall element needs to be cleaned up! */
            fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                    i,
                    chunk.value_buf[0*9+i],
                    chunk.value_buf[1*9+i],
                    chunk.value_buf[2*9+i],
                    chunk.value_buf[3*9+i],
                    chunk.value_buf[4*9+i],
                    chunk.value_buf[5*9+i],
                    chunk.value_buf[6*9+i]
              );
          }
          break;
        }
#endif /*D_WALL1*/

#ifdef D_BRICK1
        case el_brick1:
        {
          INT i;

          INT stresstype;
          INT ngauss;

          /*
           * get the element's stress type and translate in to
           * internal */
          stresstype = field->ele_param.size_buf[brick1_variables.ep_size_stresstyp];
          stresstype = c1_stresstype.table[stresstype];

          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                BRICK1\n", Id, j);
          fprintf(out,"\n");

          switch (stresstype)
          {
          case c1_gpxyz:
            fprintf(out,"INT.point   x-coord.     y-coord.     z-coord.     stress-xx    stress-yy    stress-zz    stress-xy    stress-xz    stress-yz\n");
            for (i=0; i<ngauss; i++)
            {
              /* There are 27 values to one gauss point */
              fprintf(out,"  %-6d %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                      i,
                      chunk.value_buf[24+27*i],
                      chunk.value_buf[25+27*i],
                      chunk.value_buf[26+27*i],
                      chunk.value_buf[ 6+27*i],
                      chunk.value_buf[ 7+27*i],
                      chunk.value_buf[ 8+27*i],
                      chunk.value_buf[ 9+27*i],
                      chunk.value_buf[10+27*i],
                      chunk.value_buf[11+27*i]
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
                      chunk.value_buf[24+27*i],
                      chunk.value_buf[25+27*i],
                      chunk.value_buf[26+27*i],
                      chunk.value_buf[ 0+27*i],
                      chunk.value_buf[ 1+27*i],
                      chunk.value_buf[ 2+27*i],
                      chunk.value_buf[ 3+27*i],
                      chunk.value_buf[ 4+27*i],
                      chunk.value_buf[ 5+27*i]
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
                      chunk.value_buf[12+27*i],
                      chunk.value_buf[13+27*i],
                      chunk.value_buf[14+27*i],
                      chunk.value_buf[15+27*i],
                      chunk.value_buf[16+27*i],
                      chunk.value_buf[17+27*i],
                      chunk.value_buf[18+27*i],
                      chunk.value_buf[19+27*i],
                      chunk.value_buf[20+27*i],
                      chunk.value_buf[21+27*i],
                      chunk.value_buf[22+27*i],
                      chunk.value_buf[23+27*i]
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
        case el_beam3:
        {
          INT i;
          INT ngauss;

          ngauss = field->ele_param.size_buf[beam3_variables.ep_size_nGP0];

          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                BEAM3\n", Id, j);
          fprintf(out,"\n");
          fprintf(out,"Gaussian         Nx           Vy           Vz           Mx           My           Mz\n");
          for (i=0; i<ngauss; i++) {
            /* There are 6 values to one gauss point */
            fprintf(out,"Gauss %d       %12.3E %12.3E %12.3E %12.3E %12.3E %12.3E \n",
                    i,
                    chunk.value_buf[0+6*i],
                    chunk.value_buf[1+6*i],
                    chunk.value_buf[2+6*i],
                    chunk.value_buf[3+6*i],
                    chunk.value_buf[4+6*i],
                    chunk.value_buf[5+6*i]
              );
          }
          break;
        }
#endif /*D_BEAM3*/

#ifdef D_INTERF
        case el_interf:
        {
          INT i;
          INT ngauss;
          INT stresstype;

          ngauss = field->ele_param.size_buf[interf_variables.ep_size_nGP];

          /*
           * get the element's stress type and translate in to
           * internal */
          stresstype = field->ele_param.size_buf[interf_variables.ep_size_stresstyp];
          stresstype = interf_stresstype.table[stresstype];

          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                INTERF\n", Id, j);

          switch (stresstype)
          {
          case if_xy:
            fprintf(out,"Gaussian     stresses-xx     stresses-yy    stresses-xy\n");
            break;
          case if_tn:
            fprintf(out,"Gaussian     stresses-tangential     stresses-normal\n");
            break;
          default:
            dserror("stress type %d unknown", stresstype);
          }

          for (i=0; i<ngauss; i++) {
            /* There are 5 values to one gauss point */
            fprintf(out,"Gauss %d   %12.3E %12.3E %12.3E\n",
                    i,
                    chunk.value_buf[0+5*i],
                    chunk.value_buf[1+5*i],
                    chunk.value_buf[2+5*i] );
          }
          break;
        }
#endif /*D_INTERF*/

#ifdef D_WALLGE
        case el_wallge:
        {
          INT i;
          INT nGP[2];
          INT ngauss;
          INT stresstype;

          nGP[0] = field->ele_param.size_buf[wallge_variables.ep_size_nGP0];
          nGP[1] = field->ele_param.size_buf[wallge_variables.ep_size_nGP1];
          ngauss = nGP[0]*nGP[1];

          fprintf(out,UNDERLINE);
          fprintf(out,"Element glob_Id %d loc_Id %d                WALL1\n", Id, j);
          fprintf(out,"\n");

          /*
           * get the element's stress type and translate in to
           * internal */
          stresstype = field->ele_param.size_buf[wallge_variables.ep_size_stresstyp];
          stresstype = wallge_stresstype.table[stresstype];

          switch (stresstype)
          {
          case wge_xy:
            fprintf(out,"Gaussian     Stress-xx    Stress-yy    Stress-xy    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
            break;
          case wge_rs:
            fprintf(out,"Gaussian     Stress-rr    Stress-ss    Stress-rs    Stress-zz    Max. P.S.    Min. P.S.    Angle\n");
            break;
          default:
            dserror("stress type %d unknown", stresstype);
          }

          for (i=0; i<ngauss; i++)
          {
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

      destroy_chunk_data(&chunk);

#ifdef D_WALLGE
      destroy_translation_table(&wallge_stresstype);
#endif
#ifdef D_INTERF
      destroy_translation_table(&interf_stresstype);
#endif
#ifdef D_BRICK1
      destroy_translation_table(&c1_stresstype);
#endif
#ifdef D_WALL1
      destroy_translation_table(&w1_stresstype);
#endif
    }
  }

end:

  fprintf(out,"\n");
  fprintf(out,"\n");

#ifdef DEBUG
  dstrc_exit();
#endif
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

#ifdef DEBUG
  dstrc_enter("post_out_fsi");
#endif

  fprintf(out,DBLLINE);
  fprintf(out,"FSI node connectivity global Ids:\n");
  fprintf(out,DBLLINE);
  fprintf(out,"\n");
  fprintf(out,"FLUID          ALE          STRUCTURE\n");

  /* Find the corresponding discretizations. We don't rely on any order. */
  for (i=0; i<problem->num_discr; ++i)
  {
    if (problem->discr[i].type == structure)
    {
      struct_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == fluid)
    {
      fluid_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == ale)
    {
      ale_field = &(problem->discr[i]);
    }
    else
    {
      dserror("unknown field type %d", problem->discr[i].type);
    }
  }

  post_find_fsi_coupling(problem,
                         struct_field, fluid_field, ale_field,
                         &fluid_struct_connect, &fluid_ale_connect);

  for (i=0; i<fluid_field->numnp; i++)
  {
    chunk_read_size_entry(&(fluid_field->ele_param), i);
    if ((fluid_ale_connect[i]!=-1) && (fluid_struct_connect[i]==-1))
    {
      chunk_read_size_entry(&(ale_field->ele_param),    fluid_ale_connect[i]);
      fprintf(out,"%-6d         %-6d       ------\n",
              fluid_field->ele_param.size_buf[element_variables.ep_size_Id],
              ale_field->ele_param.size_buf[element_variables.ep_size_Id]);
    }
    else if ((fluid_ale_connect[i]==-1) && (fluid_struct_connect[i]!=-1))
    {
      chunk_read_size_entry(&(struct_field->ele_param), fluid_struct_connect[i]);
      fprintf(out,"%-6d         ------       %-6d\n",
              fluid_field->ele_param.size_buf[element_variables.ep_size_Id],
              struct_field->ele_param.size_buf[element_variables.ep_size_Id]);
    }
    else if ((fluid_ale_connect[i]!=-1) && (fluid_struct_connect[i]!=-1))
    {
      chunk_read_size_entry(&(ale_field->ele_param),    fluid_ale_connect[i]);
      chunk_read_size_entry(&(struct_field->ele_param), fluid_struct_connect[i]);
      fprintf(out,"%-6d         %-6d       %-6d\n",
              fluid_field->ele_param.size_buf[element_variables.ep_size_Id],
              ale_field->ele_param.size_buf[element_variables.ep_size_Id],
              struct_field->ele_param.size_buf[element_variables.ep_size_Id]);
    }
    else
      fprintf(out,"%-6d         ------       ------\n",
              fluid_field->ele_param.size_buf[element_variables.ep_size_Id]);

  }
  fprintf(out,UNDERLINE);

#if 0
  /* Not worth it (?) */
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

#ifdef DEBUG
  dstrc_exit();
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
  CHAR filename[100];
  PROBLEM_DATA problem;
  FILE* f;
  INT i;

  init_problem_data(&problem, argc, argv);

  sprintf(filename, "%s.out", problem.basename);
  f = fopen(filename, "wb");

  post_out_general(f, &problem);

#ifdef D_FSI
  /* Output the coupling information if it's a fsi problem. It's a
   * little harder to discover a fluid-ale problem, but the function
   * is supposed to work the same way. Maybe we could find a new
   * problemtype for fluid-ale? */
  if (strcmp(map_read_string(&(problem.control_table), "problem_type"), "fsi")==0) {
    post_out_fsi(f, &problem);
  }
#endif

  /* Iterate all discretizations. */
  for (i=0; i<problem.num_discr; ++i) {
    FIELD_DATA* field;
    RESULT_DATA result;

    field = &(problem.discr[i]);

    /* Iterate all results. */
    init_result_data(field, &result);

    while (next_result(&result))
    {
      post_out_sol(f, &result);
    }

    destroy_result_data(&result);
  }

  fclose(f);
  post_log(4, "Done.\n");
  return 0;
}
