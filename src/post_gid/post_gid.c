/*!
\file
\brief Postprocessing utility that takes ccarat output and produces a
GiD input file.

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
\date 09/04

*/

#include <assert.h>

#include "gid_out.h"

#include "post_gid.h"

#include "post_ale2.h"
#include "post_ale3.h"
#include "post_axishell.h"
#include "post_beam3.h"
#include "post_brick1.h"
#include "post_fluid2.h"
#include "post_fluid2_pro.h"
#include "post_fluid3.h"
#include "post_fluid3_fast.h"
#include "post_interf.h"
#include "post_shell8.h"
#include "post_shell9.h"
#include "post_wall1.h"
#include "post_wallge.h"

#include "../output/gid.h"



/*----------------------------------------------------------------------*/
/*!
  \brief Convert an array of field local node ids to global node ids
  in GiD (fortran) style.

  This is an utility function called when meshes are written. There we
  read the mesh chunk that contains the field local node ids. But we
  have to give global node ids to GiD. Furthermore the internal ids
  are counted from zero but GiD requires ids to be counted from one.

  \param field      (i) the field we write
  \param local_ids  (i) array of known field local ids
  \param global_ids (o) array of gid ids to be found
  \param numnp      (i) number of node ids in these arrays

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void get_gid_node_ids(FIELD_DATA* field, INT* local_ids, INT* global_ids, INT numnp)
{
  INT k;

#ifdef DEBUG
  dstrc_enter("get_gid_node_ids");
#endif

  for (k=0; k<numnp; ++k)
  {
    INT node_id;
    node_id = local_ids[k];

    chunk_read_size_entry(&(field->coords), node_id);
    global_ids[k] = field->coords.size_buf[node_variables.coords_size_Id] + 1;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write coordinates.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void write_coords(FIELD_DATA* field, GIDSET* gid)
{
  INT i;
  INT j;
  DOUBLE x[3];

#ifdef DEBUG
  dstrc_enter("write_coords");
#endif

  GiD_BeginCoordinates();

  /* This is a shortcut. We'll always write three coordinates. That's
   * what the gid library will do anyway. So in case we are reading a
   * two dimensional problem just set the last component to zero. */
  x[2] = 0;

  for (i=0; i<field->numnp; ++i)
  {
    INT Id;

    chunk_read_value_entry(&(field->coords), i);
    for (j=0; j<field->problem->ndim; ++j)
    {
      x[j] = field->coords.value_buf[j];
    }

    chunk_read_size_entry(&(field->coords), i);
    Id = field->coords.size_buf[node_variables.coords_size_Id];

    /* All GiD ids are one based. */
    GiD_WriteCoordinates(Id+1, x[0], x[1], x[2]);
  }

  GiD_EndCoordinates();

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the element domains.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_domain(FIELD_DATA *field, GIDSET* gid)
{
  CHUNK_DATA chunk;

#ifdef DEBUG
  dstrc_enter("write_domain");
#endif

  post_log(3, "%s: Write domain\n", field->name);
  init_chunk_data(&(field->head), &chunk, "domain");

#ifdef D_SHELL8
  shell8_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_SHELL9
  shell9_write_domain(field, gid, &chunk);
#endif /*D_SHELL9*/

  /*--------------------------------------------------------------------*/

#ifdef D_BRICK1
  brick1_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_FLUID2
  fluid2_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_FLUID3
  fluid3_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_FLUID3_F
  fluid3_fast_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_ALE
  ale2_write_domain(field, gid, &chunk);
  ale3_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_WALL1
  wall1_write_domain(field, gid, &chunk);
#endif

  /*--------------------------------------------------------------------*/

#ifdef D_BEAM3
  beam3_write_domain(field, gid, &chunk);
#endif

  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal displacements for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_displacement(FIELD_DATA *field, RESULT_DATA* result)
{
#ifdef DEBUG
  dstrc_enter("write_displacement");
#endif

#ifdef D_SHELL8
  if (field->is_shell8_problem)
  {
    shell8_write_displacement(field, result);
  }
  else
#endif
#ifdef D_SHELL9
  if (field->is_shell9_problem)
  {
    shell9_write_displacement(field, result);
  }
  else
#endif
  {
    INT k;
    DOUBLE x[3];
    CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };
    DOUBLE time;
    INT step;
    CHAR buf[100];
    CHUNK_DATA chunk;

    init_chunk_data(result, &chunk, "displacement");

    time = map_read_real(result->group, "time");
    step = map_read_int(result->group, "step");

    post_log(3, "%s: Write displacement of step %d\n", field->name, step);

    sprintf(buf, "%s_displacement", fieldnames[field->type]);
    GiD_BeginResult(buf, "ccarat", step, GiD_Vector, GiD_OnNodes,
                    NULL, NULL, chunk.value_entry_length, componentnames);

    /* In case this is a 2d problem. */
    x[2] = 0;

    /* Iterate all nodes. */
    for (k = 0; k < field->numnp; ++k)
    {
      INT i;
      chunk_read_value_entry(&chunk, k);
      for (i=0; i<chunk.value_entry_length; ++i)
      {
        x[i] = chunk.value_buf[i];
      }

      chunk_read_size_entry(&(field->coords), k);

      GiD_WriteVector(field->coords.size_buf[node_variables.coords_size_Id]+1,
                      x[0], x[1], x[2]);
    }

    GiD_EndResult();
    destroy_chunk_data(&chunk);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal velocities for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_velocity(FIELD_DATA *field, RESULT_DATA* result)
{
  INT k;
  DOUBLE velocity[3];
  CHAR* componentnames[] = { "x-vel", "y-vel", "z-vel" };
  DOUBLE time;
  INT step;
  CHAR buf[100];
  CHUNK_DATA chunk;

#ifdef DEBUG
  dstrc_enter("write_velocity");
#endif

  init_chunk_data(result, &chunk, "velocity");

  time = map_read_real(result->group, "time");
  step = map_read_int(result->group, "step");

  post_log(3, "%s: Write velocity of step %d\n", field->name, step);

  sprintf(buf, "%s_velocity", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Vector, GiD_OnNodes,
                  NULL, NULL, chunk.value_entry_length, componentnames);

  /* In case this is a 2d problem. */
  velocity[2] = 0;

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k)
  {
    INT i;
    chunk_read_value_entry(&chunk, k);
    for (i=0; i<chunk.value_entry_length; ++i)
    {
      velocity[i] = chunk.value_buf[i];
    }

    chunk_read_size_entry(&(field->coords), k);

    GiD_WriteVector(field->coords.size_buf[node_variables.coords_size_Id]+1,
                    velocity[0], velocity[1], velocity[2]);
  }

  GiD_EndResult();
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal pressure for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_pressure(FIELD_DATA *field, RESULT_DATA* result)
{
  INT k;
  CHAR* componentnames[] = { "pressure" };
  DOUBLE time;
  INT step;
  CHAR buf[100];
  CHUNK_DATA chunk;

#ifdef DEBUG
  dstrc_enter("write_pressure");
#endif

  init_chunk_data(result, &chunk, "pressure");

  time = map_read_real(result->group, "time");
  step = map_read_int(result->group, "step");

  post_log(3, "%s: Write pressure of step %d\n", field->name, step);

  sprintf(buf, "%s_pressure", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Scalar, GiD_OnNodes,
                  NULL, NULL, chunk.value_entry_length, componentnames);

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k)
  {
    chunk_read_value_entry(&chunk, k);
    chunk_read_size_entry(&(field->coords), k);
    GiD_WriteScalar(field->coords.size_buf[node_variables.coords_size_Id]+1,
                    chunk.value_buf[0]);
  }

  GiD_EndResult();
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the gauss point stresses for one time step.

  This is element specific again.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_stress(FIELD_DATA *field, GIDSET* gid, RESULT_DATA* result)
{
  DOUBLE time;
  INT step;
  CHUNK_DATA chunk;

#ifdef DEBUG
  dstrc_enter("write_stress");
#endif

  init_chunk_data(result, &chunk, "stress");

  time = map_read_real(result->group, "time");
  step = map_read_int(result->group, "step");

  post_log(3, "%s: Write stress of step %d\n", field->name, step);

#ifdef D_SHELL8
  shell8_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_SHELL9
  shell9_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_WALL1
  wall1_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_BRICK1
  brick1_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_FLUID3
  fluid3_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_FLUID3_F
  fluid3_fast_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_AXISHELL
  axishell_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_BEAM3
  beam3_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_INTERF
  interf_write_stress(field, gid, &chunk, time, step);
#endif
#ifdef D_WALLGE
  wallge_write_stress(field, gid, &chunk, time, step);
#endif

  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find out what types of elements there are in this discretization.

  This is the most primitive way of handling element variants. So it's
  probably the best.

  This function is based on \a out_gid_sol_init of the original ccarat
  output.

  \author u.kue
  \date 11/04
  \sa out_gid_sol_init
*/
/*----------------------------------------------------------------------*/
static void setup_gid_flags(FIELD_DATA* field, GIDSET* gid)
{
  INT i;
#ifdef DEBUG
  dstrc_enter("setup_gid_flags");
#endif

  memset(gid, 0, sizeof(GIDSET));
  for (i=0; i<field->numele; ++i)
  {
    INT numnp;
    INT dis;
    INT el_type;
    INT Id;
    INT* size_ptr;

    /* read the element's data */
    get_element_params(field, i, &Id, &el_type, &dis, &numnp);

    size_ptr = field->ele_param.size_buf;

    switch (el_type)
    {
#ifdef D_SHELL8
    case el_shell8:
    {
      if (numnp==4)
      {
        gid->is_shell8_4_22   = 1;
        gid->shell8_4_22_name = "shell8_4_22";
      }
      if (numnp==8 || numnp==9)
      {
        gid->is_shell8_8_33   = 1;
        gid->shell8_8_33_name = "shell8_8_33";
      }
      break;
    }
#endif
#ifdef D_SHELL9
    case el_shell9:
    {
      INT nGP[2];

      nGP[0] = size_ptr[shell9_variables.ep_size_nGP0];
      nGP[1] = size_ptr[shell9_variables.ep_size_nGP1];

      if (numnp==4)
      {
        if (nGP[0]==2 && nGP[1]==2)
        {
          gid->is_shell9_4_22   = 1;
          gid->shell9_4_22_name = "shell9_4_22";
        }
        if (nGP[0]==3 && nGP[1]==3)
        {
          gid->is_shell9_4_33   = 1;
          gid->shell9_4_33_name = "shell9_4_33";
        }
      }
      if (numnp==8)
      {
        if (nGP[0]==2 && nGP[1]==2)
        {
          gid->is_shell9_8_22   = 1;
          gid->shell9_8_22_name = "shell9_8_22";
        }
        if (nGP[0]==3 && nGP[1]==3)
        {
          gid->is_shell9_8_33   = 1;
          gid->shell9_8_33_name = "shell9_8_33";
        }
      }
      if (numnp==9)
      {
        if (nGP[0]==2 && nGP[1]==2)
        {
          gid->is_shell9_9_22   = 1;
          gid->shell9_9_22_name = "shell9_9_22";
        }
        if (nGP[0]==3 && nGP[1]==3)
        {
          gid->is_shell9_9_33   = 1;
          gid->shell9_9_33_name = "shell9_9_33";
        }
      }
      break;
    }
#endif
#ifdef D_BRICK1
    case el_brick1:
    {
      if (numnp==8)
      {
        gid->is_brick1_222   = 1;
        gid->brick1_222_name = "brick1_222";
      }
      if (numnp==20 || numnp==27)
      {
        gid->is_brick1_333   = 1;
        gid->brick1_333_name = "brick1_333";
      }
      break;
    }
#endif
#ifdef D_FLUID2
    case el_fluid2:
    {
      if (numnp==4)
      {
        gid->is_fluid2_22    = 1;
        gid->fluid2_22_name  = "fluid2_22";
      }
      if (numnp==8  || numnp==9)
      {
        gid->is_fluid2_33    = 1;
        gid->fluid2_33_name  = "fluid2_33";
      }

      /*
       * There are many ways to implement triangle elements. But we
       * are only concerned about the number of nodes. */
      if (numnp==3)
      {
        gid->is_fluid2_tri3    = 1;
        gid->fluid2_tri3_name  = "fluid2_3";
      }
      if (numnp==6)
      {
        gid->is_fluid2_6    = 1;
        gid->fluid2_6_name  = "fluid2_6";
      }
      break;
    }
#endif
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      if (numnp==4)
      {
        gid->is_fluid2_pro_22    = 1;
        gid->fluid2_pro_22_name  = "fluid2_pro_22";
      }
      if (numnp==8  || numnp==9)
      {
        gid->is_fluid2_pro_33    = 1;
        gid->fluid2_pro_33_name  = "fluid2_pro_33";
      }
      break;
    }
#endif
#ifdef D_FLUID2TU
    case el_fluid2_tu:
    {
      dserror("no output of el_fluid2_tu yet");
      break;
    }
#endif
#ifdef D_FLUID3
    case el_fluid3:
    {
      if (numnp==8)
      {
        gid->is_fluid3_222   = 1;
        gid->fluid3_222_name = "fluid3_222";
      }
      if (numnp==20 || numnp==27)
      {
        gid->is_fluid3_333   = 1;
        gid->fluid3_333_name = "fluid3_333";
      }
      break;
    }
#endif
#ifdef D_FLUID3_F
    case el_fluid3_fast:
    {
      if (numnp==8)
      {
        gid->is_f3f_8_222   = 1;
        gid->f3f_8_222_name = "f3f_8_222";
      }
      if (numnp==20 || numnp==27)
      {
        gid->is_f3f_20_333   = 1;
        gid->f3f_20_333_name = "f3f_333";
      }
      break;
    }
#endif
#ifdef D_ALE
    case el_ale2:
    {
      INT nGP[2];

      nGP[0] = size_ptr[ale2_variables.ep_size_nGP0];
      nGP[1] = size_ptr[ale2_variables.ep_size_nGP1];

      if (numnp==4)
      {
        if (nGP[0]==1 && nGP[1]==1)
        {
          gid->is_ale_11    = 1;
          gid->ale_11_name  = "ale_11";
        }
        else if (nGP[0]==2 && nGP[1]==2 )
        {
          gid->is_ale_22    = 1;
          gid->ale_22_name  = "ale_22";
        }
      }
      if (numnp==3)
      {
        if (nGP[0] == 1)
        {
          gid->is_ale_tri_1    = 1;
          gid->ale_tri_1_name  = "ale_tri_1";
        }
        else if (nGP[0] == 3)
        {
          gid->is_ale_tri_3    = 1;
          gid->ale_tri_3_name  = "ale_tri_3";
        }
      }
      break;
    }
#endif
#ifdef D_ALE
    case el_ale3:
    {
      INT nGP[3];

      nGP[0] = size_ptr[ale3_variables.ep_size_nGP0];
      nGP[1] = size_ptr[ale3_variables.ep_size_nGP1];
      nGP[2] = size_ptr[ale3_variables.ep_size_nGP2];

      if (numnp==8)
      {
        if (nGP[0]==1 && nGP[1]==1 && nGP[2]==1 )
        {
          gid->is_ale_8_111   = 1;
          gid->ale_8_111_name = "ale_8_111";
        }
        else if (nGP[0]==2 && nGP[1]==2 && nGP[2]==2 )
        {
          gid->is_ale_8_222   = 1;
          gid->ale_8_222_name = "ale_222";
        }
      }
      if (numnp==4)
      {
        if (nGP[0] == 1)
        {
          gid->is_ale_tet_1    = 1;
          gid->ale_tet_1_name  = "ale_tet_1";
        }
        else if (nGP[0] == 4)
        {
          gid->is_ale_tet_4    = 1;
          gid->ale_tet_4_name  = "ale_tet_4";
        }
      }
      break;
    }
#endif
#ifdef D_WALL1
    case el_wall1:
    {
      INT nGP[2];

      nGP[0] = size_ptr[wall1_variables.ep_size_nGP0];
      nGP[1] = size_ptr[wall1_variables.ep_size_nGP1];

      if (nGP[0]==1)
      {
        gid->is_wall1_11    = 1;
        gid->wall1_11_name  = "wall1_11";
      }
      if (nGP[0]==2)
      {
        gid->is_wall1_22    = 1;
        gid->wall1_22_name  = "wall1_22";
      }
      if (nGP[0]==3)
      {
        gid->is_wall1_33   = 1;
        gid->wall1_33_name = "wall1_33";
      }
      break;
    }
#endif
#ifdef D_BEAM3
    case el_beam3:
    {
      INT nGP[1];

      nGP[0] = size_ptr[beam3_variables.ep_size_nGP0];

      if (numnp==2)
      {
        if (nGP[0]==1)
        {
          gid->is_beam3_21 = 1;
          gid->beam3_21_name  = "beam3_21";
        }
        if (nGP[0]==2)
        {
          gid->is_beam3_22 = 1;
          gid->beam3_22_name  = "beam3_22";
        }
      }
      if (numnp==3)
      {
        if (nGP[0]==2)
        {
          gid->is_beam3_32 = 1;
          gid->beam3_32_name  = "beam3_32";
        }
        if (nGP[0]==3)
        {
          gid->is_beam3_33 = 1;
          gid->beam3_33_name  = "beam3_33";
        }
      }
      break;
    }
#endif
#ifdef D_AXISHELL
    case el_axishell:
    {
      gid->is_axishell    = 1;
      gid->axishell_name  = "axishell";
      break;
    }
#endif
#ifdef D_INTERF
    case el_interf:
    {
      INT nGP;

      nGP = size_ptr[interf_variables.ep_size_nGP];

      if (nGP==2)
      {
        gid->is_interf_22    = 1;
        gid->interf_22_name  = "interf_22";
      }
      if (nGP==3)
      {
        gid->is_interf_33    = 1;
        gid->interf_33_name  = "interf_33";
      }
      break;
    }
#endif
#ifdef D_WALLGE
    case el_wallge:
    {
      INT nGP[2];

      nGP[0] = size_ptr[wallge_variables.ep_size_nGP0];
      nGP[1] = size_ptr[wallge_variables.ep_size_nGP1];

      if (nGP[0]==2)
      {
        gid->is_wallge_22    = 1;
        gid->wallge_22_name  = "wallge_22";
      }
      if (nGP[0]==3)
      {
        gid->is_wallge_33    = 1;
        gid->wallge_33_name  = "wallge_33";
      }
      break;
    }
#endif
    default:
      dserror_args(__FILE__, __LINE__, "element type %d unknown", el_type);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write connectivity and node coordinates.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_mesh(FIELD_DATA* field, GIDSET* gid)
{
  INT first_mesh = 1;

#ifdef DEBUG
  dstrc_enter("write_mesh");
#endif

  /* Ok, we have the most primitive description possible. No way to
   * loop this data. We have to do each case by hand. (Not such a big
   * improvement after all.) */

#ifdef D_SHELL8
  shell8_write_gauss(gid);
#endif
#ifdef D_SHELL9
  shell9_write_gauss(gid);
#endif
#ifdef D_BRICK1
  brick1_write_gauss(gid);
#endif
#ifdef D_FLUID2
  fluid2_write_gauss(gid);
#endif
#ifdef D_FLUID2_PRO
  fluid2_pro_write_gauss(gid);
#endif
#ifdef D_FLUID3
  fluid3_write_gauss(gid);
#endif
#ifdef D_FLUID3_F
  fluid3_fast_write_gauss(gid);
#endif
#ifdef D_ALE
  ale2_write_gauss(gid);
  ale3_write_gauss(gid);
#endif
#ifdef D_WALL1
  wall1_write_gauss(gid);
#endif
#ifdef D_BEAM3
  beam3_write_gauss(gid);
#endif
#ifdef D_AXISHELL
  axishell_write_gauss(gid);
#endif
#ifdef D_INTERF
  interf_write_gauss(gid);
#endif
#ifdef D_WALLGE
  wallge_write_gauss(gid);
#endif

  /*--------------------------------------------------------------------*/
  /* Write a mesh for each element variant */

  if (field->mesh.size_entry_length > MAXNOD)
  {
    dserror("MAXNOD overflow");
  }

#ifdef D_SHELL8
  shell8_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_SHELL9
  shell9_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_BRICK1
  brick1_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_FLUID2
  fluid2_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_FLUID2_PRO
  fluid2_pro_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_FLUID3
  fluid3_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_FLUID3_F
  fluid3_fast_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_ALE
  ale2_write_mesh(field, gid, &first_mesh);
  ale3_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_WALL1
  wall1_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_BEAM3
  beam3_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_AXISHELL
  axishell_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_INTERF
  interf_write_mesh(field, gid, &first_mesh);
#endif
#ifdef D_WALLGE
  wallge_write_mesh(field, gid, &first_mesh);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a whole field.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_field(PROBLEM_DATA* problem, INT num)
{
  FIELD_DATA* field;
  GIDSET gid;
  RESULT_DATA result;

#ifdef DEBUG
  dstrc_enter("write_field");
#endif

  field = &(problem->discr[num]);

  /* find the element variants that are used in this discretization */
  setup_gid_flags(field, &gid);

  /*--------------------------------------------------------------------*/
  /* Write connectivity and node coordinates. */

  write_mesh(field, &gid);

  /*--------------------------------------------------------------------*/
  /* Now we could define a few range tables. Maybe it's nice to have
   * per discretization range tables. But they are highly specific. */
#if 0
  GiD_BeginRangeTable( char * name );
  GiD_WriteMinRange( double max, char * name );
  GiD_WriteRange( double min, double max, char * name );
  GiD_WriteMaxRange( double min, char * name );
  GiD_EndRangeTable();
#endif

  /*--------------------------------------------------------------------*/
  /* Optionally there is one domain table per field. It's the first result. */

#if 1
  if (map_has_map(field->group, "domain"))
  {
    write_domain(field, &gid);
  }
#endif

  /*--------------------------------------------------------------------*/
  /* Now it's time to write the time dependend results. */

  /* Iterate all results. */
  init_result_data(field, &result);

  while (next_result(&result))
  {
    /* nodal result are written independent of the elements in this
     * discretization */
    if (map_has_map(result.group, "displacement"))
    {
      write_displacement(field, &result);
    }
    if (map_has_map(result.group, "velocity"))
    {
      write_velocity(field, &result);
    }
    if (map_has_map(result.group, "pressure"))
    {
      write_pressure(field, &result);
    }

    /* element dependend results */
    if (map_has_map(result.group, "stress"))
    {
      write_stress(field, &gid, &result);
    }
  }
  destroy_result_data(&result);

#ifdef DEBUG
  dstrc_exit();
#endif
}


#ifdef D_FSI

/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing all two fsi fields.

  The ale field is not written on it's own. Instead the ale
  displacements are written to the fluid nodes. The structure field is
  handled by standard methods, but only if it's there. There are
  fluid-ale problems, too.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void write_fsi(PROBLEM_DATA* problem)
{
  INT i;
  FIELD_DATA* struct_field = NULL;
  FIELD_DATA* fluid_field;
  FIELD_DATA* ale_field;

  INT* fluid_ale_connect;

  GIDSET struct_gid;
  GIDSET fluid_gid;

  RESULT_DATA result;
  RESULT_DATA struct_result;

#ifdef DEBUG
  dstrc_enter("write_fsi");
#endif

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
      dserror_args(__FILE__, __LINE__,
          "unknown field type %d", problem->discr[i].type);
    }
  }

  /*--------------------------------------------------------------------*/
  /* Write connectivity and node coordinates. */

  /* find the element variants that are used in the fluid discretization */
  setup_gid_flags(fluid_field, &fluid_gid);

  write_mesh(fluid_field, &fluid_gid);
  if (struct_field != NULL)
  {
    /* find the element variants that are used in the struct discretization */
    setup_gid_flags(struct_field, &struct_gid);

    write_mesh(struct_field, &struct_gid);
  }

  /*--------------------------------------------------------------------*/
  /* Optionally there is one domain table per field. It's the first result. */

#if 1
  if (map_has_map(fluid_field->group, "domain"))
  {
    write_domain(fluid_field, &fluid_gid);
  }
  if (struct_field != NULL)
  {
    if (map_has_map(struct_field->group, "domain"))
    {
      write_domain(struct_field, &struct_gid);
    }
  }
#endif

  /*--------------------------------------------------------------------*/
  /* We need to find the connection between ale and fluid nodes. */

  {
    INT *fluid_struct_connect;
    post_find_fsi_coupling(problem,
                           struct_field, fluid_field, ale_field,
                           &fluid_struct_connect, &fluid_ale_connect);

    /*
     * We don't need the connection to the structure field here. Free
     * it immediately. */
    if (struct_field != NULL)
    {
      CCAFREE(fluid_struct_connect);
    }
  }

  /*--------------------------------------------------------------------*/
  /* Now it's time to write the time dependend results. */

  /* Do the fluid results. */
  init_result_data(fluid_field, &result);
  while (next_result(&result))
  {
    if (map_has_map(result.group, "velocity"))
    {
      write_velocity(fluid_field, &result);
    }
    if (map_has_map(result.group, "pressure"))
    {
      write_pressure(fluid_field, &result);
    }

#if 0
    /* No fluid stresses yet? */

    /* element dependend results */
    if (map_has_map(result_group, "stress")) {
      write_stress(fluid_field, result_group);
    }
#endif
  }
  destroy_result_data(&result);


  /* Do the structure results. */
  if (struct_field != NULL)
  {
    init_result_data(struct_field, &result);
    while (next_result(&result))
    {
      /*
      if (map_has_map(result.group, "displacement"))
      {
        write_displacement(struct_field, &result);
      }
      */
      if (map_has_map(result.group, "velocity"))
      {
        write_velocity(struct_field, &result);
      }

      if (map_has_map(result.group, "stress"))
      {
        write_stress(struct_field, &struct_gid, &result);
      }
    }
    destroy_result_data(&result);
  }

  /* Do the ale results. */
  init_result_data(ale_field, &result);

  /* We expect struct and ale results for the same time steps! No
   * discrepancies allowed! */
  if (struct_field != NULL)
    init_result_data(struct_field, &struct_result);

  while (next_result(&result) &&
         ((struct_field == NULL) || next_result(&struct_result)))
  {
    /* We read the displacements of the ale elements and write fluid
     * displacments. That's why this is special and not treated by the
     * normal function. */
    if (map_has_map(result.group, "displacement"))
    {
      INT k;
      DOUBLE x[3];
      CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };
      DOUBLE time;
      INT step;
      CHUNK_DATA chunk;
      CHUNK_DATA struct_chunk;

      init_chunk_data(&result, &chunk, "displacement");

      time = map_read_real(result.group, "time");
      step = map_read_int(result.group, "step");

      if (struct_field != NULL)
      {
        if (!map_has_map(struct_result.group, "displacement"))
        {
          dserror("structure displacement required for fsi");
        }

        init_chunk_data(&struct_result, &struct_chunk, "displacement");
        if ((fabs(time - map_read_real(struct_result.group, "time")) > 1e-6) ||
            (step != map_read_int(struct_result.group, "step")))
        {
          dserror("struct -- ale result discrepancy");
        }
      }

      post_log(3, "%s: Write displacement of step %d\n", ale_field->name, step);

      GiD_BeginResult("displacement", "ccarat", step, GiD_Vector, GiD_OnNodes,
                      NULL, NULL, chunk.value_entry_length, componentnames);

      /* In case this is a 2d problem. */
      x[2] = 0;

      /* Iterate all fluid nodes. */
      for (k = 0; k < fluid_field->numnp; ++k)
      {
        chunk_read_size_entry(&(fluid_field->coords), k);
        if (fluid_ale_connect[k] != -1)
        {
          INT i;
          chunk_read_value_entry(&chunk, fluid_ale_connect[k]);
          for (i=0; i<chunk.value_entry_length; ++i)
          {
            x[i] = chunk.value_buf[i];
          }

          GiD_WriteVector(fluid_field->coords.size_buf[node_variables.coords_size_Id]+1,
                          x[0], x[1], x[2]);
        }
        else
        {
          GiD_WriteVector(fluid_field->coords.size_buf[node_variables.coords_size_Id]+1,
                          0, 0, 0);
        }
      }

      /*
       * All displacements are to be shown together. So they must be
       * in one result set. */
      if (struct_field != NULL)
      {
        for (k = 0; k < struct_field->numnp; ++k)
        {
          INT i;
          chunk_read_size_entry(&(struct_field->coords), k);
          chunk_read_value_entry(&struct_chunk, k);
          for (i=0; i<struct_chunk.value_entry_length; ++i)
          {
            x[i] = struct_chunk.value_buf[i];
          }

          GiD_WriteVector(struct_field->coords.size_buf[node_variables.coords_size_Id]+1,
                          x[0], x[1], x[2]);
        }
      }

      GiD_EndResult();
      if (struct_field != NULL)
        destroy_chunk_data(&struct_chunk);
      destroy_chunk_data(&chunk);
    }
  }
  destroy_result_data(&result);

  CCAFREE(fluid_ale_connect);

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief The filter's main function.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  CHAR filename[100];
  PROBLEM_DATA problem;
  INT i;

  init_problem_data(&problem, argc, argv);

  sprintf(filename, "%s.flavia.res", problem.basename);
  GiD_OpenPostResultFile(filename);

  /* Some problem types have to be handled in a special way. The
   * default is to simply write the results as they are. */
  switch (problem.type)
  {
#ifdef D_FSI
  case prb_fsi:
    /*
     * We know there is an ale, a fluid and a structure field. We need
     * to find the fluid nodes that match the ale nodes to be able to
     * output fluid node displacements. */
    dsassert(problem.num_discr==3, "three fields expected for fsi");
    write_fsi(&problem);
    break;
  case prb_fluid:
    if (problem.num_discr==2)
    {
      /*
       * So we have two discretizations. If these are two fields we
       * have a fluid-ale problem. Otherwise it might be the infamous
       * projection method or something. */
      if (map_symbol_count(&(problem.control_table), "field") == 2)
      {
        write_fsi(&problem);
        break;
      }
    }
#endif
    /* No break here! The simple fluid case is done by the default
     * code. */
  default:
    for (i=0; i<problem.num_discr; ++i)
    {
#if 0
#ifdef D_SHELL9
      if (problem.discr[i].is_shell9_problem)
      {
        write_shell9_field(&problem, i);
      }
      else
#endif
#endif
        write_field(&problem, i);
    }
  }

  GiD_ClosePostResultFile();

  post_log(4, "Done.\n");
  return 0;
}
