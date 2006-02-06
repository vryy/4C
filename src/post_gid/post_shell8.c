
#ifdef D_SHELL8

#include "post_shell8.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("shell8_write_domain");
#endif

  if (gid->is_shell8_4_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell8_4_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell8 || numnp != 4) continue;

      chunk_read_size_entry(chunk, i);

#if 0
      for (j=0; j<4; j++)/* quadrilateral version */
#endif
        for (j=0; j<8; j++)/* hexahedra version */
          GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_shell8_8_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell8_8_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell8 || numnp != 9) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<9; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("shell8_write_stress");
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("shell8_write_gauss");
#endif

  if (gid->is_shell8_4_22)
  {
    /*GiD_BeginGaussPoint(gid->shell8_22_name, GiD_Quadrilateral, gid->shell8_22_name, 4, 0, 1);*/

    /* this is the shell visualization using Hexahedra */
    GiD_BeginGaussPoint(gid->shell8_4_22_name, GiD_Hexahedra, gid->shell8_4_22_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell8_8_33)
  {
    GiD_BeginGaussPoint(gid->shell8_8_33_name, GiD_Quadrilateral, gid->shell8_8_33_name, 9, 0, 1);
    GiD_EndGaussPoint();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read the director vector of on particular node.

  Directors are element values. That is each shell8 element has an
  director in each node. Thus we need to find the elements attached
  to the node and read the directors there.

  We just use the first element's director.

  \param field          (i) the field everything belongs to
  \param node           (i) the node number (field local id)
  \param director_chunk (i) a properly initialized chunk containing
                            the director values
  \param d              (o) the director of this node, scaled
  \param thick          (o) the thick of this node

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void shell8_read_director(FIELD_DATA* field,
                                 INT node,
                                 CHUNK_DATA* director_chunk,
                                 DOUBLE d[3],
                                 DOUBLE* thick)
{
  INT j;
  INT ele_id;
  INT numnp;
  DOUBLE scal = 1.0;

#ifdef DEBUG
  dstrc_enter("shell8_read_director");
#endif

  chunk_read_size_entry(&(field->coords), node);

  /* We do it the ccarat way and take the director of the first
   * element connected here. This requires that all directors at one
   * node are the same. If that's not true an average might be
   * better. Ask someone who knows the shell8 element. */

  /*
   * Please note that the first entry is the node Id, the second the
   * number of elements here. So we have to access the third. */
  ele_id = field->coords.size_buf[node_variables.coords_size_eleid];

  /* find the position of this node in the mesh */
  chunk_read_size_entry(&(field->ele_param), ele_id);
  chunk_read_size_entry(&(field->mesh), ele_id);

  numnp = field->ele_param.size_buf[element_variables.ep_size_numnp];
  for (j=0; j<numnp; ++j)
  {
    if (field->mesh.size_buf[j] == node)
    {
      break;
    }
  }

  if (j==numnp)
  {
    dserror("node not found: connectivity breakdown");
  }

  /*
   * shell8_director is an element chunk. So we access its entries
   * by element id (field local). */
  chunk_read_value_entry(director_chunk, ele_id);

  /* read the values, scale them properly */
  *thick = director_chunk->value_buf[4*j+3];
  d[0]   = director_chunk->value_buf[4*j+0]*(*thick)*scal/2;
  d[1]   = director_chunk->value_buf[4*j+1]*(*thick)*scal/2;
  d[2]   = director_chunk->value_buf[4*j+2]*(*thick)*scal/2;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write coordinates for shell8 elements as bricks.

  The normal functions are sufficient to visualize the middle
  layer. However, there is a third dimension.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void shell8_write_coords(FIELD_DATA* field, GIDSET* gid)
{
  INT i;
  INT j;
  DOUBLE x[3];
#if defined(S8_HEX8)
  DOUBLE d[3];
  DOUBLE thick;
#endif
  CHUNK_DATA director_chunk;

#ifdef DEBUG
  dstrc_enter("shell8_write_coords");
#endif

  init_chunk_data(&(field->head), &director_chunk, "shell8_director");
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

#if defined(S8_HEX8)

    shell8_read_director(field, i, &director_chunk, d, &thick);

    /* All GiD ids are one based. */
    /* This numbering scheme assumes that there is only one
     * discretization in this problem and that there are only shell8
     * elements used. */
    GiD_WriteCoordinates(Id+1,
                         x[0]-d[0], x[1]-d[1], x[2]-d[2]);
    GiD_WriteCoordinates(Id+1 + field->problem->numnp,
                         x[0]+d[0], x[1]+d[1], x[2]+d[2]);

#else
    GiD_WriteCoordinates(Id+1, x[0], x[1], x[2]);
#endif
  }

  GiD_EndCoordinates();
  destroy_chunk_data(&director_chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("shell8_write_mesh");
#endif

  if (gid->is_shell8_4_22)
  {
    INT i;

#if defined(S8_HEX8)
    PROBLEM_DATA* problem = field->problem;

    /* this is the hexahedra version */
    GiD_BeginMesh(gid->shell8_4_22_name, 3, GiD_Hexahedra, 8);
    if (*first_mesh)
    {
      *first_mesh = 0;
      shell8_write_coords(field, gid);
    }

    GiD_BeginElements();

    for (i=0; i<field->numele; ++i)
    {
      INT k;
      INT numnp;
      INT el_type;
      INT Id;
      INT mesh_entry[MAXNOD];
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell8 || numnp !=4) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);

      /* There is a second row of nodes. These are invented by the
       * filter. We relay on beeing the only discretization of this
       * problem. */
      /* ok. This is a particularly nasty hack. */
      get_gid_node_ids(field, field->mesh.size_buf, &(mesh_entry[numnp]), numnp);
      for (k=0; k<numnp; ++k)
      {
        mesh_entry[numnp+k] += problem->numnp;
      }

      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();

#else

    /* this is the surface version */
    GiD_BeginMesh(gid->shell8_4_22_name, 3, GiD_Quadrilateral, 4);
    if (*first_mesh)
    {
      *first_mesh = 0;
      shell8_write_coords(field, gid);
    }

    GiD_BeginElements();

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT mesh_entry[MAXNOD];
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell8 || numnp !=4) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);

      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();

#endif
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_shell8_8_33)
  {
    INT i;

    GiD_BeginMesh(gid->shell8_8_33_name, 3, GiD_Quadrilateral, 9);
    if (first_mesh)
    {
      first_mesh = 0;
      write_coords(field, gid);
    }

    GiD_BeginElements();

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT mesh_entry[MAXNOD];
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell8 || numnp !=9) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write displacements for shell8 elements as bricks.

  The normal functions are sufficient to visualize the middle
  layer. However, there is a third dimension.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell8_write_displacement(FIELD_DATA *field, RESULT_DATA* result)
{
  INT k;
  DOUBLE x[6];
  CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };
  CHUNK_DATA chunk;
  DOUBLE time;
  INT step;
  CHAR buf[100];
  DOUBLE scal = 1.0;
  DOUBLE sdc;

#ifdef DEBUG
  dstrc_enter("shell8_write_displacement");
#endif

  init_chunk_data(result, &chunk, "displacement");

  time = map_read_real(result->group, "time");
  step = map_read_int(result->group, "step");

  post_log(3, "%s: Write displacement of step %d\n",
          field->name, step);

  /* Again we look at the first element only. */
  chunk_read_value_entry(&(field->ele_param), 0);
  sdc = field->ele_param.value_buf[shell8_variables.ep_value_sdc];

  sprintf(buf, "%s_displacement", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Vector, GiD_OnNodes,
                  NULL, NULL, 3, componentnames);

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k)
  {
    INT i;
    INT Id;

    chunk_read_value_entry(&chunk, k);
    for (i=0; i<chunk.value_entry_length; ++i)
    {
      x[i] = chunk.value_buf[i];
    }

    chunk_read_size_entry(&(field->coords), k);
    Id = field->coords.size_buf[node_variables.coords_size_Id];

    GiD_WriteVector(Id+1,
                    x[0]-x[3]*scal/sdc, x[1]-x[4]*scal/sdc, x[2]-x[5]*scal/sdc);
    GiD_WriteVector(Id+1 + field->problem->numnp,
                    x[0]+x[3]*scal/sdc, x[1]+x[4]*scal/sdc, x[2]+x[5]*scal/sdc);
  }

  GiD_EndResult();
  destroy_chunk_data(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
