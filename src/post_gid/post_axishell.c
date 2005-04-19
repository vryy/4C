
#ifdef D_AXISHELL

#include "post_axishell.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void axishell_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("axishell_write_stress");
#endif

  if (gid->is_axishell)
  {
    CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
    INT i;

    GiD_BeginResult("axishell_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->axishell_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_axishell || numnp !=2) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      GiD_Write3DMatrix(Id+1,
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], 0);
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
void axishell_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("axishell_write_gauss");
#endif

  if (gid->is_axishell)
  {
    GiD_BeginGaussPoint(gid->axishell_name, GiD_Linear, gid->axishell_name, 1, 0, 1);
    GiD_EndGaussPoint();
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
void axishell_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("axishell_write_mesh");
#endif

  if (gid->is_axishell)
  {
    INT i;

    GiD_BeginMesh(gid->axishell_name, 2, GiD_Linear, 2);

    if (*first_mesh)
    {
      *first_mesh = 0;
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

      if (el_type != el_axishell || numnp != 2) continue;

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


#endif
