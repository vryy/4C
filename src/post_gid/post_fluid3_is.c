
#ifdef D_FLUID3_IS

#include "post_fluid3_is.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 01/06
*/
/*----------------------------------------------------------------------*/
void fluid3_is_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("fluid3_is_write_gauss");
#endif

  if (gid->is_fluid3_is_222)
  {
    GiD_BeginGaussPoint(gid->fluid3_is_222_name, GiD_Hexahedra, gid->fluid3_is_222_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_fluid3_is_333)
  {
    GiD_BeginGaussPoint(gid->fluid3_is_333_name, GiD_Hexahedra, gid->fluid3_is_333_name, 27, 0, 1);
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
  \date 01/06
*/
/*----------------------------------------------------------------------*/
void fluid3_is_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("fluid3_is_write_mesh");
#endif

  if (gid->is_fluid3_is_222)
  {
    INT i;

    GiD_BeginMesh(gid->fluid3_is_222_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_fluid3_is || numnp != 8) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid3_is_333)
  {
    INT i;

    GiD_BeginMesh(gid->fluid3_is_333_name, 3, GiD_Hexahedra, 27);

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

      if (el_type != el_fluid3_is || numnp != 27) continue;

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
