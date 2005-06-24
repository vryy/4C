
#ifdef D_FLUID2

#include "post_fluid2.h"
#include "gid_out.h"
#include "post_gid.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void fluid2_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("fluid2_write_gauss");
#endif

  if (gid->is_fluid2_22)
  {
    GiD_BeginGaussPoint(gid->fluid2_22_name, GiD_Quadrilateral, gid->fluid2_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_fluid2_33)
  {
    GiD_BeginGaussPoint(gid->fluid2_33_name, GiD_Quadrilateral, gid->fluid2_33_name, 9, 0, 1);
    GiD_EndGaussPoint();
  }

  if (gid->is_fluid2_3)
  {
    /* The number of gauss points here is just guessing. It's not
     * important as long as there are no element values. */
    GiD_BeginGaussPoint(gid->fluid2_3_name, GiD_Triangle, gid->fluid2_3_name, 3, 0, 1);
    GiD_EndGaussPoint();
  }

  if (gid->is_fluid2_6)
  {
    /* The number of gauss points here is just guessing. It's not
     * important as long as there are no element values. */
    GiD_BeginGaussPoint(gid->fluid2_6_name, GiD_Triangle, gid->fluid2_6_name, 6, 0, 1);
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
void fluid2_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("fluid2_write_mesh");
#endif

  if (gid->is_fluid2_22)
  {
    INT i;

    GiD_BeginMesh(gid->fluid2_22_name, 2, GiD_Quadrilateral, 4);

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

      if (el_type != el_fluid2 || numnp != 4) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid2_33)
  {
    INT i;

    GiD_BeginMesh(gid->fluid2_33_name, 2, GiD_Quadrilateral, 9);

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

      if (el_type != el_fluid2 || numnp != 9) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid2_3)
  {
    INT i;

    GiD_BeginMesh(gid->fluid2_3_name, 2, GiD_Triangle, 3);

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

      if (el_type != el_fluid2 || numnp != 3) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid2_6)
  {
    INT i;

    GiD_BeginMesh(gid->fluid2_6_name, 2, GiD_Triangle, 6);

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

      if (el_type != el_fluid2 || numnp != 6) continue;

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
