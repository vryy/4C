
#ifdef D_FLUID3

#include "post_fluid3.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void fluid3_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("fluid3_write_domain");
#endif

  if (gid->is_fluid3_222)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->fluid3_222_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_fluid3 || numnp != 8) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<8; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid3_333)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->fluid3_333_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_fluid3 || (numnp != 20 && numnp != 27)) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<numnp; j++)
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
void fluid3_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("fluid3_write_stress");
#endif

  /* Do we want to output fluid stresses? */

  if (gid->is_fluid3_222)
  {
  }

  if (gid->is_fluid3_333)
  {
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
void fluid3_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("fluid3_write_gauss");
#endif

  if (gid->is_fluid3_222)
  {
    GiD_BeginGaussPoint(gid->fluid3_222_name, GiD_Hexahedra, gid->fluid3_222_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_fluid3_333)
  {
    GiD_BeginGaussPoint(gid->fluid3_333_name, GiD_Hexahedra, gid->fluid3_333_name, 27, 0, 1);
    GiD_EndGaussPoint();
  }

  if (gid->is_fluid3_tet4)
  {
    GiD_BeginGaussPoint(gid->fluid3_tet4_name, GiD_Tetrahedra, gid->fluid3_tet4_name, 1, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_fluid3_tet10)
  {
    GiD_BeginGaussPoint(gid->fluid3_tet10_name, GiD_Tetrahedra, gid->fluid3_tet10_name, 4, 0, 1);
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
void fluid3_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("fluid3_write_mesh");
#endif

  if (gid->is_fluid3_222)
  {
    INT i;

    GiD_BeginMesh(gid->fluid3_222_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_fluid3 || numnp != 8) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid3_333)
  {
    INT i;

    /*GiD_BeginMesh(gid->fluid3_333_name, 3, GiD_Hexahedra, 27);*/
    GiD_BeginMesh(gid->fluid3_333_name, 3, GiD_Hexahedra, 20);

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

      if (el_type != el_fluid3 || (numnp != 20 && numnp != 27)) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid3_tet4)
  {
    INT i;

    GiD_BeginMesh(gid->fluid3_tet4_name, 3, GiD_Tetrahedra, 4);

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

      if (el_type != el_fluid3 || numnp != 4) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_fluid3_tet10)
  {
    INT i;

    GiD_BeginMesh(gid->fluid3_tet10_name, 3, GiD_Tetrahedra, 10);

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

      if (el_type != el_fluid3 || numnp != 10) continue;

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
