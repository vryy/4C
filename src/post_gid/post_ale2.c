
#ifdef D_ALE

#include "post_ale2.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void ale2_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("ale2_write_domain");
#endif

  if (gid->is_ale_11)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->ale_11_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_ale2 || numnp != 4) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<4; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_ale_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->ale_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_ale2 || numnp != 4) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<4; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_ale_tri_1)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->ale_tri_1_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_ale2 || numnp != 3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<3; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_ale_tri_3)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->ale_tri_3_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_ale2 || numnp != 3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<3; j++)
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
void ale2_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("ale2_write_gauss");
#endif

  if (gid->is_ale_11)
  {
    GiD_BeginGaussPoint(gid->ale_11_name, GiD_Quadrilateral, gid->ale_11_name, 1, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_ale_22)
  {
    GiD_BeginGaussPoint(gid->ale_22_name, GiD_Quadrilateral, gid->ale_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_ale_tri_1)
  {
    GiD_BeginGaussPoint(gid->ale_tri_1_name, GiD_Triangle, gid->ale_tri_1_name, 1, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_ale_tri_3)
  {
    GiD_BeginGaussPoint(gid->ale_tri_3_name, GiD_Triangle, gid->ale_tri_3_name, 3, 0, 1);
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
void ale2_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("ale2_write_mesh");
#endif

  /* In fsi calculations we don't want to show the ale elements */
  if (field->problem->type != prb_fsi)
  {
    if (gid->is_ale_11)
    {
      INT i;

      GiD_BeginMesh(gid->ale_11_name, 2, GiD_Quadrilateral, 4);

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

        if (el_type != el_ale2 || numnp != 4) continue;

        chunk_read_size_entry(&(field->mesh), i);
        get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
        GiD_WriteElement(Id+1, mesh_entry);
      }

      GiD_EndElements();
      GiD_EndMesh();
    }

    /*--------------------------------------------------------------------*/

    if (gid->is_ale_22)
    {
      INT i;

      GiD_BeginMesh(gid->ale_22_name, 2, GiD_Quadrilateral, 4);

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

        if (el_type != el_ale2 || numnp != 4) continue;

        chunk_read_size_entry(&(field->mesh), i);
        get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
        GiD_WriteElement(Id+1, mesh_entry);
      }

      GiD_EndElements();
      GiD_EndMesh();
    }

    /*--------------------------------------------------------------------*/

    if (gid->is_ale_tri_1)
    {
      INT i;

      GiD_BeginMesh(gid->ale_tri_1_name, 2, GiD_Triangle, 3);

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

        if (el_type != el_ale2 || numnp != 3) continue;

        chunk_read_size_entry(&(field->mesh), i);
        get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
        GiD_WriteElement(Id+1, mesh_entry);
      }

      GiD_EndElements();
      GiD_EndMesh();
    }

    /*--------------------------------------------------------------------*/

    if (gid->is_ale_tri_3)
    {
      INT i;

      GiD_BeginMesh(gid->ale_tri_3_name, 2, GiD_Triangle, 3);

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

        if (el_type != el_ale2 || numnp != 3) continue;

        chunk_read_size_entry(&(field->mesh), i);
        get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
        GiD_WriteElement(Id+1, mesh_entry);
      }

      GiD_EndElements();
      GiD_EndMesh();
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
