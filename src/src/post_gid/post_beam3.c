
#ifdef D_BEAM3

#include "post_beam3.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void beam3_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("beam3_write_domain");
#endif

  if (gid->is_beam3_21)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->beam3_21_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<1; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->beam3_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<2; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_32)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->beam3_32_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<2; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->beam3_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

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
void beam3_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("beam3_write_stress");
#endif

  if (gid->is_beam3_21)
  {
    CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
    INT i;

    GiD_BeginResult("beam3_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->beam3_21_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      GiD_Write3DMatrix(Id+1,
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], stress[5]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_22)
  {
    CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
    INT i;

    GiD_BeginResult("beam3_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->beam3_22_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      GiD_Write3DMatrix(Id+1,
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], stress[5]);
      GiD_Write3DMatrix(Id+1,
                        stress[6], stress[7], stress[8],
                        stress[9], stress[10], stress[11]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_32)
  {
    CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
    INT i;

    GiD_BeginResult("beam3_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->beam3_32_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      GiD_Write3DMatrix(Id+1,
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], stress[5]);
      GiD_Write3DMatrix(Id+1,
                        stress[6], stress[7], stress[8],
                        stress[9], stress[10], stress[11]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_33)
  {
    CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
    INT i;

    GiD_BeginResult("beam3_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->beam3_33_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_beam3) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      GiD_Write3DMatrix(Id+1,
                        stress[0], stress[1], stress[2],
                        stress[3], stress[4], stress[5]);
      GiD_Write3DMatrix(Id+1,
                        stress[6], stress[7], stress[8],
                        stress[9], stress[10], stress[11]);
      GiD_Write3DMatrix(Id+1,
                        stress[12], stress[13], stress[14],
                        stress[15], stress[16], stress[17]);
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
void beam3_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("beam3_write_gauss");
#endif

  if (gid->is_beam3_21)
  {
    GiD_BeginGaussPoint(gid->beam3_21_name, GiD_Linear, gid->beam3_21_name, 1, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_beam3_22)
  {
    GiD_BeginGaussPoint(gid->beam3_22_name, GiD_Linear, gid->beam3_22_name, 2, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_beam3_32)
  {
    GiD_BeginGaussPoint(gid->beam3_32_name, GiD_Linear, gid->beam3_32_name, 2, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_beam3_33)
  {
    GiD_BeginGaussPoint(gid->beam3_33_name, GiD_Linear, gid->beam3_33_name, 3, 0, 1);
    GiD_EndGaussPoint();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


void beam3_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("beam3_write_mesh");
#endif

  if (gid->is_beam3_21)
  {
    INT i;
    GiD_BeginMesh(gid->beam3_21_name, 3, GiD_Linear, 2);

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

      if (el_type != el_beam3 || numnp != 2) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_22)
  {
    INT i;

    GiD_BeginMesh(gid->beam3_22_name, 3, GiD_Linear, 2);

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

      if (el_type != el_beam3 || numnp != 2) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_32)
  {
    INT i;

    GiD_BeginMesh(gid->beam3_32_name, 3, GiD_Linear, 3);

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

      if (el_type != el_beam3 || numnp != 3) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_beam3_33)
  {
    INT i;

    GiD_BeginMesh(gid->beam3_33_name, 3, GiD_Linear, 3);

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

      if (el_type != el_beam3 || numnp != 3) continue;

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
