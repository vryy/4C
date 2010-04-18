
#ifdef D_SHELL9

#include "post_shell9.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell9_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("shell9_write_domain");
#endif

  /*4-noded*/
  if (gid->is_shell9_4_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_4_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 4) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<4; j++) */    /* quadrilateral version */
      for (j=0; j<8; j++)       /* hexahedra version */
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_shell9_4_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_4_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 4) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<9; j++) */     /* quadrilateral version */
      for (j=0; j<27; j++)       /* hexahedra version */
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  /*8-noded*/
  if (gid->is_shell9_8_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_8_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 8) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<4; j++) */     /* quadrilateral version */
      for (j=0; j<8; j++)        /* hexahedra version */
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_shell9_8_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_8_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 8) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<9; j++) */    /* quadrilateral version */
      for (j=0; j<27; j++)      /* hexahedra version */
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  /*9-noded*/
  if (gid->is_shell9_9_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_9_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 9) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<4; j++) */     /* quadrilateral version */
      for (j=0; j<8; j++)        /* hexahedra version */
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_shell9_9_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->shell9_9_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_shell9 || numnp != 9) continue;

      chunk_read_size_entry(chunk, i);

      /* for (j=0; j<9; j++) */      /* quadrilateral version */
      for (j=0; j<27; j++)      /* hexahedra version */
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
void shell9_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("shell9_write_stress");
#endif

  if (gid->is_shell9_4_22 || gid->is_shell9_8_22 || gid->is_shell9_9_22 ||
      gid->is_shell9_4_33 || gid->is_shell9_8_33 || gid->is_shell9_9_33)
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
void shell9_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("shell9_write_gauss");
#endif

  /* this is the shell9 visualization using Hexahedra */
  if (gid->is_shell9_4_22)
  {
    GiD_BeginGaussPoint(gid->shell9_4_22_name, GiD_Hexahedra, gid->shell9_4_22_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell9_4_33)
  {
    GiD_BeginGaussPoint(gid->shell9_4_33_name, GiD_Hexahedra, gid->shell9_4_33_name, 27, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell9_8_22)
  {
    GiD_BeginGaussPoint(gid->shell9_8_22_name, GiD_Hexahedra, gid->shell9_8_22_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell9_8_33)
  {
    GiD_BeginGaussPoint(gid->shell9_8_33_name, GiD_Hexahedra, gid->shell9_8_33_name, 27, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell9_9_22)
  {
    GiD_BeginGaussPoint(gid->shell9_9_22_name, GiD_Hexahedra, gid->shell9_9_22_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_shell9_9_33)
  {
    GiD_BeginGaussPoint(gid->shell9_9_33_name, GiD_Hexahedra, gid->shell9_9_33_name, 27, 0, 1);
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
static void shell9_write_coords(FIELD_DATA* field, GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("shell9_write_coords");
#endif

  dserror("not yet");

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
void shell9_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("shell9_write_mesh");
#endif

  /* this is the shell9 visualization using Hexahedra */

  /* 4-noded shell9 element -> Hex8 */
  if (gid->is_shell9_4_22 || gid->is_shell9_4_33)
  {
    INT i;

    if (gid->is_shell9_4_22)
    {
      GiD_BeginMesh(gid->shell9_4_22_name, 3, GiD_Hexahedra, 8);
    }
    else
    {
      GiD_BeginMesh(gid->shell9_4_33_name, 3, GiD_Hexahedra, 8);
    }

    if (*first_mesh)
    {
      *first_mesh = 0;
      shell9_write_coords(field, gid);
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

      if (el_type != el_shell9 || numnp != 4) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }
    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  /* 8-noded shell9 element -> Hex20 */
  if (gid->is_shell9_8_22 || gid->is_shell9_8_33)
  {
    INT i;

    if (gid->is_shell9_8_22)
    {
      GiD_BeginMesh(gid->shell9_8_22_name, 3, GiD_Hexahedra, 20);
    }
    else
    {
      GiD_BeginMesh(gid->shell9_8_33_name, 3, GiD_Hexahedra, 20);
    }

    if (*first_mesh)
    {
      *first_mesh = 0;
      shell9_write_coords(field, gid);
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

      if (el_type != el_shell9 || numnp != 8) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  /* 9-noded shell9 element -> Hex27 */
  if (gid->is_shell9_9_22 || gid->is_shell9_9_33)
  {
    INT i;

    if (gid->is_shell9_9_22)
    {
      GiD_BeginMesh(gid->shell9_9_22_name, 3, GiD_Hexahedra, 27);
    }
    else
    {
      GiD_BeginMesh(gid->shell9_9_33_name, 3, GiD_Hexahedra, 27);
    }

    if (*first_mesh)
    {
      *first_mesh = 0;
      shell9_write_coords(field, gid);
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

      if (el_type != el_shell9 || numnp != 9) continue;

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
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void shell9_write_displacement(FIELD_DATA *field, RESULT_DATA* result)
{
#ifdef DEBUG
  dstrc_enter("shell9_write_displacement");
#endif

  dserror("not yet");

#ifdef DEBUG
  dstrc_exit();
#endif
}


#endif
