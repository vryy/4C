
#ifdef D_WALL1

#include "post_wall1.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void wall1_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("wall1_write_domain");
#endif

  if (gid->is_wall1_11)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->wall1_11_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1 || numnp != 3) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<1; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wall1_22)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->wall1_22_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<4; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wall1_33)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->wall1_33_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

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
void wall1_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
/*
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/

  INT           gaussperm4[4] = {3,1,0,2};
  /*INT           gaussperm8[8] = {0,4,2,6,1,5,3,7};*/
  INT           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
  /*INT           gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};*/

#ifdef DEBUG
  dstrc_enter("wall1_write_stress");
#endif

  /* wall1 element with 2x2 gausspoints */
  if (gid->is_wall1_22)
  {
    INT ngauss = 4;
    INT i;

    /* 2D plane stress - plane strain element */
    CHAR* componentnames[] = { "Stress-xx", "Stress-yy", "Stress-xy",
                               "Stress-zz", "damage", "||u||" };
    GiD_BeginResult("wall1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->wall1_22_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT l;
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      for (l=0; l<ngauss; ++l)
      {
        INT p;
        p = gaussperm4[l];
        GiD_Write3DMatrix(Id+1,
                          stress[9*p+0], stress[9*p+1], stress[9*p+2],
                          stress[9*p+3], stress[9*p+7], stress[9*p+8]);
      }
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  /* wall1 triangle-element with 1x1 gausspoint */
  if (gid->is_wall1_11)
  {
    INT ngauss = 1;
    INT i;

    /* 2D plane stress - plane strain element */
    CHAR* componentnames[] = { "Stress-xx", "Stress-yy", "Stress-xy",
                               "Stress-zz", "damage", "||u||" };
    GiD_BeginResult("wall1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->wall1_11_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT l;
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      for (l=0; l<ngauss; ++l)
      {
        /* stupid, but similar to the other cases */
        INT p;
        p = l;
        GiD_Write3DMatrix(Id+1,
                          stress[9*p+0], stress[9*p+1], stress[9*p+2],
                          stress[9*p+3], 0, 0);
      }
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  /* wall1 element with 3x3 gausspoints */
  if (gid->is_wall1_33)
  {
    INT ngauss = 9;
    INT i;

    /* 2D plane stress - plane strain element */
    CHAR* componentnames[] = { "Stress-xx", "Stress-yy", "Stress-xy",
                               "Stress-zz", "damage", "||u||" };
    GiD_BeginResult("wall1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->wall1_33_name, NULL, 6, componentnames);

    for (i=0; i<field->numele; ++i)
    {
      INT l;
      INT numnp;
      INT el_type;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      for (l=0; l<ngauss; ++l)
      {
        INT p;
        p = gaussperm9[l];
        GiD_Write3DMatrix(Id+1,
                          stress[9*p+0], stress[9*p+1], stress[9*p+2],
                          stress[9*p+3], stress[9*p+7], stress[9*p+8]);
      }
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
void wall1_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("wall1_write_gauss");
#endif

  if (gid->is_wall1_11)
  {
    GiD_BeginGaussPoint(gid->wall1_11_name, GiD_Triangle, gid->wall1_11_name, 1, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_wall1_22)
  {
    GiD_BeginGaussPoint(gid->wall1_22_name, GiD_Quadrilateral, gid->wall1_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_wall1_33)
  {
    GiD_BeginGaussPoint(gid->wall1_33_name, GiD_Quadrilateral, gid->wall1_33_name, 9, 0, 1);
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
void wall1_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("wall1_write_mesh");
#endif

  if (gid->is_wall1_11)
  {
    INT i;

    GiD_BeginMesh(gid->wall1_11_name, 2, GiD_Triangle, 3);

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

      if (el_type != el_wall1 || numnp != 3) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wall1_22)
  {
    INT numnp;
    INT i;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->wall1_22_name, 2, GiD_Quadrilateral, numnp);

    if (*first_mesh)
    {
      *first_mesh = 0;
      write_coords(field, gid);
    }

    GiD_BeginElements();

    for (i=0; i<field->numele; ++i)
    {
      INT el_type;
      INT Id;
      INT mesh_entry[MAXNOD];
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wall1_33)
  {
    INT i;
    INT numnp;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->wall1_33_name, 2, GiD_Quadrilateral, numnp);

    if (*first_mesh)
    {
      *first_mesh = 0;
      write_coords(field, gid);
    }

    GiD_BeginElements();

    for (i=0; i<field->numele; ++i)
    {
      INT el_type;
      INT Id;
      INT mesh_entry[MAXNOD];
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wall1 || numnp < 8) continue;

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
