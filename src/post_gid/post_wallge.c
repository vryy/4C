
#ifdef D_WALLGE

#include "post_wallge.h"
#include "post_gid.h"
#include "gid_out.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Write stress.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void wallge_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
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
  dstrc_enter("wallge_write_stress");
#endif

  /* gradient enhanced wall element */
  if (gid->is_wallge_22)
  {
    INT k;
    INT ngauss=4;
    CHAR* componentnames[] = { "Stress-xx","Stress-yy","Stress-xy","damage","loc_aequiv_strain","nonloc_aequiv_strain" };
    GiD_BeginResult("wallge_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->wallge_22_name, NULL, 6, componentnames);

    for (k=0; k<field->numele; ++k)
    {
      INT numnp;
      INT Id;
      INT l;
      INT el_type;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, k, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wallge) continue;

      chunk_read_value_entry(chunk, k);
      stress = chunk->value_buf;

      for (l=0; l<ngauss; ++l)
      {
        INT p;
        p = gaussperm4[l];
        GiD_Write3DMatrix(Id+1,
                          stress[6*p+0], stress[6*p+1], stress[6*p+2],
                          stress[6*p+3], stress[6*p+4], stress[6*p+5]);
      }
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wallge_33)
  {
    INT k;
    INT ngauss=9;
    CHAR* componentnames[] = { "Stress-xx","Stress-yy","Stress-xy","damage","loc_aequiv_strain","nonloc_aequiv_strain" };
    GiD_BeginResult("wallge_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                    gid->wallge_33_name, NULL, 6, componentnames);

    for (k=0; k<field->numele; ++k)
    {
      INT numnp;
      INT Id;
      INT l;
      INT el_type;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, k, &Id, &el_type, &dis, &numnp);

      if (el_type != el_wallge) continue;

      chunk_read_value_entry(chunk, k);
      stress = chunk->value_buf;

      for (l=0; l<ngauss; ++l)
      {
        INT p;
        p = gaussperm9[l];
        GiD_Write3DMatrix(Id+1,
                          stress[6*p+0], stress[6*p+1], stress[6*p+2],
                          stress[6*p+3], stress[6*p+4], stress[6*p+5]);
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
  \brief Write gauss point info.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void wallge_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("wallge_write_gauss");
#endif

  if (gid->is_wallge_22)
  {
    GiD_BeginGaussPoint(gid->wallge_22_name, GiD_Quadrilateral, gid->wallge_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_wallge_33)
  {
    GiD_BeginGaussPoint(gid->wallge_33_name, GiD_Quadrilateral, gid->wallge_33_name, 9, 0, 1);
    GiD_EndGaussPoint();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write mesh connectivity.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void wallge_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("wallge_write_mesh");
#endif

  if (gid->is_wallge_22)
  {
    INT i;
    INT numnp;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->wallge_22_name, 2, GiD_Quadrilateral, numnp);

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

      if (el_type != el_wallge) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_wallge_33)
  {
    INT i;
    INT numnp;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->wallge_33_name, 2, GiD_Quadrilateral, numnp);

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

      if (el_type != el_wallge) continue;

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
