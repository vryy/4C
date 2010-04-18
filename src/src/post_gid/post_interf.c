
#ifdef D_INTERF

#include "post_interf.h"
#include "post_gid.h"
#include "gid_out.h"

#include "../interf/interf.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
static void interf_write_22_stress(FIELD_DATA* field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE step)
{
  CHAR* interf_stresstypenames[] = INTERF_STRESSTYPE;
  TRANSLATION_TABLE interf_stresstype;
  INT stresstype;
  INT el_type;
  INT k;

  CHAR* tn_componentnames[] = { "stress-tang","stress-normal","dummy","D-nomal","D-tang","dummy" };
  CHAR* xy_componentnames[] = { "stress-sxx","stress-syy","stress-sxy","D-nomal","D-tang","dummy" };
  CHAR** componentnames;

#ifdef DEBUG
  dstrc_enter("interf_write_22_stress");
#endif

  init_translation_table(&interf_stresstype,
                         map_read_map(field->group, "interf_stresstypes"),
                         interf_stresstypenames);

  /* look at the first interface element */
  for (k=0; k<field->numele; ++k)
  {
    INT Id;
    INT numnp;
    INT dis;

    /* read the element's data */
    get_element_params(field, k, &Id, &el_type, &dis, &numnp);

    if (el_type == el_interf)
      break;
  }

  if (k==field->numele)
  {
    dserror("no interface elements found");
  }

  stresstype = field->ele_param.size_buf[interf_variables.ep_size_stresstyp];
  stresstype = interf_stresstype.table[stresstype];

  switch (stresstype)
  {
  case if_xy:
    componentnames = xy_componentnames;
    break;
  case if_tn:
    componentnames = tn_componentnames;
    break;
  default:
    dserror("stress type %d unknown", stresstype);
  }

  GiD_BeginResult("interface_stresses", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                  gid->interf_22_name, NULL, 6, componentnames);

  for (k=0; k<field->numele; ++k)
  {
    INT numnp;
    INT Id;
    INT l;
    DOUBLE* stress;
    INT dis;

    /* read the element's data */
    get_element_params(field, k, &Id, &el_type, &dis, &numnp);

    if (el_type != el_interf) continue;

    chunk_read_value_entry(chunk, k);
    stress = chunk->value_buf;

    l = 0;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 1;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 0;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
  }

  GiD_EndResult();
  destroy_translation_table(&interf_stresstype);

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
void interf_write_33_stress(FIELD_DATA* field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE step)
{
  CHAR* interf_stresstypenames[] = INTERF_STRESSTYPE;
  TRANSLATION_TABLE interf_stresstype;
  INT stresstype;
  INT el_type;
  INT k;

  CHAR* tn_componentnames[] = { "stress-tang","stress-normal","dummy","D-nomal","D-tang","dummy" };
  CHAR* xy_componentnames[] = { "stress-sxx","stress-syy","stress-sxy","D-nomal","D-tang","dummy" };
  CHAR** componentnames;

#ifdef DEBUG
  dstrc_enter("interf_write_33_stress");
#endif

  init_translation_table(&interf_stresstype,
                         map_read_map(field->group, "interf_stresstypes"),
                         interf_stresstypenames);

  /* look at the first interface element */
  for (k=0; k<field->numele; ++k)
  {
    chunk_read_size_entry(&(field->ele_param), k);

    el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];
    el_type = field->problem->element_type.table[el_type];
    if (el_type == el_interf)
      break;
  }

  if (k==field->numele)
  {
    dserror("no interface elements found");
  }

  stresstype = field->ele_param.size_buf[interf_variables.ep_size_stresstyp];
  stresstype = interf_stresstype.table[stresstype];

  switch (stresstype)
  {
  case if_xy:
    componentnames = xy_componentnames;
    break;
  case if_tn:
    componentnames = tn_componentnames;
    break;
  default:
    dserror("stress type %d unknown", stresstype);
  }

  GiD_BeginResult("interface_stresses", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                  gid->interf_22_name, NULL, 6, componentnames);

  for (k=0; k<field->numele; ++k)
  {
    INT numnp;
    INT Id;
    INT l;
    DOUBLE* stress;
    INT dis;

    /* read the element's data */
    get_element_params(field, k, &Id, &el_type, &dis, &numnp);

    if (el_type != el_interf) continue;

    chunk_read_value_entry(chunk, k);
    stress = chunk->value_buf;

    l = 0;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 2;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 0;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 1;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 2;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 1;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 0;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
    l = 1;
    GiD_Write3DMatrix(Id+1,
                      stress[5*l+0], stress[5*l+1], stress[5*l+2],
                      stress[5*l+3], stress[5*l+4], 0);
  }

  GiD_EndResult();
  destroy_translation_table(&interf_stresstype);

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
void interf_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("interf_write_stress");
#endif

  /* 1D interface element (combination only with wall) */
  if (gid->is_interf_22)
  {
    interf_write_22_stress(field, gid, chunk, step);
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_interf_33)
  {
    interf_write_33_stress(field, gid, chunk, step);
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
void interf_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("interf_write_gauss");
#endif

  if (gid->is_interf_22)
  {
    GiD_BeginGaussPoint(gid->interf_22_name, GiD_Quadrilateral, gid->interf_22_name, 4, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_interf_33)
  {
    GiD_BeginGaussPoint(gid->interf_33_name, GiD_Quadrilateral, gid->interf_33_name, 9, 0, 1);
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
void interf_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("interf_write_mesh");
#endif

  if (gid->is_interf_22)
  {
    INT i;
    INT numnp;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->interf_22_name, 2, GiD_Quadrilateral, numnp);

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

      if (el_type != el_interf) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_interf_33)
  {
    INT i;
    INT numnp;

    /*- NNODE is only checked for first element,if you use different wall-types this will be a problem -*/
    chunk_read_size_entry(&(field->ele_param), 0);
    numnp   = field->ele_param.size_buf[element_variables.ep_size_numnp];

    GiD_BeginMesh(gid->interf_33_name, 2, GiD_Quadrilateral, numnp);

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

      if (el_type != el_interf) continue;

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
