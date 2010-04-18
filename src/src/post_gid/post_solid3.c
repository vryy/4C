
#ifdef D_SOLID3

#include "post_solid3.h"
#include "../solid3/solid3.h"
#include "gid_out.h"
#include "post_gid.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author bborn
  \date 01/07
*/
/*----------------------------------------------------------------------*/
void solid3_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("solid3_write_domain");
#endif

  if (gid->is_solid3_h8_222)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->solid3_h8_222_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_solid3 || numnp != 8) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<8; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_solid3_h20_333)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->solid3_h20_333_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_solid3 || numnp != 20 ) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<27; j++)
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

  \author bborn
  \date 01/07
*/
/*----------------------------------------------------------------------*/
void solid3_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
  /* Gauss point permutations */
  INT gperm_h_222[8] = {0,4,6,2,1,5,7,3};
  INT gperm_h_333[27] = {0,10,24,6,2,20,26,8,
                         9,21,15,3,1,19,25,7,11,23,17,5,
                         12,10,22,16,4,14,13};
  /* INT gperm_t_4[4] = {0,1,2,3}; */
  /* INT gperm_t_10[10] = {0,1,2,3,4,5,6,7,8,9}; */ /* check this */
  INT *gperm;  /* pointer to current permutation */

#ifdef DEBUG
  dstrc_enter("solid3_write_stress");
#endif

  /* bricks have 6 stress - use 3D matrix */
  if (gid->is_solid3_h8_222 || gid->is_solid3_h20_333)
  {
    INT i;
    INT el_type;
    INT stresstype;

    GiD_ResultType resulttype;

    CHAR* so3_stresstypenames[] = SOLID3_STRESSTYPE;
    TRANSLATION_TABLE so3_stresstype;

    INT components;

    CHAR* eqv_componentnames[] = { "Stress-equiv" };

    CHAR* xyz_componentnames[] = { "Stress-xx", "Stress-yy", "Stress-zz",
                                   "Stress-xy", "Stress-yz", "Stress-zx" };

    CHAR* rst_componentnames[] = { "Stress-rr", "Stress-ss", "Stress-tt",
                                   "Stress-rs", "Stress-st", "Stress-tr" };

    CHAR* i23_componentnames[] = { "Stress-11", "Stress-22", "Stress-33",
                                   "Angle-x1", "Angle-y1", "Angle-z1",
                                   "Angle-x2", "Angle-y2", "Angle-z2",
                                   "Angle-x3", "Angle-y3", "Angle-z3" };

    CHAR** componentnames;

    INT ngauss = 0;  /* number of Gauss points in element domain */

    init_translation_table(&so3_stresstype,
                           map_read_map(field->group, "so3_stresstypes"),
                           so3_stresstypenames);

    /* It's ridiculous! We store the flags to all elements and now we
     * only check the first one! */
    chunk_read_size_entry(&(field->ele_param), 0);

    el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];
    el_type = field->problem->element_type.table[el_type];
    if (el_type != el_solid3)
    {
      dserror("currently only all-solid3-meshes supported");
    }

    stresstype = field->ele_param.size_buf[solid3_variables.ep_size_stresstype];
    stresstype = so3_stresstype.table[stresstype];

    switch (stresstype)
    {
    case so3_stress_none:  /* default */
      break;
    case so3_stress_gpxyz:  /* gp: global system on element level */
      resulttype = GiD_Matrix;
      components = 6;
      componentnames = xyz_componentnames;
      break;
    case so3_stress_gprst:  /* gp: local  system on element level */
      resulttype = GiD_Matrix;
      components = 6;
      componentnames = rst_componentnames;
      break;
    case so3_stress_gp123:  /* gp: principal-stresses and directions */
      resulttype = GiD_MainMatrix;
      components = 12;
      componentnames = i23_componentnames;
      break;
    case so3_stress_ndxyz:  /* extrapolation to nodes: global system on element level    */
    case so3_stress_ndrst:  /* extrapolation to nodes: local  system on element level    */
    case so3_stress_nd123:  /* extrapolation to nodes: principal-stresses                */
    case so3_stress_ndeqv:  /* extrapolation to nodes: equivalent stress                 */
    default:
      dserror("solid3 stress %d needs to be done", stresstype);
    }

    /* So we assume all elements to be of the same type. */
    if (gid->is_solid3_h8_222)
    {
      GiD_BeginResult("solid3_stress", "ccarat", step, 
                      resulttype, GiD_OnGaussPoints,
                      gid->solid3_h8_222_name, NULL, 
                      components, componentnames);
      ngauss = 8;
      gperm = &(gperm_h_222[0]);
    }
    else /* (gid->is_solid3_h20_333) */
    {
      GiD_BeginResult("solid3_stress", "ccarat", step, 
                      resulttype, GiD_OnGaussPoints,
                      gid->solid3_h20_333_name, NULL, 
                      components, componentnames);
      ngauss = 27;
      gperm = &(gperm_h_333[0]);
    }

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT Id;
      DOUBLE* stress;
      INT dis;
      INT l;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_solid3) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      switch (stresstype)
      {
      case so3_stress_none:
        break;
      case so3_stress_gpxyz:  /* gp: global system on element level */

        /*
         * Take the right stress values from the chunk entry and write
         * them to gid.
         *
         * How are stresses ordered on chunk? 
         * Probably defined in out_pack_stress.
         * However, only stresses 'gpxyz' are packed there,
         * thus 'gprst', 'gp123', etc cannot be written */

        for (l=0; l<ngauss; ++l)
        {
          INT p;
          p = gperm[l];
          GiD_Write3DMatrix(Id+1,
                            stress[6*p+0], stress[6*p+1], stress[6*p+2],
                            stress[6*p+3], stress[6*p+4], stress[6*p+5]);
        }

        break;
      case so3_stress_gprst:  /* gp: local  system on element level */
        break;
      case so3_stress_gp123:  /* gp: principal-stresses */
        /* will need GiD_WriteMainMatrix */
        break;
      default:
        dserror("stress type %d not supported", stresstype);
      }
    }
    GiD_EndResult();

    destroy_translation_table(&so3_stresstype);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author bborn
  \date 01/07
*/
/*----------------------------------------------------------------------*/
void solid3_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("solid3_write_gauss");
#endif

  if (gid->is_solid3_h8_222)
  {
    GiD_BeginGaussPoint(gid->solid3_h8_222_name, GiD_Hexahedra, gid->solid3_h8_222_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_solid3_h20_333)
  {
    GiD_BeginGaussPoint(gid->solid3_h20_333_name, GiD_Hexahedra, gid->solid3_h20_333_name, 27, 0, 1);
    GiD_EndGaussPoint();
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author bborn
  \date 01/07
*/
/*----------------------------------------------------------------------*/
void solid3_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("solid3_write_mesh");
#endif

  if (gid->is_solid3_h8_222)
  {
    INT i;

    GiD_BeginMesh(gid->solid3_h8_222_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_solid3 || numnp != 8) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_solid3_h20_333)
  {
    INT i;

    /* we only do the first eight nodes */
    GiD_BeginMesh(gid->solid3_h20_333_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_solid3 || numnp != 20) continue;

      chunk_read_size_entry(&(field->mesh), i);

      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, 8);
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
