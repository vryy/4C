
#ifdef D_BRICK1

#include "post_brick1.h"
#include "../brick1/brick1.h"
#include "gid_out.h"
#include "post_gid.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Element specific service function.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void brick1_write_domain(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("brick1_write_domain");
#endif

  if (gid->is_brick1_222)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->brick1_222_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_brick1 || numnp != 8) continue;

      chunk_read_size_entry(chunk, i);

      for (j=0; j<8; j++)
        GiD_WriteScalar(Id+1, chunk->size_buf[0]);
    }
    GiD_EndResult();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_brick1_333)
  {
    GiD_BeginResult("Domains", "ccarat", 0, GiD_Scalar, GiD_OnGaussPoints,
                    gid->brick1_333_name, NULL, 0, NULL);

    for (i=0; i<field->numele; i++)
    {
      INT numnp;
      INT el_type;
      INT Id;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_brick1 || (numnp != 20 || numnp != 27)) continue;

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

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void brick1_write_stress(FIELD_DATA *field, GIDSET* gid, CHUNK_DATA* chunk, DOUBLE time, INT step)
{
#ifdef DEBUG
  dstrc_enter("brick1_write_stress");
#endif

  /* bricks have 6 stress - use 3D matrix */
  if (gid->is_brick1_222 || gid->is_brick1_333)
  {
    INT i;
    INT el_type;
    INT stresstype;

    CHAR* c1_stresstypenames[] = BRICK1_STRESSTYPE;
    TRANSLATION_TABLE c1_stresstype;

    CHAR* equi_componentnames[] = { "equivStress" };

    CHAR* xyz_componentnames[] = { "Stress-xx", "Stress-yy", "Stress-zz",
                                   "Stress-xy", "Stress-xz", "Stress-yz" };

    CHAR* rst_componentnames[] = { "Stress-rr", "Stress-ss", "Stress-tt",
                                   "Stress-rs", "Stress-st", "Stress-tr" };

    CHAR** componentnames;

    init_translation_table(&c1_stresstype,
                           map_read_map(field->group, "c1_stresstypes"),
                           c1_stresstypenames);

    /* It's ridiculous! We store the flags to all elements and now we
     * only check the first one! */
    chunk_read_size_entry(&(field->ele_param), 0);

    el_type = field->ele_param.size_buf[element_variables.ep_size_eltyp];
    el_type = field->problem->element_type.table[el_type];
    if (el_type != el_brick1)
    {
      dserror("currently only all-brick1-meshes supported");
    }

    stresstype = field->ele_param.size_buf[brick1_variables.ep_size_stresstyp];
    stresstype = c1_stresstype.table[stresstype];

    switch (stresstype)
    {
    case c1_nostr:              /* default                           */
      dserror("stress output without stresses. odd.");
      break;
    case c1_gpxyz:              /* gp: global system on element level    */
      componentnames = xyz_componentnames;
      break;
    case c1_gprst:              /* gp: local  system on element level    */
      componentnames = equi_componentnames;
      break;
    case c1_gp123:              /* gp: principal-stresses                */
      componentnames = rst_componentnames;
      break;
    case c1_npxyz:              /* extrapolation to nodes: global system on element level    */
    case c1_nprst:              /* extrapolation to nodes: local  system on element level    */
    case c1_np123:              /* extrapolation to nodes: principal-stresses                */
    case c1_npeqs:              /* extrapolation to nodes: equivalent stress                 */
    default:
      dserror("brick1 stress needs to be done");
    }

    /* So we assume all elements to be of the same type. */
    if (gid->is_brick1_222)
    {
      GiD_BeginResult("brick1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                      gid->brick1_222_name, NULL, 6, componentnames);
    }
    else /* (gid->is_brick1_333) */
    {
      GiD_BeginResult("brick1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                      gid->brick1_333_name, NULL, 6, componentnames);
    }

    for (i=0; i<field->numele; ++i)
    {
      INT numnp;
      INT Id;
      DOUBLE* stress;
      INT dis;

      /* read the element's data */
      get_element_params(field, i, &Id, &el_type, &dis, &numnp);

      if (el_type != el_brick1) continue;

      chunk_read_value_entry(chunk, i);
      stress = chunk->value_buf;

      switch (stresstype)
      {
      case c1_gpxyz:              /* gp: global system on element level    */

        /*
         * Take the right stress values from the chunk entry and write
         * them to gid.
         *
         * Yet to be done. */

        /*
        for (l=0; l<ngauss; ++l)
        {
          INT p;
          p = gaussperm9[l];
          GiD_Write3DMatrix(Id+1,
                            stress[9*p+0], stress[9*p+1], stress[9*p+2],
                            stress[9*p+3], stress[9*p+7], stress[9*p+8]);
        }
        */
        break;
      case c1_gprst:              /* gp: local  system on element level    */
        break;
      case c1_gp123:              /* gp: principal-stresses                */
        break;
      default:
        dserror("stress type %d not supported", stresstype);
      }
    }
    GiD_EndResult();

    destroy_translation_table(&c1_stresstype);
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
void brick1_write_gauss(GIDSET* gid)
{
#ifdef DEBUG
  dstrc_enter("brick1_write_gauss");
#endif

  if (gid->is_brick1_222)
  {
    GiD_BeginGaussPoint(gid->brick1_222_name, GiD_Hexahedra, gid->brick1_222_name, 8, 0, 1);
    GiD_EndGaussPoint();
  }
  if (gid->is_brick1_333)
  {
    GiD_BeginGaussPoint(gid->brick1_333_name, GiD_Hexahedra, gid->brick1_333_name, 27, 0, 1);
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
void brick1_write_mesh(FIELD_DATA *field, GIDSET* gid, INT* first_mesh)
{
#ifdef DEBUG
  dstrc_enter("brick1_write_mesh");
#endif

  if (gid->is_brick1_222)
  {
    INT i;

    GiD_BeginMesh(gid->brick1_222_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_brick1 || numnp != 8) continue;

      chunk_read_size_entry(&(field->mesh), i);
      get_gid_node_ids(field, field->mesh.size_buf, mesh_entry, numnp);
      GiD_WriteElement(Id+1, mesh_entry);
    }

    GiD_EndElements();
    GiD_EndMesh();
  }

  /*--------------------------------------------------------------------*/

  if (gid->is_brick1_333)
  {
    INT i;

    /* we only do the first eight nodes */
    GiD_BeginMesh(gid->brick1_333_name, 3, GiD_Hexahedra, 8);

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

      if (el_type != el_brick1 || numnp != 20) continue;

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
