/*!
\file
\brief Postprocessing utility that takes ccarat output and produces a
GiD input file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Filters like this one are special inhabitants of the ccarat
world. They are always single processor applications yet they share
some code with ccarat and are closely linked to ccarat internals.

The general idea is that we cannot load the whole result data into
memory at once.

Filters are independent programs, thus they have their own main
function. This filter's main function is the last function in this
file (to keep the number of needed function prototypes as small as
possible). You might want to start reading the file from there.

\author u.kue
\date 09/04

*/

#include <assert.h>

#include "gid_out.h"

#include "post_gid.h"

/* map from DIS_TYP to GiD_ElementType */
static GiD_ElementType gid_type_map[] = {
  GiD_NoElement,                /* dis_none,       -- unknown dis type */
  GiD_Quadrilateral,            /* quad4,          -- 4 noded quadrilateral */
  GiD_Quadrilateral,            /* quad8,          -- 8 noded quadrilateral */
  GiD_Quadrilateral,            /* quad9,          -- 9 noded quadrilateral */
  GiD_Triangle,                 /* tri3,           -- 3 noded triangle */
  GiD_Triangle,                 /* tri6,           -- 6 noded triangle */
  GiD_Hexahedra,                /* hex8,           -- 8 noded hexahedra */
  GiD_Hexahedra,                /* hex20,          -- 20 noded hexahedra */
  GiD_Hexahedra,                /* hex27,          -- 27 noded hexahedra */
  GiD_Tetrahedra,               /* tet4,           -- 4 noded tetrahedra */
  GiD_Tetrahedra,               /* tet10,          -- 4 noded tetrahedra */
  GiD_Linear,                   /* line2,          -- 2 noded line */
  GiD_Linear                    /* line3           -- 3 noded line */
};

/*----------------------------------------------------------------------*/
/*!
  \brief Write coordinates.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_coords(PROBLEM_DATA* problem, INT num)
{
  INT i;
  FIELD_DATA* field;
  MAP* coords_group;
  INT value_entry_length;
  INT value_offset;
  DOUBLE x[3];

  field = &(problem->discr[num]);

  /* figure out the elements used in this discretization */
  fprintf(allfiles.out_err, "%s: Write node coordinates\n", field->name);
  coords_group = map_read_map(field->table, "coords");

  value_entry_length = map_read_int(coords_group, "value_entry_length");
  value_offset = map_read_int(coords_group, "value_offset");

  dsassert(value_entry_length == problem->ndim, "dimension mismatch");

  GiD_BeginCoordinates();

  /* This is a shortcut. We'll always write three coordinates. That's
   * what the gid library will do anyway. So in case we are reading a
   * two dimensional problem just set the last component to zero. */
  x[2] = 0;

  fseek(field->value_file, value_offset, SEEK_SET);
  for (i=0; i<field->numnp; ++i) {

    if (fread(x, sizeof(DOUBLE), value_entry_length, field->value_file)!=value_entry_length) {
      dserror("reading value file of discretization %s failed", field->name);
    }

    /* All GiD ids are one based. */
    GiD_WriteCoordinates(field->node_ids[i]+1, x[0], x[1], x[2]);
  }

  GiD_EndCoordinates();
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write elements.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_elements(PROBLEM_DATA* problem, INT num, INT element, INT nodes)
{
  INT i;
  FIELD_DATA* field;
  MAP* mesh_group;
  INT size_entry_length;
  INT size_offset;
  INT* mesh_entry;

  field = &(problem->discr[num]);

  fprintf(allfiles.out_err, "%s: Write element %s with %d nodes\n",
         field->name, element_info[element].name, nodes);

  mesh_group = map_read_map(field->table, "mesh");
  size_entry_length = map_read_int(mesh_group, "size_entry_length");
  size_offset = map_read_int(mesh_group, "size_offset");

  dsassert(size_entry_length >= nodes+3, "node count mismatch");

  mesh_entry = (INT*)CCACALLOC(size_entry_length, sizeof(INT));

  GiD_BeginElements();

  fseek(field->size_file, size_offset, SEEK_SET);
  for (i=0; i<field->numele; ++i) {
    INT major;
    INT minor;

    if (fread(mesh_entry, sizeof(INT), size_entry_length, field->size_file)!=size_entry_length) {
      dserror("reading size file of discretization %s failed", field->name);
    }

    major = mesh_entry[1];
    minor = mesh_entry[2];

    /* convert from file major/minor to internal major/minor */
    minor = field->internal_minors[major][minor];
    major = field->internal_majors[major];

    /* Write the element if major number and numbers of nodes match. */
    if ((element == major) &&
        (element_info[element].variant[minor].node_number == nodes)) {
      INT j;
      /* All GiD ids are one based. */
      /* Additionally we have to convert discretization local ids to
       * global ids. */
      for (j=0; j<nodes; ++j) {
        mesh_entry[j+3] = field->node_ids[mesh_entry[j+3]] + 1;
      }
      GiD_WriteElement(field->element_type[i].Id+1, mesh_entry+3);
    }
  }

  GiD_EndElements();

  CCAFREE(mesh_entry);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the element domains.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_domain(FIELD_DATA *field, MAP* group)
{
  INT j;
  INT k;
  GiD_ElementType EType;
  INT NNode;
  INT el;

  fprintf(allfiles.out_err, "%s: Write domain\n", field->name);

  dsassert(map_read_int(group, "size_entry_length") == 1, "domain size mismatch");

  /* Find all types of elements in this field. */
  for (el=0; el<el_count; ++el) {
    ELEMENT_INFO* info;
    info = &(element_info[el]);

    /* Check the minor numbers. */
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (field->element_flags[el][j]) {
        INT GP_number;
        CHAR GP_name[30];

        NNode = info->variant[j].node_number;
        EType = gid_type_map[info->variant[j].dis_type];
        GP_number = info->variant[j].gauss_number;

        sprintf(GP_name, "%s_%d", info->name, GP_number);
        GiD_BeginResult("Domain", "ccarat", 0.0, GiD_Scalar, GiD_OnGaussPoints,
                        GP_name, NULL, 0, NULL);

        fseek(field->size_file, map_read_int(group, "size_offset"), SEEK_SET);

        /* Iterate all elements and write those with the current type. */
        for (k = 0; k < field->numele; ++k) {
          INT domain;

          /* Read the data even if this element has the wrong
           * type. We rely on the file pointer beeing set. */
          if (fread(&domain, sizeof(INT), 1, field->size_file) != 1) {
            dserror("failed to read domain number of element %d", k);
          }

          if ((field->element_type[k].major == el) &&
              (field->element_type[k].minor == j)) {
            INT i;

            for (i = 0; i < GP_number; ++i) {
              GiD_WriteScalar(field->element_type[k].Id+1, domain);
            }
          }
        }

        GiD_EndResult();
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal displacements for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_displacement(FIELD_DATA *field, MAP* group)
{
  INT k;
  MAP* disp_group;
  INT value_entry_length;
  DOUBLE x[3];
  CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };
  DOUBLE time;
  INT step;
  CHAR buf[100];

  disp_group = map_read_map(group, "displacement");

  dsassert(map_read_int(disp_group, "size_entry_length") == 0,
           "displacement size mismatch");

  /* This is supposed to equal the number of dimensions. */
  value_entry_length = map_read_int(disp_group, "value_entry_length");
  dsassert((value_entry_length==2) || (value_entry_length==3),
           "displacement vector length corrupt");

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  fprintf(allfiles.out_err, "%s: Write displacement of step %d\n", field->name, step);

  sprintf(buf, "%s_displacement", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Vector, GiD_OnNodes,
                  NULL, NULL, value_entry_length, componentnames);

  fseek(field->value_file, map_read_int(disp_group, "value_offset"), SEEK_SET);

  /* In case this is a 2d problem. */
  x[2] = 0;

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k) {

    if (fread(x, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
      dserror("failed to read displacement of node %d", k);
    }

    GiD_WriteVector(field->node_ids[k]+1, x[0], x[1], x[2]);
  }

  GiD_EndResult();
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal velocities for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_velocity(FIELD_DATA *field, MAP* group)
{
  INT k;
  MAP* velocity_group;
  INT value_entry_length;
  DOUBLE velocity[3];
  CHAR* componentnames[] = { "x-vel", "y-vel", "z-vel" };
  DOUBLE time;
  INT step;
  CHAR buf[100];

  velocity_group = map_read_map(group, "velocity");

  dsassert(map_read_int(velocity_group, "size_entry_length") == 0,
           "velocity size mismatch");

  /* This is supposed to equal the number of dimensions. */
  value_entry_length = map_read_int(velocity_group, "value_entry_length");
  dsassert((value_entry_length==2) || (value_entry_length==3),
           "velocity vector length corrupt");

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  fprintf(allfiles.out_err, "%s: Write velocity of step %d\n", field->name, step);

  sprintf(buf, "%s_velocity", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Vector, GiD_OnNodes,
                  NULL, NULL, value_entry_length, componentnames);

  fseek(field->value_file, map_read_int(velocity_group, "value_offset"), SEEK_SET);

  /* In case this is a 2d problem. */
  velocity[2] = 0;

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k) {

    if (fread(velocity, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
      dserror("failed to read velocity of node %d", k);
    }

    GiD_WriteVector(field->node_ids[k]+1, velocity[0], velocity[1], velocity[2]);
  }

  GiD_EndResult();
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the nodal pressure for one time step.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_pressure(FIELD_DATA *field, MAP* group)
{
  INT k;
  MAP* pressure_group;
  INT value_entry_length;
  DOUBLE pressure;
  CHAR* componentnames[] = { "pressure" };
  DOUBLE time;
  INT step;
  CHAR buf[100];

  pressure_group = map_read_map(group, "pressure");

  dsassert(map_read_int(pressure_group, "size_entry_length") == 0,
           "pressure size mismatch");

  value_entry_length = map_read_int(pressure_group, "value_entry_length");
  dsassert(value_entry_length==1, "pressure item length corrupt");

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  fprintf(allfiles.out_err, "%s: Write pressure of step %d\n", field->name, step);

  sprintf(buf, "%s_pressure", fieldnames[field->type]);
  GiD_BeginResult(buf, "ccarat", step, GiD_Scalar, GiD_OnNodes,
                  NULL, NULL, value_entry_length, componentnames);

  fseek(field->value_file, map_read_int(pressure_group, "value_offset"), SEEK_SET);

  /* Iterate all nodes. */
  for (k = 0; k < field->numnp; ++k) {

    if (fread(&pressure, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
      dserror("failed to read pressure of node %d", k);
    }

    GiD_WriteScalar(field->node_ids[k]+1, pressure);
  }

  GiD_EndResult();
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the gauss point stresses for one time step.

  This is element specific again.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_stress(FIELD_DATA *field, MAP* group)
{
  INT i;
  INT j;
  INT k;
  MAP* stress_group;
  INT value_entry_length;
  DOUBLE* stress;
  DOUBLE time;
  INT step;

/*
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/

  INT           gaussperm4[4] = {3,1,0,2};
  /*INT           gaussperm8[8] = {0,4,2,6,1,5,3,7};*/
  INT           gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
  /*INT           gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};*/

  stress_group = map_read_map(group, "stress");

  dsassert(map_read_int(stress_group, "size_entry_length") == 0,
           "stress size mismatch");

  value_entry_length = map_read_int(stress_group, "value_entry_length");
  stress = (DOUBLE*)CCACALLOC(value_entry_length, sizeof(DOUBLE));

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  fprintf(allfiles.out_err, "%s: Write stress of step %d\n", field->name, step);

  for (i=0; i<el_count; ++i) {
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (field->element_flags[i][j]) {

        ELEMENT_INFO* info;
        INT GP_number;
        CHAR GP_name[30];
        GiD_ElementType EType;
        INT NNode;

        info = &(element_info[i]);

        NNode = info->variant[j].node_number;
        EType = gid_type_map[info->variant[j].dis_type];
        GP_number = info->variant[j].gauss_number;
        sprintf(GP_name, "%s_%d", info->name, GP_number);

        dsassert(value_entry_length >= element_info[i].variant[j].stress_matrix_size,
                 "stress too short for element type");

        switch (i) {
#ifdef D_WALL1
        case el_wall1: {        /* 2D plane stress - plane strain element */
          CHAR* componentnames[] = { "Stress-xx", "Stress-yy", "Stress-xy",
                                     "Stress-zz", "damage", "||u||" };
          GiD_BeginResult("wall1_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                          NULL, NULL, 6, componentnames);
          fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);

          switch (j) {
#if 0                           /* not yet... */
          case MINOR_WALL1_11: { /* 3-noded wall1 1x1 GP */

            /* Iterate all elements, choose the ones with matching type */
            for (k=0; k<field->numele; ++k) {
              if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
                dserror("failed to read stress of element %d", k);
              }
              if ((field->element_type[k].major == i) &&
                  (field->element_type[k].minor == j)) {
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[0], stress[1], stress[2],
                                  stress[3], 0, 0);
              }
            }
            break;
          }
#endif
          case MINOR_WALL1_22:  /* 4-noded wall1 2x2 GP */

#if 0                           /* not yet... */
          case MINOR_WALL1_8_33: /* 8-noded wall1 3x3 GP */
          case MINOR_WALL1_9_33: /* 9-noded wall1 3x3 GP */
#endif

            /* Iterate all elements, choose the ones with matching type */
            for (k=0; k<field->numele; ++k) {
              if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
                dserror("failed to read stress of element %d", k);
              }
              if ((field->element_type[k].major == i) &&
                  (field->element_type[k].minor == j)) {
                INT l;
                for (l=0; l<element_info[i].variant[j].gauss_number; ++l) {
                  INT p;
                  p = gaussperm4[l];
                  GiD_Write3DMatrix(field->node_ids[k]+1,
                                    stress[9*p+0], stress[9*p+1], stress[9*p+2],
                                    stress[9*p+3], stress[9*p+7], stress[9*p+8]);
                }
              }
            }
            break;
          default:
            dserror("unknown minor type %d for wall element", j);
          }

          GiD_EndResult();
          break;
        }
#endif
#ifdef D_BEAM3
        case el_beam3: {        /* structural 3D-beam element */
          CHAR* componentnames[] = { "N-x", "V-y", "V-z", "M-x", "M-y", "M-z" };
          GiD_BeginResult("beam3_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                          NULL, NULL, 6, componentnames);
          fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);

          /* Iterate all elements, choose the ones with matching type */
          for (k=0; k<field->numele; ++k) {
            if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
              dserror("failed to read stress of element %d", k);
            }
            if ((field->element_type[k].major == i) &&
                (field->element_type[k].minor == j)) {
              INT l;
              for (l=0; l<element_info[i].variant[j].gauss_number; ++l) {
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[6*l+0], stress[6*l+1], stress[6*l+2],
                                  stress[6*l+3], stress[6*l+4], stress[6*l+5]);
              }
            }
          }

          GiD_EndResult();
          break;
        }
#endif
#ifdef D_INTERF
        case el_interf: {       /* 1D interface element (combination only with wall) */
          CHAR* tn_componentnames[] = { "stress-tang","stress-normal","dummy","D-nomal","D-tang","dummy" };
          CHAR* xy_componentnames[] = { "stress-sxx","stress-syy","stress-sxy","D-nomal","D-tang","dummy" };
          CHAR** componentnames;

          if (map_has_string(stress_group, "interf_orient", "global")) {
            componentnames = xy_componentnames;
          }
          else {
            componentnames = tn_componentnames;
          }

          GiD_BeginResult("interface_stresses", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                          NULL, NULL, 6, componentnames);
          fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);

          /* Iterate all elements, choose the ones with matching type */
          for (k=0; k<field->numele; ++k) {
            if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
              dserror("failed to read stress of element %d", k);
            }
            if ((field->element_type[k].major == i) &&
                (field->element_type[k].minor == j)) {
              INT l;
              /*
               * This is fake. The interface element pretends to be of
               * higher dimension. */
              if (j == MINOR_INTERF_22) {
                l = 0;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 1;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 0;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
              }
              else if (j == MINOR_INTERF_33) {
                l = 0;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 2;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 0;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 1;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 2;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 1;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 0;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
                l = 1;
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[5*l+0], stress[5*l+1], stress[5*l+2],
                                  stress[5*l+3], stress[5*l+4], 0);
              }
              else {
                dserror("ups. impossible: j=%d", j);
              }
            }
          }

          GiD_EndResult();
          break;
        }
#endif
#ifdef D_WALLGE
        case el_wallge: {       /* gradient enhanced wall element */
          CHAR* componentnames[] = { "Stress-xx","Stress-yy","Stress-xy","damage","loc_aequiv_strain","nonloc_aequiv_strain" };
          GiD_BeginResult("wallge_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                          NULL, NULL, 6, componentnames);
          fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);

          /* Iterate all elements, choose the ones with matching type */
          for (k=0; k<field->numele; ++k) {
            if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
              dserror("failed to read stress of element %d", k);
            }
            if ((field->element_type[k].major == i) &&
                (field->element_type[k].minor == j)) {
              INT l;
              for (l=0; l<element_info[i].variant[j].gauss_number; ++l) {
                INT p;
                p = gaussperm4[l];
                GiD_Write3DMatrix(field->node_ids[k]+1,
                                  stress[6*p+0], stress[6*p+1], stress[6*p+2],
                                  stress[6*p+3], stress[6*p+4], stress[6*p+5]);
              }
            }
          }

          GiD_EndResult();
          break;
        }
#endif
#ifdef D_AXISHELL
        case el_axishell: {     /* 1D axisymmetrical shell element */
          CHAR* componentnames[] = { "n_s","n_theta","m_s","m_theta","q_s","unused" };
          GiD_BeginResult("axishell_forces", "ccarat", step, GiD_Matrix, GiD_OnGaussPoints,
                          NULL, NULL, 6, componentnames);
          fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);

          /* Iterate all elements, choose the ones with matching type */
          for (k=0; k<field->numele; ++k) {
            if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
              dserror("failed to read stress of element %d", k);
            }
            if ((field->element_type[k].major == i) &&
                (field->element_type[k].minor == j)) {
              GiD_Write3DMatrix(field->node_ids[k]+1,
                                stress[0], stress[1], stress[2],
                                stress[3], stress[4], 0);
            }
          }

          GiD_EndResult();
          break;
        }
#endif
        default:
          dserror("stress output not supported for element type %d", i);
        }
      }
    }
  }

  CCAFREE(stress);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write connectivity and node coordinates.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_mesh(PROBLEM_DATA *problem, INT num)
{
  INT i, j, k;
  GiD_ElementType EType;
  INT NNode;
  CHAR mesh_name[100];
  FIELD_DATA* field;
  INT first_mesh = 1;

  field = &(problem->discr[num]);

  /* There are different meshes for different element types. However,
   * the number of gauss points is of no importance here. So elements
   * with different minor numbers but the same number of nodes live in
   * one mesh. */

  for (i=0; i<el_count; ++i) {
    ELEMENT_INFO* info;

    /* We keep all node numbers that we've handled in an array. */
    INT node_numbers[MAX_EL_MINOR];
    GiD_ElementType gid_type[MAX_EL_MINOR];
    for (j=0; j<MAX_EL_MINOR; ++j) {
      node_numbers[j] = 0;
      gid_type[j] = GiD_NoElement;
    }

    info = &(element_info[i]);

    /* Check the minor numbers. */
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (field->element_flags[i][j]) {
        NNode = info->variant[j].node_number;
        EType = gid_type_map[info->variant[j].dis_type];

        dsassert((NNode > 0) && (EType != GiD_NoElement),
                 "non-existing element variant selected");

        for (k=0; k<MAX_EL_MINOR; ++k) {
          if (node_numbers[k] == NNode) {
            /* This type of element with this node number has already
             * been done. */
            dsassert(gid_type[k] == EType,
                     "node numbers match but it's not the same type of element");
            break;
          }
          if (node_numbers[k] == 0) {
            /* Store all these elements */

            /* encode all information in the mesh name */
            sprintf(mesh_name, "%s_%s_%d", field->name, info->name, NNode);

            /* do the writing */
            GiD_BeginMesh(mesh_name, problem->ndim, EType, NNode);
            if (first_mesh) {
              first_mesh = 0;
              write_coords(problem, num);
            }
            write_elements(problem, num, i, NNode);
            GiD_EndMesh();

            /* mark node number as done */
            node_numbers[k] = NNode;
            gid_type[k] = EType;
            break;
          }
        }
      }
    }
  }

  /*--------------------------------------------------------------------*/
  /* Tell GiD about our element's gauss points. */

  for (i=0; i<el_count; ++i) {
    ELEMENT_INFO* info;
    info = &(element_info[i]);

    /* Check the minor numbers. */
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (field->element_flags[i][j]) {
        INT GP_number;
        CHAR GP_name[30];

        NNode = info->variant[j].node_number;
        EType = gid_type_map[info->variant[j].dis_type];
        GP_number = info->variant[j].gauss_number;

        /* The gauss point definition is meant just for this
         * discretizations meshes.  */
        sprintf(mesh_name, "%s_%s_%d", field->name, info->name, NNode);

        sprintf(GP_name, "%s_%d", info->name, GP_number);

        GiD_BeginGaussPoint(GP_name, EType, mesh_name, GP_number, 0, 1);
        GiD_EndGaussPoint();
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a whole field.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_field(PROBLEM_DATA* problem, INT num)
{
  INT i;
  FIELD_DATA* field;

  field = &(problem->discr[num]);

  /*--------------------------------------------------------------------*/
  /* Write connectivity and node coordinates. */

  write_mesh(problem, num);

  /*--------------------------------------------------------------------*/
  /* Now we could define a few range tables. Maybe it's nice to have
   * per discretization range tables. But they are highly specific. */
#if 0
  GiD_BeginRangeTable( char * name );
  GiD_WriteMinRange( double max, char * name );
  GiD_WriteRange( double min, double max, char * name );
  GiD_WriteMaxRange( double min, char * name );
  GiD_EndRangeTable();
#endif

  /*--------------------------------------------------------------------*/
  /* Optionally there is one domain table per field. It's the first result. */

  if (map_has_map(field->table, "domain")) {
    write_domain(field, map_read_map(field->table, "domain"));
  }

  /*--------------------------------------------------------------------*/
  /* Now it's time to write the time dependend results. */

  for (i=0; i<problem->num_results; ++i) {
    MAP* result_group;
    result_group = problem->result_group[i];

    /* We iterate the list of all results. Here we are interested in
     * the results of this discretization. */
    if (match_field_result(field, result_group)) {

      /* nodal result are written independent of the elements in this
       * discretization */
      if (map_has_map(result_group, "displacement")) {
        write_displacement(field, result_group);
      }
      if (map_has_map(result_group, "velocity")) {
        write_velocity(field, result_group);
      }
      if (map_has_map(result_group, "pressure")) {
        write_pressure(field, result_group);
      }

      /* element dependend results */
      if (map_has_map(result_group, "stress")) {
        write_stress(field, result_group);
      }
    }
  }
}


#ifdef D_FSI

/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing all two fsi fields.

  The ale field is not written on it's own. Instead the ale
  displacements are written to the fluid nodes. The structure field is
  handled by standard methods, but only if it's there. There are
  fluid-ale problems, too.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
static void write_fsi(PROBLEM_DATA* problem)
{
  INT i;
  FIELD_DATA* struct_field = NULL;
  FIELD_DATA* fluid_field;
  FIELD_DATA* ale_field;

  INT* fluid_ale_connect;

  /* Find the corresponding discretizations. We don't rely on any order. */
  for (i=0; i<problem->num_discr; ++i) {
    if (problem->discr[i].type == structure) {
      struct_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == fluid) {
      fluid_field = &(problem->discr[i]);
    }
    else if (problem->discr[i].type == ale) {
      ale_field = &(problem->discr[i]);
    }
    else {
      dserror("unknown field type %d", problem->discr[i].type);
    }
  }

  /*--------------------------------------------------------------------*/
  /* Write connectivity and node coordinates. */

  write_mesh(problem, fluid_field - problem->discr);
  if (struct_field != NULL) {
    write_mesh(problem, struct_field - problem->discr);
  }

  /*--------------------------------------------------------------------*/
  /* Optionally there is one domain table per field. It's the first result. */

#if 0
  if (map_has_map(fluid_field->table, "domain")) {
    write_domain(fluid_field, map_read_map(fluid_field->table, "domain"));
  }
  if (struct_field != NULL) {
    if (map_has_map(struct_field->table, "domain")) {
      write_domain(struct_field, map_read_map(struct_field->table, "domain"));
    }
  }
#endif

  /*--------------------------------------------------------------------*/
  /* We need to find the connection between ale and fluid nodes. */

  {
    INT *fluid_struct_connect;
    post_find_fsi_coupling(problem,
                           struct_field, fluid_field, ale_field,
                           &fluid_struct_connect, &fluid_ale_connect);

    /*
     * We don't need the connection to the structure field here. Free
     * it immediately. */
    if (struct_field != NULL) {
      CCAFREE(fluid_struct_connect);
    }
  }

  /*--------------------------------------------------------------------*/
  /* Now it's time to write the time dependend results. */

  /* We iterate the list of all results. */
  for (i=0; i<problem->num_results; ++i) {
    MAP* result_group;
    result_group = problem->result_group[i];

    /* The ale field.
     *
     * We are interested in the displacements. But we'll write them
     * to the fluid nodes. */
    if (match_field_result(ale_field, result_group)) {

      if (map_has_map(result_group, "displacement")) {
        INT k;
        MAP* disp_group;
        INT value_entry_length;
        INT value_offset;
        DOUBLE x[3];
        CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };
        DOUBLE time;
        INT step;

        disp_group = map_read_map(result_group, "displacement");

        dsassert(map_read_int(disp_group, "size_entry_length") == 0,
                 "displacement size mismatch");

        value_offset = map_read_int(disp_group, "value_offset");
        value_entry_length = map_read_int(disp_group, "value_entry_length");
        dsassert(value_entry_length==problem->ndim,
                 "displacement vector length corrupt");

        time = map_read_real(result_group, "time");
        step = map_read_int(result_group, "step");

        fprintf(allfiles.out_err, "%s: Write displacement of step %d\n", ale_field->name, step);

        GiD_BeginResult("fluid_displacement", "ccarat", step, GiD_Vector, GiD_OnNodes,
                        NULL, NULL, value_entry_length, componentnames);

        /* In case this is a 2d problem. */
        x[2] = 0;

        /* Iterate all fluid nodes. */
        for (k = 0; k < fluid_field->numnp; ++k) {

          if (fluid_ale_connect[k] != -1) {

            /* We have to seek each time. This is not efficient. */
            fseek(ale_field->value_file,
                  value_offset + sizeof(DOUBLE)*fluid_ale_connect[k]*value_entry_length,
                  SEEK_SET);
            if (fread(x, sizeof(DOUBLE), value_entry_length, ale_field->value_file) != value_entry_length) {
              dserror("failed to read displacement of node %d", k);
            }

            GiD_WriteVector(fluid_field->node_ids[k]+1, x[0], x[1], x[2]);
          }
          else {
            GiD_WriteVector(fluid_field->node_ids[k]+1, 0, 0, 0);
          }
        }

        GiD_EndResult();
      }
    }

    /* The fluid field */
    if (match_field_result(fluid_field, result_group)) {

      if (map_has_map(result_group, "velocity")) {
        write_velocity(fluid_field, result_group);
      }
      if (map_has_map(result_group, "pressure")) {
        write_pressure(fluid_field, result_group);
      }

#if 0
      /* No fluid stresses yet? */

      /* element dependend results */
      if (map_has_map(result_group, "stress")) {
        write_stress(fluid_field, result_group);
      }
#endif
    }

    /* The structure field */
    if (struct_field != NULL) {
      if (match_field_result(struct_field, result_group)) {

        if (map_has_map(result_group, "displacement")) {
          write_displacement(struct_field, result_group);
        }
        if (map_has_map(result_group, "velocity")) {
          write_velocity(struct_field, result_group);
        }
        if (map_has_map(result_group, "pressure")) {
          write_pressure(struct_field, result_group);
        }

        if (map_has_map(result_group, "stress")) {
          write_stress(struct_field, result_group);
        }
      }
    }
  }

  CCAFREE(fluid_ale_connect);
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief The filter's main function.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
int main(int argc, char** argv)
{
  CHAR basename[100];
  CHAR filename[100];
  MAP control_table;
  PROBLEM_DATA problem;
  INT i;

  if (argc != 2) {
    printf("usage: %s control-file\n", argv[0]);
    return 1;
  }

  setup_filter(argv[1], &control_table, basename);

  dsassert(map_has_string(&control_table, "version", "0.1"),
           "expect version 0.1 control file");

  /* Debug output */
  /*map_print(stdout, &control_table, 0);*/

  init_problem_data(&problem, &control_table);

  sprintf(filename, "%s.flavia.res", basename);
  GiD_OpenPostResultFile(filename);

  /* Some problem types have to be handled in a special way. The
   * default is to simply write the results as they are. */
  switch (problem.type) {
#ifdef D_FSI
  case prb_fsi:
    /*
     * We know there is an ale, a fluid and a structure field. We need
     * to find the fluid nodes that match the ale nodes to be able to
     * output fluid node displacements. */
    dsassert(problem.num_discr==3, "three fields expected for fsi");
    write_fsi(&problem);
    break;
  case prb_fluid:
    if (problem.num_discr==2) {
      /*
       * So we have two discretizations. If these are two fields we
       * have a fluid-ale problem. Otherwise it might be the infamous
       * projection method or something. */
      if (map_symbol_count(&control_table, "field") == 2) {
        write_fsi(&problem);
        break;
      }
    }
#endif
    /* No break here! The simple fluid case is done by the default
     * code. */
  default:
    for (i=0; i<problem.num_discr; ++i) {
#ifdef D_SHELL9
      if (problem.discr[i].is_shell9_problem) {
        write_shell9_field(&problem, i);
      }
      else
#endif
        write_field(&problem, i);
    }
  }

  GiD_ClosePostResultFile();

  fprintf(allfiles.out_err, "Done.\n");
  return 0;
}
