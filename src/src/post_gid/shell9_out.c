/*!
\file
\brief The special case of a shell9 problem.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Discretizations that consist of shell9 elements (if there are shell9
elements no other elements are allowed, there must be just one shell9
variant, too.) are special. The ccarat output contains node
information for every layer of the shell. It's up to the filter to
produce meaningful results that can be visualized.

ccarat used to contain two output versions for shell9 elements. In
both of them brick elements where used. The first version made up a
mesh of bricks and gave averaged values at the nodes. The second did
not connect the brick to each other and was thus able to output the
results as calculated by the shell9 element. Both versions are
reimplemented here.

The stresses are extrapolated to the nodes because GiD had
difficulties to visualize gauss point stresses for bricks.

\author u.kue
\date 09/04

*/

#ifdef D_SHELL9

#include "post_gid.h"

#include "gid_out.h"

#include "../shell9/shell9.h"


extern CHAR* fieldnames[];


#if 0
/* In order to do some postprocessing calculations (extrapolating the
 * stresses to the nodes) we need to know the mesh. Thus we store all
 * information that is gathered while writing the mesh to GiD. Please
 * note that in general we cannot assume our meshes to be small
 * enougth be fit in one processors memory. That's the reason we care
 * so much about binary output. But shell9 meshes are special
 * cases. There are not millions of elements, and if there are one
 * day, we'll be able to abandon this variable. */
static DOUBLE* coords;
#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a shell9 discretization the
  smoothed way.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_shell9_field_smoothed(PROBLEM_DATA* problem, INT num)
{
  FIELD_DATA* field;

  field = &(problem->discr[num]);

  dserror("write_shell9_field_smoothed not implemented yet");
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write a mesh that consists of brick elements.

  Danger! If there are other discretizations this code will fail!
  That's because we don't know which node and element ids are assigned
  in those discretizations. In the normal case the value
  field->first_node_id tells us where our range starts. But here nodes
  and elements are generated that ccarat never knew about.

  This problem could be solved inside the filter. For now it's not
  worth the trouble.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_shell9_mesh_unsmoothed(PROBLEM_DATA* problem, INT num)
{
  FIELD_DATA* field;
  CHAR mesh_name[100];
  GiD_ElementType EType;
  INT NNode;
  INT value_entry_length;
  INT value_offset;
  INT el;
  INT nodes_per_layer;
  INT nodes_per_element;
  INT numnp;
  DOUBLE* x;
  INT GP_number;
  CHAR GP_name[30];

  MAP* coords_group;

  field = &(problem->discr[num]);

  /*--------------------------------------------------------------------*/
  /* Write connectivity and node coordinates. */

  switch (field->s9_minor) {
  case MINOR_SHELL9_4_22:       /* 4-noded shell9 2x2 GP */
  case MINOR_SHELL9_4_33:       /* 4-noded shell9 3x3 GP */

    /*
     * Minor number specific values. These are not the values from the
     * element_info table. These values do not describe the real
     * appearance of the shell9 element. Instead they describe the
     * brick elements that are generated, but it still helps to have
     * independent values. */
    NNode = 8;
    EType = GiD_Hexahedra;

    /* On the other hand we know the real number of nodes of this
     * variant. */
    numnp = 4;

    sprintf(mesh_name, "%s_%s_%d", field->name, "shell9", NNode);

    GiD_BeginMesh(mesh_name, 3, EType, NNode);
    GiD_BeginCoordinates();

    printf("%s: Write node coordinates\n", field->name);
    coords_group = map_read_map(field->table, "shell9_coords");

    value_entry_length = map_read_int(coords_group, "value_entry_length");
    value_offset = map_read_int(coords_group, "value_offset");

    dsassert(value_entry_length == 3*numnp*(field->s9_layers+1),
             "invalid shell9 coords entry length");
    x = (DOUBLE*)CCACALLOC(value_entry_length, sizeof(DOUBLE));

    fseek(field->value_file, value_offset, SEEK_SET);

    /* We number the nodes per element. Inside an element the nodes of
     * each layer are numbered consecutively. This deviates from the
     * way the ccarat output is organized. There the nodes are grouped
     * together and each node group contains all layers. */
    nodes_per_layer = 2*numnp;
    nodes_per_element = field->s9_layers*nodes_per_layer;

#if 0
    coords = (DOUBLE*)CCACALLOC(3*field->numele*nodes_per_element, sizeof(DOUBLE));
#endif

    for (el=0; el<field->numele; ++el) {
      INT n;

      if (fread(x, sizeof(DOUBLE), value_entry_length, field->value_file)!=value_entry_length) {
        dserror("reading value file of discretization %s failed", field->name);
      }

      for (n=0; n<numnp; ++n) {
        INT lay;

        for (lay=0; lay<field->s9_layers; ++lay) {
          INT id;
          INT item;

          /* All GiD ids are one based. */
          id = el*nodes_per_element + lay*nodes_per_layer + n + 1;

          /* ccarat wrote each node coordinate just once */
          /*
           * There is one more set of nodes than layers because we
           * have nodes at the top and at the bottom. */
          item = 3*(n*(field->s9_layers+1) + lay);

          GiD_WriteCoordinates(id, x[item+0], x[item+1], x[item+2]);
#if 0
          coords[3*id+0] = x[item+0];
          coords[3*id+1] = x[item+1];
          coords[3*id+2] = x[item+2];
#endif

          id += numnp;
          item += 3;
          GiD_WriteCoordinates(id, x[item+0], x[item+1], x[item+2]);
#if 0
          coords[3*id+0] = x[item+0];
          coords[3*id+1] = x[item+1];
          coords[3*id+2] = x[item+2];
#endif
        }
      }
    }

    GiD_EndCoordinates();
    GiD_BeginElements();

    /* Now that the coordinates are done create the elements. There is
     * no ccarat output to be read. The nodes determine the elements
     * that we need. */

    for (el=0; el<field->numele; ++el) {
      INT lay;
      INT connect[8];

      for (lay=0; lay<field->s9_layers; ++lay) {
        INT id;
        INT n;

        /* All GiD ids are one based. */
        id = el*field->s9_layers + lay + 1;

        for (n=0; n<numnp; ++n) {

          /* lower layer */
          connect[n] = el*nodes_per_element + lay*nodes_per_layer + n + 1;

          /* upper layer */
          connect[n+numnp] = connect[n] + numnp;
        }

        GiD_WriteElement(id, connect);
      }
    }

    GiD_EndElements();
    GiD_EndMesh();

    CCAFREE(x);

    /*--------------------------------------------------------------------*/
    /* Tell GiD about our element's gauss points.
     * This is not so important here because there are no results at
     * the gauss points (yet). But maybe we'll need them someday... */

    if (field->s9_minor == MINOR_SHELL9_4_22) {
      GP_number = 8;
    }
    else {
      GP_number = 27;
    }
    sprintf(GP_name, "shell9_%d", GP_number);

    GiD_BeginGaussPoint(GP_name, EType, mesh_name, GP_number, 0, 1);
    GiD_EndGaussPoint();

    break;
  case MINOR_SHELL9_8_22:       /* 8-noded shell9 2x2 GP */
  case MINOR_SHELL9_8_33:       /* 8-noded shell9 3x3 GP */
  case MINOR_SHELL9_9_22:       /* 9-noded shell9 2x2 GP */
  case MINOR_SHELL9_9_33:       /* 9-noded shell9 3x3 GP */

    /*
     * Eight and nine node shell9 elements are treated very similar. In
     * both cases there is a set of nodes in the middle of the new
     * brick element. The hex20 element (8 noded shell9) doesn't use
     * the nodes in the center of the sides, but these nodes are there
     * nevertheless. */

    if ((field->s9_minor == MINOR_SHELL9_9_22) ||
        (field->s9_minor == MINOR_SHELL9_9_33)) {
      NNode = 27;
      EType = GiD_Hexahedra;
      numnp = 9;
    }
    else {
      NNode = 20;
      EType = GiD_Hexahedra;
      numnp = 8;
    }

    sprintf(mesh_name, "%s_%s_%d", field->name, "shell9", NNode);

    GiD_BeginMesh(mesh_name, 3, EType, NNode);
    GiD_BeginCoordinates();

    printf("%s: Write node coordinates\n", field->name);
    coords_group = map_read_map(field->table, "shell9_coords");

    value_entry_length = map_read_int(coords_group, "value_entry_length");
    value_offset = map_read_int(coords_group, "value_offset");

    dsassert(value_entry_length == 3*numnp*(2*field->s9_layers+1),
             "invalid shell9 coords entry length");
    x = (DOUBLE*)CCACALLOC(value_entry_length, sizeof(DOUBLE));

    fseek(field->value_file, value_offset, SEEK_SET);

    /* We number the nodes per element. Inside an element the nodes of
     * each layer are numbered consecutively. This deviates from the
     * way the ccarat output is organized. There the nodes are grouped
     * together and each node group contains all layers. */
    nodes_per_layer = 3*numnp;
    nodes_per_element = field->s9_layers*nodes_per_layer;

#if 0
    coords = (DOUBLE*)CCACALLOC(3*field->numele*nodes_per_element, sizeof(DOUBLE));
#endif

    for (el=0; el<field->numele; ++el) {
      INT n;

      if (fread(x, sizeof(DOUBLE), value_entry_length, field->value_file)!=value_entry_length) {
        dserror("reading value file of discretization %s failed", field->name);
      }

      for (n=0; n<numnp; ++n) {
        INT lay;

        for (lay=0; lay<field->s9_layers; ++lay) {
          INT id;
          INT item;

          /* All GiD ids are one based. */
          id = el*nodes_per_element + lay*nodes_per_layer + n + 1;

          /* ccarat wrote each node coordinate just once */
          /* There are two coordinate pairs per layer. */
          item = 3*2*(n*field->s9_layers + lay);

          GiD_WriteCoordinates(id, x[item+0], x[item+1], x[item+2]);
#if 0
          coords[3*id+0] = x[item+0];
          coords[3*id+1] = x[item+1];
          coords[3*id+2] = x[item+2];
#endif

          id += numnp;
          item += 3;
          GiD_WriteCoordinates(id, x[item+0], x[item+1], x[item+2]);
#if 0
          coords[3*id+0] = x[item+0];
          coords[3*id+1] = x[item+1];
          coords[3*id+2] = x[item+2];
#endif

          id += numnp;
          item += 3;
          GiD_WriteCoordinates(id, x[item+0], x[item+1], x[item+2]);
#if 0
          coords[3*id+0] = x[item+0];
          coords[3*id+1] = x[item+1];
          coords[3*id+2] = x[item+2];
#endif
        }
      }
    }

    GiD_EndCoordinates();

    /* Now that the coordinates are done create the elements. There is
     * no ccarat output to be read. The nodes determine the elements
     * that we need. */

    for (el=0; el<field->numele; ++el) {
      INT lay;
      INT connect[27];

      for (lay=0; lay<field->s9_layers; ++lay) {
        INT id;
        INT k;
        INT *ptr;

        /* All GiD ids are one based. */
        id = el*field->s9_layers + lay + 1;

        ptr = connect;

        /*untere Ebene -> Eckknoten*/
        for (k=0; k<4; k++) {
          *ptr++ = el*nodes_per_element + lay*3*numnp + k + 1;
        }
        /*obere Ebene -> Eckknoten*/
        for (k=0; k<4; k++) {
          *ptr++ = el*nodes_per_element + (lay*3+2)*numnp + k + 1;
        }
        /*untere Ebene -> Seitenmittelknoten*/
        for (k=4; k<8; k++) {
          *ptr++ = el*nodes_per_element + (lay*3)*numnp + k + 1;
        }
        /*mittlere Ebene -> Eckknoten*/
        for (k=0; k<4; k++) {
          *ptr++ = el*nodes_per_element + (lay*3+1)*numnp + k + 1;
        }
        /*obere Ebene -> Seitenmittelknoten*/
        for (k=4; k<8; k++) {
          *ptr++ = el*nodes_per_element + (lay*3+2)*numnp + k + 1;
        }

        /* write additional node numbers for quad9 */
        if ((field->s9_minor == MINOR_SHELL9_9_22) ||
            (field->s9_minor == MINOR_SHELL9_9_33)) {

          /*untere Ebene -> Mittelknoten*/
          for (k=8; k<9; k++) {
            *ptr++ = el*nodes_per_element + (lay*3)*numnp + k + 1;
          }
          /*mittlere Ebene -> Seitenmittelknoten*/
          for (k=4; k<8; k++) {
            *ptr++ = el*nodes_per_element + (lay*3+1)*numnp + k + 1;
          }
          /*obere Ebene -> Mittelknoten*/
          for (k=8; k<9; k++) {
            *ptr++ = el*nodes_per_element + (lay*3+2)*numnp + k + 1;
          }
          /*mittlere Ebene -> Mittelknoten*/
          for (k=8; k<9; k++) {
            *ptr++ = el*nodes_per_element + (lay*3+1)*numnp + k + 1;
          }
        }

        GiD_WriteElement(id, connect);
      }
    }

    GiD_EndMesh();

    CCAFREE(x);

    /*--------------------------------------------------------------------*/
    /* Tell GiD about our element's gauss points.
     * This is not so important here because there are no results at
     * the gauss points (yet). But maybe we'll need them someday... */

    if ((field->s9_minor == MINOR_SHELL9_8_22) ||
        (field->s9_minor == MINOR_SHELL9_9_22)) {
      GP_number = 8;
    }
    else {
      GP_number = 27;
    }
    sprintf(GP_name, "shell9_%d", GP_number);

    GiD_BeginGaussPoint(GP_name, EType, mesh_name, GP_number, 0, 1);
    GiD_EndGaussPoint();

    break;
  default:
    dserror("unknown minor number for shell9 element", field->s9_minor);
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the displacements of the unconnected brick mesh.

  Danger! If there are other discretizations this code will fail!
  That's because we don't know which node and element ids are assigned
  in those discretizations. In the normal case the value
  field->first_node_id tells us where our range starts. But here nodes
  and elements are generated that ccarat never knew about.

  This problem could be solved inside the filter. For now it's not
  worth the trouble.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_shell9_displacement_unsmoothed(FIELD_DATA *field, MAP* group)
{
  INT el;
  MAP* disp_group;
  INT value_entry_length;
  DOUBLE* disp;
  DOUBLE time;
  INT step;
  INT numnp;
  INT nodes_per_element;

  CHAR* componentnames[] = { "x-displ", "y-displ", "z-displ" };

  disp_group = map_read_map(group, "shell9_displacement");

  dsassert(map_read_int(disp_group, "size_entry_length") == 0,
           "stress size mismatch");

  value_entry_length = map_read_int(disp_group, "value_entry_length");
  disp = (DOUBLE*)CCACALLOC(value_entry_length, sizeof(DOUBLE));

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  printf("%s: Write displacement of step %d\n", field->name, step);

  GiD_BeginResult("displacement", "ccarat", step, GiD_Vector, GiD_OnNodes,
                  NULL, NULL, 3, componentnames);
  fseek(field->value_file, map_read_int(disp_group, "value_offset"), SEEK_SET);

  switch (field->s9_minor) {
  case MINOR_SHELL9_4_22:       /* 4-noded shell9 2x2 GP */
  case MINOR_SHELL9_4_33:       /* 4-noded shell9 3x3 GP */

    numnp = 4;

    dsassert(value_entry_length == 3*numnp*(field->s9_layers+1),
             "invalid shell9 displacement entry length");

    /* We number the nodes per element. Inside an element the nodes of
     * each layer are numbered consecutively. This deviates from the
     * way the ccarat output is organized. There the nodes are grouped
     * together and each node group contains all layers. */
    nodes_per_element = 2*field->s9_layers*numnp;

    /* Iterate all elements */
    for (el=0; el<field->numele; ++el) {
      INT lay;

      if (fread(disp, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
        dserror("failed to read displacement of element %d", el);
      }

      for (lay=0; lay<field->s9_layers; ++lay) {
        INT n;

        for (n=0; n<numnp; ++n) {
          INT id;
          INT item;

          /* lower layer */
          /* All GiD ids are one based. */
          id = el*nodes_per_element + 2*lay*numnp + n + 1;
          item = 3*(n*(field->s9_layers+1) + lay);

          GiD_WriteVector(id, disp[item+0], disp[item+1], disp[item+2]);

          /* upper layer */
          id += numnp;
          item += 3;
          GiD_WriteVector(id, disp[item+0], disp[item+1], disp[item+2]);
        }
      }
    }

    break;
  case MINOR_SHELL9_8_22:       /* 8-noded shell9 2x2 GP */
  case MINOR_SHELL9_8_33:       /* 8-noded shell9 3x3 GP */
  case MINOR_SHELL9_9_22:       /* 9-noded shell9 2x2 GP */
  case MINOR_SHELL9_9_33:       /* 9-noded shell9 3x3 GP */

    if ((field->s9_minor == MINOR_SHELL9_9_22) ||
        (field->s9_minor == MINOR_SHELL9_9_33)) {
      numnp = 9;
    }
    else {
      numnp = 8;
    }

    dsassert(value_entry_length == 3*numnp*(2*field->s9_layers+1),
             "invalid shell9 displacement entry length");

    /* We number the nodes per element. Inside an element the nodes of
     * each layer are numbered consecutively. This deviates from the
     * way the ccarat output is organized. There the nodes are grouped
     * together and each node group contains all layers. */
    nodes_per_element = 3*field->s9_layers*numnp;

    for (el=0; el<field->numele; ++el) {
      INT n;

      if (fread(disp, sizeof(DOUBLE), value_entry_length, field->value_file)!=value_entry_length) {
        dserror("reading value file of discretization %s failed", field->name);
      }

      for (n=0; n<numnp; ++n) {
        INT lay;

        for (lay=0; lay<field->s9_layers; ++lay) {
          INT id;
          INT item;

          /* All GiD ids are one based. */
          id = el*nodes_per_element + 3*lay*numnp + n + 1;

          /* ccarat wrote each node coordinate just once */
          /* There are two coordinate pairs per layer. */
          item = 3*2*(n*field->s9_layers + lay);

          GiD_WriteVector(id, disp[item+0], disp[item+1], disp[item+2]);

          id += numnp;
          item += 3;
          GiD_WriteVector(id, disp[item+0], disp[item+1], disp[item+2]);

          id += numnp;
          item += 3;
          GiD_WriteVector(id, disp[item+0], disp[item+1], disp[item+2]);
        }
      }
    }

    break;
  default:
    dserror("unknown minor number for shell9 element", field->s9_minor);
  }

  GiD_EndResult();

  CCAFREE(disp);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write extrapolated stresses of the unconnected brick mesh.

  Here we really have to do some calculation. The values written by
  ccarat are the gauss point values. However, we waht to have values
  at the nodes (because GiD is known to handle these) and thus the
  extrapolation needs to be done. This used to be done by ccarat, but
  this way it's much more flexible. (Maybe GiD will understand the
  gauss point values some day.)

  Things are getting messy here. The original code is more restrictive
  that at other places, the strange element variants (number of nodes
  and gauss points) are not supported.

  This is the second half of ``s9_stress``.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_shell9_stress_unsmoothed(FIELD_DATA *field, MAP* group)
{
  INT gaussperm4[4] = {3,1,0,2};
  INT gaussperm9[9] = {8,2,0,6,5,1,3,7,4};

  INT el;
  MAP* stress_group;
  INT value_entry_length;
  DOUBLE* stress;
  DOUBLE time;
  INT step;
  S9_DATA data;

  /* The number of gauss points used for extrapolation. They coincide
   * with the number of nodes in r and s dimension. */
  INT nir_x;
  INT nis_x;
  INT nit_x;

  INT stress_row_length;

  INT numnp;
  INT nodes_per_layer;
  INT nodes_per_element;

  DIS_TYP distyp;

  /* shape functions */
  ARRAY funct_a;
  DOUBLE *funct;

  /* derivatives of shape functions */
  ARRAY deriv_a;
  DOUBLE **deriv;

  /* element array for stresses on nodal points */
  ARRAY strK_a;
  DOUBLE **gp_strK;

  /* interpolated values in plane before interpolation in thickness */
  DOUBLE gp_u[6][9], gp_o[6][9];

  CHAR* xyz_componentnames[] = { "Stress-xx","Stress-xy","Stress-yy","Stress-xz","Stress-yz","Stress-zz" };
  CHAR* rst_componentnames[] = { "Stress-rr","Stress-rs","Stress-ss","Stress-rt","Stress-st","Stress-tt" };
  CHAR** componentnames;

  if (strcmp(field->s9_forcetype, "xyz")==0) {
    componentnames = xyz_componentnames;
  }
  else {
    componentnames = rst_componentnames;
  }

  /* read the control data */
  stress_group = map_read_map(group, "stress");

  dsassert(map_read_int(stress_group, "size_entry_length") == 0,
           "stress size mismatch");

  value_entry_length = map_read_int(stress_group, "value_entry_length");
  stress = (DOUBLE*)CCACALLOC(value_entry_length, sizeof(DOUBLE));

  time = map_read_real(group, "time");
  step = map_read_int(group, "step");

  printf("%s: Write stress of step %d\n", field->name, step);

  GiD_BeginResult("shell9_stresses", "ccarat", step, GiD_Matrix, GiD_OnNodes,
                  NULL, NULL, 6, componentnames);
  fseek(field->value_file, map_read_int(stress_group, "value_offset"), SEEK_SET);


  /* integration parameters for stress extrapolation */
  switch (field->s9_minor) {
  case MINOR_SHELL9_4_22:       /* 4-noded shell9 2x2 GP */
    nir_x = 2;
    nis_x = 2;
    nit_x = 2;
    distyp = quad4;
    s9intg_str(&data,4);

    numnp = 4;
    nodes_per_layer = 2*numnp;
    nodes_per_element = field->s9_layers*nodes_per_layer;
    break;
  case MINOR_SHELL9_4_33:       /* 4-noded shell9 3x3 GP */
  case MINOR_SHELL9_8_22:       /* 8-noded shell9 2x2 GP */
  case MINOR_SHELL9_9_22:       /* 9-noded shell9 2x2 GP */
    dserror("minor type %d not supported", field->s9_minor);
    break;
  case MINOR_SHELL9_8_33:       /* 8-noded shell9 3x3 GP */
  case MINOR_SHELL9_9_33:       /* 9-noded shell9 3x3 GP */
    nir_x = 3;
    nis_x = 3;
    nit_x = 2;
    distyp = quad9;
    s9intg_str(&data,9);

    /* we need to make a difference when writing */
    numnp = (field->s9_minor == MINOR_SHELL9_8_33) ? 8 : 9;
    nodes_per_layer = 3*numnp;
    nodes_per_element = field->s9_layers*nodes_per_layer;
    break;
  default:
    dserror("unknown minor number for shell9 element", field->s9_minor);
  }

  /* setup supporting arrays */
  funct = amdef("funct", &funct_a,nir_x*nis_x,1,"DV");
  deriv = amdef("deriv", &deriv_a,2,nir_x*nis_x,"DA");
  gp_strK = amdef("gp_strK", &strK_a, 6, nir_x*nis_x*nit_x*field->s9_layers,"DA");

  /* The length of one row in the stress array. That's the number of
   * gauss points per shell9 element (per layer). */
  stress_row_length = nir_x*nis_x*nit_x*field->s9_layers;

  dsassert(value_entry_length == 6*stress_row_length,
           "stress value entry mismatch");

  /* write the stresses for all shell9 elements */
  for (el=0; el<field->numele; ++el) {
    INT lay;

    if (fread(stress, sizeof(DOUBLE), value_entry_length, field->value_file) != value_entry_length) {
      dserror("failed to read stress of element %d", el);
    }

    /* Each layer in the shell9 element corresponds to one brick
     * element in the output. */
    for (lay=0; lay<field->s9_layers; lay++) {
      INT j;
      INT lt;
      for (lt=0; lt<nit_x; lt++) {
        INT lr;
        INT ngauss = 0;
        for (lr=0; lr<nir_x; lr++) {
          INT ls;
          DOUBLE e1 = data.xgpr[lr];
          if (e1 != 0.0)
            e1 = 1./e1;

          for (ls=0; ls<nis_x; ls++) {
            DOUBLE e2 = data.xgps[ls];
            if (e2 != 0.0)
              e2 = 1./e2;

            /* shape functions at gp e1,e2 on mid surface */
            s9_funct_deriv(funct,deriv,e1,e2,distyp,0);

            /* interpolate in this plane to nodal points */
            for (j=0; j<6; j++) {
              INT k;

              /* stress at GPs in one plane */
              DOUBLE* strGP;
              DOUBLE strK = 0.0;

              /*
               * There had been a loop. But it copied consecutive
               * values to another array. We can avoid that. */
              strGP = &(stress[stress_row_length*j + (2*lay + lt) * (nir_x*nis_x)]);

              if (distyp == quad4)
                for (k=0; k<4; k++)
                  strK += funct[k]*strGP[gaussperm4[k]];

              /*
               * The midpoint value for quad8 is not written but we
               * handle that later. */
              else if (distyp == quad9)
                for (k=0; k<9; k++)
                  strK += funct[k]*strGP[gaussperm9[k]];

              /* store value to finally interpolate in thickness direction */
              /* lower part */
              if (lt == 0)  gp_u[j][ngauss] = strK;
              /* upper part */
              if (lt == 1)  gp_o[j][ngauss] = strK;
            }
            ngauss++;
          }
        }
      }

      /* now interpolate the gp_u and gp_o values to the nodes in thickness direction */
      /*interpolate the 6-stress komponents*/
      for (j=0; j<6; j++) {
        INT k;
        for (k=0; k<(nir_x*nis_x); k++) {
          DOUBLE P_u, P_o;
          DOUBLE N_u, N_o;
          DOUBLE m;
          INT ID_stress;

          P_u = gp_u[j][k];
          P_o = gp_o[j][k];
          m   = (P_o - P_u) / (2./sqrt(3.));
          N_u = P_u - m * ( 1. - 1./sqrt(3.));
          N_o = P_o + m * ( 1. - 1./sqrt(3.));

          /* write to gp_strK array */
          ID_stress = k + (2*lay + 0) * (nir_x*nis_x);
          gp_strK[j][ID_stress] = N_u;
          ID_stress = k + (2*lay + 1) * (nir_x*nis_x);
          gp_strK[j][ID_stress] = N_o;
        }
      }
    }

    /* now write the stresses at the nodes gp_strK back to ele_stress array */
    for (lay=0; lay<field->s9_layers; lay++) {
      INT j;
      INT lt;

      /*
       * We have to collect the stresses for upper and lower nodes in
       * one brick element because the middle nodes' stresses (in case
       * of quad8 or quad9) are just the average of upper and lower
       * ones. It's that bad.
       *
       * The lower nodes are in the first half, the upper ones come
       * after them. To each node there are six consecutive stress
       * values. */
      DOUBLE node_stress[6*9*2];

      for (lt=0; lt<nit_x; lt++) {
        INT lr;

        /*
         * This gives the number of the gauss point we are looking
         * at. But it's the node number as well. */
        INT ngauss = 0;

        /* quad8 does not have a ninth node. */
        if (ngauss < numnp) {

          /*
           * Here we loop the gauss points in r and s dimension but keep
           * in mind that we can use the indices for node numbers, too. */
          for (lr=0; lr<nir_x; lr++) {
            INT ls;
            for (ls=0; ls<nis_x; ls++) {
              INT ID_stress_perm;

              if (distyp == quad4)
                ID_stress_perm = gaussperm4[ngauss] + (2*lay + lt) * (nir_x*nis_x);
              else if (distyp == quad9)
                ID_stress_perm = gaussperm9[ngauss] + (2*lay + lt) * (nir_x*nis_x);

              /* loop the stress components */
              for (j=0; j<6; j++) {
                node_stress[6*lt*numnp + 6*ngauss + j] = gp_strK[j][ID_stress_perm];
              }

              ngauss++;
            }
          }
        }
      }

      /* Write the extrapolated values for the current brick element. */
      switch (numnp) {
      case 4:
        for (j=0; j<numnp; ++j) {
          INT id;
          INT item;

          id = el*nodes_per_element + lay*nodes_per_layer + j + 1;
          item = 6*j;
          GiD_Write3DMatrix(id,
                            node_stress[item+0], node_stress[item+1], node_stress[item+2],
                            node_stress[item+3], node_stress[item+4], node_stress[item+5]);

          id += numnp;
          item += 6*numnp;
          GiD_Write3DMatrix(id,
                            node_stress[item+0], node_stress[item+1], node_stress[item+2],
                            node_stress[item+3], node_stress[item+4], node_stress[item+5]);
        }
        break;
      case 8:
      case 9:
        for (j=0; j<numnp; ++j) {
          INT id;
          INT item;

          /* lower nodes */
          id = el*nodes_per_element + lay*3*numnp + j + 1;
          item = 6*j;
          GiD_Write3DMatrix(id,
                            node_stress[item+0], node_stress[item+1], node_stress[item+2],
                            node_stress[item+3], node_stress[item+4], node_stress[item+5]);

          /* upper nodes */
          id += 2*numnp;
          item += 6*numnp;
          GiD_Write3DMatrix(id,
                            node_stress[item+0], node_stress[item+1], node_stress[item+2],
                            node_stress[item+3], node_stress[item+4], node_stress[item+5]);
        }
        for (j=0; j<numnp; ++j) {
          INT id;
          INT item;

          /* middle nodes: interpolate */
          id = el*nodes_per_element + (lay*3+1)*numnp + j + 1;
          item = 6*j;
          GiD_Write3DMatrix(id,
                            .5*(node_stress[item+0]+node_stress[6*numnp+item+0]),
                            .5*(node_stress[item+1]+node_stress[6*numnp+item+1]),
                            .5*(node_stress[item+2]+node_stress[6*numnp+item+2]),
                            .5*(node_stress[item+3]+node_stress[6*numnp+item+3]),
                            .5*(node_stress[item+4]+node_stress[6*numnp+item+4]),
                            .5*(node_stress[item+5]+node_stress[6*numnp+item+5]));
        }
        break;
      default:
        dserror("Ups!");
      }
    }
  }

  GiD_EndResult();

  CCAFREE(stress);
  amdel(&strK_a);
  amdel(&funct_a);
  amdel(&deriv_a);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a shell9 discretization without
  smoothing.

  Danger! If there are other discretizations this code will fail!
  That's because we don't know which node and element ids are assigned
  in those discretizations. In the normal case the value
  field->first_node_id tells us where our range starts. But here nodes
  and elements are generated that ccarat never knew about.

  This problem could be solved inside the filter. For now it's not
  worth the trouble.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void write_shell9_field_unsmoothed(PROBLEM_DATA* problem, INT num)
{
  FIELD_DATA* field;
  INT i;

  field = &(problem->discr[num]);

  /*--------------------------------------------------------------------*/
  /* Have the mesh written. It's going to be a "mesh" of unconnected
   * brick elements, one per shell9 element and layer. */

  write_shell9_mesh_unsmoothed(problem, num);

  /*--------------------------------------------------------------------*/
  /* Now it's time to write the time dependend results. */

  for (i=0; i<problem->num_results; ++i) {
    MAP* result_group;
    result_group = problem->result_group[i];

    /* We iterate the list of all results. Here we are interested in
     * the results of this discretization. */
    if ((strcmp(map_read_string(result_group, "field"),
                fieldnames[field->type]) == 0) &&
        (map_read_int(result_group, "discretization") == field->disnum)) {

      if (map_has_map(result_group, "shell9_displacement")) {
        write_shell9_displacement_unsmoothed(field, result_group);
      }
      if (map_has_map(result_group, "stress")) {
        write_shell9_stress_unsmoothed(field, result_group);
      }
    }
  }

#if 0
  CCAFREE(coords);
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Starting point for writing a shell9 discretization.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void write_shell9_field(PROBLEM_DATA* problem, INT num)
{
  FIELD_DATA* field;

  field = &(problem->discr[num]);
  if (field->s9_smooth_results) {
    write_shell9_field_smoothed(problem, num);
  }
  else {
    write_shell9_field_unsmoothed(problem, num);
  }
}


#endif
