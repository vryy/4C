/*!
\file
\brief Support functions for binary IO to a single file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Here are the element specific functions that are called during binary
input and output. Most of these functions are at the bottom of the io
module and are not called by an algorithm directly. They provide the
link between the specific ccarat data structures and the generic io
module, so they need to be changed when changes to the elements are
made.

These functions and definitions are collected here intentionally. It's
one place to go in chase of changes. Of course much of this code could
be put to the individual elements as well. But we'd loose the strong
connection to the io module, overload the element routines with
support stuff and might loose even more performance due to per element
function calls.

\author u.kue
\date 08/04

*/

#ifdef BINIO

#include "io_packing.h"
#include "io_singlefile.h"

#include "../pss_full/pss_table.h"

#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../beam3/beam3.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../ale2/ale2.h"
#include "../ale3/ale3.h"
#include "../axishell/axishell.h"
#include "../interf/interf.h"
#include "../wallge/wallge.h"
#include "../ls/ls.h"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
struct _MATERIAL     *mat;

/*!----------------------------------------------------------------------
\brief vector of multilayer material law

<pre>                                                            sh 10/02
This structure struct _MULTIMAT  *multimat is defined in global_control.c
and the type is in materials.h
It holds all information about the layered section data
</pre>
*----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | struct _ALLDYNA       *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA *alldyn;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for output.

  This structure needs to be initialized at startup. The whole output
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_OUT_MAIN bin_out_main;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for input.

  This structure needs to be initialized at startup. The whole input
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_IN_MAIN bin_in_main;


/*----------------------------------------------------------------------*/
/*!
  \brief All fields names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern CHAR* fieldnames[];



/*----------------------------------------------------------------------*/
/*!
  \brief Mark all element types that are used in this discretization.

  Here the context's element_flag table is set up. That is all
  elements are checked and the corresponding entry (selected by
  major/minor type number) in the flag table gets a tic.

  The result is allreduced.

  \param context  (i/o) pointer to an output context in the setup phase
  \param actintra (i)   the communicator
  \param actpdis  (i)   the partitions discretization

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_find_element_types(struct _BIN_OUT_FIELD *context,
                            INTRA *actintra,
                            PARTDISCRET *actpdis)
{
  INT i;
  ELEMENT_FLAGS element_flag;

#ifdef DEBUG
  dstrc_enter("out_find_element_types");
#endif

  memset(element_flag, 0, ELEMENT_FLAGS_SIZE);

  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];

    switch (actele->eltyp) {
    case el_shell8:
      element_flag[el_shell8][MINOR_SHELL8(actele)] = 1;
      break;
#ifdef D_SHELL9
    case el_shell9:
      element_flag[el_shell9][MINOR_SHELL9(actele)] = 1;
      break;
#endif
#ifdef D_BRICK1
    case el_brick1:
      element_flag[el_brick1][MINOR_BRICK1(actele)] = 1;
      break;
#endif
#ifdef D_FLUID2
    case el_fluid2:
      element_flag[el_fluid2][MINOR_FLUID2(actele)] = 1;
      break;
#endif
#ifdef D_XFEM
    case el_fluid2_xfem:
      element_flag[el_fluid2_xfem][MINOR_FLUID2_XFEM(actele)] = 1;
      break;
#endif
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      element_flag[el_fluid2_pro][MINOR_FLUID2_PRO(actele)] = 1;
      break;
#endif
#ifdef D_FLUID3
    case el_fluid3:
      element_flag[el_fluid3][MINOR_FLUID3(actele)] = 1;
      break;
#endif
#ifdef D_ALE
    case el_ale2:
      element_flag[el_ale2][MINOR_ALE2(actele)] = 1;
      break;
    case el_ale3:
      element_flag[el_ale3][MINOR_ALE3(actele)] = 1;
      break;
#endif
#ifdef D_WALL1
    case el_wall1:
      element_flag[el_wall1][MINOR_WALL1(actele)] = 1;
      break;
#endif
#ifdef D_BEAM3
    case el_beam3:
      element_flag[el_beam3][MINOR_BEAM3(actele)] = 1;
      break;
#endif
#ifdef D_AXISHELL
    case el_axishell:
      element_flag[el_axishell][MINOR_AXISHELL] = 1;
      break;
#endif
#ifdef D_INTERF
    case el_interf:
      element_flag[el_interf][MINOR_INTERF(actele)] = 1;
      break;
#endif
#ifdef D_WALLGE
    case el_wallge:
      element_flag[el_wallge][MINOR_WALLGE(actele)] = 1;
      break;
#endif
#ifdef D_LS
    case el_ls2:
      element_flag[el_ls2][MINOR_LS2(actele)] = 1;
      break;
#endif

    default:
      dserror("Unknown type %d of element", actele->eltyp);
      break;
    }
  }

  /* memcpy takes the buffer length in bytes. MPI asks for the number
   * of entries. */
#ifdef PARALLEL
  MPI_Allreduce(element_flag, context->element_flag, el_count*MAX_EL_MINOR, MPI_INT,
                MPI_MAX, actintra->MPI_INTRA_COMM);
#else
  memcpy(context->element_flag, element_flag, ELEMENT_FLAGS_SIZE);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/* Here we have lots of functions that do nothing but packing
 * something (node arrays, element stuff, distributed vectors) into
 * the output buffers. These functions share a common structure. See
 * ``out_pack_items`` for an explanation of their parameters.
 *
 * These are the functions that need to be changed according to
 * changes in the elements. */
/*----------------------------------------------------------------------*/



/*----------------------------------------------------------------------*/
/*!
  \brief Pack the values of the indicated node array in some buffer to
  send them to their writing processor.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_node_arrays(BIN_OUT_CHUNK *chunk,
                                 NODE_ARRAY array,
                                 PARTDISCRET *actpdis,
                                 DOUBLE *send_buf,
                                 INT send_count,
                                 INT *send_size_buf,
                                 INT dst_first_id,
                                 INT dst_num)
{
  INT j, counter;

#ifdef DEBUG
  dstrc_enter("out_pack_node_arrays");
#endif

#ifdef PARALLEL
#define pack_send_size_buf(node_array)                          \
      send_size_buf[3*counter  ] = actnode->Id_loc;             \
      send_size_buf[3*counter+1] = actnode->node_array.fdim;    \
      send_size_buf[3*counter+2] = actnode->node_array.sdim;
#else
#define pack_send_size_buf(node_array)                          \
      send_size_buf[2*counter  ] = actnode->node_array.fdim;    \
      send_size_buf[2*counter+1] = actnode->node_array.sdim;
#endif

  /* gather the nodes to be send to processor i */
  counter = 0;
#define call_gather_nodes(node_array)                                   \
  dsassert(chunk->size_entry_length == 2, "invalid size entry length"); \
  for (j=0; j<actpdis->numnp; ++j) {                                    \
    NODE* actnode = actpdis->node[j];                                   \
    if ((actnode->Id_loc >= dst_first_id) &&                            \
        (actnode->Id_loc < dst_first_id+dst_num)) {                     \
      DOUBLE *src_ptr;                                                  \
      DOUBLE *dst_ptr;                                                  \
      INT k;                                                            \
      INT size = actnode->node_array.fdim*actnode->node_array.sdim;     \
      dsassert(size <= chunk->field->max_size[node_array_ ## node_array], \
               "outdated size calculation. panic.");                    \
      src_ptr = actnode->node_array.a.da[0];                            \
      dst_ptr = &(send_buf[chunk->value_entry_length*counter]);         \
      pack_send_size_buf(node_array);                                   \
      for (k=0; k<size; ++k) {                                          \
        *dst_ptr++ = *src_ptr++;                                        \
      }                                                                 \
      counter += 1;                                                     \
    }                                                                   \
  }                                                                     \
  dsassert(counter*chunk->value_entry_length == send_count,             \
           "node pack count mismatch");


  switch (array) {
  case node_array_sol:
    call_gather_nodes(sol);
    break;
  case node_array_sol_increment:
    call_gather_nodes(sol_increment);
    break;
  case node_array_sol_residual:
    call_gather_nodes(sol_residual);
    break;
  case node_array_sol_mf:
    call_gather_nodes(sol_mf);
    break;
  default:
    dserror("Node array %d unknown", array);
  }

#undef call_gather_nodes
#undef pack_send_size_buf

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the mesh (element id, element type and connectivity) in
  some buffer to send them to their writing processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_mesh(BIN_OUT_CHUNK *chunk,
                          PARTDISCRET *actpdis,
                          DOUBLE *send_buf,
                          INT send_count,
                          INT *send_size_buf,
                          INT dst_first_id,
                          INT dst_num)
{
  INT i;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_mesh");
#endif

  dsassert(chunk->size_entry_length > 3, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  len = chunk->size_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      INT* ptr;
      INT j;

#ifdef PARALLEL
      send_size_buf[(len+1)*counter] = actele->Id_loc;
      ptr = &(send_size_buf[(len+1)*counter+1]);
#else
      ptr = &(send_size_buf[len*counter]);
#endif

      *ptr++ = actele->Id;
      *ptr++ = actele->eltyp;
      *ptr++ = FIND_MINOR(actele);

      dsassert(actele->numnp+3 <= len, "size entry too short");

      for (j=0; j<actele->numnp; ++j) {
        /* We have to store the global Id here. */
        *ptr++ = actele->node[j]->Id;
      }

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the node id and its coordinates in some buffer to send
  them to their writing processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_coords(BIN_OUT_CHUNK *chunk,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_coords");
#endif

  dsassert(chunk->size_entry_length == 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == genprob.ndim, "invalid value entry length");

  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i) {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num)) {
      DOUBLE* ptr;
      INT j;

#ifdef PARALLEL
      send_size_buf[2*counter  ] = actnode->Id_loc;
      send_size_buf[2*counter+1] = actnode->Id;
#else
      send_size_buf[counter] = actnode->Id;
#endif

      ptr = &(send_buf[len*counter]);
      for (j=0; j<len; ++j) {
        *ptr++ = actnode->x[j];
      }

      counter += 1;
    }
  }

  dsassert(counter*len == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the displacements to send them to their writing processor.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_displacement(BIN_OUT_CHUNK *chunk,
                                  INT place,
                                  PARTDISCRET *actpdis,
                                  DOUBLE *send_buf,
                                  INT send_count,
                                  INT *send_size_buf,
                                  INT dst_first_id,
                                  INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_displacement");
#endif

  ndim = chunk->value_entry_length;

#if 0
#ifdef D_AXISHELL
  /* An ugly hack! The old out_gid_sol used to change genprob.ndim
   * itself which is very much worse. We won't do this here. However,
   * this way we might miss some side effects the old version was
   * relying on. We'll see... */
  ndim = MAX(ndim, 3);
#endif
#endif

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num)) {
      INT k;
      DOUBLE *ptr = actnode->sol.a.da[place];

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim >= ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      for (k=0; k<ndim; ++k) {
        send_buf[ndim*counter+k] = *ptr++;
      }
      counter += 1;
    }
  }
  dsassert(counter*ndim == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


#ifdef D_SHELL8

/*----------------------------------------------------------------------*/
/*!
  \brief Pack the directors.

  There are three director coordinates per node per
  element. Additionally we store one thickness value at each node.

  Using these values it's possible to visualize shell8 elements using
  bricks.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell8_director(BIN_OUT_CHUNK *chunk,
                                     PARTDISCRET *actpdis,
                                     DOUBLE *send_buf,
                                     INT send_count,
                                     INT *send_size_buf,
                                     INT dst_first_id,
                                     INT dst_num)
{
  INT i;
  INT k;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_shell8_director");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      SHELL8* s8;
      INT minor;
      DOUBLE *dst_ptr;
      INT numnp;

      s8 = actele->e.s8;
      minor = FIND_MINOR(actele);
      numnp = element_info[el_shell9].variant[minor].node_number;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      for (k=0; k<numnp; ++k) {
        *dst_ptr++ = s8->a3ref.a.da[0][k];
        *dst_ptr++ = s8->a3ref.a.da[1][k];
        *dst_ptr++ = s8->a3ref.a.da[2][k];
      }
      for (k=0; k<numnp; ++k) {
        /*
         * We could read the per node thickness here as well but this
         * is what the old ccarat output does. */
        *dst_ptr++ = s8->thick;
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


#ifdef D_SHELL9


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the node coordinates of the shell9 element's inner nodes.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell9_coords(BIN_OUT_CHUNK *chunk,
                                   PARTDISCRET *actpdis,
                                   DOUBLE *send_buf,
                                   INT send_count,
                                   INT *send_size_buf,
                                   INT dst_first_id,
                                   INT dst_num)
{
  INT i;
  INT k;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_shell9_coords");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      SHELL9* s9;
      INT minor;
      DOUBLE *dst_ptr;

      s9 = actele->e.s9;
      minor = FIND_MINOR(actele);

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      /* this is ``s9_out_gid_allcoords_unsmo`` in disguise */

      switch (minor) {
      case MINOR_SHELL9_4_22:   /* 4-noded shell9 2x2 GP */
      case MINOR_SHELL9_4_33:   /* 4-noded shell9 3x3 GP */

        for (k=0; k<actele->numnp; k++) {
          NODE *actnode;
          INT klay;
          INT jlay;
          INT l;

          /* director a3 of current layer in ref. config */
          DOUBLE a3ref_l[3];

          /* coordinats */
          DOUBLE x[3],x_u[3],x_o[3];

          actnode = actele->node[k];

          /* get director of node
             NOTE: in the undeformed geometrie, the directors of the
             different kinematic layers are equal to the director in the
             total reference layer */
          a3ref_l[0]= s9->a3ref.a.da[0][k];
          a3ref_l[1]= s9->a3ref.a.da[1][k];
          a3ref_l[2]= s9->a3ref.a.da[2][k];

          /* calculate the coordinates of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the coordinates to values of middle surface */
          for (l=0; l<3; l++)
            x[l] = actnode->x[l];

          /* continuity matrix for current kin layer at bot */
          for (jlay=0; jlay<s9->num_klay; jlay++) {
            DOUBLE klayhgt;
            DOUBLE hl;
            DOUBLE e3;
            DOUBLE zeta;

            klayhgt = s9->klayhgt[jlay];
            hl      = (klayhgt/100)*s9->thick;
            /* hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

            /* calculate coordinates */
            for (l=0; l<3; l++)
              x[l] += zeta*hl*a3ref_l[l];
          }

          /* write coordinates of first kinematic layer at bottom */
          *dst_ptr++ = x[0];
          *dst_ptr++ = x[1];
          *dst_ptr++ = x[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            x_u[l] = x[l];

          /* calculate the coordinates of the top surface of the each kin layer */
          for (klay=0; klay<s9->num_klay; klay++) {
            INT num_mlay;
            DOUBLE* mlayhgt;
            DOUBLE sum_hgt;
            INT mlay;

            /* initialize the coordinates to values of middle surface */
            for (l=0; l<3; l++)
              x[l] = actnode->x[l];

            /* continuity matrix for aktual kin layer at top */
            for (jlay=0; jlay<s9->num_klay; jlay++) {
              DOUBLE klayhgt;
              DOUBLE hl;
              DOUBLE e3;
              DOUBLE zeta;

              klayhgt = s9->klayhgt[jlay];
              hl      = (klayhgt/100)*s9->thick;
              /*hl      = 0.5*hl; */  /*norm of director is half of kinematic layer hight*/
              hl      = A3FAC_SHELL9*hl;
              e3 = +1.0; /*top*/
              zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

              /* calculate coordinates */
              for (l=0; l<3; l++)
                x[l] += zeta*hl*a3ref_l[l];
            }

            /* calculate the coordinates of the material layers of
             * aktual kinematic layer */
            num_mlay = s9->kinlay[klay].num_mlay;
            mlayhgt  = s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save coordinates for interpolation for the material layers */
            for (l=0; l<3; l++)
              x_o[l] = x[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt  = 0.0; /*initialize the total hgt*/

            /* only if more than one material layer to this kinematic layer */
            for (mlay=1; mlay<num_mlay; mlay++) {
              sum_hgt += mlayhgt[mlay-1];
              for (l=0; l<3; l++)
                x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);

              /* write the bottom coordinates of this material layer */
              *dst_ptr++ = x[0];
              *dst_ptr++ = x[1];
              *dst_ptr++ = x[2];
            }

            /* write the top coordinates of this kinematic layer */
            *dst_ptr++ = x_o[0];
            *dst_ptr++ = x_o[1];
            *dst_ptr++ = x_o[2];

            for (l=0; l<3; l++)
              x_u[l] = x_o[l];
          }
        }
        break;
      case MINOR_SHELL9_8_22:   /* 8-noded shell9 2x2 GP */
      case MINOR_SHELL9_8_33:   /* 8-noded shell9 3x3 GP */
      case MINOR_SHELL9_9_22:   /* 9-noded shell9 2x2 GP */
      case MINOR_SHELL9_9_33:   /* 9-noded shell9 3x3 GP */

        for (k=0; k<actele->numnp; k++) {
          NODE *actnode;
          INT klay;
          INT jlay;
          INT l;

          /* director a3 of current layer in ref. config */
          DOUBLE a3ref_l[3];

          /* coordinats */
          DOUBLE x[3],x_u[3],x_o[3];

          actnode = actele->node[k];

          /* get director of node
             NOTE: in the undeformed geometrie, the directors of the
             different kinematic layers are equal to the director in the
             total reference layer */
          a3ref_l[0]= s9->a3ref.a.da[0][k];
          a3ref_l[1]= s9->a3ref.a.da[1][k];
          a3ref_l[2]= s9->a3ref.a.da[2][k];

          /* calculate the coordinates of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the coordinates to values of middle surface */
          for (l=0; l<3; l++)
            x[l] = actnode->x[l];

          /* continuity matrix for aktual kin layer at bot */
          for (jlay=0; jlay<s9->num_klay; jlay++) {
            DOUBLE klayhgt;
            DOUBLE hl;
            DOUBLE e3;
            DOUBLE zeta;

            klayhgt = s9->klayhgt[jlay];
            hl      = (klayhgt/100)*s9->thick;
            /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

            /* calculate coordinates */
            for (l=0; l<3; l++)
              x[l] += zeta*hl*a3ref_l[l];
          }

          /* write coordinates of first kinematic layer at bottom */
          *dst_ptr++ = x[0];
          *dst_ptr++ = x[1];
          *dst_ptr++ = x[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            x_u[l] = x[l];

          /* calculate the coordinates of the top surface of the each
             kin layer */
          for (klay=0; klay<s9->num_klay; klay++) {
            INT num_mlay;
            DOUBLE* mlayhgt;
            DOUBLE sum_hgt;
            DOUBLE sum_hgt_mid;
            INT mlay;

            /* initialize the coordinates to values of middle surface */
            for (l=0; l<3; l++)
              x[l] = actnode->x[l];

            /* continuity matrix for aktual kin layer at top */
            for (jlay=0; jlay<s9->num_klay; jlay++) {
              DOUBLE klayhgt;
              DOUBLE hl;
              DOUBLE e3;
              DOUBLE zeta;

              klayhgt = s9->klayhgt[jlay];
              hl      = (klayhgt/100)*s9->thick;
              /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
              hl      = A3FAC_SHELL9*hl;   /*norm of director is half of kinematic layer hight*/
              e3 = +1.0; /*top*/
              zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

              /* calculate coordinates */
              for (l=0; l<3; l++)
                x[l] += zeta*hl*a3ref_l[l];
            }

            /* calculate the coordinates of the material layers of
             * aktual kinematic layer */
            num_mlay = actele->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save coordinates for interpolation for the material layers */
            for (l=0; l<3; l++)
              x_o[l] = x[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt     = 0.0; /*initialize the total hgt*/
            sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
            for (mlay=0; mlay<num_mlay; mlay++) {

              /* if more than one material layer to this kinematic layer */
              if (mlay > 0) {
                sum_hgt += mlayhgt[mlay-1];
                for (l=0; l<3; l++)
                  x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);

                /* write the bottom coordinates of this material layer */
                *dst_ptr++ = x[0];
                *dst_ptr++ = x[1];
                *dst_ptr++ = x[2];
              }

              /* a midnode has to be calculated */
              sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
              for (l=0; l<3; l++)
                x[l] = x_u[l] + (sum_hgt_mid/100.)*(x_o[l]-x_u[l]);

              /* write the middle coordinates of this material layer */
              *dst_ptr++ = x[0];
              *dst_ptr++ = x[1];
              *dst_ptr++ = x[2];
            }

            *dst_ptr++ = x_o[0];
            *dst_ptr++ = x_o[1];
            *dst_ptr++ = x_o[2];

            for (l=0; l<3; l++)
              x_u[l] = x_o[l];
          }
        }

        break;
      default:
        dserror("unknown minor type %d for shell9 element", minor);
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the displacements of all the layers of a shell9 element.

  Normally the displacements belong to nodes. But these are internal
  element displacements that might be put to the nodes of a

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell9_displacement(BIN_OUT_CHUNK *chunk,
                                         INT place,
                                         PARTDISCRET *actpdis,
                                         DOUBLE *send_buf,
                                         INT send_count,
                                         INT *send_size_buf,
                                         INT dst_first_id,
                                         INT dst_num)
{
  INT i, counter, len;

#ifdef DEBUG
  dstrc_enter("out_pack_shell9_displacement");
#endif

  len = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      SHELL9* s9;
      INT minor;
      DOUBLE *dst_ptr;
      DOUBLE dis[3], dis_u[3], dis_o[3]; /* displacement x,y,z */
      INT j;
      INT k;
      INT l;
      INT klay;

      s9 = actele->e.s9;
      minor = FIND_MINOR(actele);

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      /* this is ``s9_out_gid_sol_dis_unsmo`` in disguise */

      /* The purpose is to calculate all layer nodes'
       * displacements. That is each node reappears in each layer and
       * needs three displacement values. The difference between the
       * four node versions and the others is that with the four node
       * versions no nodes in the middle of a layer are generated.
       *
       * Please note that we're talking about node values and yet we
       * store element data. Each element has its inner nodes. It's up
       * to the filter to bring inner nodes of consecutive elements
       * together (smoothed view) or not (unconnected view). */

      switch (minor) {
      case MINOR_SHELL9_4_22:   /* 4-noded shell9 2x2 GP */
      case MINOR_SHELL9_4_33:   /* 4-noded shell9 3x3 GP */

        for (k=0; k<actele->numnp; k++) {
          NODE *actnode;

          actnode = actele->node[k];

          /* calculate the displacements of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the displacements of this layer to zero */
          for (l=0; l<3; l++)
            dis[l]=0.0;

          /* calculate the relative displacements of this layer */
          for (j=0; j<s9->num_klay; j++) {
            DOUBLE e3 = -1.0; /*bottom*/
            DOUBLE zeta = s9con(e3,s9->num_klay,klay,j,1.0);
            for (l=0; l<3; l++)
              dis[l] += zeta*actnode->sol.a.da[place][l+3*(j+1)];
          }

          /* add the displacement of reference layer to relative
           * displacement */
          for (l=0; l<3; l++)
            dis[l] += actnode->sol.a.da[place][l];

          /* write the displacements of first kinematic layer at bottom */
          *dst_ptr++ = dis[0];
          *dst_ptr++ = dis[1];
          *dst_ptr++ = dis[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            dis_u[l] = dis[l];

          /* calculate the displacements of the top surface of the
           * each kin layer */
          for (klay=0; klay<s9->num_klay; klay++) {
            DOUBLE sum_hgt;
            INT mlay;
            INT num_mlay;
            DOUBLE* mlayhgt;

            /* initialize the displacements of this layer to zero */
            for (l=0; l<3; l++)
              dis[l]=0.0;

            /* calculate the relative displacements of this layer */
            for (j=0; j<s9->num_klay; j++) {
              DOUBLE e3 = +1.0; /*top*/
              DOUBLE zeta = s9con(e3,s9->num_klay,klay,j,1.0);
              for (l=0; l<3; l++)
                dis[l] += zeta*actnode->sol.a.da[place][l+3*(j+1)];
            }

            /* add the displacement of reference layer to relative
             * displacement */
            for (l=0; l<3; l++)
              dis[l] += actnode->sol.a.da[place][l];

            /* calculate the displacements of the material layers of
             * aktual kinematic layer */
            num_mlay = s9->kinlay[klay].num_mlay;
            mlayhgt  = s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save displacements for interpolation for the material layers */
            for (l=0; l<3; l++)
              dis_o[l] = dis[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */

            /* initialize the total hgt */
            sum_hgt = 0.0;

            /* only if more than one material layer to this kinematic layer */
            for (mlay=1; mlay<num_mlay; mlay++) {
              sum_hgt += mlayhgt[mlay-1];
              for (l=0; l<3; l++)
                dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);

              /* write the bottom displacement of this material layer */
              *dst_ptr++ = dis[0];
              *dst_ptr++ = dis[1];
              *dst_ptr++ = dis[2];
            }

            /* write the top displacement of this kinematic layer */
            *dst_ptr++ = dis_o[0];
            *dst_ptr++ = dis_o[1];
            *dst_ptr++ = dis_o[2];

            for (l=0; l<3; l++)
              dis_u[l] = dis_o[l];
          }
        }
        break;
      case MINOR_SHELL9_8_22:   /* 8-noded shell9 2x2 GP */
      case MINOR_SHELL9_8_33:   /* 8-noded shell9 3x3 GP */
      case MINOR_SHELL9_9_22:   /* 9-noded shell9 2x2 GP */
      case MINOR_SHELL9_9_33:   /* 9-noded shell9 3x3 GP */

        for (k=0; k<actele->numnp; k++) {
          NODE *actnode;
          INT jlay;

          actnode = actele->node[k];

          /* calculate the displacements of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the displacements of this layer to zero */
          for (l=0; l<3; l++)
            dis[l]=0.0;

          /* calculate the relative displacements of this layer */
          for (jlay=0; jlay<s9->num_klay; jlay++) {
            DOUBLE e3   = -1.0; /*bottom*/
            DOUBLE zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);
            for (l=0; l<3; l++)
              dis[l] += zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
          }

          /* add the displacement of reference layer to relative
           * displacement */
          for (l=0; l<3; l++)
            dis[l] += actnode->sol.a.da[place][l];

          /* write the displacements of first kinematic layer at bottom */
          *dst_ptr++ = dis[0];
          *dst_ptr++ = dis[1];
          *dst_ptr++ = dis[2];

          /* save displacements for interpolation for the material layers */
          for (l=0; l<3; l++)
            dis_u[l] = dis[l];

          /* calculate the displacements of the top surface of the
           * each kin layer */
          for (klay=0; klay<s9->num_klay; klay++) {
            DOUBLE sum_hgt;
            DOUBLE sum_hgt_mid;
            INT mlay;
            INT num_mlay;
            DOUBLE* mlayhgt;

            /* initialize the displacements of this layer to zero */
            for(l=0; l<3; l++)
              dis[l]=0.0;

            /* calculate the relative displacements of this layer */
            for (jlay=0; jlay<s9->num_klay; jlay++) {
              DOUBLE e3   = +1.0; /*top*/
              DOUBLE zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);
              for (l=0; l<3; l++)
                dis[l] += zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
            }

            /* add the displacement of reference layer to relative
             * displacement */
            for (l=0; l<3; l++)
              dis[l] += actnode->sol.a.da[place][l];

            /* calculate the displacements of the material layers of
             * aktual kinematic layer */
            num_mlay = actele->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save displacements for interpolation for the material layers */
            for (l=0; l<3; l++)
              dis_o[l] = dis[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt     = 0.0; /*initialize the total hgt*/
            sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/

            for (mlay=0; mlay<num_mlay; mlay++) {

              /* if more than one material layer to this kinematic layer */
              if (mlay > 0) {
                 sum_hgt += mlayhgt[mlay-1];
                 for (l=0; l<3; l++)
                   dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);

                 /* write the bottom displacement of this material layer */
                 *dst_ptr++ = dis[0];
                 *dst_ptr++ = dis[1];
                 *dst_ptr++ = dis[2];
              }

              /* the displacement at midnode has to be calculated */
              sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
              for (l=0; l<3; l++)
                dis[l] = dis_u[l] + (sum_hgt_mid/100.)*(dis_o[l]-dis_u[l]);

              /* write the mid displacement of this material layer  */
              *dst_ptr++ = dis[0];
              *dst_ptr++ = dis[1];
              *dst_ptr++ = dis[2];
            }

            /* write the top displacement of this kinematic layer */
            *dst_ptr++ = dis_o[0];
            *dst_ptr++ = dis_o[1];
            *dst_ptr++ = dis_o[2];

            for (l=0; l<3; l++)
              dis_u[l] = dis_o[l];
          }
        }

        break;
      default:
          dserror("unknown minor type %d for shell9 element", minor);
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the velocities in some buffer to send them to their
  writing processor.

  There are two or three velocity values per node depending on the
  number of dimensions.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_velocity(BIN_OUT_CHUNK *chunk,
                              INT place,
                              PARTDISCRET *actpdis,
                              DOUBLE *send_buf,
                              INT send_count,
                              INT *send_size_buf,
                              INT dst_first_id,
                              INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_velocity");
#endif

  ndim = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num)) {
      INT k;
      DOUBLE *ptr = actnode->sol.a.da[place];

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim >= ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      for (k=0; k<ndim; ++k) {
        send_buf[ndim*counter+k] = *ptr++;
      }
      counter += 1;
    }
  }
  dsassert(counter*ndim == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the pressure in some buffer to send them to their
  writing processor.

  The point is that there is just one pressure value per node. It
  lifes in the node's sol array in the same row as the velocity
  values. That's why we need to get the current problem's number of
  dimensions.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_pressure(BIN_OUT_CHUNK *chunk,
                              INT place,
                              PARTDISCRET *actpdis,
                              DOUBLE *send_buf,
                              INT send_count,
                              INT *send_size_buf,
                              INT dst_first_id,
                              INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_pressure");
#endif

  ndim = genprob.ndim;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num)) {

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim > ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      send_buf[counter] = actnode->sol.a.da[place][ndim];

      counter += 1;
    }
  }
  dsassert(counter == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the element's stresses.

  Here the element based stresses are handled.

  This is highly element specific. That is the filter must know what
  these numbers mean. Each element is different.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_stress(BIN_OUT_CHUNK *chunk,
                            INT place,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT j;
  INT k;
  INT len;
  INT counter;

/*
   gausspoint permutation :
   On the Gausspoint number i in Gid, the results of Carats GP number gausspermn[i]
   have to be written
*/

  static INT gaussperm4[4] = {3,1,0,2};
  /*INT gaussperm8[8] = {0,4,2,6,1,5,3,7};*/
  static INT gaussperm9[9] = {8,2,0,6,5,1,3,7,4};
  /*INT gaussperm27[27] = {0,9,18,3,12,21,6,15,24,1,10,19,4,13,22,7,16,25,2,11,20,5,14,23,8,17,26};*/

#ifdef DEBUG
  dstrc_enter("out_pack_stress");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
/*   dsassert(chunk->value_entry_length == 0, "invalid value entry length"); */

  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      INT minor;
      DOUBLE **stress;
      DOUBLE *dst_ptr;
      INT rows;
      INT cols;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      minor = FIND_MINOR(actele);

      /* What needs to be done depends on major and minor element
       * number. */

      switch (actele->eltyp) {
#ifdef D_WALL1
      case el_wall1:            /* 2D plane stress - plane strain element */

        /* there is just one source for all variants of wall elements */
        stress = actele->e.w1->stress_GP.a.d3[place];

        /*
         * The triangle version needs special treatment. And no
         * thought has been spend about the D_SSI issue. */

        switch (minor) {
        case MINOR_WALL1_11:    /* 3-noded wall1 1x1 GP */
          rows = 1;
          cols = 4;
          for (k=0; k<cols; ++k) {
            *dst_ptr++ = stress[k][0];
          }
          break;
        case MINOR_WALL1_22:    /* 4-noded wall1 2x2 GP */
        case MINOR_WALL1_8_33:  /* 8-noded wall1 3x3 GP */
        case MINOR_WALL1_9_33:  /* 9-noded wall1 3x3 GP */
        case MINOR_WALL1_8_22:  /* 8-noded wall1 2x2 GP */
          rows = element_info[el_wall1].variant[minor].gauss_number;
          cols = 6;
          if (actele->e.w1->elewa != NULL) {
            for (j=0; j<rows; j++) {
              *dst_ptr++ = stress[0][gaussperm4[j]];
              *dst_ptr++ = stress[1][gaussperm4[j]];
              *dst_ptr++ = stress[2][gaussperm4[j]];
              *dst_ptr++ = stress[3][gaussperm4[j]];
              *dst_ptr++ = actele->e.w1->elewa[0].ipwa[gaussperm4[j]].damage;
              *dst_ptr++ = actele->e.w1->elewa[0].ipwa[gaussperm4[j]].aequistrain;
            }
          }
          else {
            for (j=0; j<rows; j++) {
              *dst_ptr++ = stress[0][gaussperm4[j]];
              *dst_ptr++ = stress[1][gaussperm4[j]];
              *dst_ptr++ = stress[2][gaussperm4[j]];
              *dst_ptr++ = stress[3][gaussperm4[j]];
              *dst_ptr++ = 0;
              *dst_ptr++ = 0;
            }
          }
          break;
        default:
          dserror("unknown minor type %d for wall element", minor);
        }
        break;
#endif
#ifdef D_BEAM3
      case el_beam3:            /* structural 3D-beam element */

        stress = actele->e.b3->force_GP.a.d3[place];
        cols = 6;

        switch (minor) {
        case MINOR_BEAM3_21:    /* 2-noded beam3 1 GP */
          rows = 1;
          break;
        case MINOR_BEAM3_22:    /* 2-noded beam3 2 GP */
        case MINOR_BEAM3_32:    /* 3-noded beam3 2 GP */
          rows = 2;
          break;
        case MINOR_BEAM3_33:    /* 3-noded beam3 3 GP */
          rows = 3;
          break;
        default:
          dserror("unknown minor type %d for beam3 element", minor);
        }

        /* copy the values */
        /* It's unusual that the later index varies less than the
         * first. However that's how the stresses are stored and I
         * don't want to deviate from the way the other elements do it. */
        for (j=0; j<rows; j++) {
          for (k=0; k<cols; ++k) {
            *dst_ptr++ = stress[k][j];
          }
        }

        break;
#endif
#ifdef D_INTERF
      case el_interf:           /* 1D interface element (combination only with wall) */

        /*
         * The first element on proc 0 decides what kind of stresses
         * we get. Of course it's undefined in a parallel setting
         * which elements will be on proc 0. Thus this is only valid
         * if all interf elements are the same. But this flag just
         * determines how these stresses are labeled anyway.
         *
         * What do we do when there happen to be no interf elements on
         * proc 0? */
        if (map_symbol_count(&chunk->group, "interf_orient") == 0) {
          if (actele->e.interf->stresstyp == if_tn) {
            map_insert_string_cpy(&chunk->group, "interf_orient", "local");
          }
          else {
            map_insert_string_cpy(&chunk->group, "interf_orient", "global");
          }
        }

        stress = actele->e.interf->stress_GP.a.d3[place];
        cols = 5;
        switch (minor) {
        case MINOR_INTERF_22:   /* interface 2x2 GP */
          rows = 2;
          break;
        case MINOR_INTERF_33:   /* interface 3x3 GP */
          rows = 3;
          break;
        default:
          dserror("unknown minor type %d for interface element", minor);
        }

        /* copy the values */
        for (j=0; j<rows; j++) {
          *dst_ptr++ = stress[0][j];
          *dst_ptr++ = stress[1][j];
          *dst_ptr++ = stress[2][j];
          *dst_ptr++ = actele->e.interf->elewa[0].ipwa[j].dn;
          *dst_ptr++ = actele->e.interf->elewa[0].ipwa[j].dt;
        }

        break;
#endif
#ifdef D_WALLGE
      case el_wallge:           /* gradient enhanced wall element */
        /*stress = actele->e.wallge->stress_GP.a.d3[place];*/

        switch (minor) {
        case MINOR_WALLGE_22:   /* gradient enhanced wall 2x2 GP */
          rows = 4;
          break;
        case MINOR_WALLGE_8_33: /* gradient enhanced wall 3x3 GP */
        case MINOR_WALLGE_9_33: /* gradient enhanced wall 3x3 GP */
          rows = 9;
          break;
        default:
          dserror("unknown minor type %d for wallge element", minor);
        }

        /* copy the values */
        for (j=0; j<rows; j++) {
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[0];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[1];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].sig[2];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].damage;
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].aequistrain;
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[gaussperm4[j]].aequistrain_nl;
        }

        break;
#endif
#ifdef D_AXISHELL
      case el_axishell:         /* 1D axisymmetrical shell element */
        /*
         * There were additional (?) displacements. But these are
         * nodal values. They are written at another place. */

        stress = actele->e.saxi->stress_GP.a.d3[place];
        rows = 5;

        for (j=0; j<rows; j++) {
          *dst_ptr++ = stress[j][0];
        }
        break;
#endif
#ifdef D_FLUID3
      case el_fluid3:           /* 3D fluid element */
        /* node based. ignore here. */
        break;
#endif
#ifdef D_SHELL9
      case el_shell9: {         /* multi layer shell element */
        SHELL9* s9 = actele->e.s9;

        stress = s9->gp_stress.a.da;

        /* In general send_count can be bigger that this element's
         * contribution. In a homogeneous mesh these sizes must be
         * equal, of course, but lets not rely on the homogeneousness
         * when we don't need to. */
        rows = s9->gp_stress.fdim;
        cols = s9->gp_stress.sdim;

        dsassert(rows*cols <= send_count, "shell9 entry too small");

        /* copy the values */
        for (j=0; j<rows; j++) {
          for (k=0; k<cols; k++) {
            *dst_ptr++ = stress[j][k];
          }
        }
        break;
      }
#endif
      case el_shell8:           /* 7 parameter shell element */
      case el_brick1:           /* structural brick element */
      case el_fluid2:           /* 2D fluid element */
      case el_fluid2_pro:       /* 2D fluid element */
      case el_fluid2_tu:        /* 2D fluid element for turbulence */
      case el_ale2:             /* 2D pseudo structural ale element */
      case el_ale3:             /* 3D pseudo structural ale element */
      case el_ls2:              /* 2D element for level set calculations */
      case el_fluid2_xfem:      /* 2D element for x-fem calculations */
      default:
        dserror("element based stress output not supported for element type %d", actele->eltyp);
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the stresses.

  Here the node based stresses are handled. I'd like to avoid
  this. Better have the filter do the interpolation.

  This is highly element specific. That is the filter must know what
  these numbers mean. Each element is different.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_nd_stress(BIN_OUT_CHUNK *chunk,
                               INT place,
                               PARTDISCRET *actpdis,
                               DOUBLE *send_buf,
                               INT send_count,
                               INT *send_size_buf,
                               INT dst_first_id,
                               INT dst_num)
{
  INT len;
  INT counter;
  INT j;

#ifdef DEBUG
  dstrc_enter("out_pack_nd_stress");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  len = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num)) {
      INT k;
      INT count;
      INT numele;
      ELEMENT *actele;
      INT minor;
      DOUBLE invcount;
      DOUBLE *dst_ptr;
      DOUBLE **stress;

      dst_ptr = &(send_buf[len*counter]);

      numele = actnode->numele;

      /* The approach demands that there is just one type of element
       * in the discretization. So we take the first element and
       * decide depending on its type what to do. */
      actele = actnode->element[0];

      minor = FIND_MINOR(actele);

      switch (actele->eltyp) {
#ifdef D_FLUID3
      case el_fluid3:           /* 3D fluid element */
        dsassert(minor == MINOR_FLUID3_222, "unsupported f3 variant");
        dsassert(len >= 12, "stress entry too small");

        /* This is the function f3_out_gid_sol_str in disguise. */
        count = 0;
        for (k=0; k<12; ++k) {
          dst_ptr[k] = 0.0;
        }
        for (k=0; k<numele; ++k) {
          INT i;
          INT l;
          actele = actnode->element[j];

          if (actele->eltyp != el_fluid3 || actele->numnp !=8) {
            dserror("uniform mesh expected");
          }

          count++;
          stress = actele->e.f3->stress_ND.a.da;
          for (l=0; l<8; l++)
            if (actele->node[l] == actnode) break;
          for (i=0; i<12; i++) {
            dst_ptr[i] += stress[l][i];
          }
        }

        invcount = 1.0/count;
        for (j=0; j<12; j++) {
          dst_ptr[j] *= invcount;
        }

        break;
#endif

      default:
        dserror("node based stress output not supported for element type %d", actele->eltyp);
      }


#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

/*       for (k=0; k<ndim; ++k) { */
/*         send_buf[ndim*counter+k] = *ptr++; */
/*       } */
      counter += 1;
    }
  }
  dsassert(counter*len == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the element's processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_domain(BIN_OUT_CHUNK *chunk,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_domain");
#endif

  dsassert(chunk->size_entry_length == 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  counter = 0;
  for (i=0; i<actpdis->numele; ++i) {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {

#ifdef PARALLEL
      send_size_buf[2*counter  ] = actele->Id_loc;
      send_size_buf[2*counter+1] = actele->proc;
#else
      send_size_buf[counter] = actele->proc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack an array of distributed vectors.

  The distributed vector is supposed to contain one value per dof. We
  store the dofs in ascending order, the order in the distributed
  vector at hand, however, is solver dependent.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dist_vector(BIN_OUT_CHUNK *chunk,
                                 PARTDISCRET *actpdis,
                                 DOUBLE *send_buf,
                                 INT send_count,
                                 INT *send_size_buf,
                                 INT dst_first_id,
                                 INT dst_num)
{
  INT i;
  INT len;
  INT counter;
  SPARSE_TYP sysarray_typ;
  SPARSE_ARRAY sysarray;

#ifdef DEBUG
  dstrc_enter("out_pack_dist_vector");
#endif

  /* for convenience */
  sysarray_typ = *chunk->field->sysarray_typ;
  sysarray     = *chunk->field->sysarray;

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  len = chunk->value_entry_length;

#ifdef PARALLEL
#define pack_send_size_buf send_size_buf[counter] = dof;
#else
#define pack_send_size_buf
#endif

#define boilerplate_copying_code                        \
  if ((dof >= dst_first_id) &&                          \
      (dof < dst_first_id+dst_num)) {                   \
    INT j;                                              \
    DOUBLE* dst_ptr;                                    \
    dst_ptr = &(send_buf[len*counter]);                 \
    for (j=0; j<len; ++j) {                             \
      *dst_ptr++ = chunk->vectors[j].vec.a.dv[i];       \
    }                                                   \
    pack_send_size_buf;                                 \
    counter += 1;                                       \
  }

  counter = 0;

  /* Here we have the solver dependency. This is
   * ``solserv_reddistvec`` in disguise. */
  switch (sysarray_typ) {

#ifdef AZTEC_PACKAGE
  case msr:
    for (i=0; i<sysarray.msr->numeq; ++i) {
      INT dof = sysarray.msr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef HYPRE_PACKAGE
  case parcsr: {
    INT rank = chunk->field->actintra->intra_rank;
    for (i=0; i<sysarray.parcsr->numeq; ++i) {
      INT dof = sysarray.parcsr->update.a.ia[rank][i];
      boilerplate_copying_code;
    }
    break;
  }
#endif

#ifdef PARSUPERLU_PACKAGE
  case ucchb:
    for (i=0; i<sysarray.ucchb->numeq; ++i) {
      INT dof = sysarray.ucchb->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case dense:
    for (i=0; i<sysarray.dense->numeq; ++i) {
      INT dof = sysarray.dense->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

#ifdef MLIB_PACKAGE
  case mds:
    for (i=0; i<sysarray.mds->numeq; ++i) {
      INT dof = i;
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef MUMPS_PACKAGE
  case rc_ptr:
    for (i=0; i<sysarray.rc_ptr->numeq; ++i) {
      INT dof = sysarray.rc_ptr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef SPOOLES_PACKAGE
  case spoolmatrix:
    for (i=0; i<sysarray.spo->numeq; ++i) {
      INT dof = sysarray.spo->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef UMFPACK
  case ccf:
    for (i=0; i<sysarray.ccf->numeq; ++i) {
      INT dof = sysarray.ccf->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case skymatrix:
    for (i=0; i<sysarray.sky->numeq; ++i) {
      INT dof = sysarray.sky->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

#ifdef MLPCG
  case bdcsr:
    for (i=0; i<sysarray.bdcsr->numeq; ++i) {
      INT dof = sysarray.bdcsr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case oll:
    for (i=0; i<sysarray.oll->numeq; ++i) {
      INT dof = sysarray.oll->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

  default:
    dserror("Unknown type %d of system matrix", sysarray_typ);
    break;
  }

#undef boilerplate_copying_code

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the restart data of all the elements in this
  discretization.

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_restart_element(BIN_OUT_CHUNK *chunk,
                                     PARTDISCRET *actpdis,
                                     DOUBLE *send_buf,
                                     INT send_count,
                                     INT *send_size_buf,
                                     INT dst_first_id,
                                     INT dst_num)
{
  INT el;
  INT len;
  INT slen;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_restart_element");
#endif

  len = chunk->value_entry_length;
  slen = chunk->size_entry_length;

  counter = 0;
  for (el=0; el<actpdis->numele; ++el) {
    ELEMENT* actele = actpdis->element[el];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num)) {
      DOUBLE* dst_ptr;
      INT* size_dst_ptr;
      INT j;
      INT minor;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actele->Id_loc;
      size_dst_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_dst_ptr = &(send_size_buf[slen*counter]);
#endif

      dst_ptr = &(send_buf[len*counter]);

      /* Now here we have all the different elements. */
      switch (actele->eltyp) {

#ifdef D_SHELL8
      case el_shell8: {
        INT arr_length;
        SHELL8* s8;
        DOUBLE* src_ptr;

        s8 = actele->e.s8;
        if (s8->nhyb) {

          arr_length = s8->alfa.fdim * s8->alfa.sdim;
          src_ptr = s8->alfa.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Dtildinv.fdim * s8->Dtildinv.sdim;
          src_ptr = s8->Dtildinv.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Lt.fdim * s8->Lt.sdim;
          src_ptr = s8->Lt.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Rtilde.fdim * s8->Rtilde.sdim;
          src_ptr = s8->Rtilde.a.dv;
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }
        }
        if (mat[actele->mat-1].mattyp==m_viscohyper) {
          ARRAY4D *a;
          a = s8->his1;
          arr_length = a->fdim * a->sdim * a->tdim * a->fodim;
          src_ptr = a->a.d4[0][0][0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }
        }

        break;
      }
#endif

#ifdef D_SHELL9
      case el_shell9: {
        INT arr_length;
        SHELL9* s9;
        DOUBLE* src_ptr;

        s9 = actele->e.s9;
        if (s9->nhyb) {

          arr_length = s9->alfa.fdim * s9->alfa.sdim;
          src_ptr = s9->alfa.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Dtildinv.fdim * s9->Dtildinv.sdim;
          src_ptr = s9->Dtildinv.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->L.fdim * s9->L.sdim;
          src_ptr = s9->L.a.da[0];
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Rtilde.fdim * s9->Rtilde.sdim;
          src_ptr = s9->Rtilde.a.dv;
          for (j=0; j<arr_length; ++j) {
            *dst_ptr++ = *src_ptr++;
          }
        }

        if (s9->elewa != NULL) {
          INT kl;

          for (kl=0; kl<s9->num_klay; kl++) {
            INT num_mlay;
            INT ml;
            INT actlay = 0;

            num_mlay = s9->kinlay[kl].num_mlay;
            for (ml=0; ml<num_mlay; ml++) {
              S9_IP_WA *ipwa = s9->elewa[actlay].ipwa;

              /* check if there is a ipwa for this layer */
              if (ipwa != NULL) {
                INT k;
                INT ngauss;
                MULTIMAT *actmultimat;

                actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

                /* number of gausspoints in one layer*/
                ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

                switch (actmultimat->mattyp) {
                case m_pl_mises:
                case m_pl_dp:
                  for (k=0; k<ngauss; k++) {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].epstn;
                    for (j=0; j<6; j++) {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                      *dst_ptr++ = ipwa[k].qn[j];
                    }
                  }
                  break;
                case m_pl_epc:
                  for (k=0; k<ngauss; k++) {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].kappa_t;
                    *dst_ptr++ = ipwa[k].kappa_c;
                    for (j=0; j<6; j++) {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                    }
                  }
                  break;
                case m_pl_hoff:
                  for (k=0; k<ngauss; k++) {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].dhard;
                    for (j=0; j<6; j++) {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                      *dst_ptr++ = ipwa[k].dkappa[j];
                      *dst_ptr++ = ipwa[k].gamma[j];
                    }
                    for (j=0; j<9; j++) {
                      *dst_ptr++ = ipwa[k].rkappa[j];
                    }
                  }
                  break;
                default:
                  printf("unsupported material %d for shell9 element",
                         actmultimat->mattyp);
                }
              }
              actlay += 1;
            }
          }
        }

        break;
      }
#endif

#ifdef D_WALL1
      case el_wall1: {
        W1_ELE_WA* elewa = actele->e.w1->elewa;
        if (elewa != NULL) {
          INT k;
          INT ngauss;

          minor = MINOR_WALL1(actele);
          ngauss = element_info[el_wall1].variant[minor].gauss_number;

          /* This is w1_write_restart. Disguised as usual. */
          switch (mat[actele->mat-1].mattyp) {
          case m_pl_mises:
            dsassert(len >= ngauss*(1+4*3), "value entry too short");
            for (k=0; k<ngauss; ++k) {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].epstn;
              for (j=0; j<4; j++) {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].qn[j];
              }
            }
            break;
          case m_pl_mises_3D:
            dsassert(len >= ngauss*(1+4*7), "value entry too short");
            for (k=0; k<ngauss; k++) {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].epstn;
              for (j=0; j<4; j++) {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].qn[j];
                *dst_ptr++ = elewa->ipwa[k].sigc[j];
                *dst_ptr++ = elewa->ipwa[k].sigi[j];
                *dst_ptr++ = elewa->ipwa[k].epsi[j];
                *dst_ptr++ = elewa->ipwa[k].di[j];
              }
            }
            break;
          case m_pl_epc3D: /* ???
                            * (mat[actele->mat-1].mattyp == m_pl_epc3D ) */
            dsassert(len >= ngauss*(2+4*7), "value entry too short");
            for (k=0; k<ngauss; k++) {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].dlam[0];
              *dst_ptr++ = elewa->ipwa[k].dlam[1];
              for (j=0; j<4; j++) {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].sigc[j];
                *dst_ptr++ = elewa->ipwa[k].grad[j];
                *dst_ptr++ = elewa->ipwa[k].sigi[j];
                *dst_ptr++ = elewa->ipwa[k].epsi[j];
                *dst_ptr++ = elewa->ipwa[k].di[  j];
              }
            }
            break;
          default:
            printf("WARNING: Unknown mattyp %d in wall1 restart\n",
                   mat[actele->mat-1].mattyp);
          }
        }
        break;
      }
#endif

#ifdef D_BEAM3
      case el_beam3:
        if (mat[actele->mat-1].mattyp == m_pl_mises ||
            mat[actele->mat-1].mattyp == m_pl_dp ||
            mat[actele->mat-1].mattyp == m_pl_epc) {
          INT j;
          INT size;
          DOUBLE* src_ptr;

          /* Copy the values from the element's array to the chunk. */
          src_ptr = actele->e.b3->elewa.a.da[0];
          size = actele->e.b3->elewa.fdim*actele->e.b3->elewa.sdim;
          dsassert(size <= len, "value entry too small");
          for (j=0; j<size; ++j) {
            *dst_ptr++ = *src_ptr++;
          }
        }
        break;
#endif

        /* There's nothing to be saved for interface and wallge
         * elements, right? */
      case el_interf:
        break;
      case el_wallge:
        break;

#ifdef D_FLUID2
      case el_fluid2: {
        FLUID_DYNAMIC *fdyn;
        DOUBLE* src_ptr;

        /* only valid for single field problem (?) */
        /* I don't think so. */
        fdyn = alldyn[genprob.numff].fdyn;

        dsassert(len >= 3, "value item too small");

        src_ptr = &(actele->e.f2->tau_old.a.dv[0]);
        for (j=0; j<3; ++j) {
          *dst_ptr++ = *src_ptr++;
        }
        if (fdyn->surftens > 0) {
          dsassert(len >= 3+2*actele->numnp, "value item too small");
          src_ptr = &(actele->e.f2->kappa_ND.a.dv[0]);
          if (actele->e.f2->fs_on>0) {
            for (j=0; j<2*actele->numnp; ++j) {
              *dst_ptr++ = *src_ptr++;
            }
          }
        }
        break;
      }
#endif

#ifdef D_FLUID3
      case el_fluid3: {
        DOUBLE* src_ptr;

        dsassert(len >= 3, "value item too small");
        src_ptr = &(actele->e.f3->tau_old.a.dv[0]);
        for (j=0; j<3; ++j) {
          *dst_ptr++ = *src_ptr++;
        }
        break;
      }
#endif

      default: {
        static CHAR warning[el_count];
        if (!warning[actele->eltyp]) {
          warning[actele->eltyp] = 1;
          printf(RED_LIGHT
                 "restart for element type '"
                 GREEN_LIGHT
                 "%s"
                 RED_LIGHT
                 "' not supported"
                 END_COLOR
                 "\n",
                 element_info[actele->eltyp].name);
        }
      }
      }

      counter += 1;
    }
  }

  dsassert(counter*len == send_count, "element pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the values to be saved in some buffer to send them to
  their writing processor.

  This is the hook function called by ``out_gather_values`` in order
  to collect certain values from nodes, elements or distributed
  vectors. This function or the ones called here must be changed in
  case an element is updated.

  There is a loop that calls this function for each destination
  processor in turn. That is upon one call it has to collect all items
  that go to one processor. These are always consecutive items, \a
  dst_first_id gives the first item's id, \a dst_num gives the number
  of items to be collected.

  \param *chunk          (i)  the chunk that's going to be written
  \param  type           (i)  what kind of values are to be collected
  \param  array          (i)  the row of the node array that
                              interesting; is ignored for some types
  \param *actpdis        (i)  the partition's discretization
  \param *send_buf       (o)  buffer to collect double values
  \param  send_count     (i)  number of double values to be collected
  \param *send_size_buf  (o)  buffer to collect integer values
  \param  dst_first_id   (i)  the id (Id_loc) of the first item to be written
  \param  dst_num        (i)  the number of (consecutive) items to be written

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_pack_items(struct _BIN_OUT_CHUNK *chunk,
                    CHUNK_CONTENT_TYPE type,
                    INT array,
                    PARTDISCRET *actpdis,
                    DOUBLE *send_buf,
                    INT send_count,
                    INT *send_size_buf,
                    INT dst_first_id,
                    INT dst_num)
{
#ifdef DEBUG
  dstrc_enter("out_pack_items");
#endif

  switch (type) {
  case cc_node_array:
    out_pack_node_arrays(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_mesh:
    /* array flag unused here */
    out_pack_mesh(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_coords:
    /* array flag unused here */
    out_pack_coords(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_displacement:
    out_pack_displacement(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_velocity:
    out_pack_velocity(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_pressure:
    out_pack_pressure(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_el_stress:
    out_pack_stress(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_nd_stress:
    out_pack_nd_stress(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_domain:
    /* array flag unused here */
    out_pack_domain(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_dist_vector:
    /* array flag unused here */
    out_pack_dist_vector(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#ifdef D_SHELL8
  case cc_shell8_director:
    /* array flag unused here */
    out_pack_shell8_director(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#endif
#ifdef D_SHELL9
  case cc_shell9_coords:
    /* array flag unused here */
    out_pack_shell9_coords(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_shell9_displacement:
    out_pack_shell9_displacement(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#endif
  case cc_restart_element:
    /* array flag unused here */
    out_pack_restart_element(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  default:
    dserror("unsupported chunk type %d", type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/* Here we have lots of functions that do nothing but unpacking input
 * buffers and put the values somewhere (node arrays, element stuff,
 * distributed vectors). These functions share a common structure. See
 * ``in_unpack_items`` for an explanation of their parameters.
 *
 * These are the functions that need to be changed according to
 * changes in the elements. */
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack the node arrays.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_node_arrays(BIN_IN_FIELD *context,
                                  BIN_IN_CHUNK *chunk,
                                  NODE_ARRAY array,
                                  DOUBLE *recv_buf,
                                  INT recv_count,
                                  INT *recv_size_buf,
                                  INT recv_size_count,
                                  INT src)
{
  INT j;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("in_unpack_node_arrays");
#endif

  /* There are only very small differences between the sequential
   * and the parallel version here. */

#ifdef PARALLEL
#define unpack_recv_size_buf                                            \
  dsassert(3*j+3 <= recv_size_count,                                    \
           "receive size buffer overrun");                              \
  dsassert(recv_size_buf[3*j] == actnode->Id_loc, "Id_loc mismatch");   \
  fdim = recv_size_buf[3*j+1];                                          \
  sdim = recv_size_buf[3*j+2];
#else
#define unpack_recv_size_buf                    \
  dsassert(2*j+2 <= recv_size_count,            \
           "receive size buffer overrun");      \
  fdim = recv_size_buf[2*j  ];                  \
  sdim = recv_size_buf[2*j+1];
#endif

#ifdef PARALLEL
#define get_id_part(j) (context->recv_node_ids[src][2*j+1])
#define get_recv_numnp (context->recv_numnp[src])
#else
#define get_id_part(j) (j)
#define get_recv_numnp (actpdis->numnp)
#endif

  dsassert(chunk->size_entry_length == 2, "size entry length mismatch");

  /* scatter received data to the nodes */
#define call_scatter_values(node_array)                                 \
  for (j=0; j<get_recv_numnp; ++j) {                                    \
    INT k;                                                              \
    INT Id_part;                                                        \
    NODE* actnode;                                                      \
    INT fdim;                                                           \
    INT sdim;                                                           \
    DOUBLE *ptr;                                                        \
    Id_part = get_id_part(j);                                           \
    dsassert((Id_part >= 0) && (Id_part < actpdis->numnp),              \
             "illegal node id");                                        \
    actnode = actpdis->node[Id_part];                                   \
    /* Id_part only in the more parallel version. */                    \
    /*dsassert(Id_part == actnode->Id_part, "Id_part mismatch");*/      \
                                                                        \
    unpack_recv_size_buf;                                               \
                                                                        \
    dsassert(sdim == actnode->node_array.sdim,                          \
             "node array sdim mismatch");                               \
    if (fdim > actnode->node_array.fdim) {                              \
      amredef(&(actnode->node_array),fdim,sdim,"DA");                   \
    }                                                                   \
    ptr = actnode->node_array.a.da[0];                                  \
    dsassert(chunk->value_entry_length*j+fdim*sdim <= recv_count,       \
             "receive buffer overrun");                                 \
    for (k=0; k<fdim*sdim; ++k) {                                       \
      *ptr++ = recv_buf[chunk->value_entry_length*j+k];                 \
    }                                                                   \
  }

  switch (array) {
  case node_array_sol:
    call_scatter_values(sol);
    break;
  case node_array_sol_increment:
    call_scatter_values(sol_increment);
    break;
  case node_array_sol_residual:
    call_scatter_values(sol_residual);
    break;
  case node_array_sol_mf:
    call_scatter_values(sol_mf);
    break;
  default:
    dserror("Node array %d unknown", array);
  }

#undef call_scatter_values
#undef get_id_part
#undef get_recv_numnp
#undef unpack_recv_size_buf

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack an array of distributed vectors.

  This is not solver dependent because we already know the locations
  where the dof values are to got to. The field initialization found
  this (solver dependent) information.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_dist_vector(BIN_IN_FIELD *context,
                                  BIN_IN_CHUNK *chunk,
                                  DOUBLE *recv_buf,
                                  INT recv_count,
                                  INT *recv_size_buf,
                                  INT recv_size_count,
                                  INT src)
{
  INT j;
  INT len;

#ifdef DEBUG
  dstrc_enter("in_unpack_dist_vector");
#endif

#ifdef PARALLEL
#define get_id_part(j)  (context->recv_dof_ids[src][2*j+1])
#define get_recv_numdof (context->recv_numdof[src])
#else
#define get_id_part(j)  (j)
#define get_recv_numdof (chunk->vectors[0].numeq)
#endif

  dsassert(chunk->size_entry_length == 0, "size entry length mismatch");

  len = chunk->value_entry_length;

  /* scatter received data to the nodes */
  for (j=0; j<get_recv_numdof; ++j) {
    INT k;
    INT Id_part;
    DOUBLE* src_ptr;

    Id_part = get_id_part(j);

    src_ptr = &(recv_buf[len*j]);
    for (k=0; k<len; ++k) {
      chunk->vectors[k].vec.a.dv[Id_part] = *src_ptr++;
    }
  }

#undef get_id_part
#undef get_recv_numele

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack what we received.

  This is called during restart, the elements are already set up. The
  arrays are there. Just the values are missing.

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_restart_element(BIN_IN_FIELD *context,
                                      BIN_IN_CHUNK *chunk,
                                      DOUBLE *recv_buf,
                                      INT recv_count,
                                      INT *recv_size_buf,
                                      INT recv_size_count,
                                      INT src)
{
  INT el;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
  INT len;
  INT slen;

#ifdef DEBUG
  dstrc_enter("in_unpack_restart_element");
#endif

#ifdef PARALLEL
#define get_id_part(j)  (context->recv_element_ids[src][2*j+1])
#define get_recv_numele (context->recv_numele[src])
#else
#define get_id_part(j)  (j)
#define get_recv_numele (actpdis->numele)
#endif

  len = chunk->value_entry_length;
  slen = chunk->size_entry_length;

  /* scatter received data to the nodes */
  for (el=0; el<get_recv_numele; ++el) {
    INT Id_part;
    ELEMENT* actele;
    DOUBLE* src_ptr;
    INT* size_src_ptr;

    Id_part = get_id_part(el);
    dsassert((Id_part >= 0) && (Id_part < actpdis->numele),
             "illegal element id");
    actele = actpdis->element[Id_part];
    /*dsassert(Id_part == actele->Id_part, "Id_part mismatch");*/

#ifdef PARALLEL
    dsassert(recv_size_buf[(slen+1)*el] == actele->Id_loc, "Id_loc mismatch");
    size_src_ptr = &(recv_size_buf[(slen+1)*el+1]);
#else
    size_src_ptr = &(recv_size_buf[slen*el]);
#endif

    src_ptr = &(recv_buf[len*el]);

    switch (actele->eltyp) {

#ifdef D_SHELL8
      case el_shell8: {
        INT arr_length;
        SHELL8* s8;
        DOUBLE* dst_ptr;

        s8 = actele->e.s8;
        if (s8->nhyb) {

          arr_length = s8->alfa.fdim * s8->alfa.sdim;
          dst_ptr = s8->alfa.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Dtildinv.fdim * s8->Dtildinv.sdim;
          dst_ptr = s8->Dtildinv.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Lt.fdim * s8->Lt.sdim;
          dst_ptr = s8->Lt.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Rtilde.fdim * s8->Rtilde.sdim;
          dst_ptr = s8->Rtilde.a.dv;
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }
        }
        if (mat[actele->mat-1].mattyp==m_viscohyper) {
          ARRAY4D *a;
          a = s8->his1;
          arr_length = a->fdim * a->sdim * a->tdim * a->fodim;
          dst_ptr = a->a.d4[0][0][0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }
        }

        break;
      }
#endif

#ifdef D_SHELL9
      case el_shell9: {
        INT arr_length;
        SHELL9* s9;
        DOUBLE* dst_ptr;

        s9 = actele->e.s9;
        if (s9->nhyb) {

          arr_length = s9->alfa.fdim * s9->alfa.sdim;
          dst_ptr = s9->alfa.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Dtildinv.fdim * s9->Dtildinv.sdim;
          dst_ptr = s9->Dtildinv.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->L.fdim * s9->L.sdim;
          dst_ptr = s9->L.a.da[0];
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Rtilde.fdim * s9->Rtilde.sdim;
          dst_ptr = s9->Rtilde.a.dv;
          for (el=0; el<arr_length; ++el) {
            *dst_ptr++ = *src_ptr++;
          }
        }

        if (s9->elewa != NULL) {
          INT kl;

          for (kl=0; kl<s9->num_klay; kl++) {
            INT num_mlay;
            INT ml;
            INT actlay = 0;

            num_mlay = s9->kinlay[kl].num_mlay;
            for (ml=0; ml<num_mlay; ml++) {
              S9_IP_WA *ipwa = s9->elewa[actlay].ipwa;

              /* check if there is a ipwa for this layer */
              if (ipwa != NULL) {
                INT k;
                INT ngauss;
                MULTIMAT *actmultimat;

                actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

                /* number of gausspoints in one layer*/
                ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

                switch (actmultimat->mattyp) {
                case m_pl_mises:
                case m_pl_dp:
                  for (k=0; k<ngauss; k++) {
                    ipwa[k].yip = *size_src_ptr++;
                    ipwa[k].epstn = *src_ptr++;
                    for (el=0; el<6; el++) {
                      ipwa[k].sig[el] = *src_ptr++;
                      ipwa[k].eps[el] = *src_ptr++;
                      ipwa[k].qn[el] = *src_ptr++;
                    }
                  }
                  break;
                case m_pl_epc:
                  for (k=0; k<ngauss; k++) {
                    ipwa[k].yip = *size_src_ptr++;
                    ipwa[k].kappa_t = *src_ptr++;
                    ipwa[k].kappa_c = *src_ptr++;
                    for (el=0; el<6; el++) {
                      ipwa[k].sig[el] = *src_ptr++;
                      ipwa[k].eps[el] = *src_ptr++;
                    }
                  }
                  break;
                case m_pl_hoff:
                  for (k=0; k<ngauss; k++) {
                    ipwa[k].yip = *size_src_ptr++;
                    ipwa[k].dhard = *src_ptr++;
                    for (el=0; el<6; el++) {
                      ipwa[k].sig[el] = *src_ptr++;
                      ipwa[k].eps[el] = *src_ptr++;
                      ipwa[k].dkappa[el] = *src_ptr++;
                      ipwa[k].gamma[el] = *src_ptr++;
                    }
                    for (el=0; el<9; el++) {
                      ipwa[k].rkappa[el] = *src_ptr++;
                    }
                  }
                  break;
                default:
                  printf("unsupported material %d for shell9 element",
                         actmultimat->mattyp);
                }
              }
              actlay += 1;
            }
          }
        }
        break;
      }
#endif

#ifdef D_WALL1
    case el_wall1: {
      W1_ELE_WA* elewa = actele->e.w1->elewa;
      if (elewa != NULL) {
        INT k;
        INT ngauss;
        INT minor;

        minor = MINOR_WALL1(actele);
        ngauss = element_info[el_wall1].variant[minor].gauss_number;

        /* This is w1_write_restart. Disguised as usual. */
        switch (mat[actele->mat-1].mattyp) {
        case m_pl_mises:
          dsassert(len >= ngauss*(1+4*3), "value entry too short");
          for (k=0; k<ngauss; ++k) {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].epstn = *src_ptr++;
            for (el=0; el<4; el++) {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].qn[el] = *src_ptr++;
            }
          }
          break;
        case m_pl_mises_3D:
          dsassert(len >= ngauss*(1+4*7), "value entry too short");
          for (k=0; k<ngauss; k++) {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].epstn = *src_ptr++;
            for (el=0; el<4; el++) {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].qn[el] = *src_ptr++;
              elewa->ipwa[k].sigc[el] = *src_ptr++;
              elewa->ipwa[k].sigi[el] = *src_ptr++;
              elewa->ipwa[k].epsi[el] = *src_ptr++;
              elewa->ipwa[k].di[el] = *src_ptr++;
            }
          }
          break;
        case m_pl_epc3D: /* ???
                          * (mat[actele->mat-1].mattyp == m_pl_epc3D ) */
          dsassert(len >= ngauss*(2+4*7), "value entry too short");
          for (k=0; k<ngauss; k++) {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].dlam[0] = *src_ptr++;
            elewa->ipwa[k].dlam[1] = *src_ptr++;
            for (el=0; el<4; el++) {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].sigc[el] = *src_ptr++;
              elewa->ipwa[k].grad[el] = *src_ptr++;
              elewa->ipwa[k].sigi[el] = *src_ptr++;
              elewa->ipwa[k].epsi[el] = *src_ptr++;
              elewa->ipwa[k].di[  el] = *src_ptr++;
            }
          }
          break;
        default:
          printf("WARNING: Unknown mattyp %d in wall1 restart\n",
                 mat[actele->mat-1].mattyp);
        }
      }
      break;
    }
#endif

#ifdef D_BEAM3
    case el_beam3:
      if (mat[actele->mat-1].mattyp == m_pl_mises ||
          mat[actele->mat-1].mattyp == m_pl_dp ||
          mat[actele->mat-1].mattyp == m_pl_epc) {
        INT j;
        INT size;
        DOUBLE* dst_ptr;

        /* Copy the values from the chunk to the element's array. */
        dst_ptr = actele->e.b3->elewa.a.da[0];
        size = actele->e.b3->elewa.fdim*actele->e.b3->elewa.sdim;
        dsassert(size <= len, "value entry too small");
        for (j=0; j<size; ++j) {
          *dst_ptr++ = *src_ptr++;
        }
      }
      break;
#endif

      /* There's nothing to be restored for interface and wallge
       * elements, right? */
    case el_interf:
      break;
    case el_wallge:
      break;

#ifdef D_FLUID2
    case el_fluid2: {
      FLUID_DYNAMIC *fdyn;
      DOUBLE* dst_ptr;
      INT k;

      /* only valid for single field problem (?) */
      /* I don't think so. */
      fdyn = alldyn[genprob.numff].fdyn;

      dsassert(len >= 3, "value item too small");
      dsassert(actele->e.f2->tau_old.Typ == cca_DV &&
               actele->e.f2->tau_old.fdim == 3, "uninitialized array");
      dst_ptr = &(actele->e.f2->tau_old.a.dv[0]);
      for (k=0; k<3; ++k) {
        *dst_ptr++ = *src_ptr++;
      }
      if (fdyn->surftens > 0) {
        dsassert(len >= 3+2*actele->numnp, "value item too small");
        dsassert(actele->e.f2->kappa_ND.Typ == cca_DA, "uninitialized array");
        dst_ptr = &(actele->e.f2->kappa_ND.a.dv[0]);
        if (actele->e.f2->fs_on>0) {
          for (k=0; k<2*actele->numnp; ++k) {
            *dst_ptr++ = *src_ptr++;
          }
        }
      }
      break;
    }
#endif

#ifdef D_FLUID3
    case el_fluid3: {
      DOUBLE* dst_ptr;
      INT k;

      dsassert(len >= 3, "value item too small");
      dsassert(actele->e.f3->tau_old.Typ == cca_DV &&
               actele->e.f3->tau_old.fdim == 3, "uninitialized array");
      dst_ptr = &(actele->e.f3->tau_old.a.dv[0]);
      for (k=0; k<3; ++k) {
        *dst_ptr++ = *src_ptr++;
      }
      break;
    }
#endif

    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

#undef get_id_part
#undef get_recv_numele

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack what we received.

  This is the hook function called by ``in_scatter_chunk`` in order
  to distribute certain values to nodes, elements or distributed
  vectors. This function or the ones called here must be changed in
  case an element is updated.

  There is a loop that calls this function for each source processor
  in turn. That is upon one call it has to distribute all items that
  come from one processor. The places these items go to are not
  consecutive. For this reason the \a context knows the local ids
  of these items. This has been figured out during initialization. See
  \a init_bin_in_field .

  \param *context         (i) the discretization the chunk belongs to
  \param *chunk           (i) the chunk that's read
  \param  type            (i) what kind of values are to be spread
  \param  array           (i) the row of the node array the value are
                              to go to; is ignored for some types
  \param *recv_buf        (i) buffer of double values
  \param  recv_count      (i) number of double values to be spread
  \param *recv_size_buf   (i) buffer of integer values
  \param *recv_size_count (i) buffer of integer values
  \param  src             (i) the rank of the source processor

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void in_unpack_items(struct _BIN_IN_FIELD *context,
                     struct _BIN_IN_CHUNK *chunk,
                     CHUNK_CONTENT_TYPE type,
                     NODE_ARRAY array,
                     DOUBLE *recv_buf,
                     INT recv_count,
                     INT *recv_size_buf,
                     INT recv_size_count,
                     INT src)
{
  INT value_length;
  INT size_length;

#ifdef DEBUG
  dstrc_enter("in_unpack_items");
#endif

  value_length = chunk->value_entry_length;
  size_length  = chunk->size_entry_length;

  switch (type) {
  case cc_node_array:
    in_unpack_node_arrays(context, chunk, array, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  case cc_dist_vector:
    in_unpack_dist_vector(context, chunk, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  case cc_restart_element:
    in_unpack_restart_element(context, chunk, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  default:
    dserror("unsupported chunk type %d", type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of different element types in this field.

  This simply counts the different types. That's particularly easy
  (and quick) because each field context knows what kinds of elements
  there are.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT count_element_types(struct _BIN_OUT_FIELD* field)
{
  INT i;
  INT j;
  INT counter = 0;

#ifdef DEBUG
  dstrc_enter("count_element_types");
#endif

  for (i=1; i<el_count; ++i) {
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (field->element_flag[i][j]) {
        counter++;
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return counter;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of different element variants to the given
  type (major number) in this field.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT count_element_variants(struct _BIN_OUT_FIELD* field, ELEMENT_TYP type)
{
  INT j;
  INT counter = 0;

#ifdef DEBUG
  dstrc_enter("count_element_variants");
#endif

  for (j=0; j<MAX_EL_MINOR; ++j) {
    if (field->element_flag[type][j]) {
      counter++;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return counter;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the mesh connectivity.

  The mesh connectivity consists of all the ids of those nodes that
  are connected to one particular element. So the main task here is to
  find the maximum number of nodes per element. This a done without
  looping all elements, instead the used element types are looped. To
  each element type (major/minor) the number of nodes are known,
  thanks to the global \a element_info .

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_mesh_item_length(struct _BIN_OUT_FIELD* context,
                           INT* value_length,
                           INT* size_length)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("find_mesh_item_length");
#endif

  *value_length = 0;
  *size_length = 0;

  /* We simply need to find the maximum element node number in the
   * discretization. */
  for (i=1; i<el_count; ++i) {
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (context->element_flag[i][j]) {
        *size_length = MAX(*size_length, element_info[i].variant[j].node_number);
      }
    }
  }

  /* And there are three more integer values: The element id and both
   * major and minor element type. */
  *size_length += 2 + 1;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element stresses.

  We want to store the real stress array for each
  element. Postprocessing can be done later. However the real array is
  not always available. Each element is different. Some elements do
  even extrapolate their stresses to the nodes. This should better be
  done by the filter, but for now we'll have to support this. Thus we
  need to find the size of the node based stress array, too.

  Of course extrapolating to the nodes requires that all elements at
  the node agree on the nodal stress. In effect this demands that only
  one type of element is used in the discretization.

  For many element types it's sufficient to lookup the \a element_info
  table to find the stress array's size. But dynamic elements like
  shell9 have to be treated specially.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_stress_item_length(struct _BIN_OUT_FIELD* context,
                             INT* el_value_length,
                             INT* el_size_length,
                             INT* nd_value_length,
                             INT* nd_size_length)
{
  INT i;
  INT j;

#ifdef DEBUG
  dstrc_enter("find_stress_item_length");
#endif

  *el_value_length = 0;
  *el_size_length = 0;
  *nd_value_length = 0;
  *nd_size_length = 0;

  /* We need to find the maximum size of the stress array. */
  for (i=1; i<el_count; ++i) {
    for (j=0; j<MAX_EL_MINOR; ++j) {
      if (context->element_flag[i][j]) {

        /* element based stresses */
        if (element_info[i].variant[j].el_stress_matrix_size != -1) {
          *el_value_length = MAX(*el_value_length,
                                 element_info[i].variant[j].el_stress_matrix_size);
        }
        else {
          switch (i) {
#ifdef D_SHELL9
          case el_shell9: {     /* multi layer shell element */
            ELEMENT* actele;
            SHELL9* s9;

            /*
             * If there are shell9 elements there are only shell9
             * elements. Additionally we expect all elements to have
             * the same dimension parameters. We take the values we
             * find in the first element.
             *
             * (We know that there's at least one element on each
             * processor.) */
            actele = context->actpart->pdis[context->disnum].element[0];
            dsassert(actele->eltyp == el_shell9, "shell9 expected");

            s9 = actele->e.s9;

            /* There are 6 stress components. */
            *el_value_length = MAX(*el_value_length,
                                   6*context->s9_layers*s9->nGP[0]*s9->nGP[1]*s9->nGP[2]);
            break;
          }
#endif
          default:
            dserror("no specific treatment of element type %d", i);
          }
        }

        /* stresses extrapolated to the nodes */
        if (element_info[i].variant[j].nd_stress_matrix_size != -1) {
          *nd_value_length = MAX(*nd_value_length,
                                 element_info[i].variant[j].nd_stress_matrix_size);
        }
        else {
          switch (i) {
          default:
            dserror("no specific treatment of element type %d", i);
          }
        }
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element's restart information.

  There could be an entry in the element info table that covers the
  item lengths required for restart. But most elements have flags that
  determine the required sizes. So we'll need special code anyway.

  This is a collective call, that is there's communication involved
  here. This way we can loop all the elements and find the sizes we
  need.

  Normally we'll try to avoid such element loops. That's what the
  element_flag variable is about. But here we have to. We can neither
  encode all possible element variants in one gigantic array nor can
  we stick to the most space consuming version for lack of knowledge.

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \param context      (i) pointer to an already set up output context
  \param value_length (o) the number of double to store per element
  \param size_length  (o) the number of integer to store per element

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_restart_item_length(struct _BIN_OUT_FIELD* context,
                              INT* value_length,
                              INT* size_length)
{
  INT local_lengths[2];
#ifdef PARALLEL
  INT global_lengths[2];
#else
  INT* global_lengths = local_lengths;
#endif

#ifdef DEBUG
  dstrc_enter("find_restart_item_length");
#endif

  /* Start at zero. */
  *value_length = 0;
  *size_length  = 0;

#ifdef D_SHELL8
  if (context->is_shell8_problem) {
    INT i;
    PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
    for (i=0; i<actpdis->numele; ++i) {
      ELEMENT* actele;

      actele = actpdis->element[i];
      if (actele->eltyp == el_shell8) {
        INT len = 0;
        SHELL8* s8;

        s8 = actele->e.s8;

        /* collect the array sizes */
        if (s8->nhyb) {
          len += s8->alfa.fdim * s8->alfa.sdim;
          len += s8->Dtildinv.fdim * s8->Dtildinv.sdim;
          len += s8->Lt.fdim * s8->Lt.sdim;
          len += s8->Rtilde.fdim * s8->Rtilde.sdim;
        }
        if (mat[actele->mat-1].mattyp==m_viscohyper) {
          ARRAY4D *a;
          a = s8->his1;
          len += a->fdim * a->sdim * a->tdim * a->fodim;
        }

        *value_length = MAX(*value_length, len);
      }
    }

    /* It's a shell8 problem. There are no other elements. */
    goto end;
  }
#endif

#ifdef D_SHELL9
  /*count_element_variants(context, el_shell9);*/
  if (context->is_shell9_problem) {
    INT i;
    PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
    for (i=0; i<actpdis->numele; ++i) {
      ELEMENT* actele;

      actele = actpdis->element[i];
      if (actele->eltyp == el_shell9) {
        INT len = 0;
        INT slen = 0;
        SHELL9* s9;

        s9 = actele->e.s9;

        /* collect the array sizes */
        if (s9->nhyb) {
          len += s9->alfa.fdim * s9->alfa.sdim;
          len += s9->Dtildinv.fdim * s9->Dtildinv.sdim;
          len += s9->L.fdim * s9->L.sdim;
          len += s9->Rtilde.fdim * s9->Rtilde.sdim;
        }
        if (s9->elewa != NULL) {
          INT kl;
          INT actlay = 0;

          for (kl=0; kl<s9->num_klay; kl++) {
            INT num_mlay;
            INT ml;

            num_mlay = s9->kinlay[kl].num_mlay;
            for (ml=0; ml<num_mlay; ml++) {

              /* check if there is a ipwa for this layer */
              if (s9->elewa[actlay].ipwa != NULL) {
                INT ngauss;
                MULTIMAT *actmultimat;

                actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

                /* number of gausspoints in one layer*/
                ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

                switch (actmultimat->mattyp) {
                case m_pl_mises:
                case m_pl_dp:
                  len += ngauss*(1+6*3);
                  slen += ngauss;
                  break;
                case m_pl_epc:
                  len += ngauss*(2+6*2);
                  slen += ngauss;
                  break;
                case m_pl_hoff:
                  len += ngauss*(1+6*4+9);
                  slen += ngauss;
                  break;
                default:
                  printf("unsupported material %d for shell9 element",
                         actmultimat->mattyp);
                }
              }
              actlay += 1;
            }
          }
        }
        *value_length = MAX(*value_length, len);
        *size_length  = MAX(*size_length, slen);
      }
    }

    /* It's a shell9 problem. There are no other elements. */
    goto end;
  }
#endif

#ifdef D_WALL1
  /* if there are wall1 elements */
  /*
   * There are huge differences depending on the material. This forces
   * us to iterate the elements. */
  if (count_element_variants(context, el_wall1) > 0) {
    INT i;
    PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
    for (i=0; i<actpdis->numele; ++i) {
      ELEMENT* actele;

      actele = actpdis->element[i];
      if (actele->eltyp == el_wall1) {
        WALL1* w1;

        w1 = actele->e.w1;
        if (w1->elewa != NULL) {
          INT ngauss;
          ngauss = w1->nGP[0]*w1->nGP[1];
          switch (mat[actele->mat-1].mattyp) {
          case m_pl_mises:
            *value_length = MAX(*value_length, ngauss*(1+4*3));
            *size_length  = MAX(*size_length,  ngauss*(1));
            break;
          case m_pl_mises_3D:
            *value_length = MAX(*value_length, ngauss*(1+4*7));
            *size_length  = MAX(*size_length,  ngauss*(1));
            break;
          case m_pl_epc3D:
            *value_length = MAX(*value_length, ngauss*(2+4*7));
            *size_length  = MAX(*size_length,  ngauss*(1));
            break;
          default:
            /*
             * This is strange. Do we want these messages when we
             * restart wall meshed with linear material? */
            printf("unsupported material %d for wall1 element",
                   mat[actele->mat-1].mattyp);
          }
        }
      }
    }
  }
#endif

#ifdef D_BEAM3
  /* if there are beam3 elements */
  /* Again the material matters. */
  if (count_element_variants(context, el_beam3) > 0) {
    INT i;

    PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
    for (i=0; i<actpdis->numele; ++i) {
      ELEMENT* actele;

      actele = actpdis->element[i];
      if (actele->eltyp == el_beam3) {
        if(mat[actele->mat-1].mattyp == m_pl_mises ||
           mat[actele->mat-1].mattyp == m_pl_dp ||
           mat[actele->mat-1].mattyp == m_pl_epc) {
          *value_length = MAX(*value_length, actele->e.b3->nGP[0]*5*5*40);
        }
      }
    }
  }
#endif

  /* There's nothing to be saved for interface and wallge elements,
   * right? */

#ifdef D_FLUID2
  if (count_element_variants(context, el_fluid2) > 0) {
    INT j;
    INT len = 3;
    FLUID_DYNAMIC *fdyn;

    /* only valid for single field problem (?) */
    /* I don't think so. */
    fdyn = alldyn[genprob.numff].fdyn;

    /* We have the stabilization parameter history (tau) and maybe kappa. */
    if (fdyn->surftens > 0) {
      for (j=0; j<MAX_EL_MINOR; ++j) {
        if (context->element_flag[el_fluid2][j]) {
          INT numnp;
          numnp = element_info[el_fluid2].variant[j].node_number;
          len = MAX(len, 3+2*numnp);
        }
      }
    }
    *value_length = MAX(*value_length, len);
  }
#endif

#ifdef D_FLUID3
  if (count_element_variants(context, el_fluid3) > 0) {
    /* tau. The stabilization parameter history. */
    *value_length = MAX(*value_length, 3);
  }
#endif


  /* Finally communicate what we've found. We need to do this because
   * there's no guarantee that the most demanding element type is used
   * on all processors. */
  local_lengths[0] = *value_length;
  local_lengths[1] = *size_length;
#ifdef PARALLEL
  MPI_Allreduce(local_lengths, global_lengths, 2, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#endif
  *value_length = global_lengths[0];
  *size_length  = global_lengths[1];

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Get the position of this field in the global field array.

  We need to store the field's position along with the field's type
  and the discretization number in the control file in order to
  identify the discretization. (Here the distinction between field and
  discretization matters.) The reason are problemtypes like SSI that
  contain more than one field of one type, the field type is no unique
  criterion.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT get_field_position(struct _BIN_OUT_FIELD  *context)
{
  INT field_pos;

#ifdef DEBUG
  dstrc_enter("get_field_position");
#endif

  for (field_pos=0; field_pos<genprob.numfld; ++field_pos) {
    if (context->actfield == &(field[field_pos])) {
      break;
    }
  }

  if (field_pos==genprob.numfld) {
    dserror("unregistered field object");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return field_pos;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the start of a group that describes one aspect of a
  discretization.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_main_group_head(struct _BIN_OUT_FIELD  *context, CHAR* name)
{
  INT field_pos;
#ifdef DEBUG
  dstrc_enter("out_main_group_head");
#endif

  field_pos = get_field_position(context);

  fprintf(bin_out_main.control_file, "%s:\n"
          "    field = \"%s\"\n"
          "    field_pos = %d\n"
          "    discretization = %d\n"
          "",
          name,
          fieldnames[context->actfield->fieldtyp],
          field_pos,
          context->disnum);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write all results for one step.

  Write results for potprocessing. All algorithms call this function
  (if they support binary output).

  This function can be called many times in a row per time step. But
  be careful not to mix calls of this function with calls to output
  restart data.

  \param context  pointer to an already set up output context
  \param time     current time
  \param step     current step count
  \param place    node array row that contains the results
  \param flags    the type of output needed; flags might be or'ed together

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_results(struct _BIN_OUT_FIELD* context,
                 DOUBLE time,
                 INT step,
                 INT place,
                 OUT_FLAGS flags)
{
  INT rank = context->actintra->intra_rank;

#ifdef DEBUG
  dstrc_enter("out_results");
#endif

  dsassert(flags != 0, "no output flags; I've got nothing to write");

  if (rank == 0) {
    static DOUBLE last_time = -1;
    static INT last_step = -1;
    static INT last_discr = -1;
    static INT last_field_pos = -1;

    INT field_pos;

    field_pos = get_field_position(context);

    /* Only start a new group if we have a new time step. This way
     * out_results can be called more often that once per step. Be
     * careful, however, not to mix it with restart output
     * calls. Restart starts a new group as well and everything will
     * be added to the latest group no matter who created it. */

    if ((last_step != step) || (last_time != time) ||
        (last_field_pos != field_pos) || (last_discr != context->disnum)) {
      out_main_group_head(context, "result");
      fprintf(bin_out_main.control_file,
              "    time = %f\n"
              "    step = %d\n"
              "\n",
              time, step);
      last_step = step;
      last_time = time;
      last_field_pos = field_pos;
      last_discr = context->disnum;
    }
  }

  if (flags & OUTPUT_DISPLACEMENT) {

#ifdef D_SHELL8
    /* shell8 elements have three additional dofs, all six living in
     * one node array's row. Thus we only need one call to output
     * them. */
    if (context->is_shell8_problem) {
      out_node_chunk(context, "displacement", cc_displacement, 6, 0, place);
    }
    else
#endif

      /* Displacement of normal nodes. */
      out_node_chunk(context, "displacement", cc_displacement, genprob.ndim, 0, place);

#ifdef D_SHELL9
    /* shell9 elements have layers with their own displacements. */
    if (context->is_shell9_problem) {
      INT numnp;
      INT layer_count;

      numnp = element_info[el_shell9].variant[context->s9_minor].node_number;

      if ((context->s9_minor == MINOR_SHELL9_4_22) ||
          (context->s9_minor == MINOR_SHELL9_4_33)) {
        layer_count = 1;
      }
      else {
        layer_count = 2;
      }

      out_element_chunk(context, "shell9_displacement", cc_shell9_displacement,
                        3*numnp*(layer_count*context->s9_layers+1), 0,
                        place);
    }
#endif
  }

  if (flags & OUTPUT_VELOCITY) {
    out_node_chunk(context, "velocity", cc_velocity, genprob.ndim, 0, place);
  }

  if (flags & OUTPUT_PRESSURE) {
    out_node_chunk(context, "pressure", cc_pressure, 1, 0, place);
  }

  if (flags & OUTPUT_STRESS) {
    INT el_value_length;
    INT nd_value_length;
    INT el_size_length;
    INT nd_size_length;

    find_stress_item_length(context,
                            &el_value_length, &el_size_length,
                            &nd_value_length, &nd_size_length);
    if ((el_value_length > 0) || (el_size_length > 0)) {
      /*
       * use the general name "stress" for the most natural gauss
       * point stresses. */
      out_element_chunk(context, "stress", cc_el_stress, el_value_length, el_size_length, 0);
    }
    if ((nd_value_length > 0) || (nd_size_length > 0)) {
      out_element_chunk(context, "nd_stress", cc_nd_stress, nd_value_length, nd_size_length, 0);
    }
  }

  if (flags & OUTPUT_CONTACT) {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "contact" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_EIGENMODES) {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "eigenmodes" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_THICKNESS) {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "thickness" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_AXI_LOADS) {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "axishell" RED_LIGHT " loads not yet supported\n" END_COLOR);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
