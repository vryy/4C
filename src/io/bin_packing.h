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

Two very important types are defined in this head file:

OUT_FLAGS -- The coarse grained types of postprocessor output ccarat
supports. These types are meant to be or'ed together and handed to the
public out_results function.

CHUNK_CONTENT_TYPE -- The fine grained internal output type. This type
specifies how large the chunk's entries are (or how that's figured
out), how the values are collected and how they are put pack.

\author u.kue
\date 08/04

*/

#ifndef BIN_PACKING_H
#define BIN_PACKING_H

#include "../headers/standardtypes.h"
#include "../global_full/element_info.h"



/*
 * Structures defined in another place. We need to introduce them here
 * because we have functions that expect pointers to them as
 * arguments. */
struct _BIN_OUT_FIELD;
struct _BIN_OUT_CHUNK;

struct _BIN_IN_FIELD;
struct _BIN_IN_CHUNK;


/*----------------------------------------------------------------------*/
/*!
  \brief The results that are to be written are indicated by a flag.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef INT OUT_FLAGS;

/* possible flags */
#define OUTPUT_DISPLACEMENT 0x0001
#define OUTPUT_VELOCITY     0x0002
#define OUTPUT_PRESSURE     0x0004
#define OUTPUT_STRESS       0x0008
#define OUTPUT_CONTACT      0x0010
#define OUTPUT_EIGENMODES   0x0020
#define OUTPUT_THICKNESS    0x0040
#define OUTPUT_AXI_LOADS    0x0080
#define OUTPUT_ACCELERATION 0x0100 /* not implemented */

/*----------------------------------------------------------------------*/
/*!
  \brief Internal type to distinguish the information stored in a
  chunk.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef enum _CHUNK_CONTENT_TYPE {
  cc_none,

  /* specific results */
  cc_displacement,
  cc_velocity,
  cc_pressure,
  cc_el_stress,
  cc_nd_stress,
  cc_domain,

  cc_dist_vector,

  /* element specific output */
  cc_restart_element,
  cc_shell8_director,
  cc_shell9_coords,
  cc_shell9_displacement,

  /* general output */
  cc_mesh,
  cc_coords,
  cc_node_array
} CHUNK_CONTENT_TYPE;


/*----------------------------------------------------------------------*/
/*!
  \brief Mark all element types that are used in this discretization.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_find_element_types(struct _BIN_OUT_FIELD *context,
                            INTRA *actintra,
                            PARTDISCRET *actpdis);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of different element types in this field.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT count_element_types(struct _BIN_OUT_FIELD* field);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of different element variants to the given
  type (major number) in this field.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT count_element_variants(struct _BIN_OUT_FIELD* field, ELEMENT_TYP type);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the mesh connectivity.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_mesh_item_length(struct _BIN_OUT_FIELD* context,
                           INT* value_length,
                           INT* size_length);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element stresses.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_stress_item_length(struct _BIN_OUT_FIELD* context,
                             INT* el_value_length,
                             INT* el_size_length,
                             INT* nd_value_length,
                             INT* nd_size_length);


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

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_restart_item_length(struct _BIN_OUT_FIELD* context,
                              INT* value_length,
                              INT* size_length);


/*----------------------------------------------------------------------*/
/*!
  \brief Get the position of this field in the global field array.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT get_field_position(struct _BIN_OUT_FIELD  *context);


/*----------------------------------------------------------------------*/
/*!
  \brief Write the start of a group that describes one aspect of a
  discretization.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_main_group_head(struct _BIN_OUT_FIELD  *context, CHAR* name);


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the values to be saved in some buffer to send them to
  their writing processor.

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
                    INT dst_num);


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack what we received.

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
                     INT src);


/*----------------------------------------------------------------------*/
/*!
  \brief Write all results for one step.

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
                 OUT_FLAGS flags);


#endif
