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

#ifdef BINIO

/*!
\addtogroup IO
*//*! @{ (documentation module open)*/

#ifndef BIN_PACKING_H
#define BIN_PACKING_H

#include "../headers/standardtypes.h"
#include "../global_full/global_element_info.h"



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
  cc_stress,
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


#ifdef D_SHELL8

/*----------------------------------------------------------------------*/
/*!
  \brief Shell8 specific setup of the discretization output object.

  Shell8 has six displacement dofs per node. The three normal ones and
  the three steaming from the director vector. These additional three
  must be handled here.

  Currently a discretization must contains just one type of shell8
  elements if it contains shell8 elements at all.

  \param context      (i) pointer to a context variable during setup.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_shell8_setup(struct _BIN_OUT_FIELD *context);

#endif

#ifdef D_SHELL9

/*----------------------------------------------------------------------*/
/*!
  \brief Shell9 specific setup of the discretization output object.

  Shell9 is a layered element with many displacement dofs per
  node. These demand a special threatment.

  Currently a discretization must contains just one type of shell9
  elements if it contains shell9 elements at all.

  \param context      (i) pointer to a context variable during setup.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_shell9_setup(struct _BIN_OUT_FIELD *context);

#endif

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
                           INT* size_length);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element stresses.

  We want to store the real stress array for each
  element. Postprocessing can be done later.

  For many element types it's sufficient to lookup the \a element_info
  table to find the stress array's size. But dynamic elements like
  shell9 have to be treated specially.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_stress_item_length(struct _BIN_OUT_FIELD* context,
                             INT* value_length,
                             INT* size_length);


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

  \param context      (i) pointer to an already set up output context
  \param value_length (o) the number of double to store per element
  \param size_length  (o) the number of integer to store per element

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
                    INT dst_num);


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
                     INT src);


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
                 OUT_FLAGS flags);


#endif

/*! @} (documentation module close)*/
#endif
