/*!
\file
\brief Functions for binary IO to a single file.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

The new aggressiv parallelization approach demands a parallel io
facility because the data is scattered amoung processors. However, one
wants to be independent of plattform and number of processors. That is
the output has to go to one big file (maybe a few files) in a way that
can be read on any hardware. The entries in this files are to be
sorted. But at the same time we cannot afford (at least that's the
idea, we'll have to measure whether it's true) to store each single
entry (node, element) individually. The answer is to assign each
processor a consecutive set of entries. That demands a lot of
communication to get the values to be stored to the processors that
are to do it. The same happens while reading such a file. Nevertheless
this approach allows to write an independent file with a minimum of
write calls.

To achieve highes performance each processor would have to store its
data to its own file. That's the approach when there are local file
systems and we really need highest performance. Postprocessing and
restart would be much more difficult. But that's another
approach. It's not treated here.

Output and input both follow the same pattern. Both are organized in a
hierarchical way. At the top there is one global object that contains
the overall information (BIN_OUT_MAIN vs. BIN_IN_MAIN). Below that
there are the discretization specific objects (BIN_OUT_FIELD and
BIN_IN_FIELD). These are constructed for each discretization that's
going to be read or written. The setup of these structures is quite
demanding so we do it just once per discretization. At the bottom
finally there are the chunks of data (BIN_OUT_CHUNK and
BIN_IN_CHUNK). These represent one piece of information, collected
from all the nodes or element inside the discretization. We build and
destroy them whenever needed because they consume a lot of memory and
are very specific to the io type.

All io is done using chunks.

\author u.kue
\date 08/04

*/

#ifdef BINIO

/*!
\addtogroup IO
*//*! @{ (documentation module open)*/

#ifndef BIN_SINGLEFILE_H
#define BIN_SINGLEFILE_H

#include "io.h"
#include "io_packing.h"


/*----------------------------------------------------------------------*/
/*!
  \brief The structure used for writing one chuck of data inside a field.

  Such a chunk consists of many items. Each item might have some
  double values and some integer values. The number of these values
  are the same for all items.

  There are two possible origins for those items: nodes and
  elements.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_OUT_CHUNK {

  BIN_OUT_FIELD* field;         /* the field we belong to */

  CHUNK_TYPE type;              /* the kind of entities we store */

  INT num;                      /* number of nodes or elements */
  INT fullarrays;               /* Number of processors that have the
                                 * maximal number of items. All others
                                 * have one less. */
  INT first_id;                 /* the first id on this proc (Id_loc) */

  INT value_entry_length;       /* number of doubles per item */
  INT size_entry_length;        /* number of ints per item */

  INT value_count;              /* total number of doubles on this proc */
  INT size_count;               /* total number of ints on this proc */

  DOUBLE* out_values;           /* all double values we write */
  INT* out_sizes;               /* all int values we write */

  MAP group;                    /* the group that describes this chunk */

  DIST_VECTOR* vectors;         /* in case we store distributed vectors */

} BIN_OUT_CHUNK;


/*----------------------------------------------------------------------*/
/*!
  \brief The structure used for reading one chuck of data inside a field.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_IN_CHUNK {

  CHUNK_TYPE type;              /* whether we handle nodes or elements */

  BIN_IN_FIELD* field;          /* the discretization we belong to */

  INT num;
  INT fullarrays;
  INT first_id;

  INT value_entry_length;       /* number of doubles per item */
  INT size_entry_length;        /* number of ints per item */

  INT value_count;              /* total number of doubles */
  INT size_count;               /* total number of ints */

  DOUBLE* in_values;
  INT* in_sizes;

  MAP* group_info;

  DIST_VECTOR* vectors;         /* In case we load distributed vectors. */
  SPARSE_ARRAY sysarray;        /* To be able to interpret a dist. vector */
  SPARSE_TYP sysarray_typ;      /* we need to know its solver. */

} BIN_IN_CHUNK;


/*----------------------------------------------------------------------*/
/*!
  \brief Open the data files for output.

  \param context      (i) pointer to a context variable during setup.
  \param files        (o) the pair of files to be opened.
  \param variant      (i) tag to distinguish different pairs. Might be
                          NULL if files is already initialized.
  \param step         (i) first time step to be written to this file

  \warning \a variant must be a static string. A pointer to it is
  kept.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_open_data_files(BIN_OUT_FIELD *context,
                         BIN_DATA_FILES *files,
                         CHAR *variant,
                         INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Init the data structure to write one set of results.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_out_chunk(BIN_OUT_FIELD* context,
                        BIN_OUT_CHUNK* chunk,
                        CHUNK_TYPE type,
                        INT value_entry_length,
                        INT size_entry_length);


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The output chunk's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_out_chunk(BIN_OUT_CHUNK* data);


/*----------------------------------------------------------------------*/
/*!
  \brief Send node values to the processors where they are written.

  Here is the algorithm that gathers the values on the processors
  where they life, packs them and sends them to their writing
  processors. We have to be careful to stay independent of the actual
  values transfered. The flag indicated which values are to be
  gathered. The actual process of gathering depends on which values
  are to be handled and is done in appropriate functions.

  The argument \a array has different meanings depending of the
  flag. If it's equal to OUTPUT_NODE_ARRAY, that is we are going to
  write node arrays, the value \a array indicates which one is going
  to be written. However when we are about to write a result chunk \a
  array indicates the place (line) in the solution array that's
  written.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_gather_values(BIN_OUT_FIELD* context,
                       BIN_OUT_CHUNK* chunk,
                       CHUNK_CONTENT_TYPE type,
                       INT array);


/*----------------------------------------------------------------------*/
/*!
  \brief Switch to the result files.

  This sets the output context to write anything to the result
  files. The setting remains active until it's changed explicitly.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void out_activate_result(BIN_OUT_FIELD *context);


/*----------------------------------------------------------------------*/
/*!
  \brief Switch to the restart files.

  This sets the output context to write anything to the restart
  files. The setting remains active until it's changed explicitly.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void out_activate_restart(BIN_OUT_FIELD *context);


/*----------------------------------------------------------------------*/
/*!
  \brief Write the gathered data.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_write_chunk(BIN_OUT_FIELD *context,
                     BIN_OUT_CHUNK* array_data,
                     CHAR* entry_name,
                     INT restart);


/*----------------------------------------------------------------------*/
/*!
  \brief Write one node array.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_node_arrays(BIN_OUT_FIELD* context,
                     NODE_ARRAY array);


/*----------------------------------------------------------------------*/
/*!
  \brief Write one node array.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_node_chunk(BIN_OUT_FIELD* context,
                    CHAR* chunk_name,
                    CHUNK_CONTENT_TYPE type,
                    INT value_length,
                    INT size_length,
                    INT array);


/*----------------------------------------------------------------------*/
/*!
  \brief Write one element array.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_element_chunk(BIN_OUT_FIELD* context,
                       CHAR* chunk_name,
                       CHUNK_CONTENT_TYPE type,
                       INT value_length,
                       INT size_length,
                       INT array);


/*----------------------------------------------------------------------*/
/*!
  \brief Write a distributed vector.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_distvec_chunk(BIN_OUT_FIELD* context,
                       CHAR* chunk_name,
                       INT length,
                       DIST_VECTOR* vectors);


/*----------------------------------------------------------------------*/
/*!
  \brief Read.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void in_read_chunk(BIN_IN_FIELD *context,
                   BIN_IN_CHUNK *array_data);


/*----------------------------------------------------------------------*/
/*!
  \brief Read on node array of the given discretization and distribute
  it to the nodes.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void in_node_arrays(BIN_IN_FIELD* context,
                    MAP* result_info,
                    NODE_ARRAY array);


/*----------------------------------------------------------------------*/
/*!
  \brief Read node values.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void in_node_chunk(BIN_IN_FIELD* context,
                   MAP* result_info,
                   CHAR* group_name,
                   NODE_ARRAY array);


/*----------------------------------------------------------------------*/
/*!
  \brief Read element data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void in_element_chunk(BIN_IN_FIELD* context,
                      MAP* result_info,
                      CHAR* group_name,
                      CHUNK_CONTENT_TYPE type);


/*----------------------------------------------------------------------*/
/*!
  \brief Read distributed vectors.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void in_distvec_chunk(BIN_IN_FIELD* context,
                      MAP* result_info,
                      CHAR* group_name,
                      INT length,
                      DIST_VECTOR* vectors);


/*----------------------------------------------------------------------*/
/*!
  \brief Find the description of the restart data in the control file

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
MAP *in_find_restart_group(BIN_IN_FIELD* context, INT disnum, INT step);

#endif

/*! @} (documentation module close)*/
#endif
