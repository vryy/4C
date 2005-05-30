/*!
\file
\brief Code that is common to all filters.

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

All filters need to read the ccarat output. Thus there are some common
functions.

\author u.kue
\date 10/04

*/


#ifndef POST_COMMON_H
#define POST_COMMON_H


/*
 * Preprocessor commands that make it possible to use the same sources
 * and configuration like the ccarat binary. */

#define FILTER

#ifndef BINIO
#define BINIO
#warning BINIO assumed for filters
#endif

#ifdef PARALLEL
#undef PARALLEL
#endif

#ifdef SPOOLES_PACKAGE
#undef SPOOLES_PACKAGE
#endif

#ifdef PERF
#undef PERF
#endif

#ifdef NO_TEXT_OUTPUT
#undef NO_TEXT_OUTPUT
#endif

#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"
#include "../io/io_packing.h"
#include "../io/io_elements.h"


extern struct _FILES           allfiles;

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
  \brief The data necessary to access one result.

  The point is that each result might define it's own input files. If
  it doesn't do so it inherits it's files from the previous result of
  the discretization.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
typedef struct _RESULT_DATA {
  FILE* value_file;
  FILE* size_file;

  /* position of this result in the global result array */
  INT pos;

  /* the field this result belongs to */
  struct _FIELD_DATA* field;

  MAP* group;

} RESULT_DATA;


/*----------------------------------------------------------------------*/
/*!
  \brief The data necessary to access one chunk.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
typedef struct _CHUNK_DATA {

  /* the result this chunk belongs to */
  RESULT_DATA* result;

  MAP* group;

  INT value_entry_length;
  INT value_offset;
  INT size_entry_length;
  INT size_offset;

  /* Buffers that contain the data of the current entry. */
  /* Never write to these buffers! */
  INT* size_buf;
  DOUBLE* value_buf;

#ifndef LOWMEM
  /* Internal memory that holds the whole chunk. Ignore it. Never
   * access it directly. It's a performance hack. It's for internal
   * use only. If the chunks don't fit into the memory we do it
   * without these buffers. You have been warned. */
  INT* size_data;
  DOUBLE* value_data;
#endif
} CHUNK_DATA;


/*----------------------------------------------------------------------*/
/*!
  \brief The translation from external numbers to internal enums.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
typedef struct _TRANSLATION_TABLE {
  INT* table;
  MAP* group;
  INT length;
} TRANSLATION_TABLE;


/*----------------------------------------------------------------------*/
/*!
  \brief The meta data that belongs to one discretization.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef struct _FIELD_DATA {
  struct _PROBLEM_DATA* problem;

  FIELDTYP type;

  INT field_pos;
  INT disnum;
  INT numele;
  INT numnp;
  INT numdf;

  MAP* group;
  CHAR* name;

  /* the mesh files */
  /* a fake to be able to use the chunk data structure */
  RESULT_DATA head;

  /* Support structures to read the field chunks. */
  /* This makes it easy to access the mesh data. */
  CHUNK_DATA ele_param;
  CHUNK_DATA mesh;
  CHUNK_DATA coords;

#ifdef D_SHELL8
  INT is_shell8_problem;
#endif

#ifdef D_SHELL9
  INT is_shell9_problem;
  INT s9_smooth_results;
  INT s9_layers;
#endif
} FIELD_DATA;


/*----------------------------------------------------------------------*/
/*!
  \brief The problems most general data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef struct _PROBLEM_DATA {
  PROBLEM_TYP type;
  INT ndim;
  INT num_discr;
  FIELD_DATA* discr;

  /* total number of items over all fields */
  /* Some elements introduce artificial nodes (shell8,shell9) in the
   * postprocessing step. We assume that nodes are consecutively
   * numbered from 0 to numnp-1 and assign numbers starting with numnp
   * to these new ones. */
  INT numele;
  INT numnp;

  CHAR basename[100];
  MAP control_table;

  /* start, stop and step numbers. a python like slice. */
  /* We don't have to read each result. This is set by command line
   * arguments. */
  INT start;
  INT end;
  INT step;

  INT num_results;
  MAP** result_group;

  /* translation tables from external to internal numbers */
  TRANSLATION_TABLE element_type;
  TRANSLATION_TABLE distype;

  CHAR input_dir[100];

} PROBLEM_DATA;



/*----------------------------------------------------------------------*/
/*!
  \brief This is ccarat setup in a hurry.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void setup_filter(CHAR* output_name, MAP* control_table, CHAR* basename);


/*----------------------------------------------------------------------*/
/*!
  \brief Filter log output

  Writes the message to the log file. If the given level is small
  enough it prints the message on screen as well.

  \param level (i) level of this message
  \param msg   (i) any text followed by printf like arguments

  \author u.kue
  \date 01/05
*/
/*----------------------------------------------------------------------*/
void post_log(INT level, CHAR* msg, ...);


/*----------------------------------------------------------------------*/
/*!
  \brief Init a translation table.

  The purpose of these tables is to find the internal enum value to an
  external number read from the data file. We need to translate each
  enum value we read. Otherwise the number written would be an
  implementation detail, thus we could not change our implementation
  once we had some binary files.

  \param table (o) table to be initialized
  \param group (i) control file group that defines the external values
  \param names (i) list of names ordered by internal number, NULL terminated

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void init_translation_table(TRANSLATION_TABLE* table,
                            MAP* group,
                            CHAR** names);


/*----------------------------------------------------------------------*/
/*!
  \brief Clear the memory occupied by the translation table.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_translation_table(TRANSLATION_TABLE* table);


/*----------------------------------------------------------------------*/
/*!
  \brief Extract one discretization's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_field_data(PROBLEM_DATA* problem, FIELD_DATA* field, MAP* field_info);


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the problem's data.

  \param problem       (o) problem object to be initialized
  \param argc          (i) number of command line arguments
  \param argv          (i) command line arguments

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_problem_data(PROBLEM_DATA* problem, INT argc, CHAR** argv);


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether a given result group belongs to this field.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT match_field_result(FIELD_DATA *field, MAP *result_group);


/*----------------------------------------------------------------------*/
/*!
  \brief Initialize the result data.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void init_result_data(FIELD_DATA* field, RESULT_DATA* result);


/*----------------------------------------------------------------------*/
/*!
  \brief Cleanup result data.

  Please note that there must not be any chunk data on this result
  after this function has been called.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_result_data(RESULT_DATA* result);


/*----------------------------------------------------------------------*/
/*!
  \brief Go to the next result of this discretization.

  \param result (i/o) on input the current result data

  \return zero if there is no further result. non-zero otherwise.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
INT next_result(RESULT_DATA* result);


/*----------------------------------------------------------------------*/
/*!
  \brief Set up the chunk structure to iterate the chunk's entries.

  \param result (i) on input the current result data
  \param chunk  (o) the chunk to be initialized

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void init_chunk_data(RESULT_DATA* result, CHUNK_DATA* chunk, CHAR* name);


/*----------------------------------------------------------------------*/
/*!
  \brief Free the chunk data.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void destroy_chunk_data(CHUNK_DATA* chunk);


/*----------------------------------------------------------------------*/
/*!
  \brief Read one size entry form the file and store it to this
  chunk_data's internal buffer.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void chunk_read_size_entry(CHUNK_DATA* chunk, INT id);


/*----------------------------------------------------------------------*/
/*!
  \brief Read one value entry form the file and store it to this
  chunk_data's internal buffer.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void chunk_read_value_entry(CHUNK_DATA* chunk, INT id);


/*----------------------------------------------------------------------*/
/*!
  \brief Read the element parameters common to all elements.

  This is just a service to make life easier and yet to do extensive
  error checking.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void get_element_params(FIELD_DATA* field,
                        INT i,
                        INT* Id, INT* el_type, INT* dis, INT* numnp);


#ifdef D_FSI

/*----------------------------------------------------------------------*/
/*!
  \brief Find the connection between ale and fluid nodes.

  \param problem       (i) The problem data
  \param struct_field  (i) struct field data; might be NULL
  \param fluid_field   (i) fluid field data
  \param ale_field     (i) ale field data
  \param _fluid_struct_connect  (o) per fluid node array that gives
                                    the corresponding structure nodes
  \param _fluid_ale_connect     (o) gives the corresponding ale nodes

  \warning The connection array are allocated here but must be freed
  by the caller.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void post_find_fsi_coupling(PROBLEM_DATA *problem,
                            FIELD_DATA *struct_field,
                            FIELD_DATA *fluid_field,
                            FIELD_DATA *ale_field,
                            INT **_fluid_struct_connect,
                            INT **_fluid_ale_connect);

#endif



/*----------------------------------------------------------------------*/
/*!
  \brief Nodes and elements of one discretization.

  There are filters that have to read data at once because they need
  to work on nodes and elements. (like post_visual) These filters
  benefit from the real node and element structure so here it is. We
  use the ordinary ccarat nodes and elements, but be careful. Only the
  basic things will be set up. Most values are meaningless here and
  won't be initialized.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
typedef struct _POST_DISCRETIZATION {

  /* general information about this discretization */
  FIELD_DATA* field;

  /* element array */
  ELEMENT* element;

  /* node array */
  NODE* node;

} POST_DISCRETIZATION;


/*----------------------------------------------------------------------*/
/*!
  \brief Set up a (fake) discretization.

  Create the node and element arrays, read node coordinates and mesh
  connectivity.

  \param discret       (o) Uninitialized discretization object
  \param problem       (i) problem data
  \param field         (i) general field data

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void init_post_discretization(POST_DISCRETIZATION* discret,
                              PROBLEM_DATA* problem,
                              FIELD_DATA* field);


#endif
