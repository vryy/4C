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

#ifdef PARALLEL
#undef PARALLEL
#endif

#ifdef SPOOLES_PACKAGE
#undef SPOOLES_PACKAGE
#endif

#ifdef PERF
#undef PERF
#endif

#include "../headers/standardtypes.h"
#include "../pss_full/pss_table.h"
#include "../io/io_packing.h"
#include "../global_full/global_element_info.h"


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
  \brief Element type specification.

  This struct stores the internal major and minor numbers, that is the
  converted ones not the ones read from the file.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef struct _ELEMENT_TYPE {
  INT Id;
  INT major;
  INT minor;
} ELEMENT_TYPE;


/*----------------------------------------------------------------------*/
/*!
  \brief The meta data that belongs to one discretization.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef struct _FIELD_DATA {
  FIELDTYP type;

  INT field_pos;
  INT disnum;
  INT numele;
  INT numnp;
  INT numdf;

  FILE* value_file;
  FILE* size_file;

  MAP* table;
  CHAR name[100];

  /* The internal major type numbers that correspond to the major
   * numbers found in the file. We rely on the ccarat never to loose
   * elements. That is ``el_count`` must never shrink.
   *
   * We need this and the next variable in order to maintain
   * readability of our binary files. */
  INT internal_majors[el_count];

  /* The internal minor type numbers that correspond to the major
   * numbers found in the file. */
  ELEMENT_FLAGS internal_minors;

  /* The knowledge which types of elements are used. This table is
   * accessed by internal major/minor numbers. */
  ELEMENT_FLAGS element_flags;

  /*
   * The types of all elements in this discretization. This can became
   * big. But it frees us from reading the mesh over and over again. */
  ELEMENT_TYPE* element_type;

  /* The real ids of all nodes. Another big and convenient array. */
  INT* node_ids;

#ifdef D_SHELL9
  INT is_shell9_problem;
  INT s9_smooth_results;
  INT s9_minor;
  INT s9_layers;
  CHAR* s9_forcetype;
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

  INT num_results;
  MAP** result_group;
} PROBLEM_DATA;



/*----------------------------------------------------------------------*/
/*!
  \brief The data necessary to access one chunk.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
typedef struct _CHUNK_DATA {
  INT value_entry_length;
  INT value_offset;
  INT size_entry_length;
  INT size_offset;
  MAP* group;
} CHUNK_DATA;


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
  \brief Extract one discretization's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_field_data(FIELD_DATA* field, MAP* field_info);


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the problem's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_problem_data(PROBLEM_DATA* problem, MAP* control_table);


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
  \brief Read the control data of one chunk.

  \param chunk         (o) The structure that is filled
  \param result_group  (i) The control file's group that contains the
                           group to be read.
  \param name          (i) The name of the group to be read

  \return Whether a group with the name asked for was found.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT read_chunk_group(CHUNK_DATA* chunk, MAP* result_group, CHAR* name);


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
