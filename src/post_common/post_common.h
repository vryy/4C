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
#include "../io/bin_packing.h"

extern CHAR* fieldnames[];

/*----------------------------------------------------------------------*/
/*!
  \brief Element type specification.

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

  /* The knowledge which types of elements are used */
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


#endif
