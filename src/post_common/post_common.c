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

#include <assert.h>

#include "post_common.h"


/* There are some global variables in ccarat that are needed by the
 * service functions. We need to specify them here and set them up
 * properly. */
struct _FILES           allfiles;
struct _PAR     par;

#ifdef DEBUG
struct _TRACE         trace;
#endif

CHAR* fieldnames[] = FIELDNAMES;


/*----------------------------------------------------------------------*/
/*!
  \brief This is ccarat setup in a hurry.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void setup_filter(CHAR* output_name, MAP* control_table, CHAR* basename)
{
  INT length;
  CHAR control_file_name[100];

  par.myrank = 0;
  par.nprocs = 1;

  /* Do we want to have tracing? We could. */
#ifdef DEBUG
  dsinit();

  /* We need to take two steps back. dsinit is too close to ccarat. */
  trace.deepness -= 2;
  trace.trace_on=0;
#endif

  /* The warning system is not set up. It's rather stupid anyway. */

  /* We need to open the error output file. The other ones are not
   * important. */
  length = strlen(output_name);
  if ((length > 8) && (strcmp(output_name+length-8, ".control")==0)) {
    /* dsassert isn't working yet. */
    assert(length-8+13 < 100);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length-8, ".post_gid.log");
    strcpy(control_file_name, output_name);
    strcpy(basename, output_name);
    basename[length-8] = '\0';
  }
  else {
    /* dsassert isn't working yet. */
    assert(length+13 < 100);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length, ".post_gid.log");
    strcpy(control_file_name, output_name);
    strcpy(control_file_name+length, ".control");
    strcpy(basename, output_name);
  }
  allfiles.out_err = fopen(allfiles.outputfile_name, "w");

  parse_control_file(control_table, control_file_name);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract one discretization's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_field_data(FIELD_DATA* field, MAP* field_info)
{
  INT i;
  CHAR* filename;
  CHAR* type;
  MAP* mesh_group;
  MAP* coords_group;
  INT size_entry_length;
  INT size_offset;
  INT* mesh_entry;

  field->table = field_info;
  field->field_pos = map_read_int(field_info, "field_pos");
  field->disnum = map_read_int(field_info, "discretization");
  field->numele = map_read_int(field_info, "numele");
  field->numnp = map_read_int(field_info, "numnp");
  field->numdf = map_read_int(field_info, "numdof");

  field->element_type = (ELEMENT_TYPE*)CCACALLOC(field->numele, sizeof(ELEMENT_TYPE));
  field->node_ids = (INT*)CCACALLOC(field->numnp, sizeof(INT));

  type = map_read_string(field_info, "field");
  for (i=0; fieldnames[i]!=NULL; ++i) {
    if (strcmp(fieldnames[i], type)==0) {
      field->type = i;
      break;
    }
  }
  if (fieldnames[i]==NULL) {
    dserror("unknown field type '%s'", type);
  }

  /* open the data files */
  filename = map_read_string(field_info, "value_file");
  field->value_file = fopen(filename, "r");
  filename = map_read_string(field_info, "size_file");
  field->size_file = fopen(filename, "r");

  /* find the discretization's name */
  filename = map_read_string(field_info, "value_file");
  i = strlen(filename);
  strcpy(field->name, filename);
  if ((i > 7) && (strcmp(&(filename[i-7]), ".values")==0)) {
    field->name[i-7] = '\0';
  }

  /* figure out the elements used in this discretization */
  fprintf(allfiles.out_err, "%s: Read element types\n", field->name);
  mesh_group = map_read_map(field_info, "mesh");
  size_entry_length = map_read_int(mesh_group, "size_entry_length");
  size_offset = map_read_int(mesh_group, "size_offset");
  mesh_entry = (INT*)CCACALLOC(size_entry_length, sizeof(INT));

  fseek(field->size_file, size_offset, SEEK_SET);
  for (i=0; i<field->numele; ++i) {
    INT major;
    INT minor;

    if (fread(mesh_entry, sizeof(INT), size_entry_length, field->size_file)!=size_entry_length) {
      dserror("reading mesh of discretization %s failed", field->name);
    }

    major = mesh_entry[1];
    minor = mesh_entry[2];

    dsassert((major >= 0) && (major < el_count), "major element number out of range");
    dsassert((minor >= 0) && (minor < MAX_EL_MINOR), "minor element number out of range");

    field->element_flags[major][minor] = 1;
    field->element_type[i].Id = mesh_entry[0];
    field->element_type[i].major = major;
    field->element_type[i].minor = minor;
  }

  CCAFREE(mesh_entry);

  /* read the node ids */
  fprintf(allfiles.out_err, "%s: Read node ids\n", field->name);
  coords_group = map_read_map(field_info, "coords");
  size_entry_length = map_read_int(coords_group, "size_entry_length");
  dsassert(size_entry_length == 1, "there must be one id to each node");
  size_offset = map_read_int(coords_group, "size_offset");

  fseek(field->size_file, size_offset, SEEK_SET);

  if (fread(field->node_ids,
            sizeof(INT),
            field->numnp*size_entry_length,
            field->size_file) != field->numnp*size_entry_length) {
    dserror("reading node ids of discretization %s failed", field->name);
  }


  /* special problems demand special attention. */

#ifdef D_SHELL9
  if (map_symbol_count(field_info, "shell9_smoothed") > 0) {
    /* This is a shell9 problem. There is guaranteed to be just one
     * type of element. The element_type flags are ignored. */

    field->is_shell9_problem = 1;

    field->s9_smooth_results = map_has_string(field_info, "shell9_smoothed", "yes");
    field->s9_minor = map_read_int(field_info, "shell9_minor");
    field->s9_layers = map_read_int(field_info, "shell9_layers");
    field->s9_forcetype = map_read_string(field_info, "shell9_forcetype");
  }
  else {
    field->is_shell9_problem = 0;
  }
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the problem's data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void init_problem_data(PROBLEM_DATA* problem, MAP* control_table)
{
  CHAR* problem_names[] = PROBLEMNAMES;
  CHAR* type;
  INT i;
  SYMBOL* symbol;

  problem->ndim = map_read_int(control_table, "ndim");
  dsassert((problem->ndim == 2) || (problem->ndim == 3), "illegal dimension");

  type = map_read_string(control_table, "problem_type");
  for (i=0; problem_names[i] != NULL; ++i) {
    if (strcmp(type, problem_names[i])==0) {
      problem->type = i;
      break;
    }
  }
  if (problem_names[i] == NULL) {
    dserror("unknown problem type '%s'", type);
  }

  /* We need to output each field separately. */
  problem->num_discr = map_symbol_count(control_table, "field");
  if (problem->num_discr==0) {
    dserror("no field group found");
  }
  problem->discr = (FIELD_DATA*)CCACALLOC(problem->num_discr, sizeof(FIELD_DATA));

  /* find the first field (the last one that has been written) */
  symbol = map_find_symbol(control_table, "field");

  /* read all fields headers, open the data files */
  for (i=0; i<problem->num_discr; ++i) {
    if (!symbol_is_map(symbol)) {
      dserror("failed to get field group");
    }

    init_field_data(&(problem->discr[i]), symbol_map(symbol));

    symbol = symbol->next;
  }

  /*--------------------------------------------------------------------*/
  /* Now collect all result groups */

  problem->num_results = map_symbol_count(control_table, "result");
  if (problem->num_results == 0) {
    dserror("no results found");
  }
  problem->result_group = (MAP**)CCACALLOC(problem->num_results, sizeof(MAP*));

  /* find the first result group */
  symbol = map_find_symbol(control_table, "result");

  /* We rely on the fact that groups are linked in reverse order. */
  /* That is results are written ordered by time step. */
  for (i=problem->num_results-1; i>=0; --i) {
    if (!symbol_is_map(symbol)) {
      dserror("failed to get result group");
    }

    problem->result_group[i] = symbol_map(symbol);

    symbol = symbol->next;
  }
}
