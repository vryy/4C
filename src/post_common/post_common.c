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
struct _CCA_TRACE         trace;
#endif


/*----------------------------------------------------------------------*/
/*!
  \brief All fields names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
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
    strcpy(allfiles.outputfile_name+length-8, ".post.log");
    strcpy(control_file_name, output_name);
    strcpy(basename, output_name);
    basename[length-8] = '\0';
  }
  else {
    /* dsassert isn't working yet. */
    assert(length+13 < 100);
    strcpy(allfiles.outputfile_name, output_name);
    strcpy(allfiles.outputfile_name+length, ".post.log");
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
  SYMBOL* variant_symbol;
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

  /*--------------------------------------------------------------------*/
  /* open the data files */

  filename = map_read_string(field_info, "value_file");
  field->value_file = fopen(filename, "rb");
  filename = map_read_string(field_info, "size_file");
  field->size_file = fopen(filename, "rb");

  /*--------------------------------------------------------------------*/
  /* find the discretization's name */

  filename = map_read_string(field_info, "value_file");
  i = strlen(filename);
  strcpy(field->name, filename);
  if ((i > 7) && (strcmp(&(filename[i-7]), ".values")==0)) {
    field->name[i-7] = '\0';
  }

  /*--------------------------------------------------------------------*/
  /* read the major/minor conversion */

  /* invalide all entries */
  for (i=0; i<el_count; ++i) {
    INT j;
    for (j=0; j<MAX_EL_MINOR; ++j) {
      field->internal_minors[i][j] = -1;
    }
  }

  variant_symbol = map_find_symbol(field_info, "element_variant");
  if (variant_symbol == NULL) {
    dserror("There must be at least one variant group");
  }

  while (variant_symbol != NULL) {
    INT file_major;
    INT file_minor;
    INT internal_major;
    INT internal_minor;
    MAP* variant_group;

    variant_group = symbol_map(variant_symbol);
    file_major = map_read_int(variant_group, "major");
    file_minor = map_read_int(variant_group, "minor");

    dsassert(file_major < el_count, "detected illegal el_count shrinkage");
    dsassert(file_minor < MAX_EL_MINOR, "detected illegal MAX_EL_MINOR shrinkage");

    restore_element_type(variant_group, &internal_major, &internal_minor);

    field->internal_majors[file_major] = internal_major;
    field->internal_minors[file_major][file_minor] = internal_minor;

    /* on to the next symbol */
    variant_symbol = variant_symbol->next;
  }

  /*--------------------------------------------------------------------*/
  /* figure out the elements used in this discretization */

  fprintf(allfiles.out_err, "%s: Read element types\n", field->name);
  mesh_group = map_read_map(field_info, "mesh");
  size_entry_length = map_read_int(mesh_group, "size_entry_length");
  size_offset = map_read_int(mesh_group, "size_offset");
  mesh_entry = (INT*)CCACALLOC(size_entry_length, sizeof(INT));

  fseek(field->size_file, size_offset, SEEK_SET);
  for (i=0; i<field->numele; ++i) {
    INT file_major;
    INT file_minor;
    INT internal_major;
    INT internal_minor;

    if (fread(mesh_entry, sizeof(INT), size_entry_length, field->size_file)!=size_entry_length) {
      dserror("reading mesh of discretization %s failed", field->name);
    }

    file_major = mesh_entry[1];
    file_minor = mesh_entry[2];

    if ((file_major <= 0) && (file_major >= el_count)) {
      dserror("file major element number out of range");
    }
    if ((file_minor < 0) && (file_minor >= MAX_EL_MINOR)) {
      dserror("file minor element number out of range");
    }

    internal_major = field->internal_majors[file_major];
    internal_minor = field->internal_minors[file_major][file_minor];

    if ((internal_major <= 0) && (internal_major >= el_count)) {
      dserror("internal major element number out of range");
    }
    if ((internal_minor < 0) && (internal_minor >= MAX_EL_MINOR)) {
      dserror("internal minor element number out of range");
    }

    field->element_flags[internal_major][internal_minor] = 1;
    field->element_type[i].Id = mesh_entry[0];
    field->element_type[i].major = internal_major;
    field->element_type[i].minor = internal_minor;
  }

  CCAFREE(mesh_entry);

  /*--------------------------------------------------------------------*/
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


  /*--------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether a given result group belongs to this field.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
INT match_field_result(FIELD_DATA *field, MAP *result_group)
{
  return (strcmp(map_read_string(result_group, "field"),
                 fieldnames[field->type]) == 0) &&
    (map_read_int(result_group, "field_pos") == field->field_pos) &&
    (map_read_int(result_group, "discretization") == field->disnum);
}


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
INT read_chunk_group(CHUNK_DATA* chunk, MAP* result_group, CHAR* name)
{
  if (map_has_map(result_group, name)) {
    MAP* chunk_group;

    chunk_group = map_read_map(result_group, name);
    chunk->group = chunk_group;

    chunk->value_entry_length = map_read_int(chunk_group, "value_entry_length");
    chunk->value_offset = map_read_int(chunk_group, "value_offset");

    chunk->size_entry_length = map_read_int(chunk_group, "size_entry_length");
    chunk->size_offset = map_read_int(chunk_group, "size_offset");

    return 1;
  }
  else {
    return 0;
  }
}


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
                            INT **_fluid_ale_connect)
{
  INT i;
  DOUBLE *struct_coords;
  DOUBLE *ale_coords;
  CHUNK_DATA ale_coords_chunk;
  CHUNK_DATA fluid_coords_chunk;
  CHUNK_DATA struct_coords_chunk;

  INT *fluid_struct_connect = NULL;
  INT *fluid_ale_connect;
  DOUBLE coords[3];

  /* read all struct node coordinates */
  if (struct_field != NULL) {
    if (!read_chunk_group(&struct_coords_chunk, struct_field->table, "coords")) {
      dserror("no coords in struct discretization");
    }
    dsassert(struct_coords_chunk.value_entry_length == problem->ndim,
             "wrong dimension number in struct field");

    struct_coords = (DOUBLE*)CCACALLOC(struct_coords_chunk.value_entry_length*struct_field->numnp,
                                       sizeof(DOUBLE));

    fseek(struct_field->value_file, struct_coords_chunk.value_offset, SEEK_SET);
    if (fread(struct_coords, sizeof(DOUBLE),
              struct_coords_chunk.value_entry_length*struct_field->numnp,
              struct_field->value_file) != struct_coords_chunk.value_entry_length*struct_field->numnp) {
      dserror("reading coordinates of struct field failed");
    }
  }

  /* read all ale node coordinates */
  if (!read_chunk_group(&ale_coords_chunk, ale_field->table, "coords")) {
    dserror("no coords in ale discretization");
  }
  dsassert(ale_coords_chunk.value_entry_length == problem->ndim,
           "wrong dimension number in ale field");

  ale_coords = (DOUBLE*)CCACALLOC(ale_coords_chunk.value_entry_length*ale_field->numnp,
                                  sizeof(DOUBLE));

  fseek(ale_field->value_file, ale_coords_chunk.value_offset, SEEK_SET);
  if (fread(ale_coords, sizeof(DOUBLE),
            ale_coords_chunk.value_entry_length*ale_field->numnp,
            ale_field->value_file) != ale_coords_chunk.value_entry_length*ale_field->numnp) {
    dserror("reading coordinates of ale field failed");
  }

  /* read the fluid node coordinates one by one and search for
   * matching struct and ale nodes. */
  if (!read_chunk_group(&fluid_coords_chunk, fluid_field->table, "coords")) {
    dserror("no coords in fluid discretization");
  }
  dsassert(fluid_coords_chunk.value_entry_length == problem->ndim,
           "wrong dimension number in fluid field");

  fluid_ale_connect = (INT*)CCACALLOC(fluid_field->numnp, sizeof(INT));
  if (struct_field != NULL) {
    fluid_struct_connect = (INT*)CCACALLOC(fluid_field->numnp, sizeof(INT));
  }

  /* This is a quadratic loop. If it turns out to be too slow one
   * could implement some quad- or octtree algorithm. */
  fseek(fluid_field->value_file, fluid_coords_chunk.value_offset, SEEK_SET);
  for (i=0; i<fluid_field->numnp; ++i) {
    INT n_ale;
    if (fread(coords, sizeof(DOUBLE),
              fluid_coords_chunk.value_entry_length,
              fluid_field->value_file)!=fluid_coords_chunk.value_entry_length) {
      dserror("reading coordinates of fluid field failed");
    }

    /* search the structure nodes */
    if (struct_field != NULL) {
      INT n_struct;

      /* no corresponding struct node by default */
      fluid_struct_connect[i] = -1;

      /* search the struct node */
      for (n_struct=0; n_struct<struct_field->numnp; ++n_struct) {
        INT k;
        DOUBLE diff = 0;

        /* quadratic error norm */
        for (k=0; k<problem->ndim; ++k) {
          DOUBLE d = coords[k] - struct_coords[problem->ndim*n_struct+k];
          diff += d*d;
        }

        /*
         * If the difference is very small we've found the corresponding
         * node. The tolerance here might be too big for very fine
         * meshes. I don't know. (ccarat uses the same tolerance but
         * applies it to the sqare root. I don't want to waste the time
         * to calculate it...)
         *
         * In fluid_struct_connect we store the local indices, that is
         * no real ids. */
        if (diff < 1e-10) {
          fluid_struct_connect[i] = n_struct;
          break;
        }
      }
    }

    /* search the ale nodes */

    /* no corresponding ale node by default */
    fluid_ale_connect[i] = -1;

    /* search the ale node */
    for (n_ale=0; n_ale<ale_field->numnp; ++n_ale) {
      INT k;
      DOUBLE diff = 0;

      /* quadratic error norm */
      for (k=0; k<problem->ndim; ++k) {
        DOUBLE d = coords[k] - ale_coords[problem->ndim*n_ale+k];
        diff += d*d;
      }

      /*
       * If the difference is very small we've found the corresponding
       * node. The tolerance here might be too big for very fine
       * meshes. I don't know. (ccarat uses the same tolerance but
       * applies it to the sqare root. I don't want to waste the time
       * to calculate it...)
       *
       * In fluid_ale_connect we store the local indices, that is no
       * real ids. */
      if (diff < 1e-10) {
        fluid_ale_connect[i] = n_ale;
        break;
      }
    }
  }

  *_fluid_struct_connect = fluid_struct_connect;
  *_fluid_ale_connect = fluid_ale_connect;

  CCAFREE(ale_coords);
  if (struct_field != NULL) {
    CCAFREE(struct_coords);
  }
}

#endif



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
                              FIELD_DATA* field)
{
  CHUNK_DATA mesh_chunk;
  CHUNK_DATA coords_chunk;
  INT i;
  INT *mesh_entry;

  discret->field = field;

  discret->node = (NODE*)CCACALLOC(field->numnp, sizeof(NODE));
  discret->element = (ELEMENT*)CCACALLOC(field->numele, sizeof(ELEMENT));

  /*--------------------------------------------------------------------*/
  /* read the node coordinates */

  if (!read_chunk_group(&coords_chunk, field->table, "coords")) {
    dserror("no coords chunk found");
  }

  fseek(field->value_file, coords_chunk.value_offset, SEEK_SET);
  for (i=0; i<field->numnp; ++i) {

    discret->node[i].Id_loc = i;
    discret->node[i].proc = 0;

    if (fread(discret->node[i].x,
              sizeof(DOUBLE),
              coords_chunk.value_entry_length,
              field->value_file) != coords_chunk.value_entry_length) {
      dserror("reading value file of discretization %s failed", field->name);
    }
  }

  /* the global Id is read from the integer part */
  fseek(field->size_file, coords_chunk.size_offset, SEEK_SET);
  for (i=0; i<field->numnp; ++i) {

    if (fread(&(discret->node[i].Id),
              sizeof(INT),
              coords_chunk.size_entry_length,
              field->size_file) != coords_chunk.size_entry_length) {
      dserror("reading size file of discretization %s failed", field->name);
    }
  }

  /*--------------------------------------------------------------------*/
  /* read the mesh */

  if (!read_chunk_group(&mesh_chunk, field->table, "mesh")) {
    dserror("no mesh chunk found");
  }

  mesh_entry = (INT*)CCACALLOC(mesh_chunk.size_entry_length, sizeof(INT));

  fseek(field->size_file, mesh_chunk.size_offset, SEEK_SET);
  for (i=0; i<field->numele; ++i) {
    INT major;
    INT minor;
    INT numnp;
    INT j;

    if (fread(mesh_entry,
              sizeof(INT),
              mesh_chunk.size_entry_length,
              field->size_file) != mesh_chunk.size_entry_length) {
      dserror("reading size file of discretization %s failed", field->name);
    }

    major = mesh_entry[1];
    minor = mesh_entry[2];

    /* convert from file major/minor to internal major/minor */
    minor = field->internal_minors[major][minor];
    major = field->internal_majors[major];

    /* Now setup the connectivity again. This is element type
     * dependent. */
    numnp = element_info[major].variant[minor].node_number;
    dsassert(mesh_chunk.size_entry_length >= numnp+3, "mesh chunk size too short");

    discret->element[i].Id = mesh_entry[0];
    discret->element[i].Id_loc = i;
    discret->element[i].proc = 0;
    discret->element[i].node = (NODE**)CCACALLOC(numnp, sizeof(NODE*));
    discret->element[i].numnp = numnp;
    discret->element[i].eltyp = major;
    discret->element[i].distyp = element_info[major].variant[minor].dis_type;
    for (j=0; j<numnp; ++j) {
      dsassert(mesh_entry[j+3] < field->numnp, "mesh connectivity corrupt");
      discret->element[i].node[j] = &(discret->node[mesh_entry[j+3]]);
    }
  }

  CCAFREE(mesh_entry);
}
