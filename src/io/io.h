/*!
\file
\brief General head file of the io module.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Include this file if you want to use the io module. It contains the
public interface. That is here are the main structures and the
functions that must be called from elsewhere.

\author u.kue
\date 10/04

*/

#ifdef BINIO

#ifndef IO_H
#define IO_H

#include "../pss_full/pss_table.h"


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
  \brief The central structure used for output.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_OUT_MAIN {
  CHAR name[256];
  FILE* control_file;

  /* number of time steps to write to one file */
  INT steps_per_file;

  /* central registration for output field contexts */
  struct _BIN_OUT_FIELD* fields[MAXFIELD];

} BIN_OUT_MAIN;


/*----------------------------------------------------------------------*/
/*!
  \brief The central structure used for input.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_IN_MAIN {
  MAP table;
} BIN_IN_MAIN;


/*----------------------------------------------------------------------*/
/*!
  \brief A pair of data files.

  To each chunk there is a double and an integer part. Thus there are
  two files needed. Each file comes with an offset, that's the point
  where the next chunk is going to be written.

  The point is that we want to write restart and result information to
  different files but still want to use one mechanism.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_DATA_FILES {

/* the files to write to */
#ifdef PARALLEL
  MPI_File value_file;
  MPI_File size_file;
#else
  FILE* value_file;
  FILE* size_file;
#endif

  /* the position where to write new chunks */
  INT value_file_offset;
  INT size_file_offset;

} BIN_DATA_FILES;


/*----------------------------------------------------------------------*/
/*!
  \brief The structure that is used to write one field (discretization).

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_OUT_FIELD {

  /* The solver attached to this discretization. This allows us to
   * store distributed vectors. */
  SPARSE_TYP* sysarray_typ;
  SPARSE_ARRAY* sysarray;

  /* identify me */
  FIELD *actfield;
  PARTITION *actpart;
  INTRA *actintra;
  INT disnum;

  INT result_count;             /* number of results written so far */
  INT restart_count;            /* number of restarts written so far */

#ifdef PARALLEL
  /* the numbers of nodes to send and receive in this field */
  INT* send_numnp;
  INT* recv_numnp;

  /* the numbers of elements to send and receive in this field */
  INT* send_numele;
  INT* recv_numele;

  /* the numbers of dofs to send and receive in this field */
  INT* send_numdof;
  INT* recv_numdof;
#endif

  /* flags that show what kinds of elements there are */
  INT element_flag[el_count];

  /* the biggest node array sizes in this discretization */
  INT max_size[4];

  BIN_DATA_FILES out_result;
  BIN_DATA_FILES out_restart;
  BIN_DATA_FILES* out;

#ifdef D_SHELL8
  INT is_shell8_problem;
#endif

#ifdef D_SHELL9
  INT is_shell9_problem;
  INT s9_layers;                /* number of layers per element */
  INT s9_numnp;                 /* number of nodes per element */
#endif

} BIN_OUT_FIELD;



/*----------------------------------------------------------------------*/
/*!
  \brief The type this chunk stores.

  We can store values that live in nodes or values that live in
  elements. A third possibility is the storage of distributed vectors.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
typedef enum _CHUNK_TYPE {
  chunk_node,
  chunk_element,
  chunk_dist_vec
} CHUNK_TYPE;

#define CHUNK_TYPE_NAMES { "node","element","dist_vec", NULL }



/*----------------------------------------------------------------------*/
/*!
  \brief The structure that is used to read one field (discretization).

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef struct _BIN_IN_FIELD {

  SPARSE_TYP *sysarray_typ;
  SPARSE_ARRAY *sysarray;

  /* identify me */
  FIELD *actfield;
  PARTITION *actpart;
  INTRA *actintra;
  INT disnum;

#ifdef PARALLEL
  /* For all types of chunks we can read we need to know which
   * items have to go to which processors. Therefore we need the
   * number of items to be send and the item ids themselves. */

  INT* send_numnp;
  INT* recv_numnp;

  INT* send_numele;
  INT* recv_numele;

  INT* send_numdof;
  INT* recv_numdof;

  INT** send_node_ids;
  INT** recv_node_ids;

  INT** send_element_ids;
  INT** recv_element_ids;

  /* dof numbers are ids by themselves. But we stick to the naming
   * scheme. */
  INT** send_dof_ids;
  INT** recv_dof_ids;
#endif

#ifdef PARALLEL
  MPI_File value_file;
  MPI_File size_file;
#else
  FILE* value_file;
  FILE* size_file;
#endif

  MAP* field_info;
} BIN_IN_FIELD;



/*----------------------------------------------------------------------*/
/*!
  \brief Init the main (static) data structure that is needed for
  writing.

  Open the control file and write the first lines. To open the file
  its name must be known. In case of restart we have to adjust it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_out_main(CHAR* outputname);


/*----------------------------------------------------------------------*/
/*!
  \brief Init the main (static) data structure that is needed for
  reading.

  Read the control file and put its content in a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_in_main(CHAR* inputname);


/*----------------------------------------------------------------------*/
/*!
  \brief Init the data structure to write the arrays of one
  disctretization.

  Initialize a \a context variable that can be used to write restart
  data and results of one discretization. Write some general
  information about the discretization as well as node coordinates and
  connectivity.

  You need to call this once before you can output any results,
  however when you call it the discretization must be set up
  already. The node arrays in particular must have their final sizes.

  \param context      (o) pointer to a unoccupied context variable.
  \param sysarray_typ (i) type of system matrix. might be NULL.
  \param sysarray     (i) the matrx itself. might be NULL, too.
  \param actfield     (i) the field
  \param actpart      (i) the partition
  \param actintra     (i) the communicator
  \param disnum       (i) the discretization number

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_out_field(BIN_OUT_FIELD* context,
                        SPARSE_TYP* sysarray_typ,
                        SPARSE_ARRAY* sysarray,
                        FIELD *actfield,
                        PARTITION *actpart,
                        INTRA *actintra,
                        INT disnum);


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The output field context's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_out_field(BIN_OUT_FIELD* context);


/*----------------------------------------------------------------------*/
/*!
  \brief Init the data struture to read node arrays from one
  particular field.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_in_field(BIN_IN_FIELD* context,
                       SPARSE_TYP *sysarray_typ,
                       SPARSE_ARRAY *sysarray,
                       FIELD *actfield,
                       PARTITION *actpart,
                       INTRA *actintra,
                       INT disnum);


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The input field context's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_in_field(BIN_IN_FIELD* data);



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


/*----------------------------------------------------------------------*/
/*!
  \brief Write structure restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_nlnstructdyn(struct _BIN_OUT_FIELD *context,
                                    STRUCT_DYNAMIC  *sdyn,
                                    STRUCT_DYN_CALC *dynvar,
                                    INT nrhs,  DIST_VECTOR *rhs,
                                    INT nsol,  DIST_VECTOR *sol,
                                    INT ndis,  DIST_VECTOR *dispi,
                                    INT nvel,  DIST_VECTOR *vel,
                                    INT nacc,  DIST_VECTOR *acc,
                                    INT nfie,  DIST_VECTOR *fie,
                                    INT nwork, DIST_VECTOR *work);


/*----------------------------------------------------------------------*/
/*!
  \brief Read structure restart data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_nlnstructdyn(STRUCT_DYNAMIC  *sdyn,
                                   STRUCT_DYN_CALC *dynvar,
                                   SPARSE_TYP      *sysarray_typ,
                                   SPARSE_ARRAY    *sysarray,
                                   FIELD           *actfield,
                                   PARTITION       *actpart,
                                   INT             disnum,
                                   INTRA           *actintra,
                                   INT nrhs,  DIST_VECTOR *rhs,
                                   INT nsol,  DIST_VECTOR *sol,
                                   INT ndis,  DIST_VECTOR *dispi,
                                   INT nvel,  DIST_VECTOR *vel,
                                   INT nacc,  DIST_VECTOR *acc,
                                   INT nfie,  DIST_VECTOR *fie,
                                   INT nwork, DIST_VECTOR *work,
                                   INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write static structure restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_nlnstructstat(struct _BIN_OUT_FIELD *context,
                                     STATIC_VAR      *statvar,
                                     STANLN          *nln_data,
                                     INT kstep,
                                     INT nrhs,  DIST_VECTOR *rhs,
                                     INT nsol,  DIST_VECTOR *sol,
                                     INT ndis,  DIST_VECTOR *dispi);


/*----------------------------------------------------------------------*/
/*!
  \brief Read static structure restart data.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_nlnstructstat(STATIC_VAR      *statvar,
                                    STANLN          *nln_data,
                                    SPARSE_TYP      *sysarray_typ,
                                    SPARSE_ARRAY    *sysarray,
                                    FIELD           *actfield,
                                    PARTITION       *actpart,
                                    INT              disnum,
                                    INTRA           *actintra,
                                    INT nrhs,  DIST_VECTOR *rhs,
                                    INT nsol,  DIST_VECTOR *sol,
                                    INT ndis,  DIST_VECTOR *dispi,
                                    INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write fluid restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fluiddyn(struct _BIN_OUT_FIELD   *context,
                                FLUID_DYNAMIC   *fdyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read fluid restart data.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_fluiddyn(FLUID_DYNAMIC   *fdyn,
                               SPARSE_TYP      *sysarray_typ,
                               SPARSE_ARRAY    *sysarray,
                               FIELD           *actfield,
                               PARTITION       *actpart,
                               INT              disnum,
                               INTRA           *actintra,
                               INT              step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write ale restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_aledyn(struct _BIN_OUT_FIELD *context,
                              ALE_DYNAMIC *adyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read ale restart data.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_aledyn(ALE_DYNAMIC     *adyn,
                             SPARSE_TYP      *sysarray_typ,
                             SPARSE_ARRAY    *sysarray,
                             FIELD	     *actfield,
                             PARTITION       *actpart,
                             INT              disnum,
                             INTRA	     *actintra,
                             INT              step);


/*----------------------------------------------------------------------*/
/*!
  \brief Write fsi restart data.

  This is additionally to the individual fields data. So only very
  little is stored here. No field at all.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fsidyn(FSI_DYNAMIC *fsidyn);


/*----------------------------------------------------------------------*/
/*!
  \brief Read fsi restart data.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_fsidyn(FSI_DYNAMIC *fsidyn, INT step);


/*----------------------------------------------------------------------*/
/*!
  \brief Close all output files currently open.

  We don't bother whether there are still writing chunks. This is
  supposed to be called by dserror only.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void io_emergency_close_files();


#endif
#endif
