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

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "../headers/standardtypes.h"

#include "io_singlefile.h"
#include "io_packing.h"

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


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES           allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;


/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for output.

  This structure needs to be initialized at startup. The whole output
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
BIN_OUT_MAIN bin_out_main;


/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for input.

  This structure needs to be initialized at startup. The whole input
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
BIN_IN_MAIN bin_in_main;


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
  \brief All node arrays names.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static CHAR* nodearraynames[] = NODEARRAYNAMES;


/*----------------------------------------------------------------------*/
/*!
  \brief Init the main (static) data structure that is needed for
  writing.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_out_main(CHAR* outputname)
{
#ifdef DEBUG
  dstrc_enter("init_bin_out_main");
#endif

  /* If this is a restarted run we don't want to overwrite the results
   * already written. Instead we try to find a new file name. */
  if (genprob.restart) {
    if (par.myrank == 0) {
      INT id;
      CHAR tmpbuf[256];
      INT len;
      INT number = 0;
      strcpy(tmpbuf, outputname);

      /* remove any trailing number */
      /* Hmm. This might hurt some people. Do we want this? */
      for (len = strlen(tmpbuf)-1; (len > 0) && isdigit(tmpbuf[len]); --len) {}
      if ((len < strlen(tmpbuf)-1) && (len > 0) && (tmpbuf[len] == '-')) {
        tmpbuf[len] = '\0';
        number = atoi(&(tmpbuf[len+1]));
      }

      for (;;) {
        number += 1;
        sprintf(bin_out_main.name, "%s-%d.control", tmpbuf, number);
        if ((id = open(bin_out_main.name, O_CREAT|O_EXCL, 0660)) != -1) {
          close(id);
          bin_out_main.name[strlen(bin_out_main.name)-8] = '\0';
          break;
        }
      }
    }
#ifdef PARALLEL
    if (par.nprocs > 1) {
      INT err;
      err = MPI_Bcast(bin_out_main.name, sizeof(bin_out_main.name), MPI_CHAR, 0, MPI_COMM_WORLD);
      if (err != 0) {
        dserror("MPI_Bcast failed: %d", err);
      }
    }
#endif
  }
  else {
    strcpy(bin_out_main.name, outputname);
  }

  if (par.myrank == 0) {
    static CHAR* problem_names[] = PROBLEMNAMES;
    CHAR name[256];
    sprintf(name, "%s.control", bin_out_main.name);
    bin_out_main.control_file = fopen(name, "w");

    fprintf(bin_out_main.control_file, "# ccarat output control file\n"
            "version = \"0.1\"\n"
            "problem_type = \"%s\"\n"
            "ndim = %d\n"
            "\n",
            problem_names[genprob.probtyp],
            genprob.ndim);

    fflush(bin_out_main.control_file);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Init the main (static) data structure that is needed for
  reading.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_in_main(CHAR* inputname)
{
  CHAR name[256];

#ifdef DEBUG
  dstrc_enter("init_bin_in_main");
#endif

  /* The control file is read on all processors. */
  sprintf(name, "%s.control", inputname);
  parse_control_file(&(bin_in_main.table), name);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The output field context's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_out_field(BIN_OUT_FIELD* context)
{
#ifdef PARALLEL
  INT err;
#endif

#ifdef DEBUG
  dstrc_enter("destroy_bin_out_field");
#endif

#ifdef PARALLEL

  CCAFREE(context->send_numnp);
  CCAFREE(context->recv_numnp);

  CCAFREE(context->send_numele);
  CCAFREE(context->recv_numele);

  if (context->sysarray_typ != NULL) {
    CCAFREE(context->send_numdof);
    CCAFREE(context->recv_numdof);
  }

  err = MPI_File_close(&(context->value_file));
  if (err != 0) {
    dserror("MPI_File_close failed: %d", err);
  }

  err = MPI_File_close(&(context->size_file));
  if (err != 0) {
    dserror("MPI_File_close failed: %d", err);
  }
#else
  fclose(context->value_file);
  fclose(context->size_file);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The output chunk destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_out_chunk(BIN_OUT_CHUNK* data)
{
#ifdef DEBUG
  dstrc_enter("destroy_bin_out_chunck");
#endif

  if (data->out_values != NULL)
    CCAFREE(data->out_values);

  if (data->out_sizes != NULL)
    CCAFREE(data->out_sizes);

  destroy_map(&(data->group));

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The input field context's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_in_field(BIN_IN_FIELD* context)
{
#ifdef PARALLEL
  INT i;
  INT err;
  INT nproc = context->actintra->intra_nprocs;
#endif

#ifdef DEBUG
  dstrc_enter("destroy_bin_in_field");
#endif

#ifdef PARALLEL
  for (i=0; i<nproc; ++i) {
    CCAFREE(context->send_node_ids[i]);
    CCAFREE(context->recv_node_ids[i]);

    CCAFREE(context->send_element_ids[i]);
    CCAFREE(context->recv_element_ids[i]);
  }

  CCAFREE(context->send_node_ids);
  CCAFREE(context->recv_node_ids);

  CCAFREE(context->send_element_ids);
  CCAFREE(context->recv_element_ids);

  CCAFREE(context->send_numnp);
  CCAFREE(context->recv_numnp);

  CCAFREE(context->send_numele);
  CCAFREE(context->recv_numele);

  if (context->sysarray_typ != NULL) {
    CCAFREE(context->send_dof_ids);
    CCAFREE(context->recv_dof_ids);

    CCAFREE(context->send_numdof);
    CCAFREE(context->recv_numdof);
  }

  err = MPI_File_close(&(context->value_file));
  if (err != 0) {
    dserror("MPI_File_close failed: %d", err);
  }

  err = MPI_File_close(&(context->size_file));
  if (err != 0) {
    dserror("MPI_File_close failed: %d", err);
  }
#else
  fclose(context->value_file);
  fclose(context->size_file);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  The input chunk's destructor.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_bin_in_chunk(BIN_IN_CHUNK* chunk)
{
#ifdef DEBUG
  dstrc_enter("destroy_bin_in_chunck");
#endif

  if (chunk->value_count > 0)
    CCAFREE(chunk->in_values);

  if (chunk->size_count > 0)
    CCAFREE(chunk->in_sizes);

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief Find the number dofs that need to be send and received.

  This is solver dependent.

  The main duty of this function is needed in the parallel version
  only. However it also finds the number of dofs, and we want to know
  it in the sequential version, too.

  \return the number of dofs.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static INT out_setup_dof_transfer(BIN_OUT_FIELD *context,
                                  SPARSE_TYP *sysarray_typ,
                                  SPARSE_ARRAY *sysarray,
                                  INTRA *actintra,
                                  INT nprocs,
                                  INT rank)
{
  INT i, max_num, fullarrays;
  INT tag_base = 0xABCD;
  INT ndof = 0;

#ifdef DEBUG
  dstrc_enter("out_setup_dof_transfer");
#endif


#define calc_max_number(numeq_total)                            \
  ndof = numeq_total;                                           \
  max_num = ((numeq_total) + nprocs - 1) / nprocs;              \
  fullarrays = nprocs - (nprocs*max_num - (numeq_total));

#ifdef PARALLEL
#define boilerplate_code                                        \
  {                                                             \
    INT proc;                                                   \
                                                                \
    /* Find the processor this dof has to be sent to. */        \
    if (dof < max_num*fullarrays) {                             \
      proc = dof / max_num;                                     \
    }                                                           \
    else {                                                      \
      INT Id = dof - max_num*fullarrays;                        \
      proc = Id / (max_num-1) + fullarrays;                     \
    }                                                           \
                                                                \
    /* Count for sending. */                                    \
    context->send_numdof[proc]++;                               \
  }

#else
#define boilerplate_code
#endif

  /* Find out how many dofs we have to send to which processor. */
  switch (*sysarray_typ) {

#ifdef AZTEC_PACKAGE
  case msr:
    calc_max_number(sysarray->msr->numeq_total);
    for (i=0; i<sysarray->msr->numeq; ++i) {
      INT dof = sysarray->msr->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

#ifdef HYPRE_PACKAGE
  case parcsr: {
    INT rank = chunk->field->actintra->intra_rank;
    calc_max_number(sysarray->parcsr->numeq_total);
    for (i=0; i<sysarray->parcsr->numeq; ++i) {
      INT dof = sysarray->parcsr->update.a.ia[rank][i];
      boilerplate_code;
    }
    break;
  }
#endif

#ifdef PARSUPERLU_PACKAGE
  case ucchb:
    calc_max_number(sysarray->ucchb->numeq_total);
    for (i=0; i<sysarray->ucchb->numeq; ++i) {
      INT dof = sysarray->ucchb->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

  case dense:
    calc_max_number(sysarray->dense->numeq_total);
    for (i=0; i<sysarray->dense->numeq; ++i) {
      INT dof = sysarray->dense->update.a.iv[i];
      boilerplate_code;
    }
    break;

#ifdef MLIB_PACKAGE
  case mds:
    calc_max_number(sysarray->mds->numeq_total);
    for (i=0; i<sysarray->mds->numeq; ++i) {
      INT dof = i;
      boilerplate_code;
    }
    break;
#endif

#ifdef MUMPS_PACKAGE
  case rc_ptr:
    calc_max_number(sysarray->rc_ptr->numeq_total);
    for (i=0; i<sysarray->rc_ptr->numeq; ++i) {
      INT dof = sysarray->rc_ptr->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

#ifdef SPOOLES_PACKAGE
  case spoolmatrix:
    calc_max_number(sysarray->spo->numeq_total);
    for (i=0; i<sysarray->spo->numeq; ++i) {
      INT dof = sysarray->spo->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

#ifdef UMFPACK
  case ccf:
    calc_max_number(sysarray->ccf->numeq_total);
    for (i=0; i<sysarray->ccf->numeq; ++i) {
      INT dof = sysarray->ccf->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

  case skymatrix:
    calc_max_number(sysarray->sky->numeq_total);
    for (i=0; i<sysarray->sky->numeq; ++i) {
      INT dof = sysarray->sky->update.a.iv[i];
      boilerplate_code;
    }
    break;

#ifdef MLPCG
  case bdcsr:
    calc_max_number(sysarray->bdcsr->numeq_total);
    for (i=0; i<sysarray->bdcsr->numeq; ++i) {
      INT dof = sysarray->bdcsr->update.a.iv[i];
      boilerplate_code;
    }
    break;
#endif

  case oll:
    calc_max_number(sysarray->oll->numeq_total);
    for (i=0; i<sysarray->oll->numeq; ++i) {
      INT dof = sysarray->oll->update.a.iv[i];
      boilerplate_code;
    }
    break;

  default:
    dserror("Unknown type %d of system matrix", *sysarray_typ);
    break;
  }

#undef boilerplate_code
#undef calc_max_number

#ifdef PARALLEL

  /* The information how many dofs we have to receive can only be
   * found by communication. */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    err = MPI_Sendrecv(&(context->send_numdof[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       &(context->recv_numdof[src]), 1, MPI_INT, src, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
  return ndof;
}


#ifdef PARALLEL


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number elements that need to be send and received.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void out_setup_element_transfer(BIN_OUT_FIELD *context,
                                       INTRA *actintra,
                                       INT nprocs,
                                       INT rank,
                                       DISCRET *actdis,
                                       PARTDISCRET *actpdis)
{
  INT i, max_num, fullarrays, j;
  INT tag_base = 0xABCD;

#ifdef DEBUG
  dstrc_enter("out_setup_element_transfer");
#endif

  max_num = (actdis->numele + nprocs - 1) / nprocs;
  fullarrays = nprocs - (nprocs*max_num - actdis->numele);

  /* Find out how many elements we have to send to which processor. */
  for (j=0; j<actpdis->numele; ++j) {
    ELEMENT* actele = actpdis->element[j];

    INT proc;

    /* Find the processor this element has to be sent to. */
    if (actele->Id_loc < max_num*fullarrays) {
      proc = actele->Id_loc / max_num;
    }
    else {
      INT Id = actele->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }

    /* Count for sending. */
    context->send_numele[proc]++;
  }

  /* The information how many elements we have to receive can only be
   * found by communication. */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    err = MPI_Sendrecv(&(context->send_numele[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       &(context->recv_numele[src]), 1, MPI_INT, src, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number nodes that need to be send and received.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void out_setup_node_transfer(BIN_OUT_FIELD *context,
                                    INTRA *actintra,
                                    INT nprocs,
                                    INT rank,
                                    DISCRET *actdis,
                                    PARTDISCRET *actpdis)
{
  INT i, max_num, fullarrays, j;
  INT tag_base = 0xABCD;

#ifdef DEBUG
  dstrc_enter("out_setup_node_transfer");
#endif

  max_num = (actdis->numnp + nprocs - 1) / nprocs;
  fullarrays = nprocs - (nprocs*max_num - actdis->numnp);

  /* Find out how many nodes we have to send to which processor. */
  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];

    INT proc;

    /* Find the processor this node has to be sent to. */
    if (actnode->Id_loc < max_num*fullarrays) {
      proc = actnode->Id_loc / max_num;
    }
    else {
      INT Id = actnode->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }

    /* Count for sending. */
    context->send_numnp[proc]++;
  }

  /* The information how many nodes we have to receive can only be
   * found by communication.
   *
   * We do it in a chattering fashion here. However, it could be done
   * in a round robin style as well because the amount of data is
   * rather small. If this loop turns out to be a bottleneck we can
   * change it (and hope for improvements.) */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    err = MPI_Sendrecv(&(context->send_numnp[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       &(context->recv_numnp[src]), 1, MPI_INT, src, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


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
                        INT disnum)
{
#ifdef PARALLEL
  INT err;
#endif
  INT nprocs = actintra->intra_nprocs;
  INT rank = actintra->intra_rank;
  INT j;
  DISCRET* actdis = &(actfield->dis[disnum]);
  PARTDISCRET* actpdis = &(actpart->pdis[disnum]);
  INT value_length;
  INT size_length;
  INT ndof = 0;

#ifdef PARALLEL
  INT max_size[4];
#else
  INT* max_size = context->max_size;
#endif

  INT field_pos;
  CHAR filename[256];

#ifdef DEBUG
  dstrc_enter("init_bin_out_field");
#endif

  /* Remember the variables of this field. */
  context->sysarray_typ = sysarray_typ;
  context->sysarray = sysarray;

  context->actfield = actfield;
  context->actpart = actpart;
  context->actintra = actintra;
  context->disnum = disnum;

  context->value_file_offset = 0;
  context->size_file_offset = 0;

  /*--------------------------------------------------------------------*/
  /* add new field to control file */

  if (rank == 0) {
    out_main_group_head(context, "field");
  }

  /*--------------------------------------------------------------------*/
  /* open file to write this field */

  /*
   * The file names must be unique even in (future) problems with
   * different field of the type and many discretizations. Thus the
   * names get a little weird. */

  field_pos = get_field_position(context);

  sprintf(filename, "%s.%s.%d.%d.values", bin_out_main.name,
          fieldnames[actfield->fieldtyp], field_pos, disnum);
#ifdef PARALLEL
  err = MPI_File_open(actintra->MPI_INTRA_COMM, filename,
                      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                      &(context->value_file));
  if (err != 0) {
    dserror("MPI_File_open failed for file '%s': %d", filename, err);
  }
#else
  context->value_file = fopen(filename, "w");
#endif
  if (rank == 0) {
    fprintf(bin_out_main.control_file, "    value_file = \"%s\"\n", filename);
  }

  sprintf(filename, "%s.%s.%d.%d.sizes", bin_out_main.name,
          fieldnames[actfield->fieldtyp], field_pos, disnum);
#ifdef PARALLEL
  err = MPI_File_open(actintra->MPI_INTRA_COMM, filename,
                      MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                      &(context->size_file));
  if (err != 0) {
    dserror("MPI_File_open failed for file '%s': %d", filename, err);
  }
#else
  context->size_file = fopen(filename, "w");
#endif
  if (rank == 0) {
    fprintf(bin_out_main.control_file, "    size_file = \"%s\"\n", filename);
  }

#ifdef PARALLEL

  /*--------------------------------------------------------------------*/
  /* setup communication variables */

  /* the number of nodes we need to send and receive per processor */
  context->send_numnp = (INT*)CCACALLOC(nprocs, sizeof(INT));
  context->recv_numnp = (INT*)CCACALLOC(nprocs, sizeof(INT));

  /* the number of elements we need to send and receive per processor */
  context->send_numele = (INT*)CCACALLOC(nprocs, sizeof(INT));
  context->recv_numele = (INT*)CCACALLOC(nprocs, sizeof(INT));

  if (sysarray_typ != NULL) {
    /* the numbers of dofs we need to send and receive per processor */
    context->send_numdof = (INT*)CCACALLOC(nprocs, sizeof(INT));
    context->recv_numdof = (INT*)CCACALLOC(nprocs, sizeof(INT));
  }

#endif

  max_size[0] = 0;
  max_size[1] = 0;
  max_size[2] = 0;
  max_size[3] = 0;

  for (j=0; j<actpdis->numnp; ++j) {
    NODE* actnode = actpdis->node[j];

    max_size[0] = MAX(max_size[0], actnode->sol.fdim*actnode->sol.sdim);
    max_size[1] = MAX(max_size[1], actnode->sol_increment.fdim*actnode->sol_increment.sdim);
    max_size[2] = MAX(max_size[2], actnode->sol_residual.fdim*actnode->sol_residual.sdim);
    max_size[3] = MAX(max_size[3], actnode->sol_mf.fdim*actnode->sol_mf.sdim);
  }

#ifdef PARALLEL

  /* We store the maximum node array size, too. The idea is to
   * allocate the maximum size for every array in order to be able to
   * access each node array easily. Of couse this results in wasted
   * memory if the variation in size is too big. However, we calculate
   * the sizes for every node array in every discretization. So the
   * waste shouldn't be too devastating. */
  MPI_Allreduce(max_size, context->max_size, 4, MPI_INT, MPI_MAX, actintra->MPI_INTRA_COMM);

  /*--------------------------------------------------------------------*/
  /* find number of nodes that need to be communicated */

  out_setup_node_transfer(context, actintra, nprocs, rank, actdis, actpdis);

  /*--------------------------------------------------------------------*/
  /* find number of elements that need to be communicated */

  out_setup_element_transfer(context, actintra, nprocs, rank, actdis, actpdis);

#endif

  /*--------------------------------------------------------------------*/
  /* find number of dofs that need to be communicated */

  /*
   * this is done in the sequential version, too, to get the number of
   * dofs. */
  if (sysarray_typ != NULL) {
    ndof = out_setup_dof_transfer(context, sysarray_typ, sysarray, actintra, nprocs, rank);
  }

  /*--------------------------------------------------------------------*/
  /* Now check for different types of elements. */

  out_find_element_types(context, actintra, actpdis);

  /*--------------------------------------------------------------------*/
  /* store additional size information */

  if (rank == 0) {
    fprintf(bin_out_main.control_file,
            "    numnp = %d\n"
            "    numele = %d\n"
            "    numdof = %d\n"
            "\n",
            actdis->numnp, actdis->numele, ndof);
  }

  /*
   * Some element specific code has to go here because some elements
   * require special treatment. */
#ifdef D_SHELL8
  {
    INT shell8_count = count_element_variants(context, el_shell8);
    if (shell8_count > 0) {
      ELEMENT* actele;
      SHELL8* s8;

      if (shell8_count != count_element_types(context)) {
        dserror("shell8 is supposed to be the only element in the mesh");
      }
      if (shell8_count != 1) {
        dserror("there must be only one shell8 variant");
      }

      /*
       * OK. We have a shell8 problem.  */
      context->is_shell8_problem = 1;

      actele = context->actpart->pdis[context->disnum].element[0];
      dsassert(actele->eltyp == el_shell8, "shell8 expected");

      s8 = actele->e.s8;
      context->s8_minor = MINOR_SHELL8(actele);

      /*
       * Mark this discretization. The filter will know that it is a
       * special one. */
      if (rank == 0) {
        fprintf(bin_out_main.control_file,
                "    shell8_minor = %d\n"
                "    shell8_scal = %f\n"
                "    shell8_sdc = %f\n"
                "\n",
                context->s8_minor, 1.0, s8->sdc);
      }
    }
    else {
      context->is_shell8_problem = 0;
    }
  }
#endif
#ifdef D_SHELL9
  {
    INT shell9_count = count_element_variants(context, el_shell9);
    if (shell9_count > 0) {
      ELEMENT* actele;
      SHELL9* s9;
      INT klay;

      if (shell9_count != count_element_types(context)) {
        dserror("shell9 is supposed to be the only element in the mesh");
      }
      if (shell9_count != 1) {
        dserror("there must be only one shell9 variant");
      }

      /*
       * OK. We have a shell9 problem.  */
      context->is_shell9_problem = 1;

      actele = context->actpart->pdis[context->disnum].element[0];
      dsassert(actele->eltyp == el_shell9, "shell9 expected");

      s9 = actele->e.s9;
      context->s9_minor = MINOR_SHELL9(actele);

      context->s9_layers = 0;
      for (klay=0; klay<s9->num_klay; klay++) {
        context->s9_layers += s9->kinlay[klay].num_mlay;
      }

      /*
       * Mark this discretization. The filter will know that it is a
       * special one. */
      if (rank == 0) {
        fprintf(bin_out_main.control_file,
                "    shell9_smoothed = \"%s\"\n"
                "    shell9_minor = %d\n"
                "    shell9_layers = %d\n"
                "    shell9_forcetype = \"%s\"\n"
                "\n",
                (ioflags.struct_stress_gid_smo ? "yes" : "no"),
                context->s9_minor, context->s9_layers,
                (s9->forcetyp==s9_xyz) ? "xyz" : ((s9->forcetyp==s9_rst) ? "rst" : "rst_ortho"));
      }
    }
    else {
      context->is_shell9_problem = 0;
    }
  }
#endif

  /*--------------------------------------------------------------------*/
  /* store the mesh connectivity and node coordinates */

  find_mesh_item_length(context, &value_length, &size_length);
  out_element_chunk(context, "mesh", cc_mesh, value_length, size_length, 0);
  out_node_chunk(context, "coords", cc_coords, genprob.ndim, 1, 0);

#ifdef D_SHELL8
  if (context->is_shell8_problem) {
    INT numnp;
    numnp = element_info[el_shell8].variant[context->s8_minor].node_number;

    /* There is a director at each node. Furthermore we need the
     * element thickness. Let's play it save and output one thickness
     * value per node. This is more that is currently supported by
     * ccarat's input though. */
    out_element_chunk(context, "shell8_director", cc_shell8_director, (3+1)*numnp, 0, 0);
  }
#endif

#ifdef D_SHELL9
  if (context->is_shell9_problem) {
    INT numnp;
    INT layer_count;

    /* We need space for three values per artificial node. */

    numnp = element_info[el_shell9].variant[context->s9_minor].node_number;

    if ((context->s9_minor == MINOR_SHELL9_4_22) ||
        (context->s9_minor == MINOR_SHELL9_4_33)) {
      layer_count = 1;
    }
    else {
      layer_count = 2;
    }

    out_element_chunk(context, "shell9_coords", cc_shell9_coords,
                      3*numnp*(layer_count*context->s9_layers+1), 0, 0);
  }
#endif

#ifdef PARALLEL
  /* Have the element distribution there as well. */
  out_element_chunk(context, "domain", cc_domain, 0, 1, 0);
#endif

  if (rank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                        INT size_entry_length)
{
  INT nprocs = context->actintra->intra_nprocs;
  INT rank = context->actintra->intra_rank;
  DISCRET* actdis = &(context->actfield->dis[context->disnum]);
  INT num;

#ifdef DEBUG
  dstrc_enter("init_bin_out_chunk");
#endif

  chunk->field = context;
  init_map(&(chunk->group));

  switch (type) {
  case chunk_node:
    num = actdis->numnp;
    chunk->type = chunk_node;
    break;
  case chunk_element:
    num = actdis->numele;
    chunk->type = chunk_element;
    break;
  case chunk_dist_vec:
    /* We are able to handle an array of distributed vectors. But
     * there must be a first one and all of them must have the same
     * size. */
    num = chunk->vectors[0].numeq_total;
    chunk->type = chunk_dist_vec;
    break;
  default:
    dserror("unknown chunk type %d", type);
  }

  /* Figure out how many entities are written by each processor. */
  chunk->num = (num + nprocs - 1) / nprocs;

  /* The number of processors that write the full count. The other
   * ones have one less to write. */
  chunk->fullarrays = nprocs - (nprocs*chunk->num - num);

  /* The Id_loc of the first node/element of this processor. */
  chunk->first_id = chunk->num*rank;

  if (rank >= chunk->fullarrays) {
    chunk->num -= 1;
    chunk->first_id -= rank - chunk->fullarrays;
  }

  /* remember the items size */

  map_insert_int_cpy(&chunk->group, value_entry_length, "value_entry_length");
  map_insert_int_cpy(&chunk->group, size_entry_length, "size_entry_length");

  /* allocate the buffers that are to be written */

  chunk->value_entry_length = value_entry_length;
  chunk->value_count = value_entry_length*chunk->num;
  if (chunk->value_count > 0) {
    chunk->out_values = (DOUBLE*)CCACALLOC(chunk->value_count, sizeof(DOUBLE));
  }
  else {
    chunk->out_values = NULL;
  }

  chunk->size_entry_length = size_entry_length;
  chunk->size_count = size_entry_length*chunk->num;
  if (chunk->size_count > 0) {
    chunk->out_sizes = (INT*)CCACALLOC(chunk->size_count, sizeof(INT));
  }
  else {
    chunk->out_sizes = NULL;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



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
                       INT array)
{
#ifdef PARALLEL
  INT tag_base = 0xABCD;
#endif
  INT nprocs = context->actintra->intra_nprocs;
  INT rank = context->actintra->intra_rank;
  INT i;
  DISCRET* actdis = &(context->actfield->dis[context->disnum]);
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
  INT num;

#ifdef DEBUG
  dstrc_enter("out_gather_values");
#endif

  switch (chunk->type) {
  case chunk_node:
    num = actdis->numnp;
    break;
  case chunk_element:
    num = actdis->numele;
    break;
  case chunk_dist_vec:
    /* We are able to handle an array of distributed vectors. But
     * there must be a first one and all of them must have the same
     * size. */
    num = chunk->vectors[0].numeq_total;
    break;
  default:
    dserror("unknown chunk type %d", type);
  }

  for (i=0; i<nprocs; ++i) {
#ifdef PARALLEL
    MPI_Status status;
    DOUBLE* recv_buf;
    INT recv_count = 0;
    INT* recv_size_buf;
    INT recv_size_count = 0;
    INT err;
    INT j;
    INT size_len;
    INT value_len;
    INT send_num;
    INT recv_num;
#endif
    INT dst;
    INT src;
    DOUBLE* send_buf;
    INT send_count = 0;
    INT* send_size_buf;
    INT send_size_count = 0;
    INT dst_first_id;
    INT dst_num;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    /* we need to known what nodes our destination expects */
    dst_num = (num + nprocs - 1) / nprocs;
    dst_first_id = dst_num*dst;
    if (dst >= chunk->fullarrays) {
      dst_num -= 1;
      dst_first_id -= dst - chunk->fullarrays;
    }

#ifdef PARALLEL
    switch (chunk->type) {
    case chunk_node:
      send_num = context->send_numnp[dst];
      recv_num = context->recv_numnp[src];
      break;
    case chunk_element:
      send_num = context->send_numele[dst];
      recv_num = context->recv_numele[src];
      break;
    case chunk_dist_vec:
      send_num = context->send_numdof[dst];
      recv_num = context->recv_numdof[src];
      break;
    default:
      dserror("unknown chunk type %d", type);
    }

    if (chunk->value_entry_length > 0) {
      send_count = send_num*chunk->value_entry_length;
      recv_count = recv_num*chunk->value_entry_length;

      send_buf = (DOUBLE*)CCACALLOC(send_count, sizeof(DOUBLE));
      recv_buf = (DOUBLE*)CCACALLOC(recv_count, sizeof(DOUBLE));
    }

    /*
     * We always need to send the item's id in order to identify it on
     * the receiving side. This id is never written, however. */
    send_size_count = (chunk->size_entry_length+1)*send_num;
    recv_size_count = (chunk->size_entry_length+1)*recv_num;

    send_size_buf = (INT*)CCACALLOC(send_size_count, sizeof(INT));
    recv_size_buf = (INT*)CCACALLOC(recv_size_count, sizeof(INT));
#else
    /* In the sequential version we take a shortcut. Here we gather
     * values into the array that's used for writing. Thus no
     * additional copying. However this requires us to send exactly
     * the same data that's written afterwards. In particular we must
     * not fill the sizes array with those ids that are never
     * written. */
    send_count = num*chunk->value_entry_length;
    send_buf = chunk->out_values;

    send_size_count = num*chunk->size_entry_length;
    send_size_buf = chunk->out_sizes;
#endif

    out_pack_items(chunk, type, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);

#ifdef PARALLEL

    /* communication */
    err = MPI_Sendrecv(send_size_buf, send_size_count, MPI_INT, dst, tag_base-i-1,
                       recv_size_buf, recv_size_count, MPI_INT, src, tag_base-i-1,
                       context->actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }

    if (chunk->value_entry_length > 0) {
      err = MPI_Sendrecv(send_buf, send_count, MPI_DOUBLE, dst, tag_base+i,
                         recv_buf, recv_count, MPI_DOUBLE, src, tag_base+i,
                         context->actintra->MPI_INTRA_COMM, &status);
      if (err != 0) {
        dserror("mpi sendrecv error %d", err);
      }
    }

    size_len = chunk->size_entry_length;
    value_len = chunk->value_entry_length;

    /* copy received values into the output buffers */

    for (j=0; j<recv_num; ++j) {
      INT k;
      INT Id;
      DOUBLE *src_ptr;
      DOUBLE *dst_ptr;
      INT *src_int_ptr;
      INT *dst_int_ptr;

      Id = recv_size_buf[(size_len+1)*j];
      dsassert((Id>=chunk->first_id) &&
               (Id<chunk->first_id+chunk->num), "inadequate id");

      /* we don't store the id which is always transfered first */
      src_int_ptr = &(recv_size_buf[(size_len+1)*j+1]);
      dst_int_ptr = &(chunk->out_sizes[size_len*(Id-chunk->first_id)]);
      for (k=0; k<size_len; ++k) {
        *dst_int_ptr++ = *src_int_ptr++;
      }

      dst_ptr = &(chunk->out_values[value_len*(Id-chunk->first_id)]);
      src_ptr = &(recv_buf[value_len*j]);
      for (k=0; k<value_len; ++k) {
        *dst_ptr++ = *src_ptr++;
      }
    }


    CCAFREE(recv_size_buf);
    CCAFREE(send_size_buf);

    if (chunk->value_entry_length > 0) {
      CCAFREE(recv_buf);
      CCAFREE(send_buf);
    }
#endif
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the gathered data.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_write_chunk(BIN_OUT_FIELD *context,
                     BIN_OUT_CHUNK* chunk,
                     CHAR* entry_name)
{
#ifdef PARALLEL
  MPI_Status status;
  INT err;
#endif

  INT local_sizes[2];
#ifdef PARALLEL
  INT global_sizes[2];
#else
  INT* global_sizes = local_sizes;
#endif

  INT rank = context->actintra->intra_rank;

#ifdef DEBUG
  dstrc_enter("out_write_chunck");
#endif

  /*MPI_File_get_position(fh, &my_current_offset);*/

#ifdef PARALLEL
  err = MPI_File_seek(context->value_file,
                      context->value_file_offset +
                      chunk->first_id*chunk->value_entry_length*sizeof(DOUBLE),
                      MPI_SEEK_SET);
  if (err != 0) {
    dserror("MPI_File_seek failed: %d", err);
  }
  err = MPI_File_write(context->value_file, chunk->out_values,
                       chunk->value_count, MPI_DOUBLE, &status);
  if (err != 0) {
    dserror("MPI_File_write failed: %d", err);
  }

  if (chunk->size_entry_length > 0) {
    err = MPI_File_seek(context->size_file,
                        context->size_file_offset +
                        chunk->first_id*chunk->size_entry_length*sizeof(INT),
                        MPI_SEEK_SET);
    if (err != 0) {
      dserror("MPI_File_seek failed: %d", err);
    }
    err = MPI_File_write(context->size_file, chunk->out_sizes,
                         chunk->size_count, MPI_INT, &status);
    if (err != 0) {
      dserror("MPI_File_write failed: %d", err);
    }
  }
#else
  fseek(context->value_file, context->value_file_offset, SEEK_SET);
  if (fwrite(chunk->out_values,
             sizeof(DOUBLE),
             chunk->value_count,
             context->value_file) != chunk->value_count) {
    dserror("failed to write value file");
  }
  fseek(context->size_file, context->size_file_offset, SEEK_SET);
  if (fwrite(chunk->out_sizes,
             sizeof(INT),
             chunk->size_count,
             context->size_file) != chunk->size_count) {
    dserror("failed to write size file");
  }
#endif

  if (rank == 0) {
    CHAR* names[] = CHUNK_TYPE_NAMES;
    fprintf(bin_out_main.control_file, "    %s:\n", entry_name);
    map_insert_int_cpy(&chunk->group, context->value_file_offset, "value_offset");
    map_insert_int_cpy(&chunk->group, context->size_file_offset, "size_offset");
    map_insert_string_cpy(&chunk->group, names[chunk->type], "type");
    map_print(bin_out_main.control_file, &chunk->group, 8);
  }

  /* get new file offsets */
  /* in the sequential version local_sizes and global_sizes point to
   * the same memory. Thus evil code duplication is avoided. ;) */
  local_sizes[0] = chunk->value_count*sizeof(DOUBLE);
  local_sizes[1] = chunk->size_count*sizeof(INT);
#ifdef PARALLEL
  MPI_Allreduce(local_sizes, global_sizes, 2, MPI_INT, MPI_SUM, context->actintra->MPI_INTRA_COMM);
#endif
  context->value_file_offset += global_sizes[0];
  context->size_file_offset += global_sizes[1];

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief Write one node array.

  This is just for convenience. But it's an important part of the
  module's interface. This function is used a lot in restart code.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void out_node_arrays(BIN_OUT_FIELD* context,
                     NODE_ARRAY array)
{
  INT entry_size;

#ifdef DEBUG
  dstrc_enter("out_node_arrays");
#endif

  entry_size = context->max_size[array];
  out_node_chunk(context, nodearraynames[array], cc_node_array, entry_size, 2, array);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                    INT array)
{
  BIN_OUT_CHUNK array_data;

#ifdef DEBUG
  dstrc_enter("out_node_chunk");
#endif

  init_bin_out_chunk(context, &array_data, chunk_node, value_length, size_length);
  out_gather_values(context, &array_data, type, array);
  out_write_chunk(context, &array_data, chunk_name);

  destroy_bin_out_chunk(&array_data);

#ifdef DEBUG
  dstrc_exit();
#endif
}

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
                       INT array)
{
  BIN_OUT_CHUNK chunk;

#ifdef DEBUG
  dstrc_enter("out_element_chunk");
#endif

  init_bin_out_chunk(context, &chunk, chunk_element, value_length, size_length);
  out_gather_values(context, &chunk, type, array);
  out_write_chunk(context, &chunk, chunk_name);

  destroy_bin_out_chunk(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                       DIST_VECTOR* vectors)
{
  BIN_OUT_CHUNK chunk;

#ifdef DEBUG
  dstrc_enter("out_distvec_chunk");
#endif

  if (context->sysarray_typ == NULL) {
    dserror("cannot write distributed vectors without solver information");
  }

  chunk.vectors = vectors;

  init_bin_out_chunk(context, &chunk, chunk_dist_vec, length, 0);
  out_gather_values(context, &chunk, cc_dist_vector, 0);
  out_write_chunk(context, &chunk, chunk_name);

  destroy_bin_out_chunk(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*--------------------------------------------------------------------*/
/*
  Input

  There is no use out output anything if it's not going to be read
  back. And as ccarat attempts to write in its own format we'll have
  to read it ourselves to do something useful with that data. Yet
  there are different ways to write data. One way is to dump ccarat's
  internal memory to disk. This enables us to read it back and go on
  with our calculation provided the input file can be read again as we
  don't want to dump information we have read before. Another approach
  is to write some meaningful data. This is required for post
  processing but most of the time we're not going to read it back into
  ccarat. In preparing meaningful numbers some information is omitted
  that ccarat needs to go on in its operation. Thus there is no
  attempt to read such values here.
*/
/*--------------------------------------------------------------------*/


#ifdef PARALLEL


/*----------------------------------------------------------------------*/
/*!
  \brief Find the node ids of those nodes that need to be send and
  received.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void in_setup_node_transfer(BIN_IN_FIELD *context,
                                   INT rank,
                                   INT nprocs,
                                   INT *counters,
                                   PARTDISCRET *actpdis,
                                   DISCRET *actdis)
{
  INT i, j, max_num, fullarrays;
  INT tag_base = 0xABCDEF;
  INTRA* actintra = context->actintra;

#ifdef DEBUG
  dstrc_enter("in_setup_node_transfer");
#endif

  max_num = (actdis->numnp + nprocs - 1) / nprocs;
  fullarrays = nprocs - (nprocs*max_num - actdis->numnp);

  for (j=0; j<actpdis->numnp; ++j) {
    INT proc;
    NODE* actnode = actpdis->node[j];

    /* Find the processor this node has to be received from. */
    if (actnode->Id_loc < max_num*fullarrays) {
      proc = actnode->Id_loc / max_num;
    }
    else {
      INT Id = actnode->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }
    context->recv_numnp[proc]++;
  }

  /* allocate memory for the node numbers we have to receive */

  for (i=0; i<nprocs; ++i) {
    context->recv_node_ids[i] = (INT*)CCACALLOC(2*context->recv_numnp[i], sizeof(INT));
  }

  /*
   * Collect the ids of all nodes to be received. Actually we need two
   * ids per node. We need the per discretization id because this one
   * is used to figure out to which processor the node belongs. And we
   * need the per partition id that allows to access the receiving
   * node inside its partition. */
  memset(counters, 0, nprocs*sizeof(INT));
  for (j=0; j<actpdis->numnp; ++j) {
    INT proc;
    NODE* actnode = actpdis->node[j];

    /* Find the processor this node has to be received from. */
    if (actnode->Id_loc < max_num*fullarrays) {
      proc = actnode->Id_loc / max_num;
    }
    else {
      INT Id = actnode->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }

    context->recv_node_ids[proc][2*counters[proc]  ] = actnode->Id_loc;
    context->recv_node_ids[proc][2*counters[proc]+1] = j; /* actnode->Id_part; */
    counters[proc]++;
  }
  for (j=0; j<nprocs; ++j) {
    dsassert(counters[j] == context->recv_numnp[j], "node send count mismatch");
  }

  /* The information how many nodes we have to send can only be
   * found by communication. */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    /* Tell the source how many items we expect and receive from our
     * destination how many items it wants from us. */
    err = MPI_Sendrecv(&(context->recv_numnp[src]), 1, MPI_INT, src, tag_base-i-1,
                       &(context->send_numnp[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }

    context->send_node_ids[dst] = (INT*)CCACALLOC(2*context->send_numnp[dst], sizeof(INT));

    /* Tell the source which items we expect. */
    err = MPI_Sendrecv(context->recv_node_ids[src], 2*context->recv_numnp[src], MPI_INT, src, tag_base+i,
                       context->send_node_ids[dst], 2*context->send_numnp[dst], MPI_INT, dst, tag_base+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the element ids of those elements that need to be send
  and received.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void in_setup_element_transfer(BIN_IN_FIELD *context,
                                      INT rank,
                                      INT nprocs,
                                      INT *counters,
                                      PARTDISCRET *actpdis,
                                      DISCRET *actdis)
{
  INT i, j, max_num, fullarrays;
  INT tag_base = 0xBCDEF;
  INTRA* actintra = context->actintra;

#ifdef DEBUG
  dstrc_enter("in_setup_element_transfer");
#endif

  max_num = (actdis->numele + nprocs - 1) / nprocs;
  fullarrays = nprocs - (nprocs*max_num - actdis->numele);

  for (j=0; j<actpdis->numele; ++j) {
    INT proc;
    ELEMENT* actele = actpdis->element[j];

    /* Find the processor this element has to be received from. */
    if (actele->Id_loc < max_num*fullarrays) {
      proc = actele->Id_loc / max_num;
    }
    else {
      INT Id = actele->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }
    context->recv_numele[proc]++;
  }

  /* allocate memory for the element numbers we have to receive */

  for (i=0; i<nprocs; ++i) {
    context->recv_element_ids[i] = (INT*)CCACALLOC(2*context->recv_numele[i], sizeof(INT));
  }

  /* collect the elements */
  memset(counters, 0, nprocs*sizeof(INT));
  for (j=0; j<actpdis->numele; ++j) {
    INT proc;
    ELEMENT* actele = actpdis->element[j];

    /* Find the processor this element has to be received from. */
    if (actele->Id_loc < max_num*fullarrays) {
      proc = actele->Id_loc / max_num;
    }
    else {
      INT Id = actele->Id_loc - max_num*fullarrays;
      proc = Id / (max_num-1) + fullarrays;
    }

    context->recv_element_ids[proc][2*counters[proc]  ] = actele->Id_loc;
    context->recv_element_ids[proc][2*counters[proc]+1] = j;
    counters[proc]++;
  }
  for (j=0; j<nprocs; ++j) {
    dsassert(counters[j] == context->recv_numele[j], "element send count mismatch");
  }

  /* the element send counts have to be communicated */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    /* Tell the source how many items we expect and receive from our
     * destination how many items it wants from us. */
    err = MPI_Sendrecv(&(context->recv_numele[src]), 1, MPI_INT, src, tag_base-i-1,
                       &(context->send_numele[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }

    context->send_element_ids[dst] = (INT*)CCACALLOC(2*context->send_numele[dst], sizeof(INT));

    /* Tell the source which items we expect. */
    err = MPI_Sendrecv(context->recv_element_ids[src], 2*context->recv_numele[src], MPI_INT, src, tag_base+i,
                       context->send_element_ids[dst], 2*context->send_numele[dst], MPI_INT, dst, tag_base+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the dofs that need to be send and received.

  This is solver specific.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static void in_setup_dof_transfer(BIN_IN_FIELD *context,
                                  INT rank,
                                  INT nprocs,
                                  INT *counters,
                                  PARTDISCRET *actpdis,
                                  DISCRET *actdis)
{
  INT i, j, max_num, fullarrays;
  INT tag_base = 0xCDEF;
  SPARSE_TYP *sysarray_typ;
  SPARSE_ARRAY *sysarray;
  INTRA* actintra = context->actintra;

#ifdef DEBUG
  dstrc_enter("in_setup_dof_transfer");
#endif

  sysarray_typ = context->sysarray_typ;
  sysarray     = context->sysarray;

#define calc_max_number(numeq_total)                            \
  max_num = ((numeq_total) + nprocs - 1) / nprocs;              \
  fullarrays = nprocs - (nprocs*max_num - (numeq_total));

#define boilerplate_code(dof)                                   \
  {                                                             \
    INT proc;                                                   \
                                                                \
    /* Find the processor this dof has to be sent to. */        \
    if ((dof) < max_num*fullarrays) {                           \
      proc = (dof) / max_num;                                   \
    }                                                           \
    else {                                                      \
      INT Id = (dof) - max_num*fullarrays;                      \
      proc = Id / (max_num-1) + fullarrays;                     \
    }                                                           \
                                                                \
    context->recv_numdof[proc]++;                               \
  }

  /* Find out how many dofs we have to send to which processor. */
  switch (*sysarray_typ) {

#ifdef AZTEC_PACKAGE
  case msr:
    calc_max_number(sysarray->msr->numeq_total);
    for (i=0; i<sysarray->msr->numeq; ++i) {
      INT dof = sysarray->msr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef HYPRE_PACKAGE
  case parcsr: {
    INT rank = context->actintra->intra_rank;
    calc_max_number(sysarray->parcsr->numeq_total);
    for (i=0; i<sysarray->parcsr->numeq; ++i) {
      INT dof = sysarray->parcsr->update.a.ia[rank][i];
      boilerplate_code(dof);
    }
    break;
  }
#endif

#ifdef PARSUPERLU_PACKAGE
  case ucchb:
    calc_max_number(sysarray->ucchb->numeq_total);
    for (i=0; i<sysarray->ucchb->numeq; ++i) {
      INT dof = sysarray->ucchb->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case dense:
    calc_max_number(sysarray->dense->numeq_total);
    for (i=0; i<sysarray->dense->numeq; ++i) {
      INT dof = sysarray->dense->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

#ifdef MLIB_PACKAGE
  case mds:
    calc_max_number(sysarray->mds->numeq_total);
    for (i=0; i<sysarray->mds->numeq; ++i) {
      INT dof = i;
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef MUMPS_PACKAGE
  case rc_ptr:
    calc_max_number(sysarray->rc_ptr->numeq_total);
    for (i=0; i<sysarray->rc_ptr->numeq; ++i) {
      INT dof = sysarray->rc_ptr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef SPOOLES_PACKAGE
  case spoolmatrix:
    calc_max_number(sysarray->spo->numeq_total);
    for (i=0; i<sysarray->spo->numeq; ++i) {
      INT dof = sysarray->spo->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef UMFPACK
  case ccf:
    calc_max_number(sysarray->ccf->numeq_total);
    for (i=0; i<sysarray->ccf->numeq; ++i) {
      INT dof = sysarray->ccf->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case skymatrix:
    calc_max_number(sysarray->sky->numeq_total);
    for (i=0; i<sysarray->sky->numeq; ++i) {
      INT dof = sysarray->sky->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

#ifdef MLPCG
  case bdcsr:
    calc_max_number(sysarray->bdcsr->numeq_total);
    for (i=0; i<sysarray->bdcsr->numeq; ++i) {
      INT dof = sysarray->bdcsr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case oll:
    calc_max_number(sysarray->oll->numeq_total);
    for (i=0; i<sysarray->oll->numeq; ++i) {
      INT dof = sysarray->oll->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

  default:
    dserror("Unknown type %d of system matrix", *sysarray_typ);
    break;
  }

#undef boilerplate_code
#undef calc_max_number

  /* allocate memory for the dofs we have to receive */

  for (i=0; i<nprocs; ++i) {
    context->recv_dof_ids[i] = (INT*)CCACALLOC(2*context->recv_numdof[i], sizeof(INT));
  }

  /* collect the dofs */
  memset(counters, 0, nprocs*sizeof(INT));

#define boilerplate_code(dof)                                   \
  {                                                             \
    INT proc;                                                   \
                                                                \
    /* Find the processor this dof has to be received from. */  \
    if ((dof) < max_num*fullarrays) {                           \
      proc = (dof) / max_num;                                   \
    }                                                           \
    else {                                                      \
      INT Id = (dof) - max_num*fullarrays;                      \
      proc = Id / (max_num-1) + fullarrays;                     \
    }                                                           \
                                                                \
    context->recv_dof_ids[proc][2*counters[proc]  ] = (dof);    \
    context->recv_dof_ids[proc][2*counters[proc]+1] = i;        \
    counters[proc]++;                                           \
  }

  switch (*sysarray_typ) {

#ifdef AZTEC_PACKAGE
  case msr:
    for (i=0; i<sysarray->msr->numeq; ++i) {
      INT dof = sysarray->msr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef HYPRE_PACKAGE
  case parcsr: {
    INT rank = context->actintra->intra_rank;
    for (i=0; i<sysarray->parcsr->numeq; ++i) {
      INT dof = sysarray->parcsr->update.a.ia[rank][i];
      boilerplate_code(dof);
    }
    break;
  }
#endif

#ifdef PARSUPERLU_PACKAGE
  case ucchb:
    for (i=0; i<sysarray->ucchb->numeq; ++i) {
      INT dof = sysarray->ucchb->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case dense:
    for (i=0; i<sysarray->dense->numeq; ++i) {
      INT dof = sysarray->dense->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

#ifdef MLIB_PACKAGE
  case mds:
    for (i=0; i<sysarray->mds->numeq; ++i) {
      INT dof = i;
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef MUMPS_PACKAGE
  case rc_ptr:
    for (i=0; i<sysarray->rc_ptr->numeq; ++i) {
      INT dof = sysarray->rc_ptr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef SPOOLES_PACKAGE
  case spoolmatrix:
    for (i=0; i<sysarray->spo->numeq; ++i) {
      INT dof = sysarray->spo->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

#ifdef UMFPACK
  case ccf:
    for (i=0; i<sysarray->ccf->numeq; ++i) {
      INT dof = sysarray->ccf->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case skymatrix:
    for (i=0; i<sysarray->sky->numeq; ++i) {
      INT dof = sysarray->sky->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

#ifdef MLPCG
  case bdcsr:
    for (i=0; i<sysarray->bdcsr->numeq; ++i) {
      INT dof = sysarray->bdcsr->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;
#endif

  case oll:
    for (i=0; i<sysarray->oll->numeq; ++i) {
      INT dof = sysarray->oll->update.a.iv[i];
      boilerplate_code(dof);
    }
    break;

  default:
    dserror("Unknown type %d of system matrix", *sysarray_typ);
    break;
  }

#undef boilerplate_code

  for (j=0; j<nprocs; ++j) {
    dsassert(counters[j] == context->recv_numdof[j], "dof send count mismatch");
  }

  /* the dof send counts have to be communicated */
  for (i=0; i<nprocs; ++i) {
    MPI_Status status;
    INT dst;
    INT src;
    INT err;

    dst = (rank + i + 1) % nprocs;
    src = (nprocs + rank - i - 1) % nprocs;

    /* Tell the source how many items we expect and receive from our
     * destination how many items it wants from us. */
    err = MPI_Sendrecv(&(context->recv_numdof[src]), 1, MPI_INT, src, tag_base-i-1,
                       &(context->send_numdof[dst]), 1, MPI_INT, dst, tag_base-i-1,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }

    context->send_dof_ids[dst] = (INT*)CCACALLOC(2*context->send_numdof[dst], sizeof(INT));

    /* Tell the source which items we expect. */
    err = MPI_Sendrecv(context->recv_dof_ids[src], 2*context->recv_numdof[src], MPI_INT, src, tag_base+i,
                       context->send_dof_ids[dst], 2*context->send_numdof[dst], MPI_INT, dst, tag_base+i,
                       actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


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
                       INT disnum)
{
#ifdef PARALLEL
  INT nprocs = actintra->intra_nprocs;
  INT rank = actintra->intra_rank;
  INT err;
  PARTDISCRET* actpdis = &(actpart->pdis[disnum]);
  INT* counters;
#endif
  DISCRET* actdis = &(actfield->dis[disnum]);
  SYMBOL* dir;
  CHAR* filename;

#ifdef DEBUG
  dstrc_enter("init_bin_in_field");
#endif

  /* Remember the variables of this field. */
  context->sysarray_typ = sysarray_typ;
  context->sysarray = sysarray;

  context->actfield = actfield;
  context->actpart = actpart;
  context->actintra = actintra;
  context->disnum = disnum;

  /*--------------------------------------------------------------------*/
  /* find the describtion of this field in the control file */
  dir = map_find_symbol(&(bin_in_main.table), "field");
  while (dir != NULL) {
    if (symbol_is_map(dir)) {
      MAP* map;
      symbol_get_map(dir, &map);
      if (map_has_string(map, "type", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "discretization", disnum)) {
        /* The group entry that contains the distinguishing
         * definitions. This group contains information we need to know. */
        context->field_info = map;
        break;
      }
    }
    dir = dir->next;
  }
  if (dir == NULL) {
    dserror("No field entry in symbol table. Control file corrupt?");
  }

  /*--------------------------------------------------------------------*/
  /* open file to read */

  filename = map_read_string(context->field_info, "value_file");
#ifdef PARALLEL
  err = MPI_File_open(actintra->MPI_INTRA_COMM, filename,
                      MPI_MODE_RDONLY, MPI_INFO_NULL,
                      &(context->value_file));
  if (err != 0) {
    dserror("MPI_File_open failed for file '%s': %d", filename, err);
  }
#else
  context->value_file = fopen(filename, "r");
#endif

  filename = map_read_string(context->field_info, "size_file");
#ifdef PARALLEL
  err = MPI_File_open(actintra->MPI_INTRA_COMM, filename,
                      MPI_MODE_RDONLY, MPI_INFO_NULL,
                      &(context->size_file));
  if (err != 0) {
    dserror("MPI_File_open failed for file '%s': %d", filename, err);
  }
#else
  context->size_file = fopen(filename, "r");
#endif

  {
    INT numnp;
    INT numele;

    numnp = map_read_int(context->field_info, "numnp");
    numele = map_read_int(context->field_info, "numele");

    dsassert(numnp == actdis->numnp, "wrong node number");
    dsassert(numele == actdis->numele, "wrong element number");
  }

#ifdef PARALLEL

  /*--------------------------------------------------------------------*/
  /* setup communication variables */

  context->send_numnp = (INT*)CCACALLOC(nprocs, sizeof(INT));
  context->recv_numnp = (INT*)CCACALLOC(nprocs, sizeof(INT));

  context->send_numele = (INT*)CCACALLOC(nprocs, sizeof(INT));
  context->recv_numele = (INT*)CCACALLOC(nprocs, sizeof(INT));

  context->send_node_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));
  context->recv_node_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));

  context->send_element_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));
  context->recv_element_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));

  if (sysarray_typ != NULL) {
    context->send_numdof = (INT*)CCACALLOC(nprocs, sizeof(INT));
    context->recv_numdof = (INT*)CCACALLOC(nprocs, sizeof(INT));

    context->send_dof_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));
    context->recv_dof_ids = (INT**)CCACALLOC(nprocs, sizeof(INT*));
  }

  /*--------------------------------------------------------------------*/
  /* Find out how many nodes we have to receive here. */

  counters = (INT*)CCACALLOC(nprocs, sizeof(INT));

  in_setup_node_transfer(context, rank, nprocs, counters, actpdis, actdis);

  /*--------------------------------------------------------------------*/
  /* Find out how many elements we have to receive here. */

  in_setup_element_transfer(context, rank, nprocs, counters, actpdis, actdis);

  /*--------------------------------------------------------------------*/
  /* Find out how many dofs we have to receive here. */

  /* If there is no solver we won't be able to load distributed
   * vectors. */
  if (sysarray_typ != NULL) {
    in_setup_dof_transfer(context, rank, nprocs, counters, actpdis, actdis);
  }

  CCAFREE(counters);

  /*--------------------------------------------------------------------*/

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Init the data structure to read and distribute on
  particular node array.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_bin_in_chunk(BIN_IN_FIELD* context,
                       BIN_IN_CHUNK* chunk,
                       MAP* result_info,
                       CHAR* group_name,
                       CHUNK_TYPE type,
                       NODE_ARRAY array)
{
  INT nprocs = context->actintra->intra_nprocs;
  INT rank = context->actintra->intra_rank;
  DISCRET* actdis = &(context->actfield->dis[context->disnum]);
  INT num;

#ifdef DEBUG
  dstrc_enter("init_bin_in_chunk");
#endif

  chunk->field = context;

  switch (type) {
  case chunk_node:
    num = actdis->numnp;
    chunk->type = chunk_node;
    break;
  case chunk_element:
    num = actdis->numele;
    chunk->type = chunk_element;
    break;
  case chunk_dist_vec:
    /* We are able to handle an array of distributed vectors. But
     * there must be a first one and all of them must have the same
     * size. */
    num = chunk->vectors[0].numeq_total;
    chunk->type = chunk_dist_vec;
    break;
  default:
    dserror("unknown chunk type %d", type);
  }

  /*--------------------------------------------------------------------*/
  /* find the describtion of this array in the control file */

  chunk->group_info = map_read_map(result_info, group_name);

  chunk->value_entry_length = map_read_int(chunk->group_info, "value_entry_length");
  chunk->size_entry_length  = map_read_int(chunk->group_info, "size_entry_length");

  if ((chunk->value_entry_length < 0) || (chunk->size_entry_length < 0)) {
    dserror("illegal item sites: %d, %d",
            chunk->value_entry_length, chunk->size_entry_length);
  }

  /*--------------------------------------------------------------------*/

  /* Figure out how many nodes are read by each processor. */
  chunk->num = (num + nprocs - 1) / nprocs;

  /* The number of processors that read the full count. The other
   * ones have one less to read. */
  chunk->fullarrays = nprocs - (nprocs*chunk->num - num);

  /* The Id_loc of the first node of this processor. */
  chunk->first_id = chunk->num*rank;

  if (rank >= chunk->fullarrays) {
    chunk->num -= 1;
    chunk->first_id -= rank - chunk->fullarrays;
  }

  chunk->value_count = chunk->value_entry_length*chunk->num;
  chunk->size_count = chunk->size_entry_length*chunk->num;
  if (chunk->value_count > 0) {
    chunk->in_values = (DOUBLE*)CCACALLOC(chunk->value_count, sizeof(DOUBLE));
  }
  else {
    chunk->in_values = NULL;
  }

  if (chunk->size_count > 0) {
    chunk->in_sizes = (INT*)CCACALLOC(chunk->size_count, sizeof(INT));
  }
  else {
    chunk->in_sizes = NULL;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Distribute the values read to the nodes.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void in_scatter_chunk(BIN_IN_FIELD* context,
                      BIN_IN_CHUNK* chunk,
                      CHUNK_CONTENT_TYPE type,
                      NODE_ARRAY array)
{
#ifdef PARALLEL
  INT tag_base = 0xABCD;
  INT rank = context->actintra->intra_rank;
#endif
  INT nprocs = context->actintra->intra_nprocs;
  INT proc;
  INT num;

#ifdef DEBUG
  dstrc_enter("in_node_arrays");
#endif

  switch (chunk->type) {
  case chunk_node:
    num = context->actfield->dis[context->disnum].numnp;
    break;
  case chunk_element:
    num = context->actfield->dis[context->disnum].numele;
    break;
  case chunk_dist_vec:
    /* We are able to handle an array of distributed vectors. But
     * there must be a first one and all of them must have the same
     * size. */
    num = chunk->vectors[0].numeq_total;
    break;
  default:
    dserror("unknown chunk type %d", chunk->type);
  }

  /* distribute */
  for (proc=0; proc<nprocs; ++proc) {
    INT src = 0;
#ifdef PARALLEL
    INT dst;
    INT j;
    MPI_Status status;
    DOUBLE* send_buf;
    DOUBLE* recv_buf;
    INT send_count = 0;
    INT recv_count = 0;
    INT* send_size_buf;
    INT* recv_size_buf;
    INT send_size_count = 0;
    INT recv_size_count = 0;
    INT err;
    INT size_len;
    INT value_len;
    INT send_num;
    INT recv_num;
    INT* send_ids;
    /*INT* recv_ids; */

    dst = (rank + proc + 1) % nprocs;
    src = (nprocs + rank - proc - 1) % nprocs;

    switch (chunk->type) {
    case chunk_node:
      send_num = context->send_numnp[dst];
      recv_num = context->recv_numnp[src];
      send_ids = context->send_node_ids[dst];
      /*recv_ids = context->recv_node_ids[src];*/
      break;
    case chunk_element:
      send_num = context->send_numele[dst];
      recv_num = context->recv_numele[src];
      send_ids = context->send_element_ids[dst];
      /*recv_ids = context->recv_element_ids[src];*/
      break;
    case chunk_dist_vec:
      send_num = context->send_numdof[dst];
      recv_num = context->recv_numdof[src];
      send_ids = context->send_dof_ids[dst];
      break;
    default:
      dserror("unknown chunk type %d", chunk->type);
    }

    if (chunk->value_entry_length > 0) {
      send_count = send_num*chunk->value_entry_length;
      recv_count = recv_num*chunk->value_entry_length;

      send_buf = (DOUBLE*)CCACALLOC(send_count, sizeof(DOUBLE));
      recv_buf = (DOUBLE*)CCACALLOC(recv_count, sizeof(DOUBLE));
    }

    send_size_count = send_num*(chunk->size_entry_length+1);
    recv_size_count = recv_num*(chunk->size_entry_length+1);

    send_size_buf = (INT*)CCACALLOC(send_size_count, sizeof(INT));
    recv_size_buf = (INT*)CCACALLOC(recv_size_count, sizeof(INT));

    /* abbreviations */
    size_len = chunk->size_entry_length;
    value_len = chunk->value_entry_length;

    /* pack data to be sent */
    for (j=0; j<send_num; ++j) {
      INT k;
      INT Id;
      DOUBLE *src_ptr;
      DOUBLE *dst_ptr;
      INT *src_int_ptr;
      INT *dst_int_ptr;

      Id = send_ids[2*j];

      /* make sure we read this node */
      dsassert((Id>=chunk->first_id) &&
               (Id<chunk->first_id+chunk->num), "inadequate id");

      send_size_buf[(size_len+1)*j] = Id;

      /* we don't store the id which is always transfered first */
      src_int_ptr = &(chunk->in_sizes[size_len*(Id-chunk->first_id)]);
      dst_int_ptr = &(send_size_buf[(size_len+1)*j+1]);
      for (k=0; k<size_len; ++k) {
        *dst_int_ptr++ = *src_int_ptr++;
      }

      dst_ptr = &(send_buf[value_len*j]);
      src_ptr = &(chunk->in_values[value_len*(Id-chunk->first_id)]);
      for (k=0; k<value_len; ++k) {
        *dst_ptr++ = *src_ptr++;
      }
    }

    /* communication */
    err = MPI_Sendrecv(send_size_buf, send_size_count, MPI_INT, dst, tag_base-proc-1,
                       recv_size_buf, recv_size_count, MPI_INT, src, tag_base-proc-1,
                       context->actintra->MPI_INTRA_COMM, &status);
    if (err != 0) {
      dserror("mpi sendrecv error %d", err);
    }

    if (chunk->value_entry_length > 0) {
      err = MPI_Sendrecv(send_buf, send_count, MPI_DOUBLE, dst, tag_base+proc,
                         recv_buf, recv_count, MPI_DOUBLE, src, tag_base+proc,
                         context->actintra->MPI_INTRA_COMM, &status);
      if (err != 0) {
        dserror("mpi sendrecv error %d", err);
      }
    }
#else
    DOUBLE* recv_buf = chunk->in_values;
    INT* recv_size_buf = chunk->in_sizes;

    INT recv_count = chunk->value_count;
    INT recv_size_count = chunk->size_count;

#endif

    /* unpack */
    in_unpack_items(context, chunk, type, array,
                    recv_buf, recv_count, recv_size_buf, recv_size_count, src);

#ifdef PARALLEL

    if (chunk->value_entry_length > 0) {
      CCAFREE(recv_buf);
      CCAFREE(send_buf);
    }

    CCAFREE(recv_size_buf);
    CCAFREE(send_size_buf);

#endif

  }

#ifdef PARALLEL

#if 0
  /* No yet. This is a parallelization function. */

  /* Copy the new values to the external nodes, too. */
  solserv_dist_node_arrays_any(context->actfield, context->actpart,
                               context->actintra, context->disnum, array);
#endif

#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void in_read_chunk(BIN_IN_FIELD *context,
                   BIN_IN_CHUNK *chunk)
{
  INT offset;
#ifdef PARALLEL
  INT err;
  MPI_Status status;
#endif

#ifdef DEBUG
  dstrc_enter("in_read_chunk");
#endif

  /* read the values */
  offset = map_read_int(chunk->group_info, "value_offset");
#ifdef PARALLEL
  err = MPI_File_seek(context->value_file,
                      offset+chunk->first_id*sizeof(DOUBLE)*chunk->value_entry_length,
                      MPI_SEEK_SET);
  if (err != 0) {
    dserror("MPI_File_seek failed: %d", err);
  }
  err = MPI_File_read(context->value_file, chunk->in_values,
                      chunk->value_count, MPI_DOUBLE, &status);
  if (err != 0) {
    dserror("MPI_File_read failed: %d", err);
  }
#else
  fseek(context->value_file, offset, SEEK_SET);
  if (fread(chunk->in_values,
            sizeof(DOUBLE),
            chunk->value_count,
            context->value_file) != chunk->value_count) {
    dserror("failed to read value file");
  }
#endif

  /* read the array dimensions */
  offset = map_read_int(chunk->group_info, "size_offset");
#ifdef PARALLEL
  err = MPI_File_seek(context->size_file,
                      offset+chunk->first_id*sizeof(INT)*chunk->size_entry_length,
                      MPI_SEEK_SET);
  if (err != 0) {
    dserror("MPI_File_seek failed: %d", err);
  }
  err = MPI_File_read(context->size_file, chunk->in_sizes,
                      chunk->size_count, MPI_INT, &status);
  if (err != 0) {
    dserror("MPI_File_read failed: %d", err);
  }
#else
  fseek(context->size_file, offset, SEEK_SET);
  if (fread(chunk->in_sizes,
            sizeof(INT),
            chunk->size_count,
            context->size_file) != chunk->size_count) {
    dserror("failed to read size file");
  }
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                    NODE_ARRAY array)
{
#ifdef DEBUG
  dstrc_enter("in_node_arrays");
#endif

  in_node_chunk(context, result_info, nodearraynames[array], array);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                   NODE_ARRAY array)
{
  BIN_IN_CHUNK chunk;

#ifdef DEBUG
  dstrc_enter("in_node_chunk");
#endif

  init_bin_in_chunk(context, &chunk, result_info, group_name, chunk_node, array);
  in_read_chunk(context, &chunk);
  in_scatter_chunk(context, &chunk, cc_node_array, array);

  destroy_bin_in_chunk(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                      CHUNK_CONTENT_TYPE type)
{
  BIN_IN_CHUNK chunk;

#ifdef DEBUG
  dstrc_enter("in_element_chunk");
#endif

  init_bin_in_chunk(context, &chunk, result_info, group_name, chunk_element, 0);
  in_read_chunk(context, &chunk);
  in_scatter_chunk(context, &chunk, type, 0);

  destroy_bin_in_chunk(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                      DIST_VECTOR* vectors)
{
  BIN_IN_CHUNK chunk;

#ifdef DEBUG
  dstrc_enter("in_distvec_chunk");
#endif

  if (context->sysarray_typ == NULL) {
    dserror("no solver set, cannot load distributed vector");
  }

  chunk.vectors = vectors;

  init_bin_in_chunk(context, &chunk, result_info, group_name, chunk_dist_vec, 0);
  in_read_chunk(context, &chunk);
  in_scatter_chunk(context, &chunk, cc_dist_vector, 0);

  destroy_bin_in_chunk(&chunk);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the description of the restart data in the control file

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
MAP *in_find_restart_group(FIELD *actfield, INT disnum, INT step)
{
  MAP *result_info = NULL;
  SYMBOL *symbol;
  INT field_pos;

#ifdef DEBUG
  dstrc_enter("in_find_restart_group");
#endif


  for (field_pos=0; field_pos<genprob.numfld; ++field_pos) {
    if (actfield == &(field[field_pos])) {
      break;
    }
  }

  if (field_pos==genprob.numfld) {
    dserror("unregistered field object");
  }

  symbol = map_find_symbol(&(bin_in_main.table), "restart");
  while (symbol != NULL) {
    if (symbol_is_map(symbol)) {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", fieldnames[actfield->fieldtyp]) &&
          map_has_int(map, "field_pos", field_pos) &&
          map_has_int(map, "discretization", disnum) &&
          map_has_int(map, "step", step)) {
        result_info = map;
        break;
      }
    }
    symbol = symbol->next;
  }
  if (symbol == NULL) {
    dserror("No restart entry for step %d in symbol table. Control file corrupt?", step);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return result_info;
}

#endif
