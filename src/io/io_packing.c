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

\author u.kue
\date 08/04

*/

#ifdef BINIO

/*!
  \addtogroup IO
*//*! @{ (documentation module open)*/

#include "io_packing.h"
#include "io_singlefile.h"
#include "io_elements.h"

#include "../pss_full/pss_table.h"
#include "../pss_full/pss_set.h"

#include "../shell8/shell8.h"
#include "../shell9/shell9.h"
#include "../wall1/wall1.h"
#include "../brick1/brick1.h"
#include "../beam3/beam3.h"
#include "../fluid2/fluid2.h"
#include "../fluid2_pro/fluid2pro.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#include "../fluid3_pro/fluid3pro.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_pro/fluid3pro_prototypes.h"
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

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
  *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN       *design;

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
  | structure of flags to control output                                 |
  | defined in out_global.c                                              |
  *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of material laws                                              |
  | defined in global_control.c
  *----------------------------------------------------------------------*/
extern struct _MATERIAL     *mat;

/*!----------------------------------------------------------------------
  \brief vector of multilayer material law

  <pre>                                                            sh 10/02
  This structure struct _MULTIMAT  *multimat is defined in global_control.c
  and the type is in materials.h
  It holds all information about the layered section data
  </pre>
  *----------------------------------------------------------------------*/
extern struct _MULTIMAT  *multimat;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | pointer to allocate dynamic variables if needed                      |
  | dedfined in global_control.c                                         |
  | struct _ALLDYNA       *alldyn;                                       |
  *----------------------------------------------------------------------*/
extern ALLDYNA *alldyn;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for output.

  This structure needs to be initialized at startup. The whole output
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_OUT_MAIN bin_out_main;

/*----------------------------------------------------------------------*/
/*!
  \brief The static variables used for input.

  This structure needs to be initialized at startup. The whole input
  mechanism is based on it.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
extern BIN_IN_MAIN bin_in_main;


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
  \brief Mark all element types that are used in this discretization.

  Here the context's element_flag table is set up. That is all
  elements are checked and the corresponding entry in the flag table
  gets a tic.

  The result is allreduced.

  \param context  (i/o) pointer to an output context in the setup phase
  \param actintra (i)   the communicator
  \param actpdis  (i)   the partitions discretization

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void out_find_element_types(struct _BIN_OUT_FIELD *context,
                            INTRA *actintra,
                            PARTDISCRET *actpdis)
{
  INT i;
  INT element_flag[el_count];

#ifdef DEBUG
  dstrc_enter("out_find_element_types");
#endif

  memset(element_flag, 0, el_count*sizeof(INT));

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    element_flag[actele->eltyp] = 1;
  }

  /* memcpy takes the buffer length in bytes. MPI asks for the number
   * of entries. */
#ifdef PARALLEL
  MPI_Allreduce(element_flag, context->element_flag, el_count, MPI_INT,
                MPI_MAX, actintra->MPI_INTRA_COMM);
#else
  memcpy(context->element_flag, element_flag, el_count*sizeof(INT));
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write element specific control information to the given file.

  This is only called if rank==0.

  The purpose is to write all kinds of information that are needed to
  interpret the elements' results. Everything goes to the group that
  describes the discretization the elements belong to.

  \author u.kue
  \date 11/04
*/
/*----------------------------------------------------------------------*/
void out_element_control(struct _BIN_OUT_FIELD *context,
                         FILE* file)
{
  INT i;

#ifdef DEBUG
  dstrc_enter("out_element_control");
#endif

  for (i=0; i<el_count; ++i)
  {
    if (context->element_flag[i])
    {
      switch (i)
      {
#ifdef D_SHELL8
      case el_shell8:
        break;
#endif
#ifdef D_SHELL9
      case el_shell9:
        break;
#endif
#ifdef D_BRICK1
      case el_brick1:
      {
        INT j;
        CHAR* c1_stresstype[] = BRICK1_STRESSTYPE;

        fprintf(file, "    c1_stresstypes:\n");
        for (j=0; c1_stresstype[j] != NULL; ++j)
        {
          fprintf(file, "        %s = %d\n", c1_stresstype[j], j);
        }
        fprintf(file, "\n");

        break;
      }
#endif
#ifdef D_FLUID2
      case el_fluid2:
        break;
#endif
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
        break;
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
        break;
#endif
#ifdef D_FLUID2TU
      case el_fluid2_tu:
        break;
#endif
#ifdef D_FLUID3
      case el_fluid3:
        break;
#endif
#ifdef D_FLUID3_F
      case el_fluid3_fast:
        break;
#endif
#ifdef D_ALE
      case el_ale2:
        break;
#endif
#ifdef D_ALE
      case el_ale3:
        break;
#endif
#ifdef D_WALL1
      case el_wall1:
      {
        INT j;
        CHAR* w1_stresstype[] = WALL1_STRESSTYPE;

        fprintf(file, "    w1_stresstypes:\n");
        for (j=0; w1_stresstype[j] != NULL; ++j)
        {
          fprintf(file, "        %s = %d\n", w1_stresstype[j], j);
        }
        fprintf(file, "\n");

        break;
      }
#endif
#ifdef D_BEAM3
      case el_beam3:
        break;
#endif
#ifdef D_AXISHELL
      case el_axishell:
        break;
#endif
#ifdef D_INTERF
      case el_interf:
      {
        INT j;
        CHAR* interf_stresstype[] = INTERF_STRESSTYPE;

        fprintf(file, "    interf_stresstypes:\n");
        for (j=0; interf_stresstype[j] != NULL; ++j)
        {
          fprintf(file, "        %s = %d\n", interf_stresstype[j], j);
        }
        fprintf(file, "\n");

        break;
      }
#endif
#ifdef D_WALLGE
      case el_wallge:
      {
        INT j;
        CHAR* wallge_stresstype[] = WALLGE_STRESSTYPE;

        fprintf(file, "    wallge_stresstypes:\n");
        for (j=0; wallge_stresstype[j] != NULL; ++j)
        {
          fprintf(file, "        %s = %d\n", wallge_stresstype[j], j);
        }
        fprintf(file, "\n");

        break;
      }
#endif
      default:
        dserror("element type %d unsupported", i);
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*======================================================================*/
/* There are special elements that demand a special setup for the
 * discretization output object if they are used. The following
 * functions provide just that. They are called by init_bin_out_field,
 * the discretization output constructor. */
/*======================================================================*/


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
void out_shell8_setup(struct _BIN_OUT_FIELD *context)
{
  INT rank = context->actintra->intra_rank;
  ELEMENT* actele;

#ifdef DEBUG
  dstrc_enter("out_shell8_setup");
#endif

  actele = context->actpart->pdis[context->disnum].element[0];
  if (actele->eltyp == el_shell8)
  {
    SHELL8* s8;
    INT numnp;

    /*
     * OK. We have a shell8 problem.  */
    context->is_shell8_problem = 1;
    s8 = actele->e.s8;

    /*
     * Mark this discretization. The filter will know that it is a
     * special one. */
    if (rank == 0)
    {
      INT j;
      CHAR* forcetype[] = S8_FORCETYPE;
      fprintf(bin_out_main.control_file,
              "    shell8_problem = \"yes\"\n\n"
              "    s8_forcetypes:\n");
      for (j=0; forcetype[j] != NULL; ++j)
      {
        fprintf(bin_out_main.control_file,
                "        %s = %d\n", forcetype[j], j);
      }
      fprintf(bin_out_main.control_file, "\n");
    }

    numnp = actele->numnp;

    /* There is a director at each node. Furthermore we need the
     * element thickness. Let's play it save and output one thickness
     * value per node. This is more that is currently supported by
     * ccarat's input though. */
    out_element_chunk(context, "shell8_director", cc_shell8_director, (3+1)*numnp, 0, 0);
  }
  else
  {
    context->is_shell8_problem = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

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
void out_shell9_setup(struct _BIN_OUT_FIELD *context)
{
  INT rank = context->actintra->intra_rank;
  ELEMENT* actele;

#ifdef DEBUG
  dstrc_enter("out_shell9_setup");
#endif

  actele = context->actpart->pdis[context->disnum].element[0];
  if (actele->eltyp == el_shell9)
  {
    SHELL9* s9;
    INT klay;
    INT layer_count;

    /*
     * OK. We have a shell9 problem.  */
    context->is_shell9_problem = 1;

    s9 = actele->e.s9;

    context->s9_layers = 0;
    for (klay=0; klay<s9->num_klay; klay++)
    {
      context->s9_layers += s9->kinlay[klay].num_mlay;
    }

    /*
     * Mark this discretization. The filter will know that it is a
     * special one. */
    if (rank == 0)
    {
      INT j;
      CHAR* forcetype[] = S9_FORCETYPE;
      fprintf(bin_out_main.control_file,
              "    shell9_problem = \"yes\"\n"
              "    shell9_smoothed = \"%s\"\n"
              "    shell9_layers = %d\n"
              "\n"
              "    s9_forcetypes:\n",
              (ioflags.struct_stress_smo ? "yes" : "no"),
              context->s9_layers);
      for (j=0; forcetype[j] != NULL; ++j)
      {
        fprintf(bin_out_main.control_file,
                "        %s = %d\n", forcetype[j], j);
      }
      fprintf(bin_out_main.control_file, "\n");
    }

    /* We need space for three values per artificial node. */

    context->s9_numnp = actele->numnp;

    if (context->s9_numnp == 4)
    {
      layer_count = 1;
    }
    else
    {
      layer_count = 2;
    }

    out_element_chunk(context, "shell9_coords", cc_shell9_coords,
                      3*context->s9_numnp*(layer_count*context->s9_layers+1), 0, 0);
  }
  else
  {
    context->is_shell9_problem = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*======================================================================*/
/* Here we have lots of functions that do nothing but packing
 * something (node arrays, element stuff, distributed vectors) into
 * the output buffers. These functions share a common structure. See
 * ``out_pack_items`` for an explanation of their parameters.
 *
 * These are the functions that need to be changed according to
 * changes in the elements. */
/*======================================================================*/



/*----------------------------------------------------------------------*/
/*!
  \brief Pack the values of the indicated node array in some buffer to
  send them to their writing processor.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_node_arrays(BIN_OUT_CHUNK *chunk,
                                 NODE_ARRAY array,
                                 PARTDISCRET *actpdis,
                                 DOUBLE *send_buf,
                                 INT send_count,
                                 INT *send_size_buf,
                                 INT dst_first_id,
                                 INT dst_num)
{
  INT j, counter;

#ifdef DEBUG
  dstrc_enter("out_pack_node_arrays");
#endif

#ifdef PARALLEL
#define pack_send_size_buf(node_array)                          \
  send_size_buf[3*counter  ] = actnode->Id_loc;                 \
  send_size_buf[3*counter+1] = actnode->node_array.fdim;        \
  send_size_buf[3*counter+2] = actnode->node_array.sdim;
#else
#define pack_send_size_buf(node_array)                          \
  send_size_buf[2*counter  ] = actnode->node_array.fdim;        \
  send_size_buf[2*counter+1] = actnode->node_array.sdim;
#endif

  /* gather the nodes to be send to processor i */
  counter = 0;
#define call_gather_nodes(node_array)                                   \
  dsassert(chunk->size_entry_length == 2, "invalid size entry length"); \
  for (j=0; j<actpdis->numnp; ++j)                                      \
  {                                                                     \
    NODE* actnode = actpdis->node[j];                                   \
    if ((actnode->Id_loc >= dst_first_id) &&                            \
        (actnode->Id_loc < dst_first_id+dst_num))                       \
    {                                                                   \
      DOUBLE *src_ptr;                                                  \
      DOUBLE *dst_ptr;                                                  \
      INT k;                                                            \
      INT size = actnode->node_array.fdim*actnode->node_array.sdim;     \
      dsassert(size <= chunk->field->max_size[node_array_ ## node_array], \
               "outdated size calculation. panic.");                    \
      src_ptr = actnode->node_array.a.da[0];                            \
      dst_ptr = &(send_buf[chunk->value_entry_length*counter]);         \
      pack_send_size_buf(node_array);                                   \
      for (k=0; k<size; ++k)                                            \
      {                                                                 \
        *dst_ptr++ = *src_ptr++;                                        \
      }                                                                 \
      counter += 1;                                                     \
    }                                                                   \
  }                                                                     \
  dsassert(counter*chunk->value_entry_length == send_count,             \
           "node pack count mismatch");


  switch (array)
  {
  case node_array_sol:
    call_gather_nodes(sol);
    break;
  case node_array_sol_increment:
    call_gather_nodes(sol_increment);
    break;
  case node_array_sol_residual:
    call_gather_nodes(sol_residual);
    break;
  case node_array_sol_mf:
    call_gather_nodes(sol_mf);
    break;
  default:
    dserror("Node array %d unknown", array);
  }

#undef call_gather_nodes
#undef pack_send_size_buf

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the mesh (element id, element type and connectivity) in
  some buffer to send them to their writing processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_mesh(BIN_OUT_CHUNK *chunk,
                          PARTDISCRET *actpdis,
                          DOUBLE *send_buf,
                          INT send_count,
                          INT *send_size_buf,
                          INT dst_first_id,
                          INT dst_num)
{
  INT i;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_mesh");
#endif

  dsassert(chunk->size_entry_length > 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  len = chunk->size_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      INT* ptr;
      INT j;

#ifdef PARALLEL
      send_size_buf[(len+1)*counter] = actele->Id_loc;
      ptr = &(send_size_buf[(len+1)*counter+1]);
#else
      ptr = &(send_size_buf[len*counter]);
#endif

      dsassert(actele->numnp <= len, "size entry too short");

      for (j=0; j<actele->numnp; ++j)
      {
        /* We have to store the local Id here, of course. */
        *ptr++ = actele->node[j]->Id_loc;
      }

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the node id and its coordinates in some buffer to send
  them to their writing processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_coords(BIN_OUT_CHUNK *chunk,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT len;
  INT slen;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_coords");
#endif

  dsassert(chunk->size_entry_length > 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == genprob.ndim, "invalid value entry length");

  len = chunk->value_entry_length;
  slen = chunk->size_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      DOUBLE* ptr;
      INT* size_ptr;
      INT j;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actnode->Id_loc;
      size_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[slen*counter]);
#endif

      /* first store the global node Id */
      size_ptr[node_variables.coords_size_Id] = actnode->Id;

      /* then store the number of elements to this node */
      size_ptr[node_variables.coords_size_numele] = actnode->numele;

      /* afterwards the local Ids of the elements connected to this node */
      for (j=0; j<actnode->numele; ++j)
      {
        size_ptr[j+node_variables.coords_size_eleid] = actnode->element[j]->Id_loc;
      }

      /* and finally store the node coordinates */
      ptr = &(send_buf[len*counter]);
      for (j=0; j<len; ++j)
      {
        *ptr++ = actnode->x[j];
      }

      counter += 1;
    }
  }

  dsassert(counter*len == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack any element specific parameters.

  We wanted to avoid to number our element variants. Thus there is not
  way to compare to elements and group them. Each element is an object
  and has to store all element parameters.

  The first entry in the integer chunk is the (global) element id. The
  second one is the ordinary ccarat element type. Then there is the
  dis type number and the number of nodes. All remaining entries are
  element type dependent.

  \author u.kue
  \date 11/04
  \sa out_pack_items
  \sa find_ele_param_item_length
*/
/*----------------------------------------------------------------------*/
static void out_pack_ele_params(BIN_OUT_CHUNK *chunk,
                                PARTDISCRET *actpdis,
                                DOUBLE *send_buf,
                                INT send_count,
                                INT *send_size_buf,
                                INT dst_first_id,
                                INT dst_num)
{
  INT i;
  INT size_len;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_ele_params");
#endif

  dsassert(chunk->size_entry_length >= 4, "invalid size entry length");
  dsassert(chunk->value_entry_length >= 0, "invalid value entry length");

  len = chunk->value_entry_length;
  size_len = chunk->size_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      DOUBLE* ptr;
      INT* size_ptr;
      ELEMENT_VARIABLES* el_vars = &element_variables;

#ifdef PARALLEL
      send_size_buf[(size_len+1)*counter] = actele->Id_loc;
      size_ptr = &(send_size_buf[(size_len+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[size_len*counter]);
#endif
      ptr = &(send_buf[len*counter]);

      /* This is commen to all elements */
      size_ptr[el_vars->ep_size_Id]     = actele->Id;
      size_ptr[el_vars->ep_size_eltyp]  = actele->eltyp;
      size_ptr[el_vars->ep_size_distyp] = actele->distyp;
      size_ptr[el_vars->ep_size_numnp]  = actele->numnp;

      switch (actele->eltyp)
      {
#ifdef D_SHELL8
      case el_shell8:
      {
        SHELL8* s8 = actele->e.s8;
        SHELL8_VARIABLES* vars = &shell8_variables;

        size_ptr[vars->ep_size_nGP0] = s8->nGP[0];
        size_ptr[vars->ep_size_nGP1] = s8->nGP[1];
        size_ptr[vars->ep_size_nGP2] = s8->nGP[2];
        size_ptr[vars->ep_size_nGP_tri] = s8->nGP_tri;
        size_ptr[vars->ep_size_forcetyp] = s8->forcetyp;

        ptr[vars->ep_value_sdc] = s8->sdc;

        break;
      }
#endif
#ifdef D_SHELL9
      case el_shell9:
      {
        SHELL9* s9 = actele->e.s9;
        SHELL9_VARIABLES* vars = &shell9_variables;

        size_ptr[vars->ep_size_nGP0] = s9->nGP[0];
        size_ptr[vars->ep_size_nGP1] = s9->nGP[1];
        size_ptr[vars->ep_size_nGP2] = s9->nGP[2];
        size_ptr[vars->ep_size_nGP_tri] = s9->nGP_tri;
        size_ptr[vars->ep_size_forcetyp] = s9->forcetyp;

        break;
      }
#endif
#ifdef D_BRICK1
      case el_brick1:
      {
        BRICK1* c1 = actele->e.c1;
        BRICK1_VARIABLES* vars = &brick1_variables;

        size_ptr[vars->ep_size_nGP0] = c1->nGP[0];
        size_ptr[vars->ep_size_nGP1] = c1->nGP[1];
        size_ptr[vars->ep_size_nGP2] = c1->nGP[2];
        size_ptr[vars->ep_size_stresstyp] = c1->stresstyp;

        break;
      }
#endif
#ifdef D_FLUID2
      case el_fluid2:
      {
        FLUID2* f2 = actele->e.f2;
        FLUID2_VARIABLES* vars = &fluid2_variables;

        size_ptr[vars->ep_size_nGP0] = f2->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f2->nGP[1];

        break;
      }
#endif
#ifdef D_FLUID2_PRO
      case el_fluid2_pro:
      {
        FLUID2_PRO* f2pro = actele->e.f2pro;
        FLUID2_PRO_VARIABLES* vars = &fluid2_pro_variables;

        size_ptr[vars->ep_size_nGP0] = f2pro->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f2pro->nGP[1];

        break;
      }
#endif
#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
      {
        FLUID3_PRO* f3pro = actele->e.f3pro;
        FLUID3_PRO_VARIABLES* vars = &fluid3_pro_variables;

        size_ptr[vars->ep_size_nGP0] = f3pro->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f3pro->nGP[1];
        size_ptr[vars->ep_size_nGP2] = f3pro->nGP[2];

        break;
      }
#endif
#ifdef D_FLUID2TU
      case el_fluid2_tu:
      {
        FLUID2* f2 = actele->e.f2;
        FLUID2TU_VARIABLES* vars = &fluid2tu_variables;

        size_ptr[vars->ep_size_nGP0] = f2->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f2->nGP[1];

        break;
      }
#endif
#ifdef D_FLUID3
      case el_fluid3:
      {
        FLUID3* f3 = actele->e.f3;
        FLUID3_VARIABLES* vars = &fluid3_variables;

        size_ptr[vars->ep_size_nGP0] = f3->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f3->nGP[1];
        size_ptr[vars->ep_size_nGP2] = f3->nGP[2];

        break;
      }
#endif
#ifdef D_FLUID3_F
      case el_fluid3_fast:
      {
        FLUID3* f3 = actele->e.f3;
        FLUID3_FAST_VARIABLES* vars = &fluid3_fast_variables;

        size_ptr[vars->ep_size_nGP0] = f3->nGP[0];
        size_ptr[vars->ep_size_nGP1] = f3->nGP[1];
        size_ptr[vars->ep_size_nGP2] = f3->nGP[2];

        break;
      }
#endif
#ifdef D_ALE
      case el_ale2:
      {
        ALE2* ale2 = actele->e.ale2;
        ALE2_VARIABLES* vars = &ale2_variables;

        size_ptr[vars->ep_size_nGP0] = ale2->nGP[0];
        size_ptr[vars->ep_size_nGP1] = ale2->nGP[1];

        break;
      }
#endif
#ifdef D_ALE
      case el_ale3:
      {
        ALE3* ale3 = actele->e.ale3;
        ALE3_VARIABLES* vars = &ale3_variables;

        size_ptr[vars->ep_size_nGP0] = ale3->nGP[0];
        size_ptr[vars->ep_size_nGP1] = ale3->nGP[1];
        size_ptr[vars->ep_size_nGP2] = ale3->nGP[2];

        break;
      }
#endif
#ifdef D_WALL1
      case el_wall1:
      {
        WALL1* w1 = actele->e.w1;
        WALL1_VARIABLES* vars = &wall1_variables;

        size_ptr[vars->ep_size_nGP0] = w1->nGP[0];
        size_ptr[vars->ep_size_nGP1] = w1->nGP[1];
        size_ptr[vars->ep_size_nGP2] = w1->nGP[2];
        size_ptr[vars->ep_size_nGP3] = w1->nGP[3];
        size_ptr[vars->ep_size_stresstyp] = w1->stresstyp;

        ptr[vars->ep_value_thick] = w1->thick;

        break;
      }
#endif
#ifdef D_BEAM3
      case el_beam3:
      {
        BEAM3* b3 = actele->e.b3;
        BEAM3_VARIABLES* vars = &beam3_variables;

        size_ptr[vars->ep_size_nGP0] = b3->nGP[0];

        break;
      }
#endif
#ifdef D_AXISHELL
      case el_axishell:
      {
        AXISHELL* saxi = actele->e.saxi;
        AXISHELL_VARIABLES* vars = &axishell_variables;

        /* Is this an element parameter? Not really... */
        ptr[vars->ep_value_thick0] = saxi->thick[0];
        ptr[vars->ep_value_thick1] = saxi->thick[1];

        break;
      }
#endif
#ifdef D_INTERF
      case el_interf:
      {
        INTERF* interf = actele->e.interf;
        INTERF_VARIABLES* vars = &interf_variables;

        size_ptr[vars->ep_size_nGP] = interf->nGP;
        size_ptr[vars->ep_size_stresstyp] = interf->stresstyp;

        ptr[vars->ep_value_thick] = interf->thick;

        break;
      }
#endif
#ifdef D_WALLGE
      case el_wallge:
      {
        WALLGE* wallge = actele->e.wallge;
        WALLGE_VARIABLES* vars = &wallge_variables;

        size_ptr[vars->ep_size_nGP0] = wallge->nGP[0];
        size_ptr[vars->ep_size_nGP1] = wallge->nGP[1];
        size_ptr[vars->ep_size_nGP2] = wallge->nGP[2];
        size_ptr[vars->ep_size_nGP3] = wallge->nGP[3];
        size_ptr[vars->ep_size_stresstyp] = wallge->stresstyp;

        ptr[vars->ep_value_thick] = wallge->thick;

        break;
      }
#endif
      default:
        dserror("element type %d unsupported", actele->eltyp);
      }

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the displacements to send them to their writing processor.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_displacement(BIN_OUT_CHUNK *chunk,
                                  INT place,
                                  PARTDISCRET *actpdis,
                                  DOUBLE *send_buf,
                                  INT send_count,
                                  INT *send_size_buf,
                                  INT dst_first_id,
                                  INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_displacement");
#endif

  ndim = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j)
  {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT k;
      DOUBLE *ptr = actnode->sol.a.da[place];

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim >= ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      for (k=0; k<ndim; ++k)
      {
        send_buf[ndim*counter+k] = *ptr++;
      }
      counter += 1;
    }
  }
  dsassert(counter*ndim == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


#ifdef D_SHELL8

/*----------------------------------------------------------------------*/
/*!
  \brief Pack the directors.

  There are three director coordinates per node per
  element. Additionally we store one thickness value at each node.

  Using these values it's possible to visualize shell8 elements using
  bricks.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell8_director(BIN_OUT_CHUNK *chunk,
                                     PARTDISCRET *actpdis,
                                     DOUBLE *send_buf,
                                     INT send_count,
                                     INT *send_size_buf,
                                     INT dst_first_id,
                                     INT dst_num)
{
  INT i;
  INT k;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_shell8_director");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      SHELL8* s8;
      DOUBLE *dst_ptr;
      INT numnp;

      s8 = actele->e.s8;
      numnp = actele->numnp;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      for (k=0; k<numnp; ++k)
      {
        *dst_ptr++ = s8->a3ref.a.da[0][k];
        *dst_ptr++ = s8->a3ref.a.da[1][k];
        *dst_ptr++ = s8->a3ref.a.da[2][k];

        /*
         * We could read the per node thickness here as well but this
         * is what the old ccarat output does. */
        *dst_ptr++ = s8->thick;
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


#ifdef D_SHELL9


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the node coordinates of the shell9 element's inner nodes.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell9_coords(BIN_OUT_CHUNK *chunk,
                                   PARTDISCRET *actpdis,
                                   DOUBLE *send_buf,
                                   INT send_count,
                                   INT *send_size_buf,
                                   INT dst_first_id,
                                   INT dst_num)
{
  INT i;
  INT k;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_shell9_coords");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      SHELL9* s9;
      DOUBLE *dst_ptr;

      s9 = actele->e.s9;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      /* this is ``s9_out_gid_allcoords_unsmo`` in disguise */

      if (actele->numnp == 4)
      {
        for (k=0; k<actele->numnp; k++)
        {
          NODE *actnode;
          INT klay;
          INT jlay;
          INT l;

          /* director a3 of current layer in ref. config */
          DOUBLE a3ref_l[3];

          /* coordinats */
          DOUBLE x[3],x_u[3],x_o[3];

          actnode = actele->node[k];

          /* get director of node
             NOTE: in the undeformed geometrie, the directors of the
             different kinematic layers are equal to the director in the
             total reference layer */
          a3ref_l[0]= s9->a3ref.a.da[0][k];
          a3ref_l[1]= s9->a3ref.a.da[1][k];
          a3ref_l[2]= s9->a3ref.a.da[2][k];

          /* calculate the coordinates of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the coordinates to values of middle surface */
          for (l=0; l<3; l++)
            x[l] = actnode->x[l];

          /* continuity matrix for current kin layer at bot */
          for (jlay=0; jlay<s9->num_klay; jlay++)
          {
            DOUBLE klayhgt;
            DOUBLE hl;
            DOUBLE e3;
            DOUBLE zeta;

            klayhgt = s9->klayhgt[jlay];
            hl      = (klayhgt/100)*s9->thick;
            /* hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

            /* calculate coordinates */
            for (l=0; l<3; l++)
              x[l] += zeta*hl*a3ref_l[l];
          }

          /* write coordinates of first kinematic layer at bottom */
          *dst_ptr++ = x[0];
          *dst_ptr++ = x[1];
          *dst_ptr++ = x[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            x_u[l] = x[l];

          /* calculate the coordinates of the top surface of the each kin layer */
          for (klay=0; klay<s9->num_klay; klay++)
          {
            INT num_mlay;
            DOUBLE* mlayhgt;
            DOUBLE sum_hgt;
            INT mlay;

            /* initialize the coordinates to values of middle surface */
            for (l=0; l<3; l++)
              x[l] = actnode->x[l];

            /* continuity matrix for aktual kin layer at top */
            for (jlay=0; jlay<s9->num_klay; jlay++)
            {
              DOUBLE klayhgt;
              DOUBLE hl;
              DOUBLE e3;
              DOUBLE zeta;

              klayhgt = s9->klayhgt[jlay];
              hl      = (klayhgt/100)*s9->thick;
              /*hl      = 0.5*hl; */  /*norm of director is half of kinematic layer hight*/
              hl      = A3FAC_SHELL9*hl;
              e3 = +1.0; /*top*/
              zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

              /* calculate coordinates */
              for (l=0; l<3; l++)
                x[l] += zeta*hl*a3ref_l[l];
            }

            /* calculate the coordinates of the material layers of
             * aktual kinematic layer */
            num_mlay = s9->kinlay[klay].num_mlay;
            mlayhgt  = s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save coordinates for interpolation for the material layers */
            for (l=0; l<3; l++)
              x_o[l] = x[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt  = 0.0; /*initialize the total hgt*/

            /* only if more than one material layer to this kinematic layer */
            for (mlay=1; mlay<num_mlay; mlay++)
            {
              sum_hgt += mlayhgt[mlay-1];
              for (l=0; l<3; l++)
                x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);

              /* write the bottom coordinates of this material layer */
              *dst_ptr++ = x[0];
              *dst_ptr++ = x[1];
              *dst_ptr++ = x[2];
            }

            /* write the top coordinates of this kinematic layer */
            *dst_ptr++ = x_o[0];
            *dst_ptr++ = x_o[1];
            *dst_ptr++ = x_o[2];

            for (l=0; l<3; l++)
              x_u[l] = x_o[l];
          }
        }
      }
      else if ((actele->numnp == 8) || (actele->numnp == 9))
      {
        for (k=0; k<actele->numnp; k++)
        {
          NODE *actnode;
          INT klay;
          INT jlay;
          INT l;

          /* director a3 of current layer in ref. config */
          DOUBLE a3ref_l[3];

          /* coordinats */
          DOUBLE x[3],x_u[3],x_o[3];

          actnode = actele->node[k];

          /* get director of node
             NOTE: in the undeformed geometrie, the directors of the
             different kinematic layers are equal to the director in the
             total reference layer */
          a3ref_l[0]= s9->a3ref.a.da[0][k];
          a3ref_l[1]= s9->a3ref.a.da[1][k];
          a3ref_l[2]= s9->a3ref.a.da[2][k];

          /* calculate the coordinates of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the coordinates to values of middle surface */
          for (l=0; l<3; l++)
            x[l] = actnode->x[l];

          /* continuity matrix for aktual kin layer at bot */
          for (jlay=0; jlay<s9->num_klay; jlay++)
          {
            DOUBLE klayhgt;
            DOUBLE hl;
            DOUBLE e3;
            DOUBLE zeta;

            klayhgt = s9->klayhgt[jlay];
            hl      = (klayhgt/100)*s9->thick;
            /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

            /* calculate coordinates */
            for (l=0; l<3; l++)
              x[l] += zeta*hl*a3ref_l[l];
          }

          /* write coordinates of first kinematic layer at bottom */
          *dst_ptr++ = x[0];
          *dst_ptr++ = x[1];
          *dst_ptr++ = x[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            x_u[l] = x[l];

          /* calculate the coordinates of the top surface of the each
             kin layer */
          for (klay=0; klay<s9->num_klay; klay++)
          {
            INT num_mlay;
            DOUBLE* mlayhgt;
            DOUBLE sum_hgt;
            DOUBLE sum_hgt_mid;
            INT mlay;

            /* initialize the coordinates to values of middle surface */
            for (l=0; l<3; l++)
              x[l] = actnode->x[l];

            /* continuity matrix for aktual kin layer at top */
            for (jlay=0; jlay<s9->num_klay; jlay++)
            {
              DOUBLE klayhgt;
              DOUBLE hl;
              DOUBLE e3;
              DOUBLE zeta;

              klayhgt = s9->klayhgt[jlay];
              hl      = (klayhgt/100)*s9->thick;
              /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
              hl      = A3FAC_SHELL9*hl;   /*norm of director is half of kinematic layer hight*/
              e3 = +1.0; /*top*/
              zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);

              /* calculate coordinates */
              for (l=0; l<3; l++)
                x[l] += zeta*hl*a3ref_l[l];
            }

            /* calculate the coordinates of the material layers of
             * aktual kinematic layer */
            num_mlay = actele->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save coordinates for interpolation for the material layers */
            for (l=0; l<3; l++)
              x_o[l] = x[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt     = 0.0; /*initialize the total hgt*/
            sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
            for (mlay=0; mlay<num_mlay; mlay++)
            {

              /* if more than one material layer to this kinematic layer */
              if (mlay > 0)
              {
                sum_hgt += mlayhgt[mlay-1];
                for (l=0; l<3; l++)
                  x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);

                /* write the bottom coordinates of this material layer */
                *dst_ptr++ = x[0];
                *dst_ptr++ = x[1];
                *dst_ptr++ = x[2];
              }

              /* a midnode has to be calculated */
              sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
              for (l=0; l<3; l++)
                x[l] = x_u[l] + (sum_hgt_mid/100.)*(x_o[l]-x_u[l]);

              /* write the middle coordinates of this material layer */
              *dst_ptr++ = x[0];
              *dst_ptr++ = x[1];
              *dst_ptr++ = x[2];
            }

            *dst_ptr++ = x_o[0];
            *dst_ptr++ = x_o[1];
            *dst_ptr++ = x_o[2];

            for (l=0; l<3; l++)
              x_u[l] = x_o[l];
          }
        }
      }
      else {
        dserror("unknown element type for shell9 element");
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the displacements of all the layers of a shell9 element.

  Normally the displacements belong to nodes. But these are internal
  element displacements that might be put to the nodes of a

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_shell9_displacement(BIN_OUT_CHUNK *chunk,
                                         INT place,
                                         PARTDISCRET *actpdis,
                                         DOUBLE *send_buf,
                                         INT send_count,
                                         INT *send_size_buf,
                                         INT dst_first_id,
                                         INT dst_num)
{
  INT i, counter, len;

#ifdef DEBUG
  dstrc_enter("out_pack_shell9_displacement");
#endif

  len = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      SHELL9* s9;
      DOUBLE *dst_ptr;
      DOUBLE dis[3], dis_u[3], dis_o[3]; /* displacement x,y,z */
      INT j;
      INT k;
      INT l;
      INT klay;

      s9 = actele->e.s9;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      /* this is ``s9_out_gid_sol_dis_unsmo`` in disguise */

      /* The purpose is to calculate all layer nodes'
       * displacements. That is each node reappears in each layer and
       * needs three displacement values. The difference between the
       * four node versions and the others is that with the four node
       * versions no nodes in the middle of a layer are generated.
       *
       * Please note that we're talking about node values and yet we
       * store element data. Each element has its inner nodes. It's up
       * to the filter to bring inner nodes of consecutive elements
       * together (smoothed view) or not (unconnected view). */

      if (actele->numnp == 4)
      {
        for (k=0; k<actele->numnp; k++)
        {
          NODE *actnode;

          actnode = actele->node[k];

          /* calculate the displacements of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the displacements of this layer to zero */
          for (l=0; l<3; l++)
            dis[l]=0.0;

          /* calculate the relative displacements of this layer */
          for (j=0; j<s9->num_klay; j++)
          {
            DOUBLE e3 = -1.0; /*bottom*/
            DOUBLE zeta = s9con(e3,s9->num_klay,klay,j,1.0);
            for (l=0; l<3; l++)
              dis[l] += zeta*actnode->sol.a.da[place][l+3*(j+1)];
          }

          /* add the displacement of reference layer to relative
           * displacement */
          for (l=0; l<3; l++)
            dis[l] += actnode->sol.a.da[place][l];

          /* write the displacements of first kinematic layer at bottom */
          *dst_ptr++ = dis[0];
          *dst_ptr++ = dis[1];
          *dst_ptr++ = dis[2];

          /* save coordinates for interpolation for the material layers */
          for (l=0; l<3; l++)
            dis_u[l] = dis[l];

          /* calculate the displacements of the top surface of the
           * each kin layer */
          for (klay=0; klay<s9->num_klay; klay++)
          {
            DOUBLE sum_hgt;
            INT mlay;
            INT num_mlay;
            DOUBLE* mlayhgt;

            /* initialize the displacements of this layer to zero */
            for (l=0; l<3; l++)
              dis[l]=0.0;

            /* calculate the relative displacements of this layer */
            for (j=0; j<s9->num_klay; j++)
            {
              DOUBLE e3 = +1.0; /*top*/
              DOUBLE zeta = s9con(e3,s9->num_klay,klay,j,1.0);
              for (l=0; l<3; l++)
                dis[l] += zeta*actnode->sol.a.da[place][l+3*(j+1)];
            }

            /* add the displacement of reference layer to relative
             * displacement */
            for (l=0; l<3; l++)
              dis[l] += actnode->sol.a.da[place][l];

            /* calculate the displacements of the material layers of
             * aktual kinematic layer */
            num_mlay = s9->kinlay[klay].num_mlay;
            mlayhgt  = s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save displacements for interpolation for the material layers */
            for (l=0; l<3; l++)
              dis_o[l] = dis[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */

            /* initialize the total hgt */
            sum_hgt = 0.0;

            /* only if more than one material layer to this kinematic layer */
            for (mlay=1; mlay<num_mlay; mlay++)
            {
              sum_hgt += mlayhgt[mlay-1];
              for (l=0; l<3; l++)
                dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);

              /* write the bottom displacement of this material layer */
              *dst_ptr++ = dis[0];
              *dst_ptr++ = dis[1];
              *dst_ptr++ = dis[2];
            }

            /* write the top displacement of this kinematic layer */
            *dst_ptr++ = dis_o[0];
            *dst_ptr++ = dis_o[1];
            *dst_ptr++ = dis_o[2];

            for (l=0; l<3; l++)
              dis_u[l] = dis_o[l];
          }
        }
      }
      else if ((actele->numnp == 8) || (actele->numnp == 9))
      {
        for (k=0; k<actele->numnp; k++)
        {
          NODE *actnode;
          INT jlay;

          actnode = actele->node[k];

          /* calculate the displacements of the bottom surface of the
           * first kinematic layer */
          klay = 0;

          /* initialize the displacements of this layer to zero */
          for (l=0; l<3; l++)
            dis[l]=0.0;

          /* calculate the relative displacements of this layer */
          for (jlay=0; jlay<s9->num_klay; jlay++)
          {
            DOUBLE e3   = -1.0; /*bottom*/
            DOUBLE zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);
            for (l=0; l<3; l++)
              dis[l] += zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
          }

          /* add the displacement of reference layer to relative
           * displacement */
          for (l=0; l<3; l++)
            dis[l] += actnode->sol.a.da[place][l];

          /* write the displacements of first kinematic layer at bottom */
          *dst_ptr++ = dis[0];
          *dst_ptr++ = dis[1];
          *dst_ptr++ = dis[2];

          /* save displacements for interpolation for the material layers */
          for (l=0; l<3; l++)
            dis_u[l] = dis[l];

          /* calculate the displacements of the top surface of the
           * each kin layer */
          for (klay=0; klay<s9->num_klay; klay++)
          {
            DOUBLE sum_hgt;
            DOUBLE sum_hgt_mid;
            INT mlay;
            INT num_mlay;
            DOUBLE* mlayhgt;

            /* initialize the displacements of this layer to zero */
            for(l=0; l<3; l++)
              dis[l]=0.0;

            /* calculate the relative displacements of this layer */
            for (jlay=0; jlay<s9->num_klay; jlay++)
            {
              DOUBLE e3   = +1.0; /*top*/
              DOUBLE zeta = s9con(e3,s9->num_klay,klay,jlay,1.0);
              for (l=0; l<3; l++)
                dis[l] += zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
            }

            /* add the displacement of reference layer to relative
             * displacement */
            for (l=0; l<3; l++)
              dis[l] += actnode->sol.a.da[place][l];

            /* calculate the displacements of the material layers of
             * aktual kinematic layer */
            num_mlay = actele->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */

            /* save displacements for interpolation for the material layers */
            for (l=0; l<3; l++)
              dis_o[l] = dis[l];

            /* interpolate for the material layers within the actual
             * kinematic layer */
            sum_hgt     = 0.0; /*initialize the total hgt*/
            sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/

            for (mlay=0; mlay<num_mlay; mlay++)
            {

              /* if more than one material layer to this kinematic layer */
              if (mlay > 0)
              {
                sum_hgt += mlayhgt[mlay-1];
                for (l=0; l<3; l++)
                  dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);

                /* write the bottom displacement of this material layer */
                *dst_ptr++ = dis[0];
                *dst_ptr++ = dis[1];
                *dst_ptr++ = dis[2];
              }

              /* the displacement at midnode has to be calculated */
              sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
              for (l=0; l<3; l++)
                dis[l] = dis_u[l] + (sum_hgt_mid/100.)*(dis_o[l]-dis_u[l]);

              /* write the mid displacement of this material layer  */
              *dst_ptr++ = dis[0];
              *dst_ptr++ = dis[1];
              *dst_ptr++ = dis[2];
            }

            /* write the top displacement of this kinematic layer */
            *dst_ptr++ = dis_o[0];
            *dst_ptr++ = dis_o[1];
            *dst_ptr++ = dis_o[2];

            for (l=0; l<3; l++)
              dis_u[l] = dis_o[l];
          }
        }
      }
      else
      {
        dserror("unknown element type for shell9 element");
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the velocities in some buffer to send them to their
  writing processor.

  There are two or three velocity values per node depending on the
  number of dimensions.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_velocity(BIN_OUT_CHUNK *chunk,
                              INT place,
                              PARTDISCRET *actpdis,
                              DOUBLE *send_buf,
                              INT send_count,
                              INT *send_size_buf,
                              INT dst_first_id,
                              INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_velocity");
#endif

  ndim = chunk->value_entry_length;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j)
  {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT k;
      DOUBLE *ptr = actnode->sol.a.da[place];

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim >= ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      for (k=0; k<ndim; ++k)
      {
        send_buf[ndim*counter+k] = *ptr++;
      }
      counter += 1;
    }
  }
  dsassert(counter*ndim == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the pressure in some buffer to send them to their
  writing processor.

  The point is that there is just one pressure value per node. It
  lifes in the node's sol array in the same row as the velocity
  values. That's why we need to get the current problem's number of
  dimensions.

  \author u.kue
  \date 08/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_pressure(BIN_OUT_CHUNK *chunk,
                              INT place,
                              PARTDISCRET *actpdis,
                              DOUBLE *send_buf,
                              INT send_count,
                              INT *send_size_buf,
                              INT dst_first_id,
                              INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_pressure");
#endif

  ndim = genprob.ndim;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  for (j=0; j<actpdis->numnp; ++j)
  {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {

      dsassert((actnode->sol.fdim > place) &&
               (actnode->sol.sdim > ndim), "sol array too small");

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      send_buf[counter] = actnode->sol.a.da[place][ndim];

      counter += 1;
    }
  }
  dsassert(counter == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the averaged pressure in some buffer to send them to their
  writing processor.

  The projection method uses discontinuous pressures that are
  difficult to plot. So here we extrapolate the pressure to the nodes
  and average them.

  \author u.kue
  \date 01/06
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_average_pressure(BIN_OUT_CHUNK *chunk,
                                      INT place,
                                      PARTDISCRET *actpdis,
                                      DOUBLE *send_buf,
                                      INT send_count,
                                      INT *send_size_buf,
                                      INT dst_first_id,
                                      INT dst_num)
{
  INT j, counter, ndim;

#ifdef DEBUG
  dstrc_enter("out_pack_average_pressure");
#endif

  ndim = genprob.ndim;

  /* gather the nodes to be send to processor i */
  counter = 0;
  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  for (j=0; j<actpdis->numnp; ++j)
  {
    NODE* actnode = actpdis->node[j];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT k;

#ifdef PARALLEL
      send_size_buf[counter] = actnode->Id_loc;
#endif

      send_buf[counter] = 0;

      for (k=0; k<actnode->numele; ++k)
      {
        INT l;
        DOUBLE el_press=0;
        ELEMENT* actele;
        actele = actnode->element[k];

        /* find element local index of node */
        for (l=0; l<actele->numnp; ++l)
        {
          NODE* node;
          node = actele->node[l];
          if (actnode==node)
          {
            break;
          }
        }
        if (l==actele->numnp)
          dserror("Node %d not in element %d. Suspect!", actnode->Id, actele->Id);

        switch (actele->eltyp)
        {
#ifdef D_FLUID2_PRO
        case el_fluid2_pro:
          f2pro_addnodepressure(actele, l, &el_press);
          break;
#endif
#ifdef D_FLUID3_PRO
        case el_fluid3_pro:
          f3pro_addnodepressure(actele, l, &el_press);
          break;
#endif
        default:
          dserror("Unsupported element type %d for averaged pressure", actele->eltyp);
        }

        send_buf[counter] += el_press;
      }
      counter += 1;
    }
  }
  dsassert(counter == send_count, "node pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the element's stresses.

  Here the element based stresses are handled, that is the stresses at
  the gauss points. We visit all elements and copy its stress array to
  the correct send_buf place. We copy in such a way that all stress
  values that belong to one gauss point will be close to each other in
  the send_buf. The internal ccarat arrays do it the other way
  round. There each gauss point corresponds to a column and each row
  has its meaning. Thus the C memory layout will bring all values with
  the same meaning together.

  This is highly element specific. That is the filter must know what
  these numbers mean. Each element is different.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_stress(BIN_OUT_CHUNK *chunk,
                            INT place,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT j;
  INT k;
  INT len;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_stress");
#endif

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");

  len = chunk->value_entry_length;

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      DOUBLE **stress;
      DOUBLE *dst_ptr;
      INT rows;
      INT cols;

      /* Where this element's values are to go. */
      dst_ptr = &(send_buf[len*counter]);

      switch (actele->eltyp)
      {
#ifdef D_WALL1
      case el_wall1:            /* 2D plane stress - plane strain element */

        /* there is just one source for all variants of wall elements */
        stress = actele->e.w1->stress_GP.a.d3[place];

        if (actele->distyp == tri3)
        {
          *dst_ptr++ = stress[0][0];
          *dst_ptr++ = stress[1][0];
          *dst_ptr++ = stress[2][0];
          *dst_ptr++ = stress[3][0];
        }
        else if (actele->distyp == tri6)
        {
          dserror("yet to be done");
        }
        else if ((actele->distyp == quad4) ||
                 (actele->distyp == quad8) ||
                 (actele->distyp == quad9))
        {
          rows = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];
          cols = 9;

          /*
           * This is bad. Only the first four values are stresses
           * (only a four gauss point version is supported here). The
           * remaining values would better be put to their own
           * chunk. */
          if (actele->e.w1->elewa != NULL)
          {
            for (j=0; j<rows; j++)
            {
              *dst_ptr++ = stress[0][j];
              *dst_ptr++ = stress[1][j];
              *dst_ptr++ = stress[2][j];
              *dst_ptr++ = stress[3][j];
              *dst_ptr++ = stress[4][j];
              *dst_ptr++ = stress[5][j];
              *dst_ptr++ = stress[6][j];
              *dst_ptr++ = actele->e.w1->elewa[0].ipwa[j].damage;
              *dst_ptr++ = actele->e.w1->elewa[0].ipwa[j].aequistrain;
            }
          }
          else
          {
            for (j=0; j<rows; j++)
            {
              *dst_ptr++ = stress[0][j];
              *dst_ptr++ = stress[1][j];
              *dst_ptr++ = stress[2][j];
              *dst_ptr++ = stress[3][j];
              *dst_ptr++ = stress[4][j];
              *dst_ptr++ = stress[5][j];
              *dst_ptr++ = stress[6][j];
              *dst_ptr++ = 0;
              *dst_ptr++ = 0;
            }
          }
        }
        else
        {
          dserror("distyp %d unsupported", actele->distyp);
        }

        break;
#endif
#ifdef D_BEAM3
      case el_beam3:            /* structural 3D-beam element */

        stress = actele->e.b3->force_GP.a.d3[place];
        cols = 6;

        rows = actele->e.b3->nGP[0];

        /* copy the values */
        for (j=0; j<rows; j++)
        {
          *dst_ptr++ = stress[0][j];
          *dst_ptr++ = stress[1][j];
          *dst_ptr++ = stress[2][j];
          *dst_ptr++ = stress[3][j];
          *dst_ptr++ = stress[4][j];
          *dst_ptr++ = stress[5][j];
        }

        break;
#endif
#ifdef D_INTERF
      case el_interf:           /* 1D interface element (combination only with wall) */

        stress = actele->e.interf->stress_GP.a.d3[place];
        cols = 5;
        rows = actele->e.interf->nGP;

        /* copy the values */
        for (j=0; j<rows; j++)
        {
          *dst_ptr++ = stress[0][j];
          *dst_ptr++ = stress[1][j];
          *dst_ptr++ = stress[2][j];
          *dst_ptr++ = actele->e.interf->elewa[0].ipwa[j].dn;
          *dst_ptr++ = actele->e.interf->elewa[0].ipwa[j].dt;
        }

        break;
#endif
#ifdef D_WALLGE
      case el_wallge:           /* gradient enhanced wall element */
        /*stress = actele->e.wallge->stress_GP.a.d3[place];*/

        if ((actele->distyp == tri3) || (actele->distyp == tri6))
        {
          dserror("yet to be done");
        }
        else if ((actele->distyp == quad4) ||
                 (actele->distyp == quad8) ||
                 (actele->distyp == quad9))
        {
          rows = actele->e.wallge->nGP[0]*actele->e.wallge->nGP[1];
        }
        else
        {
          dserror("distyp %d unsupported", actele->distyp);
        }

        /* copy the values */
        for (j=0; j<rows; j++)
        {
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].sig[0];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].sig[1];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].sig[2];
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].damage;
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].aequistrain;
          *dst_ptr++ = actele->e.wallge->elwa[0].iptwa[j].aequistrain_nl;
        }

        break;
#endif
#ifdef D_AXISHELL
      case el_axishell:         /* 1D axisymmetrical shell element */
        /*
         * There were additional (?) displacements. But these are
         * nodal values. They are written at another place. */

        stress = actele->e.saxi->stress_GP.a.d3[place];
        rows = 5;

        for (j=0; j<rows; j++)
        {
          *dst_ptr++ = stress[j][0];
        }
        break;
#endif
#ifdef D_FLUID3
      case el_fluid3:           /* 3D fluid element */
        /* node based. ignore here. */
        break;
#endif
#ifdef D_SHELL9
      case el_shell9:
      {         /* multi layer shell element */
        SHELL9* s9 = actele->e.s9;

        stress = s9->gp_stress.a.da;

        /* In general send_count can be bigger that this element's
         * contribution. In a homogeneous mesh these sizes must be
         * equal, of course, but lets not rely on the homogeneousness
         * when we don't need to. */
        rows = s9->gp_stress.fdim;
        cols = s9->gp_stress.sdim;

        dsassert(rows*cols <= send_count, "shell9 entry too small");

        /* copy the values */
        /*
         * This is special! We keep the memory layout of the internal
         * ccarat array. The filters seem to be easier that way. */
        for (j=0; j<rows; j++)
        {
          for (k=0; k<cols; k++)
          {
            *dst_ptr++ = stress[j][k];
          }
        }
        break;
      }
#endif
#ifdef D_SHELL8
      case el_shell8:
      {         /* 7 parameter shell element */
        SHELL8* s8 = actele->e.s8;

        /*
         * We rely on the fact that the first index is not used. It
         * steams back to the old idea of keeping all results in
         * memory... */
        stress = s8->forces.a.d3[0];

        /* thus the second dimension gives the number of rows. */
        rows = s8->forces.sdim;

        /*
         * Are there triangular versions those number of gauss points
         * cannot be calculated like this? */
        cols = s8->nGP[0]*s8->nGP[1]*s8->nGP[2];

        dsassert(rows == 18, "shell8 changed but binary output not updated");
        dsassert(rows*cols <= send_count, "shell8 entry too small");

        /* copy the values */
        /*
         * The inner loop must loop the rows in order to keep the
         * values of one gauss point together. */
        for (k=0; k<cols; k++)
        {
          for (j=0; j<rows; j++)
          {
            *dst_ptr++ = stress[j][k];
          }
        }
        break;
      }
#endif
#ifdef D_BRICK1
      case el_brick1:
      {         /* structural brick element */
        BRICK1* c1 = actele->e.c1;

        /*
         * We rely on the fact that the first index is not used. It
         * steams back to the old idea of keeping all results in
         * memory... */
        stress = c1->stress_GP.a.d3[0];

        /* thus the second dimension gives the number of rows. */
        rows = c1->stress_GP.sdim;

        /*
         * Are there triangular versions those number of gauss points
         * cannot be calculated like this? */
        cols = c1->nGP[0]*c1->nGP[1]*c1->nGP[2];

        dsassert(rows == 27, "brick1 changed but binary output not updated");
        dsassert(rows*cols <= send_count, "brick1 entry too small");

        /* copy the values */
        /*
         * The inner loop must loop the rows in order to keep the
         * values of one gauss point together. */
        for (k=0; k<cols; k++)
        {
          for (j=0; j<rows; j++)
          {
            *dst_ptr++ = stress[j][k];
          }
        }
        break;
      }
#endif
      case el_fluid2:           /* 2D fluid element */
      case el_fluid2_pro:       /* 2D fluid element */
      case el_fluid2_tu:        /* 2D fluid element for turbulence */
      case el_ale2:             /* 2D pseudo structural ale element */
      case el_ale3:             /* 3D pseudo structural ale element */
      default:
        dserror("element based stress output not supported for element type %d", actele->eltyp);
      }

#ifdef PARALLEL
      send_size_buf[counter] = actele->Id_loc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the element's processor.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_domain(BIN_OUT_CHUNK *chunk,
                            PARTDISCRET *actpdis,
                            DOUBLE *send_buf,
                            INT send_count,
                            INT *send_size_buf,
                            INT dst_first_id,
                            INT dst_num)
{
  INT i;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_domain");
#endif

  dsassert(chunk->size_entry_length == 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  counter = 0;
  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {

#ifdef PARALLEL
      send_size_buf[2*counter  ] = actele->Id_loc;
      send_size_buf[2*counter+1] = actele->proc;
#else
      send_size_buf[counter] = actele->proc;
#endif

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack an array of distributed vectors.

  The distributed vector is supposed to contain one value per dof. We
  store the dofs in ascending order, the order in the distributed
  vector at hand, however, is solver dependent.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dist_vector(BIN_OUT_CHUNK *chunk,
                                 PARTDISCRET *actpdis,
                                 DOUBLE *send_buf,
                                 INT send_count,
                                 INT *send_size_buf,
                                 INT dst_first_id,
                                 INT dst_num)
{
  INT i;
  INT len;
  INT counter;
  SPARSE_TYP sysarray_typ;
  SPARSE_ARRAY sysarray;

#ifdef DEBUG
  dstrc_enter("out_pack_dist_vector");
#endif

  /* for convenience */
  sysarray_typ = *chunk->field->sysarray_typ;
  sysarray     = *chunk->field->sysarray;

  dsassert(chunk->size_entry_length == 0, "invalid size entry length");
  len = chunk->value_entry_length;

#ifdef PARALLEL
#define pack_send_size_buf send_size_buf[counter] = dof;
#else
#define pack_send_size_buf
#endif

#define boilerplate_copying_code                        \
  if ((dof >= dst_first_id) &&                          \
      (dof < dst_first_id+dst_num))                     \
  {                                                     \
    INT j;                                              \
    DOUBLE* dst_ptr;                                    \
    dst_ptr = &(send_buf[len*counter]);                 \
    for (j=0; j<len; ++j)                               \
    {                                                   \
      *dst_ptr++ = chunk->vectors[j].vec.a.dv[i];       \
    }                                                   \
    pack_send_size_buf;                                 \
    counter += 1;                                       \
  }

  counter = 0;

  /* Here we have the solver dependency. This is
   * ``solserv_reddistvec`` in disguise. */
  switch (sysarray_typ)
  {

#ifdef TRILINOS_PACKAGE
  case trilinos:
    {
      const int numeq = sysarray.trilinos->numeq;
      int* update = sysarray.trilinos->update.a.iv;
      for (i=0; i<numeq; ++i)
      {
        INT dof = update[i];
        boilerplate_copying_code;
      }
    }
    break;
#endif

#ifdef AZTEC_PACKAGE
  case msr:
    for (i=0; i<sysarray.msr->numeq; ++i)
    {
      INT dof = sysarray.msr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef HYPRE_PACKAGE
  case parcsr:
  {
    INT rank = chunk->field->actintra->intra_rank;
    for (i=0; i<sysarray.parcsr->numeq; ++i)
    {
      INT dof = sysarray.parcsr->update.a.ia[rank][i];
      boilerplate_copying_code;
    }
    break;
  }
#endif

#ifdef PARSUPERLU_PACKAGE
  case ucchb:
    for (i=0; i<sysarray.ucchb->numeq; ++i)
    {
      INT dof = sysarray.ucchb->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case dense:
    for (i=0; i<sysarray.dense->numeq; ++i)
    {
      INT dof = sysarray.dense->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

#ifdef MLIB_PACKAGE
  case mds:
    for (i=0; i<sysarray.mds->numeq; ++i)
    {
      INT dof = i;
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef MUMPS_PACKAGE
  case rc_ptr:
    for (i=0; i<sysarray.rc_ptr->numeq; ++i)
    {
      INT dof = sysarray.rc_ptr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef SPOOLES_PACKAGE
  case spoolmatrix:
    for (i=0; i<sysarray.spo->numeq; ++i)
    {
      INT dof = sysarray.spo->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

#ifdef UMFPACK
  case ccf:
    for (i=0; i<sysarray.ccf->numeq; ++i)
    {
      INT dof = sysarray.ccf->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case skymatrix:
    for (i=0; i<sysarray.sky->numeq; ++i)
    {
      INT dof = sysarray.sky->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

#ifdef MLPCG
  case bdcsr:
    for (i=0; i<sysarray.bdcsr->numeq; ++i)
    {
      INT dof = sysarray.bdcsr->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;
#endif

  case oll:
    for (i=0; i<sysarray.oll->numeq; ++i)
    {
      INT dof = sysarray.oll->update.a.iv[i];
      boilerplate_copying_code;
    }
    break;

  default:
    dserror("Unknown type %d of system matrix", sysarray_typ);
    break;
  }

#undef boilerplate_copying_code

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the restart data of all the elements in this
  discretization.

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \author u.kue
  \date 09/04
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_restart_element(BIN_OUT_CHUNK *chunk,
                                     PARTDISCRET *actpdis,
                                     DOUBLE *send_buf,
                                     INT send_count,
                                     INT *send_size_buf,
                                     INT dst_first_id,
                                     INT dst_num)
{
  INT el;
  INT len;
  INT slen;
  INT counter;

#ifdef DEBUG
  dstrc_enter("out_pack_restart_element");
#endif

  len = chunk->value_entry_length;
  slen = chunk->size_entry_length;

  counter = 0;
  for (el=0; el<actpdis->numele; ++el)
  {
    ELEMENT* actele = actpdis->element[el];
    if ((actele->Id_loc >= dst_first_id) &&
        (actele->Id_loc < dst_first_id+dst_num))
    {
      DOUBLE* dst_ptr;
      INT* size_dst_ptr;
      INT j;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actele->Id_loc;
      size_dst_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_dst_ptr = &(send_size_buf[slen*counter]);
#endif

      dst_ptr = &(send_buf[len*counter]);

      /* Now here we have all the different elements. */
      switch (actele->eltyp)
      {

#ifdef D_SHELL8
      case el_shell8:
      {
        INT arr_length;
        SHELL8* s8;
        DOUBLE* src_ptr;

        s8 = actele->e.s8;
        if (s8->nhyb)
        {

          arr_length = s8->alfa.fdim * s8->alfa.sdim;
          src_ptr = s8->alfa.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Dtildinv.fdim * s8->Dtildinv.sdim;
          src_ptr = s8->Dtildinv.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Lt.fdim * s8->Lt.sdim;
          src_ptr = s8->Lt.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s8->Rtilde.fdim * s8->Rtilde.sdim;
          src_ptr = s8->Rtilde.a.dv;
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }
        }
        if (mat[actele->mat-1].mattyp==m_viscohyper)
        {
          ARRAY4D *a;
          a = s8->his1;
          arr_length = a->fdim * a->sdim * a->tdim * a->fodim;
          src_ptr = a->a.d4[0][0][0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }
        }

        break;
      }
#endif

#ifdef D_SHELL9
      case el_shell9:
      {
        INT arr_length;
        SHELL9* s9;
        DOUBLE* src_ptr;

        s9 = actele->e.s9;
        if (s9->nhyb)
        {

          arr_length = s9->alfa.fdim * s9->alfa.sdim;
          src_ptr = s9->alfa.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Dtildinv.fdim * s9->Dtildinv.sdim;
          src_ptr = s9->Dtildinv.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->L.fdim * s9->L.sdim;
          src_ptr = s9->L.a.da[0];
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }

          arr_length = s9->Rtilde.fdim * s9->Rtilde.sdim;
          src_ptr = s9->Rtilde.a.dv;
          for (j=0; j<arr_length; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }
        }

        if (s9->elewa != NULL)
        {
          INT kl;

          for (kl=0; kl<s9->num_klay; kl++)
          {
            INT num_mlay;
            INT ml;
            INT actlay = 0;

            num_mlay = s9->kinlay[kl].num_mlay;
            for (ml=0; ml<num_mlay; ml++)
            {
              S9_IP_WA *ipwa = s9->elewa[actlay].ipwa;

              /* check if there is a ipwa for this layer */
              if (ipwa != NULL)
              {
                INT k;
                INT ngauss;
                MULTIMAT *actmultimat;

                actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

                /* number of gausspoints in one layer*/
                ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

                switch (actmultimat->mattyp)
                {
                case m_pl_mises:
                case m_pl_dp:
                  for (k=0; k<ngauss; k++)
                  {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].epstn;
                    for (j=0; j<6; j++)
                    {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                      *dst_ptr++ = ipwa[k].qn[j];
                    }
                  }
                  break;
                case m_pl_epc:
                  for (k=0; k<ngauss; k++)
                  {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].kappa_t;
                    *dst_ptr++ = ipwa[k].kappa_c;
                    for (j=0; j<6; j++)
                    {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                    }
                  }
                  break;
                case m_pl_hoff:
                  for (k=0; k<ngauss; k++)
                  {
                    *size_dst_ptr++ = ipwa[k].yip;
                    *dst_ptr++ = ipwa[k].dhard;
                    for (j=0; j<6; j++)
                    {
                      *dst_ptr++ = ipwa[k].sig[j];
                      *dst_ptr++ = ipwa[k].eps[j];
                      *dst_ptr++ = ipwa[k].dkappa[j];
                      *dst_ptr++ = ipwa[k].gamma[j];
                    }
                    for (j=0; j<9; j++)
                    {
                      *dst_ptr++ = ipwa[k].rkappa[j];
                    }
                  }
                  break;
                default:
                  printf("unsupported material %d for shell9 element",
                         actmultimat->mattyp);
                }
              }
              actlay += 1;
            }
          }
        }

        break;
      }
#endif

#ifdef D_WALL1
      case el_wall1:
      {
        W1_ELE_WA* elewa = actele->e.w1->elewa;
        if (elewa != NULL)
        {
          INT k;
          INT ngauss;

          dsassert((actele->distyp == quad4) ||
                   (actele->distyp == quad8) ||
                   (actele->distyp == quad9), "only quads supported up to now");

          ngauss = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];

          /* This is w1_write_restart. Disguised as usual. */
          switch (mat[actele->mat-1].mattyp)
          {
          case m_pl_mises:
            dsassert(len >= ngauss*(1+4*3), "value entry too short");
            for (k=0; k<ngauss; ++k)
            {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].epstn;
              for (j=0; j<4; j++)
              {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].qn[j];
              }
            }
            break;
          case m_pl_mises_3D:
            dsassert(len >= ngauss*(1+4*7), "value entry too short");
            for (k=0; k<ngauss; k++)
            {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].epstn;
              for (j=0; j<4; j++)
              {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].qn[j];
                *dst_ptr++ = elewa->ipwa[k].sigc[j];
                *dst_ptr++ = elewa->ipwa[k].sigi[j];
                *dst_ptr++ = elewa->ipwa[k].epsi[j];
                *dst_ptr++ = elewa->ipwa[k].di[j];
              }
            }
            break;
          case m_pl_epc3D: /* ???
                            * (mat[actele->mat-1].mattyp == m_pl_epc3D ) */
            dsassert(len >= ngauss*(2+4*7), "value entry too short");
            for (k=0; k<ngauss; k++)
            {
              *size_dst_ptr++ = elewa->ipwa[k].yip;
              *dst_ptr++ = elewa->ipwa[k].dlam[0];
              *dst_ptr++ = elewa->ipwa[k].dlam[1];
              for (j=0; j<4; j++)
              {
                *dst_ptr++ = elewa->ipwa[k].sig[j];
                *dst_ptr++ = elewa->ipwa[k].eps[j];
                *dst_ptr++ = elewa->ipwa[k].sigc[j];
                *dst_ptr++ = elewa->ipwa[k].grad[j];
                *dst_ptr++ = elewa->ipwa[k].sigi[j];
                *dst_ptr++ = elewa->ipwa[k].epsi[j];
                *dst_ptr++ = elewa->ipwa[k].di[  j];
              }
            }
            break;
          default:
            printf("WARNING: Unknown mattyp %d in wall1 restart\n",
                   mat[actele->mat-1].mattyp);
          }
        }
        break;
      }
#endif

#ifdef D_BEAM3
      case el_beam3:
        if (mat[actele->mat-1].mattyp == m_pl_mises ||
            mat[actele->mat-1].mattyp == m_pl_dp ||
            mat[actele->mat-1].mattyp == m_pl_epc)
        {
          INT j;
          INT size;
          DOUBLE* src_ptr;

          /* Copy the values from the element's array to the chunk. */
          src_ptr = actele->e.b3->elewa.a.da[0];
          size = actele->e.b3->elewa.fdim*actele->e.b3->elewa.sdim;
          dsassert(size <= len, "value entry too small");
          for (j=0; j<size; ++j)
          {
            *dst_ptr++ = *src_ptr++;
          }
        }
        break;
#endif

#ifdef D_INTERF
      case el_interf:
      {
        IF_ELE_WA* elewa = actele->e.interf->elewa;
        if (elewa != NULL)
        {
          INT k;
          INT ngauss;

          ngauss = actele->e.interf->nGP;

          dsassert(len >= ngauss*10, "value entry too short");
          for (k=0; k<ngauss; ++k)
          {
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Tt;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Tn;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].dt;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].dn;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].yip;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].jump_ut_pl;
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Q[0][0];
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Q[0][1];
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Q[1][0];
            *dst_ptr++ = actele->e.interf->elewa[0].ipwa[k].Q[1][1];
          }
        }
        break;
      }
#endif

#ifdef D_WALLGE
      /* There's nothing to be saved for wallge elements, right? */
      case el_wallge:
        break;
#endif

#ifdef D_FLUID2
      case el_fluid2:
      {
        FLUID_DYNAMIC *fdyn;
        DOUBLE* src_ptr;

        /* only valid for single field problem (?) */
        /* I don't think so. */
        fdyn = alldyn[genprob.numff].fdyn;

        dsassert(len >= 3, "value item too small");

        src_ptr = &(actele->e.f2->tau_old.a.dv[0]);
        for (j=0; j<3; ++j)
        {
          *dst_ptr++ = *src_ptr++;
        }
        if (fdyn->surftens > 0)
        {
          dsassert(len >= 3+2*actele->numnp, "value item too small");
          src_ptr = &(actele->e.f2->kappa_ND.a.dv[0]);
          if (actele->e.f2->fs_on>0)
          {
            for (j=0; j<2*actele->numnp; ++j)
            {
              *dst_ptr++ = *src_ptr++;
            }
          }
        }
        break;
      }
#endif

#ifdef D_FLUID3
      case el_fluid3:
      {
        DOUBLE* src_ptr;

        dsassert(len >= 3, "value item too small");
        src_ptr = &(actele->e.f3->tau_old.a.dv[0]);
        for (j=0; j<3; ++j)
        {
          *dst_ptr++ = *src_ptr++;
        }
        break;
      }
#endif

#ifdef D_FLUID3_F
      case el_fluid3_fast:
      {
        DOUBLE* src_ptr;

        dsassert(len >= 3, "value item too small");
        src_ptr = &(actele->e.f3->tau_old.a.dv[0]);
        for (j=0; j<3; ++j)
        {
          *dst_ptr++ = *src_ptr++;
        }
        break;
      }
#endif

#ifdef D_FLUID2_PRO

      /* Note that we do not reconstruct all pressure values on all
       * elements but only on those elements that are needed for
       * calculation. You have to calculate one time step before there
       * are again all pressures on each processor. */

      case el_fluid2_pro:
      {
        DOUBLE* src_ptr;

	switch (actele->e.f2pro->dm)
	{
	case dm_q2pm1:
	  dsassert(len >= 6, "value item too small");
	  src_ptr = actele->e.f2pro->press;
	  for (j=0; j<3; ++j)
	  {
	    *dst_ptr++ = *src_ptr++;
	  }
	  src_ptr = actele->e.f2pro->phi;
	  for (j=0; j<3; ++j)
	  {
	    *dst_ptr++ = *src_ptr++;
	  }
	  break;
	case dm_q1p0:
	  dsassert(len >= 2, "value item too small");

          src_ptr = actele->e.f2pro->press;
          *dst_ptr++ = *src_ptr++;

          src_ptr = actele->e.f2pro->phi;
          *dst_ptr++ = *src_ptr++;
	  break;
	default:
	  dserror("unsupported discretization mode %d", actele->e.f2pro->dm);
	}
        break;
      }
#endif

#ifdef D_FLUID3_PRO
      case el_fluid3_pro:
      {
        DOUBLE* src_ptr;

	switch (actele->e.f3pro->dm)
	{
	case dm_q2pm1:
	  dsassert(len >= 8, "value item too small");
	  src_ptr = actele->e.f3pro->press;
	  for (j=0; j<4; ++j)
	  {
	    *dst_ptr++ = *src_ptr++;
	  }
	  src_ptr = actele->e.f3pro->phi;
	  for (j=0; j<4; ++j)
	  {
	    *dst_ptr++ = *src_ptr++;
	  }
	  break;
	case dm_q1p0:
	  dsassert(len >= 2, "value item too small");

          src_ptr = actele->e.f3pro->press;
          *dst_ptr++ = *src_ptr++;

          src_ptr = actele->e.f3pro->phi;
          *dst_ptr++ = *src_ptr++;
	  break;
	default:
	  dserror("unsupported discretization mode %d", actele->e.f3pro->dm);
	}
        break;
      }
#endif

#ifdef D_BRICK1
      case el_brick1:
        /* Nothing to do?! */
        break;
#endif

      default:
      {
        static CHAR warning[el_count];
        if (!warning[actele->eltyp])
        {
          warning[actele->eltyp] = 1;
          printf(RED_LIGHT
                 "restart for element type '"
                 GREEN_LIGHT
                 "%d"
                 RED_LIGHT
                 "' not supported"
                 END_COLOR
                 "\n",
                 actele->eltyp);
        }
      }
      }

      counter += 1;
    }
  }

  dsassert(counter*len == send_count, "element pack count mismatch");

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack the design node this node is connected with (if any).

  \author u.kue
  \date 12/05
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dnode(BIN_OUT_CHUNK *chunk,
                           PARTDISCRET *actpdis,
                           DOUBLE *send_buf,
                           INT send_count,
                           INT *send_size_buf,
                           INT send_size_count,
                           INT dst_first_id,
                           INT dst_num)
{
  INT i;
  INT counter;
  INT slen;

#ifdef DEBUG
  dstrc_enter("out_pack_dline");
#endif

  slen = chunk->size_entry_length;

  dsassert(slen == 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  for (i=0; i<send_size_count; ++i)
  {
    send_size_buf[i] = -1;
  }

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      DNODE* dnode;
      INT* size_ptr;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actnode->Id_loc;
      size_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[slen*counter]);
#endif

      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        dnode = actnode->gnode->d.dnode;
        *size_ptr++ = dnode->Id;
        break;
      case ondline:
      case ondsurf:
      case ondvol:
        /* no dnode on this node */
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      counter += 1;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Pack all design lines this node is connected with.

  That is a little troublesome because we need to take indirect
  connections into account.

  \author u.kue
  \date 12/05
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dline(BIN_OUT_CHUNK *chunk,
                           PARTDISCRET *actpdis,
                           DOUBLE *send_buf,
                           INT send_count,
                           INT *send_size_buf,
                           INT send_size_count,
                           INT dst_first_id,
                           INT dst_num)
{
  INT i;
  INT counter;
  INT slen;
  INTSET set;

#ifdef DEBUG
  dstrc_enter("out_pack_dline");
#endif

  slen = chunk->size_entry_length;

  dsassert(slen >= 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  for (i=0; i<send_size_count; ++i)
  {
    send_size_buf[i] = -1;
  }

  intset_init(&set, 10);

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT j=0;
      INT ndnode=0;
      DNODE* dnode=NULL;
      DLINE** pdline=NULL;
      INT ndline=0;
      DLINE* dline;

      INT* size_ptr;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actnode->Id_loc;
      size_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[slen*counter]);
#endif

      intset_clear(&set);

      /* determine starting point */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        dnode = actnode->gnode->d.dnode;
        ndnode = 1;
        break;
      case ondline:
        pdline = &(actnode->gnode->d.dline);
        ndline = 1;
        break;
      case ondsurf:
        /* Nothing to do :) */
        break;
      case ondvol:
        /* Nothing to do :) */
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      /* count dnode, dline, dsurf, dvol to one gnode */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        ndline = dnode->ndline;
        pdline = dnode->dline;
      case ondline:
        for (j=0; j<ndline; ++j)
        {
          dline = pdline[j];
	  if (!intset_contains(&set, dline->Id))
	  {
	    intset_add(&set, dline->Id);
	  }
        }
        break;
      case ondsurf:
        /* Nothing to do :) */
        break;
      case ondvol:
        /* Nothing to do :) */
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      for (j=0; j<set.count; ++j)
      {
        *size_ptr++ = set.value[j];
      }

      counter += 1;
    }
  }

  intset_destroy(&set);

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief Pack all design surfaces this node is connected with.

  That is a little troublesome because we need to take indirect
  connections into account.

  \author u.kue
  \date 12/05
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dsurf(BIN_OUT_CHUNK *chunk,
                           PARTDISCRET *actpdis,
                           DOUBLE *send_buf,
                           INT send_count,
                           INT *send_size_buf,
                           INT send_size_count,
                           INT dst_first_id,
                           INT dst_num)
{
  INT i;
  INT counter;
  INT slen;
  INTSET set;

#ifdef DEBUG
  dstrc_enter("out_pack_dsurf");
#endif

  slen = chunk->size_entry_length;

  dsassert(slen >= 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  for (i=0; i<send_size_count; ++i)
  {
    send_size_buf[i] = -1;
  }

  intset_init(&set, 10);

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT j=0;
      INT ndnode=0;
      DNODE* dnode=NULL;
      DLINE** pdline=NULL;
      INT k=0;
      INT ndline=0;
      DLINE* dline;
      DSURF** pdsurf=NULL;
      INT ndsurf=0;
      DSURF* dsurf;

      INT* size_ptr;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actnode->Id_loc;
      size_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[slen*counter]);
#endif

      intset_clear(&set);

      /* determine starting point */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        dnode = actnode->gnode->d.dnode;
        ndnode = 1;
        break;
      case ondline:
        pdline = &(actnode->gnode->d.dline);
        ndline = 1;
        break;
      case ondsurf:
        pdsurf = &(actnode->gnode->d.dsurf);
        ndsurf = 1;
        break;
      case ondvol:
        /* Nothing to do :) */
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      /* count dnode, dline, dsurf, dvol to one gnode */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        ndline = dnode->ndline;
        pdline = dnode->dline;
      case ondline:
        for (j=0; j<ndline; ++j)
        {
          dline = pdline[j];
          ndsurf = dline->ndsurf;
          pdsurf = dline->dsurf;
        case ondsurf:
          for (k=0; k<ndsurf; ++k)
          {
            dsurf = pdsurf[k];
	    if (!intset_contains(&set, dsurf->Id))
	    {
	      intset_add(&set, dsurf->Id);
	    }
          }
        }
        break;
      case ondvol:
        /* Nothing to do :) */
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      for (j=0; j<set.count; ++j)
      {
        *size_ptr++ = set.value[j];
      }

      counter += 1;
    }
  }

  intset_destroy(&set);

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief Pack all design volumes this node is connected with.

  That is a little troublesome because we need to take indirect
  connections into account.

  \author u.kue
  \date 12/05
  \sa out_pack_items
*/
/*----------------------------------------------------------------------*/
static void out_pack_dvol(BIN_OUT_CHUNK *chunk,
                          PARTDISCRET *actpdis,
                          DOUBLE *send_buf,
                          INT send_count,
                          INT *send_size_buf,
                          INT send_size_count,
                          INT dst_first_id,
                          INT dst_num)
{
  INT i;
  INT counter;
  INT slen;
  INTSET set;

#ifdef DEBUG
  dstrc_enter("out_pack_dvol");
#endif

  slen = chunk->size_entry_length;

  dsassert(slen >= 1, "invalid size entry length");
  dsassert(chunk->value_entry_length == 0, "invalid value entry length");

  for (i=0; i<send_size_count; ++i)
  {
    send_size_buf[i] = -1;
  }

  intset_init(&set, 10);

  counter = 0;
  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode = actpdis->node[i];
    if ((actnode->Id_loc >= dst_first_id) &&
        (actnode->Id_loc < dst_first_id+dst_num))
    {
      INT j=0;
      INT ndnode=0;
      DNODE* dnode=NULL;
      DLINE** pdline=NULL;
      INT k=0;
      INT ndline=0;
      DLINE* dline;
      DSURF** pdsurf=NULL;
      INT ndsurf=0;
      DSURF* dsurf;
      INT l=0;
      INT ndvol=0;
      DVOL* dvol;
      DVOL** pdvol=NULL;

      INT* size_ptr;

#ifdef PARALLEL
      send_size_buf[(slen+1)*counter] = actnode->Id_loc;
      size_ptr = &(send_size_buf[(slen+1)*counter+1]);
#else
      size_ptr = &(send_size_buf[slen*counter]);
#endif

      intset_clear(&set);

      /* determine starting point */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        dnode = actnode->gnode->d.dnode;
        ndnode = 1;
        break;
      case ondline:
        pdline = &(actnode->gnode->d.dline);
        ndline = 1;
        break;
      case ondsurf:
        pdsurf = &(actnode->gnode->d.dsurf);
        ndsurf = 1;
        break;
      case ondvol:
        pdvol = &(actnode->gnode->d.dvol);
        ndvol = 1;
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      /* count dnode, dline, dsurf, dvol to one gnode */
      switch (actnode->gnode->ondesigntyp)
      {
      case ondnode:
        ndline = dnode->ndline;
        pdline = dnode->dline;
      case ondline:
        for (j=0; j<ndline; ++j)
        {
          dline = pdline[j];
          ndsurf = dline->ndsurf;
          pdsurf = dline->dsurf;
        case ondsurf:
          for (k=0; k<ndsurf; ++k)
          {
            dsurf = pdsurf[k];
            ndvol = dsurf->ndvol;
            pdvol = dsurf->dvol;
          case ondvol:
            for (l=0; l<ndvol; ++l)
            {
              dvol = pdvol[l];
              if (!intset_contains(&set, dvol->Id))
              {
                intset_add(&set, dvol->Id);
              }
            }
          }
        }
        break;
      default:
        dserror("invalid design type %d", actnode->gnode->ondesigntyp);
      }

      for (j=0; j<set.count; ++j)
      {
        *size_ptr++ = set.value[j];
      }

      counter += 1;
    }
  }

  intset_destroy(&set);

#ifdef DEBUG
  dstrc_exit();
#endif
}



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
  \param  dst_num        (i)  the number of (consecutive) items to be
  written on the receiving processor in total

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
                    INT send_size_count,
                    INT dst_first_id,
                    INT dst_num)
{
#ifdef DEBUG
  dstrc_enter("out_pack_items");
#endif

  switch (type)
  {
  case cc_node_array:
    out_pack_node_arrays(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_mesh:
    /* array flag unused here */
    out_pack_mesh(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_coords:
    /* array flag unused here */
    out_pack_coords(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_ele_params:
    /* array flag unused here */
    out_pack_ele_params(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_displacement:
    out_pack_displacement(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_velocity:
    out_pack_velocity(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_pressure:
    out_pack_pressure(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_av_pressure:
    out_pack_average_pressure(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_stress:
    out_pack_stress(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_domain:
    /* array flag unused here */
    out_pack_domain(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_dist_vector:
    /* array flag unused here */
    out_pack_dist_vector(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#ifdef D_SHELL8
  case cc_shell8_director:
    /* array flag unused here */
    out_pack_shell8_director(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#endif
#ifdef D_SHELL9
  case cc_shell9_coords:
    /* array flag unused here */
    out_pack_shell9_coords(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_shell9_displacement:
    out_pack_shell9_displacement(chunk, array, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
#endif
  case cc_restart_element:
    /* array flag unused here */
    out_pack_restart_element(chunk, actpdis, send_buf, send_count, send_size_buf, dst_first_id, dst_num);
    break;
  case cc_dnode:
    /* array flag unused here */
    out_pack_dnode(chunk, actpdis, send_buf, send_count, send_size_buf, send_size_count, dst_first_id, dst_num);
    break;
  case cc_dline:
    /* array flag unused here */
    out_pack_dline(chunk, actpdis, send_buf, send_count, send_size_buf, send_size_count, dst_first_id, dst_num);
    break;
  case cc_dsurf:
    /* array flag unused here */
    out_pack_dsurf(chunk, actpdis, send_buf, send_count, send_size_buf, send_size_count, dst_first_id, dst_num);
    break;
  case cc_dvol:
    /* array flag unused here */
    out_pack_dvol(chunk, actpdis, send_buf, send_count, send_size_buf, send_size_count, dst_first_id, dst_num);
    break;
  default:
    dserror("unsupported chunk type %d", type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*======================================================================*/
/* Here we have lots of functions that do nothing but unpacking input
 * buffers and put the values somewhere (node arrays, element stuff,
 * distributed vectors). These functions share a common structure. See
 * ``in_unpack_items`` for an explanation of their parameters.
 *
 * These are the functions that need to be changed according to
 * changes in the elements. */
/*======================================================================*/


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack the node arrays.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_node_arrays(BIN_IN_FIELD *context,
                                  BIN_IN_CHUNK *chunk,
                                  NODE_ARRAY array,
                                  DOUBLE *recv_buf,
                                  INT recv_count,
                                  INT *recv_size_buf,
                                  INT recv_size_count,
                                  INT src)
{
  INT j;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("in_unpack_node_arrays");
#endif

  /* There are only very small differences between the sequential
   * and the parallel version here. */

#ifdef PARALLEL
#define unpack_recv_size_buf                                            \
  dsassert(3*j+3 <= recv_size_count,                                    \
           "receive size buffer overrun");                              \
  dsassert(recv_size_buf[3*j] == actnode->Id_loc, "Id_loc mismatch");   \
  fdim = recv_size_buf[3*j+1];                                          \
  sdim = recv_size_buf[3*j+2];
#else
#define unpack_recv_size_buf                    \
  dsassert(2*j+2 <= recv_size_count,            \
           "receive size buffer overrun");      \
  fdim = recv_size_buf[2*j  ];                  \
  sdim = recv_size_buf[2*j+1];
#endif

#ifdef PARALLEL
#define get_id_part(j) (context->recv_node_ids[src][2*j+1])
#define get_recv_numnp (context->recv_numnp[src])
#else
#define get_id_part(j) (j)
#define get_recv_numnp (actpdis->numnp)
#endif

  dsassert(chunk->size_entry_length == 2, "size entry length mismatch");

  /* scatter received data to the nodes */
#define call_scatter_values(node_array)                                 \
  for (j=0; j<get_recv_numnp; ++j)                                      \
  {                                                                     \
    INT k;                                                              \
    INT Id_part;                                                        \
    NODE* actnode;                                                      \
    INT fdim;                                                           \
    INT sdim;                                                           \
    DOUBLE *ptr;                                                        \
    Id_part = get_id_part(j);                                           \
    dsassert((Id_part >= 0) && (Id_part < actpdis->numnp),              \
             "illegal node id");                                        \
    actnode = actpdis->node[Id_part];                                   \
    /* Id_part only in the more parallel version. */                    \
    /*dsassert(Id_part == actnode->Id_part, "Id_part mismatch");*/      \
                                                                        \
    unpack_recv_size_buf;                                               \
                                                                        \
    dsassert(sdim == actnode->node_array.sdim,                          \
             "node array sdim mismatch");                               \
    if (fdim > actnode->node_array.fdim)                                \
    {                                                                   \
      amredef(&(actnode->node_array),fdim,sdim,"DA");                   \
    }                                                                   \
    ptr = actnode->node_array.a.da[0];                                  \
    dsassert(chunk->value_entry_length*j+fdim*sdim <= recv_count,       \
             "receive buffer overrun");                                 \
    for (k=0; k<fdim*sdim; ++k)                                         \
    {                                                                   \
      *ptr++ = recv_buf[chunk->value_entry_length*j+k];                 \
    }                                                                   \
  }

  switch (array)
  {
  case node_array_sol:
    call_scatter_values(sol);
    break;
  case node_array_sol_increment:
    call_scatter_values(sol_increment);
    break;
  case node_array_sol_residual:
    call_scatter_values(sol_residual);
    break;
  case node_array_sol_mf:
    call_scatter_values(sol_mf);
    break;
  default:
    dserror("Node array %d unknown", array);
  }

#undef call_scatter_values
#undef get_id_part
#undef get_recv_numnp
#undef unpack_recv_size_buf

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack an array of distributed vectors.

  This is not solver dependent because we already know the locations
  where the dof values are to got to. The field initialization found
  this (solver dependent) information.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_dist_vector(BIN_IN_FIELD *context,
                                  BIN_IN_CHUNK *chunk,
                                  DOUBLE *recv_buf,
                                  INT recv_count,
                                  INT *recv_size_buf,
                                  INT recv_size_count,
                                  INT src)
{
  INT j;
  INT len;

#ifdef DEBUG
  dstrc_enter("in_unpack_dist_vector");
#endif

#ifdef PARALLEL
#define get_id_part(j)  (context->recv_dof_ids[src][2*j+1])
#define get_recv_numdof (context->recv_numdof[src])
#else
#define get_id_part(j)  (j)
#define get_recv_numdof (chunk->vectors[0].numeq)
#endif

  dsassert(chunk->size_entry_length == 0, "size entry length mismatch");

  len = chunk->value_entry_length;

  /* scatter received data to the nodes */
  for (j=0; j<get_recv_numdof; ++j)
  {
    INT k;
    INT Id_part;
    DOUBLE* src_ptr;

    Id_part = get_id_part(j);

    src_ptr = &(recv_buf[len*j]);
    for (k=0; k<len; ++k)
    {
      chunk->vectors[k].vec.a.dv[Id_part] = *src_ptr++;
    }
  }

#undef get_id_part
#undef get_recv_numele

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Unpack what we received.

  This is called during restart, the elements are already set up. The
  arrays are there. Just the values are missing.

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \author u.kue
  \date 09/04
  \sa in_unpack_items
*/
/*----------------------------------------------------------------------*/
static void in_unpack_restart_element(BIN_IN_FIELD *context,
                                      BIN_IN_CHUNK *chunk,
                                      DOUBLE *recv_buf,
                                      INT recv_count,
                                      INT *recv_size_buf,
                                      INT recv_size_count,
                                      INT src)
{
  INT el;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);
  INT len;
  INT slen;

#ifdef DEBUG
  dstrc_enter("in_unpack_restart_element");
#endif

#ifdef PARALLEL
#define get_id_part(j)  (context->recv_element_ids[src][2*j+1])
#define get_recv_numele (context->recv_numele[src])
#else
#define get_id_part(j)  (j)
#define get_recv_numele (actpdis->numele)
#endif

  len = chunk->value_entry_length;
  slen = chunk->size_entry_length;

  /* scatter received data to the nodes */
  for (el=0; el<get_recv_numele; ++el)
  {
    INT Id_part;
    ELEMENT* actele;
    DOUBLE* src_ptr;
    INT* size_src_ptr;

    Id_part = get_id_part(el);
    dsassert((Id_part >= 0) && (Id_part < actpdis->numele),
             "illegal element id");
    actele = actpdis->element[Id_part];
    /*dsassert(Id_part == actele->Id_part, "Id_part mismatch");*/

#ifdef PARALLEL
    dsassert(recv_size_buf[(slen+1)*el] == actele->Id_loc, "Id_loc mismatch");
    size_src_ptr = &(recv_size_buf[(slen+1)*el+1]);
#else
    size_src_ptr = &(recv_size_buf[slen*el]);
#endif

    src_ptr = &(recv_buf[len*el]);

    switch (actele->eltyp)
    {

#ifdef D_SHELL8
    case el_shell8:
    {
      INT el;
      INT arr_length;
      SHELL8* s8;
      DOUBLE* dst_ptr;

      s8 = actele->e.s8;
      if (s8->nhyb)
      {

        arr_length = s8->alfa.fdim * s8->alfa.sdim;
        dst_ptr = s8->alfa.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s8->Dtildinv.fdim * s8->Dtildinv.sdim;
        dst_ptr = s8->Dtildinv.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s8->Lt.fdim * s8->Lt.sdim;
        dst_ptr = s8->Lt.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s8->Rtilde.fdim * s8->Rtilde.sdim;
        dst_ptr = s8->Rtilde.a.dv;
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }
      }
      if (mat[actele->mat-1].mattyp==m_viscohyper)
      {
        ARRAY4D *a;
        a = s8->his1;
        arr_length = a->fdim * a->sdim * a->tdim * a->fodim;
        dst_ptr = a->a.d4[0][0][0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }
      }

      break;
    }
#endif

#ifdef D_SHELL9
    case el_shell9:
    {
      INT el;
      INT arr_length;
      SHELL9* s9;
      DOUBLE* dst_ptr;

      s9 = actele->e.s9;
      if (s9->nhyb)
      {

        arr_length = s9->alfa.fdim * s9->alfa.sdim;
        dst_ptr = s9->alfa.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s9->Dtildinv.fdim * s9->Dtildinv.sdim;
        dst_ptr = s9->Dtildinv.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s9->L.fdim * s9->L.sdim;
        dst_ptr = s9->L.a.da[0];
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }

        arr_length = s9->Rtilde.fdim * s9->Rtilde.sdim;
        dst_ptr = s9->Rtilde.a.dv;
        for (el=0; el<arr_length; ++el)
        {
          *dst_ptr++ = *src_ptr++;
        }
      }

      if (s9->elewa != NULL)
      {
        INT kl;

        for (kl=0; kl<s9->num_klay; kl++)
        {
          INT num_mlay;
          INT ml;
          INT actlay = 0;

          num_mlay = s9->kinlay[kl].num_mlay;
          for (ml=0; ml<num_mlay; ml++)
          {
            S9_IP_WA *ipwa = s9->elewa[actlay].ipwa;

            /* check if there is a ipwa for this layer */
            if (ipwa != NULL)
            {
              INT k;
              INT ngauss;
              MULTIMAT *actmultimat;

              actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

              /* number of gausspoints in one layer*/
              ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

              switch (actmultimat->mattyp)
              {
              case m_pl_mises:
              case m_pl_dp:
                for (k=0; k<ngauss; k++)
                {
                  INT el;
                  ipwa[k].yip = *size_src_ptr++;
                  ipwa[k].epstn = *src_ptr++;
                  for (el=0; el<6; el++)
                  {
                    ipwa[k].sig[el] = *src_ptr++;
                    ipwa[k].eps[el] = *src_ptr++;
                    ipwa[k].qn[el] = *src_ptr++;
                  }
                }
                break;
              case m_pl_epc:
                for (k=0; k<ngauss; k++)
                {
                  INT el;
                  ipwa[k].yip = *size_src_ptr++;
                  ipwa[k].kappa_t = *src_ptr++;
                  ipwa[k].kappa_c = *src_ptr++;
                  for (el=0; el<6; el++)
                  {
                    ipwa[k].sig[el] = *src_ptr++;
                    ipwa[k].eps[el] = *src_ptr++;
                  }
                }
                break;
              case m_pl_hoff:
                for (k=0; k<ngauss; k++)
                {
                  INT el;
                  ipwa[k].yip = *size_src_ptr++;
                  ipwa[k].dhard = *src_ptr++;
                  for (el=0; el<6; el++)
                  {
                    ipwa[k].sig[el] = *src_ptr++;
                    ipwa[k].eps[el] = *src_ptr++;
                    ipwa[k].dkappa[el] = *src_ptr++;
                    ipwa[k].gamma[el] = *src_ptr++;
                  }
                  for (el=0; el<9; el++)
                  {
                    ipwa[k].rkappa[el] = *src_ptr++;
                  }
                }
                break;
              default:
                printf("unsupported material %d for shell9 element",
                       actmultimat->mattyp);
              }
            }
            actlay += 1;
          }
        }
      }
      break;
    }
#endif

#ifdef D_WALL1
    case el_wall1:
    {
      W1_ELE_WA* elewa = actele->e.w1->elewa;
      if (elewa != NULL)
      {
        INT k;
        INT ngauss;

        dsassert((actele->distyp == quad4) ||
                 (actele->distyp == quad8) ||
                 (actele->distyp == quad9), "only quads supported up to now");

        ngauss = actele->e.w1->nGP[0]*actele->e.w1->nGP[1];

        /* This is w1_write_restart. Disguised as usual. */
        switch (mat[actele->mat-1].mattyp)
        {
        case m_pl_mises:
          dsassert(len >= ngauss*(1+4*3), "value entry too short");
          for (k=0; k<ngauss; ++k)
          {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].epstn = *src_ptr++;
            for (el=0; el<4; el++)
            {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].qn[el] = *src_ptr++;
            }
          }
          break;
        case m_pl_mises_3D:
          dsassert(len >= ngauss*(1+4*7), "value entry too short");
          for (k=0; k<ngauss; k++)
          {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].epstn = *src_ptr++;
            for (el=0; el<4; el++)
            {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].qn[el] = *src_ptr++;
              elewa->ipwa[k].sigc[el] = *src_ptr++;
              elewa->ipwa[k].sigi[el] = *src_ptr++;
              elewa->ipwa[k].epsi[el] = *src_ptr++;
              elewa->ipwa[k].di[el] = *src_ptr++;
            }
          }
          break;
        case m_pl_epc3D: /* ???
                          * (mat[actele->mat-1].mattyp == m_pl_epc3D ) */
          dsassert(len >= ngauss*(2+4*7), "value entry too short");
          for (k=0; k<ngauss; k++)
          {
            elewa->ipwa[k].yip = *size_src_ptr++;
            elewa->ipwa[k].dlam[0] = *src_ptr++;
            elewa->ipwa[k].dlam[1] = *src_ptr++;
            for (el=0; el<4; el++)
            {
              elewa->ipwa[k].sig[el] = *src_ptr++;
              elewa->ipwa[k].eps[el] = *src_ptr++;
              elewa->ipwa[k].sigc[el] = *src_ptr++;
              elewa->ipwa[k].grad[el] = *src_ptr++;
              elewa->ipwa[k].sigi[el] = *src_ptr++;
              elewa->ipwa[k].epsi[el] = *src_ptr++;
              elewa->ipwa[k].di[  el] = *src_ptr++;
            }
          }
          break;
        default:
          printf("WARNING: Unknown mattyp %d in wall1 restart\n",
                 mat[actele->mat-1].mattyp);
        }
      }
      break;
    }
#endif

#ifdef D_BEAM3
    case el_beam3:
      if (mat[actele->mat-1].mattyp == m_pl_mises ||
          mat[actele->mat-1].mattyp == m_pl_dp ||
          mat[actele->mat-1].mattyp == m_pl_epc)
      {
        INT j;
        INT size;
        DOUBLE* dst_ptr;

        /* Copy the values from the chunk to the element's array. */
        dst_ptr = actele->e.b3->elewa.a.da[0];
        size = actele->e.b3->elewa.fdim*actele->e.b3->elewa.sdim;
        dsassert(size <= len, "value entry too small");
        for (j=0; j<size; ++j)
        {
          *dst_ptr++ = *src_ptr++;
        }
      }
      break;
#endif

#ifdef D_INTERF
    case el_interf:
    {
      IF_ELE_WA* elewa = actele->e.interf->elewa;
      if (elewa != NULL)
      {
        INT k;
        INT ngauss;

        ngauss = actele->e.interf->nGP;

        dsassert(len >= ngauss*10, "value entry too short");
        for (k=0; k<ngauss; ++k)
        {
          actele->e.interf->elewa[0].ipwa[k].Tt = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].Tn = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].dt = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].dn = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].yip = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].jump_ut_pl = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].Q[0][0] = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].Q[0][1] = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].Q[1][0] = *src_ptr++;
          actele->e.interf->elewa[0].ipwa[k].Q[1][1] = *src_ptr++;
        }
      }
      break;
    }
#endif

#ifdef D_WALLGE
    /* There's nothing to be restored for wallge elements, right? */
    case el_wallge:
      break;
#endif

#ifdef D_FLUID2
    case el_fluid2:
    {
      FLUID_DYNAMIC *fdyn;
      DOUBLE* dst_ptr;
      INT k;

      /* only valid for single field problem (?) */
      /* I don't think so. */
      fdyn = alldyn[genprob.numff].fdyn;

      dsassert(len >= 3, "value item too small");
      dsassert(actele->e.f2->tau_old.Typ == cca_DV &&
               actele->e.f2->tau_old.fdim == 3, "uninitialized array");
      dst_ptr = &(actele->e.f2->tau_old.a.dv[0]);
      for (k=0; k<3; ++k)
      {
        *dst_ptr++ = *src_ptr++;
      }
      if (fdyn->surftens > 0)
      {
        dsassert(len >= 3+2*actele->numnp, "value item too small");
        dsassert(actele->e.f2->kappa_ND.Typ == cca_DA, "uninitialized array");
        dst_ptr = &(actele->e.f2->kappa_ND.a.dv[0]);
        if (actele->e.f2->fs_on>0)
        {
          for (k=0; k<2*actele->numnp; ++k)
          {
            *dst_ptr++ = *src_ptr++;
          }
        }
      }
      break;
    }
#endif

#ifdef D_FLUID3
    case el_fluid3:
    {
      DOUBLE* dst_ptr;
      INT k;

      dsassert(len >= 3, "value item too small");
      dsassert(actele->e.f3->tau_old.Typ == cca_DV &&
               actele->e.f3->tau_old.fdim == 3, "uninitialized array");
      dst_ptr = &(actele->e.f3->tau_old.a.dv[0]);
      for (k=0; k<3; ++k)
      {
        *dst_ptr++ = *src_ptr++;
      }
      break;
    }
#endif

#ifdef D_FLUID3_F
    case el_fluid3_fast:
    {
      DOUBLE* dst_ptr;
      INT k;

      dsassert(len >= 3, "value item too small");
      dsassert(actele->e.f3->tau_old.Typ == cca_DV &&
               actele->e.f3->tau_old.fdim == 3, "uninitialized array");
      dst_ptr = &(actele->e.f3->tau_old.a.dv[0]);
      for (k=0; k<3; ++k)
      {
        *dst_ptr++ = *src_ptr++;
      }
      break;
    }
#endif

#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      DOUBLE* dst_ptr;
      INT k;

      switch (actele->e.f2pro->dm)
      {
      case dm_q2pm1:
	dsassert(len >= 6, "value item too small");
	dst_ptr = actele->e.f2pro->press;
	for (k=0; k<3; ++k)
	{
	  *dst_ptr++ = *src_ptr++;
	}
	src_ptr = actele->e.f2pro->phi;
	for (k=0; k<3; ++k)
	{
	  *dst_ptr++ = *src_ptr++;
	}
	break;
      case dm_q1p0:
        dsassert(len >= 2, "value item too small");

	dst_ptr = actele->e.f2pro->press;
        *dst_ptr++ = *src_ptr++;

        src_ptr = actele->e.f2pro->phi;
        *dst_ptr++ = *src_ptr++;
        break;
      default:
	dserror("unsupported discretization mode %d", actele->e.f2pro->dm);
      }
      break;
    }
#endif

#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      DOUBLE* dst_ptr;
      INT k;

      switch (actele->e.f3pro->dm)
      {
      case dm_q2pm1:
	dsassert(len >= 8, "value item too small");
	dst_ptr = actele->e.f3pro->press;
	for (k=0; k<4; ++k)
	{
	  *dst_ptr++ = *src_ptr++;
	}
	src_ptr = actele->e.f3pro->phi;
	for (k=0; k<4; ++k)
	{
	  *dst_ptr++ = *src_ptr++;
	}
	break;
      case dm_q1p0:
        dsassert(len >= 2, "value item too small");

	dst_ptr = actele->e.f3pro->press;
        *dst_ptr++ = *src_ptr++;

        src_ptr = actele->e.f3pro->phi;
        *dst_ptr++ = *src_ptr++;
        break;
      default:
	dserror("unsupported discretization mode %d", actele->e.f3pro->dm);
      }
      break;
    }
#endif

#ifdef D_BRICK1
    case el_brick1:
      /* Nothing to do?! */
      break;
#endif

    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

#undef get_id_part
#undef get_recv_numele

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                     INT src)
{
  INT value_length;
  INT size_length;

#ifdef DEBUG
  dstrc_enter("in_unpack_items");
#endif

  value_length = chunk->value_entry_length;
  size_length  = chunk->size_entry_length;

  switch (type)
  {
  case cc_node_array:
    in_unpack_node_arrays(context, chunk, array, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  case cc_dist_vector:
    in_unpack_dist_vector(context, chunk, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  case cc_restart_element:
    in_unpack_restart_element(context, chunk, recv_buf, recv_count, recv_size_buf, recv_size_count, src);
    break;
  default:
    dserror("unsupported chunk type %d", type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*======================================================================*/
/* After packing and unpacking there is a third element dependent
 * task. To tell the required sizes. Most of the time these functions
 * are needed to provide the space to pack the elements for
 * writing. When we read elements we the control file tells the
 * required sizes. */
/*======================================================================*/


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the design object ids for each node.

  \author u.kue
  \date 12/05
  \sa out_pack_ele_params
*/
/*----------------------------------------------------------------------*/
void find_design_item_length(struct _BIN_OUT_FIELD* context,
                             INT* dnodemax,
                             INT* dlinemax,
                             INT* dsurfmax,
                             INT* dvolmax)
{
  INT i;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

  INT dobjectmax[4] = { 0, 0, 0, 0 };
  INTSET dline_set;
  INTSET dsurf_set;
  INTSET dvol_set;

#ifdef PARALLEL
  INT recv_length[4];
#endif

#ifdef DEBUG
  dstrc_enter("find_design_item_length");
#endif

  intset_init(&dline_set,10);
  intset_init(&dsurf_set,10);
  intset_init(&dvol_set,10);

  for (i=0; i<actpdis->numnp; ++i)
  {
    INT dnodecount = 0;

    INT j=0;
    INT ndnode=0;
    DNODE* dnode=NULL;
    DLINE** pdline=NULL;
    INT k=0;
    INT ndline=0;
    DLINE* dline;
    DSURF** pdsurf=NULL;
    INT ndsurf=0;
    DSURF* dsurf;
    INT l=0;
    INT ndvol=0;
    DVOL* dvol;
    DVOL** pdvol=NULL;

    NODE* actnode;
    actnode = actpdis->node[i];

    intset_clear(&dline_set);
    intset_clear(&dsurf_set);
    intset_clear(&dvol_set);

    /* determine starting point */
    switch (actnode->gnode->ondesigntyp)
    {
    case ondnode:
      dnode = actnode->gnode->d.dnode;
      ndnode = 1;
      break;
    case ondline:
      pdline = &(actnode->gnode->d.dline);
      ndline = 1;
      break;
    case ondsurf:
      pdsurf = &(actnode->gnode->d.dsurf);
      ndsurf = 1;
      break;
    case ondvol:
      pdvol = &(actnode->gnode->d.dvol);
      ndvol = 1;
      break;
    default:
      dserror("invalid design type %d", actnode->gnode->ondesigntyp);
    }

    /* count dnode, dline, dsurf, dvol to one gnode */
    /* note the missing breaks! this switch is nothing but a
     * structured goto construct. */
    switch (actnode->gnode->ondesigntyp)
    {
    case ondnode:
      dnodecount += ndnode;
      ndline = dnode->ndline;
      pdline = dnode->dline;
    case ondline:
      for (j=0; j<ndline; ++j)
      {
        dline = pdline[j];
        ndsurf = dline->ndsurf;
        pdsurf = dline->dsurf;
	if (!intset_contains(&dline_set, dline->Id))
	{
	  intset_add(&dline_set, dline->Id);
	}
      case ondsurf:
        for (k=0; k<ndsurf; ++k)
        {
          dsurf = pdsurf[k];
          ndvol = dsurf->ndvol;
	  pdvol = dsurf->dvol;
	  if (!intset_contains(&dsurf_set, dsurf->Id))
	  {
	    intset_add(&dsurf_set, dsurf->Id);
	  }
        case ondvol:
	  for (l=0; l<ndvol; ++l)
	  {
	    dvol = pdvol[l];
	    if (!intset_contains(&dvol_set, dvol->Id))
	    {
	      intset_add(&dvol_set, dvol->Id);
	    }
	  }
        }
      }
      break;
    default:
      dserror("invalid design type %d", actnode->gnode->ondesigntyp);
    }

    dobjectmax[0] = MAX(dobjectmax[0], dnodecount);
    dobjectmax[1] = MAX(dobjectmax[1], dline_set.count);
    dobjectmax[2] = MAX(dobjectmax[2], dsurf_set.count);
    dobjectmax[3] = MAX(dobjectmax[3], dvol_set.count);
  }

#ifdef PARALLEL
  MPI_Allreduce(dobjectmax, recv_length, 4, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
  *dnodemax = recv_length[0];
  *dlinemax = recv_length[1];
  *dsurfmax = recv_length[2];
  *dvolmax  = recv_length[3];
#else
  *dnodemax = dobjectmax[0];
  *dlinemax = dobjectmax[1];
  *dsurfmax = dobjectmax[2];
  *dvolmax  = dobjectmax[3];
#endif

  intset_destroy(&dline_set);
  intset_destroy(&dsurf_set);
  intset_destroy(&dvol_set);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element parameters.

  \author u.kue
  \date 11/04
  \sa out_pack_ele_params
*/
/*----------------------------------------------------------------------*/
void find_ele_param_item_length(struct _BIN_OUT_FIELD* context,
                                INT* value_length,
                                INT* size_length)
{
  INT vlen = 0;
  INT slen = 0;
  INT i;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("find_ele_param_item_length");
#endif

  /* element specific parameter length */

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele = actpdis->element[i];

    switch (actele->eltyp)
    {
#ifdef D_SHELL8
    case el_shell8:
      slen = MAX(slen, shell8_variables.ep_size_length);
      vlen = MAX(vlen, shell8_variables.ep_value_length);
      break;
#endif
#ifdef D_SHELL9
    case el_shell9:
      slen = MAX(slen, shell9_variables.ep_size_length);
      vlen = MAX(vlen, shell9_variables.ep_value_length);
      break;
#endif
#ifdef D_BRICK1
    case el_brick1:
      slen = MAX(slen, brick1_variables.ep_size_length);
      vlen = MAX(vlen, brick1_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID2
    case el_fluid2:
      slen = MAX(slen, fluid2_variables.ep_size_length);
      vlen = MAX(vlen, fluid2_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      slen = MAX(slen, fluid2_pro_variables.ep_size_length);
      vlen = MAX(vlen, fluid2_pro_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
      slen = MAX(slen, fluid3_pro_variables.ep_size_length);
      vlen = MAX(vlen, fluid3_pro_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID2TU
    case el_fluid2_tu:
      slen = MAX(slen, fluid2tu_variables.ep_size_length);
      vlen = MAX(vlen, fluid2tu_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID3
    case el_fluid3:
      slen = MAX(slen, fluid3_variables.ep_size_length);
      vlen = MAX(vlen, fluid3_variables.ep_value_length);
      break;
#endif
#ifdef D_FLUID3_F
    case el_fluid3_fast:
      slen = MAX(slen, fluid3_fast_variables.ep_size_length);
      vlen = MAX(vlen, fluid3_fast_variables.ep_value_length);
      break;
#endif
#ifdef D_ALE
    case el_ale2:
      slen = MAX(slen, ale2_variables.ep_size_length);
      vlen = MAX(vlen, ale2_variables.ep_value_length);
      break;
#endif
#ifdef D_ALE
    case el_ale3:
      slen = MAX(slen, ale3_variables.ep_size_length);
      vlen = MAX(vlen, ale3_variables.ep_value_length);
      break;
#endif
#ifdef D_WALL1
    case el_wall1:
      slen = MAX(slen, wall1_variables.ep_size_length);
      vlen = MAX(vlen, wall1_variables.ep_value_length);
      break;
#endif
#ifdef D_BEAM3
    case el_beam3:
      slen = MAX(slen, beam3_variables.ep_size_length);
      vlen = MAX(vlen, beam3_variables.ep_value_length);
      break;
#endif
#ifdef D_AXISHELL
    case el_axishell:
      slen = MAX(slen, axishell_variables.ep_size_length);
      vlen = MAX(vlen, axishell_variables.ep_value_length);
      break;
#endif
#ifdef D_INTERF
    case el_interf:
      slen = MAX(slen, interf_variables.ep_size_length);
      vlen = MAX(vlen, interf_variables.ep_value_length);
      break;
#endif
#ifdef D_WALLGE
    case el_wallge:
      slen = MAX(slen, wallge_variables.ep_size_length);
      vlen = MAX(vlen, wallge_variables.ep_value_length);
      break;
#endif
    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(&slen, size_length, 1, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
  MPI_Allreduce(&vlen, value_length, 1, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#else
  *size_length = slen;
  *value_length = vlen;
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the mesh connectivity.

  The mesh connectivity consists of all the ids of those nodes that
  are connected to one particular element. So the task here is to find
  the maximum number of nodes per element. This a done looping all
  elements.

  In the parallel case we need to communicate.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_mesh_item_length(struct _BIN_OUT_FIELD* context,
                           INT* value_length,
                           INT* size_length)
{
  INT i;
  INT numnp = 0;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("find_mesh_item_length");
#endif

  *value_length = 0;
  *size_length = 0;

  /* We simply need to find the maximum element node number in the
   * discretization. */

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele;

    actele = actpdis->element[i];
    numnp = MAX(numnp, actele->numnp);
  }

#ifdef PARALLEL
  MPI_Allreduce(&numnp, size_length, 1, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#else
  *size_length = numnp;
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the node coordinates.

  The number of coordinates per node equals the number of
  dimensions. The value entry length is therefore easy to
  determine. But the size entry contains (a) the node's global Id and
  (b) the number of elements connected to this node and (c) their
  local Ids.

  In the parallel case we need to communicate.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void find_coords_item_length(struct _BIN_OUT_FIELD* context,
                             INT* value_length,
                             INT* size_length)
{
  INT i;
  INT numele = 0;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("find_coords_item_length");
#endif

  *value_length = genprob.ndim;
  *size_length = 0;

  for (i=0; i<actpdis->numnp; ++i)
  {
    NODE* actnode;

    actnode = actpdis->node[i];
    numele = MAX(numele, actnode->numele);
  }

#ifdef PARALLEL
  MPI_Allreduce(&numele, size_length, 1, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#else
  *size_length = numele;
#endif

  /* plus one for the global node Id and one for the number of elements */
  *size_length += 2;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the number of double and integer values that are needed
  to store the element stresses.

  We want to store the real stress array for each
  element. Postprocessing can be done later.

  In the parallel case we need to communicate.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_stress_item_length(struct _BIN_OUT_FIELD* context,
                             INT* value_length,
                             INT* size_length)
{
  INT i;
  INT length = 0;
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("find_stress_item_length");
#endif

  *value_length = 0;
  *size_length = 0;

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele;
    actele = actpdis->element[i];

    switch (actele->eltyp)
    {
#ifdef D_SHELL8
    case el_shell8:
    {
      INT rows;
      INT cols;
      SHELL8* s8;
      s8 = actele->e.s8;

      rows = s8->forces.sdim;
      cols = s8->nGP[0]*s8->nGP[1]*s8->nGP[2];

      length = MAX(length, rows*cols);
      break;
    }
#endif
#ifdef D_SHELL9
    case el_shell9:             /* multi layer shell element */
    {
      SHELL9* s9;
      s9 = actele->e.s9;

      /* There are 6 stress components. */
      length = MAX(length, 6*context->s9_layers*s9->nGP[0]*s9->nGP[1]*s9->nGP[2]);
      break;
    }
#endif
#ifdef D_BRICK1
    case el_brick1:
    {
      INT rows;
      INT cols;
      BRICK1* c1 = actele->e.c1;

      rows = c1->stress_GP.sdim;
      cols = c1->nGP[0]*c1->nGP[1]*c1->nGP[2];
      length = MAX(length, rows*cols);
      break;
    }
#endif
#ifdef D_FLUID2
    case el_fluid2:
      break;
#endif
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
      break;
#endif
#ifdef D_FLUID2TU
    case el_fluid2_tu:
      break;
#endif
#ifdef D_FLUID3
    case el_fluid3:
      break;
#endif
#ifdef D_FLUID3_F
    case el_fluid3_fast:
      break;
#endif
#ifdef D_ALE
    case el_ale2:
      break;
#endif
#ifdef D_ALE
    case el_ale3:
      break;
#endif
#ifdef D_WALL1
    case el_wall1:
    {
      WALL1* w1 = actele->e.w1;
      if (actele->distyp == tri3)
      {
        length = MAX(length, 4);
      }
      else if (actele->distyp == tri6)
      {
        dserror("yet to be done");
      }
      else if ((actele->distyp == quad4) ||
               (actele->distyp == quad8) ||
               (actele->distyp == quad9))
      {
        length = MAX(length, 9*w1->nGP[0]*w1->nGP[1]);
      }
      else
      {
        dserror("distyp %d unsupported", actele->distyp);
      }
      break;
    }
#endif
#ifdef D_BEAM3
    case el_beam3:
    {
      BEAM3* b3 = actele->e.b3;
      length = MAX(length, 6*b3->nGP[0]);
      break;
    }
#endif
#ifdef D_AXISHELL
    case el_axishell:
      break;
#endif
#ifdef D_INTERF
    case el_interf:
    {
      INTERF* interf = actele->e.interf;
      length = MAX(length, 5*interf->nGP);
      break;
    }
#endif
#ifdef D_WALLGE
    case el_wallge:
    {
      WALLGE* wallge = actele->e.wallge;
      if ((actele->distyp == tri3) || (actele->distyp == tri6))
      {
        dserror("yet to be done");
      }
      else if ((actele->distyp == quad4) ||
               (actele->distyp == quad8) ||
               (actele->distyp == quad9))
      {
        length = MAX(length, 6*wallge->nGP[0]*wallge->nGP[1]);
      }
      else
      {
        dserror("distyp %d unsupported", actele->distyp);
      }
      break;
    }
#endif
    default:
      dserror("element type %d not supported", actele->eltyp);
    }
  }

#ifdef PARALLEL
  MPI_Allreduce(&length, value_length, 1, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#else
  *value_length = length;
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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

  There are three functions that do the element handling during
  restart. These are ``out_pack_restart_element``,
  ``in_unpack_restart_element`` and ``find_restart_item_length``. They
  are responsible for gathering the element data, scattering the
  element data and finding the element data's size respectively. If
  you change one of those you probably want to change the others
  too. These functions are very closely coupled, they share one
  structure.

  And these are the functions you have to change when you want to add
  or change ccarat's elements.

  \param context      (i) pointer to an already set up output context
  \param value_length (o) the number of double to store per element
  \param size_length  (o) the number of integer to store per element

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
void find_restart_item_length(struct _BIN_OUT_FIELD* context,
                              INT* value_length,
                              INT* size_length)
{
  INT i;
  INT local_lengths[2];
#ifdef PARALLEL
  INT global_lengths[2];
#else
  INT* global_lengths = local_lengths;
#endif
  PARTDISCRET* actpdis = &(context->actpart->pdis[context->disnum]);

#ifdef DEBUG
  dstrc_enter("find_restart_item_length");
#endif

  /* Start at zero. */
  *value_length = 0;
  *size_length  = 0;

  for (i=0; i<actpdis->numele; ++i)
  {
    ELEMENT* actele;
    actele = actpdis->element[i];

    switch (actele->eltyp)
    {
#ifdef D_SHELL8
    case el_shell8:
    {
      INT len = 0;
      SHELL8* s8;

      s8 = actele->e.s8;

      /* collect the array sizes */
      if (s8->nhyb)
      {
        len += s8->alfa.fdim * s8->alfa.sdim;
        len += s8->Dtildinv.fdim * s8->Dtildinv.sdim;
        len += s8->Lt.fdim * s8->Lt.sdim;
        len += s8->Rtilde.fdim * s8->Rtilde.sdim;
      }
      if (mat[actele->mat-1].mattyp==m_viscohyper)
      {
        ARRAY4D *a;
        a = s8->his1;
        len += a->fdim * a->sdim * a->tdim * a->fodim;
      }

      *value_length = MAX(*value_length, len);
      break;
    }
#endif

#ifdef D_SHELL9
    case el_shell9:             /* multi layer shell element */
    {
      INT len = 0;
      INT slen = 0;
      SHELL9* s9;

      s9 = actele->e.s9;

      /* collect the array sizes */
      if (s9->nhyb)
      {
        len += s9->alfa.fdim * s9->alfa.sdim;
        len += s9->Dtildinv.fdim * s9->Dtildinv.sdim;
        len += s9->L.fdim * s9->L.sdim;
        len += s9->Rtilde.fdim * s9->Rtilde.sdim;
      }
      if (s9->elewa != NULL)
      {
        INT kl;
        INT actlay = 0;

        for (kl=0; kl<s9->num_klay; kl++)
        {
          INT num_mlay;
          INT ml;

          num_mlay = s9->kinlay[kl].num_mlay;
          for (ml=0; ml<num_mlay; ml++)
          {

            /* check if there is a ipwa for this layer */
            if (s9->elewa[actlay].ipwa != NULL)
            {
              INT ngauss;
              MULTIMAT *actmultimat;

              actmultimat = &(multimat[s9->kinlay[kl].mmatID[ml]-1]);

              /* number of gausspoints in one layer*/
              ngauss = s9->nGP[0]*s9->nGP[1]*s9->nGP[2];

              switch (actmultimat->mattyp)
              {
              case m_pl_mises:
              case m_pl_dp:
                len += ngauss*(1+6*3);
                slen += ngauss;
                break;
              case m_pl_epc:
                len += ngauss*(2+6*2);
                slen += ngauss;
                break;
              case m_pl_hoff:
                len += ngauss*(1+6*4+9);
                slen += ngauss;
                break;
              default:
                printf("unsupported material %d for shell9 element",
                       actmultimat->mattyp);
              }
            }
            actlay += 1;
          }
        }
      }
      *value_length = MAX(*value_length, len);
      *size_length  = MAX(*size_length, slen);
      break;
    }
#endif

#ifdef D_BRICK1
    case el_brick1:
    {
      /* nothing? */
      break;
    }
#endif

#ifdef D_FLUID2
    case el_fluid2:
    {
      INT len = 3;
      FLUID_DYNAMIC *fdyn;

      /* only valid for single field problem (?) */
      /* I don't think so. */
      fdyn = alldyn[genprob.numff].fdyn;

      /* We have the stabilization parameter history (tau) and maybe kappa. */
      if (fdyn->surftens > 0)
      {
        len = MAX(len, 3+2*actele->numnp);
      }

      *value_length = MAX(*value_length, len);
      break;
    }
#endif

#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      INT numpdof;
      switch (actele->e.f2pro->dm)
      {
      case dm_q2pm1:
	numpdof = 3;
	break;
      case dm_q1p0:
        numpdof = 1;
        break;
      default:
	dserror("unsupported discretization mode %d", actele->e.f2pro->dm);
      }

      *value_length = MAX(*value_length, 2*numpdof);
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      INT numpdof;
      switch (actele->e.f3pro->dm)
      {
      case dm_q2pm1:
	numpdof = 4;
	break;
      case dm_q1p0:
        numpdof = 1;
        break;
      default:
	dserror("unsupported discretization mode %d", actele->e.f3pro->dm);
      }

      *value_length = MAX(*value_length, 2*numpdof);
      break;
    }
#endif
#ifdef D_FLUID2TU
    case el_fluid2_tu:
      break;
#endif
#ifdef D_FLUID3
    case el_fluid3:
      *value_length = MAX(*value_length, 3);
      break;
#endif
#ifdef D_FLUID3_F
    case el_fluid3_fast:
      *value_length = MAX(*value_length, 3);
      break;
#endif
#ifdef D_ALE
    case el_ale2:
      break;
#endif
#ifdef D_ALE
    case el_ale3:
      break;
#endif

#ifdef D_WALL1
    case el_wall1:
    {
      WALL1* w1 = actele->e.w1;

      /* There are huge differences depending on the material. */

      if (w1->elewa != NULL)
      {
        INT ngauss;
        ngauss = w1->nGP[0]*w1->nGP[1];
        switch (mat[actele->mat-1].mattyp)
        {
        case m_pl_mises:
          *value_length = MAX(*value_length, ngauss*(1+4*3));
          *size_length  = MAX(*size_length,  ngauss*(1));
          break;
        case m_pl_mises_3D:
          *value_length = MAX(*value_length, ngauss*(1+4*7));
          *size_length  = MAX(*size_length,  ngauss*(1));
          break;
        case m_pl_epc3D:
          *value_length = MAX(*value_length, ngauss*(2+4*7));
          *size_length  = MAX(*size_length,  ngauss*(1));
          break;
        default:
          /*
           * This is strange. Do we want these messages when we
           * restart wall meshed with linear material? */
          printf("unsupported material %d for wall1 element",
                 mat[actele->mat-1].mattyp);
        }
      }
      break;
    }
#endif

#ifdef D_BEAM3
    case el_beam3:
    {
      BEAM3* b3 = actele->e.b3;
      if(mat[actele->mat-1].mattyp == m_pl_mises ||
         mat[actele->mat-1].mattyp == m_pl_dp ||
         mat[actele->mat-1].mattyp == m_pl_epc)
      {
        *value_length = MAX(*value_length, b3->nGP[0]*5*5*40);
      }
      break;
    }
#endif
#ifdef D_AXISHELL
    case el_axishell:
      break;
#endif

#ifdef D_INTERF
    case el_interf:
    {
      INTERF* interf = actele->e.interf;
      INT ngauss = interf->nGP;
      *value_length = MAX(*value_length, ngauss*10);
      break;
    }
#endif

#ifdef D_WALLGE
    case el_wallge:
    {
      /*WALLGE* wallge = actele->e.wallge;*/
      /* There's nothing to be saved for wallge elements, right? */
      break;
    }
#endif

    default:
      dserror("element type %d unsupported", actele->eltyp);
    }
  }

  /* Finally communicate what we've found. We need to do this because
   * there's no guarantee that the most demanding element type is used
   * on all processors. */
  local_lengths[0] = *value_length;
  local_lengths[1] = *size_length;
#ifdef PARALLEL
  MPI_Allreduce(local_lengths, global_lengths, 2, MPI_INT, MPI_MAX, context->actintra->MPI_INTRA_COMM);
#endif
  *value_length = global_lengths[0];
  *size_length  = global_lengths[1];

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*======================================================================*/
/* Largely unrelated utility functions */
/*======================================================================*/


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
INT get_field_position(struct _BIN_OUT_FIELD  *context)
{
  INT field_pos;

#ifdef DEBUG
  dstrc_enter("get_field_position");
#endif

  for (field_pos=0; field_pos<genprob.numfld; ++field_pos)
  {
    if (context->actfield == &(field[field_pos]))
    {
      break;
    }
  }

  if (field_pos==genprob.numfld)
  {
    dserror("unregistered field object");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return field_pos;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the start of a group that describes one aspect of a
  discretization.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void out_main_group_head(struct _BIN_OUT_FIELD  *context, CHAR* name)
{
  INT field_pos;
#ifdef DEBUG
  dstrc_enter("out_main_group_head");
#endif

  field_pos = get_field_position(context);

  fprintf(bin_out_main.control_file, "%s:\n"
          "    field = \"%s\"\n"
          "    field_pos = %d\n"
          "    discretization = %d\n"
          "",
          name,
          fieldnames[context->actfield->fieldtyp],
          field_pos,
          context->disnum);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*! @} (documentation module close)*/
#endif
