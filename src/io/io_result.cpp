/*!
\file
\brief The public result output function.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Everybody who wants to output results calls the function defined here.

\author u.kue
\date 11/04

*/

#ifdef BINIO

/*!
  \addtogroup IO
*//*! @{ (documentation module open)*/

extern "C" {

#include "../headers/standardtypes.h"

}

#include "io.h"
#include "io_packing.h"
#include "io_singlefile.h"


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

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
  *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
  \brief ranks and communicators

  <pre>                                                         m.gee 8/00
  This structure struct _PAR par; is defined in main_ccarat.c
  and the type is in partition.h
  </pre>

  *----------------------------------------------------------------------*/
extern struct _PAR   par;


/*----------------------------------------------------------------------*/
/*!
  \brief Write all results for one step.

  Write results for postprocessing. All algorithms call this function
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
                 OUT_FLAGS flags)
{
  INT rank = context->actintra->intra_rank;

  /* We remember the step for with we are called. This way we notice
   * when a new time step comes and output the appropriate header. */
  static DOUBLE last_time = -1;
  static INT last_step = -1;
  static INT last_discr = -1;
  static INT last_field_pos = -1;

  INT field_pos;

#ifdef DEBUG
  dstrc_enter("out_results");
#endif

#ifdef PERF
  perf_begin(71);
#endif

  dsassert(flags != 0, "no output flags; I've got nothing to write");

  /* It's crucial to make sure everything that's written will go to
   * the result files (not the restart ones). */
  out_activate_result(context);

  field_pos = get_field_position(context);

  /* Only start a new group if we have a new time step. This way
   * out_results can be called more often that once per step. Be
   * careful, however, not to mix it with restart output
   * calls. Restart starts a new group as well and everything will
   * be added to the latest group no matter who created it. */

  if ((last_step != step) || (last_time != time) ||
      (last_field_pos != field_pos) || (last_discr != context->disnum))
  {
    last_step = step;
    last_time = time;
    last_field_pos = field_pos;
    last_discr = context->disnum;

    if (rank == 0)
    {
      out_main_group_head(context, "result");
      fprintf(bin_out_main.control_file,
              "    time = %f\n"
              "    step = %d\n"
              "\n",
              time, step);
    }

    /* If we've reached the number of steps to put into one file
     * open a new one. */
    if ((context->result_count % bin_out_main.steps_per_file) == 0)
    {
      out_open_data_files(context, context->out, "result", context->result_count);
    }

    context->result_count += 1;
  }

  /*--------------------------------------------------------------------*/
  /* Now let's see what we have to write. */

  if (flags & OUTPUT_DISPLACEMENT)
  {

#ifdef D_SHELL8
    /* shell8 elements have three additional dofs, all six living in
     * one node array's row. Thus we only need one call to output
     * them. */
    if (context->is_shell8_problem)
    {
      out_node_chunk(context, "displacement", cc_displacement, 6, 0, place);
    }
    else
#endif

      /* Displacement of normal nodes. */
      out_node_chunk(context, "displacement", cc_displacement, genprob.ndim, 0, place);

#ifdef D_SHELL9
    /* shell9 elements have layers with their own displacements. */
    if (context->is_shell9_problem)
    {
      INT layer_count;

      if (context->s9_numnp == 4)
      {
        layer_count = 1;
      }
      else
      {
        layer_count = 2;
      }

      out_element_chunk(context, "shell9_displacement", cc_shell9_displacement,
                        3*context->s9_numnp*(layer_count*context->s9_layers+1), 0,
                        place);
    }
#endif
  }

  if (flags & OUTPUT_VELOCITY)
  {
    out_node_chunk(context, "velocity", cc_velocity, genprob.ndim, 0, place);
  }

  if (flags & OUTPUT_PRESSURE)
  {
    out_node_chunk(context, "pressure", cc_pressure, 1, 0, place);
  }

  if (flags & OUTPUT_AV_PRESSURE)
  {
    out_node_chunk(context, "average_pressure", cc_av_pressure, 1, 0, place);
  }

  if (flags & OUTPUT_STRESS)
  {
    INT value_length;
    INT size_length;

    find_stress_item_length(context, &value_length, &size_length);
    if ((value_length > 0) || (size_length > 0))
    {
      out_element_chunk(context, "stress", cc_stress, value_length, size_length, 0);
    }
  }

  if (flags & OUTPUT_CONTACT)
  {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "contact" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_EIGENMODES)
  {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "eigenmodes" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_THICKNESS)
  {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "thickness" RED_LIGHT " not yet supported\n" END_COLOR);
  }

  if (flags & OUTPUT_AXI_LOADS)
  {
    printf(RED_LIGHT "binary output of " GREEN_LIGHT "axishell" RED_LIGHT " loads not yet supported\n" END_COLOR);
  }

  if (par.myrank == 0)
  {
    fflush(bin_out_main.control_file);
  }

#ifdef PARALLEL
  /* How to flush in parallel? */
#else
  fflush(context->out->value_file);
  fflush(context->out->size_file);
#endif

#ifdef PERF
  perf_end(71);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*! @} (documentation module close)*/
#endif
