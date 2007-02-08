/*!
\file
\brief Restart functions based on the binary io.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Here the restart functions for all algorithms are collected. These
functions are supposed to be rather high level and easy to
understand. Feel free to add your own.

\author u.kue
\date 08/04

*/

#ifdef BINIO

/*!
\addtogroup IO
*//*! @{ (documentation module open)*/

extern "C" {

#include "../headers/standardtypes.h"

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

}

#include "io_packing.h"
#include "io_singlefile.h"


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
                                    INT nwork, DIST_VECTOR *work)
{
  INT rank = context->actintra->intra_rank;
  MAP map;
  INT i;
  INT value_length;
  INT size_length;

#ifdef DEBUG
  dstrc_enter("restart_write_bin_nlnstructdyn");
#endif

#ifdef PERF
  perf_begin(72);
#endif

  /* It's crucial to make sure everything that's written will go to
   * the restart files (not the result ones). */
  out_activate_restart(context);

  if (rank == 0) {
    out_main_group_head(context, "restart");
    fprintf(bin_out_main.control_file,
            "    step = %d\n"
            "    time = %20.20f\n"
            "\n",
            sdyn->step, sdyn->time);
  }

  /* If we've reached the number of steps to put into one file
   * open a new one. */
  if ((context->restart_count % bin_out_main.steps_per_file) == 0)
  {
    out_open_data_files(context, context->out, "restart", context->restart_count);
  }

  context->restart_count += 1;

  /*--------------------------------------------------------------------*/
  /* Write the control structures. */

  /* This is a little too explicit. However the old ccarat restart is
   * that verbose (for whatever reason) and this explicitness is the
   * only way to store the structures in a portable fashion. An
   * alternative would be to write a big hex dump. But it would be
   * tied to the (compiler dependent) memory layout of the
   * structures.
   *
   * The drawback of this explicit scheme is that we have to adjust
   * whenever the structures change. */

  if (rank == 0) {
    /* I use an extra map here to keep name and type together. */
    init_map(&map);

    map_insert_int_cpy(&map, sdyn->Typ, "type");
    map_insert_int_cpy(&map, sdyn->updevry_disp, "updevry_disp");
    map_insert_int_cpy(&map, sdyn->updevry_stress, "updevry_stress");
    map_insert_int_cpy(&map, sdyn->res_write_evry, "res_write_evry");
    map_insert_int_cpy(&map, sdyn->nstep, "nstep");
    map_insert_int_cpy(&map, sdyn->step, "step");
    map_insert_int_cpy(&map, sdyn->damp, "damp");
    map_insert_int_cpy(&map, sdyn->iter, "iter");
    map_insert_int_cpy(&map, sdyn->maxiter, "maxiter");
    map_insert_int_cpy(&map, sdyn->eigen, "eigen");
    map_insert_int_cpy(&map, sdyn->contact, "contact");
    map_insert_real_cpy(&map, sdyn->toldisp, "toldisp");
    map_insert_real_cpy(&map, sdyn->dt, "dt");
    map_insert_real_cpy(&map, sdyn->maxtime, "maxtime");
    map_insert_real_cpy(&map, sdyn->time, "time");
    map_insert_real_cpy(&map, sdyn->beta, "beta");
    map_insert_real_cpy(&map, sdyn->gamma, "gamma");
    map_insert_real_cpy(&map, sdyn->alpha_m, "alpha_m");
    map_insert_real_cpy(&map, sdyn->alpha_f, "alpha_f");
#ifdef GEMM
    map_insert_real_cpy(&map, sdyn->xsi, "xsi");
#endif
    map_insert_real_cpy(&map, sdyn->m_damp, "m_damp");
    map_insert_real_cpy(&map, sdyn->k_damp, "k_damp");

    map_insert_int_cpy(&map, sdyn->timeadapt, "timeadapt");
    map_insert_int_cpy(&map, sdyn->itwant, "itwant");
    map_insert_real_cpy(&map, sdyn->maxdt, "maxdt");
    map_insert_real_cpy(&map, sdyn->resultdt, "resultdt");
    map_insert_int_cpy(&map, sdyn->writecounter, "writecounter");

    /* output to the control file */
    fprintf(bin_out_main.control_file, "    sdyn:\n");
    fprintf(bin_out_main.control_file,
            "        # Maybe we don't need to store all of these variables.\n"
            "        # Please check.\n");
    map_print(bin_out_main.control_file, &map, 8);

    destroy_map(&map);

    /* second structure: dynvar*/
    init_map(&map);

    map_insert_real_cpy(&map, dynvar->rldfac, "rldfac");
    map_insert_real_cpy(&map, dynvar->rnorm, "rnorm");

    for (i=0; i<20; ++i) {
      CHAR buf[15];
      sprintf(buf, "constants_%d", i);
      map_insert_real_cpy(&map, dynvar->constants[i], buf);
    }

    map_insert_real_cpy(&map, dynvar->epot, "epot");
    map_insert_real_cpy(&map, dynvar->eout, "eout");
    map_insert_real_cpy(&map, dynvar->etot, "etot");
    map_insert_real_cpy(&map, dynvar->ekin, "ekin");

#ifdef D_WALL1
    map_insert_real_cpy(&map, dynvar->total_linmom[0], "total_linmom_0");
    map_insert_real_cpy(&map, dynvar->total_linmom[1], "total_linmom_1");
    map_insert_real_cpy(&map, dynvar->total_angular_momentum, "total_angular_momentum");
    map_insert_real_cpy(&map, dynvar->total_strain_energy, "total_strain_energy");
    map_insert_real_cpy(&map, dynvar->total_kinetic_energy, "total_kinetic_energy");
    map_insert_real_cpy(&map, dynvar->local_linmom[0], "local_linmom_0");
    map_insert_real_cpy(&map, dynvar->local_linmom[1], "local_linmom_1");
    map_insert_real_cpy(&map, dynvar->local_angular_momentum, "local_angular_momentum");
    map_insert_real_cpy(&map, dynvar->local_strain_energy, "local_strain_energy");
    map_insert_real_cpy(&map, dynvar->local_kinetic_energy, "local_kinetic_energy");
#endif
    map_insert_real_cpy(&map, dynvar->dinorm, "dinorm");
    map_insert_real_cpy(&map, dynvar->dnorm, "dnorm");

    /* output to the control file */
    fprintf(bin_out_main.control_file, "    dynvar:\n");
    fprintf(bin_out_main.control_file,
            "        # Maybe we don't need to store all of these variables.\n"
            "        # Please check.\n");
    map_print(bin_out_main.control_file, &map, 8);

    destroy_map(&map);
  }

  /*--------------------------------------------------------------------*/
  /* Write the distributed vectors. */

  out_distvec_chunk(context, "rhs_vec",  nrhs,  rhs);
  out_distvec_chunk(context, "sol_vec",  nsol,  sol);
  out_distvec_chunk(context, "disp_vec", ndis,  dispi);
  out_distvec_chunk(context, "vel_vec",  nvel,  vel);
  out_distvec_chunk(context, "acc_vec",  nacc,  acc);
  out_distvec_chunk(context, "fie_vec",  nfie,  fie);
  out_distvec_chunk(context, "work_vec", nwork, work);

  /*--------------------------------------------------------------------*/
  /* Write the node arrays. */

  out_node_arrays(context, node_array_sol);
  out_node_arrays(context, node_array_sol_increment);
  out_node_arrays(context, node_array_sol_residual);
  if (context->max_size[node_array_sol_mf] > 0) {
    out_node_arrays(context, node_array_sol_mf);
  }

  /*--------------------------------------------------------------------*/
  /* Write the element data. */

  /* This call looks very easy. But you need to go down and change
   * functions at the bottom of the io module if you want to add new
   * elements. */

  find_restart_item_length(context, &value_length, &size_length);
  out_element_chunk(context, "element_data", cc_restart_element,
                    value_length, size_length, 0);

  if (rank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifndef PARALLEL
  fflush(context->out_restart.value_file);
  fflush(context->out_restart.size_file);
#endif

#ifdef PERF
  perf_end(72);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                                   INT              disnum,
                                   INTRA           *actintra,
                                   INT nrhs,  DIST_VECTOR *rhs,
                                   INT nsol,  DIST_VECTOR *sol,
                                   INT ndis,  DIST_VECTOR *dispi,
                                   INT nvel,  DIST_VECTOR *vel,
                                   INT nacc,  DIST_VECTOR *acc,
                                   INT nfie,  DIST_VECTOR *fie,
                                   INT nwork, DIST_VECTOR *work,
                                   INT step)
{
  BIN_IN_FIELD context;
  MAP* result_info;
  MAP* map;
  INT i;

#ifdef DEBUG
  dstrc_enter("restart_read_bin_nlnstructdyn");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize input of this field */

  init_bin_in_field(&context, sysarray_typ, sysarray, actfield, actpart, actintra, disnum);

  /*--------------------------------------------------------------------*/
  /* find the describtion of the restart data in the control file */

  result_info = in_find_restart_group(&context, disnum, step);

  /*--------------------------------------------------------------------*/
  /* additional values that must be read */

  /* Get a pointer to the first group. We don't own the
   * group. Therefore we can savely forget the pointer afterwards. */
  map = map_read_map(result_info, "sdyn");

  // no means to cast to an anonymous enum. Lets rely on the normal
  // input to set the right type.
  //sdyn->Typ = map_read_int(map, "type");
  sdyn->updevry_disp = map_read_int(map, "updevry_disp");
  sdyn->updevry_stress = map_read_int(map, "updevry_stress");
  sdyn->res_write_evry = map_read_int(map, "res_write_evry");
  sdyn->nstep = map_read_int(map, "nstep");
  sdyn->step = map_read_int(map, "step");
  sdyn->damp = map_read_int(map, "damp");
  sdyn->iter = map_read_int(map, "iter");
  sdyn->maxiter = map_read_int(map, "maxiter");
  sdyn->eigen = map_read_int(map, "eigen");
  sdyn->contact = map_read_int(map, "contact");
  sdyn->toldisp = map_read_real(map, "toldisp");
  sdyn->dt = map_read_real(map, "dt");
  sdyn->maxtime = map_read_real(map, "maxtime");
  sdyn->time = map_read_real(map, "time");
  sdyn->beta = map_read_real(map, "beta");
  sdyn->gamma = map_read_real(map, "gamma");
  sdyn->alpha_m = map_read_real(map, "alpha_m");
  sdyn->alpha_f = map_read_real(map, "alpha_f");
#ifdef GEMM
  sdyn->xsi = map_read_real(map, "xsi");
#endif
  sdyn->m_damp = map_read_real(map, "m_damp");
  sdyn->k_damp = map_read_real(map, "k_damp");

  sdyn->timeadapt = map_read_int(map, "timeadapt");
  sdyn->itwant = map_read_int(map, "itwant");
  sdyn->maxdt = map_read_real(map, "maxdt");
  sdyn->resultdt = map_read_real(map, "resultdt");
  sdyn->writecounter = map_read_int(map, "writecounter");

  /* Get a pointer to the second group. */
  map = map_read_map(result_info, "dynvar");

  dynvar->rldfac = map_read_real(map, "rldfac");
  dynvar->rnorm = map_read_real(map, "rnorm");

  for (i=0; i<20; ++i) {
    CHAR buf[15];
    sprintf(buf, "constants_%d", i);
    dynvar->constants[i] = map_read_real(map, buf);
  }

  dynvar->epot = map_read_real(map, "epot");
  dynvar->eout = map_read_real(map, "eout");
  dynvar->etot = map_read_real(map, "etot");
  dynvar->ekin = map_read_real(map, "ekin");

#ifdef D_WALL1
  dynvar->total_linmom[0] = map_read_real(map, "total_linmom_0");
  dynvar->total_linmom[1] = map_read_real(map, "total_linmom_1");
  dynvar->total_angular_momentum = map_read_real(map, "total_angular_momentum");
  dynvar->total_strain_energy = map_read_real(map, "total_strain_energy");
  dynvar->total_kinetic_energy = map_read_real(map, "total_kinetic_energy");
  dynvar->local_linmom[0] = map_read_real(map, "local_linmom_0");
  dynvar->local_linmom[1] = map_read_real(map, "local_linmom_1");
  dynvar->local_angular_momentum = map_read_real(map, "local_angular_momentum");
  dynvar->local_strain_energy = map_read_real(map, "local_strain_energy");
  dynvar->local_kinetic_energy = map_read_real(map, "local_kinetic_energy");
#endif
  dynvar->dinorm = map_read_real(map, "dinorm");
  dynvar->dnorm = map_read_real(map, "dnorm");

  /*--------------------------------------------------------------------*/
  /* read the distributed vectors */

  in_distvec_chunk(&context, result_info, "rhs_vec",  nrhs,  rhs);
  in_distvec_chunk(&context, result_info, "sol_vec",  nsol,  sol);
  in_distvec_chunk(&context, result_info, "disp_vec", ndis,  dispi);
  in_distvec_chunk(&context, result_info, "vel_vec",  nvel,  vel);
  in_distvec_chunk(&context, result_info, "acc_vec",  nacc,  acc);
  in_distvec_chunk(&context, result_info, "fie_vec",  nfie,  fie);
  in_distvec_chunk(&context, result_info, "work_vec", nwork, work);

  /*--------------------------------------------------------------------*/
  /* read the individual node arrays */

  in_node_arrays(&context, result_info, node_array_sol);
  in_node_arrays(&context, result_info, node_array_sol_increment);
  in_node_arrays(&context, result_info, node_array_sol_residual);
  if (map_has_map(result_info, "sol_mf")) {
    in_node_arrays(&context, result_info, node_array_sol_mf);
  }

  /*--------------------------------------------------------------------*/
  /* read element data */

  /* This call looks very easy. But you need to go down and change
   * functions at the bottom of the io module if you want to add new
   * elements. */

  in_element_chunk(&context, result_info, "element_data", cc_restart_element);

  /*--------------------------------------------------------------------*/
  /* done with this field */

  destroy_bin_in_field(&context);

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                                     INT ndis,  DIST_VECTOR *dispi)
{
  INT rank = context->actintra->intra_rank;
  MAP map;
  INT i;
  INT value_length;
  INT size_length;

#ifdef DEBUG
  dstrc_enter("restart_write_bin_nlnstructstat");
#endif

#ifdef PERF
  perf_begin(72);
#endif

  /* It's crucial to make sure everything that's written will go to
   * the restart files (not the result ones). */
  out_activate_restart(context);

  if (rank == 0) {
    out_main_group_head(context, "restart");
    fprintf(bin_out_main.control_file,
            "    step = %d\n"
            "    time = %20.20f\n"
            "\n",
            kstep, 0.0);
  }

  /* If we've reached the number of steps to put into one file
   * open a new one. */
  if ((context->restart_count % bin_out_main.steps_per_file) == 0)
  {
    out_open_data_files(context, context->out, "restart", context->restart_count);
  }

  context->restart_count += 1;

  /*--------------------------------------------------------------------*/
  /* Write the control structures. */

  if (rank == 0) {
    /* I use an extra map here to keep name and type together. */
    init_map(&map);

    map_insert_int_cpy(&map, statvar->linear, "linear");
    map_insert_int_cpy(&map, statvar->nonlinear, "nonlinear");
    map_insert_int_cpy(&map, statvar->kintyp, "kintyp");
    map_insert_int_cpy(&map, statvar->nr_controltyp, "nr_controltyp");
    map_insert_int_cpy(&map, statvar->nstep, "nstep");
    map_insert_int_cpy(&map, statvar->maxiter, "maxiter");
    map_insert_real_cpy(&map, statvar->tolresid, "tolresid");
    map_insert_real_cpy(&map, statvar->toldisp, "toldisp");
    map_insert_real_cpy(&map, statvar->stepsize, "stepsize");
    map_insert_int_cpy(&map, statvar->iarc, "iarc");
    map_insert_real_cpy(&map, statvar->arcscl, "arcscl");
    map_insert_int_cpy(&map, statvar->signchcsp, "signchcsp");
    map_insert_int_cpy(&map, statvar->resevry_disp, "resevry_disp");
    map_insert_int_cpy(&map, statvar->resevry_stress, "resevry_stress");
    map_insert_int_cpy(&map, statvar->resevery_restart, "resevery_restart");
    map_insert_int_cpy(&map, statvar->graderw, "graderw");

    map_insert_int_cpy(&map, statvar->control_node_global, "control_node_global");
    map_insert_int_cpy(&map, statvar->control_dof, "control_dof");

    map_insert_int_cpy(&map, statvar->isrelstepsize, "isrelstepsize");
    for (i=0; i<20; ++i) {
      CHAR buf[15];
      sprintf(buf, "actstep_%d", i);
      map_insert_int_cpy(&map, statvar->actstep[i], buf);
    }
    for (i=0; i<20; ++i) {
      CHAR buf[20];
      sprintf(buf, "actstepsize_%d", i);
      map_insert_real_cpy(&map, statvar->actstepsize[i], buf);
    }
    map_insert_int_cpy(&map, statvar->numcurve, "numcurve");

    for (i=0; i<6; ++i) {
      CHAR buf[20];
      sprintf(buf, "reldisnode_ID_%d", i);
      map_insert_int_cpy(&map, statvar->reldisnode_ID[i], buf);
    }
    for (i=0; i<6; ++i) {
      CHAR buf[15];
      sprintf(buf, "reldis_dof_%d", i);
      map_insert_int_cpy(&map, statvar->reldis_dof[i], buf);
    }

    /* output to the control file */
    fprintf(bin_out_main.control_file, "    statvar:\n");
    fprintf(bin_out_main.control_file,
            "        # Maybe we don't need to store all of these variables.\n"
            "        # Please check.\n");
    map_print(bin_out_main.control_file, &map, 8);

    destroy_map(&map);

    /* second structure: dynvar*/
    init_map(&map);

    map_insert_real_cpy(&map, nln_data->sp1, "sp1");
    map_insert_real_cpy(&map, nln_data->csp, "csp");
    map_insert_real_cpy(&map, nln_data->rlold, "rlold");
    map_insert_real_cpy(&map, nln_data->rlnew, "rlnew");
    map_insert_real_cpy(&map, nln_data->rlpre, "rlpre");

    map_insert_real_cpy(&map, nln_data->renorm, "renorm");
    map_insert_real_cpy(&map, nln_data->rinorm, "rinorm");
    map_insert_real_cpy(&map, nln_data->rrnorm, "rrnorm");

    map_insert_real_cpy(&map, nln_data->renergy, "renergy");

#if 0
    /* This is especially crude. But I think we need to do
     * this... even though the old restart code does nonsense with
     * this array. */
    map_insert_int_cpy(&map, nln_data->arcfac.fdim, "arcfac_fdim");
    map_insert_int_cpy(&map, nln_data->arcfac.sdim, "arcfac_sdim");
    for (i=0; i<nln_data->arcfac.fdim*nln_data->arcfac.sdim; ++i) {
      CHAR buf[15];
      sprintf(buf, "arcfac_%d", i);
      map_insert_real_cpy(&map, nln_data->arcfac.a.da[0][i], buf);
    }
#endif

    /* output to the control file */
    fprintf(bin_out_main.control_file, "    nln_data:\n");
    fprintf(bin_out_main.control_file,
            "        # Maybe we don't need to store all of these variables.\n"
            "        # Please check.\n");
    map_print(bin_out_main.control_file, &map, 8);

    destroy_map(&map);
  }

  /*--------------------------------------------------------------------*/
  /* Write the distributed vectors. */

  out_distvec_chunk(context, "rhs_vec",  nrhs,  rhs);
  out_distvec_chunk(context, "sol_vec",  nsol,  sol);
  out_distvec_chunk(context, "disp_vec", ndis,  dispi);

  /*--------------------------------------------------------------------*/
  /* Write the node arrays. */

  out_node_arrays(context, node_array_sol);
  out_node_arrays(context, node_array_sol_increment);
  out_node_arrays(context, node_array_sol_residual);

  /*--------------------------------------------------------------------*/
  /* Write the element data. */

  /* There are certain elements we expect here. */
  find_restart_item_length(context, &value_length, &size_length);
  out_element_chunk(context, "element_data", cc_restart_element,
                    value_length, size_length, 0);

  if (rank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifndef PARALLEL
  fflush(context->out_restart.value_file);
  fflush(context->out_restart.size_file);
#endif

#ifdef PERF
  perf_end(72);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                                    INT step)
{
  BIN_IN_FIELD context;
  MAP* result_info;
  MAP* map;
  INT i;

#ifdef DEBUG
  dstrc_enter("restart_read_bin_nlnstructstat");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize input of this field */

  init_bin_in_field(&context, sysarray_typ, sysarray, actfield, actpart, actintra, disnum);

  /*--------------------------------------------------------------------*/
  /* find the describtion of the restart data in the control file */

  result_info = in_find_restart_group(&context, disnum, step);

  /*--------------------------------------------------------------------*/
  /* additional values that must be read */

  /* Get a pointer to the first group. We don't own the
   * group. Therefore we can savely forget the pointer afterwards. */
  map = map_read_map(result_info, "statvar");

  statvar->linear = map_read_int(map, "linear");
  statvar->nonlinear = map_read_int(map, "nonlinear");
  statvar->kintyp = static_cast<KINTYP>(map_read_int(map, "kintyp"));
  statvar->nr_controltyp = static_cast<NR_CONTROLTYP>(map_read_int(map, "nr_controltyp"));
  statvar->nstep = map_read_int(map, "nstep");
  statvar->maxiter = map_read_int(map, "maxiter");
  statvar->tolresid = map_read_real(map, "tolresid");
  statvar->toldisp = map_read_real(map, "toldisp");
  statvar->stepsize = map_read_real(map, "stepsize");
  statvar->iarc = map_read_int(map, "iarc");
  statvar->arcscl = map_read_real(map, "arcscl");
  statvar->signchcsp = map_read_int(map, "signchcsp");
  statvar->resevry_disp = map_read_int(map, "resevry_disp");
  statvar->resevry_stress = map_read_int(map, "resevry_stress");
  statvar->resevery_restart = map_read_int(map, "resevery_restart");
  statvar->graderw = map_read_int(map, "graderw");

  statvar->control_node_global = map_read_int(map, "control_node_global");
  statvar->control_dof = map_read_int(map, "control_dof");

  /* We need to reestablish the control node pointer. */
  {
    INT cdof;
    calstatserv_findcontroldof(actfield,
                               statvar->control_node_global,
                               statvar->control_dof,
                               &(statvar->controlnode),
                               &cdof);
  }

  statvar->isrelstepsize = map_read_int(map, "isrelstepsize");
  for (i=0; i<20; ++i) {
    CHAR buf[15];
    sprintf(buf, "actstep_%d", i);
    statvar->actstep[20] = map_read_int(map, buf);
  }
  for (i=0; i<20; ++i) {
    CHAR buf[20];
    sprintf(buf, "actstepsize_%d", i);
    statvar->actstepsize[20] = map_read_real(map, buf);
  }
  statvar->numcurve = map_read_int(map, "numcurve");

  for (i=0; i<6; ++i) {
    CHAR buf[20];
    sprintf(buf, "reldisnode_ID_%d", i);
    statvar->reldisnode_ID[6] = map_read_int(map, buf);
  }
  for (i=0; i<6; ++i) {
    CHAR buf[15];
    sprintf(buf, "reldis_dof_%d", i);
    statvar->reldis_dof[6] = map_read_int(map, buf);
  }

  map = map_read_map(result_info, "nln_data");

  nln_data->sp1 = map_read_real(map, "sp1");
  nln_data->csp = map_read_real(map, "csp");
  nln_data->rlold = map_read_real(map, "rlold");
  nln_data->rlnew = map_read_real(map, "rlnew");
  nln_data->rlpre = map_read_real(map, "rlpre");

  nln_data->renorm = map_read_real(map, "renorm");
  nln_data->rinorm = map_read_real(map, "rinorm");
  nln_data->rrnorm = map_read_real(map, "rrnorm");

  nln_data->renergy = map_read_real(map, "renergy");

#if 0
  if ((map_read_int(map, "arcfac_fdim") != nln_data->arcfac.fdim) ||
      (map_read_int(map, "arcfac_sdim") != nln_data->arcfac.sdim)) {
    dserror("arcfac dimension mismatch");
  }
  for (i=0; i<nln_data->arcfac.fdim*nln_data->arcfac.sdim; ++i) {
    CHAR buf[15];
    sprintf(buf, "arcfac_%d", i);
    nln_data->arcfac.a.da[0][i] = map_read_real(map, buf);
  }
#endif

  /*--------------------------------------------------------------------*/
  /* read the distributed vectors */

  in_distvec_chunk(&context, result_info, "rhs_vec",  nrhs,  rhs);
  in_distvec_chunk(&context, result_info, "sol_vec",  nsol,  sol);
  in_distvec_chunk(&context, result_info, "disp_vec", ndis,  dispi);

  /*--------------------------------------------------------------------*/
  /* read the individual node arrays */

  in_node_arrays(&context, result_info, node_array_sol);
  in_node_arrays(&context, result_info, node_array_sol_increment);
  in_node_arrays(&context, result_info, node_array_sol_residual);

  /*--------------------------------------------------------------------*/
  /* read element data */

  in_element_chunk(&context, result_info, "element_data", cc_restart_element);

  /*--------------------------------------------------------------------*/
  /* done with this field */

  destroy_bin_in_field(&context);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write fluid restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fluiddyn(struct _BIN_OUT_FIELD   *context,
                                FLUID_DYNAMIC   *fdyn)
{
  INT rank = context->actintra->intra_rank;
  INT value_length;
  INT size_length;

#ifdef DEBUG
  dstrc_enter("restart_write_bin_fluiddyn");
#endif

#ifdef PERF
  perf_begin(72);
#endif

  /* It's crucial to make sure everything that's written will go to
   * the restart files (not the result ones). */
  out_activate_restart(context);

  if (rank == 0) {
    out_main_group_head(context, "restart");
    fprintf(bin_out_main.control_file,
            "    step = %d\n"
            "    time = %20.20f\n"
            "\n",
            fdyn->step, fdyn->acttime);
  }

  /* If we've reached the number of steps to put into one file
   * open a new one. */
  if ((context->restart_count % bin_out_main.steps_per_file) == 0)
  {
    out_open_data_files(context, context->out, "restart", context->restart_count);
  }

  context->restart_count += 1;

  out_node_arrays(context, node_array_sol);
  out_node_arrays(context, node_array_sol_increment);
  if (context->max_size[node_array_sol_mf] > 0) {
    out_node_arrays(context, node_array_sol_mf);
  }

  /*--------------------------------------------------------------------*/
  /* Write the element data. */

  /* This call looks very easy. But you need to go down and change
   * functions at the bottom of the io module if you want to add new
   * elements. */

  find_restart_item_length(context, &value_length, &size_length);
  out_element_chunk(context, "element_data", cc_restart_element,
                    value_length, size_length, 0);

  if (rank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifndef PARALLEL
  fflush(context->out_restart.value_file);
  fflush(context->out_restart.size_file);
#endif

#ifdef PERF
  perf_end(72);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                               INT              step)
{
  BIN_IN_FIELD context;
  MAP* result_info;

#ifdef DEBUG
  dstrc_enter("restart_read_bin_fluiddyn");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize input of this field */

  init_bin_in_field(&context, sysarray_typ, sysarray, actfield, actpart, actintra, disnum);

  /*--------------------------------------------------------------------*/
  /* find the describtion of the restart data in the control file */

  result_info = in_find_restart_group(&context, disnum, step);

  /*--------------------------------------------------------------------*/
  /* additional values that must be read */

  fdyn->step = step;
  fdyn->acttime = map_read_real(result_info, "time");

  /*--------------------------------------------------------------------*/
  /* read the individual node arrays */

  in_node_arrays(&context, result_info, node_array_sol);
  in_node_arrays(&context, result_info, node_array_sol_increment);
  if (map_has_map(result_info, "sol_mf")) {
    in_node_arrays(&context, result_info, node_array_sol_mf);
  }

  /*--------------------------------------------------------------------*/
  /* read element data */

  in_element_chunk(&context, result_info, "element_data", cc_restart_element);

  /*--------------------------------------------------------------------*/
  /* done with this field */

  destroy_bin_in_field(&context);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write ale restart data.

  \param context  pointer to an already set up output context

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_aledyn(struct _BIN_OUT_FIELD *context,
                              ALE_DYNAMIC *adyn)
{
  INT rank = context->actintra->intra_rank;

#ifdef DEBUG
  dstrc_enter("restart_write_bin_aledyn");
#endif

#ifdef PERF
  perf_begin(72);
#endif

  /* It's crucial to make sure everything that's written will go to
   * the restart files (not the result ones). */
  out_activate_restart(context);

  if (rank == 0) {
    out_main_group_head(context, "restart");
    fprintf(bin_out_main.control_file,
            "    step = %d\n"
            "    time = %20.20f\n"
            "\n",
            adyn->step, adyn->time);
  }

  /* If we've reached the number of steps to put into one file
   * open a new one. */
  if ((context->restart_count % bin_out_main.steps_per_file) == 0)
  {
    out_open_data_files(context, context->out, "restart", context->restart_count);
  }

  context->restart_count += 1;

  out_node_arrays(context, node_array_sol);
  out_node_arrays(context, node_array_sol_increment);
  if (context->max_size[node_array_sol_mf] > 0) {
    out_node_arrays(context, node_array_sol_mf);
  }

  /* No element data needs to be stored :) */

  if (rank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifndef PARALLEL
  fflush(context->out_restart.value_file);
  fflush(context->out_restart.size_file);
#endif

#ifdef PERF
  perf_end(72);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


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
                             INT              step)
{
  BIN_IN_FIELD context;
  MAP* result_info;

#ifdef DEBUG
  dstrc_enter("restart_read_bin_aledyn");
#endif

  /*--------------------------------------------------------------------*/
  /* initialize input of this field */

  init_bin_in_field(&context, sysarray_typ, sysarray, actfield, actpart, actintra, disnum);

  /*--------------------------------------------------------------------*/
  /* find the describtion of the restart data in the control file */

  result_info = in_find_restart_group(&context, disnum, step);

  /*--------------------------------------------------------------------*/
  /* additional values that must be read */

  adyn->step = step;
  adyn->time = map_read_real(result_info, "time");

  /*--------------------------------------------------------------------*/
  /* read the individual node arrays */

  in_node_arrays(&context, result_info, node_array_sol);
  in_node_arrays(&context, result_info, node_array_sol_increment);
  if (map_has_map(result_info, "sol_mf")) {
    in_node_arrays(&context, result_info, node_array_sol_mf);
  }

  /*--------------------------------------------------------------------*/
  /* done with this field */

  destroy_bin_in_field(&context);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write fsi restart data.

  This is additionally to the individual fields data. So only very
  little is stored here. No field at all.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_write_bin_fsidyn(FSI_DYNAMIC *fsidyn)
{
#ifdef DEBUG
  dstrc_enter("restart_write_bin_fsidyn");
#endif

#ifdef PERF
  perf_begin(72);
#endif

  /* It's crucial to make sure everything that's written will go to
   * the restart files (not the result ones). */
  /*out_activate_restart(context);*/

  /*
   * We take the topmost rank here because there is no field to ask
   * for a local rank. Actually the current implementation silently
   * assumes that all partitions are spread over all processors and
   * the rank zero processor is the same in all partitions. Else we
   * would need a control file per partition. */

  if (par.myrank == 0) {
    fprintf(bin_out_main.control_file, "restart:\n"
            "    field = \"fsi\"\n"
            "    step = %d\n"
            "    time = %20.20f\n"
            "    relax = %20.20f\n"
            "\n",
            fsidyn->step, fsidyn->time, fsidyn->relax);
  }

  /* No binary data is written here. Don't count this. */
#if 0
  /* If we've reached the number of steps to put into one file
   * open a new one. */
  if ((context->restart_count % bin_out_main.steps_per_file) == 0)
  {
    out_open_data_files(context, context->out, "restart", context->restart_count);
  }

  context->restart_count += 1;
#endif

  if (par.myrank == 0) {
    fflush(bin_out_main.control_file);
  }

#ifdef PERF
  perf_end(72);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Read fsi restart data.

  \author u.kue
  \date 10/04
*/
/*----------------------------------------------------------------------*/
void restart_read_bin_fsidyn(FSI_DYNAMIC *fsidyn, INT step)
{
  MAP *result_info = NULL;
  SYMBOL *symbol;

#ifdef DEBUG
  dstrc_enter("restart_read_bin_fsidyn");
#endif

  symbol = map_find_symbol(&(bin_in_main.table), "restart");
  while (symbol != NULL) {
    if (symbol_is_map(symbol)) {
      MAP* map;
      symbol_get_map(symbol, &map);
      if (map_has_string(map, "field", "fsi") &&
          map_has_int(map, "step", step)) {
        result_info = map;
        break;
      }
    }
    symbol = symbol->next;
  }
  if (symbol == NULL) {
    dserror(
        "No restart entry for step %d in symbol table. Control file corrupt?", step);
  }

  fsidyn->step  = step;
  fsidyn->time  = map_read_real(result_info, "time");
  fsidyn->relax = map_read_real(result_info, "relax");

#ifdef DEBUG
  dstrc_exit();
#endif
}

/*! @} (documentation module close)*/
#endif
