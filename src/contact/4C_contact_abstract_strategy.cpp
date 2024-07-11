/*---------------------------------------------------------------------*/
/*! \file
\brief Main abstract class for contact solution strategies

\level 2


*/
/*---------------------------------------------------------------------*/
#include "4C_contact_abstract_strategy.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_noxinterface.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_contact_utils_parallel.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_defines.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Epetra_SerialComm.h>
#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AbstractStratDataContainer::AbstractStratDataContainer()
    : glmdofrowmap_(Teuchos::null),
      gsnoderowmap_(Teuchos::null),
      gmnoderowmap_(Teuchos::null),
      gsdofrowmap_(Teuchos::null),
      gmdofrowmap_(Teuchos::null),
      gndofrowmap_(Teuchos::null),
      gsmdofrowmap_(Teuchos::null),
      gdisprowmap_(Teuchos::null),
      gactivenodes_(Teuchos::null),
      gactivedofs_(Teuchos::null),
      ginactivenodes_(Teuchos::null),
      ginactivedofs_(Teuchos::null),
      gactiven_(Teuchos::null),
      gactivet_(Teuchos::null),
      gslipnodes_(Teuchos::null),
      gslipdofs_(Teuchos::null),
      gslipt_(Teuchos::null),
      gsdof_vertex_(Teuchos::null),
      gsdof_edge_(Teuchos::null),
      gsdof_surf_(Teuchos::null),
      unbalance_evaluation_time_(0),
      unbalance_num_slave_elements_(0),
      pglmdofrowmap_(Teuchos::null),
      pgsdofrowmap_(Teuchos::null),
      pgmdofrowmap_(Teuchos::null),
      pgsmdofrowmap_(Teuchos::null),
      pgsdirichtoggle_(Teuchos::null),
      partype_(Inpar::Mortar::ParallelRedist::redist_none),
      initial_elecolmap_(Teuchos::null),
      dmatrix_(Teuchos::null),
      mmatrix_(Teuchos::null),
      wgap_(Teuchos::null),
      tangrhs_(Teuchos::null),
      inactiverhs_(Teuchos::null),
      str_contact_rhs_ptr_(Teuchos::null),
      constrrhs_(Teuchos::null),
      lindmatrix_(Teuchos::null),
      linmmatrix_(Teuchos::null),
      kteffnew_(Teuchos::null),
      dold_(Teuchos::null),
      mold_(Teuchos::null),
      z_(Teuchos::null),
      zold_(Teuchos::null),
      zincr_(Teuchos::null),
      zuzawa_(Teuchos::null),
      stressnormal_(Teuchos::null),
      stresstangential_(Teuchos::null),
      forcenormal_(Teuchos::null),
      forcetangential_(Teuchos::null),
      stepnp_(-1),
      iter_(-1),
      isincontact_(false),
      wasincontact_(false),
      wasincontactlts_(false),
      isselfcontact_(false),
      friction_(false),
      non_smooth_contact_(false),
      regularized_(false),
      dualquadslavetrafo_(false),
      trafo_(Teuchos::null),
      invtrafo_(Teuchos::null),
      dmatrixmod_(Teuchos::null),
      doldmod_(Teuchos::null),
      inttime_(0.0),
      ivel_(0),
      stype_(Inpar::CONTACT::solution_vague),
      constr_direction_(Inpar::CONTACT::constr_vague)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CONTACT::AbstractStrategy::AbstractStrategy(
    const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params_in, const int spatialDim,
    const Teuchos::RCP<const Epetra_Comm>& comm, const double alphaf, const int maxdof)
    : Mortar::StrategyBase(
          data_ptr, dof_row_map, NodeRowMap, params_in, spatialDim, comm, alphaf, maxdof),
      glmdofrowmap_(data_ptr->global_lm_dof_row_map_ptr()),
      gsnoderowmap_(data_ptr->global_slave_node_row_map_ptr()),
      gmnoderowmap_(data_ptr->global_master_node_row_map_ptr()),
      gsdofrowmap_(data_ptr->global_slave_dof_row_map_ptr()),
      gmdofrowmap_(data_ptr->global_master_dof_row_map_ptr()),
      gndofrowmap_(data_ptr->global_internal_dof_row_map_ptr()),
      gsmdofrowmap_(data_ptr->global_slave_master_dof_row_map_ptr()),
      gdisprowmap_(data_ptr->global_disp_dof_row_map_ptr()),
      gactivenodes_(data_ptr->global_active_node_row_map_ptr()),
      gactivedofs_(data_ptr->global_active_dof_row_map_ptr()),
      ginactivenodes_(data_ptr->global_inactive_node_row_map_ptr()),
      ginactivedofs_(data_ptr->global_inactive_dof_row_map_ptr()),
      gactiven_(data_ptr->global_active_n_dof_row_map_ptr()),
      gactivet_(data_ptr->global_active_t_dof_row_map_ptr()),
      gslipnodes_(data_ptr->global_slip_node_row_map_ptr()),
      gslipdofs_(data_ptr->global_slip_dof_row_map_ptr()),
      gslipt_(data_ptr->global_slip_t_dof_row_map_ptr()),
      gsdofVertex_(data_ptr->global_slave_dof_vertex_row_map_ptr()),
      gsdofEdge_(data_ptr->global_slave_dof_edge_row_map_ptr()),
      gsdofSurf_(data_ptr->global_slave_dof_surface_row_map_ptr()),
      unbalanceEvaluationTime_(data_ptr->unbalance_time_factors()),
      unbalanceNumSlaveElements_(data_ptr->unbalance_element_factors()),
      pglmdofrowmap_(data_ptr->pg_lm_dof_row_map_ptr()),
      pgsdofrowmap_(data_ptr->pg_sl_dof_row_map_ptr()),
      pgmdofrowmap_(data_ptr->pg_ma_dof_row_map_ptr()),
      pgsmdofrowmap_(data_ptr->pg_sl_ma_dof_row_map_ptr()),
      pgsdirichtoggle_(data_ptr->pg_sl_dirich_toggle_dof_row_map_ptr()),
      initial_elecolmap_(data_ptr->initial_sl_ma_ele_col_map()),
      dmatrix_(data_ptr->d_matrix_ptr()),
      mmatrix_(data_ptr->m_matrix_ptr()),
      wgap_(data_ptr->w_gap_ptr()),
      tangrhs_(data_ptr->tang_rhs_ptr()),
      inactiverhs_(data_ptr->inactive_rhs_ptr()),
      strcontactrhs_(data_ptr->str_contact_rhs_ptr()),
      constrrhs_(data_ptr->constr_rhs_ptr()),
      lindmatrix_(data_ptr->d_lin_matrix_ptr()),
      linmmatrix_(data_ptr->m_lin_matrix_ptr()),
      kteffnew_(data_ptr->kteffnew_matrix_ptr()),
      dold_(data_ptr->old_d_matrix_ptr()),
      mold_(data_ptr->old_m_matrix_ptr()),
      z_(data_ptr->lm_ptr()),
      zold_(data_ptr->old_lm_ptr()),
      zincr_(data_ptr->lm_incr_ptr()),
      zuzawa_(data_ptr->lm_uzawa_ptr()),
      stressnormal_(data_ptr->stress_normal_ptr()),
      stresstangential_(data_ptr->stress_tangential_ptr()),
      forcenormal_(data_ptr->force_normal_ptr()),
      forcetangential_(data_ptr->force_tangential_ptr()),
      step_(data_ptr->step_np()),
      iter_(data_ptr->nln_iter()),
      isincontact_(data_ptr->is_in_contact()),
      wasincontact_(data_ptr->was_in_contact()),
      wasincontactlts_(data_ptr->was_in_contact_last_time_step()),
      isselfcontact_(data_ptr->is_self_contact()),
      friction_(data_ptr->is_friction()),
      nonSmoothContact_(data_ptr->is_non_smooth_contact()),
      regularized_(data_ptr->is_regularized()),
      dualquadslavetrafo_(data_ptr->is_dual_quad_slave_trafo()),
      trafo_(data_ptr->trafo_ptr()),
      invtrafo_(data_ptr->inv_trafo_ptr()),
      dmatrixmod_(data_ptr->modified_d_matrix_ptr()),
      doldmod_(data_ptr->old_modified_d_matrix_ptr()),
      inttime_(data_ptr->int_time()),
      ivel_(data_ptr->mean_interface_vels()),
      stype_(data_ptr->sol_type()),
      constr_direction_(data_ptr->constr_direction()),
      data_ptr_(data_ptr),
      noxinterface_ptr_(Teuchos::null)
{
  // set data container pointer (only PRIVATE direct access!)
  data_ptr_->sol_type() =
      Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(params_in, "STRATEGY");
  data_ptr_->constr_direction() = Core::UTILS::IntegralValue<Inpar::CONTACT::ConstraintDirection>(
      params_in, "CONSTRAINT_DIRECTIONS");
  data_ptr_->par_type() = Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(
      params_in.sublist("PARALLEL REDISTRIBUTION"), "PARALLEL_REDIST");

  Inpar::CONTACT::FrictionType ftype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::FrictionType>(params(), "FRICTION");

  // set frictional contact status
  if (ftype != Inpar::CONTACT::friction_none) friction_ = true;

  // set nonsmooth contact status
  if (Core::UTILS::IntegralValue<int>(params(), "NONSMOOTH_GEOMETRIES")) nonSmoothContact_ = true;

  if (Core::UTILS::IntegralValue<Inpar::CONTACT::Regularization>(
          params(), "CONTACT_REGULARIZATION") != Inpar::CONTACT::reg_none)
    regularized_ = true;

  // initialize storage fields for parallel redistribution
  unbalanceEvaluationTime_.clear();
  unbalanceNumSlaveElements_.clear();

  // build the NOX::Nln::CONSTRAINT::Interface::Required object
  noxinterface_ptr_ = Teuchos::rcp(new CONTACT::NoxInterface);
  noxinterface_ptr_->init(Teuchos::rcp(this, false));
  noxinterface_ptr_->setup();

  return;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::AbstractStrategy& strategy)
{
  strategy.print(os);
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::is_rebalancing_necessary(const bool first_time_step)
{
  // No rebalancing of a serial run, since it makes no sense.
  if (get_comm().NumProc() == 1) return false;

  bool perform_rebalancing = false;
  const double max_time_unbalance =
      params().sublist("PARALLEL REDISTRIBUTION").get<double>("MAX_BALANCE_EVAL_TIME");
  const double max_ele_unbalance =
      params().sublist("PARALLEL REDISTRIBUTION").get<double>("MAX_BALANCE_SLAVE_ELES");

  double time_average = 0.0;
  double elements_average = 0.0;
  if (!first_time_step)
    compute_and_reset_parallel_balance_indicators(time_average, elements_average);

  switch (which_parallel_redistribution())
  {
    case Inpar::Mortar::ParallelRedist::redist_none:
    {
      break;
    }
    case Inpar::Mortar::ParallelRedist::redist_static:
    {
      // Static redistribution: ONLY at time t=0 or after restart
      if (first_time_step)
      {
        // The user demanded to perform rebalancing, so let's do it.
        perform_rebalancing = true;
      }

      break;
    }
    case Inpar::Mortar::ParallelRedist::redist_dynamic:
    {
      // Dynamic redistribution: whenever system is out of balance

      // This is the first time step (t=0) or restart
      if (first_time_step)
      {
        // Always perform rebalancing in the first time step
        perform_rebalancing = true;
      }

      // This is a regular time step (neither t=0 nor restart)
      else
      {
        /* Decide on redistribution
         *
         * We allow a maximum value of the balance measure in the system as defined in the input
         * parameter MAX_BALANCE_EVAL_TIME, i.e. the maximum local processor workload and the
         * minimum local processor workload for mortar evaluation of all interfaces may not differ
         * by more than (MAX_BALANCE_EVAL_TIME - 1.0)*100%)
         *
         * Moreover, we redistribute if in the majority of iteration steps of the last time step
         * there has been an unbalance in element distribution.
         */
        if (time_average >= max_time_unbalance || elements_average >= max_ele_unbalance)
          perform_rebalancing = true;
      }

      break;
    }
  }

  print_parallel_balance_indicators(time_average, elements_average, max_time_unbalance);

  return perform_rebalancing;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::compute_and_reset_parallel_balance_indicators(
    double& time_average, double& elements_average)
{
  FOUR_C_ASSERT(unbalanceEvaluationTime_.size() > 0, "Vector should have length > 0.");
  FOUR_C_ASSERT(unbalanceNumSlaveElements_.size() > 0, "Vector should have length > 0.");

  // compute average balance factors of last time step
  for (const auto& time : unbalanceEvaluationTime_) time_average += time;
  time_average /= static_cast<double>(unbalanceEvaluationTime_.size());
  for (const auto& num_elements : unbalanceNumSlaveElements_) elements_average += num_elements;
  elements_average /= static_cast<double>(unbalanceNumSlaveElements_.size());

  // Reset balance factors of last time step
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::print_parallel_balance_indicators(
    double& time_average, double& elements_average, const double& max_time_unbalance) const
{
  // Screen output only on proc 0
  if (get_comm().MyPID() == 0)
  {
    std::cout << "*************** DATA OF PREVIOUS TIME STEP ***************" << std::endl;
    if (time_average > 0)
    {
      std::cout << "Parallel balance (time): " << time_average << " (limit " << max_time_unbalance
                << ")\n"
                << "Parallel balance (eles): " << elements_average << " (limit 0.5)" << std::endl;
    }
    else
      std::cout << "Parallel balance: t=0/restart" << std::endl;
    std::cout << "**********************************************************" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::is_update_of_ghosting_necessary(
    const Inpar::Mortar::ExtendGhosting& ghosting_strategy, const bool first_time_step) const
{
  bool enforce_update_of_ghosting = false;
  switch (ghosting_strategy)
  {
    case Inpar::Mortar::ExtendGhosting::redundant_all:
    case Inpar::Mortar::ExtendGhosting::redundant_master:
    {
      // this is the first time step (t=0) or restart
      if (first_time_step)
        enforce_update_of_ghosting = true;
      else
        enforce_update_of_ghosting = false;

      break;
    }
    case Inpar::Mortar::ExtendGhosting::roundrobin:
    case Inpar::Mortar::ExtendGhosting::binning:
    {
      enforce_update_of_ghosting = true;
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown strategy to extend ghosting if necessary.");
    }
  }

  return enforce_update_of_ghosting;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::redistribute_contact(
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel)
{
  bool redistributed = false;

  if (CONTACT::UTILS::UseSafeRedistributeAndGhosting(params()))
    redistributed = redistribute_with_safe_ghosting(*dis, *vel);
  else
  {
    if (get_comm().MyPID() == 0)
    {
      std::cout << "+++++++++++++++++++++++++++++++ WARNING +++++++++++++++++++++++++++++++\n"
                << "+++ You're using an outdated contact redistribution implementation, +++\n"
                << "+++ that might deliver an insufficient master-side ghosting.        +++\n"
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
                << std::endl;
    }
    redistributed = redistribute_contact_old(dis, vel);
  }

  return redistributed;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::redistribute_with_safe_ghosting(
    const Epetra_Vector& displacement, const Epetra_Vector& velocity)
{
  // time measurement
  get_comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  const Inpar::Mortar::ExtendGhosting ghosting_strategy =
      Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
          params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY");

  bool first_time_step = is_first_time_step();
  const bool perform_rebalancing = is_rebalancing_necessary(first_time_step);
  const bool enforce_ghosting_update =
      is_update_of_ghosting_necessary(ghosting_strategy, first_time_step);

  // Prepare for extending the ghosting
  ivel_.resize(interfaces().size(), 0.0);  // initialize to zero for non-binning strategies
  if (ghosting_strategy == Inpar::Mortar::ExtendGhosting::binning)
    calc_mean_velocity_for_binning(velocity);

  // Set old and current displacement state (needed for search within redistribution)
  if (perform_rebalancing)
  {
    set_state(Mortar::state_new_displacement, displacement);
    set_state(Mortar::state_old_displacement, displacement);
  }

  // Update parallel distribution and ghosting of all interfaces
  for (std::size_t i = 0; i < interfaces().size(); ++i)
    interfaces()[i]->update_parallel_layout_and_data_structures(
        perform_rebalancing, enforce_ghosting_update, maxdof_, ivel_[i]);

  // Re-setup strategy to update internal map objects
  if (perform_rebalancing) setup(true, false);

  // time measurement
  get_comm().Barrier();
  const double t_end = Teuchos::Time::wallTime() - t_start;
  if (get_comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution..............." << std::scientific
              << std::setprecision(6) << t_end << " secs\n"
              << std::endl;

  return perform_rebalancing;
}

/*----------------------------------------------------------------------*
 | parallel redistribution                                   popp 09/10 |
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::redistribute_contact_old(
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<const Epetra_Vector> vel)
{
  // decide whether redistribution should be applied or not
  bool first_time_step = is_first_time_step();
  const bool doredist = is_rebalancing_necessary(first_time_step);

  // get out of here if simulation is still in balance
  if (!doredist) return false;

  // time measurement
  get_comm().Barrier();
  const double t_start = Teuchos::Time::wallTime();

  // Prepare for extending the ghosting
  ivel_.resize(interfaces().size(), 0.0);  // initialize to zero for non-binning strategies
  if (Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
          params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY") ==
      Inpar::Mortar::ExtendGhosting::binning)
    calc_mean_velocity_for_binning(*vel);

  /* set old and current displacement state
   * (needed for search within redistribution) */
  set_state(Mortar::state_new_displacement, *dis);
  set_state(Mortar::state_old_displacement, *dis);

  // parallel redistribution of all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // redistribute optimally among procs
    interfaces()[i]->redistribute();

    // call fill complete again
    interfaces()[i]->fill_complete(true, maxdof_, ivel_[i]);

    // print new parallel distribution
    if (get_comm().MyPID() == 0)
      std::cout << "Interface parallel distribution after rebalancing:" << std::endl;
    interfaces()[i]->print_parallel_distribution();

    // re-create binary search tree
    interfaces()[i]->create_search_tree();
  }

  // re-setup strategy with redistributed=TRUE, init=FALSE
  setup(true, false);

  // time measurement
  get_comm().Barrier();
  double t_end = Teuchos::Time::wallTime() - t_start;
  if (get_comm().MyPID() == 0)
    std::cout << "\nTime for parallel redistribution..............." << std::scientific
              << std::setprecision(6) << t_end << " secs\n"
              << std::endl;

  return doredist;
}

/*----------------------------------------------------------------------*
 | setup this strategy object                                popp 08/10 |
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::setup(bool redistributed, bool init)
{
  if (init)
  {
    // set potential global self contact status
    // (this is TRUE if at least one contact interface is a self contact interface)
    bool selfcontact = false;
    for (unsigned i = 0; i < interfaces().size(); ++i)
      if (interfaces()[i]->self_contact()) selfcontact = true;

    if (selfcontact) isselfcontact_ = true;
  }

  // ------------------------------------------------------------------------
  // setup global accessible Epetra_Maps
  // ------------------------------------------------------------------------

  // make sure to remove all existing maps first
  // (do NOT remove map of non-interface dofs after redistribution)
  gsdofrowmap_ = Teuchos::null;
  gmdofrowmap_ = Teuchos::null;
  gsmdofrowmap_ = Teuchos::null;
  glmdofrowmap_ = Teuchos::null;
  gdisprowmap_ = Teuchos::null;
  gsnoderowmap_ = Teuchos::null;
  gmnoderowmap_ = Teuchos::null;
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  ginactivenodes_ = Teuchos::null;
  ginactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  if (!redistributed) gndofrowmap_ = Teuchos::null;
  if (init) initial_elecolmap_.clear();
  initial_elecolmap_.resize(0);

  if (friction_)
  {
    gslipnodes_ = Teuchos::null;
    gslipdofs_ = Teuchos::null;
    gslipt_ = Teuchos::null;
  }

  // initialize vertex, edge and surface maps for nonsmooth case
  if (Core::UTILS::IntegralValue<int>(params(), "NONSMOOTH_GEOMETRIES"))
  {
    gsdofVertex_ = Teuchos::null;
    gsdofEdge_ = Teuchos::null;
    gsdofSurf_ = Teuchos::null;
  }

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // merge interface maps to global maps
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // build Lagrange multiplier dof map
    if (is_self_contact())
    {
      if (redistributed) FOUR_C_THROW("SELF-CONTACT: Parallel redistribution is not supported!");

      Interface& inter = *interfaces()[i];
      Teuchos::RCP<const Epetra_Map> refdofrowmap = Teuchos::null;
      if (inter.self_contact())
        refdofrowmap = Core::LinAlg::MergeMap(inter.slave_row_dofs(), inter.master_row_dofs());
      else
        refdofrowmap = inter.slave_row_dofs();

      Teuchos::RCP<Epetra_Map> selfcontact_lmmap =
          interfaces()[i]->update_lag_mult_sets(offset_if, redistributed, *refdofrowmap);

      Teuchos::RCP<Epetra_Map>& gsc_refdofmap_ptr =
          data().global_self_contact_ref_dof_row_map_ptr();
      Teuchos::RCP<Epetra_Map>& gsc_lmdofmap_ptr = data().global_self_contact_lm_dof_row_map_ptr();
      gsc_lmdofmap_ptr = Core::LinAlg::MergeMap(selfcontact_lmmap, gsc_lmdofmap_ptr);
      gsc_refdofmap_ptr = Core::LinAlg::MergeMap(refdofrowmap, gsc_refdofmap_ptr);

      const int loffset_interface = selfcontact_lmmap->NumGlobalElements();
      if (loffset_interface > 0) offset_if += loffset_interface;
    }
    else
    {
      interfaces()[i]->update_lag_mult_sets(offset_if, redistributed);
      const int loffset_interface = interfaces()[i]->lag_mult_dofs()->NumGlobalElements();
      if (loffset_interface > 0) offset_if += loffset_interface;
    }

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ =
        Core::LinAlg::MergeMap(slave_row_nodes_ptr(), interfaces()[i]->slave_row_nodes());
    gmnoderowmap_ =
        Core::LinAlg::MergeMap(master_row_nodes_ptr(), interfaces()[i]->master_row_nodes());
    gsdofrowmap_ =
        Core::LinAlg::MergeMap(slave_dof_row_map_ptr(true), interfaces()[i]->slave_row_dofs());
    gmdofrowmap_ = Core::LinAlg::MergeMap(gmdofrowmap_, interfaces()[i]->master_row_dofs());

    // merge active sets and slip sets of all interfaces
    // (these maps are NOT allowed to be overlapping !!!)
    interfaces()[i]->build_active_set(init);
    gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interfaces()[i]->active_nodes(), false);
    gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interfaces()[i]->active_dofs(), false);

    ginactivenodes_ =
        Core::LinAlg::MergeMap(ginactivenodes_, interfaces()[i]->in_active_nodes(), false);
    ginactivedofs_ =
        Core::LinAlg::MergeMap(ginactivedofs_, interfaces()[i]->in_active_dofs(), false);

    gactiven_ = Core::LinAlg::MergeMap(gactiven_, interfaces()[i]->active_n_dofs(), false);
    gactivet_ = Core::LinAlg::MergeMap(gactivet_, interfaces()[i]->active_t_dofs(), false);

    // store initial element col map for binning strategy
    initial_elecolmap_.push_back(
        Teuchos::rcp<Epetra_Map>(new Epetra_Map(*interfaces()[i]->discret().element_col_map())));

    // ****************************************************
    // friction
    // ****************************************************
    if (friction_)
    {
      gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interfaces()[i]->slip_nodes(), false);
      gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interfaces()[i]->slip_dofs(), false);
      gslipt_ = Core::LinAlg::MergeMap(gslipt_, interfaces()[i]->slip_t_dofs(), false);
    }

    // define maps for nonsmooth case
    if (Core::UTILS::IntegralValue<int>(params(), "NONSMOOTH_GEOMETRIES"))
    {
      gsdofVertex_ = Core::LinAlg::MergeMap(gsdofVertex_, interfaces()[i]->sdof_vertex_rowmap());
      gsdofEdge_ = Core::LinAlg::MergeMap(gsdofEdge_, interfaces()[i]->sdof_edge_rowmap());
      gsdofSurf_ = Core::LinAlg::MergeMap(gsdofSurf_, interfaces()[i]->sdof_surf_rowmap());
    }
  }

  // create the global Lagrange multiplier DoF row map
  glmdofrowmap_ = create_deterministic_lm_dof_row_map(*gsdofrowmap_);

  // setup global non-slave-or-master dof map
  // (this is done by splitting from the discretization dof map)
  // (no need to rebuild this map after redistribution)
  if (!redistributed)
  {
    gndofrowmap_ = Core::LinAlg::SplitMap(*(problem_dofs()), slave_dof_row_map(true));
    gndofrowmap_ = Core::LinAlg::SplitMap(*gndofrowmap_, *gmdofrowmap_);
  }

  // setup combined global slave and master dof map
  // setup global displacement dof map
  gsmdofrowmap_ = Core::LinAlg::MergeMap(slave_dof_row_map(true), *gmdofrowmap_, false);
  gdisprowmap_ = Core::LinAlg::MergeMap(*gndofrowmap_, *gsmdofrowmap_, false);

  // TODO: check if necessary!
  // due to boundary modification we have to extend master map to slave dofs
  //  if(Core::UTILS::IntegralValue<int>(Params(),"NONSMOOTH_GEOMETRIES"))
  //  {
  //    gmdofrowmap_ = Core::LinAlg::MergeMap(slave_dof_row_map(true), *gmdofrowmap_, false);
  //  }

  // initialize flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // ------------------------------------------------------------------------
  // setup global accessible vectors and matrices
  // ------------------------------------------------------------------------

  // initialize vectors and matrices
  if (!redistributed)
  {
    // setup Lagrange multiplier vectors
    z_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    zincr_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    zold_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

    // setup global mortar matrices Dold and Mold
    dold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 1, true, false));
    dold_->zero();
    dold_->complete();
    mold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 1, true, false));
    mold_->zero();
    mold_->complete(*gmdofrowmap_, slave_dof_row_map(true));
  }

  // In the redistribution case, first check if the vectors and
  // matrices have already been defined, If yes, transform them
  // to the new redistributed maps. If not, initialize them.
  // Moreover, store redistributed quantities into nodes!!!
  else
  {
    // setup Lagrange multiplier vectors
    if (z_ == Teuchos::null)
    {
      z_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newz = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*z_, *newz);
      z_ = newz;
    }

    if (zincr_ == Teuchos::null)
    {
      zincr_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> newzincr =
          Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*zincr_, *newzincr);
      zincr_ = newzincr;
    }

    if (zold_ == Teuchos::null)
      zold_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzold =
          Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*zold_, *newzold);
      zold_ = newzold;
    }

    if (zuzawa_ == Teuchos::null)
      zuzawa_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    else
    {
      Teuchos::RCP<Epetra_Vector> newzuzawa =
          Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
      Core::LinAlg::export_to(*zuzawa_, *newzuzawa);
      zuzawa_ = newzuzawa;
    }

    // setup global Mortar matrices Dold and Mold
    if (dold_ == Teuchos::null)
    {
      dold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 1, true, false));
      dold_->zero();
      dold_->complete();
    }
    else if (dold_->row_map().NumGlobalElements() > 0)
      dold_ = Mortar::matrix_row_col_transform(
          dold_, slave_dof_row_map_ptr(true), slave_dof_row_map_ptr(true));

    if (mold_ == Teuchos::null)
    {
      mold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 1, true, false));
      mold_->zero();
      mold_->complete(*gmdofrowmap_, slave_dof_row_map(true));
    }
    else if (mold_->row_map().NumGlobalElements() > 0)
      mold_ = Mortar::matrix_row_col_transform(mold_, slave_dof_row_map_ptr(true), gmdofrowmap_);
  }

  // output contact stress vectors
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // These matrices need to be applied to the slave displacements
  // in the cases of dual LM interpolation for tet10/hex20 meshes
  // in 3D or for locally linear Lagrange multipliers for line3 meshes
  // in 2D. Here, the displacement basis functions have been modified
  // in order to assure positivity of the D matrix entries and at
  // the same time biorthogonality. Thus, to scale back the modified
  // discrete displacements \hat{d} to the nodal discrete displacements
  // {d}, we have to apply the transformation matrix T and vice versa
  // with the transformation matrix T^(-1).
  //----------------------------------------------------------------------
  Inpar::Mortar::ShapeFcn shapefcn =
      Core::UTILS::IntegralValue<Inpar::Mortar::ShapeFcn>(params(), "LM_SHAPEFCN");
  Inpar::Mortar::LagMultQuad lagmultquad =
      Core::UTILS::IntegralValue<Inpar::Mortar::LagMultQuad>(params(), "LM_QUAD");
  if ((shapefcn == Inpar::Mortar::shape_dual || shapefcn == Inpar::Mortar::shape_petrovgalerkin) &&
      (n_dim() == 3 || (n_dim() == 2 && lagmultquad == Inpar::Mortar::lagmult_lin)))
    for (int i = 0; i < (int)interfaces().size(); ++i)
      dualquadslavetrafo_ += interfaces()[i]->quadslave();

  //----------------------------------------------------------------------
  // IF SO, COMPUTE TRAFO MATRIX AND ITS INVERSE
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // for locally linear Lagrange multipliers, consider both slave and master DOFs,
    // and otherwise, only consider slave DOFs
    if (lagmultquad == Inpar::Mortar::lagmult_lin)
    {
      trafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsmdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsmdofrowmap_, 10));
    }
    else
    {
      trafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
      invtrafo_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*gsdofrowmap_, 10));
    }

    // set of already processed nodes
    // (in order to avoid double-assembly for N interfaces)
    std::set<int> donebefore;

    // for all interfaces
    for (int i = 0; i < (int)interfaces().size(); ++i)
      interfaces()[i]->assemble_trafo(*trafo_, *invtrafo_, donebefore);

    // fill_complete() transformation matrices
    trafo_->complete();
    invtrafo_->complete();
  }

  // transform modified old D-matrix in case of friction
  // (ony necessary after parallel redistribution)
  if (redistributed && friction_ && is_dual_quad_slave_trafo())
  {
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 1, true, false));
      doldmod_->zero();
      doldmod_->complete();
    }
    else
      doldmod_ = Mortar::matrix_row_col_transform(
          doldmod_, slave_dof_row_map_ptr(true), slave_dof_row_map_ptr(true));
  }

  if (init)
  {
    // store interface maps with parallel distribution of underlying
    // problem discretization (i.e. interface maps before parallel
    // redistribution of slave and master sides)
    if (parallel_redistribution_status())
    {
      for (std::size_t i = 0; i < interfaces().size(); ++i)
        interfaces()[i]->store_unredistributed_maps();
      if (lm_dof_row_map_ptr(true) != Teuchos::null)
        pglmdofrowmap_ = Teuchos::rcp(new Epetra_Map(lm_dof_row_map(true)));
      pgsdofrowmap_ = Teuchos::rcp(new Epetra_Map(slave_dof_row_map(true)));
      pgmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gmdofrowmap_));
      pgsmdofrowmap_ = Teuchos::rcp(new Epetra_Map(*gsmdofrowmap_));
    }
  }

  post_setup(redistributed, init);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CONTACT::AbstractStrategy::create_deterministic_lm_dof_row_map(
    const Epetra_Map& gsdofrowmap) const
{
  const unsigned num_my_sdofs = gsdofrowmap.NumMyElements();
  const int* my_sdof_gids = gsdofrowmap.MyGlobalElements();

  std::vector<int> my_lm_gids(num_my_sdofs, -1);

  for (unsigned slid = 0; slid < num_my_sdofs; ++slid)
  {
    const int sgid = my_sdof_gids[slid];

    // find slid of the interface map
    unsigned interface_id = 0;
    int interface_slid = -1;
    for (auto cit = interfaces().begin(); cit != interfaces().end(); ++cit, ++interface_id)
    {
      const Interface& interface = **cit;
      Teuchos::RCP<const Epetra_Map> sdof_map = interface.slave_row_dofs();

      interface_slid = sdof_map->LID(sgid);
      if (interface_slid != -1) break;
    }

    if (interface_slid == -1)
      FOUR_C_THROW(
          "Couldn't find the global slave dof id #%d in the local interface "
          "maps on proc #%d!",
          sgid, get_comm().MyPID());

    // get the corresponding Lagrange Multiplier GID
    const int interface_lmgid = interfaces()[interface_id]->lag_mult_dofs()->GID(interface_slid);
    if (interface_lmgid == -1)
      FOUR_C_THROW(
          "Couldn't find the corresponding Lagrange multiplier GID! "
          "Note that the UpdateLagMultSets() must be called on each interface "
          "beforehand.");

    my_lm_gids[slid] = interface_lmgid;
  }
  return Teuchos::rcp(
      new Epetra_Map(-1, static_cast<int>(my_lm_gids.size()), my_lm_gids.data(), 0, get_comm()));
}


/*----------------------------------------------------------------------*
 | global evaluation method called from time integrator      popp 06/09 |
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::apply_force_stiff_cmt(Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f,
    const int timeStep, const int nonlinearIteration, bool predictor)
{
  // update step and iteration counters
  step_ = timeStep;
  iter_ = nonlinearIteration;

  // Create timing reports?
  bool doAccurateTimeMeasurements =
      Core::UTILS::IntegralValue<bool>(data().s_contact(), "TIMING_DETAILS");

  if (doAccurateTimeMeasurements)
  {
    // mortar initialization and evaluation
    get_comm().Barrier();
    const double t_start1 = Teuchos::Time::wallTime();
    set_state(Mortar::state_new_displacement, *dis);
    get_comm().Barrier();
    const double t_end1 = Teuchos::Time::wallTime() - t_start1;

    get_comm().Barrier();
    const double t_start2 = Teuchos::Time::wallTime();
    //---------------------------------------------------------------
    // For selfcontact the master/slave sets are updated within the -
    // contact search, see SelfBinaryTree.                          -
    // Therefore, we have to initialize the mortar matrices after   -
    // interface evaluations.                                       -
    //---------------------------------------------------------------
    if (is_self_contact())
    {
      initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
      initialize_mortar();                  // initialize mortar matrices and vectors
      assemble_mortar();                    // assemble mortar terms into global matrices
    }
    else
    {
      initialize_mortar();                  // initialize mortar matrices and vectors
      initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
      assemble_mortar();                    // assemble mortar terms into global matrices
    }
    get_comm().Barrier();
    const double t_end2 = Teuchos::Time::wallTime() - t_start2;

    // evaluate relative movement for friction
    get_comm().Barrier();
    const double t_start3 = Teuchos::Time::wallTime();
    if (predictor)
      predict_relative_movement();
    else
      evaluate_relative_movement();

    // update active set
    if (!predictor) update_active_set_semi_smooth();

    get_comm().Barrier();
    const double t_end3 = Teuchos::Time::wallTime() - t_start3;

    // apply contact forces and stiffness
    get_comm().Barrier();
    const double t_start4 = Teuchos::Time::wallTime();
    initialize();           // init lin-matrices
    evaluate(kt, f, dis);   // assemble lin. matrices, condensation ...
    evaluate_constr_rhs();  // evaluate the constraint rhs (saddle-point system only)

    get_comm().Barrier();
    const double t_end4 = Teuchos::Time::wallTime() - t_start4;

    // only for debugging:
    interface_forces();

    if (get_comm().MyPID() == 0)
    {
      std::cout << "    -->setstate :\t" << t_end1 << " seconds\n";
      std::cout << "    -->interface eval. :\t" << t_end2 << " seconds\n";
      std::cout << "    -->update active set :\t" << t_end3 << " seconds\n";
      std::cout << "    -->modify global system :\t" << t_end4 << " seconds\n";
    }
  }
  else
  {
    // mortar initialization and evaluation
    set_state(Mortar::state_new_displacement, *dis);

    //---------------------------------------------------------------
    // For selfcontact the master/slave sets are updated within the -
    // contact search, see SelfBinaryTree.                          -
    // Therefore, we have to initialize the mortar matrices after   -
    // interface evaluations.                                       -
    //---------------------------------------------------------------
    if (is_self_contact())
    {
      initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
      initialize_mortar();                  // initialize mortar matrices and vectors
      assemble_mortar();                    // assemble mortar terms into global matrices
    }
    else
    {
      initialize_mortar();                  // initialize mortar matrices and vectors
      initialize_and_evaluate_interface();  // evaluate mortar terms (integrate...)
      assemble_mortar();                    // assemble mortar terms into global matrices
    }

    // evaluate relative movement for friction
    if (predictor)
      predict_relative_movement();
    else
      evaluate_relative_movement();

    // update active set
    if (!predictor) update_active_set_semi_smooth();

    // apply contact forces and stiffness
    initialize();           // init lin-matrices
    evaluate(kt, f, dis);   // assemble lin. matrices, condensation ...
    evaluate_constr_rhs();  // evaluate the constraint rhs (saddle-point system only)

    // only for debugging:
    interface_forces();
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::set_state(
    const enum Mortar::StateType& statetype, const Epetra_Vector& vec)
{
  switch (statetype)
  {
    case Mortar::state_new_displacement:
    case Mortar::state_old_displacement:
    {
      // set state on interfaces
      for (int i = 0; i < (int)interfaces().size(); ++i) interfaces()[i]->set_state(statetype, vec);
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Unsupported state type! (state type = %s)", Mortar::StateType2String(statetype).c_str());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | update global master and slave sets (public)               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::update_global_self_contact_state()
{
  if (not is_self_contact()) return;

  // reset global slave / master Epetra Maps
  gsnoderowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
  gsdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
  gmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));
  glmdofrowmap_ = Teuchos::rcp(new Epetra_Map(0, 0, get_comm()));

  // make numbering of LM dofs consecutive and unique across N interfaces
  int offset_if = 0;

  // setup global slave / master Epetra_Maps
  // (this is done by looping over all interfaces and merging)
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // build Lagrange multiplier dof map
    interfaces()[i]->update_self_contact_lag_mult_set(global_self_contact_lm_map(), *gsmdofrowmap_);

    // merge interface Lagrange multiplier dof maps to global LM dof map
    glmdofrowmap_ =
        Core::LinAlg::MergeMap(lm_dof_row_map_ptr(true), interfaces()[i]->lag_mult_dofs());
    offset_if = lm_dof_row_map(true).NumGlobalElements();
    if (offset_if < 0) offset_if = 0;

    // merge interface master, slave maps to global master, slave map
    gsnoderowmap_ =
        Core::LinAlg::MergeMap(slave_row_nodes_ptr(), interfaces()[i]->slave_row_nodes());
    gsdofrowmap_ =
        Core::LinAlg::MergeMap(slave_dof_row_map_ptr(true), interfaces()[i]->slave_row_dofs());
    gmdofrowmap_ = Core::LinAlg::MergeMap(gmdofrowmap_, interfaces()[i]->master_row_dofs());
  }

  Teuchos::RCP<Epetra_Vector> tmp_ptr = Teuchos::rcp(new Epetra_Vector(*gsdofrowmap_, true));

  {
    const int* oldgids = zincr_->Map().MyGlobalElements();
    for (int i = 0; i < zincr_->Map().NumMyElements(); ++i)
    {
      if (std::abs((*zincr_)[i]) > std::numeric_limits<double>::epsilon())
      {
        const int new_lid = gsdofrowmap_->LID(oldgids[i]);
        if (new_lid == -1)
          FOUR_C_THROW(
              "Self contact: The Lagrange multiplier increment vector "
              "could not be transferred consistently.");
        else
          (*tmp_ptr)[new_lid] = (*zincr_)[i];
      }
    }
    zincr_ = Teuchos::rcp(new Epetra_Vector(*tmp_ptr));
  }

  tmp_ptr->PutScalar(0.0);
  {
    const int* oldgids = z_->Map().MyGlobalElements();
    for (int i = 0; i < z_->Map().NumMyElements(); ++i)
    {
      if (std::abs((*z_)[i]) > std::numeric_limits<double>::epsilon())
      {
        const int new_lid = gsdofrowmap_->LID(oldgids[i]);
        if (new_lid == -1)
          FOUR_C_THROW(
              "Self contact: The Lagrange multiplier vector "
              "could not be transferred consistently.");
        else
          (*tmp_ptr)[new_lid] = (*z_)[i];
      }
    }
    z_ = tmp_ptr;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::calc_mean_velocity_for_binning(const Epetra_Vector& velocity)
{
  ivel_.clear();
  ivel_.resize(0);

  // create vector of interface velocities
  for (const auto& interface : interfaces())
  {
    Teuchos::RCP<Epetra_Vector> interfaceVelocity =
        Teuchos::rcp(new Epetra_Vector(*interface->discret().dof_row_map()));
    Core::LinAlg::export_to(velocity, *interfaceVelocity);

    double meanVelocity = 0.0;

    int err = interfaceVelocity->MeanValue(&meanVelocity);
    if (err)
      FOUR_C_THROW("Calculation of mean velocity for interface %s failed.",
          interface->discret().name().c_str());
    meanVelocity = abs(meanVelocity);

    ivel_.push_back(meanVelocity);
  }
  return;
}

/*----------------------------------------------------------------------*
 | initialize + evaluate interface for next Newton step       popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::initialize_and_evaluate_interface(
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // time measurement (on each processor)
  const double t_start = Teuchos::Time::wallTime();

  // get type of parallel strategy
  const Teuchos::ParameterList& mortarParallelRedistParams =
      params().sublist("PARALLEL REDISTRIBUTION");
  Inpar::Mortar::ExtendGhosting extendghosting =
      Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
          mortarParallelRedistParams, "GHOSTING_STRATEGY");

  // Evaluation for all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // initialize / reset interfaces
    interfaces()[i]->initialize();

    // store required integration time
    inttime_ += interfaces()[i]->inttime();

    switch (extendghosting)
    {
      case Inpar::Mortar::ExtendGhosting::roundrobin:
      {
        // first perform rrloop to detect the required ghosting
        interfaces()[i]->round_robin_detect_ghosting();

        // second step --> evaluate
        interfaces()[i]->evaluate(0, step_, iter_);
        break;
      }
      case Inpar::Mortar::ExtendGhosting::binning:
      {
        // required master elements are already ghosted (preparestepcontact) !!!
        // call evaluation
        interfaces()[i]->evaluate(0, step_, iter_);
        break;
      }
      case Inpar::Mortar::ExtendGhosting::redundant_all:
      case Inpar::Mortar::ExtendGhosting::redundant_master:
      {
        interfaces()[i]->evaluate(0, step_, iter_);
        break;
      }
    }
  }  // end interface loop

  // check the parallel distribution
  check_parallel_distribution(t_start);

  //**********************************************************************
  // OVERVIEW OF PARALLEL MORTAR COUPLING STATUS
  //**********************************************************************
#ifdef CONTACTSTATUS
  // total numbers per processor
  std::vector<int> smpairs(1);
  std::vector<int> smintpairs(1);
  std::vector<int> intcells(1);

  // add numbers of all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    smpairs[0] += Interfaces()[i]->SlaveMasterPairs();
    smintpairs[0] += Interfaces()[i]->SlaveMasterIntPairs();
    intcells[0] += Interfaces()[i]->IntegrationCells();
  }

  // vector containing all proc ids
  const int numproc = Comm().NumProc();
  std::vector<int> allproc(numproc);
  for (int i = 0; i < numproc; ++i) allproc[i] = i;

  // global numbers
  std::vector<int> gsmpairs, gsmintpairs, gintcells;
  Core::LinAlg::Gather<int>(smpairs, gsmpairs, numproc, allproc.data(), Comm());
  Core::LinAlg::Gather<int>(smintpairs, gsmintpairs, numproc, allproc.data(), Comm());
  Core::LinAlg::Gather<int>(intcells, gintcells, numproc, allproc.data(), Comm());

  // output to screen
  if (Comm().MyPID() == 0)
  {
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
    std::cout << std::setw(10) << "proc ID" << std::setw(16) << "# s/m pairs" << std::setw(16)
              << "# s/m intpairs" << std::setw(16) << "# intcells" << std::endl;
    for (int i = 0; i < numproc; ++i)
    {
      std::cout << std::setw(10) << i << std::setw(16) << gsmpairs[i] << std::setw(16)
                << gsmintpairs[i] << std::setw(16) << gintcells[i] << std::endl;
    }
    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
  }
#endif
  //**********************************************************************

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::check_parallel_distribution(const double& t_start)
{
  const double my_total_time = Teuchos::Time::wallTime() - t_start;
  update_parallel_distribution_status(my_total_time);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::update_parallel_distribution_status(const double& my_total_time)
{
  //**********************************************************************
  // PARALLEL REDISTRIBUTION
  //**********************************************************************
  // don't do this if this is a single processor (serial) job
  if (get_comm().NumProc() == 1) return;

  // collect information about participation in coupling evaluation
  // and in parallel distribution of the individual interfaces
  std::vector<int> numloadele((int)interfaces().size());
  std::vector<int> numcrowele((int)interfaces().size());
  for (int i = 0; i < (int)interfaces().size(); ++i)
    interfaces()[i]->collect_distribution_data(numloadele[i], numcrowele[i]);

  // time measurement (on each processor)
  double t_end_for_minall = my_total_time;
  double t_end_for_maxall = my_total_time;

  // restrict time measurement to procs that own at least some part
  // of the "close" slave interface section(s) on the global level,
  // i.e. restrict to procs that actually have to do some work
  int gnumloadele = 0;
  for (int i = 0; i < (int)numloadele.size(); ++i) gnumloadele += numloadele[i];

  // for non-loaded procs, set time measurement to values 0.0 / 1.0e12,
  // which do not affect the maximum and minimum identification
  if (gnumloadele == 0)
  {
    t_end_for_minall = 1.0e12;
    t_end_for_maxall = 0.0;
  }

  // store time indicator for parallel redistribution
  // (indicator is the maximum local processor time
  // divided by the minimum local processor time)
  double maxall = 0.0;
  double minall = 0.0;
  get_comm().MaxAll(&t_end_for_maxall, &maxall, 1);
  get_comm().MinAll(&t_end_for_minall, &minall, 1);

  // check for plausibility before storing
  if (maxall == 0.0 && minall == 1.0e12)
    data().unbalance_time_factors().push_back(1.0);
  else
    data().unbalance_time_factors().push_back(maxall / minall);

  // obtain info whether there is an unbalance in element distribution
  bool eleunbalance = false;
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // find out how many close slave elements in total
    int totrowele = 0;
    get_comm().SumAll(&numcrowele[i], &totrowele, 1);

    // find out how many procs have work on this interface
    int lhascrowele = 0;
    int ghascrowele = 0;
    if (numcrowele[i] > 0) lhascrowele = 1;
    get_comm().SumAll(&lhascrowele, &ghascrowele, 1);

    // minimum number of elements per proc
    int minele = params().sublist("PARALLEL REDISTRIBUTION").get<int>("MIN_ELEPROC");
    int numproc = get_comm().NumProc();

    //--------------------------------------------------------------------
    // check if there is an element unbalance
    //--------------------------------------------------------------------
    // CASE 0: if minimum number of elements per proc is zero, but
    // further procs are still available and more than numproc elements
    if ((minele == 0) && (totrowele > numproc) && (ghascrowele < numproc)) eleunbalance = true;

    // CASE 1: in total too few close slave elements but more than one
    // proc is active (otherwise, i.e. if interface small, we have no choice)
    if ((minele > 0) && (totrowele < ghascrowele * minele) && (ghascrowele > 1))
      eleunbalance = true;

    // CASE 2: in total too many close slave elements, but further procs
    // are still available for redsitribution
    if ((minele > 0) && (totrowele >= (ghascrowele + 1) * minele) && (ghascrowele < numproc))
      eleunbalance = true;
  }

  // obtain global info on element unbalance
  int geleunbalance = 0;
  int leleunbalance = (int)(eleunbalance);
  get_comm().SumAll(&leleunbalance, &geleunbalance, 1);
  if (geleunbalance > 0)
    data().unbalance_element_factors().push_back(1);
  else
    data().unbalance_element_factors().push_back(0);

  // debugging output
  // std::cout << "PROC: " << Comm().MyPID() << "\t LOADELE: " << numloadele[0] << "\t ROWELE: " <<
  // numcrowele[0]
  //     << "\t MIN: " << minall << "\t MAX: " << maxall
  //     << "\t tmin: " << t_end_for_minall << "\t tmax: " << t_end_for_maxall
  //     << "\t TUNBALANCE: " << unbalanceEvaluationTime_[(int)unbalanceEvaluationTime_.size()-1]
  //     << "\t EUNBALANCE: " <<
  //     unbalanceNumSlaveElements_[(int)unbalanceNumSlaveElements_.size()-1] << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | initialize mortar stuff for next Newton step               popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::initialize_mortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing D,M etc.
  update_global_self_contact_state();

  // initialize Dold and Mold if not done already
  if (dold_ == Teuchos::null)
  {
    dold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 10));
    dold_->zero();
    dold_->complete();
  }
  if (mold_ == Teuchos::null)
  {
    mold_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100));
    mold_->zero();
    mold_->complete(*gmdofrowmap_, slave_dof_row_map(true));
  }

  // (re)setup global Mortar Core::LinAlg::SparseMatrices and Epetra_Vectors
  dmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 10));
  mmatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100));

  if (constr_direction_ == Inpar::CONTACT::constr_xyz)
    wgap_ = Core::LinAlg::CreateVector(slave_dof_row_map(true), true);
  else if (constr_direction_ == Inpar::CONTACT::constr_ntt)
    wgap_ = Core::LinAlg::CreateVector(slave_row_nodes(), true);
  else
    FOUR_C_THROW("unknown contact constraint direction");

  // in the case of frictional dual quad 3D, also the modified D matrices are setup
  if (friction_ && is_dual_quad_slave_trafo())
  {
    // initialize Dold and Mold if not done already
    if (doldmod_ == Teuchos::null)
    {
      doldmod_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 10));
      doldmod_->zero();
      doldmod_->complete();
    }
    // setup of dmatrixmod_
    dmatrixmod_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 10));
  }

  return;
}

/*----------------------------------------------------------------------*
 | Assemble mortar stuff for next Newton step                 popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::assemble_mortar()
{
  // for all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // assemble D-, M-matrix and g-vector, store them globally
    interfaces()[i]->assemble_dm(*dmatrix_, *mmatrix_);
    interfaces()[i]->assemble_g(*wgap_);

#ifdef CONTACTFDNORMAL
    // FD check of normal derivatives
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
    //    Interfaces()[i]->fd_check_normal_deriv();
    Interfaces()[i]->fd_check_normal_cpp_deriv();
    std::cout << " -- CONTACTFDNORMAL- -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDNORMAL
#ifdef CONTACTFDMORTARD
    // FD check of Mortar matrix D derivatives
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
    dmatrix_->Complete();
    if (dmatrix_->NormOne()) Interfaces()[i]->FDCheckMortarDDeriv();
    dmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARD -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARD
#ifdef CONTACTFDMORTARM
    // FD check of Mortar matrix M derivatives
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
    mmatrix_->Complete(*gmdofrowmap_, *gsdofrowmap_);
    if (mmatrix_->NormOne()) Interfaces()[i]->FDCheckMortarMDeriv();
    mmatrix_->UnComplete();
    std::cout << " -- CONTACTFDMORTARM -----------------------------------" << std::endl;
#endif  // #ifdef CONTACTFDMORTARM
  }

  // fill_complete() global Mortar matrices
  dmatrix_->complete();
  mmatrix_->complete(*gmdofrowmap_, slave_dof_row_map(true));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_reference_state()
{
  // flag for initialization of contact with nodal gaps
  bool initcontactbygap = Core::UTILS::IntegralValue<int>(params(), "INITCONTACTBYGAP");

  // only do something for frictional case
  // or for initialization of initial contact set with nodal gap
  if (!friction_ and !initcontactbygap) return;

  // do mortar calculation
  initialize_mortar();
  initialize_and_evaluate_interface();
  assemble_mortar();

  // (1) GAP INITIALIZATION CASE
  // initialize init contact with nodal gap
  if (initcontactbygap)
  {
    // merge interface maps to global maps
    for (const auto& interface : interfaces())
    {
      // merge active sets and slip sets of all interfaces
      // (these maps are NOT allowed to be overlapping !!!)
      interface->build_active_set(true);
      gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interface->active_nodes(), false);
      gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interface->active_dofs(), false);
      gactiven_ = Core::LinAlg::MergeMap(gactiven_, interface->active_n_dofs(), false);
      gactivet_ = Core::LinAlg::MergeMap(gactivet_, interface->active_t_dofs(), false);

      if (friction_)
      {
        gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interface->slip_nodes(), false);
        gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interface->slip_dofs(), false);
        gslipt_ = Core::LinAlg::MergeMap(gslipt_, interface->slip_t_dofs(), false);
      }
    }

    // initialize flags for global contact status
    if (gactivenodes_->NumGlobalElements())
    {
      isincontact_ = true;
      wasincontact_ = true;
      wasincontactlts_ = true;
    }

    // error if no nodes are initialized to active
    if (gactivenodes_->NumGlobalElements() == 0)
      FOUR_C_THROW("No active nodes: Choose bigger value for INITCONTACTGAPVALUE!");
  }

  // (2) FRICTIONAL CONTACT CASE
  // do some friction stuff
  if (friction_)
  {
    // store contact state to contact nodes (active or inactive)
    store_nodal_quantities(Mortar::StrategyBase::activeold);

    // store D and M to old ones
    store_dm("old");

    // store nodal entries from D and M to old ones
    store_to_old(Mortar::StrategyBase::dm);

    // store nodal normals
    store_to_old(Mortar::StrategyBase::n_old);

    // transform dold_ in the case of dual quadratic 3d
    if (is_dual_quad_slave_trafo())
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix> tempold =
          Core::LinAlg::MLMultiply(*dold_, false, *invtrafo_, false, false, false, true);
      doldmod_ = tempold;
    }

    // evaluate relative movement
    // needed because it is not called in the predictor of the
    // lagrange multiplier strategy
    evaluate_relative_movement();
  }

  // reset unbalance factors for redistribution
  // (since the interface has been evaluated once above)
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate relative movement of contact bodies           gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_relative_movement()
{
  // only for fricional contact
  if (!friction_) return;

  // transformation of slave displacement dofs
  // Dmod       ---->   D * T^(-1)
  if (is_dual_quad_slave_trafo())
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrixmod_ = temp;
  }

  // vector of slave coordinates xs
  Teuchos::RCP<Epetra_Vector> xsmod = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

  for (int i = 0; i < (int)interfaces().size(); ++i) interfaces()[i]->assemble_slave_coord(xsmod);

  // in case of 3D dual quadratic case, slave coordinates xs are modified
  if (is_dual_quad_slave_trafo()) invtrafo_->Apply(*xsmod, *xsmod);

  // ATTENTION: for evaluate_relative_movement() we need the vector xsmod in
  // fully overlapping layout. Thus, export here. First, allreduce
  // slave dof row map to obtain fully overlapping slave dof map.
  Teuchos::RCP<Epetra_Map> fullsdofs = Core::LinAlg::AllreduceEMap(slave_dof_row_map(true));
  Teuchos::RCP<Epetra_Vector> xsmodfull = Teuchos::rcp(new Epetra_Vector(*fullsdofs));
  Core::LinAlg::export_to(*xsmod, *xsmodfull);
  xsmod = xsmodfull;

  // evaluation of obj. invariant slip increment
  // do the evaluation on the interface
  // loop over all slave row nodes on the current interface
  if (Core::UTILS::IntegralValue<int>(params(), "GP_SLIP_INCR") == false)
    for (int i = 0; i < (int)interfaces().size(); ++i)
      interfaces()[i]->evaluate_relative_movement(xsmod, dmatrixmod_, doldmod_);

  return;
}

/*----------------------------------------------------------------------*
 | call appropriate evaluate for contact evaluation           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
    Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis)
{
  // treat frictional and frictionless cases differently
  if (friction_)
    evaluate_friction(kteff, feff);
  else
    evaluate_contact(kteff, feff);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate matrix of normals (for VelocityUpdate)            popp 10/11|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::AbstractStrategy::evaluate_normals(
    Teuchos::RCP<Epetra_Vector> dis)
{
  // set displacement state and evaluate nodal normals
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    interfaces()[i]->set_state(Mortar::state_new_displacement, *dis);
    interfaces()[i]->evaluate_nodal_normals();
  }

  // create empty global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> normals =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_row_nodes(), 3));

  // assemble nodal normals
  for (int i = 0; i < (int)interfaces().size(); ++i) interfaces()[i]->assemble_normals(*normals);

  // complete global matrix
  // (rectangular: rows=snodes, cols=sdofs)
  normals->complete(slave_dof_row_map(true), slave_row_nodes());

  return normals;
}

/*----------------------------------------------------------------------*
 |  Store Lagrange mulitpliers and disp. jumps into CNode     popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::store_nodal_quantities(Mortar::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // get global quantity to be stored in nodes
    Teuchos::RCP<Epetra_Vector> vectorglobal = Teuchos::null;

    // start type switch
    switch (type)
    {
      case Mortar::StrategyBase::lmold:
      {
        vectorglobal = lagrange_multiplier_old();
        break;
      }
      case Mortar::StrategyBase::lmcurrent:
      case Mortar::StrategyBase::lmupdate:
      {
        vectorglobal = lagrange_multiplier();
        break;
      }
      case Mortar::StrategyBase::lmuzawa:
      {
        vectorglobal = lagrange_multiplier_uzawa();
        break;
      }
      case Mortar::StrategyBase::activeold:
      case Mortar::StrategyBase::slipold:
      {
        break;
      }
      default:
        FOUR_C_THROW("store_nodal_quantities: Unknown state std::string variable!");
        break;
    }  // switch

    // slave dof and node map of the interface
    // columnmap for current or updated LM
    // rowmap for remaining cases
    Teuchos::RCP<Epetra_Map> sdofmap, snodemap;
    if (type == Mortar::StrategyBase::lmupdate or type == Mortar::StrategyBase::lmcurrent)
    {
      sdofmap = interfaces()[i]->slave_col_dofs();
      snodemap = interfaces()[i]->slave_col_nodes();
    }
    else
    {
      sdofmap = interfaces()[i]->slave_row_dofs();
      snodemap = interfaces()[i]->slave_row_nodes();
    }

    // export global quantity to current interface slave dof map (column or row)
    Teuchos::RCP<Epetra_Vector> vectorinterface = Teuchos::null;
    vectorinterface = Teuchos::rcp(new Epetra_Vector(*sdofmap));
    if (vectorglobal != Teuchos::null)  // necessary for case "activeold" and wear
      Core::LinAlg::export_to(*vectorglobal, *vectorinterface);

    // loop over all slave nodes (column or row) on the current interface
    for (int j = 0; j < snodemap->NumMyElements(); ++j)
    {
      int gid = snodemap->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // be aware of problem dimension
      const int numdof = cnode->num_dof();
      if (n_dim() != numdof) FOUR_C_THROW("Inconsisteny Dim <-> NumDof");

      // find indices for DOFs of current node in Epetra_Vector
      // and extract this node's quantity from vectorinterface
      std::vector<int> locindex(n_dim());

      for (int dof = 0; dof < n_dim(); ++dof)
      {
        locindex[dof] = (vectorinterface->Map()).LID(cnode->dofs()[dof]);
        if (locindex[dof] < 0) FOUR_C_THROW("StoreNodalQuantites: Did not find dof in map");

        switch (type)
        {
          case Mortar::StrategyBase::lmcurrent:
          {
            cnode->mo_data().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case Mortar::StrategyBase::lmold:
          {
            cnode->mo_data().lmold()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case Mortar::StrategyBase::lmuzawa:
          {
            cnode->mo_data().lmuzawa()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case Mortar::StrategyBase::lmupdate:
          {
#ifndef CONTACTPSEUDO2D
            // throw a FOUR_C_THROW if node is Active and DBC
            if (cnode->is_dbc() && cnode->active())
              FOUR_C_THROW("Slave node %i is active AND carries D.B.C.s!", cnode->id());
#endif  // #ifndef CONTACTPSEUDO2D

            // store updated LM into node
            cnode->mo_data().lm()[dof] = (*vectorinterface)[locindex[dof]];
            break;
          }
          case Mortar::StrategyBase::activeold:
          {
            cnode->data().active_old() = cnode->active();
            break;
          }
          case Mortar::StrategyBase::slipold:
          {
            if (!friction_) FOUR_C_THROW("Slip just for friction problems!");

            FriNode* fnode = dynamic_cast<FriNode*>(cnode);
            fnode->fri_data().slip_old() = fnode->fri_data().slip();
            break;
          }
          default:
            FOUR_C_THROW("store_nodal_quantities: Unknown state std::string variable!");
            break;
        }  // switch
      }
    }  // end slave loop
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Output vector of normal/tang. contact stresses        gitterle 08/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::compute_contact_stresses()
{
  // reset contact stress class variables
  stressnormal_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  stresstangential_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // be aware of problem dimension
      const int numdof = cnode->num_dof();
      if (n_dim() != numdof) FOUR_C_THROW("Inconsisteny Dim <-> NumDof");

      double nn[3];
      double nt1[3];
      double nt2[3];
      double lmn = 0.0;
      double lmt1 = 0.0;
      double lmt2 = 0.0;

      for (int j = 0; j < 3; ++j)
      {
        nn[j] = cnode->mo_data().n()[j];
        nt1[j] = cnode->data().txi()[j];
        nt2[j] = cnode->data().teta()[j];
        lmn += nn[j] * cnode->mo_data().lm()[j];
        lmt1 += nt1[j] * cnode->mo_data().lm()[j];
        lmt2 += nt2[j] * cnode->mo_data().lm()[j];
      }

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(n_dim());

      // normal stress components
      for (int dof = 0; dof < n_dim(); ++dof)
      {
        locindex[dof] = (stressnormal_->Map()).LID(cnode->dofs()[dof]);
        (*stressnormal_)[locindex[dof]] = -lmn * nn[dof];
      }

      // tangential stress components
      for (int dof = 0; dof < n_dim(); ++dof)
      {
        locindex[dof] = (stresstangential_->Map()).LID(cnode->dofs()[dof]);
        (*stresstangential_)[locindex[dof]] = -lmt1 * nt1[dof] - lmt2 * nt2[dof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::store_dirichlet_status(
    Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmaps)
{
  // loop over all interfaces
  for (unsigned i = 0; i < interfaces().size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // check if this node's dofs are in dbcmap
      for (int k = 0; k < cnode->num_dof(); ++k)
      {
        int currdof = cnode->dofs()[k];
        int lid = (dbcmaps->cond_map())->LID(currdof);

        // store dbc status if found
        if (lid >= 0 && cnode->dbc_dofs()[k] == false) cnode->set_dbc() = true;

        // check compatibility of contact symmetry condition and displacement dirichlet conditions
        if (lid < 0 && cnode->dbc_dofs()[k] == true)
        {
          std::cout << "node " << cnode->id() << " at: " << cnode->x()[0] << " " << cnode->x()[1]
                    << " " << cnode->x()[2] << std::endl;
          std::cout << "dbcdofs: " << cnode->dbc_dofs()[0] << cnode->dbc_dofs()[1]
                    << cnode->dbc_dofs()[2] << std::endl;
          FOUR_C_THROW(
              "Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
        }
      }
    }
  }
  // create old style dirichtoggle vector (supposed to go away)
  pgsdirichtoggle_ = Core::LinAlg::CreateVector(slave_dof_row_map(true), true);
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->cond_map())));
  temp->PutScalar(1.0);
  Core::LinAlg::export_to(*temp, *pgsdirichtoggle_);

  post_store_dirichlet_status(dbcmaps);

  return;
}

/*----------------------------------------------------------------------*
 |  Store D and M last coverged step <-> current step         popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::store_dm(const std::string& state)
{
  // store Dold and Mold matrix in D and M
  if (state == "current")
  {
    dmatrix_ = dold_;
    mmatrix_ = mold_;
  }

  // store D and M matrix in Dold and Mold
  else if (state == "old")
  {
    dold_ = dmatrix_;
    mold_ = mmatrix_;
    if (friction_ && is_dual_quad_slave_trafo()) doldmod_ = dmatrixmod_;
  }

  // unknown conversion
  else
  {
    FOUR_C_THROW("StoreDM: Unknown conversion requested!");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Store nodal quant. to old ones (last conv. time step)  gitterle 02/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::store_to_old(Mortar::StrategyBase::QuantityType type)
{
  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i) interfaces()[i]->store_to_old(type);

  return;
}

/*----------------------------------------------------------------------*
 |  Update and output contact at end of time step             popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::update(Teuchos::RCP<const Epetra_Vector> dis)
{
  // store Lagrange multipliers, D and M
  // (we need this for interpolation of the next generalized mid-point)
  // in the case of self contact, the size of z may have changed
  if (is_self_contact()) zold_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

  zold_->Scale(1.0, *z_);
  store_nodal_quantities(Mortar::StrategyBase::lmold);
  store_dm("old");

  // store contact state to contact nodes (active or inactive)
  store_nodal_quantities(Mortar::StrategyBase::activeold);

  // old displacements in nodes
  // (this is NOT only needed for friction but also for calculating
  // the auxiliary positions in binarytree contact search)
  set_state(Mortar::state_old_displacement, *dis);

  // reset active set status for next time step
  reset_active_set();

  // update flag for global contact status of last time step
  if (gactivenodes_->NumGlobalElements())
  {
    wasincontact_ = true;
    wasincontactlts_ = true;
  }
  else
  {
    wasincontact_ = false;
    wasincontactlts_ = false;
  }

  //----------------------------------------friction: store history values
  // in the case of frictional contact we have to store several
  // information and quantities at the end of a time step (converged
  // state) which is needed in the next time step as history
  // information / quantities.
  if (friction_)
  {
    // store contact state to friction nodes (slip or stick)
    store_nodal_quantities(Mortar::StrategyBase::slipold);

    // store nodal entries of D and M to old ones
    store_to_old(Mortar::StrategyBase::dm);

    // store nodal entries form penalty contact tractions to old ones
    store_to_old(Mortar::StrategyBase::pentrac);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  write restart information for contact                     popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::do_write_restart(
    std::map<std::string, Teuchos::RCP<Epetra_Vector>>& restart_vectors, bool forcedrestart) const
{
  // initalize
  Teuchos::RCP<Epetra_Vector> activetoggle = Teuchos::rcp(new Epetra_Vector(slave_row_nodes()));
  Teuchos::RCP<Epetra_Vector> sliptoggle = Teuchos::null;

  // write toggle
  restart_vectors["activetoggle"] = activetoggle;
  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(slave_row_nodes()));
    restart_vectors["sliptoggle"] = sliptoggle;
  }

  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);
      int dof = (activetoggle->Map()).LID(gid);

      if (forcedrestart)
      {
        // set value active / inactive in toggle vector
        if (cnode->data().active_old()) (*activetoggle)[dof] = 1;
      }
      else
      {
        // set value active / inactive in toggle vector
        if (cnode->active()) (*activetoggle)[dof] = 1;
      }

      // set value slip / stick in the toggle vector
      if (friction_)
      {
        CONTACT::FriNode* frinode = dynamic_cast<CONTACT::FriNode*>(cnode);
        if (forcedrestart)
        {
          if (frinode->fri_data().slip_old()) (*sliptoggle)[dof] = 1;
        }
        else
        {
          if (frinode->fri_data().slip()) (*sliptoggle)[dof] = 1;
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  read restart information for contact                      popp 03/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::do_read_restart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  // check whether this is a restart with contact of a previously
  // non-contact simulation run (if yes, we have to be careful not
  // to try to read certain, in this case non-existing, vectors
  // such as the activetoggle or sliptoggle vectors, but rather
  // initialize the restart active and slip sets as being empty)
  bool restartwithcontact = Core::UTILS::IntegralValue<int>(params(), "RESTART_WITH_CONTACT");

  // set restart displacement state
  set_state(Mortar::state_new_displacement, *dis);
  set_state(Mortar::state_old_displacement, *dis);

  // evaluate interface and restart mortar quantities
  // in the case of SELF CONTACT, also re-setup master/slave maps
  initialize_mortar();
  initialize_and_evaluate_interface(cparams_ptr);
  assemble_mortar();

  //----------------------------------------------------------------------
  // CHECK IF WE NEED TRANSFORMATION MATRICES FOR SLAVE DISPLACEMENT DOFS
  //----------------------------------------------------------------------
  // Concretely, we apply the following transformations:
  // D         ---->   D * T^(-1)
  //----------------------------------------------------------------------
  if (is_dual_quad_slave_trafo())
  {
    // modify dmatrix_
    Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
        Core::LinAlg::MLMultiply(*dmatrix_, false, *invtrafo_, false, false, false, true);
    dmatrix_ = temp;
  }

  // read restart information on active set and slip set (leave sets empty
  // if this is a restart with contact of a non-contact simulation run)
  Teuchos::RCP<Epetra_Vector> activetoggle =
      Teuchos::rcp(new Epetra_Vector(slave_row_nodes(), true));
  if (!restartwithcontact) reader.read_vector(activetoggle, "activetoggle");

  // friction
  Teuchos::RCP<Epetra_Vector> sliptoggle;
  Teuchos::RCP<Epetra_Vector> weightedwear;

  if (friction_)
  {
    sliptoggle = Teuchos::rcp(new Epetra_Vector(slave_row_nodes()));
    if (!restartwithcontact) reader.read_vector(sliptoggle, "sliptoggle");
  }

  // store restart information on active set and slip set
  // into nodes, therefore first loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < (interfaces()[i]->slave_row_nodes())->NumMyElements(); ++j)
    {
      int gid = (interfaces()[i]->slave_row_nodes())->GID(j);
      int dof = (activetoggle->Map()).LID(gid);

      if ((*activetoggle)[dof] == 1)
      {
        Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
        if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
        Node* cnode = dynamic_cast<Node*>(node);

        // set value active / inactive in cnode
        cnode->active() = true;

        if (friction_)
        {
          // set value stick / slip in cnode
          if ((*sliptoggle)[dof] == 1)
            dynamic_cast<CONTACT::FriNode*>(cnode)->fri_data().slip() = true;
        }
      }
    }
  }

  // read restart information on Lagrange multipliers
  z_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  zold_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  if (!restartwithcontact)
    if (!(Global::Problem::instance()->structural_dynamic_params().get<std::string>(
              "INT_STRATEGY") == "Standard" &&
            is_penalty()))
    {
      reader.read_vector(lagrange_multiplier(), "lagrmultold");
      reader.read_vector(lagrange_multiplier_old(), "lagrmultold");
    }

  // Lagrange multiplier increment is always zero (no restart value to be read)
  zincr_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  // store restart information on Lagrange multipliers into nodes
  store_nodal_quantities(Mortar::StrategyBase::lmcurrent);
  store_nodal_quantities(Mortar::StrategyBase::lmold);

  // only for Uzawa augmented strategy
  // TODO: this should be moved to contact_penalty_strategy
  if (stype_ == Inpar::CONTACT::solution_uzawa)
  {
    zuzawa_ = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
    if (!restartwithcontact) reader.read_vector(lagrange_multiplier_uzawa(), "lagrmultold");
    store_nodal_quantities(Mortar::StrategyBase::lmuzawa);
  }

  // store restart Mortar quantities
  store_dm("old");

  if (friction_)
  {
    store_nodal_quantities(Mortar::StrategyBase::activeold);
    store_to_old(Mortar::StrategyBase::dm);
  }

  // (re)setup active global Epetra_Maps
  gactivenodes_ = Teuchos::null;
  gactivedofs_ = Teuchos::null;
  gactiven_ = Teuchos::null;
  gactivet_ = Teuchos::null;
  gslipnodes_ = Teuchos::null;
  gslipdofs_ = Teuchos::null;
  gslipt_ = Teuchos::null;

  // update active sets of all interfaces
  // (these maps are NOT allowed to be overlapping !!!)
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    interfaces()[i]->build_active_set();
    gactivenodes_ = Core::LinAlg::MergeMap(gactivenodes_, interfaces()[i]->active_nodes(), false);
    gactivedofs_ = Core::LinAlg::MergeMap(gactivedofs_, interfaces()[i]->active_dofs(), false);
    gactiven_ = Core::LinAlg::MergeMap(gactiven_, interfaces()[i]->active_n_dofs(), false);
    gactivet_ = Core::LinAlg::MergeMap(gactivet_, interfaces()[i]->active_t_dofs(), false);
    if (friction_)
    {
      gslipnodes_ = Core::LinAlg::MergeMap(gslipnodes_, interfaces()[i]->slip_nodes(), false);
      gslipdofs_ = Core::LinAlg::MergeMap(gslipdofs_, interfaces()[i]->slip_dofs(), false);
      gslipt_ = Core::LinAlg::MergeMap(gslipt_, interfaces()[i]->slip_t_dofs(), false);
    }
  }

  // update flags for global contact status
  if (gactivenodes_->NumGlobalElements())
  {
    isincontact_ = true;
    wasincontact_ = true;
    wasincontactlts_ = true;
  }

  // evaluate relative movement (jump)
  // needed because it is not called in the predictor of the
  // lagrange multiplier strategy
  evaluate_relative_movement();

  // reset unbalance factors for redistribution
  // (during restart the interface has been evaluated once)
  unbalanceEvaluationTime_.resize(0);
  unbalanceNumSlaveElements_.resize(0);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute interface forces (for debugging only)             popp 02/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::interface_forces(bool output)
{
  // check chosen output option
  Inpar::CONTACT::EmOutputType emtype =
      Core::UTILS::IntegralValue<Inpar::CONTACT::EmOutputType>(params(), "EMOUTPUT");

  // get out of here if no output wanted
  if (emtype == Inpar::CONTACT::output_none) return;

  // compute discrete slave and master interface forces
  Teuchos::RCP<Epetra_Vector> fcslavetemp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
  Teuchos::RCP<Epetra_Vector> fcmastertemp =
      Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));

  // for self contact, slave and master sets may have changed,
  // thus we have to export z to new D and M dimensions
  if (is_self_contact())
  {
    Teuchos::RCP<Epetra_Vector> zexp = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
    if (dmatrix_->row_map().NumGlobalElements()) Core::LinAlg::export_to(*z_, *zexp);
    dmatrix_->multiply(true, *zexp, *fcslavetemp);
    mmatrix_->multiply(true, *zexp, *fcmastertemp);
  }
  // if there is no self contact everything is ok
  else
  {
    dmatrix_->multiply(true, *z_, *fcslavetemp);
    mmatrix_->multiply(true, *z_, *fcmastertemp);
  }

  // export the interface forces to full dof layout
  Teuchos::RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Teuchos::RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(*problem_dofs()));
  Core::LinAlg::export_to(*fcslavetemp, *fcslave);
  Core::LinAlg::export_to(*fcmastertemp, *fcmaster);

  // contact forces and moments
  std::vector<double> gfcs(3);
  std::vector<double> ggfcs(3);
  std::vector<double> gfcm(3);
  std::vector<double> ggfcm(3);
  std::vector<double> gmcs(3);
  std::vector<double> ggmcs(3);
  std::vector<double> gmcm(3);
  std::vector<double> ggmcm(3);

  std::vector<double> gmcsnew(3);
  std::vector<double> ggmcsnew(3);
  std::vector<double> gmcmnew(3);
  std::vector<double> ggmcmnew(3);

  // weighted gap vector
  Teuchos::RCP<Epetra_Vector> gapslave = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
  Teuchos::RCP<Epetra_Vector> gapmaster = Teuchos::rcp(new Epetra_Vector(mmatrix_->domain_map()));

  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < n_dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->dofs()[d]);
        if (dofid < 0) FOUR_C_THROW("ContactForces: Did not find slave dof in map");
        nodeforce[d] = (*fcslavetemp)[dofid];
        gfcs[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcs[d] += nodemoment[d];

      // weighted gap
      Core::LinAlg::SerialDenseVector posnode(n_dim());
      std::vector<int> lm(n_dim());
      std::vector<int> lmowner(n_dim());
      for (int d = 0; d < n_dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->dofs()[d];
        lmowner[d] = cnode->owner();
      }
      Core::LinAlg::Assemble(*gapslave, posnode, lm, lmowner);
    }

    // loop over all master nodes on the current interface
    for (int j = 0; j < interfaces()[i]->master_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->master_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      std::vector<double> nodeforce(3);
      std::vector<double> position(3);

      // forces and positions
      for (int d = 0; d < n_dim(); ++d)
      {
        int dofid = (fcmastertemp->Map()).LID(cnode->dofs()[d]);
        if (dofid < 0) FOUR_C_THROW("ContactForces: Did not find master dof in map");
        nodeforce[d] = -(*fcmastertemp)[dofid];
        gfcm[d] += nodeforce[d];
        position[d] = cnode->xspatial()[d];
      }

      // moments
      std::vector<double> nodemoment(3);
      nodemoment[0] = position[1] * nodeforce[2] - position[2] * nodeforce[1];
      nodemoment[1] = position[2] * nodeforce[0] - position[0] * nodeforce[2];
      nodemoment[2] = position[0] * nodeforce[1] - position[1] * nodeforce[0];
      for (int d = 0; d < 3; ++d) gmcm[d] += nodemoment[d];

      // weighted gap
      Core::LinAlg::SerialDenseVector posnode(n_dim());
      std::vector<int> lm(n_dim());
      std::vector<int> lmowner(n_dim());
      for (int d = 0; d < n_dim(); ++d)
      {
        posnode[d] = cnode->xspatial()[d];
        lm[d] = cnode->dofs()[d];
        lmowner[d] = cnode->owner();
      }
      Core::LinAlg::Assemble(*gapmaster, posnode, lm, lmowner);
    }
  }

  // weighted gap
  Teuchos::RCP<Epetra_Vector> gapslavefinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
  Teuchos::RCP<Epetra_Vector> gapmasterfinal = Teuchos::rcp(new Epetra_Vector(mmatrix_->row_map()));
  dmatrix_->multiply(false, *gapslave, *gapslavefinal);
  mmatrix_->multiply(false, *gapmaster, *gapmasterfinal);
  Teuchos::RCP<Epetra_Vector> gapfinal = Teuchos::rcp(new Epetra_Vector(dmatrix_->row_map()));
  gapfinal->Update(1.0, *gapslavefinal, 0.0);
  gapfinal->Update(-1.0, *gapmasterfinal, 1.0);

  // again, for alternative moment: lambda x gap
  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      std::vector<double> lm(3);
      std::vector<double> nodegaps(3);
      std::vector<double> nodegapm(3);

      // LMs and gaps
      for (int d = 0; d < n_dim(); ++d)
      {
        int dofid = (fcslavetemp->Map()).LID(cnode->dofs()[d]);
        if (dofid < 0) FOUR_C_THROW("ContactForces: Did not find slave dof in map");
        nodegaps[d] = (*gapslavefinal)[dofid];
        nodegapm[d] = (*gapmasterfinal)[dofid];
        lm[d] = cnode->mo_data().lm()[d];
      }

      // moments
      std::vector<double> nodemoments(3);
      std::vector<double> nodemomentm(3);
      nodemoments[0] = nodegaps[1] * lm[2] - nodegaps[2] * lm[1];
      nodemoments[1] = nodegaps[2] * lm[0] - nodegaps[0] * lm[2];
      nodemoments[2] = nodegaps[0] * lm[1] - nodegaps[1] * lm[0];
      nodemomentm[0] = nodegapm[1] * lm[2] - nodegapm[2] * lm[1];
      nodemomentm[1] = nodegapm[2] * lm[0] - nodegapm[0] * lm[2];
      nodemomentm[2] = nodegapm[0] * lm[1] - nodegapm[1] * lm[0];
      for (int d = 0; d < 3; ++d)
      {
        gmcsnew[d] += nodemoments[d];
        gmcmnew[d] -= nodemomentm[d];
      }
    }
  }

  // summing up over all processors
  for (int i = 0; i < 3; ++i)
  {
    get_comm().SumAll(&gfcs[i], &ggfcs[i], 1);
    get_comm().SumAll(&gfcm[i], &ggfcm[i], 1);
    get_comm().SumAll(&gmcs[i], &ggmcs[i], 1);
    get_comm().SumAll(&gmcm[i], &ggmcm[i], 1);
    get_comm().SumAll(&gmcsnew[i], &ggmcsnew[i], 1);
    get_comm().SumAll(&gmcmnew[i], &ggmcmnew[i], 1);
  }

  // print interface results to file
  if (emtype == Inpar::CONTACT::output_file || emtype == Inpar::CONTACT::output_both)
  {
    // do this at end of time step only (output==true)!
    // processor 0 does all the work
    if (output && get_comm().MyPID() == 0)
    {
      FILE* MyFile = nullptr;
      std::ostringstream filename;
      const std::string filebase = Global::Problem::instance()->output_control_file()->file_name();
      filename << filebase << ".interface";
      MyFile = fopen(filename.str().c_str(), "at+");

      if (MyFile)
      {
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggfcm[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcs[i]);
        for (int i = 0; i < 3; i++) fprintf(MyFile, "%g\t", ggmcm[i]);
        // for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsfgh[i]);
        // for (int i=0; i<3; i++) fprintf(MyFile, "%g\t", gsmgh[i]);
        fprintf(MyFile, "\n");
        fclose(MyFile);
      }
      else
        FOUR_C_THROW("File for writing meshtying forces could not be opened.");
    }
  }

  // print interface results to screen
  if (emtype == Inpar::CONTACT::output_screen || emtype == Inpar::CONTACT::output_both)
  {
    // do this during Newton steps only (output==false)!
    // processor 0 does all the work
    if (!output && get_comm().MyPID() == 0)
    {
      double snorm = sqrt(ggfcs[0] * ggfcs[0] + ggfcs[1] * ggfcs[1] + ggfcs[2] * ggfcs[2]);
      double mnorm = sqrt(ggfcm[0] * ggfcm[0] + ggfcm[1] * ggfcm[1] + ggfcm[2] * ggfcm[2]);
      printf("Slave Contact Force:   % e  % e  % e \tNorm: % e\n", ggfcs[0], ggfcs[1], ggfcs[2],
          snorm);
      printf("Master Contact Force:  % e  % e  % e \tNorm: % e\n", ggfcm[0], ggfcm[1], ggfcm[2],
          mnorm);
      printf("Slave Contact Moment:  % e  % e  % e\n", ggmcs[0], ggmcs[1], ggmcs[2]);
      // printf("Slave Contact Moment:  % e  % e  % e\n",ggmcsnew[0],ggmcsnew[1],ggmcsnew[2]);
      printf("Master Contact Moment: % e  % e  % e\n", ggmcm[0], ggmcm[1], ggmcm[2]);
      // printf("Master Contact Moment: % e  % e  % e\n",ggmcmnew[0],ggmcmnew[1],ggmcmnew[2]);
      fflush(stdout);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  print interfaces (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::print(std::ostream& os) const
{
  if (get_comm().MyPID() == 0)
  {
    os << "--------------------------------- CONTACT::AbstractStrategy\n"
       << "Contact interfaces: " << (int)interfaces().size() << std::endl
       << "-------------------------------------------------------------\n";
  }
  get_comm().Barrier();
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    std::cout << *(interfaces()[i]);
  }
  get_comm().Barrier();

  return;
}

/*----------------------------------------------------------------------*
 | print active set information                               popp 06/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::print_active_set() const
{
  // output message
  get_comm().Barrier();
  if (get_comm().MyPID() == 0)
  {
    printf("\nActive contact set--------------------------------------------------------------\n");
    fflush(stdout);
  }

  //**********************************************************************
  // detailed active set output
  //**********************************************************************

#ifdef CONTACTASOUTPUT
  // create storage for local and global data
  std::vector<int> lnid, gnid;
  std::vector<double> llmn, glmn;
  std::vector<double> lgap, ggap;

  std::vector<double> Xposl, Xposg;
  std::vector<double> Yposl, Yposg;
  std::vector<double> Zposl, Zposg;

  std::vector<double> xposl, xposg;
  std::vector<double> yposl, yposg;
  std::vector<double> zposl, zposg;

  // introduce integer variable status
  // (0=inactive, 1=active, 2=slip, 3=stick)
  // (this is necessary as all data will be written by proc 0, but
  // the knowledge of the above status ONLY exists on the owner
  // processor of the respective node. Thus this information must
  // also be communicated to proc 0 in addition to the actual data!)
  std::vector<int> lsta, gsta;

  // some more storage for local and global friction data
  std::vector<double> llmt, glmt;
  std::vector<double> ljtx, gjtx;
  std::vector<double> ljte, gjte;
  std::vector<double> lwear, gwear;

  // loop over all interfaces
  for (int i = 0; i < (int)Interfaces().size(); ++i)
  {
    // loop over all slave row nodes on the current interface
    for (int j = 0; j < Interfaces()[i]->SlaveRowNodes()->NumMyElements(); ++j)
    {
      // gid of current node
      int gid = Interfaces()[i]->SlaveRowNodes()->GID(j);
      Core::Nodes::Node* node = Interfaces()[i]->Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

      //--------------------------------------------------------------------
      // FRICTIONLESS CASE
      //--------------------------------------------------------------------
      if (!friction_)
      {
        // cast to Node
        Node* cnode = dynamic_cast<Node*>(node);

        // compute weighted gap
        double wgap = (*wgap_)[wgap_->Map().LID(gid)];

        double Xpos = cnode->X()[0];
        double Ypos = cnode->X()[1];
        double Zpos = cnode->X()[2];

        double xpos = cnode->xspatial()[0];
        double ypos = cnode->xspatial()[1];
        double zpos = cnode->xspatial()[2];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k = 0; k < 3; ++k) nz += cnode->MoData().n()[k] * cnode->MoData().lm()[k];

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        Xposl.push_back(Xpos);
        Yposl.push_back(Ypos);
        Zposl.push_back(Zpos);
        xposl.push_back(xpos);
        yposl.push_back(ypos);
        zposl.push_back(zpos);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active())
          lsta.push_back(1);
        else
          lsta.push_back(0);
      }

      //--------------------------------------------------------------------
      // FRICTIONAL CASE
      //--------------------------------------------------------------------
      else
      {
        // cast to Node and FriNode
        Node* cnode = dynamic_cast<Node*>(node);
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);

        // compute weighted gap
        double wgap = (*wgap_)[wgap_->Map().LID(gid)];

        // compute normal part of Lagrange multiplier
        double nz = 0.0;
        for (int k = 0; k < 3; ++k) nz += frinode->MoData().n()[k] * frinode->MoData().lm()[k];

        // compute tangential parts of Lagrange multiplier and jumps and wear
        double txiz = 0.0;
        double tetaz = 0.0;
        double jumptxi = 0.0;
        double jumpteta = 0.0;
        double wear = 0.0;

        for (int k = 0; k < Dim(); ++k)
        {
          txiz += frinode->data().txi()[k] * frinode->MoData().lm()[k];
          tetaz += frinode->data().teta()[k] * frinode->MoData().lm()[k];
          jumptxi += frinode->data().txi()[k] * frinode->FriData().jump()[k];
          jumpteta += frinode->data().teta()[k] * frinode->FriData().jump()[k];
        }

        // total tangential component
        double tz = sqrt(txiz * txiz + tetaz * tetaz);

        // check for dimensions
        if (Dim() == 2 && abs(jumpteta) > 0.0001)
          FOUR_C_THROW("Error: Jumpteta should be zero for 2D");

        // store node id
        lnid.push_back(gid);

        // store relevant data
        llmn.push_back(nz);
        lgap.push_back(wgap);
        llmt.push_back(tz);
        ljtx.push_back(jumptxi);
        ljte.push_back(jumpteta);
        lwear.push_back(wear);

        // store status (0=inactive, 1=active, 2=slip, 3=stick)
        if (cnode->Active())
        {
          if (frinode->FriData().Slip())
            lsta.push_back(2);
          else
            lsta.push_back(3);
        }
        else
        {
          lsta.push_back(0);
        }
      }
    }
  }

  // we want to gather data from on all procs
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;

  // communicate all data to proc 0
  Core::LinAlg::Gather<int>(lnid, gnid, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(llmn, glmn, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(lgap, ggap, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<int>(lsta, gsta, (int)allproc.size(), allproc.data(), Comm());

  Core::LinAlg::Gather<double>(Xposl, Xposg, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(Yposl, Yposg, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(Zposl, Zposg, (int)allproc.size(), allproc.data(), Comm());

  Core::LinAlg::Gather<double>(xposl, xposg, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(yposl, yposg, (int)allproc.size(), allproc.data(), Comm());
  Core::LinAlg::Gather<double>(zposl, zposg, (int)allproc.size(), allproc.data(), Comm());

  // communicate some more data to proc 0 for friction
  if (friction_)
  {
    Core::LinAlg::Gather<double>(llmt, glmt, (int)allproc.size(), allproc.data(), Comm());
    Core::LinAlg::Gather<double>(ljtx, gjtx, (int)allproc.size(), allproc.data(), Comm());
    Core::LinAlg::Gather<double>(ljte, gjte, (int)allproc.size(), allproc.data(), Comm());
    Core::LinAlg::Gather<double>(lwear, gwear, (int)allproc.size(), allproc.data(), Comm());
  }

  // output is solely done by proc 0
  if (Comm().MyPID() == 0)
  {
    //--------------------------------------------------------------------
    // FRICTIONLESS CASE
    //--------------------------------------------------------------------
    if (!friction_)
    {
      // loop over all nodes
      for (int k = 0; k < (int)gnid.size(); ++k)
      {
        // print nodes of inactive set *************************************
        if (gsta[k] == 0)
        {
          printf("INACTIVE: %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",
              gnid[k], ggap[k], glmn[k], Xposg[k], Yposg[k], Zposg[k]);
          fflush(stdout);
        }

        // print nodes of active set ***************************************
        else if (gsta[k] == 1)
        {
          printf("ACTIVE:   %d \t wgap: % e \t lm: % e \t Xref: % e \t Yref: % e \t Zref: % e \n",
              gnid[k], ggap[k], glmn[k], Xposg[k], Yposg[k], Zposg[k]);
          fflush(stdout);
        }

        // invalid status **************************************************
        else
          FOUR_C_THROW("Invalid node status %i for frictionless case", gsta[k]);
      }
    }

    //--------------------------------------------------------------------
    // FRICTIONAL CASE
    //--------------------------------------------------------------------
    else
    {
#ifdef CONTACTEXPORT
      // export variables
      double sum_jumpx = 0.0;
      double sum_jumpe = 0.0;
      double sum_jumpall = 0.0;
      double iter_slip = 0.0;
#endif

      // loop over all nodes
      for (int k = 0; k < (int)gnid.size(); ++k)
      {
        // print nodes of slip set **************************************
        if (gsta[k] == 2)
        {
          printf("SLIP:  %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",
              gnid[k], glmn[k], glmt[k], gjtx[k], gjte[k], gwear[k]);
          fflush(stdout);
#ifdef CONTACTEXPORT
          // preparation for output
          sum_jumpx += gjtx[k];
          sum_jumpe += gjte[k];
          sum_jumpall += sqrt(gjtx[k] * gjtx[k] + gjte[k] * gjte[k]);
          iter_slip = iter_slip + 1.0;
#endif
        }

        // print nodes of stick set *************************************
        else if (gsta[k] == 3)
        {
          printf("STICK: %d \t lm_n: % e \t lm_t: % e \t jump1: % e \t jump2: % e \t wear: % e \n",
              gnid[k], glmn[k], glmt[k], gjtx[k], gjte[k], gwear[k]);
          fflush(stdout);
        }

        // print nodes of inactive set *************************************
        else if (gsta[k] == 0)
        {
          // do nothing
        }

        // invalid status **************************************************
        else
          FOUR_C_THROW("Invalid node status %i for frictional case", gsta[k]);
      }

#ifdef CONTACTEXPORT
      // export averaged slip increments to xxx.jump
      double sum_jumpx_final = 0.0;
      double sum_jumpe_final = 0.0;
      double sum_jumpall_final = 0.0;

      if (iter_slip > 0.0)
      {
        sum_jumpx_final = sum_jumpx / iter_slip;
        sum_jumpe_final = sum_jumpe / iter_slip;
        sum_jumpall_final = sum_jumpall / iter_slip;
      }

      FILE* MyFile = nullptr;
      std::ostringstream filename;
      const std::string filebase =
          Global::Problem::instance()->OutputControlFile()->file_name_only_prefix();
      filename << filebase << ".jump";
      MyFile = fopen(filename.str().c_str(), "at+");
      if (MyFile)
      {
        // fprintf(MyFile,valuename.c_str());
        fprintf(MyFile, "%g\t", sum_jumpx_final);
        fprintf(MyFile, "%g\t", sum_jumpe_final);
        fprintf(MyFile, "%g\n", sum_jumpall_final);
        fclose(MyFile);
      }
      else
        FOUR_C_THROW("File for Output could not be opened.");
#endif
    }
  }

#else
  //**********************************************************************
  // reduced active set output
  //**********************************************************************

  // counters
  int activenodes = 0;
  int gactivenodes = 0;
  int inactivenodes = 0;
  int ginactivenodes = 0;
  int slipnodes = 0;
  int gslipnodes = 0;

  // counters for non-smooth contact
  int surfacenodes = 0;
  int gsurfacenodes = 0;
  int edgenodes = 0;
  int gedgenodes = 0;
  int cornernodes = 0;
  int gcornernodes = 0;

  // nonsmooth contact active?
  bool nonsmooth = Core::UTILS::IntegralValue<int>(params(), "NONSMOOTH_GEOMETRIES");

  // loop over all interfaces
  for (int i = 0; i < (int)interfaces().size(); ++i)
  {
    // loop over all slave nodes on the current interface
    for (int j = 0; j < interfaces()[i]->slave_row_nodes()->NumMyElements(); ++j)
    {
      int gid = interfaces()[i]->slave_row_nodes()->GID(j);
      Core::Nodes::Node* node = interfaces()[i]->discret().g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);

      // increase active counters
      Node* cnode = dynamic_cast<Node*>(node);

      if (cnode->active())
        activenodes += 1;
      else
        inactivenodes += 1;

      // increase friction counters
      if (friction_)
      {
        FriNode* frinode = dynamic_cast<FriNode*>(cnode);
        if (cnode->active() && frinode->fri_data().slip()) slipnodes += 1;
      }

      // get nonsmooth contact states
      if (nonsmooth)
      {
        if (cnode->active() && cnode->is_on_edge() && !cnode->is_on_corner()) edgenodes += 1;
        if (cnode->active() && cnode->is_on_corner()) cornernodes += 1;
        if (cnode->active() && !cnode->is_on_corner_edge()) surfacenodes += 1;
      }
    }
  }

  // sum among all processors
  get_comm().SumAll(&activenodes, &gactivenodes, 1);
  get_comm().SumAll(&inactivenodes, &ginactivenodes, 1);
  get_comm().SumAll(&slipnodes, &gslipnodes, 1);
  get_comm().SumAll(&edgenodes, &gedgenodes, 1);
  get_comm().SumAll(&cornernodes, &gcornernodes, 1);
  get_comm().SumAll(&surfacenodes, &gsurfacenodes, 1);

  // print active set information
  if (get_comm().MyPID() == 0)
  {
    if (nonsmooth)
    {
      std::cout << "Total ACTIVE SURFACE nodes:\t" << gsurfacenodes << std::endl;
      std::cout << "Total    ACTIVE EDGE nodes:\t" << gedgenodes << std::endl;
      std::cout << "Total  ACTIVE CORNER nodes:\t" << gcornernodes << std::endl;
      std::cout << "Total       INACTIVE nodes:\t" << ginactivenodes << std::endl;
    }
    else if (friction_)
    {
      std::cout << "Total     SLIP nodes:\t" << gslipnodes << std::endl;
      std::cout << "Total    STICK nodes:\t" << gactivenodes - gslipnodes << std::endl;
      std::cout << "Total INACTIVE nodes:\t" << ginactivenodes << std::endl;
    }
    else
    {
      std::cout << "Total   ACTIVE nodes:\t" << gactivenodes << std::endl;
      std::cout << "Total INACTIVE nodes:\t" << ginactivenodes << std::endl;
    }
  }
#endif  // #ifdef CONTACTASOUTPUT
  // output line
  get_comm().Barrier();
  if (get_comm().MyPID() == 0)
  {
    printf("--------------------------------------------------------------------------------\n\n");
    fflush(stdout);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Visualization of contact segments with gmsh                popp 08/08|
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::visualize_gmsh(const int step, const int iter)
{
  // visualization with gmsh
  for (int i = 0; i < (int)interfaces().size(); ++i) interfaces()[i]->visualize_gmsh(step, iter);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::collect_maps_for_preconditioner(
    Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap,
    Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap) const
{
  InnerDofMap = gndofrowmap_;   // global internal dof row map
  ActiveDofMap = gactivedofs_;  // global active slave dof row map

  // check if parallel redistribution is used
  // if parallel redistribution is activated, then use (original) maps before redistribution
  // otherwise we use just the standard master/slave maps
  if (pgsdofrowmap_ != Teuchos::null)
    SlaveDofMap = pgsdofrowmap_;
  else
    SlaveDofMap = gsdofrowmap_;
  if (pgmdofrowmap_ != Teuchos::null)
    MasterDofMap = pgmdofrowmap_;
  else
    MasterDofMap = gmdofrowmap_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  set_state(Mortar::state_new_displacement, dispnp);
  reset_lagrange_multipliers(cparams, xnew);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate(CONTACT::ParamsInterface& cparams,
    const std::vector<Teuchos::RCP<const Epetra_Vector>>* eval_vec,
    const std::vector<Teuchos::RCP<Epetra_Vector>>* eval_vec_mutable)
{
  pre_evaluate(cparams);

  const enum Mortar::ActionType& act = cparams.get_action_type();
  switch (act)
  {
    // -------------------------------------------------------------------
    // evaluate only the contact forces / contact right hand side
    // -------------------------------------------------------------------
    case Mortar::eval_force:
    {
      evaluate_force(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // evaluate the contact matrix blocks and the rhs contributions
    // -------------------------------------------------------------------
    case Mortar::eval_force_stiff:
    {
      evaluate_force_stiff(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // run before an evaluate call in the Solid::ModelEvaluator class
    // -------------------------------------------------------------------
    case Mortar::eval_run_pre_evaluate:
    {
      run_pre_evaluate(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // run after an evaluate call in the Solid::ModelEvaluator class
    // -------------------------------------------------------------------
    case Mortar::eval_run_post_evaluate:
    {
      run_post_evaluate(cparams);
      break;
    }
    // -------------------------------------------------------------------
    // reset internal stored solution quantities (e.g.
    // displacement state, Lagrange multi.)
    // -------------------------------------------------------------------
    case Mortar::eval_reset:
    {
      if (not eval_vec) FOUR_C_THROW("Missing evaluation vectors!");
      if (eval_vec->size() != 2)
        FOUR_C_THROW(
            "The \"Mortar::eval_reset\" action expects \n"
            "exactly 2 evaluation vector pointers! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());
      const Epetra_Vector& dispnp = *((*eval_vec)[0]);
      const Epetra_Vector& xnew = *((*eval_vec)[1]);
      reset(cparams, dispnp, xnew);

      break;
    }
    // -------------------------------------------------------------------
    // recover internal stored solution quantities (e.g. Lagrange multi.)
    // -------------------------------------------------------------------
    case Mortar::eval_run_post_compute_x:
    {
      if (not eval_vec) FOUR_C_THROW("Missing evaluation vectors!");

      if (eval_vec->size() != 3)
        FOUR_C_THROW(
            "The \"Mortar::eval_recover\" action expects \n"
            "exactly 3 evaluation vector pointers! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      const Teuchos::RCP<const Epetra_Vector>& xold_ptr = (*eval_vec)[0];
      if (xold_ptr.is_null() or !xold_ptr.is_valid_ptr()) FOUR_C_THROW("xold_ptr is undefined!");

      const Teuchos::RCP<const Epetra_Vector>& dir_ptr = (*eval_vec)[1];
      if (dir_ptr.is_null() or !dir_ptr.is_valid_ptr()) FOUR_C_THROW("dir_ptr is undefined!");

      const Teuchos::RCP<const Epetra_Vector>& xnew_ptr = (*eval_vec)[2];
      if (xnew_ptr.is_null() or !xnew_ptr.is_valid_ptr()) FOUR_C_THROW("xnew_ptr is undefined!");

      run_post_compute_x(cparams, *xold_ptr, *dir_ptr, *xnew_ptr);

      break;
    }
    case Mortar::eval_run_pre_compute_x:
    {
      if (not eval_vec) FOUR_C_THROW("Missing constant evaluation vectors!");
      if (not eval_vec_mutable) FOUR_C_THROW("Missing mutable evaluation vectors!");

      if (eval_vec->size() != 1)
        FOUR_C_THROW(
            "The \"Mortar::eval_augment_direction\" action expects \n"
            "exactly 1 constant evaluation vector pointer! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      if (eval_vec_mutable->size() != 1)
        FOUR_C_THROW(
            "The \"Mortar::eval_augment_direction\" action expects \n"
            "exactly 1 mutable evaluation vector pointer! But you \n"
            "passed %i vector pointers!",
            eval_vec->size());

      const Teuchos::RCP<const Epetra_Vector>& xold_ptr = eval_vec->front();
      if (xold_ptr.is_null()) FOUR_C_THROW("Missing xold vector!");

      const Teuchos::RCP<Epetra_Vector>& dir_mutable_ptr = eval_vec_mutable->front();
      if (dir_mutable_ptr.is_null()) FOUR_C_THROW("Missing dir_mutable vector!");

      run_pre_compute_x(cparams, *xold_ptr, *dir_mutable_ptr);

      break;
    }
    case Mortar::eval_run_post_iterate:
    {
      run_post_iterate(cparams);

      break;
    }
    case Mortar::eval_run_post_apply_jacobian_inverse:
    {
      const Epetra_Vector* rhs = cparams.get<const Epetra_Vector>(0);
      Epetra_Vector* result = cparams.get<Epetra_Vector>(1);
      const Epetra_Vector* xold = cparams.get<const Epetra_Vector>(2);
      const NOX::Nln::Group* grp = cparams.get<const NOX::Nln::Group>(3);

      run_post_apply_jacobian_inverse(cparams, *rhs, *result, *xold, *grp);

      break;
    }
    case Mortar::eval_correct_parameters:
    {
      const NOX::Nln::CorrectionType* type = cparams.get<const NOX::Nln::CorrectionType>(0);

      correct_parameters(cparams, *type);

      break;
    }
    case Mortar::eval_wgap_gradient_error:
    {
      evaluate_weighted_gap_gradient_error(cparams);

      break;
    }
    case Mortar::eval_static_constraint_rhs:
    {
      evaluate_static_constraint_rhs(cparams);

      break;
    }
    case Mortar::remove_condensed_contributions_from_str_rhs:
    {
      if (not eval_vec_mutable) FOUR_C_THROW("The mutable evaluation vector is expected!");
      if (eval_vec_mutable->size() < 1)
        FOUR_C_THROW(
            "The eval_vec_mutable is supposed to have at least a length"
            " of 1!");

      Epetra_Vector& str_rhs = *eval_vec_mutable->front();
      remove_condensed_contributions_from_rhs(str_rhs);

      break;
    }
    case Mortar::eval_run_pre_solve:
    {
      if (not eval_vec) FOUR_C_THROW("The read-only evaluation vector is expected!");
      if (eval_vec->size() < 1)
        FOUR_C_THROW(
            "The eval_vec is supposed to have at least a length"
            " of 1!");

      const Teuchos::RCP<const Epetra_Vector>& curr_disp = eval_vec->front();
      run_pre_solve(curr_disp, cparams);

      break;
    }
    // -------------------------------------------------------------------
    // no suitable action could be found
    // -------------------------------------------------------------------
    default:
    {
      FOUR_C_THROW("Unsupported action type: %i | %s", act, ActionType2String(act).c_str());
      break;
    }
  }

  post_evaluate(cparams);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_force(CONTACT::ParamsInterface& cparams)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_force_stiff(CONTACT::ParamsInterface& cparams)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_static_constraint_rhs(CONTACT::ParamsInterface& cparams)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::remove_condensed_contributions_from_rhs(
    Epetra_Vector& str_rhs) const
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_pre_evaluate(CONTACT::ParamsInterface& cparams)
{
  /* Not yet implemented! See the CONTACT::AUG framework for an
   * example. */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_post_evaluate(CONTACT::ParamsInterface& cparams)
{
  /* Not yet implemented! See the CONTACT::Aug::ComboStrategy for an
   * example. */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_post_compute_x(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_pre_compute_x(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold, Epetra_Vector& dir_mutable)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_post_iterate(const CONTACT::ParamsInterface& cparams)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_pre_solve(
    const Teuchos::RCP<const Epetra_Vector>& curr_disp, const CONTACT::ParamsInterface& cparams)
{
  // do nothing
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::run_post_apply_jacobian_inverse(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& rhs, Epetra_Vector& result,
    const Epetra_Vector& xold, const NOX::Nln::Group& grp)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::evaluate_weighted_gap_gradient_error(
    CONTACT::ParamsInterface& cparams)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::reset_lagrange_multipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  FOUR_C_THROW(
      "Not yet implemented! See the CONTACT::Aug::Strategy for an "
      "example.");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::correct_parameters(
    CONTACT::ParamsInterface& cparams, const NOX::Nln::CorrectionType type)
{
  /* do nothing */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::is_saddle_point_system() const
{
  if ((stype_ == Inpar::CONTACT::solution_lagmult) and
      system_type() == Inpar::CONTACT::system_saddlepoint)
  {
    if (is_in_contact() or was_in_contact() or was_in_contact_last_time_step()) return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::is_condensed_system() const
{
  if (stype_ == Inpar::CONTACT::solution_lagmult and
      system_type() != Inpar::CONTACT::system_saddlepoint)
  {
    if (is_in_contact() or was_in_contact() or was_in_contact_last_time_step()) return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::fill_maps_for_preconditioner(
    std::vector<Teuchos::RCP<Epetra_Map>>& maps) const
{
  /* FixMe This function replaces the deprecated collect_maps_for_preconditioner(),
   * the old version can be deleted, as soon as the contact uses the new
   * structure framework. */

  if (maps.size() != 4) FOUR_C_THROW("The vector size has to be 4!");
  /* check if parallel redistribution is used
   * if parallel redistribution is activated, then use (original) maps
   * before redistribution otherwise we use just the standard master/slave
   * maps */

  // (0) masterDofMap
  if (pgmdofrowmap_ != Teuchos::null)
    maps[0] = pgmdofrowmap_;
  else
    maps[0] = gmdofrowmap_;

  // (1) slaveDofMap
  if (pgsdofrowmap_ != Teuchos::null)
    maps[1] = pgsdofrowmap_;
  else
    maps[1] = gsdofrowmap_;

  // (2) innerDofMap
  maps[2] = gndofrowmap_;

  // (3) activeDofMap
  maps[3] = gactivedofs_;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  FOUR_C_THROW("Not implemented!");
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AbstractStrategy::lagrange_multiplier_np(
    const bool& redist) const
{
  if ((redist) or not parallel_redistribution_status()) return z_;

  Teuchos::RCP<Epetra_Vector> z_unredist =
      Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(false)));
  Core::LinAlg::export_to(*z_, *z_unredist);
  return z_unredist;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::AbstractStrategy::lagrange_multiplier_n(
    const bool& redist) const
{
  if ((redist) or not parallel_redistribution_status()) return zold_;

  Teuchos::RCP<Epetra_Vector> zold_unredist =
      Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(false)));
  Core::LinAlg::export_to(*zold_, *zold_unredist);
  return zold_unredist;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AbstractStrategy::get_potential_value(
    const enum NOX::Nln::MeritFunction::MeritFctName mrt_type) const
{
  FOUR_C_THROW("The currently active strategy \"%s\" does not support this method!",
      Inpar::CONTACT::SolvingStrategy2String(type()).c_str());
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::AbstractStrategy::get_linearized_potential_value_terms(const Epetra_Vector& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName mrt_type,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  FOUR_C_THROW("The currently active strategy \"%s\" does not support this method!",
      Inpar::CONTACT::SolvingStrategy2String(type()).c_str());
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::AbstractStrategy::postprocess_quantities_per_interface(
    Teuchos::RCP<Teuchos::ParameterList> outputParams)
{
  using Teuchos::RCP;

  // Evaluate slave and master forces
  {
    RCP<Epetra_Vector> fcslave = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true), true));
    RCP<Epetra_Vector> fcmaster = Teuchos::rcp(new Epetra_Vector(master_dof_row_map(true), true));

    // Mortar matrices might not be initialized, e.g. in the initial state. If so, keep zero vector.
    if (!d_matrix().is_null()) d_matrix()->multiply(true, *zold_, *fcslave);
    if (!m_matrix().is_null()) m_matrix()->multiply(true, *zold_, *fcmaster);

    // Append data to parameter list
    outputParams->set<RCP<const Epetra_Vector>>("interface traction", zold_);
    outputParams->set<RCP<const Epetra_Vector>>("slave forces", fcslave);
    outputParams->set<RCP<const Epetra_Vector>>("master forces", fcmaster);
  }

  // Postprocess contact stresses
  {
    compute_contact_stresses();

    // Append data to parameter list
    outputParams->set<RCP<const Epetra_Vector>>("norcontactstress", stressnormal_);
    outputParams->set<RCP<const Epetra_Vector>>("tancontactstress", stresstangential_);
  }

  for (auto& interface : interfaces()) interface->postprocess_quantities(*outputParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::AbstractStrategy::is_first_time_step() const
{
  bool first_time_step = false;
  if (unbalanceEvaluationTime_.size() == 0 && unbalanceNumSlaveElements_.size() == 0)
    first_time_step = true;

  return first_time_step;
}

FOUR_C_NAMESPACE_CLOSE
