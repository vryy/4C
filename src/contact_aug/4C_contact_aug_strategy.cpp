/*----------------------------------------------------------------------*/
/*! \file
\brief Augmented contact solving strategy with standard Lagrangian
       multipliers.

\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_contact_aug_strategy.hpp"

#include "4C_contact_aug_active_set.hpp"
#include "4C_contact_aug_interface.hpp"
#include "4C_contact_aug_parallel_distribution_controller.hpp"
#include "4C_contact_aug_potential.hpp"
#include "4C_contact_aug_steepest_ascent_strategy.hpp"
#include "4C_contact_aug_timemonitor.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_lagrange_strategy.hpp"
#include "4C_contact_node.hpp"
#include "4C_contact_paramsinterface.hpp"
#include "4C_contact_utils.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_mortar_matrix_transform.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_utils_epetra_exceptions.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::DataContainer::DataContainer()
    : wasincontactlastiter_(false),
      isactivesetconverged_(false),
      printlinearconservation_(false),
      printangularconservation_(false),
      is_semi_smooth_newton_(false),
      matrix_maps_valid_(false),
      vector_maps_valid_(false),
      cn_(-1.0),
      eval_state_(Mortar::eval_none),
      ghosting_strategy_(Inpar::Mortar::ExtendGhosting::redundant_master),
      var_type_(Inpar::CONTACT::var_unknown),
      fd_check_type_(Inpar::CONTACT::FdCheck::off),
      potentialPtr_(Teuchos::null),
      mat_row_col_transformer_(Teuchos::null),
      BMatrixPtr_(Teuchos::null),
      dGLmLinMatrixPtr_(Teuchos::null),
      dGGLinMatrixPtr_(Teuchos::null),
      dLmNWGapLinMatrixPtr_(Teuchos::null),
      dLmTLmTMatrixPtr_(Teuchos::null),
      dLmTLmTLinMatrixPtr_(Teuchos::null),
      inactiveLinMatrixPtr_(Teuchos::null),
      inactiveDiagMatrixPtr_(Teuchos::null),
      aPtr_(Teuchos::null),
      kappaPtr_(Teuchos::null),
      lmNPtr_(Teuchos::null),
      aWGapPtr_(Teuchos::null),
      wGapAllPtr_(Teuchos::null),
      dLmTLmTRhsPtr_(Teuchos::null),
      slForceLmPtr_(Teuchos::null),
      slForceGPtr_(Teuchos::null),
      maForceLmPtr_(Teuchos::null),
      maForceGPtr_(Teuchos::null),
      cnPtr_(Teuchos::null),
      gsndofrowmapPtr_(Teuchos::null),
      gstdofrowmapPtr_(Teuchos::null),
      gOldActiveSlaveNodesPtr_(Teuchos::null)
{
  // empty constructor

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::DataContainer::init_matrix_row_col_transformer()
{
  Mortar::MatrixRowColTransformer::plain_block_map_pairs redistributed_row(4);
  redistributed_row[CONTACT::MatBlockType::displ_displ] = &global_slave_master_dof_row_map_ptr();
  redistributed_row[CONTACT::MatBlockType::displ_lm] = &global_disp_dof_row_map_ptr();
  redistributed_row[CONTACT::MatBlockType::lm_displ] = &global_lm_dof_row_map_ptr();
  redistributed_row[CONTACT::MatBlockType::lm_lm] = &global_lm_dof_row_map_ptr();

  Mortar::MatrixRowColTransformer::plain_block_map_pairs redistributed_col(4);
  redistributed_col[CONTACT::MatBlockType::displ_displ] = &global_slave_master_dof_row_map_ptr();
  redistributed_col[CONTACT::MatBlockType::displ_lm] = &global_lm_dof_row_map_ptr();
  redistributed_col[CONTACT::MatBlockType::lm_displ] = &global_disp_dof_row_map_ptr();
  redistributed_col[CONTACT::MatBlockType::lm_lm] = &global_lm_dof_row_map_ptr();

  Mortar::MatrixRowColTransformer::plain_block_map_pairs unredistributed_row(4);
  unredistributed_row[CONTACT::MatBlockType::displ_displ] = &PGSlMaDofRowMapPtr();
  unredistributed_row[CONTACT::MatBlockType::displ_lm] = &ProbDofsPtr();
  unredistributed_row[CONTACT::MatBlockType::lm_displ] = &PGLmDofRowMapPtr();
  unredistributed_row[CONTACT::MatBlockType::lm_lm] = &PGLmDofRowMapPtr();


  Mortar::MatrixRowColTransformer::plain_block_map_pairs unredistributed_col(4);
  unredistributed_col[CONTACT::MatBlockType::displ_displ] = &PGSlMaDofRowMapPtr();
  unredistributed_col[CONTACT::MatBlockType::displ_lm] = &PGLmDofRowMapPtr();
  unredistributed_col[CONTACT::MatBlockType::lm_displ] = &ProbDofsPtr();
  unredistributed_col[CONTACT::MatBlockType::lm_lm] = &PGLmDofRowMapPtr();

  mat_row_col_transformer_->init(
      redistributed_row, redistributed_col, unredistributed_row, unredistributed_col);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::DataContainer::init_sub_data_container(
    const Inpar::CONTACT::SolvingStrategy strat_type)
{
  switch (strat_type)
  {
    case Inpar::CONTACT::solution_steepest_ascent:
    case Inpar::CONTACT::solution_steepest_ascent_sp:
      // do not initialize it twice (combo strategy)
      if (sa_data_ptr_.is_null())
        sa_data_ptr_ = Teuchos::RCP<CONTACT::Aug::SteepestAscent::DataContainer>(
            new CONTACT::Aug::SteepestAscent::DataContainer());
      break;
    default:
      FOUR_C_THROW(
          "There is no known sub-data container for the given "
          "strategy type! ( strat_type = %s )",
          Inpar::CONTACT::SolvingStrategy2String(strat_type).c_str());
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::Strategy::Strategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data_ptr,
    const Epetra_Map* dof_row_map, const Epetra_Map* NodeRowMap,
    const Teuchos::ParameterList& params, const plain_interface_set& interfaces, int dim,
    const Teuchos::RCP<const Epetra_Comm>& comm, int maxdof)
    : CONTACT::AbstractStrategy(data_ptr, dof_row_map, NodeRowMap, params, dim, comm, 0.0, maxdof),
      aug_data_ptr_(Teuchos::rcp_dynamic_cast<CONTACT::Aug::DataContainer>(data_ptr, true)),
      aug_data_(*aug_data_ptr_)
{
  // store values of the parameter list
  const Teuchos::ParameterList& p_aug = params.sublist("AUGMENTED");

  data().print_linear_mom_conservation() =
      Core::UTILS::IntegralValue<bool>(p_aug, "PRINT_LINEAR_CONSERVATION");

  data().print_angular_mom_conservation() =
      Core::UTILS::IntegralValue<bool>(p_aug, "PRINT_ANGULAR_CONSERVATION");

  data().add_inactiv_force_contributions() =
      Core::UTILS::IntegralValue<bool>(p_aug, "ADD_INACTIVE_FORCE_CONTRIBUTIONS");

  data().set_is_semi_smooth_newton(Core::UTILS::IntegralValue<bool>(params, "SEMI_SMOOTH_NEWTON"));

  data().SetConstantCn(Params().get<double>("SEMI_SMOOTH_CN"));

  data().SetGhostingStrategy(Teuchos::getIntegralValue<Inpar::Mortar::ExtendGhosting>(
      Params().sublist("PARALLEL REDISTRIBUTION"), "GHOSTING_STRATEGY"));

  data().set_variational_approach_type(
      Core::UTILS::IntegralValue<Inpar::CONTACT::VariationalApproach>(
          p_aug, "VARIATIONAL_APPROACH"));

  if (data().add_inactiv_force_contributions() and
      data().variational_approach_type() != Inpar::CONTACT::var_complete)
    FOUR_C_THROW(
        "The \"ADD_INACTIVE_FORCE_CONTRIBUTIONS\" option is only supported "
        "by the complete variational approach.");

  data().SetPotential(Teuchos::rcp(new Potential(*this, *aug_data_ptr_)));

  data().set_matrix_row_col_transformer(Teuchos::rcp(new Mortar::MatrixRowColTransformer(4)));

  const int par_redist_interval = p_aug.get<int>("PARALLEL_REDIST_INTERVAL");
  data().SetPDController(
      Teuchos::rcp(new ParallelDistributionController(*this, data(), par_redist_interval)));

  interface_.reserve(interfaces.size());
  for (plain_interface_set::const_iterator cit = interfaces.begin(); cit != interfaces.end(); ++cit)
  {
    const Teuchos::RCP<CONTACT::Interface>& interface = *cit;
    // cast to augmented interfaces just as sanity check
    interface_.push_back(Teuchos::rcp_dynamic_cast<CONTACT::Aug::Interface>(interface, true));
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::Aug::Strategy::is_saddle_point_system() const
{
  return (is_in_contact() or was_in_contact() or was_in_contact_last_time_step());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::Aug::plain_interface_set& CONTACT::Aug::Strategy::interfaces() { return interface_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const CONTACT::Aug::plain_interface_set& CONTACT::Aug::Strategy::interfaces() const
{
  return interface_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::post_setup(bool redistributed, bool init)
{
  if (init or redistributed)
  {
    // reassemble the global slave normal/tangential dof row maps
    assemble_global_sl_nt_dof_row_maps();

    assemble_global_ele_maps();
  }

  if (init)
  {
    // initialize cn
    initialize_cn(data().ConstantCn());

    if (parallel_redistribution_status()) data().init_matrix_row_col_transformer();
  }

  // just used for the redistributed case
  if (redistributed)
  {
    // redistribute the cn-vector
    redistribute_cn();

    // redistribute the global augmented old active slave nodes map
    if ((not data().g_old_active_slave_nodes_ptr().is_null()) and
        (data().g_old_active_slave_nodes().NumGlobalElements() > 0))
    {
      data().g_old_active_slave_nodes_ptr() = Core::Rebalance::RebalanceInAccordanceWithReference(
          slave_row_nodes(), data().g_old_active_slave_nodes());
    }
  }

  // in both cases the maps change and we have to re-build all matrices
  data().SetMatrixMapsValid(false);
  data().SetVectorMapsValid(false);

  // setup the potential class with the current maps
  data().Potential().setup();

  // setup the row column transformer object
  if (parallel_redistribution_status()) data().matrix_row_col_transformer().setup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::Update(Teuchos::RCP<const Epetra_Vector> dis)
{
  CONTACT::AbstractStrategy::Update(dis);
  initialize_cn(data().ConstantCn());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::assemble_global_sl_nt_dof_row_maps()
{
  data().g_sl_normal_dof_row_map_ptr() = Teuchos::null;
  data().g_sl_tangential_dof_row_map_ptr() = Teuchos::null;

  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    data().g_sl_normal_dof_row_map_ptr() =
        Core::LinAlg::MergeMap(data().g_sl_normal_dof_row_map_ptr(), interface.SlaveRowNDofs());
    data().g_sl_tangential_dof_row_map_ptr() =
        Core::LinAlg::MergeMap(data().g_sl_tangential_dof_row_map_ptr(), interface.SlaveRowTDofs());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::assemble_global_ele_maps()
{
  data().GSeleRowMapPtr() = Teuchos::null;
  data().GSeleColMapPtr() = Teuchos::null;

  data().GMeleRowMapPtr() = Teuchos::null;
  data().GMeleColMapPtr() = Teuchos::null;

  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    data().GSeleRowMapPtr() =
        Core::LinAlg::MergeMap(data().GSeleRowMapPtr(), (*cit)->SlaveRowElements());
    data().GSeleColMapPtr() =
        Core::LinAlg::MergeMap(data().GSeleColMapPtr(), (*cit)->SlaveColElements());

    data().GMeleRowMapPtr() =
        Core::LinAlg::MergeMap(data().GMeleRowMapPtr(), (*cit)->MasterRowElements());
    data().GMeleColMapPtr() =
        Core::LinAlg::MergeMap(data().GMeleColMapPtr(), (*cit)->MasterColElements());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::initialize_cn(const double cn_init)
{
  if (cn_init <= 0.0) FOUR_C_THROW("The initial CN value must be greater than zero!");

  if (data().CnPtr().is_null() or data().Cn().GlobalLength() == 0)
    data().CnPtr() = Core::LinAlg::CreateVector(slave_row_nodes(), true);

  // set all nodal cn-values to the input value
  data().Cn().PutScalar(cn_init);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::redistribute_cn()
{
  // redistribute the cn-vector
  Teuchos::RCP<Epetra_Vector> newcn = Teuchos::rcp(new Epetra_Vector(slave_row_nodes()));
  Core::LinAlg::Export(data().Cn(), *newcn);
  data().CnPtr() = newcn;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::DoReadRestart(Core::IO::DiscretizationReader& reader,
    Teuchos::RCP<const Epetra_Vector> dis, Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  CONTACT::AbstractStrategy::DoReadRestart(reader, dis, cparams_ptr);
  post_setup(false, false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::initialize_mortar()
{
  // for self contact, slave and master sets may have changed,
  // thus we have to update them before initializing Dn, Mn etc.
  update_global_self_contact_state();

  // (re)setup global Mortar Core::LinAlg::SparseMatrices and Epetra_Vectors
  data().DMatrixPtr() = Teuchos::null;
  data().MMatrixPtr() = Teuchos::null;
  data().BMatrixPtr() = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      data().GSlNormalDofRowMap(), 100, false, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::assemble_mortar()
{
  // for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    interface.AssembleBMatrix(data().BMatrix());
  }

  data().BMatrix().Complete(slave_master_dof_row_map(true), data().GSlNormalDofRowMap());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::split_mortar()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  data().DMatrixPtr() = ExtractMatrix(data().BMatrix(), data().global_active_n_dof_row_map(),
      *data().global_slave_dof_row_map_ptr());

  data().MMatrixPtr() = ExtractMatrix(data().BMatrix(), data().global_active_n_dof_row_map(),
      *data().global_master_dof_row_map_ptr());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::zeroize_stiffness_state()
{
  // *** zeroize existing global matrices ***
  data().DGLmLinMatrix().reset();
  data().DGGLinMatrix().reset();

  data().DLmNWGapLinMatrix().reset();
  data().DLmTLmTMatrix().reset();
  data().DLmTLmTLinMatrix().reset();

  data().InactiveLinMatrix().reset();
  data().InactiveDDMatrix().reset();
  data().InactiveDiagMatrix().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::create_stiffness_state(const Epetra_Map& gAugInactiveSlaveDofs)
{
  // *** (re)setup global matrices ***
  data().DGLmLinMatrixPtr() = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      slave_master_dof_row_map(true), 100, false, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
  data().DGGLinMatrixPtr() = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      slave_master_dof_row_map(true), 100, false, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

  data().d_lm_nw_gap_lin_matrix_ptr() = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(data().global_active_n_dof_row_map(), 100, false, false));
  data().DLmTLmTMatrixPtr() = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(data().global_active_t_dof_row_map(), 100, false, false));
  data().DLmTLmTLinMatrixPtr() = Teuchos::rcp(
      new Core::LinAlg::SparseMatrix(data().global_active_t_dof_row_map(), 100, false, false));

  data().inactive_diag_matrix_ptr() = Core::LinAlg::CreateVector(gAugInactiveSlaveDofs, true);
  data().inactive_lin_matrix_ptr() =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(gAugInactiveSlaveDofs, 100, false, false));
  data().InactiveDDMatrixPtr() = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
      slave_dof_row_map(true), 100, false, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::zeroize_rhs_state()
{
  // *** zeroize existing global matrices ***
  data().LmN().PutScalar(0.0);
  data().AWGap().PutScalar(0.0);

  data().DLmTLmTRhs().PutScalar(0.0);
  data().InactiveRhs().PutScalar(0.0);

  data().AVec().PutScalar(0.0);
  data().KappaVec().PutScalar(0.0);
  data().WGap().PutScalar(0.0);
  data().WGapAllSlNodes().PutScalar(0.0);

  data().SlForceLm().PutScalar(0.0);
  data().SlForceLmInactive().PutScalar(0.0);
  data().SlForceG().PutScalar(0.0);
  data().MaForceLm().PutScalar(0.0);
  data().MaForceG().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::create_rhs_state(const Epetra_Map& gAugInactiveSlaveDofs)
{
  // *** (re)setup global augmented Epetra_Vectors ***
  data().LmNPtr() = Teuchos::rcp(new Epetra_Vector(data().global_active_n_dof_row_map(), true));
  data().AWGapPtr() = Teuchos::rcp(new Epetra_Vector(data().global_active_n_dof_row_map(), true));
  data().DLmTLmTRhsPtr() =
      Teuchos::rcp(new Epetra_Vector(data().global_active_t_dof_row_map(), true));
  data().InactiveRhsPtr() = Teuchos::rcp(new Epetra_Vector(gAugInactiveSlaveDofs, true));

  data().AVecPtr() = Teuchos::rcp(new Epetra_Vector(slave_row_nodes(), true));
  data().KappaVecPtr() = Teuchos::rcp(new Epetra_Vector(data().global_active_node_row_map(), true));
  data().WGapPtr() = Teuchos::rcp(new Epetra_Vector(data().global_active_n_dof_row_map(), true));
  data().WGapAllSlNodesPtr() = Teuchos::rcp(new Epetra_Vector(data().GSlNormalDofRowMap(), true));

  data().SlForceLmPtr() = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true), true));
  data().sl_force_lm_inactive_ptr() =
      Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true), true));
  data().SlForceGPtr() = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)), true);
  data().MaForceLmPtr() = Teuchos::rcp(new Epetra_Vector(master_dof_row_map(true), true));
  data().MaForceGPtr() = Teuchos::rcp(new Epetra_Vector(master_dof_row_map(true)), true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::initialize(enum Mortar::ActionType actiontype)
{
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs = Teuchos::null;
  Core::LinAlg::SplitMap(slave_dof_row_map(true), data().global_active_dof_row_map());

  switch (actiontype)
  {
    case Mortar::eval_force_stiff:
    {
      if (data().MatrixMapsValid())
      {
        zeroize_stiffness_state();
      }
      else
      {
        gAugInactiveSlaveDofs =
            Core::LinAlg::SplitMap(slave_dof_row_map(true), data().global_active_dof_row_map());
        create_stiffness_state(*gAugInactiveSlaveDofs);
        data().SetMatrixMapsValid(true);
      }
      [[fallthrough]];
    }
    case Mortar::eval_force:
    case Mortar::eval_static_constraint_rhs:
    {
      if (data().VectorMapsValid())
      {
        zeroize_rhs_state();
      }
      else
      {
        if (gAugInactiveSlaveDofs.is_null())
          gAugInactiveSlaveDofs =
              Core::LinAlg::SplitMap(slave_dof_row_map(true), data().global_active_dof_row_map());
        create_rhs_state(*gAugInactiveSlaveDofs);
        data().SetVectorMapsValid(true);
      }

      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Unsupported action type detected: %s", Mortar::ActionType2String(actiontype).c_str());
      exit(EXIT_FAILURE);
    }
  }

  if (data().is_friction())
    FOUR_C_THROW(
        "AugmentedLagrangeStrategy::Initialize: "
        "Frictional case is not yet considered!");

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_force_stiff(CONTACT::ParamsInterface& cparams)
{
  // call the evaluate force routine
  evaluate_force(cparams);

  // --- Assemble stiffness matrix ---------------------------------------
  assemble_contact_stiff();

  post_eval_force_stiff(cparams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::post_eval_force_stiff(CONTACT::ParamsInterface& cparams)
{
  if (not is_in_contact() and not was_in_contact() and not was_in_contact_last_time_step()) return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::initialize_and_evaluate_interface(
    Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr)
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // set variational approach
  cparams_ptr->set_variational_approach_type(data().variational_approach_type());
  data().PDController().setup(*cparams_ptr);

  // get type of parallel strategy
  Inpar::Mortar::ExtendGhosting strat = data().ghosting_strategy();

  // Evaluation for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    // initialize / reset interfaces
    interface.initialize();

    // store required integration time
    data().IntTime() += interface.Inttime();
    switch (strat)
    {
      /*----------------------------------------------------------*
       |  Fully redundant ghosting of master side                 |
       *----------------------------------------------------------*/
      case Inpar::Mortar::ExtendGhosting::redundant_all:
      case Inpar::Mortar::ExtendGhosting::redundant_master:
      {
        evaluate_interface(interface, 0, cparams_ptr);

        break;
      }
      default:
      {
        FOUR_C_THROW(
            "Augmented Lagrange strategy supports only redundant ghosting of interfaces.\".");
        break;
      }
    }
  }  // end interface loop

  data().PDController().check(*cparams_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::spread_global_sele_eval_times_to_interfaces()
{
  // Evaluation for all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    interface.StoreSeleEvalTimes(*data().GSeleEvalTimesPtr());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::check_parallel_distribution(const GlobalTimeMonitor& global_timer)
{
  const double my_total_time = global_timer.getMyTotalTime();
  update_parallel_distribution_status(my_total_time);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::Aug::Strategy::dyn_redistribute_contact(const Teuchos::RCP<const Epetra_Vector>& dis,
    Teuchos::RCP<const Epetra_Vector> vel, const int nlniter)
{
  return data().PDController().redistribute(dis, vel, nlniter);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_interface(CONTACT::Aug::Interface& interface,
    const int rriter, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
{
  const enum Mortar::ActionType atype = cparams_ptr->get_action_type();
  switch (atype)
  {
    case Mortar::eval_force:
    case Mortar::eval_force_stiff:
    {
      // evaluate averaged weighted gap
      interface.evaluate(rriter, cparams_ptr);

      // evaluate remaining entities and linearization
      interface.RedEvaluate(cparams_ptr);

      break;
    }
    case Mortar::eval_static_constraint_rhs:
    {
      interface.eval_active_contributions(rriter, cparams_ptr);
      interface.RedEvaluate(cparams_ptr);

      break;
    }
    case Mortar::eval_wgap_gradient_error:
    {
      interface.eval_active_contributions(rriter, cparams_ptr);

      break;
    }
    default:
    {
      FOUR_C_THROW("What shall be integrated? (enum=%d | \"%s\")", atype,
          Mortar::ActionType2String(atype).c_str());

      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::run_post_compute_x(const CONTACT::ParamsInterface& cparams,
    const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we use this routine just to store the Lagrange
   * multiplier increment. */
  Epetra_Vector zincr_exp(lm_dof_row_map(true));
  Core::LinAlg::Export(dir, zincr_exp);
  int err = zincr_exp.ReplaceMap(slave_dof_row_map(true));
  if (err) FOUR_C_THROW("ReplaceMap failed with error code %d.", err);

  // get the current step length
  const double stepLength = cparams.get_step_length();
  // ---------------------------------------------------------------------
  /* store the SCALED Lagrange multiplier increment in the contact
   * strategy */
  // ---------------------------------------------------------------------
  CATCH_EPETRA_ERROR(data().LmIncrPtr()->Update(stepLength, zincr_exp, 0.0));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::reset(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp, const Epetra_Vector& xnew)
{
  data().set_current_eval_state(Mortar::eval_none);
  CONTACT::AbstractStrategy::reset(cparams, dispnp, xnew);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::reset_lagrange_multipliers(
    const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew)
{
  /* Since the augmented Lagrangian strategy does not support any kind
   * of condensation, we do not have to check if it is a saddle
   * point system. */
  Epetra_Vector& z = *data().LmPtr();

  z.PutScalar(0.0);

  if (z.ReplaceMap(lm_dof_row_map(true))) FOUR_C_THROW("ReplaceMap failed!");

  Core::LinAlg::Export(xnew, z);

  if (z.ReplaceMap(slave_dof_row_map(true))) FOUR_C_THROW("ReplaceMap failed!");

  // ---------------------------------------------------------------------
  // store the new Lagrange multiplier in the nodes
  // ---------------------------------------------------------------------
  store_nodal_quantities(Mortar::StrategyBase::lmupdate);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::pre_eval_force(CONTACT::ParamsInterface& cparams)
{
  // set current evaluation action type
  set_current_eval_state(cparams);

  /*--------------------------------------------------------------*
   | For self-contact the master/slave sets are updated within the|
   | contact search, see SelfBinaryTree.                          |
   | Therefore, we have to initialize the mortar matrices after   |
   | interface evaluations.                                       |
   *--------------------------------------------------------------*/
  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcp(&cparams, false);
  if (is_self_contact())
  {
    // evaluate mortar terms (integrate...)
    initialize_and_evaluate_interface(cparams_ptr);
    // initialize mortar matrices and vectors
    initialize_mortar();
    // assemble mortar terms into global matrices
    assemble_mortar();
  }
  else
  {
    // initialize mortar matrices and vectors
    initialize_mortar();
    // evaluate mortar terms (integrate...)
    initialize_and_evaluate_interface(cparams_ptr);
    // assemble mortar terms into global matrices
    assemble_mortar();
  }

  if (cparams.is_predictor())
  {
    // evaluate relative movement for friction
    evaluate_rel_mov_predict();
  }
  else
    evaluate_relative_movement();

  // update active set
  update_active_set_semi_smooth(cparams);

  /* Split the Dn/Mn matrices to get only the active rows
   * (only necessary for the augmented Lagrangian formulation) */
  split_mortar();

  // initialize all rhs vectors and linearization matrices
  initialize(cparams.get_action_type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_force(CONTACT::ParamsInterface& cparams)
{
  // --- Prepare the evaluation and integrate all quantities ------------------
  pre_eval_force(cparams);

  // --- Assemble the gap vectors ---------------------------------------------
  assemble_gap();

  // --- compute the augmented forces -----------------------------------------
  evaluate_augmented_forces();

  // --- Assemble the right hand side terms -----------------------------------
  assemble_contact_rhs();

  /* Evaluate structural and constraint rhs. This is also necessary, if the
   * rhs did not change during the predictor step, but a redistribution was
   * executed! */
  evaluate_str_contact_rhs();  // update structural contact rhs
  evaluate_constr_rhs();       // update the constrRHS

  post_eval_force(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::post_eval_force(CONTACT::ParamsInterface& cparams)
{
  // Check linear and angular momentum conservation
  if (cparams.get_action_type() == Mortar::eval_force or  // only one per Newton
      cparams.get_nln_iter() == 0)                        // predictor
    check_conservation_laws(cparams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::assemble_gap()
{
  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;
  /*--------------------------------------------------------------------*
   | Assembly                                                           |
   *--------------------------------------------------------------------*
   | --> weighted gap                                                   |
   | --> averaged weighted gap                                          |
   *--------------------------------------------------------------------*/
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    interface.assemble_active_gap_vectors(data().AWGap(), data().WGap());
    interface.assemble_gap_vector_of_all_sl_nodes(data().WGapAllSlNodes());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CONTACT::Aug::Strategy::assemble_contact_rhs()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return false;
  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE CONTACT RIGHT HAND SIDE                   |
   *--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*
   | Assembly                                                           |
   *--------------------------------------------------------------------*
   | --> normal Lagrange multiplier                                     |
   | --> tangential constraint right hand side for the frictionless case|
   | --> normal and tangential inactive rhs                             |
   *--------------------------------------------------------------------*/
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    // --- augmented Lagrange formulation --------------------------------
    // --- FORCE BALANCE -------------------------------------------------
    interface.AssembleLmNVector(data().LmN());

    // --- CONSTRAINTS ---------------------------------------------------
    // active - normal direction
    // --> wGapRhs_

    // tributary area of inactive and active nodes
    interface.AssembleAugAVector(data().AVec(), data().KappaVec());

    // active - tangential direction
    interface.AssembleDLmTLmTRhs(data().DLmTLmTRhs());

    // inactive - all directions
    interface.assemble_aug_inactive_rhs(data().InactiveRhs(), data().Cn(), inactive_scale_factor());
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::assemble_contact_stiff()
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // check if contact contributions are present,
  // if not we can skip this routine to speed things up
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step()) return;

  /*--------------------------------------------------------------------*
   |             ASSEMBLE THE TANGENTIAL STIFFNESS MATRIX               |
   *--------------------------------------------------------------------*/
  // --- augmented Lagrange formulation ----------------------------------
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::Aug::Interface& interface = dynamic_cast<CONTACT::Aug::Interface&>(**cit);

    // --- Force Balance ------------------------------------------------
    // linearization w.r.t. displ.
    interface.assemble_dg_lm_lin_matrix(data().DGLmLinMatrix());
    interface.assemble_dgg_lin_matrix(data().DGGLinMatrix(), data().Cn());

    // --- Constraints --------------------------------------------------
    // linearization w.r.t. LM
    interface.assemble_d_lm_t_lm_t_matrix(data().DLmTLmTMatrix());
    interface.assemble_aug_inactive_diag_matrix(
        data().InactiveDiagMatrix(), data().Cn(), inactive_scale_factor());

    // linearization w.r.t. displ.
    // active - normal direction
    interface.assemble_d_lm_nw_gap_lin_matrix(data().DLmNWGapLinMatrix());

    // active - tangential direction
    interface.assemble_d_lm_t_lm_t_lin_matrix(data().DLmTLmTLinMatrix());

    /*--------------------------------------------------------------------*
     | The linearization of the nodal area w.r.t. the displ. for inactive |
     | nodes can help to prevent oscillations of the active set, because  |
     | the reduction of the inactive lm-value is decelerated.             |
     |                                                                    |
     | In general, the following relation does NOT hold anymore:          |
     |                    Delta(z_{n,i}^{k+1}) = - z_{n,i}^{k}            |
     |                          z_{n,i}^{k+1}  =   0                      |
     *--------------------------------------------------------------------*/
    interface.assemble_aug_inactive_lin_matrix(
        data().InactiveLinMatrix(), data().Cn(), inactive_scale_factor());

    interface.assemble_inactive_dd_matrix(
        data().InactiveDDMatrix(), data().Cn(), inactive_scale_factor());
  }

  // --- START - fill_complete matrices ----------------------------------
  // domainmap: columnmap | rangemap: rowmap
  // get inactive slave dofs
  Teuchos::RCP<Epetra_Map> gAugInactiveSlaveDofs =
      Core::LinAlg::SplitMap(slave_dof_row_map(true), data().global_active_dof_row_map());

  // --- Force Balance --------------------------------------------------
  // linearization w.r.t. displ.
  data().DGLmLinMatrix().Complete(slave_master_dof_row_map(true), slave_master_dof_row_map(true));
  data().DGGLinMatrix().Complete(slave_master_dof_row_map(true), slave_master_dof_row_map(true));

  // --- Constraints ----------------------------------------------------
  // linearization w.r.t. LM
  data().DLmTLmTMatrix().Complete(
      data().global_active_t_dof_row_map(), data().global_active_t_dof_row_map());

  // linearization w.r.t. displ.
  data().DLmNWGapLinMatrix().Complete(
      slave_master_dof_row_map(true), data().global_active_n_dof_row_map());
  data().DLmTLmTLinMatrix().Complete(slave_dof_row_map(true), data().global_active_t_dof_row_map());
  data().InactiveLinMatrix().Complete(slave_dof_row_map(true), *gAugInactiveSlaveDofs);
  data().InactiveDDMatrix().Complete(slave_dof_row_map(true), slave_dof_row_map(true));

  // --- END - fill_complete matrices ------------------------------------

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::update_active_set_semi_smooth(const CONTACT::ParamsInterface& cparams)
{
  ActiveSet active_set(*this);
  active_set.Update(cparams);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_str_contact_rhs()
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step())
  {
    data().StrContactRhsPtr() = Teuchos::null;
    return;
  }
  data().StrContactRhsPtr() = Teuchos::rcp(new Epetra_Vector(*ProblemDofs(), true));


  // For self contact, slave and master sets may have changed,
  if (is_self_contact())
    FOUR_C_THROW("Augmented Lagrange Formulation: Self contact is not yet considered!");

  // --- add contact force terms ----------------------------------------------
  // *** Slave side ***
  Epetra_Vector augfs(data().SlForceLm());
  CATCH_EPETRA_ERROR(augfs.Update(-1.0, data().SlForceG(), 1.0));
  if (data().add_inactiv_force_contributions())
    CATCH_EPETRA_ERROR(augfs.Update(1.0, data().SlForceLmInactive(), 1.0));

  Epetra_Vector augfs_exp(*ProblemDofs());
  Core::LinAlg::Export(augfs, augfs_exp);
  data().StrContactRhs().Scale(-1.0, augfs_exp);

  // Master side
  Epetra_Vector augfm(data().MaForceLm());
  CATCH_EPETRA_ERROR(augfm.Update(-1.0, data().MaForceG(), 1.0));

  Epetra_Vector augfm_exp(*ProblemDofs());
  Core::LinAlg::Export(augfm, augfm_exp);
  CATCH_EPETRA_ERROR(data().StrContactRhs().Update(-1.0, augfm_exp, 1.0));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_constr_rhs()
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step())
  {
    // (re)setup the vector
    data().ConstrRhsPtr() = Teuchos::null;
    return;
  }

  // initialize constraint r.h.s. (still with wrong map)
  Teuchos::RCP<Epetra_Vector> augConstrRhs =
      Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true), true));

  // We solve for the incremental Lagrange multiplier dz_. Hence,
  // we can keep the contact force terms on the right-hand side!
  add_contributions_to_constr_rhs(*augConstrRhs);

  // replace row map
  augConstrRhs->ReplaceMap(lm_dof_row_map(true));

  // export and set constraint rhs vector
  if (parallel_redistribution_status())
  {
    data().ConstrRhsPtr() = Teuchos::rcp(new Epetra_Vector(lm_dof_row_map(false)));
    Core::LinAlg::Export(*augConstrRhs, *data().ConstrRhsPtr());
  }
  else
    data().ConstrRhsPtr() = augConstrRhs;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_static_constraint_rhs(CONTACT::ParamsInterface& cparams)
{
  set_current_eval_state(cparams);
  initialize_and_evaluate_interface(cparams);

  // --- Assemble the gap vectors ---------------------------------------------
  initialize(cparams.get_action_type());
  assemble_gap();

  // --- Evaluate only the forces coming from the constraints -----------------
  evaluate_constraint_forces();

  // --- Zeroize all contributions from the LM values -------------------------
  zeroize_lm_forces();

  // --- Assemble the right hand side terms -----------------------------------
  assemble_contact_rhs();

  // --- update the rhs w.r.t. all constraint contributions -------------------
  evaluate_str_contact_rhs();
  evaluate_constr_rhs();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::add_contributions_to_constr_rhs(Epetra_Vector& augConstrRhs) const
{
  // Add active constraints in normal direction:
  Core::LinAlg::AssembleMyVector(0.0, augConstrRhs, 1.0, *data().WGapPtr());

  if (Core::IO::cout.requested_output_level() >= Core::IO::debug)
  {
    double wgap_nrm2 = 0.0;
    data().WGapPtr()->Norm2(&wgap_nrm2);
    Core::IO::cout << __FUNCTION__ << " [wgap-norm2] = " << wgap_nrm2 << Core::IO::endl;
  }

  // Add inactive constraints
  Core::LinAlg::AssembleMyVector(1.0, augConstrRhs, 1.0, *data().InactiveRhsPtr());

  // Add tangential frictionless constraints
  Core::LinAlg::AssembleMyVector(1.0, augConstrRhs, 1.0, *data().DLmTLmTRhsPtr());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::check_conservation_laws(CONTACT::ParamsInterface& cparams)
{
  // get forces due to Lagrange multipliers of the slave side
  const Teuchos::RCP<Epetra_Vector>& augfs_lm = data().SlForceLmPtr();

  // get forces due to the regularization term of the slave side
  const Teuchos::RCP<Epetra_Vector>& augfs_g = data().SlForceGPtr();

  // get forces due to Lagrange multipliers of the master side
  const Teuchos::RCP<Epetra_Vector>& augfm_lm = data().MaForceLmPtr();

  // get forces due to the regularization term of the master side
  const Teuchos::RCP<Epetra_Vector>& augfm_g = data().MaForceGPtr();

  // complete force on the slave side
  Epetra_Vector augfs(data().SlForceLm());
  CATCH_EPETRA_ERROR(augfs.Update(-1.0, data().SlForceG(), 1.0));
  if (data().add_inactiv_force_contributions())
    CATCH_EPETRA_ERROR(augfs.Update(1.0, data().SlForceLmInactive(), 1.0));

  // compelet force on the master side
  Epetra_Vector augfm(data().MaForceLm());
  CATCH_EPETRA_ERROR(augfm.Update(-1.0, data().MaForceG(), 1.0));

  /*-------------------------------*
   | LINEAR MOMENTUM CONSERVATION  |
   *-------------------------------*/
  if (data().print_linear_mom_conservation())
  {
    const int probdim = Dim();

    double lssum[3] = {0.0, 0.0, 0.0};  // local slave sum
    double gssum[3] = {0.0, 0.0, 0.0};  // global slave sum
    double lmsum[3] = {0.0, 0.0, 0.0};  // local master sum
    double gmsum[3] = {0.0, 0.0, 0.0};  // global master sum
    double gcsum[3] = {0.0, 0.0, 0.0};  // global complete sum
    // slave
    for (int i = 0; i < augfs_lm->MyLength(); ++i) lssum[i % probdim] += (*augfs_lm)[i];

    Comm().SumAll(lssum, gssum, probdim);
    // master
    for (int i = 0; i < augfm_lm->MyLength(); ++i) lmsum[i % probdim] += (*augfm_lm)[i];

    Comm().SumAll(lmsum, gmsum, probdim);
    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Linear Momentum Conservation      <<\n";
      std::cout << ">>  comp.-wise in x-, y- and z-direction  <<\n";
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        FOUR_C_THROW("Conservation of linear momentum is not fulfilled!");
    }

    // slave
    std::fill(lssum, lssum + probdim, 0.0);
    for (int i = 0; i < augfs_g->MyLength(); ++i) lssum[i % probdim] += (*augfs_g)[i];

    Comm().SumAll(lssum, gssum, probdim);

    // master
    std::fill(lmsum, lmsum + probdim, 0.0);
    for (int i = 0; i < augfm_g->MyLength(); ++i) lmsum[i % probdim] += (*augfm_g)[i];

    Comm().SumAll(lmsum, gmsum, probdim);

    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
      std::cout << "--------------------------------------------" << std::endl;
      std::cout << ">>      Complete                          <<" << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        FOUR_C_THROW("Conservation of linear momentum is not fulfilled!");
    }

    // slave
    std::fill(lssum, lssum + probdim, 0.0);
    for (int i = 0; i < augfs.MyLength(); ++i) lssum[i % probdim] += augfs[i];

    Comm().SumAll(lssum, gssum, probdim);

    // master
    std::fill(lmsum, lmsum + probdim, 0.0);
    for (int i = 0; i < augfm.MyLength(); ++i) lmsum[i % probdim] += augfm[i];

    Comm().SumAll(lmsum, gmsum, probdim);

    // complete balance check
    for (int i = 0; i < probdim; ++i) gcsum[i] = gssum[i] + gmsum[i];

    if (Comm().MyPID() == 0)
    {
      std::cout << "SLAVE:   " << std::setw(14) << gssum[0] << "," << std::setw(14) << gssum[1]
                << "," << std::setw(14) << gssum[2] << std::endl;

      std::cout << "MASTER:  " << std::setw(14) << gmsum[0] << "," << std::setw(14) << gmsum[1]
                << "," << std::setw(14) << gmsum[2] << std::endl;

      std::cout << "Balance: " << std::setw(14) << gcsum[0] << "," << std::setw(14) << gcsum[1]
                << "," << std::setw(14) << gcsum[2] << std::endl;
    }

    for (int i = 0; i < probdim; ++i)
    {
      if (std::abs(gcsum[i]) >
          1.0e-11 * std::max(1.0, std::max(std::abs(gssum[i]), std::abs(gmsum[i]))))
        FOUR_C_THROW("Conservation of linear momentum is not fulfilled!");
    }
  }
  /*-------------------------------*
   | ANGULAR MOMENTUM CONSERVATION |
   *-------------------------------*/
  if (data().print_angular_mom_conservation())
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n";
      std::cout << ">>      Angular Momentum Conservation     <<\n";
      std::cout << ">>  comp.-wise in x-, y- and z-direction  <<\n";
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
    }

    unsigned icount = 0;
    const int nln_iter = cparams.get_nln_iter();

    Core::LinAlg::SerialDenseMatrix conservation_data(18, 1, false);
    for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end();
         ++cit, ++icount)
    {
      const CONTACT::Interface& interface = **cit;

      if (Comm().MyPID() == 0)
      {
        std::cout << ">>----- Interface " << std::setw(2) << icount;
        std::cout << " ---------------------<<" << std::endl;
        std::cout << ">>      Standard terms (lm)               <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_lm, *augfm_lm, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.get_output_file_path(), "lm_terms");

      if (Comm().MyPID() == 0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Regularization terms (awgap)      <<" << std::endl;
      }
      interface.EvalResultantMoment(*augfs_g, *augfm_g, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.get_output_file_path(), "regularization_terms");

      if (Comm().MyPID() == 0)
      {
        std::cout << "--------------------------------------------" << std::endl;
        std::cout << ">>      Complete                          <<" << std::endl;
      }
      interface.EvalResultantMoment(augfs, augfm, &conservation_data);
      CONTACT::UTILS::WriteConservationDataToFile(Comm().MyPID(), icount, nln_iter,
          conservation_data, cparams.get_output_file_path(), "complete");
    }
    if (Comm().MyPID() == 0)
      std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<" << std::endl;
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_augmented_forces()
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step()) return;

  // augmented force terms
  enum ForceTerm
  {
    lm_term = 0,
    g_term = 1
  };

  const Epetra_Vector& zn = data().Potential().GetZnActive();
  Epetra_Vector awgapn(data().global_active_n_dof_row_map(), true);
  Core::LinAlg::ExtractMyVector(data().AWGap(), awgapn);

  // scale the averaged weighted gap elementwise by cn
  MultiplyElementwise(data().Cn(), data().global_active_node_row_map(), awgapn, false);

  std::array<double*, 2> values = {nullptr, nullptr};

  values[lm_term] = zn.Values();
  values[g_term] = awgapn.Values();

  const Epetra_MultiVector zn_awgapn(View, data().global_active_n_dof_row_map(), values.data(), 2);

  /*----------------- SLAVE SIDE ---------------------------------------------*/
  std::array<double*, 2> f_values = {nullptr, nullptr};
  f_values[lm_term] = data().SlForceLm().Values();
  f_values[g_term] = data().SlForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector slForces(View, slave_dof_row_map(true), f_values.data(), 2);

  data().DMatrix().Multiply(true, zn_awgapn, slForces);

  for (auto& inter : interface_)
  {
    const Interface& aug_inter = dynamic_cast<const Interface&>(*inter);
    aug_inter.Add_Var_A_GG(*slForces(g_term), data().Cn());
    aug_inter.assemble_sl_force_lm_inactive(
        data().SlForceLmInactive(), data().Cn(), inactive_scale_factor());
  }

  /*----------------- MASTER SIDE --------------------------------------------*/
  f_values[lm_term] = data().MaForceLm().Values();
  f_values[g_term] = data().MaForceG().Values();

  // interface forces on the slave side
  Epetra_MultiVector maForces(View, master_dof_row_map(true), f_values.data(), 2);

  data().MMatrix().Multiply(true, zn_awgapn, maForces);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_constraint_forces()
{
  if (!is_in_contact() and !was_in_contact() and !was_in_contact_last_time_step()) return;

  Epetra_Vector awgapn(data().global_active_n_dof_row_map(), true);
  Core::LinAlg::ExtractMyVector(data().AWGap(), awgapn);

  // scale the averaged weighted gap elementwise by cn
  MultiplyElementwise(data().Cn(), data().global_active_node_row_map(), awgapn, false);

  data().DMatrix().Multiply(true, awgapn, data().SlForceG());

  data().MMatrix().Multiply(true, awgapn, data().MaForceG());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::zeroize_lm_forces()
{
  // set the lagrange multiplier forces to zero
  data().SlForceLm().PutScalar(0.0);
  data().MaForceLm().PutScalar(0.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_augmented_forces(Epetra_Vector& augfs_lm,
    Epetra_Vector& augfs_g, Epetra_Vector& augfm_lm, Epetra_Vector& augfm_g) const
{
  if (!is_in_contact()) return;

  /****************** SLAVE SIDE ****************************************/
  // *** standard Lagrange multiplier fraction ***
  // Export
  Core::LinAlg::Export(*data().SlForceLmPtr(), augfs_lm);

  // *** regularization fraction ***
  // Export
  Core::LinAlg::Export(*data().SlForceGPtr(), augfs_g);
  augfs_g.Scale(-1.0);

  /****************** MASTER SIDE ***************************************/
  // *** standard lagrange multiplier fraction ***
  // Export
  Core::LinAlg::Export(*data().MaForceLmPtr(), augfm_lm);

  // *** regularization fraction ***
  // Export
  Core::LinAlg::Export(*data().MaForceGPtr(), augfm_g);
  augfm_g.Scale(-1.0);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::compute_contact_stresses()
{
  // reset contact stress class variables
  data().StressNormalPtr() = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));
  data().StressTangentialPtr() = Teuchos::rcp(new Epetra_Vector(slave_dof_row_map(true)));

  // loop over all interfaces
  for (plain_interface_set::const_iterator cit = interface_.begin(); cit != interface_.end(); ++cit)
  {
    const CONTACT::Interface& interface = **cit;

    // loop over all slave row nodes on the current interface
    for (int j = 0; j < interface.SlaveRowNodes()->NumMyElements(); ++j)
    {
      int gid = interface.SlaveRowNodes()->GID(j);
      Core::Nodes::Node* node = interface.Discret().gNode(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      Node* cnode = dynamic_cast<Node*>(node);

      // be aware of problem dimension
      int dim = Dim();
      int numdof = cnode->NumDof();
      if (dim != numdof) FOUR_C_THROW("Inconsisteny Dim <-> NumDof");

      // get nodal normal and tangential directions
      double* nn = cnode->MoData().n();
      double* nt1 = cnode->Data().txi();
      double* nt2 = cnode->Data().teta();
      double lmn = cnode->MoData().lm()[0];
      double lmt1 = cnode->MoData().lm()[1];
      double lmt2 = cnode->MoData().lm()[2];

      // find indices for DOFs of current node in Epetra_Vector
      // and put node values (normal and tangential stress components) at these DOFs

      std::vector<int> locindex(dim);

      // normal stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (data().StressNormalPtr()->Map()).LID(cnode->Dofs()[dof]);
        (*data().StressNormalPtr())[locindex[dof]] = -lmn * nn[dof];
      }

      // tangential stress components
      for (int dof = 0; dof < dim; ++dof)
      {
        locindex[dof] = (data().StressTangentialPtr()->Map()).LID(cnode->Dofs()[dof]);
        (*data().StressTangentialPtr())[locindex[dof]] = -lmt1 * nt1[dof] - lmt2 * nt2[dof];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::WriteOutput(Core::IO::DiscretizationWriter& writer) const
{
  Teuchos::RCP<Epetra_Vector> augfs_lm = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfs_g = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfm_lm = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));
  Teuchos::RCP<Epetra_Vector> augfm_g = Teuchos::rcp(new Epetra_Vector(*ProblemDofs()));

  // evaluate augmented contact forces
  evaluate_augmented_forces(*augfs_lm, *augfs_g, *augfm_lm, *augfm_g);

  // contact forces on slave and master side
  writer.write_vector("norslaveforcelm", augfs_lm);
  writer.write_vector("norslaveforceg", augfs_g);
  writer.write_vector("normasterforcelm", augfm_lm);
  writer.write_vector("normasterforceg", augfm_g);

  Epetra_Vector str_row_node_owners(*ProblemNodes(), false);
  str_row_node_owners.PutScalar(-1.0);

  for (auto& cinterface : interfaces())
  {
    const CONTACT::Aug::Interface& interface =
        dynamic_cast<const CONTACT::Aug::Interface&>(*cinterface);

    Teuchos::RCP<Epetra_Vector> irow_node_owners =
        interface.collect_row_node_owners(writer.get_discretization());

    Core::LinAlg::Export(*irow_node_owners, str_row_node_owners);
  }

  writer.write_vector(
      "contactowner", Teuchos::rcpFromRef(str_row_node_owners), Core::IO::nodevector);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> CONTACT::Aug::Strategy::get_rhs_block_ptr(
    const enum CONTACT::VecBlockType& bt) const
{
  // if there are no active contact contributions
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
    return Teuchos::null;

  // get the desired vector and return it (read-only)
  Teuchos::RCP<const Epetra_Vector> vec_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::VecBlockType::displ:
    {
      vec_ptr = data().StrContactRhsPtr();

      if (Core::IO::cout.requested_output_level() >= Core::IO::debug)
      {
        double vec_nrm2 = 0.0;
        vec_ptr->Norm2(&vec_nrm2);
        Core::IO::cout << __FUNCTION__ << " [CONTACT::VecBlockType::displ] = " << vec_nrm2
                       << Core::IO::endl;
      }

      break;
    }
    case CONTACT::VecBlockType::constraint:
    {
      vec_ptr = data().ConstrRhsPtr();
      if (Core::IO::cout.requested_output_level() >= Core::IO::debug)
      {
        double vec_nrm2 = 0.0;
        vec_ptr->Norm2(&vec_nrm2);
        Core::IO::cout << __FUNCTION__ << " [CONTACT::VecBlockType::constraint] = " << vec_nrm2
                       << Core::IO::endl;
      }

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown CONTACT::VecBlockType!");
      break;
    }
  }

  return vec_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> CONTACT::Aug::Strategy::get_matrix_block_ptr(
    const enum CONTACT::MatBlockType& bt, const ParamsInterface* cparams) const
{
  TEUCHOS_FUNC_TIME_MONITOR(CONTACT_FUNC_NAME);

  // if there are no active contact contributions
  if (!is_in_contact() && !was_in_contact() && !was_in_contact_last_time_step())
    return Teuchos::null;

  Teuchos::RCP<Core::LinAlg::SparseMatrix> mat_ptr = Teuchos::null;
  switch (bt)
  {
    case CONTACT::MatBlockType::displ_displ:
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix>& kdd_ptr = mat_ptr;
      kdd_ptr = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(slave_master_dof_row_map(true), 100, false, true));

      // build matrix kdd
      add_contributions_to_matrix_block_displ_displ(*kdd_ptr);
      kdd_ptr->Complete(slave_master_dof_row_map(true), slave_master_dof_row_map(true));

      break;
    }
    case CONTACT::MatBlockType::displ_lm:
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix>& kdz_ptr = mat_ptr;
      kdz_ptr = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(*data().global_disp_dof_row_map_ptr(), 100, false, true));

      // build constraint matrix kdz
      add_contributions_to_matrix_block_displ_lm(*kdz_ptr);
      kdz_ptr->Complete(slave_dof_row_map(true), *data().global_disp_dof_row_map_ptr());

      // transform constraint matrix kzd to lmdofmap (MatrixColTransform)
      static Teuchos::RCP<Epetra_Map> newcolmap = Teuchos::null;
      Mortar::ReplaceColumnAndDomainMap(*kdz_ptr, *lm_dof_row_map_ptr(true), &newcolmap);

      break;
    }
    case CONTACT::MatBlockType::lm_displ:
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix>& kzd_ptr = mat_ptr;
      kzd_ptr =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100, false, true));

      // build constraint matrix kzd
      add_contributions_to_matrix_block_lm_displ(*kzd_ptr);
      int err = kzd_ptr->ReplaceRowMap(*lm_dof_row_map_ptr(true));
      if (err) FOUR_C_THROW("ReplaceMap failed on kzd_ptr! (err = %d)", err);

      kzd_ptr->Complete(*data().global_disp_dof_row_map_ptr(), *lm_dof_row_map_ptr(true));

      break;
    }
    case CONTACT::MatBlockType::lm_lm:
    {
      Teuchos::RCP<Core::LinAlg::SparseMatrix>& kzz_ptr = mat_ptr;
      kzz_ptr =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(slave_dof_row_map(true), 100, false, true));

      // build constraint matrix kzz
      add_contributions_to_matrix_block_lm_lm(*kzz_ptr);

      // replace row map
      int err = kzz_ptr->ReplaceRowMap(*lm_dof_row_map_ptr(true));
      if (err) FOUR_C_THROW("ReplaceMap failed on kzz_ptr! (err = %d)", err);

      kzz_ptr->Complete(slave_dof_row_map(true), *lm_dof_row_map_ptr(true));

      // transform constraint matrix kzz to lmdofmap (columns, 2nd step)
      static Teuchos::RCP<Epetra_Map> newcolmap = Teuchos::null;
      Mortar::ReplaceColumnAndDomainMap(*kzz_ptr, *lm_dof_row_map_ptr(true), &newcolmap);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown Solid::MatBlockType!");
      break;
    }
  }

  // transform parallel row/column distribution
  // (only necessary in the parallel redistribution case)
  if (parallel_redistribution_status())
  {
    Mortar::MatrixRowColTransformer& transformer = data().matrix_row_col_transformer();
    mat_ptr = transformer.redistributed_to_unredistributed(bt, *mat_ptr);
  }

  return mat_ptr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::add_contributions_to_matrix_block_displ_displ(
    Core::LinAlg::SparseMatrix& kdd, const CONTACT::ParamsInterface* cparams) const
{
  kdd.Add(*data().DGLmLinMatrixPtr(), false, -1.0, 1.0);
  kdd.Add(*data().DGGLinMatrixPtr(), false, 1.0, 1.0);

  if (data().add_inactiv_force_contributions())
    kdd.Add(*data().InactiveDDMatrixPtr(), false, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::add_contributions_to_matrix_block_displ_lm(
    Core::LinAlg::SparseMatrix& kdz) const
{
  kdz.Add(*data().DMatrixPtr(), true, -1.0, 1.0);
  kdz.Add(*data().MMatrixPtr(), true, -1.0, 1.0);

  if (data().add_inactiv_force_contributions())
    kdz.Add(*data().inactive_lin_matrix_ptr(), true, -1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::add_contributions_to_matrix_block_lm_displ(
    Core::LinAlg::SparseMatrix& kzd) const
{
  kzd.Add(*data().d_lm_nw_gap_lin_matrix_ptr(), false, 1.0, 1.0);
  kzd.Add(*data().DLmTLmTLinMatrixPtr(), false, 1.0, 1.0);
  kzd.Add(*data().inactive_lin_matrix_ptr(), false, 1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::add_contributions_to_matrix_block_lm_lm(
    Core::LinAlg::SparseMatrix& kzz) const
{
  if (Core::LinAlg::InsertMyRowDiagonalIntoUnfilledMatrix(kzz, *data().inactive_diag_matrix_ptr()))
  {
    Epetra_Vector kzz_diag = Epetra_Vector(kzz.RangeMap(), true);
    // extract the diagonal and avoid to replace already set values
    kzz.ExtractDiagonalCopy(kzz_diag);
    Core::LinAlg::AssembleMyVector(1.0, kzz_diag, 1.0, *data().inactive_diag_matrix_ptr());

    // if the matrix is filled, we try to replace the diagonal
    if (kzz.replace_diagonal_values(kzz_diag)) FOUR_C_THROW("replace_diagonal_values failed!");
  }

  kzz.Add(*data().DLmTLmTMatrixPtr(), false, 1.0, 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Strategy::constraint_norm() const
{
  double nrm2 = 0.0;
  data().ConstrRhsPtr()->Norm2(&nrm2);
  return nrm2;
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
const Epetra_Vector& CONTACT::Aug::Strategy::get_weighted_gap(const enum MapType type) const
{
  switch (type)
  {
    case MapType::active_slave_nodes:
    {
      if (data().WGapPtr().is_null()) FOUR_C_THROW("The weighted gap vector is not initialized!");

      return *data().WGapPtr();
    }
    case MapType::all_slave_nodes:
    {
      if (data().WGapAllSlNodesPtr().is_null())
        FOUR_C_THROW("The weighted gap of all slave nodes is not initialized!");

      return *data().WGapAllSlNodesPtr();
    }
    default:
    {
      FOUR_C_THROW("Unknown weighted gap type! (type=%d)", type);
      exit(EXIT_FAILURE);
    }
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::SparseMatrix> CONTACT::Aug::Strategy::get_weighted_gap_gradient(
    const enum WGapGradientType grad_type, const enum MapType map_type) const
{
  switch (grad_type)
  {
    case WGapGradientType::force_balance:
    {
      switch (map_type)
      {
        case MapType::all_slave_nodes:
        {
          if (data().BMatrixPtr().is_null())
            FOUR_C_THROW("The modified weighted gap gradient is not initialized!");

          return data().BMatrixPtr();
        }
        case MapType::active_slave_nodes:
        {
          Teuchos::RCP<Core::LinAlg::SparseMatrix> wgap_grad = Teuchos::rcp(
              new Core::LinAlg::SparseMatrix(*data().global_active_n_dof_row_map_ptr(), 100));

          if (data().DMatrixPtr().is_null() or data().MMatrixPtr().is_null())
            FOUR_C_THROW("D-Matrix or/and M-Matrix are nullptr!");

          wgap_grad->Add(*data().DMatrixPtr(), false, 1.0, 0.0);
          wgap_grad->Add(*data().MMatrixPtr(), false, 1.0, 0.0);

          wgap_grad->Complete(*data().global_slave_master_dof_row_map_ptr(),
              *data().global_active_n_dof_row_map_ptr());

          return wgap_grad;
        }
        default:
          FOUR_C_THROW("Unsupported MapType!");
          exit(EXIT_FAILURE);
      }
      break;
    }
    case WGapGradientType::constraint_enforcement:
    {
      switch (map_type)
      {
        case MapType::all_slave_nodes:
        {
          // If the complete variational approach is active, there is no need
          // to re-assemble any quantities. Due to the fact that the system is
          // symmetric, the B-Matrix can be returned.
          if (data().variational_approach_type() == Inpar::CONTACT::var_complete)
          {
            return get_weighted_gap_gradient(
                WGapGradientType::force_balance, MapType::all_slave_nodes);
          }

          Teuchos::RCP<Core::LinAlg::SparseMatrix> wgap_grad =
              Teuchos::rcp(new Core::LinAlg::SparseMatrix(*data().g_sl_normal_dof_row_map_ptr(),
                  100, false, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

          for (plain_interface_set::const_iterator cit = interface_.begin();
               cit != interface_.end(); ++cit)
          {
            const CONTACT::Aug::Interface& interface =
                dynamic_cast<CONTACT::Aug::Interface&>(**cit);

            interface.assemble_d_lm_nw_gap_lin_matrix(*wgap_grad, map_type);
          }

          wgap_grad->Complete(
              *data().global_slave_master_dof_row_map_ptr(), *data().g_sl_normal_dof_row_map_ptr());

          return wgap_grad;
        }
        case MapType::active_slave_nodes:
        {
          if (data().d_lm_nw_gap_lin_matrix_ptr().is_null())
            FOUR_C_THROW("The active weighted gap gradient is not initialized!");

          return data().d_lm_nw_gap_lin_matrix_ptr();
        }
        default:
          FOUR_C_THROW("Unsupported MapType!");
          exit(EXIT_FAILURE);
      }
    }
    default:
      FOUR_C_THROW("Unsupported WGapGradientType!");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_weighted_gap_gradient_error(CONTACT::ParamsInterface& cparams)
{
  const int num_my_activenodes = data().global_active_node_row_map_ptr()->NumMyElements();
  const int num_gactivenodes = data().global_active_node_row_map_ptr()->NumGlobalElements();
  if (num_gactivenodes == 0) return;

  std::unordered_map<int, Deriv1stMap> grad_error_ma;
  grad_error_ma.reserve(num_my_activenodes);
  std::unordered_map<int, Deriv1stMap> grad_error_jac;
  grad_error_jac.reserve(num_my_activenodes);

  cparams.SetUnorderedMap(&grad_error_ma, 0);
  cparams.SetUnorderedMap(&grad_error_jac, 1);

  Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr = Teuchos::rcpFromRef(cparams);
  initialize_and_evaluate_interface(cparams_ptr);

  // --- Nodal gradient error due to the convective parametric master coordinate
  {
    std::vector<std::pair<int, double>>& grad_error_ma_per_node =
        data().nodal_gradient_error_ma_proj();

    L2ErrorNormPerNode(grad_error_ma, grad_error_ma_per_node);
    Core::IO::cout(Core::IO::standard) << "Nodal gradient error: projected parametric "
                                          "master coordinate\n";
    for (auto& ge_per_node : grad_error_ma_per_node)
      Core::IO::cout(Core::IO::standard)
          << "GID #" << ge_per_node.first << ", e = " << ge_per_node.second << Core::IO::endl;
  }

  // --- Nodal gradient error due to the slave jacobian determinant
  {
    std::vector<std::pair<int, double>>& grad_error_jac_per_node =
        data().nodal_gradient_error_jacobian();

    L2ErrorNormPerNode(grad_error_jac, grad_error_jac_per_node);
    Core::IO::cout(Core::IO::standard) << "Nodal gradient error: jacobian\n";
    for (auto& ge_per_node : grad_error_jac_per_node)
      Core::IO::cout(Core::IO::standard)
          << "GID #" << ge_per_node.first << ", e = " << ge_per_node.second << Core::IO::endl;
  }

  // --- total gradient error (master proj + slave jacobian)
  {
    std::unordered_map<int, Deriv1stMap>* error_map_vec[2] = {&grad_error_ma, &grad_error_jac};
    double my_total_error = MyTotalSquareError(error_map_vec, 2);
    double& total_error = data().total_gradient_error();
    Comm().SumAll(&my_total_error, &total_error, 1);

    total_error /= num_gactivenodes;
    total_error = std::sqrt(total_error);

    Core::IO::cout(Core::IO::standard) << "total_error = " << total_error << Core::IO::endl;
  }

  cparams.ClearAll(Core::Gen::AnyDataContainer::DataType::unordered_map);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::MyTotalSquareError(
    const std::unordered_map<int, Deriv1stMap>* const* error_map_vec, const unsigned num_vecs)
{
  double total_error = 0.0;
  for (unsigned i = 0; i < num_vecs; ++i)
    for (auto& error_map_i : (*error_map_vec[i]))
      for (auto& error_map_ij : error_map_i.second)
        total_error += error_map_ij.second * error_map_ij.second;

  return total_error;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::L2ErrorNormPerNode(const std::unordered_map<int, Deriv1stMap>& error_map,
    std::vector<std::pair<int, double>>& error_norm_per_node)
{
  error_norm_per_node.clear();
  error_norm_per_node.resize(error_map.size(), std::pair<int, double>());

  unsigned count = 0;
  for (auto& error_map_i : error_map)
  {
    std::pair<int, double>& nodal_error_norm = error_norm_per_node[count++];
    nodal_error_norm.first = error_map_i.first;
    for (auto& error_map_ij : error_map_i.second)
      nodal_error_norm.second += error_map_ij.second * error_map_ij.second;

    // compute L2-norm for each nodal contribution
    nodal_error_norm.second = std::sqrt(nodal_error_norm.second);
    //    std::cout << "nodal_error_norm of #" << nodal_error_norm.first << " = " <<
    //        nodal_error_norm.second << std::endl;
  }

  if (count != error_map.size()) FOUR_C_THROW("Size <--> count mismatch!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Strategy::get_potential_value(
    const enum NOX::Nln::MeritFunction::MeritFctName mrt_type) const
{
  // Since constness is broken in this implementation, we resort to a const-cast here.
  auto& pot = const_cast<Potential&>(data().Potential());
  pot.Compute();

  switch (mrt_type)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    {
      return data().Potential().get(POTENTIAL::Type::lagrangian, POTENTIAL::SetType::all);
    }
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    {
      return data().Potential().get(POTENTIAL::Type::lagrangian, POTENTIAL::SetType::active);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    {
      return data().Potential().get(
          POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::all);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      return data().Potential().get(
          POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::active);
    }
    default:
    {
      FOUR_C_THROW(
          "The specified merit function type is not yet supported. "
          "type = %s | %d",
          NOX::Nln::MeritFunction::MeritFuncName2String(mrt_type).c_str(), mrt_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CONTACT::Aug::Strategy::get_linearized_potential_value_terms(const Epetra_Vector& dir,
    const enum NOX::Nln::MeritFunction::MeritFctName mrt_type,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  switch (mrt_type)
  {
    case NOX::Nln::MeritFunction::mrtfct_lagrangian:
    {
      return get_linearized_potential_model_terms(
          dir, POTENTIAL::Type::lagrangian, POTENTIAL::SetType::all, linorder, lintype);
    }
    case NOX::Nln::MeritFunction::mrtfct_lagrangian_active:
    {
      return get_linearized_potential_model_terms(
          dir, POTENTIAL::Type::lagrangian, POTENTIAL::SetType::active, linorder, lintype);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm:
    {
      return get_linearized_potential_model_terms(
          dir, POTENTIAL::Type::infeasibility_measure, POTENTIAL::SetType::all, linorder, lintype);
    }
    case NOX::Nln::MeritFunction::mrtfct_infeasibility_two_norm_active:
    {
      return get_linearized_potential_model_terms(dir, POTENTIAL::Type::infeasibility_measure,
          POTENTIAL::SetType::active, linorder, lintype);
    }
    default:
    {
      FOUR_C_THROW(
          "The specified merit function type is not yet supported. "
          "type = %s | %d",
          NOX::Nln::MeritFunction::MeritFuncName2String(mrt_type).c_str(), mrt_type);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Strategy::get_linearized_potential_model_terms(const Epetra_Vector& dir,
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::Nln::MeritFunction::LinOrder linorder,
    const enum NOX::Nln::MeritFunction::LinType lintype) const
{
  double linval = 0.0;

  // Since constness is broken in this implementation, we resort to a const-cast here.
  auto& pot = const_cast<Potential&>(data().Potential());
  pot.ComputeLin(dir);

  switch (linorder)
  {
    case NOX::Nln::MeritFunction::linorder_first:
    {
      get_linearized_potential_model_terms_1st_order(pottype, potset, lintype, linval);

      break;
    }
    case NOX::Nln::MeritFunction::linorder_second:
    {
      get_linearized_potential_model_terms_2nd_order(pottype, potset, lintype, linval);

      break;
    }
    case NOX::Nln::MeritFunction::linorder_all:
    {
      get_linearized_potential_model_terms_1st_order(pottype, potset, lintype, linval);
      get_linearized_potential_model_terms_2nd_order(pottype, potset, lintype, linval);

      break;
    }
  }

  return linval;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::get_linearized_potential_model_terms_1st_order(
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::Nln::MeritFunction::LinType lintype, double& linval) const
{
  const CONTACT::Aug::Potential& pot = data().Potential();

  switch (lintype)
  {
    case NOX::Nln::MeritFunction::lin_wrt_primary_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_d);

      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_lagrange_multiplier_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_z);

      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_mixed_dofs:
    {
      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_all_dofs:
    {
      get_linearized_potential_model_terms_1st_order(
          pottype, potset, NOX::Nln::MeritFunction::lin_wrt_primary_dofs, linval);
      get_linearized_potential_model_terms_1st_order(
          pottype, potset, NOX::Nln::MeritFunction::lin_wrt_lagrange_multiplier_dofs, linval);

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::get_linearized_potential_model_terms_2nd_order(
    const enum POTENTIAL::Type pottype, const enum POTENTIAL::SetType potset,
    const enum NOX::Nln::MeritFunction::LinType lintype, double& linval) const
{
  const CONTACT::Aug::Potential& pot = data().Potential();

  switch (lintype)
  {
    case NOX::Nln::MeritFunction::lin_wrt_primary_dofs:
    {
      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_lagrange_multiplier_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_z_and_z);

      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_mixed_dofs:
    {
      linval += pot.GetLin(pottype, potset, POTENTIAL::LinTerm::wrt_d_and_z);

      break;
    }
    case NOX::Nln::MeritFunction::lin_wrt_all_dofs:
    {
      get_linearized_potential_model_terms_2nd_order(
          pottype, potset, NOX::Nln::MeritFunction::lin_wrt_primary_dofs, linval);
      get_linearized_potential_model_terms_2nd_order(
          pottype, potset, NOX::Nln::MeritFunction::lin_wrt_lagrange_multiplier_dofs, linval);
      get_linearized_potential_model_terms_2nd_order(
          pottype, potset, NOX::Nln::MeritFunction::lin_wrt_mixed_dofs, linval);

      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::split_state_vector(const Epetra_Vector& full_state,
    Teuchos::RCP<Epetra_Vector>& displ_state_slma_ptr,
    Teuchos::RCP<Epetra_Vector>& z_state_active_ptr,
    Teuchos::RCP<Epetra_Vector>& z_state_inactive_ptr) const
{
  if (displ_state_slma_ptr.is_null())
    displ_state_slma_ptr =
        Teuchos::rcp(new Epetra_Vector(*data().global_slave_master_dof_row_map_ptr()));

  if (z_state_active_ptr.is_null())
    z_state_active_ptr = Teuchos::rcp(new Epetra_Vector(*data().global_active_n_dof_row_map_ptr()));

  if (z_state_inactive_ptr.is_null())
  {
    Teuchos::RCP<Epetra_Map> ginactivendofs = Core::LinAlg::SplitMap(
        *data().g_sl_normal_dof_row_map_ptr(), *data().global_active_n_dof_row_map_ptr());
    z_state_inactive_ptr = Teuchos::rcp(new Epetra_Vector(*ginactivendofs));
  }

  split_state_vector(full_state, *displ_state_slma_ptr, *z_state_active_ptr, *z_state_inactive_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::split_state_vector(const Epetra_Vector& full_state,
    Epetra_Vector& displ_state_slma, Epetra_Vector& z_state_active,
    Epetra_Vector& z_state_inactive) const
{
  // extract slave/master part of the displacement increment
  Epetra_Vector displ_exp(*data().global_disp_dof_row_map_ptr());
  Core::LinAlg::Export(full_state, displ_exp);
  Core::LinAlg::ExtractMyVector(displ_exp, displ_state_slma);

  // extract active/inactive part of the solution vector
  Epetra_Vector z_exp(lm_dof_row_map(true));
  Core::LinAlg::Export(full_state, z_exp);
  CATCH_EPETRA_ERROR(z_exp.ReplaceMap(slave_dof_row_map(true)));

  Core::LinAlg::ExtractMyVector(z_exp, z_state_active);
  Core::LinAlg::ExtractMyVector(z_exp, z_state_inactive);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::Aug::Strategy::characteristic_interface_element_length(
    const enum CONTACT::Aug::SideType stype) const
{
  //  find the maximal characteristic interface element length
  double my_max_ih = -1.0;

  for (const Teuchos::RCP<CONTACT::Interface>& iptr : interface_)
  {
    const Interface& interface = dynamic_cast<const Interface&>(*iptr);
    my_max_ih = std::max(interface.my_characteristic_element_length(stype), my_max_ih);
  }

  double max_ih = -1;
  Comm().MaxAll(&my_max_ih, &max_ih, 1);

  return max_ih;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::set_current_eval_state(const CONTACT::ParamsInterface& cparams)
{
  data().set_current_eval_state(cparams.get_action_type());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Aug::Strategy::evaluate_reference_state()
{
  // do nothing for the augmented strategy
}

FOUR_C_NAMESPACE_CLOSE
