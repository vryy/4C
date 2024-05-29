/*----------------------------------------------------------------------*/
/*! \file
 \brief Evaluates flux between ScaTra and ScaTra on manifolds incl. coupling matrices

 \level 2


 *------------------------------------------------------------------------------------------------*/

#include "4C_ssi_manifold_utils.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_comm_utils_gid_vector.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_assemblestrategy.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_s2i.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_ssi_monolithic.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::ManifoldScaTraCoupling::ManifoldScaTraCoupling(Teuchos::RCP<DRT::Discretization> manifolddis,
    Teuchos::RCP<DRT::Discretization> scatradis, CORE::Conditions::Condition* condition_manifold,
    CORE::Conditions::Condition* condition_kinetics, const int ndof_per_node)
    : condition_kinetics_(condition_kinetics),
      condition_manifold_(condition_manifold),
      coupling_adapter_(Teuchos::rcp(new CORE::ADAPTER::Coupling())),
      inv_thickness_(1.0 / condition_manifold->parameters().Get<double>("thickness")),
      manifold_condition_id_(condition_manifold->parameters().Get<int>("ConditionID")),
      kinetics_condition_id_(condition_kinetics->parameters().Get<int>("ConditionID")),
      manifold_map_extractor_(Teuchos::null),
      master_converter_(Teuchos::null),
      scatra_map_extractor_(Teuchos::null),
      size_matrix_graph_()
{
  std::vector<int> inodegidvec_manifold;
  CORE::COMM::AddOwnedNodeGIDFromList(
      *manifolddis, *condition_manifold->GetNodes(), inodegidvec_manifold);

  std::vector<int> inodegidvec_scatra;
  CORE::COMM::AddOwnedNodeGIDFromList(
      *scatradis, *condition_kinetics->GetNodes(), inodegidvec_scatra);

  coupling_adapter_->setup_coupling(*scatradis, *manifolddis, inodegidvec_scatra,
      inodegidvec_manifold, ndof_per_node, true, 1.0e-8);
  master_converter_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMasterConverter(*coupling_adapter_));

  scatra_map_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(
      *scatradis->dof_row_map(), coupling_adapter_->MasterDofMap(), true));

  manifold_map_extractor_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(
      *manifolddis->dof_row_map(), coupling_adapter_->SlaveDofMap(), true));

  // initially, the matrices are empty
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::ManifoldScaTra, 0));
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::ManifoldStructure, 0));
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::ScaTraManifold, 0));
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::ScaTraStructure, 0));
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::SysMatManifold, 0));
  size_matrix_graph_.insert(std::make_pair(BlockMatrixType::SysMatScaTra, 0));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
bool SSI::ManifoldScaTraCoupling::check_and_set_size_of_matrix_graph(
    const BlockMatrixType block, const int size)
{
  // check, if size of matrix graph changed between last evaluation and this evaluation
  const bool changed_size = size != size_matrix_graph_.at(block);

  // update new size
  size_matrix_graph_.at(block) = size;

  return changed_size;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::ScaTraManifoldScaTraFluxEvaluator::ScaTraManifoldScaTraFluxEvaluator(
    const SSI::SsiMono& ssi_mono)
    : block_map_scatra_(ssi_mono.BlockMapScaTra()),
      block_map_scatra_manifold_(ssi_mono.block_map_sca_tra_manifold()),
      block_map_structure_(ssi_mono.BlockMapStructure()),
      do_output_(CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->SSIControlParams().sublist("MANIFOLD"), "OUTPUT_INFLOW")),
      full_map_manifold_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold))),
      full_map_scatra_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport))),
      full_map_structure_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::structure))),
      scatra_(ssi_mono.ScaTraBaseAlgorithm()),
      scatra_manifold_(ssi_mono.sca_tra_manifold_base_algorithm()),
      ssi_structure_meshtying_(ssi_mono.ssi_structure_mesh_tying())
{
  // safety check before setup of coupling
  if (ssi_mono.ScaTraField()->NumDofPerNode() != ssi_mono.ScaTraManifold()->NumDofPerNode())
    FOUR_C_THROW("Number of dofs per node of scatra field and scatra manifold field must be equal");

  std::vector<CORE::Conditions::Condition*> conditions_manifold;
  scatra_manifold_->ScaTraField()->discretization()->GetCondition(
      "SSISurfaceManifold", conditions_manifold);

  std::vector<CORE::Conditions::Condition*> conditions_manifold_kinetics_scatra;
  scatra_->ScaTraField()->discretization()->GetCondition(
      "SSISurfaceManifoldKinetics", conditions_manifold_kinetics_scatra);

  // create pair: manifold condition - kinetics condition
  for (const auto& condition_manifold : conditions_manifold)
  {
    for (const auto& condition_kinetics : conditions_manifold_kinetics_scatra)
    {
      if (condition_manifold->parameters().Get<int>("ConditionID") ==
          condition_kinetics->parameters().Get<int>("ManifoldConditionID"))
      {
        scatra_manifold_couplings_.emplace_back(Teuchos::rcp(
            new SSI::ManifoldScaTraCoupling(scatra_manifold_->ScaTraField()->discretization(),
                scatra_->ScaTraField()->discretization(), condition_manifold, condition_kinetics,
                ssi_mono.ScaTraManifold()->NumDofPerNode())));
      }
    }
  }

  rhs_manifold_ = CORE::LINALG::CreateVector(*full_map_manifold_, true);
  rhs_scatra_ = CORE::LINALG::CreateVector(*full_map_scatra_, true);

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      systemmatrix_manifold_ = SSI::UTILS::SSIMatrices::setup_block_matrix(
          block_map_scatra_manifold_, block_map_scatra_manifold_);
      systemmatrix_scatra_ =
          SSI::UTILS::SSIMatrices::setup_block_matrix(block_map_scatra_, block_map_scatra_);
      matrix_manifold_structure_ = SSI::UTILS::SSIMatrices::setup_block_matrix(
          block_map_scatra_manifold_, block_map_structure_);
      matrix_manifold_scatra_ = SSI::UTILS::SSIMatrices::setup_block_matrix(
          block_map_scatra_manifold_, block_map_scatra_);
      matrix_scatra_manifold_ = SSI::UTILS::SSIMatrices::setup_block_matrix(
          block_map_scatra_, block_map_scatra_manifold_);
      matrix_scatra_structure_ =
          SSI::UTILS::SSIMatrices::setup_block_matrix(block_map_scatra_, block_map_structure_);

      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      systemmatrix_manifold_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
      systemmatrix_scatra_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);
      matrix_manifold_structure_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
      matrix_manifold_scatra_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
      matrix_scatra_manifold_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);
      matrix_scatra_structure_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);

      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  rhs_manifold_cond_ = CORE::LINALG::CreateVector(*full_map_manifold_, true);
  rhs_scatra_cond_ = CORE::LINALG::CreateVector(*full_map_scatra_, true);

  systemmatrix_manifold_cond_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
  systemmatrix_scatra_cond_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);
  matrix_manifold_scatra_cond_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
  matrix_manifold_structure_cond_ =
      SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_manifold_);
  matrix_scatra_manifold_cond_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);
  matrix_scatra_structure_cond_ = SSI::UTILS::SSIMatrices::setup_sparse_matrix(full_map_scatra_);

  // Prepare runtime csv writer
  if (DoOutput())
  {
    runtime_csvwriter_.emplace(ssi_mono.Comm().MyPID(),
        *GLOBAL::Problem::Instance()->OutputControlFile(), "manifold_inflow");

    for (const auto& condition_manifold : conditions_manifold)
    {
      const std::string manifold_string =
          "manifold " + std::to_string(condition_manifold->parameters().Get<int>("ConditionID"));

      runtime_csvwriter_->register_data_vector("Integral of " + manifold_string, 1, 16);

      for (int k = 0; k < ssi_mono.ScaTraManifold()->NumDofPerNode(); ++k)
      {
        runtime_csvwriter_->register_data_vector(
            "Total flux of scalar " + std::to_string(k + 1) + " into " + manifold_string, 1, 16);

        runtime_csvwriter_->register_data_vector(
            "Mean flux of scalar " + std::to_string(k + 1) + " into " + manifold_string, 1, 16);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_matrix_manifold_sca_tra()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      matrix_manifold_scatra_->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      matrix_manifold_scatra_->Complete(*full_map_scatra_, *full_map_manifold_);
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_matrix_manifold_structure()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      matrix_manifold_structure_->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      matrix_manifold_structure_->Complete(*full_map_structure_, *full_map_manifold_);
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_matrix_sca_tra_manifold()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      matrix_scatra_manifold_->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      matrix_scatra_manifold_->Complete(*full_map_manifold_, *full_map_scatra_);
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_matrix_sca_tra_structure()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      matrix_scatra_structure_->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      matrix_scatra_structure_->Complete(*full_map_structure_, *full_map_scatra_);
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_system_matrix_manifold()
{
  systemmatrix_manifold_->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::complete_system_matrix_sca_tra()
{
  systemmatrix_scatra_->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::Evaluate()
{
  // clear matrices and rhs from last evaluation
  systemmatrix_manifold_->Zero();
  systemmatrix_scatra_->Zero();
  matrix_manifold_scatra_->Zero();
  matrix_manifold_structure_->Zero();
  matrix_scatra_manifold_->Zero();
  matrix_scatra_structure_->Zero();

  rhs_manifold_->PutScalar(0.0);
  rhs_scatra_->PutScalar(0.0);

  // evaluate all scatra-manifold coupling conditions
  for (const auto& scatra_manifold_coupling : scatra_manifold_couplings_)
  {
    // clear matrices and rhs from last condition. Maps are different for each condition (need for
    // UnComplete()).
    systemmatrix_manifold_cond_->Zero();
    systemmatrix_manifold_cond_->UnComplete();
    systemmatrix_scatra_cond_->Zero();
    systemmatrix_scatra_cond_->UnComplete();
    matrix_manifold_scatra_cond_->Zero();
    matrix_manifold_scatra_cond_->UnComplete();
    matrix_manifold_structure_cond_->Zero();
    matrix_manifold_structure_cond_->UnComplete();
    matrix_scatra_manifold_cond_->Zero();
    matrix_scatra_manifold_cond_->UnComplete();
    matrix_scatra_structure_cond_->Zero();
    matrix_scatra_structure_cond_->UnComplete();

    rhs_manifold_cond_->PutScalar(0.0);
    rhs_scatra_cond_->PutScalar(0.0);

    evaluate_bulk_side(scatra_manifold_coupling);

    copy_sca_tra_sca_tra_manifold_side(scatra_manifold_coupling);

    // This is needed because the graph of the matrices could change from step to step in case we
    // have zero flux (and then zero entries in the matrices)
    un_complete_matrices_if_necessary(scatra_manifold_coupling);

    add_condition_contribution();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::evaluate_bulk_side(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  // First: Set parameters to elements
  pre_evaluate(scatra_manifold_coupling);

  scatra_->ScaTraField()->discretization()->set_state("phinp", scatra_->ScaTraField()->Phinp());

  // Second: Evaluate condition
  {
    // manifold-scatra coupling matrix evaluated on scatra side
    auto matrix_scatra_manifold_cond_on_scatra_side =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));

    Teuchos::ParameterList condparams;
    CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
        "action", SCATRA::BoundaryAction::calc_s2icoupling, condparams);
    condparams.set<int>("evaluate_manifold_coupling", 1);

    scatra_->ScaTraField()->add_time_integration_specific_vectors();

    // Evaluation of RHS and scatra-manifold coupling matrices
    {
      CORE::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
          "differentiationtype", SCATRA::DifferentiationType::elch, condparams);

      // dscatra_dscatra, dscatra_dmanifold (on scatra side)
      CORE::FE::AssembleStrategy strategyscatra(0, 0, systemmatrix_scatra_cond_,
          matrix_scatra_manifold_cond_on_scatra_side, rhs_scatra_cond_, Teuchos::null,
          Teuchos::null);

      scatra_->ScaTraField()->discretization()->evaluate_condition(condparams, strategyscatra,
          "SSISurfaceManifoldKinetics", scatra_manifold_coupling->KineticsConditionID());

      systemmatrix_scatra_cond_->Complete();
      matrix_scatra_manifold_cond_on_scatra_side->Complete();

      // dscatra_dmanifold (on scatra side) -> dscatra_dmanifold
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_manifold_cond_on_scatra_side,
          *full_map_scatra_, *full_map_scatra_, 1.0, nullptr,
          &*scatra_manifold_coupling->MasterConverter(), *matrix_scatra_manifold_cond_, true, true);
      matrix_scatra_manifold_cond_->Complete(*full_map_manifold_, *full_map_scatra_);
    }

    // Evaluation of linearization w.r.t. displacement
    {
      CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
          "action", SCATRA::BoundaryAction::calc_s2icoupling_od, condparams);

      CORE::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
          "differentiationtype", SCATRA::DifferentiationType::disp, condparams);

      // dscatra_dstructure
      auto matrix_scatra_structure_cond_slave_side_disp_evaluate =
          Teuchos::rcp(new CORE::LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));

      CORE::FE::AssembleStrategy strategyscatra(0, 1,
          matrix_scatra_structure_cond_slave_side_disp_evaluate, Teuchos::null, Teuchos::null,
          Teuchos::null, Teuchos::null);

      scatra_->ScaTraField()->discretization()->evaluate_condition(condparams, strategyscatra,
          "SSISurfaceManifoldKinetics", scatra_manifold_coupling->KineticsConditionID());

      matrix_scatra_structure_cond_slave_side_disp_evaluate->Complete(
          *full_map_structure_, *full_map_scatra_);

      // "slave side" from manifold and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      auto matrix_scatra_structure_cond_slave_side_disp =
          Teuchos::rcp(new CORE::LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->slave_slave_transformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            CORE::ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        CORE::LINALG::MatrixLogicalSplitAndTransform()(
            *matrix_scatra_structure_cond_slave_side_disp_evaluate, *full_map_scatra_, *slave_map,
            1.0, nullptr, &slave_slave_converter, *matrix_scatra_structure_cond_slave_side_disp,
            true, true);
      }
      matrix_scatra_structure_cond_slave_side_disp->Complete(
          *full_map_structure_, *full_map_scatra_);

      // Add slave side disp. contributions
      matrix_scatra_structure_cond_->Add(
          *matrix_scatra_structure_cond_slave_side_disp, false, 1.0, 0.0);

      // Add master side disp. contributions
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto converter = meshtying->SlaveSideConverter();

        // assemble derivatives of x w.r.t. structure slave dofs
        CORE::LINALG::MatrixLogicalSplitAndTransform()(
            *matrix_scatra_structure_cond_slave_side_disp,
            matrix_scatra_structure_cond_slave_side_disp->RangeMap(), *cond_slave_dof_map, 1.0,
            nullptr, &(*converter), *matrix_scatra_structure_cond_, true, true);
      }
      matrix_scatra_structure_cond_->Complete(*full_map_structure_, *full_map_scatra_);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::copy_sca_tra_sca_tra_manifold_side(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const double inv_thickness = scatra_manifold_coupling->InvThickness();
  {
    auto rhs_scatra_cond_extract =
        scatra_manifold_coupling->ScaTraMapExtractor()->ExtractCondVector(rhs_scatra_cond_);

    auto rhs_manifold_cond_extract =
        scatra_manifold_coupling->CouplingAdapter()->MasterToSlave(rhs_scatra_cond_extract);

    scatra_manifold_coupling->manifold_map_extractor()->AddCondVector(
        rhs_manifold_cond_extract, rhs_manifold_cond_);
    rhs_manifold_cond_->Scale(-inv_thickness);
  }

  // dmanifold_dscatra: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_scatra_cond_, *full_map_scatra_,
      *full_map_scatra_, -inv_thickness, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *matrix_manifold_scatra_cond_, true, true);

  matrix_manifold_scatra_cond_->Complete(*full_map_scatra_, *full_map_manifold_);

  // dmanifold_dmanifold: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_manifold_cond_, *full_map_scatra_,
      *full_map_manifold_, -inv_thickness, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *systemmatrix_manifold_cond_, true, true);

  systemmatrix_manifold_cond_->Complete();

  // dmanifold_dstructure: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_structure_cond_, *full_map_scatra_,
      *full_map_structure_, -inv_thickness, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *matrix_manifold_structure_cond_, true, true);

  matrix_manifold_structure_cond_->Complete(*full_map_structure_, *full_map_manifold_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::add_condition_contribution()
{
  rhs_manifold_->Update(1.0, *rhs_manifold_cond_, 1.0);
  rhs_scatra_->Update(1.0, *rhs_scatra_cond_, 1.0);

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      auto blockmaps_manifold = *scatra_manifold_->ScaTraField()->BlockMaps();

      /*
       * m: manifold
       * s: scatra
       * d: structure/displacements
       */
      auto flux_manifold_scatra_mm_block =
          systemmatrix_manifold_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              blockmaps_manifold, blockmaps_manifold);
      auto flux_manifold_scatra_md_block =
          matrix_manifold_structure_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, blockmaps_manifold);
      auto flux_manifold_scatra_ms_block =
          matrix_manifold_scatra_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_scatra_, blockmaps_manifold);
      auto flux_manifold_scatra_sm_block =
          matrix_scatra_manifold_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              blockmaps_manifold, *block_map_scatra_);
      auto flux_manifold_scatra_sd_block =
          matrix_scatra_structure_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, *block_map_scatra_);
      auto flux_manifold_scatra_ss_block =
          systemmatrix_scatra_cond_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_scatra_, *block_map_scatra_);

      flux_manifold_scatra_mm_block->Complete();
      flux_manifold_scatra_md_block->Complete();
      flux_manifold_scatra_ms_block->Complete();
      flux_manifold_scatra_sm_block->Complete();
      flux_manifold_scatra_sd_block->Complete();
      flux_manifold_scatra_ss_block->Complete();

      systemmatrix_manifold_->Add(*flux_manifold_scatra_mm_block, false, 1.0, 1.0);
      matrix_manifold_scatra_->Add(*flux_manifold_scatra_ms_block, false, 1.0, 1.0);
      matrix_manifold_structure_->Add(*flux_manifold_scatra_md_block, false, 1.0, 1.0);
      systemmatrix_scatra_->Add(*flux_manifold_scatra_ss_block, false, 1.0, 1.0);
      matrix_scatra_manifold_->Add(*flux_manifold_scatra_sm_block, false, 1.0, 1.0);
      matrix_scatra_structure_->Add(*flux_manifold_scatra_sd_block, false, 1.0, 1.0);

      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      systemmatrix_manifold_->Add(*systemmatrix_manifold_cond_, false, 1.0, 1.0);
      matrix_manifold_structure_->Add(*matrix_manifold_structure_cond_, false, 1.0, 1.0);
      matrix_manifold_scatra_->Add(*matrix_manifold_scatra_cond_, false, 1.0, 1.0);
      systemmatrix_scatra_->Add(*systemmatrix_scatra_cond_, false, 1.0, 1.0);
      matrix_scatra_structure_->Add(*matrix_scatra_structure_cond_, false, 1.0, 1.0);
      matrix_scatra_manifold_->Add(*matrix_scatra_manifold_cond_, false, 1.0, 1.0);

      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::evaluate_sca_tra_manifold_inflow()
{
  inflow_.clear();
  domainintegral_.clear();

  scatra_->ScaTraField()->discretization()->set_state("phinp", scatra_->ScaTraField()->Phinp());

  for (const auto& scatra_manifold_coupling : scatra_manifold_couplings_)
  {
    const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

    std::vector<double> zero_scalar_vector(scatra_->ScaTraField()->NumDofPerNode(), 0.0);
    inflow_.insert(std::make_pair(kineticsID, zero_scalar_vector));

    // First: set parameters to elements
    pre_evaluate(scatra_manifold_coupling);

    // Second: evaluate condition
    evaluate_sca_tra_manifold_inflow_integral(scatra_manifold_coupling);

    // Third: evaluate domain integral
    evaluate_sca_tra_manifold_domain_integral(scatra_manifold_coupling);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::evaluate_sca_tra_manifold_domain_integral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

  // integrate only if not done so far
  if (domainintegral_.find(kineticsID) == domainintegral_.end())
  {
    Teuchos::ParameterList condparams;

    CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
        "action", SCATRA::BoundaryAction::calc_boundary_integral, condparams);

    // integrated domain of this condition
    auto domainintegral_cond = Teuchos::rcp(new CORE::LINALG::SerialDenseVector(1));

    scatra_->ScaTraField()->discretization()->EvaluateScalars(
        condparams, domainintegral_cond, "SSISurfaceManifold", kineticsID);

    domainintegral_.insert(std::make_pair(kineticsID, domainintegral_cond->values()[0]));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::evaluate_sca_tra_manifold_inflow_integral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

  Teuchos::ParameterList condparams;

  CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_flux, condparams);

  condparams.set<bool>("only_positive_fluxes", true);

  condparams.set<int>("evaluate_manifold_coupling", 1);

  // integrated scalars of this condition
  auto inflow_cond =
      Teuchos::rcp(new CORE::LINALG::SerialDenseVector(scatra_->ScaTraField()->NumDofPerNode()));

  scatra_->ScaTraField()->discretization()->EvaluateScalars(
      condparams, inflow_cond, "SSISurfaceManifold", kineticsID);

  for (int i = 0; i < inflow_cond->length(); ++i)
    inflow_.at(kineticsID).at(i) += inflow_cond->values()[i];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::pre_evaluate(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  Teuchos::ParameterList eleparams;

  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::set_scatra_ele_boundary_parameter, eleparams);

  eleparams.set<CORE::Conditions::ConditionType>(
      "condition type", CORE::Conditions::ConditionType::S2IKinetics);

  switch (scatra_manifold_coupling->ConditionKinetics()->parameters().Get<int>("kinetic model"))
  {
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_constantinterfaceresistance);
      eleparams.set<double>("resistance",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<double>("resistance"));
      eleparams.set<const std::vector<int>*>("onoff",
          &scatra_manifold_coupling->ConditionKinetics()->parameters().Get<std::vector<int>>(
              "onoff"));
      eleparams.set<int>("numelectrons",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<int>("e-"));
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_butlervolmerreduced);
      eleparams.set<int>("numscal",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<int>("numscal"));
      eleparams.set<const std::vector<int>*>("stoichiometries",
          &scatra_manifold_coupling->ConditionKinetics()->parameters().Get<std::vector<int>>(
              "stoichiometries"));
      eleparams.set<int>("numelectrons",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<int>("e-"));
      eleparams.set<double>(
          "k_r", scatra_manifold_coupling->ConditionKinetics()->parameters().Get<double>("k_r"));
      eleparams.set<double>("alpha_a",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<double>("alpha_a"));
      eleparams.set<double>("alpha_c",
          scatra_manifold_coupling->ConditionKinetics()->parameters().Get<double>("alpha_c"));
      break;
    }
    case INPAR::S2I::kinetics_nointerfaceflux:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_nointerfaceflux);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown kinetics type for manifold couplign");
      break;
    }
  }
  scatra_manifold_->ScaTraField()->discretization()->Evaluate(
      eleparams, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::Output()
{
  FOUR_C_ASSERT(runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");

  std::map<std::string, std::vector<double>> output_data;
  for (const auto& inflow_comp : inflow_)
  {
    const std::string manifold_string = "manifold " + std::to_string(inflow_comp.first);

    // find pair of domain integrals with same key for current condition and get domain integral
    const auto domainintegral_cond = domainintegral_.find(inflow_comp.first);
    const double domainint = domainintegral_cond->second;

    output_data["Integral of " + manifold_string] = {domainint};

    for (int i = 0; i < static_cast<int>(inflow_comp.second.size()); ++i)
    {
      output_data["Total flux of scalar " + std::to_string(i + 1) + " into " + manifold_string] = {
          inflow_comp.second[i]};
      output_data["Mean flux of scalar " + std::to_string(i + 1) + " into " + manifold_string] = {
          inflow_comp.second[i] / domainint};
    }
  }
  runtime_csvwriter_->WriteDataToFile(scatra_manifold_->ScaTraField()->Time(),
      scatra_manifold_->ScaTraField()->Step(), output_data);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::un_complete_matrices_if_necessary(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  // definition of lambda function to get size of graph for given matrix
  auto graph_size = [](Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix)
  { return matrix->EpetraMatrix()->Graph().NumGlobalEntries(); };

  // get size of graphs of conditions matrices
  const int size_manifold_scatra_graph_ = graph_size(matrix_manifold_scatra_cond_);
  const int size_manifold_structure_graph = graph_size(matrix_manifold_structure_cond_);
  const int size_scatra_manifold_graph = graph_size(matrix_scatra_manifold_cond_);
  const int size_scatra_structure_graph = graph_size(matrix_scatra_structure_cond_);
  const int size_manifold_sysmat_graph = graph_size(systemmatrix_manifold_cond_);
  const int size_scatra_sysmat_graph = graph_size(systemmatrix_scatra_cond_);

  // check if size of any condition matrix was updated and store new size
  bool do_uncomplete = false;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::ManifoldScaTra, size_manifold_scatra_graph_))
    do_uncomplete = true;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::ManifoldStructure, size_manifold_structure_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::ScaTraManifold, size_scatra_manifold_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::ScaTraStructure, size_scatra_structure_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::SysMatManifold, size_manifold_sysmat_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->check_and_set_size_of_matrix_graph(
          BlockMatrixType::SysMatScaTra, size_scatra_sysmat_graph))
    do_uncomplete = true;

  // uncomplete all global matrices if condition matrices have updated graph
  if (do_uncomplete)
  {
    matrix_manifold_scatra_->UnComplete();
    matrix_manifold_structure_->UnComplete();
    matrix_scatra_manifold_->UnComplete();
    matrix_scatra_structure_->UnComplete();
    systemmatrix_manifold_->UnComplete();
    systemmatrix_scatra_->UnComplete();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ManifoldMeshTyingStrategyBase::ManifoldMeshTyingStrategyBase(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps, const bool is_manifold_meshtying)
    : is_manifold_meshtying_(is_manifold_meshtying),
      condensed_dof_map_(Teuchos::null),
      ssi_maps_(std::move(ssi_maps)),
      ssi_meshtying_(Teuchos::null)
{
  if (is_manifold_meshtying_)
  {
    ssi_meshtying_ = Teuchos::rcp(
        new SSI::UTILS::SSIMeshTying("SSISurfaceManifold", scatra_manifold_dis, false, false));

    if (ssi_meshtying_->MeshTyingHandlers().empty())
    {
      FOUR_C_THROW(
          "Could not create mesh tying between manifold fields. They are not intersecting. "
          "Disable 'MESHTYING_MANIFOLD' or create intersecting manifold conditions.");
    }

    // merge slave dof maps from all mesh tying conditions
    Teuchos::RCP<Epetra_Map> slave_dof_map = Teuchos::null;
    for (const auto& meshtying : ssi_meshtying_->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      if (slave_dof_map == Teuchos::null)
        slave_dof_map = Teuchos::rcp(new Epetra_Map(*coupling_adapter->SlaveDofMap()));
      else
      {
        auto slave_dof_map_old = Teuchos::rcp(new Epetra_Map(*slave_dof_map));
        slave_dof_map = CORE::LINALG::MergeMap(slave_dof_map_old, coupling_adapter->SlaveDofMap());
      }
    }
    // exclusive interior and master dofs across all slave conditions
    condensed_dof_map_ =
        CORE::LINALG::SplitMap(*ssi_maps_->sca_tra_manifold_dof_row_map(), *slave_dof_map);
  }
  else
    condensed_dof_map_ = ssi_maps_->sca_tra_manifold_dof_row_map();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ManifoldMeshTyingStrategySparse::ManifoldMeshTyingStrategySparse(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis, Teuchos::RCP<UTILS::SSIMaps> ssi_maps,
    const bool is_manifold_meshtying)
    : ManifoldMeshTyingStrategyBase(scatra_manifold_dis, ssi_maps, is_manifold_meshtying)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ManifoldMeshTyingStrategyBlock::ManifoldMeshTyingStrategyBlock(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps, const bool is_manifold_meshtying)
    : ManifoldMeshTyingStrategyBase(scatra_manifold_dis, ssi_maps, is_manifold_meshtying),
      condensed_block_dof_map_(Teuchos::null),
      meshtying_block_handler_()
{
  // split condensed_dof_map into blocks
  std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps_condensed_block_dof_map;
  for (int i = 0; i < ssi_maps_->block_map_sca_tra_manifold()->NumMaps(); ++i)
  {
    partial_maps_condensed_block_dof_map.emplace_back(CORE::LINALG::IntersectMap(
        *condensed_dof_map_, *ssi_maps_->block_map_sca_tra_manifold()->Map(i)));
  }

  condensed_block_dof_map_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor(
      *condensed_dof_map_, partial_maps_condensed_block_dof_map));

  if (is_manifold_meshtying_)
  {
    // couple meshyting_handler_ and condensed_block_dof_map_ to meshtying_block_handler_
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto slave_dof_map = coupling_adapter->SlaveDofMap();
      auto perm_slave_dof_map = coupling_adapter->PermSlaveDofMap();
      auto master_dof_map = coupling_adapter->MasterDofMap();
      auto perm_master_dof_map = coupling_adapter->PermMasterDofMap();

      // split maps according to split of matrix blocks, i.e. block maps
      std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps_slave_block_dof_map;
      std::vector<Teuchos::RCP<CORE::ADAPTER::Coupling>> partial_block_adapters;
      for (int i = 0; i < ssi_maps_->block_map_sca_tra_manifold()->NumMaps(); ++i)
      {
        auto [slave_block_map, perm_master_block_map] =
            intersect_coupling_maps_block_map(ssi_maps_->block_map_sca_tra_manifold()->Map(i),
                slave_dof_map, perm_master_dof_map, scatra_manifold_dis->Comm());

        auto [master_block_map, perm_slave_block_map] =
            intersect_coupling_maps_block_map(ssi_maps_->block_map_sca_tra_manifold()->Map(i),
                master_dof_map, perm_slave_dof_map, scatra_manifold_dis->Comm());

        auto coupling_adapter_block = Teuchos::rcp(new CORE::ADAPTER::Coupling());
        coupling_adapter_block->setup_coupling(
            slave_block_map, perm_slave_block_map, master_block_map, perm_master_block_map);

        partial_maps_slave_block_dof_map.emplace_back(slave_block_map);
        partial_block_adapters.emplace_back(coupling_adapter_block);
      }

      meshtying_block_handler_.emplace_back(
          partial_block_adapters, partial_maps_slave_block_dof_map);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::tuple<Teuchos::RCP<const Epetra_Map>, Teuchos::RCP<const Epetra_Map>>
SSI::ManifoldMeshTyingStrategyBlock::intersect_coupling_maps_block_map(
    Teuchos::RCP<const Epetra_Map> block_map, Teuchos::RCP<const Epetra_Map> intersecting_map,
    Teuchos::RCP<const Epetra_Map> permuted_map, const Epetra_Comm& comm)
{
  std::vector<int> intersecting_map_vec, permuted_intersecting_map_vec;
  for (int slave_lid = 0; slave_lid < intersecting_map->NumMyElements(); ++slave_lid)
  {
    const int slave_gid = intersecting_map->GID(slave_lid);
    if (block_map->LID(slave_gid) != -1)
    {
      intersecting_map_vec.emplace_back(slave_gid);
      permuted_intersecting_map_vec.emplace_back(permuted_map->GID(slave_lid));
    }
  }

  auto intersected_map = Teuchos::rcp(new const Epetra_Map(
      -1, static_cast<int>(intersecting_map_vec.size()), intersecting_map_vec.data(), 0, comm));

  auto permuted_intersected_map =
      Teuchos::rcp(new const Epetra_Map(-1, static_cast<int>(permuted_intersecting_map_vec.size()),
          permuted_intersecting_map_vec.data(), 0, comm));

  return {intersected_map, permuted_intersected_map};
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBase::apply_mesh_tying_to_manifold_rhs(
    Teuchos::RCP<Epetra_Vector> rhs_manifold)
{
  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();
      auto multimap = meshtying->slave_master_extractor();

      auto slave_dofs = multimap->ExtractVector(rhs_manifold, 1);
      auto slave_to_master_dofs = coupling_adapter->SlaveToMaster(slave_dofs);
      multimap->AddVector(slave_to_master_dofs, 2, rhs_manifold);
      multimap->PutScalar(*rhs_manifold, 1, 0.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::apply_meshtying_to_manifold_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_sparse = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_sparse = CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_matrix);

  // add derivs. of interior/master dofs. w.r.t. interior/master dofs
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *condensed_dof_map_,
      *condensed_dof_map_, 1.0, nullptr, nullptr, *ssi_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = CORE::ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs. w.r.t. slave dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *cond_slave_dof_map, 1.0, &converter, &converter, *ssi_manifold_sparse, true, true);
      // add derivs. of slave dofs. w.r.t. interior/master dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *condensed_dof_map_, 1.0, &converter, nullptr, *ssi_manifold_sparse, true, true);
      // add derivs. of interior/master dofs w.r.t. slave dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *condensed_dof_map_,
          *cond_slave_dof_map, 1.0, nullptr, &converter, *ssi_manifold_sparse, true, true);
    }

    // Finalize: put 1.0 on main diag of slave dofs
    const double one = 1.0;
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();

      auto slave_dof_map = coupling_adapter->SlaveDofMap();
      for (int doflid_slave = 0; doflid_slave < slave_dof_map->NumMyElements(); ++doflid_slave)
      {
        // extract global ID of current slave-side row
        const int dofgid_slave = slave_dof_map->GID(doflid_slave);
        if (dofgid_slave < 0) FOUR_C_THROW("Local ID not found!");

        // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
        if (ssi_manifold_sparse->Filled())
        {
          const int rowlid_slave = ssi_manifold_sparse->RowMap().LID(dofgid_slave);
          if (rowlid_slave < 0) FOUR_C_THROW("Global ID not found!");
          if (ssi_manifold_sparse->EpetraMatrix()->ReplaceMyValues(
                  rowlid_slave, 1, &one, &rowlid_slave))
            FOUR_C_THROW("ReplaceMyValues failed!");
        }

        // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column
        // indices
        else if (ssi_manifold_sparse->EpetraMatrix()->InsertGlobalValues(
                     dofgid_slave, 1, &one, &dofgid_slave))
          FOUR_C_THROW("InsertGlobalValues failed!");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::apply_meshtying_to_manifold_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_block =
      CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
      // add derivs. of interior/master dofs. w.r.t. interior/master dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
          *condensed_block_dof_map_->Map(row), *condensed_block_dof_map_->Map(col), 1.0, nullptr,
          nullptr, ssi_manifold_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : mesh_tying_block_handler())
        {
          auto meshtying = block_meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter_row = CORE::ADAPTER::CouplingSlaveConverter(*meshtying[row]);
          auto converter_col = CORE::ADAPTER::CouplingSlaveConverter(*meshtying[col]);

          // add derivs. of slave dofs. w.r.t. slave dofs
          CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map[row], *cond_block_slave_dof_map[col], 1.0, &converter_row,
              &converter_col, ssi_manifold_block->Matrix(row, col), true, true);
          // add derivs. of slave dofs. w.r.t. interior/master dofs
          CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map[row], *condensed_block_dof_map_->Map(col), 1.0,
              &converter_row, nullptr, ssi_manifold_block->Matrix(row, col), true, true);
          // add derivs. of interior/master dofs w.r.t. slave dofs
          CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *condensed_block_dof_map_->Map(row), *cond_block_slave_dof_map[col], 1.0, nullptr,
              &converter_col, ssi_manifold_block->Matrix(row, col), true, true);
        }

        // Finalize: put 1.0 on main diag of slave dofs
        const double one = 1.0;
        if (row == col)
        {
          for (const auto& block_meshtying : mesh_tying_block_handler())
          {
            auto coupling_adapter = block_meshtying.first[row];
            auto slave_dof_map = coupling_adapter->SlaveDofMap();
            for (int doflid_slave = 0; doflid_slave < slave_dof_map->NumMyElements();
                 ++doflid_slave)
            {
              // extract global ID of current slave-side row
              const int dofgid_slave = slave_dof_map->GID(doflid_slave);
              if (dofgid_slave < 0) FOUR_C_THROW("Local ID not found!");

              // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column
              // indices
              if (ssi_manifold_block->Matrix(row, row).Filled())
              {
                const int rowlid_slave =
                    ssi_manifold_block->Matrix(row, row).RowMap().LID(dofgid_slave);
                if (rowlid_slave < 0) FOUR_C_THROW("Global ID not found!");
                if (ssi_manifold_block->Matrix(row, row).EpetraMatrix()->ReplaceMyValues(
                        rowlid_slave, 1, &one, &rowlid_slave))
                  FOUR_C_THROW("ReplaceMyValues failed!");
              }

              // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and
              // column indices
              else if (ssi_manifold_block->Matrix(row, row).EpetraMatrix()->InsertGlobalValues(
                           dofgid_slave, 1, &one, &dofgid_slave))
                FOUR_C_THROW("InsertGlobalValues failed!");
            }
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::apply_meshtying_to_manifold_scatra_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_scatra_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_sparse =
      CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_scatra_matrix);

  // add derivs. of interior/master dofs w.r.t. scatra dofs
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_scatra_sparse, *condensed_dof_map_,
      *ssi_maps_->ScaTraDofRowMap(), 1.0, nullptr, nullptr, *ssi_manifold_scatra_sparse, true,
      true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = CORE::ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs w.r.t. scatra dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_scatra_sparse, *cond_slave_dof_map,
          *ssi_maps_->ScaTraDofRowMap(), 1.0, &converter, nullptr, *ssi_manifold_scatra_sparse,
          true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::apply_meshtying_to_manifold_scatra_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_scatra_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_block =
      CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_scatra_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < ssi_maps_->BlockMapScaTra()->NumMaps(); ++col)
    {
      // add derivs. of interior/master dofs w.r.t. scatra dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_scatra_block->Matrix(row, col),
          *condensed_block_dof_map_->Map(row), *ssi_maps_->BlockMapScaTra()->Map(col), 1.0, nullptr,
          nullptr, ssi_manifold_scatra_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : mesh_tying_block_handler())
        {
          auto meshtying = block_meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter_row = CORE::ADAPTER::CouplingSlaveConverter(*meshtying[row]);

          // add derivs. of slave dofs w.r.t. scatra dofs
          CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_scatra_block->Matrix(row, col),
              *cond_block_slave_dof_map[row], *ssi_maps_->BlockMapScaTra()->Map(col), 1.0,
              &converter_row, nullptr, ssi_manifold_scatra_block->Matrix(row, col), true, true);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::apply_meshtying_to_manifold_structure_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_structure_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_structure_matrix,
    const bool do_uncomplete)
{
  auto ssi_manifold_structure_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_sparse =
      CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure =
      SSI::UTILS::SSIMatrices::setup_sparse_matrix(ssi_maps_->sca_tra_manifold_dof_row_map());

  // add derivs. of interior/master dofs w.r.t. structure dofs
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*ssi_manifold_structure_sparse,
      *condensed_dof_map_, *ssi_maps_->StructureDofRowMap(), 1.0, nullptr, nullptr,
      *temp_manifold_structure, true, true);
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_structure_sparse, *condensed_dof_map_,
      *ssi_maps_->StructureDofRowMap(), 1.0, nullptr, nullptr, *temp_manifold_structure, true,
      true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = CORE::ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs w.r.t. structure dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*ssi_manifold_structure_sparse,
          *cond_slave_dof_map, *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr,
          *temp_manifold_structure, true, true);
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*manifold_structure_sparse,
          *cond_slave_dof_map, *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr,
          *temp_manifold_structure, true, true);
    }
  }

  ssi_manifold_structure_sparse->Zero();
  if (do_uncomplete) ssi_manifold_structure_sparse->UnComplete();
  temp_manifold_structure->Complete(
      *ssi_maps_->StructureDofRowMap(), *ssi_maps_->sca_tra_manifold_dof_row_map());
  ssi_manifold_structure_sparse->Add(*temp_manifold_structure, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::apply_meshtying_to_manifold_structure_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_manifold_structure_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> manifold_structure_matrix,
    const bool do_uncomplete)
{
  auto ssi_manifold_structure_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_block =
      CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure = SSI::UTILS::SSIMatrices::setup_block_matrix(
      ssi_maps_->block_map_sca_tra_manifold(), ssi_maps_->BlockMapStructure());

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    // add derivs. of interior/master dofs w.r.t. structure dofs
    CORE::LINALG::MatrixLogicalSplitAndTransform()(ssi_manifold_structure_block->Matrix(row, 0),
        *condensed_block_dof_map_->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, nullptr,
        nullptr, temp_manifold_structure->Matrix(row, 0), true, true);
    CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_structure_block->Matrix(row, 0),
        *condensed_block_dof_map_->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, nullptr,
        nullptr, temp_manifold_structure->Matrix(row, 0), true, true);

    if (is_manifold_meshtying_)
    {
      for (const auto& block_meshtying : mesh_tying_block_handler())
      {
        auto meshtying = block_meshtying.first;
        auto cond_block_slave_dof_map = block_meshtying.second;
        auto converter_row = CORE::ADAPTER::CouplingSlaveConverter(*meshtying[row]);

        // add derivs. of slave dofs w.r.t. structure dofs
        CORE::LINALG::MatrixLogicalSplitAndTransform()(ssi_manifold_structure_block->Matrix(row, 0),
            *cond_block_slave_dof_map[row], *ssi_maps_->StructureDofRowMap(), 1.0, &converter_row,
            nullptr, temp_manifold_structure->Matrix(row, 0), true, true);
        CORE::LINALG::MatrixLogicalSplitAndTransform()(manifold_structure_block->Matrix(row, 0),
            *cond_block_slave_dof_map[row], *ssi_maps_->StructureDofRowMap(), 1.0, &converter_row,
            nullptr, temp_manifold_structure->Matrix(row, 0), true, true);
      }
    }
  }

  ssi_manifold_structure_block->Zero();
  if (do_uncomplete) ssi_manifold_structure_block->UnComplete();
  temp_manifold_structure->Complete();
  ssi_manifold_structure_block->Add(*temp_manifold_structure, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::apply_meshtying_to_scatra_manifold_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_scatra_manifold_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> scatra_manifold_matrix,
    const bool do_uncomplete)
{
  auto ssi_scatra_manifold_sparse =
      CORE::LINALG::CastToSparseMatrixAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_sparse =
      CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_manifold_matrix);

  if (do_uncomplete) ssi_scatra_manifold_sparse->UnComplete();

  // add derivs. of scatra w.r.t. interior/master dofs
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*scatra_manifold_sparse,
      *ssi_maps_->ScaTraDofRowMap(), *condensed_dof_map_, 1.0, nullptr, nullptr,
      *ssi_scatra_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : SSIMeshTying()->MeshTyingHandlers())
    {
      auto coupling_adapter = meshtying->SlaveMasterCoupling();

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = CORE::ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of scatra w.r.t. slave dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(*scatra_manifold_sparse,
          *ssi_maps_->ScaTraDofRowMap(), *cond_slave_dof_map, 1.0, nullptr, &converter,
          *ssi_scatra_manifold_sparse, true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::apply_meshtying_to_scatra_manifold_matrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> ssi_scatra_manifold_matrix,
    Teuchos::RCP<const CORE::LINALG::SparseOperator> scatra_manifold_matrix,
    const bool do_uncomplete)
{
  auto ssi_scatra_manifold_block =
      CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_block =
      CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_manifold_matrix);

  if (do_uncomplete) ssi_scatra_manifold_block->UnComplete();

  for (int row = 0; row < ssi_maps_->BlockMapScaTra()->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
      // add derivs. of scatra w.r.t. interior/master dofs
      CORE::LINALG::MatrixLogicalSplitAndTransform()(scatra_manifold_block->Matrix(row, col),
          *ssi_maps_->BlockMapScaTra()->Map(row), *condensed_block_dof_map_->Map(col), 1.0, nullptr,
          nullptr, ssi_scatra_manifold_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : mesh_tying_block_handler())
        {
          auto meshtying = block_meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter_col = CORE::ADAPTER::CouplingSlaveConverter(*meshtying[col]);

          // add derivs. of scatra w.r.t. slave dofs
          CORE::LINALG::MatrixLogicalSplitAndTransform()(scatra_manifold_block->Matrix(row, col),
              *ssi_maps_->BlockMapScaTra()->Map(row), *cond_block_slave_dof_map[col], 1.0, nullptr,
              &converter_col, ssi_scatra_manifold_block->Matrix(row, col), true, true);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> SSI::BuildManifoldMeshTyingStrategy(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis, Teuchos::RCP<UTILS::SSIMaps> ssi_maps,
    const bool is_manifold_meshtying, CORE::LINALG::MatrixType matrixtype_manifold)
{
  Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> meshtyingstrategy = Teuchos::null;

  switch (matrixtype_manifold)
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      meshtyingstrategy = Teuchos::rcp(new SSI::ManifoldMeshTyingStrategyBlock(
          scatra_manifold_dis, ssi_maps, is_manifold_meshtying));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      meshtyingstrategy = Teuchos::rcp(new SSI::ManifoldMeshTyingStrategySparse(
          scatra_manifold_dis, ssi_maps, is_manifold_meshtying));
      break;
    }

    default:
    {
      FOUR_C_THROW("unknown matrix type of Manifold field");
      break;
    }
  }

  return meshtyingstrategy;
}

FOUR_C_NAMESPACE_CLOSE
