/*----------------------------------------------------------------------*/
/*! \file
 \brief Evaluates flux between ScaTra and ScaTra on manifolds incl. coupling matrices

 \level 2


 *------------------------------------------------------------------------------------------------*/

#include "ssi_manifold_utils.H"

#include "ssi_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_s2i.H"
#include "../drt_inpar/inpar_ssi.H"

#include "../drt_io/runtime_csv_writer.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_matchingoctree.H"
#include "../drt_lib/drt_utils_gid_vector.H"
#include "../drt_lib/drt_utils_parameter_list.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::ManifoldScaTraCoupling::ManifoldScaTraCoupling(Teuchos::RCP<DRT::Discretization> manifolddis,
    Teuchos::RCP<DRT::Discretization> scatradis, DRT::Condition* condition_manifold,
    DRT::Condition* condition_kinetics, const int ndof_per_node)
    : condition_kinetics_(condition_kinetics),
      condition_manifold_(condition_manifold),
      coupling_adapter_(Teuchos::rcp(new ADAPTER::Coupling())),
      manifold_conditionID_(condition_manifold->GetInt("ConditionID")),
      kinetics_conditionID_(condition_kinetics->GetInt("ConditionID")),
      manifold_map_extractor_(Teuchos::null),
      scatra_map_extractor_(Teuchos::null),
      size_matrix_graph_(),
      master_converter_(Teuchos::null)
{
  std::vector<int> inodegidvec_manifold;
  DRT::UTILS::AddOwnedNodeGIDVector(
      *manifolddis, *condition_manifold->Nodes(), inodegidvec_manifold);

  std::vector<int> inodegidvec_scatra;
  DRT::UTILS::AddOwnedNodeGIDVector(*scatradis, *condition_kinetics->Nodes(), inodegidvec_scatra);

  coupling_adapter_->SetupCoupling(*scatradis, *manifolddis, inodegidvec_scatra,
      inodegidvec_manifold, ndof_per_node, true, 1.0e-8);
  master_converter_ = Teuchos::rcp(new ADAPTER::CouplingMasterConverter(*coupling_adapter_));

  scatra_map_extractor_ = Teuchos::rcp(
      new LINALG::MapExtractor(*scatradis->DofRowMap(), coupling_adapter_->MasterDofMap(), true));

  manifold_map_extractor_ = Teuchos::rcp(
      new LINALG::MapExtractor(*manifolddis->DofRowMap(), coupling_adapter_->SlaveDofMap(), true));

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
bool SSI::ManifoldScaTraCoupling::CheckAndSetSizeOfMatrixGraph(
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
    const SSI::SSIMono& ssi_mono)
    : block_map_scatra_(ssi_mono.BlockMapScaTra()),
      block_map_scatra_manifold_(ssi_mono.BlockMapScaTraManifold()),
      block_map_structure_(ssi_mono.BlockMapStructure()),
      do_output_(DRT::INPUT::IntegralValue<bool>(
          DRT::Problem::Instance()->SSIControlParams().sublist("MANIFOLD"), "OUTPUT_INFLOW")),
      full_map_manifold_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold))),
      full_map_scatra_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport))),
      full_map_structure_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::structure))),
      matrix_manifold_scatra_(Teuchos::null),
      matrix_manifold_structure_(Teuchos::null),
      matrix_scatra_manifold_(Teuchos::null),
      matrix_scatra_structure_(Teuchos::null),
      rhs_manifold_(Teuchos::null),
      rhs_scatra_(Teuchos::null),
      runtime_csvwriter_(nullptr),
      scatra_(ssi_mono.ScaTraBaseAlgorithm()),
      scatra_manifold_(ssi_mono.ScaTraManifoldBaseAlgorithm()),
      scatra_manifold_couplings_(Teuchos::null),
      ssi_structure_meshtying_(ssi_mono.SSIStructureMeshTying()),
      systemmatrix_manifold_(Teuchos::null),
      systemmatrix_scatra_(Teuchos::null)
{
  // safety check before setup of coupling
  if (ssi_mono.ScaTraField()->NumDofPerNode() != ssi_mono.ScaTraManifold()->NumDofPerNode())
    dserror("Number of dofs per node of scatra field and scatra manifold field must be equal");

  std::vector<DRT::Condition*> conditions_manifold;
  scatra_manifold_->ScaTraField()->Discretization()->GetCondition(
      "SSISurfaceManifold", conditions_manifold);

  std::vector<DRT::Condition*> conditions_manifold_kinetics_scatra;
  scatra_->ScaTraField()->Discretization()->GetCondition(
      "SSISurfaceManifoldKinetics", conditions_manifold_kinetics_scatra);

  // create pair: manifold condition - kinetics condition
  for (const auto& condition_manifold : conditions_manifold)
  {
    for (const auto& condition_kinetics : conditions_manifold_kinetics_scatra)
    {
      if (condition_manifold->GetInt("ConditionID") ==
          condition_kinetics->GetInt("ManifoldConditionID"))
      {
        scatra_manifold_couplings_.emplace_back(Teuchos::rcp(
            new SSI::ManifoldScaTraCoupling(scatra_manifold_->ScaTraField()->Discretization(),
                scatra_->ScaTraField()->Discretization(), condition_manifold, condition_kinetics,
                ssi_mono.ScaTraManifold()->NumDofPerNode())));
      }
    }
  }

  rhs_manifold_ = LINALG::CreateVector(*full_map_manifold_, true);
  rhs_scatra_ = LINALG::CreateVector(*full_map_scatra_, true);

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      systemmatrix_manifold_ = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
          block_map_scatra_manifold_, block_map_scatra_manifold_);
      systemmatrix_scatra_ =
          SSI::UTILS::SSIMatrices::SetupBlockMatrix(block_map_scatra_, block_map_scatra_);
      matrix_manifold_structure_ = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
          block_map_scatra_manifold_, block_map_structure_);
      matrix_manifold_scatra_ =
          SSI::UTILS::SSIMatrices::SetupBlockMatrix(block_map_scatra_manifold_, block_map_scatra_);
      matrix_scatra_manifold_ =
          SSI::UTILS::SSIMatrices::SetupBlockMatrix(block_map_scatra_, block_map_scatra_manifold_);
      matrix_scatra_structure_ =
          SSI::UTILS::SSIMatrices::SetupBlockMatrix(block_map_scatra_, block_map_structure_);

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      systemmatrix_manifold_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
      systemmatrix_scatra_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);
      matrix_manifold_structure_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
      matrix_manifold_scatra_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
      matrix_scatra_manifold_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);
      matrix_scatra_structure_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  rhs_manifold_cond_ = LINALG::CreateVector(*full_map_manifold_, true);
  rhs_scatra_cond_ = LINALG::CreateVector(*full_map_scatra_, true);

  systemmatrix_manifold_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
  systemmatrix_scatra_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);
  matrix_manifold_scatra_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
  matrix_manifold_structure_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_manifold_);
  matrix_scatra_manifold_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);
  matrix_scatra_structure_cond_ = SSI::UTILS::SSIMatrices::SetupSparseMatrix(full_map_scatra_);

  // Prepare runtime csv writer
  if (DoOutput())
  {
    runtime_csvwriter_ = std::make_shared<RuntimeCsvWriter>(ssi_mono.Comm().MyPID());

    runtime_csvwriter_->Init("manifold_inflow");

    for (const auto& condition_manifold : conditions_manifold)
    {
      const std::string manifold_string =
          "manifold " + std::to_string(condition_manifold->GetInt("ConditionID"));

      runtime_csvwriter_->RegisterDataVector("Integral of " + manifold_string, 1, 16);

      for (int k = 0; k < ssi_mono.ScaTraManifold()->NumDofPerNode(); ++k)
      {
        runtime_csvwriter_->RegisterDataVector(
            "Total flux of scalar " + std::to_string(k + 1) + " into " + manifold_string, 1, 16);

        runtime_csvwriter_->RegisterDataVector(
            "Mean flux of scalar " + std::to_string(k + 1) + " into " + manifold_string, 1, 16);
      }
    }

    runtime_csvwriter_->Setup();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteMatrixManifoldScaTra()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      matrix_manifold_scatra_->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      matrix_manifold_scatra_->Complete(*full_map_scatra_, *full_map_manifold_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteMatrixManifoldStructure()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      matrix_manifold_structure_->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      matrix_manifold_structure_->Complete(*full_map_structure_, *full_map_manifold_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteMatrixScaTraManifold()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      matrix_scatra_manifold_->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      matrix_scatra_manifold_->Complete(*full_map_manifold_, *full_map_scatra_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteMatrixScaTraStructure()
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      matrix_scatra_structure_->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      matrix_scatra_structure_->Complete(*full_map_structure_, *full_map_scatra_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteSystemMatrixManifold()
{
  systemmatrix_manifold_->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteSystemMatrixScaTra()
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
  for (auto scatra_manifold_coupling : scatra_manifold_couplings_)
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

    EvaluateBulkSide(scatra_manifold_coupling);

    CopyScaTraScaTraManifoldSide(scatra_manifold_coupling);

    // This is needed because the graph of the matrices could change from step to step in case we
    // have zero flux (and then zero entries in the matrices)
    UnCompleteMatricesIfNecessary(scatra_manifold_coupling);

    AddConditionContribution();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateBulkSide(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  // First: Set parameters to elements
  PreEvaluate(scatra_manifold_coupling);

  scatra_->ScaTraField()->Discretization()->ClearState();
  scatra_->ScaTraField()->Discretization()->SetState("phinp", scatra_->ScaTraField()->Phinp());

  // Second: Evaluate condition
  {
    // manifold-scatra coupling matrix evaluated on scatra side
    auto matrix_scatra_manifold_cond_on_scatra_side =
        Teuchos::rcp(new LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));

    Teuchos::ParameterList condparams;
    DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
        "action", SCATRA::BoundaryAction::calc_s2icoupling, condparams);
    condparams.set<int>("evaluate_manifold_coupling", 1);

    scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

    // Evaluation of RHS and scatra-manifold coupling matrices
    {
      DRT::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
          "differentiationtype", SCATRA::DifferentiationType::elch, condparams);

      // dscatra_dscatra, dscatra_dmanifold (on scatra side)
      DRT::AssembleStrategy strategyscatra(0, 0, systemmatrix_scatra_cond_,
          matrix_scatra_manifold_cond_on_scatra_side, rhs_scatra_cond_, Teuchos::null,
          Teuchos::null);

      scatra_->ScaTraField()->Discretization()->EvaluateCondition(condparams, strategyscatra,
          "SSISurfaceManifoldKinetics", scatra_manifold_coupling->KineticsConditionID());

      systemmatrix_scatra_cond_->Complete();
      matrix_scatra_manifold_cond_on_scatra_side->Complete();

      // dscatra_dmanifold (on scatra side) -> dscatra_dmanifold
      LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_manifold_cond_on_scatra_side,
          *full_map_scatra_, *full_map_scatra_, 1.0, nullptr,
          &*scatra_manifold_coupling->MasterConverter(), *matrix_scatra_manifold_cond_, true, true);
      matrix_scatra_manifold_cond_->Complete(*full_map_manifold_, *full_map_scatra_);
    }

    // Evaluation of linearization w.r.t. displacement
    {
      DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
          "action", SCATRA::BoundaryAction::calc_s2icoupling_od, condparams);

      DRT::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
          "differentiationtype", SCATRA::DifferentiationType::disp, condparams);

      // dscatra_dstructure
      auto matrix_scatra_structure_cond_slave_side_disp_evaluate =
          Teuchos::rcp(new LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));

      DRT::AssembleStrategy strategyscatra(0, 1,
          matrix_scatra_structure_cond_slave_side_disp_evaluate, Teuchos::null, Teuchos::null,
          Teuchos::null, Teuchos::null);

      scatra_->ScaTraField()->Discretization()->EvaluateCondition(condparams, strategyscatra,
          "SSISurfaceManifoldKinetics", scatra_manifold_coupling->KineticsConditionID());

      scatra_manifold_->ScaTraField()->Discretization()->ClearState();

      matrix_scatra_structure_cond_slave_side_disp_evaluate->Complete(
          *full_map_structure_, *full_map_scatra_);

      // "slave side" from manifold and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      auto matrix_scatra_structure_cond_slave_side_disp =
          Teuchos::rcp(new LINALG::SparseMatrix(*full_map_scatra_, 27, false, true));
      for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter = ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        LINALG::MatrixLogicalSplitAndTransform()(
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
      for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
      {
        auto cond_slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto converter = meshtying->SlaveSideConverter();

        // assemble derivatives of x w.r.t. structure slave dofs
        LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_structure_cond_slave_side_disp,
            matrix_scatra_structure_cond_slave_side_disp->RangeMap(), *cond_slave_dof_map, 1.0,
            nullptr, &(*converter), *matrix_scatra_structure_cond_, true, true);
      }
      matrix_scatra_structure_cond_->Complete(*full_map_structure_, *full_map_scatra_);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CopyScaTraScaTraManifoldSide(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  {
    auto rhs_scatra_cond_extract =
        scatra_manifold_coupling->ScaTraMapExtractor()->ExtractCondVector(rhs_scatra_cond_);

    auto rhs_manifold_cond_extract =
        scatra_manifold_coupling->CouplingAdapter()->MasterToSlave(rhs_scatra_cond_extract);

    scatra_manifold_coupling->ManifoldMapExtractor()->AddCondVector(
        rhs_manifold_cond_extract, rhs_manifold_cond_);
    rhs_manifold_cond_->Scale(-1.0);
  }

  // dmanifold_dscatra: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_scatra_cond_, *full_map_scatra_,
      *full_map_scatra_, -1.0, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *matrix_manifold_scatra_cond_, true, true);

  matrix_manifold_scatra_cond_->Complete(*full_map_scatra_, *full_map_manifold_);

  // dmanifold_dmanifold: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_manifold_cond_, *full_map_scatra_,
      *full_map_manifold_, -1.0, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *systemmatrix_manifold_cond_, true, true);

  systemmatrix_manifold_cond_->Complete();

  // dmanifold_dstructure: scatra rows are transformed to manifold side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*matrix_scatra_structure_cond_, *full_map_scatra_,
      *full_map_structure_, -1.0, &*scatra_manifold_coupling->MasterConverter(), nullptr,
      *matrix_manifold_structure_cond_, true, true);

  matrix_manifold_structure_cond_->Complete(*full_map_structure_, *full_map_manifold_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::AddConditionContribution()
{
  rhs_manifold_->Update(1.0, *rhs_manifold_cond_, 1.0);
  rhs_scatra_->Update(1.0, *rhs_scatra_cond_, 1.0);

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      auto blockmaps_manifold = scatra_manifold_->ScaTraField()->BlockMaps();

      /*
       * m: manifold
       * s: scatra
       * d: structure/displacements
       */
      auto flux_manifold_scatra_mm_block =
          systemmatrix_manifold_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
              blockmaps_manifold, blockmaps_manifold);
      auto flux_manifold_scatra_md_block =
          matrix_manifold_structure_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, blockmaps_manifold);
      auto flux_manifold_scatra_ms_block =
          matrix_manifold_scatra_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
              *block_map_scatra_, blockmaps_manifold);
      auto flux_manifold_scatra_sm_block =
          matrix_scatra_manifold_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
              blockmaps_manifold, *block_map_scatra_);
      auto flux_manifold_scatra_sd_block =
          matrix_scatra_structure_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, *block_map_scatra_);
      auto flux_manifold_scatra_ss_block =
          systemmatrix_scatra_cond_->Split<LINALG::DefaultBlockMatrixStrategy>(
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
    case LINALG::MatrixType::sparse:
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
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateScaTraManifoldInflow()
{
  inflow_.clear();
  domainintegral_.clear();

  scatra_->ScaTraField()->Discretization()->SetState("phinp", scatra_->ScaTraField()->Phinp());

  for (const auto& scatra_manifold_coupling : scatra_manifold_couplings_)
  {
    const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

    std::vector<double> zero_scalar_vector(scatra_->ScaTraField()->NumDofPerNode(), 0.0);
    inflow_.insert(std::make_pair(kineticsID, zero_scalar_vector));

    // First: set parameters to elements
    PreEvaluate(scatra_manifold_coupling);

    // Second: evaluate condition
    EvaluateScaTraManifoldInflowIntegral(scatra_manifold_coupling);

    // Third: evaluate domain integral
    EvaluateScaTraManifoldDomainIntegral(scatra_manifold_coupling);
  }
  scatra_->ScaTraField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateScaTraManifoldDomainIntegral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

  // integrate only if not done so far
  if (domainintegral_.find(kineticsID) == domainintegral_.end())
  {
    Teuchos::ParameterList condparams;

    DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
        "action", SCATRA::BoundaryAction::calc_boundary_integral, condparams);

    // integrated domain of this condition
    auto domainintegral_cond = Teuchos::rcp(new Epetra_SerialDenseVector(1));

    scatra_->ScaTraField()->Discretization()->EvaluateScalars(
        condparams, domainintegral_cond, "SSISurfaceManifold", kineticsID);

    domainintegral_.insert(std::make_pair(kineticsID, domainintegral_cond->Values()[0]));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateScaTraManifoldInflowIntegral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int kineticsID = scatra_manifold_coupling->KineticsConditionID();

  Teuchos::ParameterList condparams;

  DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_flux, condparams);

  condparams.set<int>("evaluate_manifold_coupling", 1);

  // integrated scalars of this condition
  auto inflow_cond =
      Teuchos::rcp(new Epetra_SerialDenseVector(scatra_->ScaTraField()->NumDofPerNode()));

  scatra_->ScaTraField()->Discretization()->EvaluateScalars(
      condparams, inflow_cond, "SSISurfaceManifold", kineticsID);

  inflow_cond->Print(std::cout << std::scientific << std::setprecision(16));

  for (int i = 0; i < inflow_cond->Length(); ++i)
    inflow_.at(kineticsID).at(i) += inflow_cond->Values()[i];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::PreEvaluate(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  Teuchos::ParameterList eleparams;

  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::set_scatra_ele_boundary_parameter, eleparams);

  eleparams.set<DRT::Condition::ConditionType>(
      "condition type", DRT::Condition::ConditionType::S2IKinetics);

  switch (scatra_manifold_coupling->ConditionKinetics()->GetInt("kinetic model"))
  {
    case INPAR::S2I::kinetics_constantinterfaceresistance:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_constantinterfaceresistance);
      eleparams.set<double>(
          "resistance", scatra_manifold_coupling->ConditionKinetics()->GetDouble("resistance"));
      eleparams.set<std::vector<int>*>("onoff",
          scatra_manifold_coupling->ConditionKinetics()->GetMutable<std::vector<int>>("onoff"));
      eleparams.set<int>(
          "numelectrons", scatra_manifold_coupling->ConditionKinetics()->GetInt("e-"));
      break;
    }
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_butlervolmerreduced);
      eleparams.set<int>(
          "numscal", scatra_manifold_coupling->ConditionKinetics()->GetInt("numscal"));
      eleparams.set<std::vector<int>*>("stoichiometries",
          scatra_manifold_coupling->ConditionKinetics()->GetMutable<std::vector<int>>(
              "stoichiometries"));
      eleparams.set<int>(
          "numelectrons", scatra_manifold_coupling->ConditionKinetics()->GetInt("e-"));
      eleparams.set<double>("k_r", scatra_manifold_coupling->ConditionKinetics()->GetDouble("k_r"));
      eleparams.set<double>(
          "alpha_a", scatra_manifold_coupling->ConditionKinetics()->GetDouble("alpha_a"));
      eleparams.set<double>(
          "alpha_c", scatra_manifold_coupling->ConditionKinetics()->GetDouble("alpha_c"));
      break;
    }
    case INPAR::S2I::kinetics_nointerfaceflux:
    {
      eleparams.set<int>("kinetic model", INPAR::S2I::kinetics_nointerfaceflux);
      break;
    }
    default:
    {
      dserror("Unknown kinetics type for manifold couplign");
      break;
    }
  }
  scatra_manifold_->ScaTraField()->Discretization()->Evaluate(
      eleparams, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::Output()
{
  for (const auto& inflow_comp : inflow_)
  {
    const std::string manifold_string = "manifold " + std::to_string(inflow_comp.first);

    // find pair of domain integrals with same key for current condition and get domain integral
    const auto domainintegral_cond = domainintegral_.find(inflow_comp.first);
    const double domainint = domainintegral_cond->second;

    runtime_csvwriter_->AppendDataVector("Integral of " + manifold_string, {domainint});

    for (int i = 0; i < static_cast<int>(inflow_comp.second.size()); ++i)
    {
      runtime_csvwriter_->AppendDataVector(
          "Total flux of scalar " + std::to_string(i + 1) + " into " + manifold_string,
          {inflow_comp.second[i]});
      runtime_csvwriter_->AppendDataVector(
          "Mean flux of scalar " + std::to_string(i + 1) + " into " + manifold_string,
          {inflow_comp.second[i] / domainint});
    }
  }

  runtime_csvwriter_->ResetTimeAndTimeStep(
      scatra_manifold_->ScaTraField()->Time(), scatra_manifold_->ScaTraField()->Step());

  runtime_csvwriter_->WriteFile();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::UnCompleteMatricesIfNecessary(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  // definition of lambda function to get size of graph for given matrix
  auto graph_size = [](Teuchos::RCP<LINALG::SparseMatrix> matrix)
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
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
          BlockMatrixType::ManifoldScaTra, size_manifold_scatra_graph_))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
          BlockMatrixType::ManifoldStructure, size_manifold_structure_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
          BlockMatrixType::ScaTraManifold, size_scatra_manifold_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
          BlockMatrixType::ScaTraStructure, size_scatra_structure_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
          BlockMatrixType::SysMatManifold, size_manifold_sysmat_graph))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeOfMatrixGraph(
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
      meshtying_handler_()
{
  if (is_manifold_meshtying_)
  {
    SetupMeshTyingHandler(scatra_manifold_dis, ssi_maps_);
    Teuchos::RCP<Epetra_Map> slave_dof_map = Teuchos::null;

    if (meshtying_handler_.empty())
    {
      dserror(
          "Could not create mesh tying between manifold fields. They are not intersecting. "
          "Disable 'MESHTYING_MANIFOLD' or create intersecting manifold conditions.");
    }

    // merge slave dof maps from all mesh tying conditions
    for (const auto& meshtying : meshtying_handler_)
    {
      auto coupling_adapter = meshtying.first;
      if (slave_dof_map == Teuchos::null)
        slave_dof_map = Teuchos::rcp(new Epetra_Map(*coupling_adapter->SlaveDofMap()));
      else
      {
        auto slave_dof_map_old = Teuchos::rcp(new Epetra_Map(*slave_dof_map));
        slave_dof_map = LINALG::MergeMap(slave_dof_map_old, coupling_adapter->SlaveDofMap());
      }
    }
    // exclusive interior and master dofs across all slave conditions
    condensed_dof_map_ = LINALG::SplitMap(*ssi_maps_->ScaTraManifoldDofRowMap(), *slave_dof_map);
  }
  else
    condensed_dof_map_ = ssi_maps_->ScaTraManifoldDofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBase::SetupMeshTyingHandler(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps)
{
  // find all couplings between nodes based on SSISurfaceManifold condition
  auto coupling_pairs = ConstructCouplingPairs(scatra_manifold_dis);

  // construct set of unique master GIDs within coupling_vec -> master_gids
  std::set<int> master_gids = FindMasterNodeGIDS(coupling_pairs);

  // assigning slave GIDs (unique key) to master GIDs (non unique values, a master node can have
  // multiple slave nodes) -> coupling_pair
  std::map<int, int> master_slave_pair = DefineMasterSlavePairing(coupling_pairs, master_gids);

  // get number of slave nodes per master node -> max. number gives number of needed adapters
  int my_max_adapters = 0;
  std::map<int, int> assigned_slave_to_master_nodes;
  for (auto pair : master_slave_pair)
  {
    const int master_node_gid = pair.second;
    if (assigned_slave_to_master_nodes.empty())
    {
      assigned_slave_to_master_nodes.insert(std::make_pair(master_node_gid, 1));
      my_max_adapters = 1;
    }
    else
    {
      if (assigned_slave_to_master_nodes.find(master_node_gid) !=
          assigned_slave_to_master_nodes.end())
      {
        assigned_slave_to_master_nodes[master_node_gid]++;
        if (my_max_adapters < assigned_slave_to_master_nodes[master_node_gid])
          my_max_adapters = assigned_slave_to_master_nodes[master_node_gid];
      }
      else
        assigned_slave_to_master_nodes.insert(std::make_pair(master_node_gid, 1));
    }
  }

  int glob_max_adapters = 0;
  scatra_manifold_dis->Comm().MaxAll(&my_max_adapters, &glob_max_adapters, 1);

  // setup coupling adapters
  for (int iadapter = 0; iadapter < glob_max_adapters; ++iadapter)
  {
    std::vector<int> inodegidvec_master;
    std::vector<int> inodegidvec_slave;

    if (!assigned_slave_to_master_nodes.empty())
    {
      for (auto master_gid : master_gids)
      {
        // check if this master node has iadapter + 1 slave nodes
        if (assigned_slave_to_master_nodes.at(master_gid) <= iadapter + 1)
        {
          DRT::UTILS::AddOwnedNodeGID(*scatra_manifold_dis, master_gid, inodegidvec_master);

          int counter = 0;
          for (auto pair : master_slave_pair)
          {
            const int master_gid_coupling = pair.second;

            if (master_gid_coupling == master_gid and counter == iadapter)
            {
              const int slave_gid_coupling = pair.first;

              DRT::UTILS::AddOwnedNodeGID(
                  *scatra_manifold_dis, slave_gid_coupling, inodegidvec_slave);

              counter++;
            }
          }
        }
      }
    }

    // setup coupling adapter
    auto coupling_adapter = Teuchos::rcp(new ADAPTER::Coupling());
    coupling_adapter->SetupCoupling(*scatra_manifold_dis, *scatra_manifold_dis, inodegidvec_master,
        inodegidvec_slave, DRT::Problem::Instance()->NDim() - 1, true, 1.0e-8);

    // setup multimap extractor for each coupling adapter
    auto slave_map = coupling_adapter->SlaveDofMap();
    auto master_map = coupling_adapter->MasterDofMap();
    auto interior_map = LINALG::SplitMap(
        *ssi_maps->ScaTraManifoldDofRowMap(), *LINALG::MergeMap(slave_map, master_map));

    std::vector<Teuchos::RCP<const Epetra_Map>> maps(0, Teuchos::null);
    maps.emplace_back(interior_map);
    maps.emplace_back(master_map);
    maps.emplace_back(slave_map);

    auto coupling_map_extractor =
        Teuchos::rcp(new LINALG::MultiMapExtractor(*ssi_maps->ScaTraManifoldDofRowMap(), maps));
    coupling_map_extractor->CheckForValidMapExtractor();

    // combine coupling adapter and multimap extractor
    meshtying_handler_.emplace_back(std::make_pair(coupling_adapter, coupling_map_extractor));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<std::pair<int, int>> SSI::ManifoldMeshTyingStrategyBase::ConstructCouplingPairs(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis)
{
  // data types to handle coupling between slave and master in coupling_vec
  // coupled_gid_nodes_type: key: master gid, value: slave gid and distance between both
  // condition_pair_type: key: condition a, value: condition b
  using coupled_gid_nodes_type = std::map<int, std::pair<int, double>>;
  using condition_pair_type = std::pair<int, int>;
  std::map<condition_pair_type, coupled_gid_nodes_type> condition_wise_coupling_pairs;

  std::vector<DRT::Condition*> manifold_conditions(0, nullptr);
  scatra_manifold_dis->GetCondition("SSISurfaceManifold", manifold_conditions);

  // fill coupling_vec with slave master pairing using NodeMatchingOctree to find matching nodes
  // between all SSISurfaceManifold conditions (nodes from all procs)
  for (int a = 0; a < static_cast<int>(manifold_conditions.size()); ++a)
  {
    auto* manifold_condition_a = manifold_conditions.at(a);

    // nodes of manifold_condition_a owned by this proc
    std::vector<int> inodegidvec_a;
    DRT::UTILS::AddOwnedNodeGIDVector(
        *scatra_manifold_dis, *manifold_condition_a->Nodes(), inodegidvec_a);

    DRT::UTILS::NodeMatchingOctree tree = DRT::UTILS::NodeMatchingOctree();
    tree.Init(*scatra_manifold_dis, inodegidvec_a, 150, 1.0e-8);
    tree.Setup();

    for (int b = a + 1; b < static_cast<int>(manifold_conditions.size()); ++b)
    {
      auto* manifold_condition_b = manifold_conditions.at(b);

      // nodes of manifold_condition_b owned by this proc
      std::vector<int> inodegidvec_b;
      DRT::UTILS::AddOwnedNodeGIDVector(
          *scatra_manifold_dis, *manifold_condition_b->Nodes(), inodegidvec_b);

      coupled_gid_nodes_type coupled_gid_nodes;
      condition_pair_type condition_pair = std::make_pair(a, b);
      condition_wise_coupling_pairs.emplace(condition_pair, coupled_gid_nodes);

      tree.FindMatch(
          *scatra_manifold_dis, inodegidvec_b, condition_wise_coupling_pairs.at(condition_pair));
    }
  }

  // coupled nodes on from all conditions on all procs
  // split map into vectors to be able to communicate and unite afterwards again
  std::vector<std::pair<int, int>> coupling_pairs;

  std::vector<int> my_gid_vec1;
  std::vector<int> my_gid_vec2;

  // loop over all condition pairs
  for (const auto& coupling : condition_wise_coupling_pairs)
  {
    // loop over all nodal couplings
    for (auto pair : coupling.second)
    {
      my_gid_vec1.emplace_back(pair.first);
      my_gid_vec2.emplace_back(pair.second.first);
    }
  }

  if (my_gid_vec1.size() != my_gid_vec2.size())
    dserror("Size of node GID vectors to be coupled does not match.");

  auto const& comm = scatra_manifold_dis->Comm();
  for (int iproc = 0; iproc < comm.NumProc(); ++iproc)
  {
    // size of vectors of proc iproc
    int size_1 = static_cast<int>(my_gid_vec1.size());
    int size_2 = static_cast<int>(my_gid_vec2.size());

    comm.Broadcast(&size_1, 1, iproc);
    comm.Broadcast(&size_2, 1, iproc);

    // new vectors to be filled (by this proc, if MyPID == iproc or other procs by communication)
    std::vector<int> vec_1, vec_2;
    if (iproc == comm.MyPID())
    {
      vec_1 = my_gid_vec1;
      vec_2 = my_gid_vec2;
    }
    vec_1.resize(size_1);
    vec_2.resize(size_2);
    comm.Broadcast(&vec_1[0], size_1, iproc);
    comm.Broadcast(&vec_2[0], size_2, iproc);

    // reassemble to coupling map on this proc
    for (int i = 0; i < static_cast<int>(size_1); ++i)
    {
      // do not add duplicates
      if (vec_1[i] != vec_2[i]) coupling_pairs.emplace_back(std::make_pair(vec_1[i], vec_2[i]));
    }
  }

  return coupling_pairs;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::set<int> SSI::ManifoldMeshTyingStrategyBase::FindMasterNodeGIDS(
    std::vector<std::pair<int, int>> coupling_pairs)
{
  std::set<int> master_gids;

  for (const auto& nodal_coupling_outer_loop : coupling_pairs)
  {
    const int gid1 = nodal_coupling_outer_loop.first;
    const int gid2 = nodal_coupling_outer_loop.second;

    if (master_gids.empty())
      master_gids.insert(nodal_coupling_outer_loop.first);
    else
    {
      if (master_gids.find(gid1) != master_gids.end() and
          master_gids.find(gid2) != master_gids.end())
      {
        // match in both -> erase one random entry and change coupling
        master_gids.erase(master_gids.find(gid2));

        // search, if erased GID is part of any other coupling pair and if so replace with
        for (const auto& nodal_coupling_inner_loop : coupling_pairs)
        {
          std::pair<int, int> new_pair;
          if (nodal_coupling_inner_loop.first == gid2 or nodal_coupling_inner_loop.second == gid2)
          {
            // new pair: gid1 + (gid, that is not gid2 from inner loop)
            const int other_gid = nodal_coupling_inner_loop.first == gid2
                                      ? nodal_coupling_inner_loop.second
                                      : nodal_coupling_inner_loop.first;
            new_pair = std::make_pair(gid1, other_gid);
            coupling_pairs.erase(
                std::find(coupling_pairs.begin(), coupling_pairs.end(), nodal_coupling_inner_loop));
            coupling_pairs.emplace_back(new_pair);
          }
        }
      }
      else if (master_gids.find(gid1) != master_gids.end() or
               master_gids.find(gid2) != master_gids.end())
      {
        // match in one -> do nothing
      }
      else
      {
        // no  match -> add one random entry
        master_gids.insert(gid2);
      }
    }
  }

  return master_gids;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::map<int, int> SSI::ManifoldMeshTyingStrategyBase::DefineMasterSlavePairing(
    const std::vector<std::pair<int, int>>& coupling_pairs, std::set<int> master_gids)
{
  std::map<int, int> coupling_pair;
  for (const auto& nodal_coupling : coupling_pairs)
  {
    const int gid_node1 = nodal_coupling.first;
    const int gid_node2 = nodal_coupling.second;
    if (gid_node1 != gid_node2)
    {
      if (master_gids.find(gid_node1) != master_gids.end())
        coupling_pair.insert(std::make_pair(gid_node2, gid_node1));
      else if (master_gids.find(gid_node2) != master_gids.end())
        coupling_pair.insert(std::make_pair(gid_node1, gid_node2));
      else
        dserror("Could not find master GID.");
    }
  }
  return coupling_pair;
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
  for (int i = 0; i < ssi_maps_->BlockMapScaTraManifold()->NumMaps(); ++i)
  {
    partial_maps_condensed_block_dof_map.emplace_back(
        LINALG::IntersectMap(*condensed_dof_map_, *ssi_maps_->BlockMapScaTraManifold()->Map(i)));
  }

  condensed_block_dof_map_ = Teuchos::rcp(
      new LINALG::MultiMapExtractor(*condensed_dof_map_, partial_maps_condensed_block_dof_map));

  // couple meshyting_handler_ and condensed_block_dof_map_ to meshtying_block_handler_
  for (const auto& meshtying : MeshTyingHandler())
  {
    auto coupling_adapter = meshtying.first;
    auto slave_dof_map = coupling_adapter->SlaveDofMap();

    std::vector<Teuchos::RCP<const Epetra_Map>> partial_maps_slave_block_dof_map;
    for (int i = 0; i < ssi_maps_->BlockMapScaTraManifold()->NumMaps(); ++i)
    {
      partial_maps_slave_block_dof_map.emplace_back(
          LINALG::IntersectMap(*slave_dof_map, *ssi_maps_->BlockMapScaTraManifold()->Map(i)));
    }
    auto slave_dof_block_map = Teuchos::rcp(
        new LINALG::MultiMapExtractor(*slave_dof_map, partial_maps_slave_block_dof_map));

    meshtying_block_handler_.emplace_back(std::make_pair(meshtying, slave_dof_block_map));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBase::ApplyMeshTyingToManifoldRHS(
    Teuchos::RCP<Epetra_Vector> rhs_manifold)
{
  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : meshtying_handler_)
    {
      auto multimap = meshtying.second;
      auto coupling_adapter = meshtying.first;

      auto slave_dofs = multimap->ExtractVector(rhs_manifold, 2);
      auto slave_to_master_dofs = coupling_adapter->SlaveToMaster(slave_dofs);
      multimap->AddVector(slave_to_master_dofs, 1, rhs_manifold);
      multimap->PutScalar(*rhs_manifold, 2, 0.0);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::ApplyMeshtyingToManifoldMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_matrix);

  // add derivs. of interior/master dofs. w.r.t. interior/master dofs
  LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *condensed_dof_map_,
      *condensed_dof_map_, 1.0, nullptr, nullptr, *ssi_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs. w.r.t. slave dofs
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *cond_slave_dof_map, 1.0, &converter, &converter, *ssi_manifold_sparse, true, true);
      // add derivs. of slave dofs. w.r.t. interior/master dofs
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *condensed_dof_map_, 1.0, &converter, nullptr, *ssi_manifold_sparse, true, true);
      // add derivs. of interior/master dofs w.r.t. slave dofs
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *condensed_dof_map_,
          *cond_slave_dof_map, 1.0, nullptr, &converter, *ssi_manifold_sparse, true, true);
    }

    // Finalize: put 1.0 on main diag of slave dofs
    const double one = 1.0;
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto slave_dof_map = coupling_adapter->SlaveDofMap();
      for (int doflid_slave = 0; doflid_slave < slave_dof_map->NumMyElements(); ++doflid_slave)
      {
        // extract global ID of current slave-side row
        const int dofgid_slave = slave_dof_map->GID(doflid_slave);
        if (dofgid_slave < 0) dserror("Local ID not found!");

        // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column indices
        if (ssi_manifold_sparse->Filled())
        {
          const int rowlid_slave = ssi_manifold_sparse->RowMap().LID(dofgid_slave);
          if (rowlid_slave < 0) dserror("Global ID not found!");
          if (ssi_manifold_sparse->EpetraMatrix()->ReplaceMyValues(
                  rowlid_slave, 1, &one, &rowlid_slave))
            dserror("ReplaceMyValues failed!");
        }

        // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and column
        // indices
        else if (ssi_manifold_sparse->EpetraMatrix()->InsertGlobalValues(
                     dofgid_slave, 1, &one, &dofgid_slave))
          dserror("InsertGlobalValues failed!");
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::ApplyMeshtyingToManifoldMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_block = LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
      // add derivs. of interior/master dofs. w.r.t. interior/master dofs
      LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
          *condensed_block_dof_map_->Map(row), *condensed_block_dof_map_->Map(col), 1.0, nullptr,
          nullptr, ssi_manifold_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : MeshTyingBlockHandler())
        {
          auto meshtying = block_meshtying.first;
          auto coupling_adapter = meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

          // add derivs. of slave dofs. w.r.t. slave dofs
          LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map->Map(row), *cond_block_slave_dof_map->Map(col), 1.0,
              &converter, &converter, ssi_manifold_block->Matrix(row, col), true, true);
          // add derivs. of slave dofs. w.r.t. interior/master dofs
          LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map->Map(row), *condensed_block_dof_map_->Map(col), 1.0,
              &converter, nullptr, ssi_manifold_block->Matrix(row, col), true, true);
          // add derivs. of interior/master dofs w.r.t. slave dofs
          LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *condensed_block_dof_map_->Map(row), *cond_block_slave_dof_map->Map(col), 1.0,
              nullptr, &converter, ssi_manifold_block->Matrix(row, col), true, true);
        }

        // Finalize: put 1.0 on main diag of slave dofs
        const double one = 1.0;
        if (row == col)
        {
          for (const auto& meshtying : MeshTyingHandler())
          {
            auto coupling_adapter = meshtying.first;

            auto slave_dof_map = coupling_adapter->SlaveDofMap();
            for (int doflid_slave = 0; doflid_slave < slave_dof_map->NumMyElements();
                 ++doflid_slave)
            {
              // extract global ID of current slave-side row
              const int dofgid_slave = slave_dof_map->GID(doflid_slave);
              if (dofgid_slave < 0) dserror("Local ID not found!");

              // apply pseudo Dirichlet conditions to filled matrix, i.e., to local row and column
              // indices
              if (ssi_manifold_block->Matrix(row, row).Filled())
              {
                const int rowlid_slave =
                    ssi_manifold_block->Matrix(row, row).RowMap().LID(dofgid_slave);
                if (rowlid_slave < 0) dserror("Global ID not found!");
                if (ssi_manifold_block->Matrix(row, row).EpetraMatrix()->ReplaceMyValues(
                        rowlid_slave, 1, &one, &rowlid_slave))
                  dserror("ReplaceMyValues failed!");
              }

              // apply pseudo Dirichlet conditions to unfilled matrix, i.e., to global row and
              // column indices
              else if (ssi_manifold_block->Matrix(row, row).EpetraMatrix()->InsertGlobalValues(
                           dofgid_slave, 1, &one, &dofgid_slave))
                dserror("InsertGlobalValues failed!");
            }
          }
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::ApplyMeshtyingToManifoldScatraMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_scatra_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_scatra_matrix);

  // add derivs. of interior/master dofs w.r.t. scatra dofs
  LINALG::MatrixLogicalSplitAndTransform()(*manifold_scatra_sparse, *condensed_dof_map_,
      *ssi_maps_->ScaTraDofRowMap(), 1.0, nullptr, nullptr, *ssi_manifold_scatra_sparse, true,
      true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs w.r.t. scatra dofs
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_scatra_sparse, *cond_slave_dof_map,
          *ssi_maps_->ScaTraDofRowMap(), 1.0, &converter, nullptr, *ssi_manifold_scatra_sparse,
          true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::ApplyMeshtyingToManifoldScatraMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_scatra_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_scatra_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < ssi_maps_->BlockMapScaTra()->NumMaps(); ++col)
    {
      // add derivs. of interior/master dofs w.r.t. scatra dofs
      LINALG::MatrixLogicalSplitAndTransform()(manifold_scatra_block->Matrix(row, col),
          *condensed_block_dof_map_->Map(row), *ssi_maps_->BlockMapScaTra()->Map(col), 1.0, nullptr,
          nullptr, ssi_manifold_scatra_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : MeshTyingBlockHandler())
        {
          auto meshtying = block_meshtying.first;
          auto coupling_adapter = meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

          // add derivs. of slave dofs w.r.t. scatra dofs
          LINALG::MatrixLogicalSplitAndTransform()(manifold_scatra_block->Matrix(row, col),
              *cond_block_slave_dof_map->Map(row), *ssi_maps_->BlockMapScaTra()->Map(col), 1.0,
              &converter, nullptr, ssi_manifold_scatra_block->Matrix(row, col), true, true);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::ApplyMeshtyingToManifoldStructureMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_structure_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_structure_matrix, const bool do_uncomplete)
{
  auto ssi_manifold_structure_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_maps_->ScaTraManifoldDofRowMap());

  // add derivs. of interior/master dofs w.r.t. structure dofs
  LINALG::MatrixLogicalSplitAndTransform()(*ssi_manifold_structure_sparse, *condensed_dof_map_,
      *ssi_maps_->StructureDofRowMap(), 1.0, nullptr, nullptr, *temp_manifold_structure, true,
      true);
  LINALG::MatrixLogicalSplitAndTransform()(*manifold_structure_sparse, *condensed_dof_map_,
      *ssi_maps_->StructureDofRowMap(), 1.0, nullptr, nullptr, *temp_manifold_structure, true,
      true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of slave dofs w.r.t. structure dofs
      LINALG::MatrixLogicalSplitAndTransform()(*ssi_manifold_structure_sparse, *cond_slave_dof_map,
          *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr, *temp_manifold_structure,
          true, true);
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_structure_sparse, *cond_slave_dof_map,
          *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr, *temp_manifold_structure,
          true, true);
    }
  }

  ssi_manifold_structure_sparse->Zero();
  if (do_uncomplete) ssi_manifold_structure_sparse->UnComplete();
  temp_manifold_structure->Complete(
      *ssi_maps_->StructureDofRowMap(), *ssi_maps_->ScaTraManifoldDofRowMap());
  ssi_manifold_structure_sparse->Add(*temp_manifold_structure, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::ApplyMeshtyingToManifoldStructureMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_structure_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> manifold_structure_matrix, const bool do_uncomplete)
{
  auto ssi_manifold_structure_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
      ssi_maps_->BlockMapScaTraManifold(), ssi_maps_->BlockMapStructure());

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    // add derivs. of interior/master dofs w.r.t. structure dofs
    LINALG::MatrixLogicalSplitAndTransform()(ssi_manifold_structure_block->Matrix(row, 0),
        *condensed_block_dof_map_->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, nullptr,
        nullptr, temp_manifold_structure->Matrix(row, 0), true, true);
    LINALG::MatrixLogicalSplitAndTransform()(manifold_structure_block->Matrix(row, 0),
        *condensed_block_dof_map_->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, nullptr,
        nullptr, temp_manifold_structure->Matrix(row, 0), true, true);

    if (is_manifold_meshtying_)
    {
      for (const auto& block_meshtying : MeshTyingBlockHandler())
      {
        auto meshtying = block_meshtying.first;
        auto coupling_adapter = meshtying.first;
        auto cond_block_slave_dof_map = block_meshtying.second;
        auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

        // add derivs. of slave dofs w.r.t. structure dofs
        LINALG::MatrixLogicalSplitAndTransform()(ssi_manifold_structure_block->Matrix(row, 0),
            *cond_block_slave_dof_map->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, &converter,
            nullptr, temp_manifold_structure->Matrix(row, 0), true, true);
        LINALG::MatrixLogicalSplitAndTransform()(manifold_structure_block->Matrix(row, 0),
            *cond_block_slave_dof_map->Map(row), *ssi_maps_->StructureDofRowMap(), 1.0, &converter,
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
void SSI::ManifoldMeshTyingStrategySparse::ApplyMeshtyingToScatraManifoldMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_scatra_manifold_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_manifold_matrix, const bool do_uncomplete)
{
  auto ssi_scatra_manifold_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_sparse =
      LINALG::CastToConstSparseMatrixAndCheckSuccess(scatra_manifold_matrix);

  if (do_uncomplete) ssi_scatra_manifold_sparse->UnComplete();

  // add derivs. of scatra w.r.t. interior/master dofs
  LINALG::MatrixLogicalSplitAndTransform()(*scatra_manifold_sparse, *ssi_maps_->ScaTraDofRowMap(),
      *condensed_dof_map_, 1.0, nullptr, nullptr, *ssi_scatra_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      // add derivs. of scatra w.r.t. slave dofs
      LINALG::MatrixLogicalSplitAndTransform()(*scatra_manifold_sparse,
          *ssi_maps_->ScaTraDofRowMap(), *cond_slave_dof_map, 1.0, nullptr, &converter,
          *ssi_scatra_manifold_sparse, true, true);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::ApplyMeshtyingToScatraManifoldMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_scatra_manifold_matrix,
    Teuchos::RCP<const LINALG::SparseOperator> scatra_manifold_matrix, const bool do_uncomplete)
{
  auto ssi_scatra_manifold_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_block =
      LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(scatra_manifold_matrix);

  if (do_uncomplete) ssi_scatra_manifold_block->UnComplete();

  for (int row = 0; row < ssi_maps_->BlockMapScaTra()->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
      // add derivs. of scatra w.r.t. interior/master dofs
      LINALG::MatrixLogicalSplitAndTransform()(scatra_manifold_block->Matrix(row, col),
          *ssi_maps_->BlockMapScaTra()->Map(row), *condensed_block_dof_map_->Map(col), 1.0, nullptr,
          nullptr, ssi_scatra_manifold_block->Matrix(row, col), true, true);

      if (is_manifold_meshtying_)
      {
        for (const auto& block_meshtying : MeshTyingBlockHandler())
        {
          auto meshtying = block_meshtying.first;
          auto coupling_adapter = meshtying.first;
          auto cond_block_slave_dof_map = block_meshtying.second;
          auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

          // add derivs. of scatra w.r.t. slave dofs
          LINALG::MatrixLogicalSplitAndTransform()(scatra_manifold_block->Matrix(row, col),
              *ssi_maps_->BlockMapScaTra()->Map(row), *cond_block_slave_dof_map->Map(col), 1.0,
              nullptr, &converter, ssi_scatra_manifold_block->Matrix(row, col), true, true);
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> SSI::BuildManifoldMeshTyingStrategy(
    Teuchos::RCP<DRT::Discretization> scatra_manifold_dis, Teuchos::RCP<UTILS::SSIMaps> ssi_maps,
    const bool is_manifold_meshtying, LINALG::MatrixType matrixtype_manifold)
{
  Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> meshtyingstrategy = Teuchos::null;

  switch (matrixtype_manifold)
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      meshtyingstrategy = Teuchos::rcp(new SSI::ManifoldMeshTyingStrategyBlock(
          scatra_manifold_dis, ssi_maps, is_manifold_meshtying));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      meshtyingstrategy = Teuchos::rcp(new SSI::ManifoldMeshTyingStrategySparse(
          scatra_manifold_dis, ssi_maps, is_manifold_meshtying));
      break;
    }

    default:
    {
      dserror("unknown matrix type of Manifold field");
      break;
    }
  }

  return meshtyingstrategy;
}