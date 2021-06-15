/*----------------------------------------------------------------------*/
/*! \file
 \brief Evaluates flux between ScaTra and ScaTra on manifolds incl. coupling matrices

 \level 2


 *------------------------------------------------------------------------------------------------*/

#include "ssi_manifold_flux_evaluator.H"

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
#include "../drt_lib/drt_utils_gid_vector.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
SSI::ManifoldScaTraCoupling::ManifoldScaTraCoupling(Teuchos::RCP<DRT::Discretization> manifolddis,
    Teuchos::RCP<DRT::Discretization> scatradis, DRT::Condition* condition_manifold,
    DRT::Condition* condition_kinetics, const int ndof_per_node)
    : condition_kinetics_(condition_kinetics),
      condition_manifold_(condition_manifold),
      coupling_adapter_(Teuchos::rcp(new ADAPTER::Coupling())),
      evaluate_master_side_(
          !(DRT::UTILS::HaveSameNodes(condition_manifold, condition_kinetics, false))),
      manifold_conditionID_(condition_manifold->GetInt("ConditionID")),
      manifold_map_extractor_(Teuchos::null),
      scatra_map_extractor_(Teuchos::null),
      size_matrix_graph_(),
      slave_converter_(Teuchos::null)
{
  std::vector<int> inodegidvec_manifold;
  DRT::UTILS::AddOwnedNodeGIDVector(
      manifolddis, *condition_manifold->Nodes(), inodegidvec_manifold);

  std::vector<int> inodegidvec_scatra;
  DRT::UTILS::AddOwnedNodeGIDVector(scatradis, *condition_kinetics->Nodes(), inodegidvec_scatra);

  coupling_adapter_->SetupCoupling(*scatradis, *manifolddis, inodegidvec_scatra,
      inodegidvec_manifold, ndof_per_node, true, 1.0e-8);
  slave_converter_ = Teuchos::rcp(new ADAPTER::CouplingSlaveConverter(*coupling_adapter_));

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
bool SSI::ManifoldScaTraCoupling::CheckAndSetSizeGraph(const BlockMatrixType block, const int size)
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
          DRT::Problem::Instance()->SSIControlParams().sublist("MANIFOLD"), "ADD_MANIFOLD")),
      full_map_manifold_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold))),
      full_map_scatra_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport))),
      full_map_structure_(ssi_mono.MapsSubProblems()->Map(
          UTILS::SSIMaps::GetProblemPosition(Subproblem::structure))),
      icoup_structure_(ssi_mono.SSIStructureMeshTying()->InterfaceCouplingAdapterStructure()),
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
      systemmatrix_manifold_(Teuchos::null),
      systemmatrix_scatra_(Teuchos::null)
{
  // safety check befor setup of coupling
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
  systemmatrix_scatra_->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CompleteSystemMatrixScaTra()
{
  systemmatrix_manifold_->Complete();
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

    EvaluateManifoldSide(scatra_manifold_coupling);

    CopyScaTraManifoldScaTraMasterSide(scatra_manifold_coupling);

    // This is needed because the graph of the matrices could change from step to step in case we
    // have zero flux (and then zero entries in the matrices)
    UnCompleteMatricesIfNecessary(scatra_manifold_coupling);

    AddConditionContribution();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateManifoldSide(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  // First: Set parameters to elements
  PreEvaluate(scatra_manifold_coupling);

  // Second: Evaluate condition
  {
    // manifold-scatra coupling matrix evaluated on manifold side
    auto matrix_manifold_scatra_manifold_side =
        Teuchos::rcp(new LINALG::SparseMatrix(*full_map_manifold_, 27, false, true));

    Teuchos::ParameterList condparams;

    condparams.set<int>("action", SCATRA::calc_scatra_manifold_flux);

    condparams.set<int>("ndsdisp", 1);

    scatra_manifold_->ScaTraField()->Discretization()->ClearState();

    scatra_manifold_->ScaTraField()->AddTimeIntegrationSpecificVectors();

    // Evaluation of RHS and scatra-manifold coupling matrices
    {
      condparams.set<int>(
          "differentiationtype", static_cast<int>(SCATRA::DifferentiationType::elch));

      DRT::AssembleStrategy strategymanifold(0, 0, systemmatrix_manifold_cond_,
          matrix_manifold_scatra_manifold_side, rhs_manifold_cond_, Teuchos::null, Teuchos::null);

      scatra_manifold_->ScaTraField()->Discretization()->EvaluateCondition(condparams,
          strategymanifold, "SSISurfaceManifold", scatra_manifold_coupling->ManifoldConditionID());

      systemmatrix_manifold_cond_->Complete();
      matrix_manifold_scatra_manifold_side->Complete();

      // column dofs are so far on manifold dis. They are transformed to scatra dis
      LINALG::MatrixLogicalSplitAndTransform()(*matrix_manifold_scatra_manifold_side,
          *full_map_manifold_, *full_map_scatra_, 1.0, nullptr,
          &*scatra_manifold_coupling->SlaveConverter(), *matrix_manifold_scatra_cond_, true, true);

      matrix_manifold_scatra_cond_->Complete(*full_map_scatra_, *full_map_manifold_);
    }

    // Evaluation of linearization w.r.t. displacement
    {
      condparams.set<int>(
          "differentiationtype", static_cast<int>(SCATRA::DifferentiationType::disp));

      auto flux_manifold_scatra_md_cond_slave_side_disp =
          Teuchos::rcp(new LINALG::SparseMatrix(*full_map_manifold_, 27, false, true));

      DRT::AssembleStrategy strategymanifold(0, 1, flux_manifold_scatra_md_cond_slave_side_disp,
          Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

      scatra_manifold_->ScaTraField()->Discretization()->EvaluateCondition(condparams,
          strategymanifold, "SSISurfaceManifold", scatra_manifold_coupling->ManifoldConditionID());

      flux_manifold_scatra_md_cond_slave_side_disp->Complete(
          *full_map_structure_, *full_map_manifold_);

      // Add slave side disp. contributions
      matrix_manifold_structure_cond_->Add(
          *flux_manifold_scatra_md_cond_slave_side_disp, false, 1.0, 1.0);

      // Add master side disp. contributions
      ADAPTER::CouplingSlaveConverter converter(*icoup_structure_);
      LINALG::MatrixLogicalSplitAndTransform()(*flux_manifold_scatra_md_cond_slave_side_disp,
          *full_map_manifold_, *full_map_structure_, 1.0, nullptr, &converter,
          *matrix_manifold_structure_cond_, true, true);

      matrix_manifold_structure_cond_->Complete(*full_map_structure_, *full_map_manifold_);
    }

    scatra_manifold_->ScaTraField()->Discretization()->ClearState();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::CopyScaTraManifoldScaTraMasterSide(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  {
    auto flux_manifold_scatra_m_cond_extract =
        scatra_manifold_coupling->ManifoldMapExtractor()->ExtractCondVector(rhs_manifold_cond_);

    auto flux_manifold_domain_RHS_m_cond_to_s =
        scatra_manifold_coupling->CouplingAdapter()->SlaveToMaster(
            flux_manifold_scatra_m_cond_extract);

    scatra_manifold_coupling->ScaTraMapExtractor()->AddCondVector(
        flux_manifold_domain_RHS_m_cond_to_s, rhs_scatra_cond_);
    rhs_scatra_cond_->Scale(-1.0);
  }

  // djscatra_dmanifold: manifold rows are transformed to scatra side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_manifold_cond_, *full_map_scatra_,
      *full_map_manifold_, -1.0, &*scatra_manifold_coupling->SlaveConverter(), nullptr,
      *matrix_scatra_manifold_cond_, true, true);

  // djscatra_dscatra: manifold rows are transformed to scatra side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*matrix_manifold_scatra_cond_, *full_map_scatra_,
      *full_map_scatra_, -1.0, &*scatra_manifold_coupling->SlaveConverter(), nullptr,
      *systemmatrix_scatra_cond_, true, true);

  matrix_scatra_manifold_cond_->Complete(*full_map_manifold_, *full_map_scatra_);
  systemmatrix_scatra_cond_->Complete();

  // djscatra_dstructure: manifold rows are transformed to scatra side (flux is scaled by -1.0)
  LINALG::MatrixLogicalSplitAndTransform()(*matrix_manifold_structure_cond_, *full_map_scatra_,
      *full_map_structure_, -1.0, &*scatra_manifold_coupling->SlaveConverter(), nullptr,
      *matrix_scatra_structure_cond_, true, true);

  matrix_scatra_structure_cond_->Complete(*full_map_structure_, *full_map_scatra_);
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

  scatra_manifold_->ScaTraField()->Discretization()->SetState(
      "phinp", scatra_manifold_->ScaTraField()->Phinp());

  for (const auto& scatra_manifold_coupling : scatra_manifold_couplings_)
  {
    const int manifoldID = scatra_manifold_coupling->ManifoldConditionID();

    std::vector<double> zero_scalar_vector(scatra_manifold_->ScaTraField()->NumDofPerNode(), 0.0);
    inflow_.insert(std::make_pair(manifoldID, zero_scalar_vector));

    // First: set parameters to elements
    PreEvaluate(scatra_manifold_coupling);

    // Second: evaluate condition
    EvaluateScaTraManifoldInflowIntegral(scatra_manifold_coupling);

    // Third: evaluate domain integral
    EvaluateScaTraManifoldDomainIntegral(scatra_manifold_coupling);
  }
  scatra_manifold_->ScaTraField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateScaTraManifoldDomainIntegral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int manifoldID = scatra_manifold_coupling->ManifoldConditionID();

  // integrate only if not done so far
  if (domainintegral_.find(manifoldID) == domainintegral_.end())
  {
    Teuchos::ParameterList condparams;

    condparams.set<int>("action", SCATRA::calc_domain_integral);

    condparams.set<int>("ndsdisp", 1);

    // integrated domain of this condition
    auto domainintegral_cond = Teuchos::rcp(new Epetra_SerialDenseVector(1));

    scatra_manifold_->ScaTraField()->Discretization()->EvaluateScalars(
        condparams, domainintegral_cond, "SSISurfaceManifold", manifoldID);

    domainintegral_.insert(std::make_pair(manifoldID, domainintegral_cond->Values()[0]));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::EvaluateScaTraManifoldInflowIntegral(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  const int manifoldID = scatra_manifold_coupling->ManifoldConditionID();

  Teuchos::ParameterList condparams;

  condparams.set<int>("action", SCATRA::calc_manifold_inflow);

  condparams.set<int>("ndsdisp", 1);

  // integrated scalars of this condition
  auto inflow_cond =
      Teuchos::rcp(new Epetra_SerialDenseVector(scatra_manifold_->ScaTraField()->NumDofPerNode()));

  scatra_manifold_->ScaTraField()->Discretization()->EvaluateScalars(
      condparams, inflow_cond, "SSISurfaceManifold", manifoldID);

  for (int i = 0; i < inflow_cond->Length(); ++i)
    inflow_.at(scatra_manifold_coupling->ManifoldConditionID()).at(i) += inflow_cond->Values()[i];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScaTraManifoldScaTraFluxEvaluator::PreEvaluate(
    Teuchos::RCP<ManifoldScaTraCoupling> scatra_manifold_coupling)
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", SCATRA::set_elch_scatra_manifold_parameter);

  eleparams.set<bool>("evaluate_master_side", scatra_manifold_coupling->EvaluateMasterSide());

  eleparams.set<int>(
      "kinetic_model", scatra_manifold_coupling->ConditionKinetics()->GetInt("KineticModel"));

  if (scatra_manifold_coupling->ConditionKinetics()->GetInt("KineticModel") ==
      INPAR::SSI::kinetics_constantinterfaceresistance)
  {
    eleparams.set<int>(
        "num_electrons", scatra_manifold_coupling->ConditionKinetics()->GetInt("e-"));
    eleparams.set<double>(
        "resistance", scatra_manifold_coupling->ConditionKinetics()->GetDouble("resistance"));
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
  bool do_uncomplete = false;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::ManifoldScaTra,
          matrix_manifold_scatra_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::ManifoldStructure,
          matrix_manifold_structure_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::ScaTraManifold,
          matrix_scatra_manifold_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::ScaTraStructure,
          matrix_scatra_structure_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::SysMatManifold,
          systemmatrix_manifold_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;
  if (scatra_manifold_coupling->CheckAndSetSizeGraph(BlockMatrixType::SysMatScaTra,
          systemmatrix_scatra_cond_->EpetraMatrix()->MaxNumEntries()))
    do_uncomplete = true;

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