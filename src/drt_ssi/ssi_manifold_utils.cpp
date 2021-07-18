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
  // definition of lambda function to get size of graph for given matrix
  auto graph_size = [](Teuchos::RCP<LINALG::SparseMatrix> matrix) {
    return matrix->EpetraMatrix()->Graph().NumGlobalEntries();
  };

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
  // data types to handle coupling between slave and master in coupling_vec
  // coupling_type: key: master gid, value: slave gid and distance between both
  // key_type: key: condition a, value: condition b
  using coupling_type = std::map<int, std::pair<int, double>>;
  using key_type = std::pair<int, int>;
  std::map<key_type, coupling_type> coupling_vec;

  std::vector<DRT::Condition*> manifold_conditions(0, nullptr);
  scatra_manifold_dis->GetCondition("SSISurfaceManifold", manifold_conditions);

  // fill coupling_vec with slave master pairing using NodeMatchingOctree to find matching nodes
  // between all SSISurfaceManifold conditions (nodes from all procs)
  for (int a = 0; a < static_cast<int>(manifold_conditions.size()); ++a)
  {
    auto* manifold_condition_a = manifold_conditions.at(a);

    std::vector<int> inodegidvec_a = *manifold_condition_a->Nodes();

    DRT::UTILS::NodeMatchingOctree tree = DRT::UTILS::NodeMatchingOctree();
    tree.Init(*scatra_manifold_dis, inodegidvec_a, 150, 1.0e-8);
    tree.Setup();

    for (int b = a + 1; b < static_cast<int>(manifold_conditions.size()); ++b)
    {
      auto* manifold_condition_b = manifold_conditions.at(b);

      std::vector<int> inodegidvec_b = *manifold_condition_b->Nodes();

      coupling_type val;
      key_type key = std::make_pair(a, b);
      coupling_vec.emplace(key, val);

      tree.FindMatch(*scatra_manifold_dis, inodegidvec_b, coupling_vec.at(key));
    }
  }

  // construct set of unique master GIDs within coupling_vec -> master_gids
  std::set<int> master_gids;

  for (const auto& coupling : coupling_vec)
  {
    for (const auto& pair : coupling.second)
    {
      const int gid_node1 = pair.first;
      const int gid_node2 = pair.second.first;

      if (gid_node1 == gid_node2) master_gids.insert(gid_node1);
    }
  }

  // assigning slave GIDs (unique key) to master GIDs (non unique values, a master node can have
  // multiple slave nodes) -> coupling_pair
  std::map<int, int> coupling_pair;

  for (const auto& coupling : coupling_vec)
  {
    for (const auto& pair : coupling.second)
    {
      const int gid_node1 = pair.first;
      const int gid_node2 = pair.second.first;
      if (gid_node1 != gid_node2)
      {
        if (master_gids.find(gid_node1) != master_gids.end())
          coupling_pair.insert(std::make_pair(gid_node2, gid_node1));
        else if (master_gids.find(gid_node2) != master_gids.end())
          coupling_pair.insert(std::make_pair(gid_node1, gid_node2));
      }
    }
  }

  // get number of slave nodes per master node -> max. number gives number of needed adapters
  int my_max_adapters = 0;
  std::map<int, int> master_coupling;
  for (auto pair : coupling_pair)
  {
    const int master_node_gid = pair.second;
    if (master_coupling.empty())
    {
      master_coupling.insert(std::make_pair(master_node_gid, 1));
      my_max_adapters = 1;
    }
    else
    {
      if (master_coupling.find(master_node_gid) != master_coupling.end())
      {
        master_coupling[master_node_gid]++;
        if (my_max_adapters < master_coupling[master_node_gid])
          my_max_adapters = master_coupling[master_node_gid];
      }
      else
        master_coupling.insert(std::make_pair(master_node_gid, 1));
    }
  }

  int glob_max_adapters = 0;
  scatra_manifold_dis->Comm().MaxAll(&my_max_adapters, &glob_max_adapters, 1);

  // setup coupling adapters
  for (int iadapter = 0; iadapter < glob_max_adapters; ++iadapter)
  {
    std::vector<int> inodegidvec_master;
    std::vector<int> inodegidvec_slave;

    for (auto master_gid : master_gids)
    {
      // check if this master node has iadapter + 1 slave nodes
      if (!master_coupling.empty())
      {
        if (master_coupling.at(master_gid) == iadapter + 1)
        {
          DRT::UTILS::AddOwnedNodeGID(scatra_manifold_dis, master_gid, inodegidvec_master);

          int counter = 0;
          for (auto pair : coupling_pair)
          {
            const int master_gid_coupling = pair.second;

            if (master_gid_coupling == master_gid and counter == iadapter)
            {
              if (counter == iadapter)
              {
                const int slave_gid_coupling = pair.first;

                DRT::UTILS::AddOwnedNodeGID(
                    scatra_manifold_dis, slave_gid_coupling, inodegidvec_slave);
              }
              counter++;
            }
          }
        }
      }
    }

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
    Teuchos::RCP<LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(manifold_matrix);

  LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *condensed_dof_map_,
      *condensed_dof_map_, 1.0, nullptr, nullptr, *ssi_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

      LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *cond_slave_dof_map, 1.0, &converter, &converter, *ssi_manifold_sparse, true, true);
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_sparse, *cond_slave_dof_map,
          *condensed_dof_map_, 1.0, &converter, nullptr, *ssi_manifold_sparse, true, true);
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
    Teuchos::RCP<LINALG::SparseOperator> manifold_matrix)
{
  auto ssi_manifold_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_matrix);
  auto manifold_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifold_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
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

          LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map->Map(row), *cond_block_slave_dof_map->Map(col), 1.0,
              &converter, &converter, ssi_manifold_block->Matrix(row, col), true, true);
          LINALG::MatrixLogicalSplitAndTransform()(manifold_block->Matrix(row, col),
              *cond_block_slave_dof_map->Map(row), *condensed_block_dof_map_->Map(col), 1.0,
              &converter, nullptr, ssi_manifold_block->Matrix(row, col), true, true);
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
    Teuchos::RCP<LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(manifold_scatra_matrix);

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
    Teuchos::RCP<LINALG::SparseOperator> manifold_scatra_matrix)
{
  auto ssi_manifold_scatra_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_scatra_matrix);
  auto manifold_scatra_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifold_scatra_matrix);

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
    for (int col = 0; col < ssi_maps_->BlockMapScaTra()->NumMaps(); ++col)
    {
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
    Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix)
{
  auto ssi_manifold_structure_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure =
      SSI::UTILS::SSIMatrices::SetupSparseMatrix(ssi_maps_->ScaTraManifoldDofRowMap());
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

      // manifold - structure
      LINALG::MatrixLogicalSplitAndTransform()(*ssi_manifold_structure_sparse, *cond_slave_dof_map,
          *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr, *temp_manifold_structure,
          true, true);
      LINALG::MatrixLogicalSplitAndTransform()(*manifold_structure_sparse, *cond_slave_dof_map,
          *ssi_maps_->StructureDofRowMap(), 1.0, &converter, nullptr, *temp_manifold_structure,
          true, true);
    }
  }

  ssi_manifold_structure_sparse->Zero();
  temp_manifold_structure->Complete(
      *ssi_maps_->StructureDofRowMap(), *ssi_maps_->ScaTraManifoldDofRowMap());
  ssi_manifold_structure_sparse->UnComplete();
  ssi_manifold_structure_sparse->Add(*temp_manifold_structure, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategyBlock::ApplyMeshtyingToManifoldStructureMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_manifold_structure_matrix,
    Teuchos::RCP<LINALG::SparseOperator> manifold_structure_matrix)
{
  auto ssi_manifold_structure_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_manifold_structure_matrix);
  auto manifold_structure_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(manifold_structure_matrix);

  auto temp_manifold_structure = SSI::UTILS::SSIMatrices::SetupBlockMatrix(
      ssi_maps_->BlockMapScaTraManifold(), ssi_maps_->BlockMapStructure());

  for (int row = 0; row < condensed_block_dof_map_->NumMaps(); ++row)
  {
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

        // manifold - structure
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
  temp_manifold_structure->Complete();
  ssi_manifold_structure_block->UnComplete();
  ssi_manifold_structure_block->Add(*temp_manifold_structure, false, 1.0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ManifoldMeshTyingStrategySparse::ApplyMeshtyingToScatraManifoldMatrix(
    Teuchos::RCP<LINALG::SparseOperator> ssi_scatra_manifold_matrix,
    Teuchos::RCP<LINALG::SparseOperator> scatra_manifold_matrix)
{
  auto ssi_scatra_manifold_sparse =
      LINALG::CastToSparseMatrixAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(scatra_manifold_matrix);

  ssi_scatra_manifold_sparse->UnComplete();
  LINALG::MatrixLogicalSplitAndTransform()(*scatra_manifold_sparse, *ssi_maps_->ScaTraDofRowMap(),
      *condensed_dof_map_, 1.0, nullptr, nullptr, *ssi_scatra_manifold_sparse, true, true);

  if (is_manifold_meshtying_)
  {
    for (const auto& meshtying : MeshTyingHandler())
    {
      auto coupling_adapter = meshtying.first;

      auto cond_slave_dof_map = coupling_adapter->SlaveDofMap();
      auto converter = ADAPTER::CouplingSlaveConverter(*coupling_adapter);

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
    Teuchos::RCP<LINALG::SparseOperator> scatra_manifold_matrix)
{
  auto ssi_scatra_manifold_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(ssi_scatra_manifold_matrix);
  auto scatra_manifold_block =
      LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(scatra_manifold_matrix);

  ssi_scatra_manifold_block->UnComplete();

  for (int row = 0; row < ssi_maps_->BlockMapScaTra()->NumMaps(); ++row)
  {
    for (int col = 0; col < condensed_block_dof_map_->NumMaps(); ++col)
    {
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