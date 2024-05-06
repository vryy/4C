/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSI
\level 2


 */
/*----------------------------------------------------------------------*/
#include "4C_ssi_monolithic_evaluate_OffDiag.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_lib_assemblestrategy.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_ssi_utils.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraStructureOffDiagCoupling::ScatraStructureOffDiagCoupling(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : block_map_structure_(std::move(block_map_structure)),
      full_map_structure_(std::move(full_map_structure)),
      meshtying_strategy_s2i_(std::move(meshtying_strategy_s2i)),
      scatra_(std::move(scatra)),
      structure_(std::move(structure)),
      ssi_structure_meshtying_(std::move(ssi_structure_meshtying))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraManifoldStructureOffDiagCoupling::ScatraManifoldStructureOffDiagCoupling(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : ScatraStructureOffDiagCoupling(std::move(block_map_structure), std::move(full_map_structure),
          std::move(ssi_structure_meshtying), std::move(meshtying_strategy_s2i), std::move(scatra),
          std::move(structure)),
      scatra_manifold_(std::move(scatra_manifold))
{
}


/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructureblock)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_scatra_mono_odblock_mesh, eleparams);

  // add state vectors to scalar transport discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-structure matrix block
  DRT::AssembleStrategy strategyscatrastructure(
      0,  // row assembly based on number of dofset associated with scalar transport dofs on
      // scalar transport discretization
      1,  // column assembly based on number of dofset associated with structural dofs on scalar
      // transport discretization
      scatrastructureblock,  // scatra-structure matrix block
      Teuchos::null,         // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble scatra-structure matrix block
  ScaTraField()->Discretization()->Evaluate(eleparams, strategyscatrastructure);
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraManifoldStructureDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatramanifoldstructureblock)
{
  FOUR_C_THROW("not implemented");
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraManifoldStructureOffDiagCoupling::EvaluateOffDiagBlockScatraManifoldStructureDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatramanifoldstructureblock)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_scatra_mono_odblock_mesh, eleparams);

  // add state vectors to scalar transport discretization
  scatra_manifold_->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-structure matrix block
  DRT::AssembleStrategy strategyscatrastructure(
      0,  // row assembly based on number of dofset associated with scalar transport dofs on
      // scalar transport discretization
      1,  // column assembly based on number of dofset associated with structural dofs on scalar
      // transport discretization
      scatramanifoldstructureblock,  // scatra-structure matrix block
      Teuchos::null,                 // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble scatra-structure matrix block
  scatra_manifold_->Discretization()->Evaluate(eleparams, strategyscatrastructure);
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureInterface(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrastructureinterface)
{
  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<CORE::LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      slavematrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *block_map_structure_, meshtying_strategy_s2i_->BlockMapsMaster(), 81, false, true));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true));
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // evaluate symmetric interface contributions on slave side
  EvaluateScatraStructureSymmetricInterfaceContributionsSlaveSide(slavematrix);

  // copy symmetric interface contributions from slave side to master side
  CopySlaveToMasterScatraStructureSymmetricInterfaceContributions(slavematrix, mastermatrix);

  // evaluate non-symmetric interface contributions
  EvaluateScatraStructureNonSymmetricInterfaceContributionsSlaveSide(slavematrix, mastermatrix);

  // add contributions from slave side and master side
  scatrastructureinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrastructureinterface->Add(*mastermatrix, false, 1.0, 1.0);
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockStructureScatraDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> structurescatradomain) const
{
  // create parameter list for element evaluation and fill it
  Teuchos::ParameterList eleparams;
  // set action
  eleparams.set("action", "calc_struct_stiffscalar");

  // linearization of structural residuals w.r.t. elch
  eleparams.set<int>("differentiationtype", static_cast<int>(STR::DifferentiationType::elch));

  // set time
  eleparams.set<double>("total time", structure_->Time());
  // set numscatradofspernode
  eleparams.set<int>("numscatradofspernode", ScaTraField()->NumDofPerNode());

  // remove state vectors from structure discretization
  structure_->Discretization()->ClearState();

  // set the current displacement state vector
  structure_->Discretization()->SetState("displacement", structure_->Dispnp());

  // create strategy for assembly of structure-scatra matrix block
  DRT::AssembleStrategy strategystructurescatra(
      0,  // row assembly based on number of dofset associated with structure dofs on structural
      // discretization
      1,  // column assembly based on number of dofset associated with scalar transport dofs on
      // structural discretization
      structurescatradomain,  // structure-scatra matrix block
      Teuchos::null,          // no additional matrices or vectors needed
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble structure-scatra matrix block
  structure_->Discretization()->Evaluate(eleparams, strategystructurescatra);

  // need to scale structurescatrablock_ with 'timefac' (e.g. with theta for OST-scheme) to get
  // correct implementation
  const double timeintparam = structure_->TimIntParam();
  // scale with theta
  structurescatradomain->Scale(1.0 - timeintparam);
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::
    CopySlaveToMasterScatraStructureSymmetricInterfaceContributions(
        Teuchos::RCP<const CORE::LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<CORE::LINALG::SparseOperator>& mastermatrix)
{
  mastermatrix->Zero();
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      const int numberscatrablocks = ScaTraField()->BlockMaps()->NumMaps();

      // cast master and slave matrix
      auto blockslavematrix =
          CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(slavematrix);
      auto blockmastermatrix =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      auto mastermatrixsparse = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true));

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
      {
        for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
        {
          auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
          auto slave_side_converter_struct = meshtying->SlaveSideConverter();

          auto slave_side_converter_scatra =
              CORE::ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

          CORE::LINALG::MatrixLogicalSplitAndTransform()(blockslavematrix->Matrix(iblock, 0),
              *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
              &slave_side_converter_scatra, &(*slave_side_converter_struct), *mastermatrixsparse,
              true, true);
        }
      }

      // finalize auxiliary system matrix
      mastermatrixsparse->Complete(*FullMapStructure(), *ScaTraField()->DofRowMap());

      // split auxiliary system matrix and assemble into scatra-structure matrix block
      auto mastermatrix_split = mastermatrixsparse->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *block_map_structure_, *ScaTraField()->BlockMaps());
      mastermatrix_split->Complete();
      blockmastermatrix->Add(*mastermatrix_split, false, 1.0, 1.0);

      mastermatrix->Complete();

      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      // cast master and slave matrix
      auto sparseslavematrix = CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(slavematrix);
      auto sparsemastermatrix = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(mastermatrix);

      // copy slave side values to master side and scale with minus 1. Insert into
      // scatrastructureinterface_sparse
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto slave_side_converter_struct = meshtying->SlaveSideConverter();
        auto slave_side_converter_scatra =
            CORE::ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

        CORE::LINALG::MatrixLogicalSplitAndTransform()(*sparseslavematrix,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
            &slave_side_converter_scatra, &(*slave_side_converter_struct), *sparsemastermatrix,
            true, true);
      }
      // finalize
      mastermatrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap());
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
void SSI::ScatraStructureOffDiagCoupling::
    EvaluateScatraStructureNonSymmetricInterfaceContributionsSlaveSide(
        Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix,
        Teuchos::RCP<CORE::LINALG::SparseOperator> mastermatrix)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_capacitance_od, condparams);

  // linearization of boundary flux w.r.t. displacement
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
      "differentiationtype", SCATRA::DifferentiationType::disp, condparams);

  // add state vectors to scalar transport discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

  // set up necessary matrices
  Teuchos::RCP<CORE::LINALG::SparseOperator>
      scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix;
  Teuchos::RCP<CORE::LINALG::SparseOperator>
      scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix;
  Teuchos::RCP<CORE::LINALG::SparseOperator>
      scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix;

  if (ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse)
  {
    scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
    scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
    scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
  }
  else
  {
    scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix =
        Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
            *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
    scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix =
        Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
            *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
    scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix =
        Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
            *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
  }

  // create strategy for assembly of auxiliary system matrix
  DRT::AssembleStrategy strategyscatras2istructure(
      0,  // row assembly based on number of dofset associated with scalar transport dofs on
      // scalar transport discretization
      1,  // column assembly based on number of dofset associated with structural dofs on
      // structural discretization
      scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix,
      scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix,
      // no additional vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  for (auto kinetics_slave_cond : meshtying_strategy_s2i_->KineticsConditionsMeshtyingSlaveSide())
  {
    if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") ==
        static_cast<int>(INPAR::S2I::kinetics_butlervolmerreducedcapacitance))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_s2i_->SetConditionSpecificScaTraParameters(*kinetics_slave_cond.second);
      // evaluate the condition
      ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatras2istructure, "S2IKinetics", kinetics_slave_cond.first);
    }
  }

  // finalize scatra-structure matrix block
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::sparse:
    {
      scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());
      scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

      auto scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix_sparse =
          CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(
              scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix);
      auto slavematrix_sparse = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(slavematrix);

      auto scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix_sparse =
          CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(
              scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix);
      auto scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_sparse =
          CORE::LINALG::CastToSparseMatrixAndCheckSuccess(
              scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix);
      auto mastermatrix_sparse = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(mastermatrix);

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();

        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            CORE::ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        // add slave contributions to slave matrix
        CORE::LINALG::MatrixLogicalSplitAndTransform()(
            *scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix_sparse,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter, *slavematrix_sparse, true, true);
        // convert structure slave dofs on scatra discretization to slave dofs on structure
        // discretization
        CORE::LINALG::MatrixLogicalSplitAndTransform()(
            *scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix_sparse,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter,
            *scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_sparse, true,
            true);

        scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_sparse->Complete(
            *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

        auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto slave_side_converter_struct = meshtying->SlaveSideConverter();
        auto slave_side_converter_scatra =
            CORE::ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

        CORE::LINALG::MatrixLogicalSplitAndTransform()(
            *scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_sparse,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, 1.0,
            &slave_side_converter_scatra, &(*slave_side_converter_struct), *mastermatrix_sparse,
            true, true);
      }

      break;
    }

    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix->Complete();
      scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix->Complete();

      slavematrix->UnComplete();
      mastermatrix->UnComplete();

      auto scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix_block =
          CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
              scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix);
      auto slavematrix_block =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(slavematrix);

      auto scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix_block =
          CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
              scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix);
      auto mastermatrix_block =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(mastermatrix);
      auto scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_block =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(
              scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      auto mastermatrixsparse = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true));

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            CORE::ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        for (int iblock = 0; iblock < ScaTraField()->BlockMaps()->NumMaps(); ++iblock)
        {
          auto scatra_slave_flux_structure_slave_dofs_on_scatra_slave_iblock =
              scatra_slave_flux_structure_slave_dofs_on_scatra_slave_matrix_block->Matrix(
                  iblock, 0);
          auto slave_iblock = slavematrix_block->Matrix(iblock, 0);

          auto scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_iblock =
              scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_matrix_block
                  ->Matrix(iblock, 0);
          auto scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_iblock =
              scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_matrix_block->Matrix(
                  iblock, 0);

          auto scatra_block_mapi =
              CORE::LINALG::IntersectMap(*ScaTraField()->BlockMaps()->Map(iblock),
                  *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

          CORE::LINALG::MatrixLogicalSplitAndTransform()(
              scatra_slave_flux_structure_slave_dofs_on_scatra_slave_iblock, *scatra_block_mapi,
              *slave_map, 1.0, nullptr, &slave_slave_converter, slave_iblock, true, true);
          CORE::LINALG::MatrixLogicalSplitAndTransform()(
              scatra_master_flux_on_scatra_slave_structure_slave_dofs_on_scatra_slave_iblock,
              *scatra_block_mapi, *slave_map, 1.0, nullptr, &slave_slave_converter,
              scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_iblock, true, true);

          scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_iblock.Complete();

          auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
          auto slave_side_converter_struct = meshtying->SlaveSideConverter();
          auto slave_side_converter_scatra =
              CORE::ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

          CORE::LINALG::MatrixLogicalSplitAndTransform()(
              scatra_master_flux_on_scatra_slave_dofs_structure_slave_dofs_iblock,
              *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, 1.0,
              &slave_side_converter_scatra, &(*slave_side_converter_struct), *mastermatrixsparse,
              true, true);
        }
      }

      // finalize auxiliary system matrix
      mastermatrixsparse->Complete(*FullMapStructure(), *ScaTraField()->DofRowMap());

      // split auxiliary system matrix and assemble into scatra-structure matrix block
      auto mastermatrix_split = mastermatrixsparse->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *block_map_structure_, *ScaTraField()->BlockMaps());
      mastermatrix_split->Complete();
      mastermatrix_block->Add(*mastermatrix_split, false, 1.0, 1.0);

      mastermatrix->Complete();
      slavematrix->Complete();
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
void SSI::ScatraStructureOffDiagCoupling::
    EvaluateScatraStructureSymmetricInterfaceContributionsSlaveSide(
        Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_od, condparams);

  // linearization of boundary flux w.r.t. displacement
  CORE::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
      "differentiationtype", SCATRA::DifferentiationType::disp, condparams);

  // add state vectors to scalar transport discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

  Teuchos::RCP<CORE::LINALG::SparseOperator> evaluate_matrix;
  if (ScaTraField()->MatrixType() == CORE::LINALG::MatrixType::sparse)
  {
    evaluate_matrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
        *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
  }
  else
  {
    evaluate_matrix =
        Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
            *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
  }

  // create strategy for assembly of auxiliary system matrix
  DRT::AssembleStrategy strategyscatrastructures2i(
      0,  // row assembly based on number of dofset associated with scalar transport dofs on
          // scalar transport discretization
      1,  // column assembly based on number of dofset associated with structural dofs on
          // structural discretization
      evaluate_matrix,  // auxiliary system matrix
      Teuchos::null,    // no additional matrices of vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  for (auto kinetics_slave_cond : meshtying_strategy_s2i_->KineticsConditionsMeshtyingSlaveSide())
  {
    if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") !=
        static_cast<int>(INPAR::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_s2i_->SetConditionSpecificScaTraParameters(*kinetics_slave_cond.second);
      // evaluate the condition
      ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatrastructures2i, "S2IKinetics", kinetics_slave_cond.first);
    }
  }

  // finalize scatra-structure matrix block
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      evaluate_matrix->Complete();

      auto evaluate_matrix_block =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(evaluate_matrix);
      auto slavematrix_block =
          CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(slavematrix);

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            CORE::ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        for (int iblock = 0; iblock < ScaTraField()->BlockMaps()->NumMaps(); ++iblock)
        {
          auto evaluate_iblock = evaluate_matrix_block->Matrix(iblock, 0);
          auto slave_iblock = slavematrix_block->Matrix(iblock, 0);

          auto scatra_slave_block_mapi =
              CORE::LINALG::IntersectMap(*ScaTraField()->BlockMaps()->Map(iblock),
                  *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

          CORE::LINALG::MatrixLogicalSplitAndTransform()(evaluate_iblock, *scatra_slave_block_mapi,
              *slave_map, 1.0, nullptr, &slave_slave_converter, slave_iblock, true, true);
        }
      }
      slavematrix->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      evaluate_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

      auto evaluate_matrix_sparse =
          CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(evaluate_matrix);
      auto slavematrix_sparse = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(slavematrix);

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            CORE::ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        CORE::LINALG::MatrixLogicalSplitAndTransform()(*evaluate_matrix_sparse,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter, *slavematrix_sparse, true, true);
      }
      slavematrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

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
SSI::ScatraStructureOffDiagCouplingSSTI::ScatraStructureOffDiagCouplingSSTI(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_scatra,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : ScatraStructureOffDiagCoupling(std::move(block_map_structure), std::move(full_map_structure),
          std::move(ssi_structure_meshtying), std::move(meshtying_strategy_s2i), std::move(scatra),
          std::move(structure)),
      full_map_scatra_(std::move(full_map_scatra))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCouplingSSTI::EvaluateOffDiagBlockStructureScatraDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> structurescatradomain) const
{
  ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockStructureScatraDomain(structurescatradomain);

  // finalize structure-scatra matrix block
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    case CORE::LINALG::MatrixType::block_condition_dof:
    {
      structurescatradomain->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
    {
      structurescatradomain->Complete(*full_map_scatra_, *FullMapStructure());
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
