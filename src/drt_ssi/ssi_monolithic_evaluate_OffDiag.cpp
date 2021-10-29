/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSI
\level 2


 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_evaluate_OffDiag.H"

#include "ssi_utils.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_structure_new/str_enum_lists.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraStructureOffDiagCoupling::ScatraStructureOffDiagCoupling(
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
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
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
    : ScatraStructureOffDiagCoupling(std::move(block_map_structure), std::move(full_map_structure),
          std::move(ssi_structure_meshtying), std::move(meshtying_strategy_s2i), std::move(scatra),
          std::move(structure)),
      scatra_manifold_(std::move(scatra_manifold))
{
}


/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureblock)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", SCATRA::calc_scatra_mono_odblock_mesh);

  // number of dofset associated with displacement-related dofs on scalar transport
  // discretization
  eleparams.set<int>("ndsdisp", 1);

  // number of dofset associated with velocity-related dofs on scalar transport discretization
  eleparams.set<int>("ndsvel", 1);

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();

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

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraManifoldStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructureblock)
{
  dserror("not implemented");
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraManifoldStructureOffDiagCoupling::EvaluateOffDiagBlockScatraManifoldStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> scatramanifoldstructureblock)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", SCATRA::calc_scatra_mono_odblock_mesh);

  // number of dofset associated with displacement-related dofs on scalar transport
  // discretization
  eleparams.set<int>("ndsdisp", 1);

  // number of dofset associated with velocity-related dofs on scalar transport discretization
  eleparams.set<int>("ndsvel", 1);

  // remove state vectors from scalar transport discretization
  scatra_manifold_->Discretization()->ClearState();

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

  // remove state vectors from scalar transport discretization
  scatra_manifold_->Discretization()->ClearState();
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> scatrastructureinterface)
{
  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *block_map_structure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *block_map_structure_, meshtying_strategy_s2i_->BlockMapsMaster(), 81, false, true));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true));
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // evaluate interface contributions on slave side
  EvaluateScatraStructureInterfaceSlaveSide(slavematrix);

  // copy interface contributions from slave side to master side
  CopySlaveToMasterScatraStructureInterface(slavematrix, mastermatrix);

  // add contributions from slave side and master side
  scatrastructureinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrastructureinterface->Add(*mastermatrix, false, 1.0, 1.0);
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) const
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
void SSI::ScatraStructureOffDiagCoupling::CopySlaveToMasterScatraStructureInterface(
    Teuchos::RCP<const LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<LINALG::SparseOperator>& mastermatrix)
{
  mastermatrix->Zero();
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      const int numberscatrablocks = ScaTraField()->BlockMaps().NumMaps();

      // cast scatrastructureinterfaceslaveside
      // cast master and slave matrix
      auto blockslavematrix =
          Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      LINALG::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
      {
        for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
        {
          auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
          auto slave_side_converter_struct = meshtying->SlaveSideConverter();

          auto slave_side_converter_scatra =
              ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

          LINALG::MatrixLogicalSplitAndTransform()(blockslavematrix->Matrix(iblock, 0),
              *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
              &slave_side_converter_scatra, &(*slave_side_converter_struct), mastermatrixsparse,
              true, true);
        }
      }

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*FullMapStructure(), *ScaTraField()->DofRowMap());

      // split auxiliary system matrix and assemble into scatra-structure matrix block
      auto mastermatrix_split = mastermatrixsparse.Split<LINALG::DefaultBlockMatrixStrategy>(
          *block_map_structure_, ScaTraField()->BlockMaps());
      mastermatrix_split->Complete();
      blockmastermatrix->Add(*mastermatrix_split, false, 1.0, 1.0);

      mastermatrix->Complete();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // cast scatrastructureinterfaceslaveside
      auto sparseslavematrix = Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mastermatrix);

      // copy slave side values to master side and scale with minus 1. Insert into
      // scatrastructureinterface_sparse
      for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto slave_side_converter_struct = meshtying->SlaveSideConverter();
        auto slave_side_converter_scatra =
            ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter());

        LINALG::MatrixLogicalSplitAndTransform()(*sparseslavematrix,
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
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateScatraStructureInterfaceSlaveSide(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix)
{
  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // linearization of boundary flux w.r.t. displacement
  condparams.set<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::disp));

  // number of dofset associated with displacement-related dofs on scalar transport discretization
  condparams.set<int>("ndsdisp", 1);

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

  Teuchos::RCP<LINALG::SparseOperator> evaluate_matrix;
  if (ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
  {
    evaluate_matrix = Teuchos::rcp(new LINALG::SparseMatrix(
        *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
  }
  else
  {
    evaluate_matrix =
        Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
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
    if (kinetics_slave_cond.second->GetInt("kinetic model") !=
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
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      evaluate_matrix->Complete();

      auto evaluate_matrix_block =
          LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(evaluate_matrix);
      auto slavematrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(slavematrix);

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter = ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        for (int iblock = 0; iblock < ScaTraField()->BlockMaps().NumMaps(); ++iblock)
        {
          auto evaluate_iblock = evaluate_matrix_block->Matrix(iblock, 0);
          auto slave_iblock = slavematrix_block->Matrix(iblock, 0);

          auto scatra_slave_block_mapi =
              LINALG::IntersectMap(*ScaTraField()->BlockMaps().Map(iblock),
                  *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

          LINALG::MatrixLogicalSplitAndTransform()(evaluate_iblock, *scatra_slave_block_mapi,
              *slave_map, 1.0, nullptr, &slave_slave_converter, slave_iblock, true, true);
        }
      }
      slavematrix->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      evaluate_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

      auto evaluate_matrix_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(evaluate_matrix);
      auto slavematrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(slavematrix);

      // "slave side" from scatra and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssi_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter = ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        LINALG::MatrixLogicalSplitAndTransform()(*evaluate_matrix_sparse,
            *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter, *slavematrix_sparse, true, true);
      }
      slavematrix->Complete(
          *full_map_structure_, *meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap());

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraStructureOffDiagCouplingSSTI::ScatraStructureOffDiagCouplingSSTI(
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_scatra,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssi_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
    : ScatraStructureOffDiagCoupling(std::move(block_map_structure), std::move(full_map_structure),
          std::move(ssi_structure_meshtying), std::move(meshtying_strategy_s2i), std::move(scatra),
          std::move(structure)),
      full_map_scatra_(std::move(full_map_scatra))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCouplingSSTI::EvaluateOffDiagBlockStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator> structurescatradomain) const
{
  ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockStructureScatraDomain(structurescatradomain);

  // finalize structure-scatra matrix block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      structurescatradomain->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      structurescatradomain->Complete(*full_map_scatra_, *FullMapStructure());
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}
