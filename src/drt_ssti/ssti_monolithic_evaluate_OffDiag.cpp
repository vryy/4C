/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#include "ssti_monolithic_evaluate_OffDiag.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_ssi/ssi_utils.H"

#include "../drt_structure_new/str_enum_lists.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSTI::ThermoStructureOffDiagCoupling::ThermoStructureOffDiagCoupling(
    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapstructure,
    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapthermo,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const SSI::UTILS::SSIStructureMeshTying> ssti_structure_meshtying,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : blockmapstructure_(std::move(blockmapstructure)),
      blockmapthermo_(std::move(blockmapthermo)),
      full_map_structure_(std::move(full_map_structure)),
      full_map_thermo_(std::move(full_map_thermo)),
      meshtying_strategy_thermo_(std::move(meshtying_strategy_thermo)),
      ssti_structure_meshtying_(std::move(ssti_structure_meshtying)),
      structure_(std::move(structure)),
      thermo_(std::move(thermo))
{
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::EvaluateOffDiagBlockThermoStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator> thermostructuredomain)
{
  // initialize thermo-structure matrix block
  thermostructuredomain->Zero();

  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", SCATRA::calc_scatra_mono_odblock_mesh);

  // number of dofset associated with displacement-related dofs on thermo discretization
  eleparams.set<int>("ndsdisp", 1);

  // number of dofset associated with velocity-related dofs on thermo discretization
  eleparams.set<int>("ndsvel", 1);

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of thermo-structure matrix block
  DRT::AssembleStrategy strategyscatrastructure(
      0,  // row assembly based on number of dofset associated with thermo dofs on thermo
          // discretization
      1,  // column assembly based on number of dofset associated with structural dofs on thermo
          // discretization
      thermostructuredomain,  // thermo-structure matrix block
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  thermo_->ScaTraField()->Discretization()->Evaluate(eleparams, strategyscatrastructure);

  // finalize thermo-structure matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      thermostructuredomain->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      thermostructuredomain->Complete(*full_map_structure_, *full_map_thermo_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  thermo_->ScaTraField()->Discretization()->ClearState();
}
/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::EvaluateOffDiagBlockThermoStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> thermostructureinterface)
{
  thermostructureinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsMaster(), 81, false, true));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *meshtying_strategy_thermo_->CouplingAdapter()->MasterDofMap(), 27, false, true));
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  EvaluateThermoStructureInterfaceSlaveSide(slavematrix);

  CopySlaveToMasterThermoStructureInterface(slavematrix, mastermatrix);

  thermostructureinterface->Add(*slavematrix, false, 1.0, 1.0);
  thermostructureinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize thermo-structure matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      thermostructureinterface->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      thermostructureinterface->Complete(*full_map_structure_, *full_map_thermo_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  thermo_->ScaTraField()->Discretization()->ClearState();
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::EvaluateOffDiagBlockStructureThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator> structurethermodomain)
{
  structurethermodomain->Zero();

  Teuchos::ParameterList eleparams;

  eleparams.set("action", "calc_struct_stiffscalar");

  eleparams.set<int>("differentiationtype", static_cast<int>(STR::DifferentiationType::temp));

  eleparams.set<double>("total time", structure_->Time());

  structure_->Discretization()->ClearState();

  structure_->Discretization()->SetState("displacement", structure_->Dispnp());

  // create strategy for assembly of structure-thermo matrix block
  DRT::AssembleStrategy strategystructurescatra(
      0,  // row assembly based on number of dofset associated with structure dofs on structural
          // discretization
      2,  // column assembly based on number of dofset associated with thermo dofs on structural
          // discretization
      structurethermodomain,  // structure-thermo matrix block
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  structure_->Discretization()->Evaluate(eleparams, strategystructurescatra);

  // need to scale structurethermoblock_ with 'timefac' to getcorrect implementation
  structurethermodomain->Scale(1.0 - structure_->TimIntParam());

  // finalize structure-thermo matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      structurethermodomain->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      structurethermodomain->Complete(*full_map_thermo_, *full_map_structure_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  structure_->Discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::CopySlaveToMasterThermoStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<LINALG::SparseOperator>& mastermatrix)
{
  mastermatrix->Zero();

  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      const int numberthermoblocks = thermo_->ScaTraField()->BlockMaps().NumMaps();

      auto blockslavematrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      LINALG::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_thermo_->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
      {
        for (const auto& meshtying : ssti_structure_meshtying_->MeshtyingHandlers())
        {
          auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
          auto slave_side_converter_struct = meshtying->SlaveSideConverter();

          auto slave_side_converter_thermo =
              ADAPTER::CouplingSlaveConverter(*meshtying_strategy_thermo_->CouplingAdapter());

          LINALG::MatrixLogicalSplitAndTransform()(blockslavematrix->Matrix(iblock, 0),
              *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
              &slave_side_converter_thermo, &(*slave_side_converter_struct), mastermatrixsparse,
              true, true);
        }
      }

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*full_map_structure_, *full_map_thermo_);

      // split sparse matrix to block matrix
      auto mastermatrix_split = mastermatrixsparse.Split<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, *blockmapthermo_);
      mastermatrix_split->Complete();
      blockmastermatrix->Add(*mastermatrix_split, false, 1.0, 1.0);

      mastermatrix->Complete();

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      auto sparseslavematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (const auto& meshtying : ssti_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto slave_side_converter_struct = meshtying->SlaveSideConverter();
        auto slave_side_converter_thermo =
            ADAPTER::CouplingSlaveConverter(*meshtying_strategy_thermo_->CouplingAdapter());

        LINALG::MatrixLogicalSplitAndTransform()(*sparseslavematrix,
            *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
            &slave_side_converter_thermo, &(*slave_side_converter_struct), *sparsemastermatrix,
            true, true);
      }

      mastermatrix->Complete(
          *full_map_structure_, *meshtying_strategy_thermo_->CouplingAdapter()->MasterDofMap());
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
void SSTI::ThermoStructureOffDiagCoupling::EvaluateThermoStructureInterfaceSlaveSide(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix)
{
  Teuchos::ParameterList condparams;

  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // number of dofset associated with displacement-related dofs on scalar transport discretization
  condparams.set<int>("ndsdisp", 1);

  condparams.set<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::disp));

  thermo_->ScaTraField()->Discretization()->ClearState();

  thermo_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  Teuchos::RCP<LINALG::SparseOperator> evaluate_matrix;
  if (thermo_->ScaTraField()->MatrixType() == LINALG::MatrixType::sparse)
  {
    evaluate_matrix = Teuchos::rcp(new LINALG::SparseMatrix(
        *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
  }
  else
  {
    evaluate_matrix =
        Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
            *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsSlave(), 81, false, true));
  }

  // create strategy for assembly of auxiliary system matrix
  DRT::AssembleStrategy strategyscatrastructures2i(
      0,  // row assembly based on number of dofset associated with thermo dofs on
          // thermo discretization
      1,  // column assembly based on number of dofset associated with structural dofs on
          // thermo discretization
      evaluate_matrix, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate interface coupling
  for (auto kinetics_slave_cond :
      meshtying_strategy_thermo_->KineticsConditionsMeshtyingSlaveSide())
  {
    if (kinetics_slave_cond.second->GetInt("kinetic model") !=
        static_cast<int>(INPAR::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo_->SetConditionSpecificScaTraParameters(*kinetics_slave_cond.second);
      // evaluate the condition
      thermo_->ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatrastructures2i, "S2IKinetics", kinetics_slave_cond.first);
    }
  }

  // finalize thermo-structure matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      evaluate_matrix->Complete();

      auto evaluate_matrix_block =
          LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(evaluate_matrix);
      auto slavematrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(slavematrix);

      // "slave side" from thermo and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssti_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter = ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        for (int iblock = 0; iblock < thermo_->ScaTraField()->BlockMaps().NumMaps(); ++iblock)
        {
          auto evaluate_iblock = evaluate_matrix_block->Matrix(iblock, 0);
          auto slave_iblock = slavematrix_block->Matrix(iblock, 0);

          auto scatra_slave_block_mapi =
              LINALG::IntersectMap(*thermo_->ScaTraField()->BlockMaps().Map(iblock),
                  *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

          LINALG::MatrixLogicalSplitAndTransform()(evaluate_iblock, *scatra_slave_block_mapi,
              *slave_map, 1.0, nullptr, &slave_slave_converter, slave_iblock, true, true);
        }
      }
      slavematrix->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      auto evaluate_matrix_sparse = LINALG::CastToConstSparseMatrixAndCheckSuccess(evaluate_matrix);
      auto slavematrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(slavematrix);

      evaluate_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      // "slave side" from thermo and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssti_structure_meshtying_->MeshtyingHandlers())
      {
        auto slave_slave_transformation = meshtying->SlaveSlaveTransformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter = ADAPTER::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        LINALG::MatrixLogicalSplitAndTransform()(*evaluate_matrix_sparse,
            *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter, *slavematrix_sparse, true, true);
      }
      slavematrix->Complete(
          *full_map_structure_, *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}
