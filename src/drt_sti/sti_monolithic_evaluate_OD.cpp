/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of OD blocks for monolithic STI
\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/

#include "sti_monolithic_evaluate_OD.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_adapter/adapter_coupling.H"

#include "../linalg/linalg_matrixtransform.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoODCoupling::ScatraThermoODCoupling(
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
    const Teuchos::RCP<const Epetra_Map> full_map_scatra,
    const Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : blockmaps_(blockmaps),
      blockmapthermo_(blockmapthermo),
      full_map_scatra_(full_map_scatra),
      full_map_thermo_(full_map_thermo),
      meshtying_strategy_s2i_(meshtying_strategy_s2i),
      meshtying_strategy_thermo_(meshtying_strategy_thermo),
      scatra_(scatra),
      thermo_(thermo)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCoupling::EvaluateODBlockScatraThermoDomain(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblock)
{
  // initialize scatra-thermo matrix block
  scatrathermoblock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", SCATRA::calc_scatra_mono_odblock_scatrathermo);

  // number of dofset associated with velocity-related dofs on scatra discretization
  eleparams.set<int>("ndsvel", 1);

  // remove state vectors from scatra discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scatra discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of scatra-thermo matrix block
  DRT::AssembleStrategy strategyscatrathermo(
      0,  // row assembly based on number of dofset associated with scatra dofs on scatra
          // discretization
      2,  // column assembly based on number of dofset associated with thermo dofs on scatra
          // discretization
      scatrathermoblock,  // scatra-thermo matrix block
      Teuchos::null,      // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble scatra-thermo matrix block
  scatra_->ScaTraField()->Discretization()->Evaluate(eleparams, strategyscatrathermo);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // finalize scatra-thermo block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      scatrathermoblock->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      scatrathermoblock->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCoupling::EvaluateODBlockThermoScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablock)
{
  // initialize thermo-scatra matrix block
  thermoscatrablock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", SCATRA::calc_scatra_mono_odblock_thermoscatra);

  // number of dofset associated with velocity-related dofs on thermo discretization
  eleparams.set<int>("ndsvel", 1);

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of thermo-scatra matrix block
  DRT::AssembleStrategy strategythermoscatra(
      0,  // row assembly based on number of dofset associated with thermo dofs on thermo
          // discretization
      2,  // column assembly based on number of dofset associated with scatra dofs on thermo
          // discretization
      thermoscatrablock,  // thermo-scatra matrix block
      Teuchos::null,      // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble thermo-scatra matrix block
  thermo_->ScaTraField()->Discretization()->Evaluate(eleparams, strategythermoscatra);

  // finalize thermo-scatra matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      thermoscatrablock->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      thermoscatrablock->Complete(*full_map_scatra_, *full_map_thermo_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoODCouplingMatchingNodes::ScatraThermoODCouplingMatchingNodes(
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
    const Teuchos::RCP<const Epetra_Map> full_map_scatra,
    const Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : ScatraThermoODCoupling(blockmaps, blockmapthermo, full_map_scatra, full_map_thermo,
          meshtying_strategy_s2i, meshtying_strategy_thermo, scatra, thermo)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMatchingNodes::EvaluateODBlockScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapthermo_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapthermo_, meshtying_strategy_s2i_->BlockMapsMaster(), 81, false, true));
      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
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

  // evaluate interface contibutions on slave side
  EvaluateScatraThermoInterfaceSlaveSide(slavematrix);

  // copy interface contributions from slave side to master side
  CopySlaveToMasterScatraThermoInterface(slavematrix, mastermatrix);

  // add contributions from slave side and master side
  scatrathermoblockinterface->Add(*slavematrix, false, 1.0, 0.0);
  scatrathermoblockinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize scatra-thermo matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMatchingNodes::EvaluateScatraThermoInterfaceSlaveSide(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix)
{
  // zero out slavematrtix
  slavematrix->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of auxiliary system matrix
  DRT::AssembleStrategy strategyscatrathermos2i(
      0,            // row assembly based on number of dofset associated with scatra dofs on scatra
                    // discretization
      2,            // column assembly based on number of dofset associated with thermo dofs on
                    // scatra discretization
      slavematrix,  // auxiliary system matrix
      Teuchos::null,  // no additional matrices of vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  std::vector<DRT::Condition*> conditions;
  scatra_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (auto& condition : conditions)
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_s2i_->SetConditionSpecificScaTraParameters(*condition);
      // evaluate the condition
      scatra_->ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatrathermos2i, "S2ICoupling", condition->GetInt("ConditionID"));
    }

  // finalize slave matrix
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMatchingNodes::CopySlaveToMasterScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<LINALG::SparseOperator>& mastermatrix)
{
  // zero out master matrix
  mastermatrix->Zero();

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      // cast master and slave matrix
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockslavematrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(slavematrix);
      Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra
      // fluxes w.r.t. slave-side thermo dofs
      LINALG::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into auxiliary system matrix
      for (int iblock = 0; iblock < meshtying_strategy_s2i_->BlockMapsSlave().NumMaps(); ++iblock)
        LINALG::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.0,
            ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()),
            mastermatrixsparse, true);

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap());

      // split auxiliary system matrix and assemble into scatra-thermo matrix block
      blockmastermatrix = mastermatrixsparse.Split<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapthermo_, scatra_->ScaTraField()->BlockMaps());

      // finalize master matrix
      mastermatrix->Complete();

      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      // cast master and slave matrix
      Teuchos::RCP<LINALG::SparseMatrix> sparseslavematrix =
          Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(slavematrix);
      Teuchos::RCP<LINALG::SparseMatrix> sparsemastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into scatra-thermo matrix block
      LINALG::MatrixRowTransform()(*sparseslavematrix, -1.0,
          ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()),
          *sparsemastermatrix, false);

      // finalize master matrix
      mastermatrix->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");

      break;
    }
  }
  // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed,
  // since these dofs will be condensed out later
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMatchingNodes::EvaluateODBlockThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // initialize slave and master matrix
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  meshtying_strategy_thermo_->MasterMatrix()->Zero();
  Teuchos::RCP<LINALG::SparseMatrix> mastermatrix = meshtying_strategy_thermo_->MasterMatrix();
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          meshtying_strategy_s2i_->BlockMapsSlave(), *blockmapthermo_, 81, false, true));
      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      meshtying_strategy_thermo_->SlaveMatrix()->Zero();
      slavematrix = meshtying_strategy_thermo_->SlaveMatrix();
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // create strategy for assembly of auxiliary system matrices
  DRT::AssembleStrategy strategythermoscatras2i(
      0,              // row assembly based on number of dofset associated with thermo dofs on
                      // thermo discretization
      2,              // column assembly based on number of dofset associated with scatra dofs on
                      // thermo discretization
      slavematrix,    // auxiliary system matrix for slave side
      mastermatrix,   // auxiliary system matrix for master side
      Teuchos::null,  // no additional matrices of vectors
      Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  std::vector<DRT::Condition*> conditions;
  thermo_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (auto& condition : conditions)
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo_->SetConditionSpecificScaTraParameters(*condition);
      // evaluate the condition
      thermo_->ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategythermoscatras2i, "S2ICoupling", condition->GetInt("ConditionID"));
    }

  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      // finalize auxiliary system matrices
      slavematrix->Complete();
      mastermatrix->Complete(*meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 0.0);

      // initialize temporary matrix
      LINALG::SparseMatrix ksm(
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), 27, false, true);

      // transform linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(), *mastermatrix,
          1.0, ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()), ksm,
          true, false);

      // finalize temporary matrix
      ksm.Complete(*meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(),
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      // split temporary matrix and assemble into thermo-scatra matrix block
      const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockksm(
          ksm.Split<LINALG::DefaultBlockMatrixStrategy>(
              scatra_->ScaTraField()->BlockMaps(), *blockmapthermo_));
      blockksm->Complete();
      thermoscatrablockinterface->Add(*blockksm, false, 1.0, 1.0);

      // finalize matrix
      thermoscatrablockinterface->Complete();

      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix->Complete(*meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());
      mastermatrix->Complete(*meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 0.0);

      // derive linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      // and assemble into thermo-scatra matrix block
      LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(), *mastermatrix,
          1.0, ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()),
          *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablockinterface), true, true);

      // finalize matrix
      thermoscatrablockinterface->Complete(*full_map_scatra_, *full_map_thermo_);

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoODCouplingMortarStandard::ScatraThermoODCouplingMortarStandard(
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmaps,
    Teuchos::RCP<LINALG::MultiMapExtractor> blockmapthermo,
    const Teuchos::RCP<const Epetra_Map> full_map_scatra,
    const Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : ScatraThermoODCoupling(blockmaps, blockmapthermo, full_map_scatra, full_map_thermo,
          meshtying_strategy_s2i, meshtying_strategy_thermo, scatra, thermo)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMortarStandard::EvaluateODBlockScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // initialize auxiliary system matrices for linearizations of slave-side and master-side
  // scatra fluxes w.r.t. slave-side thermo dofs
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  meshtying_strategy_s2i_->MasterMatrix()->Zero();
  Teuchos::RCP<LINALG::SparseMatrix> mastermatrix_sparse = meshtying_strategy_s2i_->MasterMatrix();
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapthermo_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix = meshtying_strategy_s2i_->SlaveMatrix();
      slavematrix->Zero();
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }


  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", INPAR::S2I::evaluate_condition_od);

  // create strategy for assembly of auxiliary system matrices
  SCATRA::MortarCellAssemblyStrategy strategyscatrathermos2i(slavematrix, INPAR::S2I::side_slave,
      INPAR::S2I::side_slave, Teuchos::null, INPAR::S2I::side_undefined, INPAR::S2I::side_undefined,
      mastermatrix_sparse, INPAR::S2I::side_master, INPAR::S2I::side_slave, Teuchos::null,
      INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
      INPAR::S2I::side_undefined, Teuchos::null, INPAR::S2I::side_undefined, 0, 1);

  // extract scatra-scatra interface coupling conditions
  std::vector<DRT::Condition*> conditions;
  scatra_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<DRT::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_s2i_->SetConditionSpecificScaTraParameters(*condition);
      // evaluate mortar integration cells
      meshtying_strategy_s2i_->EvaluateMortarCells(
          meshtying_strategy_s2i_->MortarDiscretization(condition->GetInt("ConditionID")),
          condparams, strategyscatrathermos2i);
    }
  }

  // finalize auxiliary system matrices
  mastermatrix_sparse->Complete(
      *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix->Complete();
      mastermatrix =
          meshtying_strategy_s2i_->MasterMatrix()->Split<LINALG::DefaultBlockMatrixStrategy>(
              *blockmapthermo_, meshtying_strategy_s2i_->BlockMapsMaster());
      mastermatrix->Complete();

      break;
    }
    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
      mastermatrix = mastermatrix_sparse;

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // assemble linearizations of slave-side and master-side scatra fluxes w.r.t. slave-side
  // thermo dofs into scatra-thermo matrix block
  scatrathermoblockinterface->Add(*slavematrix, false, 1.0, 0.0);
  scatrathermoblockinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed, since
  // these dofs will be condensed out later
  // finalize scatra-thermo matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(
          *thermo_->ScaTraField()->Discretization()->DofRowMap(), *full_map_scatra_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scatra discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoODCouplingMortarStandard::EvaluateODBlockThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes
  // w.r.t. slave-side and master-side scatra dofs
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          scatra_->ScaTraField()->BlockMaps(), *blockmapthermo_, 81, false, true));
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix = meshtying_strategy_thermo_->SlaveMatrix();
      slavematrix->Zero();
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", INPAR::S2I::evaluate_condition_od);

  // create strategy for assembly of auxiliary system matrix
  SCATRA::MortarCellAssemblyStrategy strategythermoscatras2i(slavematrix, INPAR::S2I::side_slave,
      INPAR::S2I::side_slave, slavematrix, INPAR::S2I::side_slave, INPAR::S2I::side_master,
      Teuchos::null, INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
      INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
      INPAR::S2I::side_undefined, Teuchos::null, INPAR::S2I::side_undefined, 0, 1);

  // extract scatra-scatra interface coupling conditions
  std::vector<DRT::Condition*> conditions;
  thermo_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<DRT::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo_->SetConditionSpecificScaTraParameters(*condition);
      // evaluate mortar integration cells
      meshtying_strategy_thermo_->EvaluateMortarCells(
          meshtying_strategy_thermo_->MortarDiscretization(condition->GetInt("ConditionID")),
          condparams, strategythermoscatras2i);
    }
  }

  // finalize auxiliary system matrix
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      slavematrix->Complete(*full_map_scatra_, *full_map_thermo_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side and master-side
  // scatra dofs into thermo-scatra matrix block
  thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 0.0);

  // linearizations of master-side thermo fluxes w.r.t. scatra dofs are not needed, since
  // thermo fluxes are source terms and thus only evaluated once on slave side

  // finalize thermo-scatra matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      thermoscatrablockinterface->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      thermoscatrablockinterface->Complete(*full_map_scatra_, *full_map_thermo_);
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->Discretization()->ClearState();

  return;
}
