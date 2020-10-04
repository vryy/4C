/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic STI
\level 2


 */
/*----------------------------------------------------------------------*/

#include "sti_monolithic_evaluate_OffDiag.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
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
STI::ScatraThermoOffDiagCoupling::ScatraThermoOffDiagCoupling(
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : block_map_thermo_(std::move(block_map_thermo)),
      block_map_thermo_interface_(std::move(block_map_thermo_interface)),
      full_map_scatra_(std::move(full_map_scatra)),
      full_map_thermo_(std::move(full_map_thermo)),
      interface_map_scatra_(std::move(interface_map_scatra)),
      interface_map_thermo_(std::move(interface_map_thermo)),
      isale_(std::move(isale)),
      meshtying_strategy_scatra_(std::move(meshtying_strategy_scatra)),
      meshtying_strategy_thermo_(std::move(meshtying_strategy_thermo)),
      scatra_(std::move(scatra)),
      thermo_(std::move(thermo))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCoupling::EvaluateOffDiagBlockScatraThermoDomain(
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
  ScaTraField()->Discretization()->ClearState();

  // add state vectors to scatra discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

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
  ScaTraField()->Discretization()->Evaluate(eleparams, strategyscatrathermo);

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();

  // finalize scatra-thermo block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      scatrathermoblock->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatrathermoblock->Complete(*FullMapThermo(), *FullMapScatra());
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
void STI::ScatraThermoOffDiagCoupling::EvaluateOffDiagBlockThermoScatraDomain(
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

  // in case of deforming mesh: set displacement
  if (IsAle()) eleparams.set<int>("ndsdisp", 1);

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  ThermoField()->AddTimeIntegrationSpecificVectors();

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
  ThermoField()->Discretization()->Evaluate(eleparams, strategythermoscatra);

  // finalize thermo-scatra matrix block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      thermoscatrablock->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      thermoscatrablock->Complete(*FullMapScatra(), *FullMapThermo());
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMatchingNodes::ScatraThermoOffDiagCouplingMatchingNodes(
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface_slave,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : ScatraThermoOffDiagCoupling(std::move(block_map_thermo),
          std::move(block_map_thermo_interface), std::move(full_map_scatra),
          std::move(full_map_thermo), std::move(interface_map_scatra),
          std::move(interface_map_thermo), std::move(isale), std::move(meshtying_strategy_scatra),
          std::move(meshtying_strategy_thermo), std::move(scatra), std::move(thermo)),
      block_map_thermo_interface_slave_(block_map_thermo_interface_slave)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateOffDiagBlockScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsSlave(), 81, false,
          true));
      mastermatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsMaster(), 81, false,
          true));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::SparseMatrix(
          *MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap(), 27, false, true));
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
  scatrathermoblockinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrathermoblockinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize scatra-thermo matrix block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(*InterfaceMapThermo(), *InterfaceMapScaTra());
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

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateScatraThermoInterfaceSlaveSide(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix)
{
  // zero out slavematrtix
  slavematrix->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // in case of deforming mesh: set displacement
  if (IsAle()) condparams.set<int>("ndsdisp", 1);

  // set type of differentiation to temperature
  condparams.set<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::temp));

  // remove state vectors from scalar transport discretization
  ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  ScaTraField()->AddTimeIntegrationSpecificVectors();

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
  ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (const auto& condition : conditions)
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyScaTra()->SetConditionSpecificScaTraParameters(*condition);
      // evaluate the condition
      ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatrathermos2i, "S2ICoupling", condition->GetInt("ConditionID"));
    }

  // finalize slave matrix
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix->Complete(
          *InterfaceMapThermo(), *MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap());
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
void STI::ScatraThermoOffDiagCouplingMatchingNodes::CopySlaveToMasterScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<LINALG::SparseOperator>& mastermatrix)
{
  // zero out master matrix
  mastermatrix->Zero();

  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      // cast master and slave matrix
      const auto blockslavematrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra
      // fluxes w.r.t. slave-side thermo dofs
      LINALG::SparseMatrix mastermatrixsparse(
          *MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into auxiliary system matrix
      for (int iblock = 0; iblock < MeshtyingStrategyScaTra()->BlockMapsSlave().NumMaps(); ++iblock)
        LINALG::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.0,
            ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
            mastermatrixsparse, true);

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap());

      // split auxiliary system matrix and assemble into scatra-thermo matrix block
      blockmastermatrix = mastermatrixsparse.Split<LINALG::DefaultBlockMatrixStrategy>(
          *BlockMapThermo(), ScaTraField()->BlockMaps());

      // finalize master matrix
      mastermatrix->Complete();

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      // cast master and slave matrix
      const auto sparseslavematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into scatra-thermo matrix block
      LINALG::MatrixRowTransform()(*sparseslavematrix, -1.0,
          ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
          *sparsemastermatrix, false);

      // finalize master matrix
      mastermatrix->Complete(
          *InterfaceMapThermo(), *MeshtyingStrategyScaTra()->InterfaceMaps()->Map(2));

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
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateOffDiagBlockThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // initialize slave and master matrix
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  MeshtyingStrategyThermo()->MasterMatrix()->Zero();
  auto mastermatrix = MeshtyingStrategyThermo()->MasterMatrix();
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          MeshtyingStrategyScaTra()->BlockMapsSlave(), *BlockMapThermoInterfaceSlave(), 81, false,
          true));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      MeshtyingStrategyThermo()->SlaveMatrix()->Zero();
      slavematrix = MeshtyingStrategyThermo()->SlaveMatrix();
      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  ThermoField()->AddTimeIntegrationSpecificVectors();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_od);

  // set differentiation type to elch
  condparams.set<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::elch));

  // in case of deforming mesh: set displacement
  if (IsAle()) condparams.set<int>("ndsdisp", 1);

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
  ThermoField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (const auto& condition : conditions)
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyThermo()->SetConditionSpecificScaTraParameters(*condition);
      // evaluate the condition
      ThermoField()->Discretization()->EvaluateCondition(
          condparams, strategythermoscatras2i, "S2ICoupling", condition->GetInt("ConditionID"));
    }

  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      // finalize auxiliary system matrices
      slavematrix->Complete();
      mastermatrix->Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

      // initialize temporary matrix
      LINALG::SparseMatrix ksm(
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap(), 27, false, true);

      // transform linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(), *mastermatrix,
          1.0, ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()), ksm,
          true, false);

      // finalize temporary matrix
      ksm.Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      // split temporary matrix and assemble into thermo-scatra matrix block
      const auto blockksm = ksm.Split<LINALG::DefaultBlockMatrixStrategy>(
          MeshtyingStrategyScaTra()->BlockMapsMaster(), *BlockMapThermoInterfaceSlave());
      blockksm->Complete();
      thermoscatrablockinterface->Add(*blockksm, false, 1.0, 1.0);

      // finalize matrix
      thermoscatrablockinterface->Complete();

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix->Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());
      mastermatrix->Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

      // derive linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      // and assemble into thermo-scatra matrix block
      LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(), *mastermatrix,
          1.0, ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
          *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(thermoscatrablockinterface), true, true);

      // finalize matrix
      thermoscatrablockinterface->Complete(
          *InterfaceMapScaTra(), *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      break;
    }
    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMortarStandard::ScatraThermoOffDiagCouplingMortarStandard(
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
    : ScatraThermoOffDiagCoupling(std::move(block_map_thermo),
          std::move(block_map_thermo_interface), std::move(full_map_scatra),
          std::move(full_map_thermo), std::move(interface_map_scatra),
          std::move(interface_map_thermo), std::move(isale), std::move(meshtying_strategy_scatra),
          std::move(meshtying_strategy_thermo), std::move(scatra), std::move(thermo))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::EvaluateOffDiagBlockScatraThermoInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // initialize auxiliary system matrices for linearizations of slave-side and master-side
  // scatra fluxes w.r.t. slave-side thermo dofs
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  MeshtyingStrategyScaTra()->MasterMatrix()->Zero();
  Teuchos::RCP<LINALG::SparseMatrix> mastermatrix_sparse =
      MeshtyingStrategyScaTra()->MasterMatrix();
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsSlave(), 81, false,
          true));
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      slavematrix = MeshtyingStrategyScaTra()->SlaveMatrix();
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
  ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<DRT::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyScaTra()->SetConditionSpecificScaTraParameters(*condition);
      // evaluate mortar integration cells
      MeshtyingStrategyScaTra()->EvaluateMortarCells(
          MeshtyingStrategyScaTra()->MortarDiscretization(condition->GetInt("ConditionID")),
          condparams, strategyscatrathermos2i);
    }
  }

  // finalize auxiliary system matrices
  mastermatrix_sparse->Complete(
      *InterfaceMapThermo(), *MeshtyingStrategyScaTra()->InterfaceMaps()->Map(2));

  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      mastermatrix =
          MeshtyingStrategyScaTra()->MasterMatrix()->Split<LINALG::DefaultBlockMatrixStrategy>(
              *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsMaster());
      mastermatrix->Complete();

      break;
    }
    case LINALG::MatrixType::sparse:
    {
      slavematrix->Complete(
          *InterfaceMapThermo(), *MeshtyingStrategyScaTra()->InterfaceMaps()->Map(1));
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
  scatrathermoblockinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrathermoblockinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed, since
  // these dofs will be condensed out later
  // finalize scatra-thermo matrix block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(*InterfaceMapThermo(), *InterfaceMapScaTra());
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scatra discretization
  ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::EvaluateOffDiagBlockThermoScatraInterface(
    Teuchos::RCP<LINALG::SparseOperator>& thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  ThermoField()->AddTimeIntegrationSpecificVectors();

  // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes
  // w.r.t. slave-side and master-side scatra dofs
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          ScaTraField()->BlockMaps(), *BlockMapThermoInterface(), 81, false, true));
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      slavematrix = MeshtyingStrategyThermo()->SlaveMatrix();
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
  ThermoField()->Discretization()->GetCondition("S2ICoupling", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<DRT::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyThermo()->SetConditionSpecificScaTraParameters(*condition);
      // evaluate mortar integration cells
      MeshtyingStrategyThermo()->EvaluateMortarCells(
          MeshtyingStrategyThermo()->MortarDiscretization(condition->GetInt("ConditionID")),
          condparams, strategythermoscatras2i);
    }
  }

  // finalize auxiliary system matrix
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      slavematrix->Complete(
          *InterfaceMapScaTra(), *MeshtyingStrategyThermo()->InterfaceMaps()->Map(1));
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
  thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

  // linearizations of master-side thermo fluxes w.r.t. scatra dofs are not needed, since
  // thermo fluxes are source terms and thus only evaluated once on slave side

  // finalize thermo-scatra matrix block
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    {
      thermoscatrablockinterface->Complete();
      break;
    }

    case LINALG::MatrixType::sparse:
    {
      thermoscatrablockinterface->Complete(*InterfaceMapScaTra(), *InterfaceMapThermo());
      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> STI::BuildScatraThermoOffDiagCoupling(
    const INPAR::S2I::CouplingType& couplingtype,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const LINALG::MultiMapExtractor> block_map_thermo_interface_slave,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo)
{
  Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> scatrathermooffdiagcoupling = Teuchos::null;

  switch (couplingtype)
  {
    case INPAR::S2I::coupling_matching_nodes:
    {
      scatrathermooffdiagcoupling = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMatchingNodes(
          block_map_thermo, block_map_thermo_interface, block_map_thermo_interface_slave,
          full_map_scatra, full_map_thermo, interface_map_scatra, interface_map_thermo, isale,
          meshtying_strategy_scatra, meshtying_strategy_thermo, scatra, thermo));
      break;
    }
    case INPAR::S2I::coupling_mortar_standard:
    {
      scatrathermooffdiagcoupling = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMortarStandard(
          block_map_thermo, block_map_thermo_interface, full_map_scatra, full_map_thermo,
          interface_map_scatra, interface_map_thermo, isale, meshtying_strategy_scatra,
          meshtying_strategy_thermo, scatra, thermo));
      break;
    }
    default:
    {
      dserror("Not supported coupling type");
      break;
    }
  }

  return scatrathermooffdiagcoupling;
}
