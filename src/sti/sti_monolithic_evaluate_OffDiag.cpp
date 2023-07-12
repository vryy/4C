/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic STI
\level 2


 */
/*----------------------------------------------------------------------*/

#include <utility>

#include "sti_monolithic_evaluate_OffDiag.H"

#include "scatra_ele_action.H"
#include "scatra_timint_implicit.H"

#include "linalg_mapextractor.H"
#include "linalg_sparseoperator.H"

#include "lib_assemblestrategy.H"
#include "lib_discret.H"
#include "lib_utils_parameter_list.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

#include "linalg_utils_sparse_algebra_create.H"

#include "coupling_adapter.H"
#include "coupling_adapter_converter.H"

#include "linalg_matrixtransform.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCoupling::ScatraThermoOffDiagCoupling(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface,
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
      isale_(isale),
      meshtying_strategy_scatra_(std::move(meshtying_strategy_scatra)),
      meshtying_strategy_thermo_(std::move(meshtying_strategy_thermo)),
      scatra_(std::move(scatra)),
      thermo_(std::move(thermo))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCoupling::EvaluateOffDiagBlockScatraThermoDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrathermoblock)
{
  // initialize scatra-thermo matrix block
  scatrathermoblock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_scatra_mono_odblock_scatrathermo, eleparams);

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
    case CORE::LINALG::MatrixType::block_condition:
    {
      scatrathermoblock->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCoupling::EvaluateOffDiagBlockThermoScatraDomain(
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermoscatrablock)
{
  // initialize thermo-scatra matrix block
  thermoscatrablock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_scatra_mono_odblock_thermoscatra, eleparams);

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
    case CORE::LINALG::MatrixType::block_condition:
    {
      thermoscatrablock->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMatchingNodes::ScatraThermoOffDiagCouplingMatchingNodes(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface_slave,
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
          std::move(interface_map_thermo), isale, std::move(meshtying_strategy_scatra),
          std::move(meshtying_strategy_thermo), std::move(scatra), std::move(thermo)),
      block_map_thermo_interface_slave_(std::move(block_map_thermo_interface_slave))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateOffDiagBlockScatraThermoInterface(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<CORE::LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsSlave(), 81, false,
              true));
      mastermatrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsMaster(), 81, false,
              true));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
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
    case CORE::LINALG::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateScatraThermoInterfaceSlaveSide(
    Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix)
{
  // zero out slavematrtix
  slavematrix->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_od, condparams);

  // set type of differentiation to temperature
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
      "differentiationtype", SCATRA::DifferentiationType::temp, condparams);

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

  // evaluate scatra-scatra interface kinetics
  std::vector<DRT::Condition*> conditions;
  for (const auto& kinetics_slave_cond :
      MeshtyingStrategyScaTra()->KineticsConditionsMeshtyingSlaveSide())
  {
    if (kinetics_slave_cond.second->GetInt("kinetic model") !=
        static_cast<int>(INPAR::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyScaTra()->SetConditionSpecificScaTraParameters(*kinetics_slave_cond.second);
      // evaluate the condition
      ScaTraField()->Discretization()->EvaluateCondition(condparams, strategyscatrathermos2i,
          "S2IKinetics", kinetics_slave_cond.second->GetInt("ConditionID"));
    }
  }

  // finalize slave matrix
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::CopySlaveToMasterScatraThermoInterface(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<CORE::LINALG::SparseOperator>& mastermatrix)
{
  // zero out master matrix
  mastermatrix->Zero();

  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      // cast master and slave matrix
      const auto blockslavematrix =
          Teuchos::rcp_dynamic_cast<const CORE::LINALG::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra
      // fluxes w.r.t. slave-side thermo dofs
      CORE::LINALG::SparseMatrix mastermatrixsparse(
          *MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into auxiliary system matrix
      for (int iblock = 0; iblock < MeshtyingStrategyScaTra()->BlockMapsSlave().NumMaps(); ++iblock)
      {
        CORE::LINALG::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.0,
            CORE::ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
            mastermatrixsparse, true);
      }

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap());

      // split auxiliary system matrix and assemble into scatra-thermo matrix block
      blockmastermatrix = mastermatrixsparse.Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *BlockMapThermo(), *ScaTraField()->BlockMaps());

      // finalize master matrix
      mastermatrix->Complete();

      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      // cast master and slave matrix
      const auto sparseslavematrix =
          Teuchos::rcp_dynamic_cast<const CORE::LINALG::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into scatra-thermo matrix block
      CORE::LINALG::MatrixRowTransform()(*sparseslavematrix, -1.0,
          CORE::ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::EvaluateOffDiagBlockThermoScatraInterface(
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // initialize slave and master matrix
  Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix(Teuchos::null);
  MeshtyingStrategyThermo()->MasterMatrix()->Zero();
  auto mastermatrix = MeshtyingStrategyThermo()->MasterMatrix();
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              MeshtyingStrategyScaTra()->BlockMapsSlave(), *BlockMapThermoInterfaceSlave(), 81,
              false, true));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
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
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
      "action", SCATRA::BoundaryAction::calc_s2icoupling_od, condparams);

  // set differentiation type to elch
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::DifferentiationType>(
      "differentiationtype", SCATRA::DifferentiationType::elch, condparams);

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

  // evaluate scatra-scatra interface kinetics
  for (const auto& kinetics_slave_cond :
      MeshtyingStrategyThermo()->KineticsConditionsMeshtyingSlaveSide())
  {
    if (kinetics_slave_cond.second->GetInt("kinetic model") !=
        static_cast<int>(INPAR::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      MeshtyingStrategyThermo()->SetConditionSpecificScaTraParameters(*kinetics_slave_cond.second);
      // evaluate the condition
      ThermoField()->Discretization()->EvaluateCondition(condparams, strategythermoscatras2i,
          "S2IKinetics", kinetics_slave_cond.second->GetInt("ConditionID"));
    }
  }

  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      // finalize auxiliary system matrices
      slavematrix->Complete();
      mastermatrix->Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->SlaveDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

      // initialize temporary matrix
      CORE::LINALG::SparseMatrix ksm(
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap(), 27, false, true);

      // transform linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      CORE::LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(),
          *mastermatrix, 1.0,
          CORE::ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()), ksm,
          true, false);

      // finalize temporary matrix
      ksm.Complete(*MeshtyingStrategyScaTra()->CouplingAdapter()->MasterDofMap(),
          *MeshtyingStrategyThermo()->CouplingAdapter()->SlaveDofMap());

      // split temporary matrix and assemble into thermo-scatra matrix block
      const auto blockksm = ksm.Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
          MeshtyingStrategyScaTra()->BlockMapsMaster(), *BlockMapThermoInterfaceSlave());
      blockksm->Complete();
      thermoscatrablockinterface->Add(*blockksm, false, 1.0, 1.0);

      // finalize matrix
      thermoscatrablockinterface->Complete();

      break;
    }
    case CORE::LINALG::MatrixType::sparse:
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
      CORE::LINALG::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(),
          *mastermatrix, 1.0,
          CORE::ADAPTER::CouplingSlaveConverter(*MeshtyingStrategyScaTra()->CouplingAdapter()),
          *Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(thermoscatrablockinterface), true,
          true);

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMortarStandard::ScatraThermoOffDiagCouplingMortarStandard(
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface,
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
          std::move(interface_map_thermo), isale, std::move(meshtying_strategy_scatra),
          std::move(meshtying_strategy_thermo), std::move(scatra), std::move(thermo))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::EvaluateOffDiagBlockScatraThermoInterface(
    Teuchos::RCP<CORE::LINALG::SparseOperator> scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // initialize auxiliary system matrices for linearizations of slave-side and master-side
  // scatra fluxes w.r.t. slave-side thermo dofs
  Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix(Teuchos::null);
  MeshtyingStrategyScaTra()->MasterMatrix()->Zero();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mastermatrix_sparse =
      MeshtyingStrategyScaTra()->MasterMatrix();
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsSlave(), 81, false,
              true));
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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

  // extract scatra-scatra interface kinetics conditions
  std::vector<DRT::Condition*> conditions;
  ScaTraField()->Discretization()->GetCondition("S2IKinetics", conditions);

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

  Teuchos::RCP<CORE::LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      mastermatrix =
          MeshtyingStrategyScaTra()
              ->MasterMatrix()
              ->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(
                  *BlockMapThermoInterface(), MeshtyingStrategyScaTra()->BlockMapsMaster());
      mastermatrix->Complete();

      break;
    }
    case CORE::LINALG::MatrixType::sparse:
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
    case CORE::LINALG::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::EvaluateOffDiagBlockThermoScatraInterface(
    Teuchos::RCP<CORE::LINALG::SparseOperator> thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // remove state vectors from thermo discretization
  ThermoField()->Discretization()->ClearState();

  // add state vectors to thermo discretization
  ThermoField()->AddTimeIntegrationSpecificVectors();

  // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes
  // w.r.t. slave-side and master-side scatra dofs
  Teuchos::RCP<CORE::LINALG::SparseOperator> slavematrix(Teuchos::null);
  switch (ScaTraField()->MatrixType())
  {
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
              *ScaTraField()->BlockMaps(), *BlockMapThermoInterface(), 81, false, true));
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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

  // extract scatra-scatra interface kinetics conditions
  std::vector<DRT::Condition*> conditions;
  ThermoField()->Discretization()->GetCondition("S2IKinetics", conditions);

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
    case CORE::LINALG::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
    case CORE::LINALG::MatrixType::block_condition:
    {
      thermoscatrablockinterface->Complete();
      break;
    }

    case CORE::LINALG::MatrixType::sparse:
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> STI::BuildScatraThermoOffDiagCoupling(
    const INPAR::S2I::CouplingType& couplingtype,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> block_map_thermo_interface_slave,
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
