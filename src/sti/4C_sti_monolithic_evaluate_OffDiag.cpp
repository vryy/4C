/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic STI
\level 2


 */
/*----------------------------------------------------------------------*/

#include "4C_sti_monolithic_evaluate_OffDiag.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_discretization_fem_general_assemblestrategy.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_utils_parameter_list.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCoupling::ScatraThermoOffDiagCoupling(
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo)
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
void STI::ScatraThermoOffDiagCoupling::evaluate_off_diag_block_scatra_thermo_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermoblock)
{
  // initialize scatra-thermo matrix block
  scatrathermoblock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_scatra_mono_odblock_scatrathermo, eleparams);

  // remove state vectors from scatra discretization
  sca_tra_field()->discretization()->ClearState();

  // add state vectors to scatra discretization
  sca_tra_field()->add_time_integration_specific_vectors();

  // create strategy for assembly of scatra-thermo matrix block
  Core::FE::AssembleStrategy strategyscatrathermo(
      0,  // row assembly based on number of dofset associated with scatra dofs on scatra
          // discretization
      2,  // column assembly based on number of dofset associated with thermo dofs on scatra
          // discretization
      scatrathermoblock,  // scatra-thermo matrix block
      Teuchos::null,      // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble scatra-thermo matrix block
  sca_tra_field()->discretization()->Evaluate(eleparams, strategyscatrathermo);

  // remove state vectors from scalar transport discretization
  sca_tra_field()->discretization()->ClearState();

  // finalize scatra-thermo block
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatrathermoblock->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      scatrathermoblock->Complete(*full_map_thermo(), *full_map_scatra());
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
void STI::ScatraThermoOffDiagCoupling::evaluate_off_diag_block_thermo_scatra_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatrablock)
{
  // initialize thermo-scatra matrix block
  thermoscatrablock->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList eleparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_scatra_mono_odblock_thermoscatra, eleparams);

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_field()->add_time_integration_specific_vectors();

  // create strategy for assembly of thermo-scatra matrix block
  Core::FE::AssembleStrategy strategythermoscatra(
      0,  // row assembly based on number of dofset associated with thermo dofs on thermo
          // discretization
      2,  // column assembly based on number of dofset associated with scatra dofs on thermo
          // discretization
      thermoscatrablock,  // thermo-scatra matrix block
      Teuchos::null,      // no additional matrices or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // assemble thermo-scatra matrix block
  thermo_field()->discretization()->Evaluate(eleparams, strategythermoscatra);

  // finalize thermo-scatra matrix block
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      thermoscatrablock->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      thermoscatrablock->Complete(*full_map_scatra(), *full_map_thermo());
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMatchingNodes::ScatraThermoOffDiagCouplingMatchingNodes(
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface_slave,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo)
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
void STI::ScatraThermoOffDiagCouplingMatchingNodes::evaluate_off_diag_block_scatra_thermo_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::SparseOperator> mastermatrix(Teuchos::null);
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *block_map_thermo_interface(), meshtying_strategy_sca_tra()->BlockMapsSlave(), 81,
              false, true));
      mastermatrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *block_map_thermo_interface(), meshtying_strategy_sca_tra()->BlockMapsMaster(), 81,
              false, true));
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *meshtying_strategy_sca_tra()->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *meshtying_strategy_sca_tra()->CouplingAdapter()->MasterDofMap(), 27, false, true));
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // evaluate interface contibutions on slave side
  evaluate_scatra_thermo_interface_slave_side(slavematrix);

  // copy interface contributions from slave side to master side
  copy_slave_to_master_scatra_thermo_interface(slavematrix, mastermatrix);

  // add contributions from slave side and master side
  scatrathermoblockinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrathermoblockinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize scatra-thermo matrix block
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(*interface_map_thermo(), *interface_map_sca_tra());
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scalar transport discretization
  sca_tra_field()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::evaluate_scatra_thermo_interface_slave_side(
    Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix)
{
  // zero out slavematrtix
  slavematrix->Zero();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_s2icoupling_od, condparams);

  // set type of differentiation to temperature
  Core::UTILS::AddEnumClassToParameterList<ScaTra::DifferentiationType>(
      "differentiationtype", ScaTra::DifferentiationType::temp, condparams);

  // remove state vectors from scalar transport discretization
  sca_tra_field()->discretization()->ClearState();

  // add state vectors to scalar transport discretization
  sca_tra_field()->add_time_integration_specific_vectors();

  // create strategy for assembly of auxiliary system matrix
  Core::FE::AssembleStrategy strategyscatrathermos2i(
      0,            // row assembly based on number of dofset associated with scatra dofs on scatra
                    // discretization
      2,            // column assembly based on number of dofset associated with thermo dofs on
                    // scatra discretization
      slavematrix,  // auxiliary system matrix
      Teuchos::null,  // no additional matrices of vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface kinetics
  std::vector<Core::Conditions::Condition*> conditions;
  for (const auto& kinetics_slave_cond :
      meshtying_strategy_sca_tra()->kinetics_conditions_meshtying_slave_side())
  {
    if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") !=
        static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_sca_tra()->set_condition_specific_sca_tra_parameters(
          *kinetics_slave_cond.second);
      // evaluate the condition
      sca_tra_field()->discretization()->evaluate_condition(condparams, strategyscatrathermos2i,
          "S2IKinetics", kinetics_slave_cond.second->parameters().Get<int>("ConditionID"));
    }
  }

  // finalize slave matrix
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix->Complete(
          *interface_map_thermo(), *meshtying_strategy_sca_tra()->CouplingAdapter()->SlaveDofMap());
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
void STI::ScatraThermoOffDiagCouplingMatchingNodes::copy_slave_to_master_scatra_thermo_interface(
    Teuchos::RCP<const Core::LinAlg::SparseOperator> slavematrix,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& mastermatrix)
{
  // zero out master matrix
  mastermatrix->Zero();

  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      // cast master and slave matrix
      const auto blockslavematrix =
          Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra
      // fluxes w.r.t. slave-side thermo dofs
      Core::LinAlg::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_sca_tra()->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into auxiliary system matrix
      for (int iblock = 0; iblock < meshtying_strategy_sca_tra()->BlockMapsSlave().NumMaps();
           ++iblock)
      {
        Core::LinAlg::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.0,
            Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_sca_tra()->CouplingAdapter()),
            mastermatrixsparse, true);
      }

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_sca_tra()->CouplingAdapter()->MasterDofMap());

      // split auxiliary system matrix and assemble into scatra-thermo matrix block
      blockmastermatrix = mastermatrixsparse.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *block_map_thermo(), *sca_tra_field()->BlockMaps());

      // finalize master matrix
      mastermatrix->Complete();

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      // cast master and slave matrix
      const auto sparseslavematrix =
          Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. slave-side thermo dofs
      // and assemble into scatra-thermo matrix block
      Core::LinAlg::MatrixRowTransform()(*sparseslavematrix, -1.0,
          Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_sca_tra()->CouplingAdapter()),
          *sparsemastermatrix, false);

      // finalize master matrix
      mastermatrix->Complete(
          *interface_map_thermo(), *meshtying_strategy_sca_tra()->InterfaceMaps()->Map(2));

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");

      break;
    }
  }
  // linearizations of scatra fluxes w.r.t. master-side thermo dofs are not needed,
  // since these dofs will be condensed out later
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMatchingNodes::evaluate_off_diag_block_thermo_scatra_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // initialize slave and master matrix
  Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix(Teuchos::null);
  meshtying_strategy_thermo()->MasterMatrix()->Zero();
  auto mastermatrix = meshtying_strategy_thermo()->MasterMatrix();
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              meshtying_strategy_sca_tra()->BlockMapsSlave(), *block_map_thermo_interface_slave(),
              81, false, true));
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      meshtying_strategy_thermo()->SlaveMatrix()->Zero();
      slavematrix = meshtying_strategy_thermo()->SlaveMatrix();
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_field()->add_time_integration_specific_vectors();

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_s2icoupling_od, condparams);

  // set differentiation type to elch
  Core::UTILS::AddEnumClassToParameterList<ScaTra::DifferentiationType>(
      "differentiationtype", ScaTra::DifferentiationType::elch, condparams);

  // create strategy for assembly of auxiliary system matrices
  Core::FE::AssembleStrategy strategythermoscatras2i(
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
      meshtying_strategy_thermo()->kinetics_conditions_meshtying_slave_side())
  {
    if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") !=
        static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo()->set_condition_specific_sca_tra_parameters(
          *kinetics_slave_cond.second);
      // evaluate the condition
      thermo_field()->discretization()->evaluate_condition(condparams, strategythermoscatras2i,
          "S2IKinetics", kinetics_slave_cond.second->parameters().Get<int>("ConditionID"));
    }
  }

  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      // finalize auxiliary system matrices
      slavematrix->Complete();
      mastermatrix->Complete(*meshtying_strategy_sca_tra()->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

      // initialize temporary matrix
      Core::LinAlg::SparseMatrix ksm(
          *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap(), 27, false, true);

      // transform linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      Core::LinAlg::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(),
          *mastermatrix, 1.0,
          Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_sca_tra()->CouplingAdapter()),
          ksm, true, false);

      // finalize temporary matrix
      ksm.Complete(*meshtying_strategy_sca_tra()->CouplingAdapter()->MasterDofMap(),
          *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap());

      // split temporary matrix and assemble into thermo-scatra matrix block
      const auto blockksm = ksm.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          meshtying_strategy_sca_tra()->BlockMapsMaster(), *block_map_thermo_interface_slave());
      blockksm->Complete();
      thermoscatrablockinterface->Add(*blockksm, false, 1.0, 1.0);

      // finalize matrix
      thermoscatrablockinterface->Complete();

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix->Complete(*meshtying_strategy_sca_tra()->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap());
      mastermatrix->Complete(*meshtying_strategy_sca_tra()->CouplingAdapter()->SlaveDofMap(),
          *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap());

      // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side scatra dofs
      // into thermo-scatra matrix block
      thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

      // derive linearizations of slave-side thermo fluxes w.r.t. master-side scatra dofs
      // and assemble into thermo-scatra matrix block
      Core::LinAlg::MatrixColTransform()(mastermatrix->RowMap(), mastermatrix->ColMap(),
          *mastermatrix, 1.0,
          Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_sca_tra()->CouplingAdapter()),
          *Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(thermoscatrablockinterface), true,
          true);

      // finalize matrix
      thermoscatrablockinterface->Complete(
          *interface_map_sca_tra(), *meshtying_strategy_thermo()->CouplingAdapter()->SlaveDofMap());

      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STI::ScatraThermoOffDiagCouplingMortarStandard::ScatraThermoOffDiagCouplingMortarStandard(
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo)
    : ScatraThermoOffDiagCoupling(std::move(block_map_thermo),
          std::move(block_map_thermo_interface), std::move(full_map_scatra),
          std::move(full_map_thermo), std::move(interface_map_scatra),
          std::move(interface_map_thermo), isale, std::move(meshtying_strategy_scatra),
          std::move(meshtying_strategy_thermo), std::move(scatra), std::move(thermo))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::
    evaluate_off_diag_block_scatra_thermo_interface(
        Teuchos::RCP<Core::LinAlg::SparseOperator> scatrathermoblockinterface)
{
  // zero out matrix
  scatrathermoblockinterface->Zero();

  // initialize auxiliary system matrices for linearizations of slave-side and master-side
  // scatra fluxes w.r.t. slave-side thermo dofs
  Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix(Teuchos::null);
  meshtying_strategy_sca_tra()->MasterMatrix()->Zero();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mastermatrix_sparse =
      meshtying_strategy_sca_tra()->MasterMatrix();
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *block_map_thermo_interface(), meshtying_strategy_sca_tra()->BlockMapsSlave(), 81,
              false, true));
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix = meshtying_strategy_sca_tra()->SlaveMatrix();
      slavematrix->Zero();
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", Inpar::S2I::evaluate_condition_od);

  // create strategy for assembly of auxiliary system matrices
  ScaTra::MortarCellAssemblyStrategy strategyscatrathermos2i(slavematrix, Inpar::S2I::side_slave,
      Inpar::S2I::side_slave, Teuchos::null, Inpar::S2I::side_undefined, Inpar::S2I::side_undefined,
      mastermatrix_sparse, Inpar::S2I::side_master, Inpar::S2I::side_slave, Teuchos::null,
      Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
      Inpar::S2I::side_undefined, Teuchos::null, Inpar::S2I::side_undefined, 0, 1);

  // extract scatra-scatra interface kinetics conditions
  std::vector<Core::Conditions::Condition*> conditions;
  sca_tra_field()->discretization()->GetCondition("S2IKinetics", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->parameters().Get<int>("interface side") == Inpar::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<Core::Conditions::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_sca_tra()->set_condition_specific_sca_tra_parameters(*condition);
      // evaluate mortar integration cells
      meshtying_strategy_sca_tra()->evaluate_mortar_cells(
          meshtying_strategy_sca_tra()->mortar_discretization(
              condition->parameters().Get<int>("ConditionID")),
          condparams, strategyscatrathermos2i);
    }
  }

  // finalize auxiliary system matrices
  mastermatrix_sparse->Complete(
      *interface_map_thermo(), *meshtying_strategy_sca_tra()->InterfaceMaps()->Map(2));

  Teuchos::RCP<Core::LinAlg::SparseOperator> mastermatrix(Teuchos::null);
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix->Complete();
      mastermatrix =
          meshtying_strategy_sca_tra()
              ->MasterMatrix()
              ->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                  *block_map_thermo_interface(), meshtying_strategy_sca_tra()->BlockMapsMaster());
      mastermatrix->Complete();

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix->Complete(
          *interface_map_thermo(), *meshtying_strategy_sca_tra()->InterfaceMaps()->Map(1));
      mastermatrix = mastermatrix_sparse;

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
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
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      scatrathermoblockinterface->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      scatrathermoblockinterface->Complete(*interface_map_thermo(), *interface_map_sca_tra());
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from scatra discretization
  sca_tra_field()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STI::ScatraThermoOffDiagCouplingMortarStandard::
    evaluate_off_diag_block_thermo_scatra_interface(
        Teuchos::RCP<Core::LinAlg::SparseOperator> thermoscatrablockinterface)
{
  // zero out matrix
  thermoscatrablockinterface->Zero();

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_field()->add_time_integration_specific_vectors();

  // initialize auxiliary system matrix for linearizations of slave-side thermo fluxes
  // w.r.t. slave-side and master-side scatra dofs
  Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix(Teuchos::null);
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *sca_tra_field()->BlockMaps(), *block_map_thermo_interface(), 81, false, true));
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix = meshtying_strategy_thermo()->SlaveMatrix();
      slavematrix->Zero();
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // create parameter list for element evaluation
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", Inpar::S2I::evaluate_condition_od);

  // create strategy for assembly of auxiliary system matrix
  ScaTra::MortarCellAssemblyStrategy strategythermoscatras2i(slavematrix, Inpar::S2I::side_slave,
      Inpar::S2I::side_slave, slavematrix, Inpar::S2I::side_slave, Inpar::S2I::side_master,
      Teuchos::null, Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
      Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
      Inpar::S2I::side_undefined, Teuchos::null, Inpar::S2I::side_undefined, 0, 1);

  // extract scatra-scatra interface kinetics conditions
  std::vector<Core::Conditions::Condition*> conditions;
  thermo_field()->discretization()->GetCondition("S2IKinetics", conditions);

  // loop over all conditions
  for (const auto& condition : conditions)
  {
    // consider conditions for slave side only
    if (condition->parameters().Get<int>("interface side") == Inpar::S2I::side_slave)
    {
      // add condition to parameter list
      condparams.set<Core::Conditions::Condition*>("condition", condition);

      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo()->set_condition_specific_sca_tra_parameters(*condition);
      // evaluate mortar integration cells
      meshtying_strategy_thermo()->evaluate_mortar_cells(
          meshtying_strategy_thermo()->mortar_discretization(
              condition->parameters().Get<int>("ConditionID")),
          condparams, strategythermoscatras2i);
    }
  }

  // finalize auxiliary system matrix
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix->Complete(
          *interface_map_sca_tra(), *meshtying_strategy_thermo()->InterfaceMaps()->Map(1));
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // assemble linearizations of slave-side thermo fluxes w.r.t. slave-side and master-side
  // scatra dofs into thermo-scatra matrix block
  thermoscatrablockinterface->Add(*slavematrix, false, 1.0, 1.0);

  // linearizations of master-side thermo fluxes w.r.t. scatra dofs are not needed, since
  // thermo fluxes are source terms and thus only evaluated once on slave side

  // finalize thermo-scatra matrix block
  switch (sca_tra_field()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      thermoscatrablockinterface->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      thermoscatrablockinterface->Complete(*interface_map_sca_tra(), *interface_map_thermo());
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // remove state vectors from thermo discretization
  thermo_field()->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> STI::BuildScatraThermoOffDiagCoupling(
    const Inpar::S2I::CouplingType& couplingtype,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_thermo_interface_slave,
    Teuchos::RCP<const Epetra_Map> full_map_scatra, Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const Epetra_Map> interface_map_thermo, bool isale,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_scatra,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo)
{
  Teuchos::RCP<STI::ScatraThermoOffDiagCoupling> scatrathermooffdiagcoupling = Teuchos::null;

  switch (couplingtype)
  {
    case Inpar::S2I::coupling_matching_nodes:
    {
      scatrathermooffdiagcoupling = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMatchingNodes(
          block_map_thermo, block_map_thermo_interface, block_map_thermo_interface_slave,
          full_map_scatra, full_map_thermo, interface_map_scatra, interface_map_thermo, isale,
          meshtying_strategy_scatra, meshtying_strategy_thermo, scatra, thermo));
      break;
    }
    case Inpar::S2I::coupling_mortar_standard:
    {
      scatrathermooffdiagcoupling = Teuchos::rcp(new STI::ScatraThermoOffDiagCouplingMortarStandard(
          block_map_thermo, block_map_thermo_interface, full_map_scatra, full_map_thermo,
          interface_map_scatra, interface_map_thermo, isale, meshtying_strategy_scatra,
          meshtying_strategy_thermo, scatra, thermo));
      break;
    }
    default:
    {
      FOUR_C_THROW("Not supported coupling type");
      break;
    }
  }

  return scatrathermooffdiagcoupling;
}

FOUR_C_NAMESPACE_CLOSE
