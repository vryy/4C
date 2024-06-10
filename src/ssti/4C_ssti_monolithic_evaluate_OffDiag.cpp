/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSTI

\level 2

*----------------------------------------------------------------------*/
#include "4C_ssti_monolithic_evaluate_OffDiag.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
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
SSTI::ThermoStructureOffDiagCoupling::ThermoStructureOffDiagCoupling(
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapstructure,
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> blockmapthermo,
    Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<const Epetra_Map> full_map_thermo,
    Teuchos::RCP<const SSI::UTILS::SSIMeshTying> ssti_structure_meshtying,
    Teuchos::RCP<const ScaTra::MeshtyingStrategyS2I> meshtying_strategy_thermo,
    Teuchos::RCP<Adapter::SSIStructureWrapper> structure,
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> thermo)
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
void SSTI::ThermoStructureOffDiagCoupling::evaluate_off_diag_block_thermo_structure_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructuredomain)
{
  // initialize thermo-structure matrix block
  thermostructuredomain->Zero();

  Teuchos::ParameterList eleparams;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_scatra_mono_odblock_mesh, eleparams);

  // remove state vectors from thermo discretization
  thermo_->ScaTraField()->discretization()->ClearState();

  // add state vectors to thermo discretization
  thermo_->ScaTraField()->add_time_integration_specific_vectors();

  // create strategy for assembly of thermo-structure matrix block
  Core::FE::AssembleStrategy strategyscatrastructure(
      0,  // row assembly based on number of dofset associated with thermo dofs on thermo
          // discretization
      1,  // column assembly based on number of dofset associated with structural dofs on thermo
          // discretization
      thermostructuredomain,  // thermo-structure matrix block
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  thermo_->ScaTraField()->discretization()->Evaluate(eleparams, strategyscatrastructure);

  thermo_->ScaTraField()->discretization()->ClearState();
}
/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::evaluate_off_diag_block_thermo_structure_interface(
    Teuchos::RCP<Core::LinAlg::SparseOperator> thermostructureinterface)
{
  thermostructureinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<Core::LinAlg::SparseOperator> mastermatrix(Teuchos::null);
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsMaster(), 81, false, true));
      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      slavematrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
      mastermatrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
          *meshtying_strategy_thermo_->CouplingAdapter()->MasterDofMap(), 27, false, true));
      break;
    }
    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  evaluate_thermo_structure_interface_slave_side(slavematrix);

  copy_slave_to_master_thermo_structure_interface(slavematrix, mastermatrix);

  thermostructureinterface->Add(*slavematrix, false, 1.0, 1.0);
  thermostructureinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize thermo-structure matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      thermostructureinterface->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      thermostructureinterface->Complete(*full_map_structure_, *full_map_thermo_);
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  thermo_->ScaTraField()->discretization()->ClearState();
}

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::evaluate_off_diag_block_structure_thermo_domain(
    Teuchos::RCP<Core::LinAlg::SparseOperator> structurethermodomain)
{
  structurethermodomain->Zero();

  Teuchos::ParameterList eleparams;

  eleparams.set("action", "calc_struct_stiffscalar");

  eleparams.set<int>("differentiationtype", static_cast<int>(STR::DifferentiationType::temp));

  eleparams.set<double>("total time", structure_->Time());

  structure_->discretization()->ClearState();

  structure_->discretization()->set_state("displacement", structure_->Dispnp());

  // create strategy for assembly of structure-thermo matrix block
  Core::FE::AssembleStrategy strategystructurescatra(
      0,  // row assembly based on number of dofset associated with structure dofs on structural
          // discretization
      2,  // column assembly based on number of dofset associated with thermo dofs on structural
          // discretization
      structurethermodomain,  // structure-thermo matrix block
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  structure_->discretization()->Evaluate(eleparams, strategystructurescatra);

  // need to scale structurethermoblock_ with 'timefac' to getcorrect implementation
  structurethermodomain->Scale(1.0 - structure_->TimIntParam());

  structure_->discretization()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::copy_slave_to_master_thermo_structure_interface(
    Teuchos::RCP<const Core::LinAlg::SparseOperator> slavematrix,
    Teuchos::RCP<Core::LinAlg::SparseOperator>& mastermatrix)
{
  mastermatrix->Zero();

  switch (thermo_->ScaTraField()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      const int numberthermoblocks = thermo_->ScaTraField()->BlockMaps()->NumMaps();

      auto blockslavematrix =
          Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      Core::LinAlg::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_thermo_->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (int iblock = 0; iblock < numberthermoblocks; ++iblock)
      {
        for (const auto& meshtying : ssti_structure_meshtying_->MeshTyingHandlers())
        {
          auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
          auto slave_side_converter_struct = meshtying->SlaveSideConverter();

          auto slave_side_converter_thermo =
              Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_thermo_->CouplingAdapter());

          Core::LinAlg::MatrixLogicalSplitAndTransform()(blockslavematrix->Matrix(iblock, 0),
              *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), *slave_dof_map, -1.0,
              &slave_side_converter_thermo, &(*slave_side_converter_struct), mastermatrixsparse,
              true, true);
        }
      }

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*full_map_structure_, *full_map_thermo_);

      // split sparse matrix to block matrix
      auto mastermatrix_split = mastermatrixsparse.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, *blockmapthermo_);
      mastermatrix_split->Complete();
      blockmastermatrix->Add(*mastermatrix_split, false, 1.0, 1.0);

      mastermatrix->Complete();

      break;
    }
    case Core::LinAlg::MatrixType::sparse:
    {
      auto sparseslavematrix =
          Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(mastermatrix);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (const auto& meshtying : ssti_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_dof_map = meshtying->SlaveMasterCoupling()->SlaveDofMap();
        auto slave_side_converter_struct = meshtying->SlaveSideConverter();
        auto slave_side_converter_thermo =
            Core::Adapter::CouplingSlaveConverter(*meshtying_strategy_thermo_->CouplingAdapter());

        Core::LinAlg::MatrixLogicalSplitAndTransform()(*sparseslavematrix,
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
      FOUR_C_THROW("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSTI::ThermoStructureOffDiagCoupling::evaluate_thermo_structure_interface_slave_side(
    Teuchos::RCP<Core::LinAlg::SparseOperator> slavematrix)
{
  Teuchos::ParameterList condparams;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_s2icoupling_od, condparams);

  Core::UTILS::AddEnumClassToParameterList<ScaTra::DifferentiationType>(
      "differentiationtype", ScaTra::DifferentiationType::disp, condparams);

  thermo_->ScaTraField()->discretization()->ClearState();

  thermo_->ScaTraField()->add_time_integration_specific_vectors();

  Teuchos::RCP<Core::LinAlg::SparseOperator> evaluate_matrix;
  if (thermo_->ScaTraField()->MatrixType() == Core::LinAlg::MatrixType::sparse)
  {
    evaluate_matrix = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
        *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), 27, false, true));
  }
  else
  {
    evaluate_matrix =
        Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
            *blockmapstructure_, meshtying_strategy_thermo_->BlockMapsSlave(), 81, false, true));
  }

  // create strategy for assembly of auxiliary system matrix
  Core::FE::AssembleStrategy strategyscatrastructures2i(
      0,  // row assembly based on number of dofset associated with thermo dofs on
          // thermo discretization
      1,  // column assembly based on number of dofset associated with structural dofs on
          // thermo discretization
      evaluate_matrix, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate interface coupling
  for (auto kinetics_slave_cond :
      meshtying_strategy_thermo_->kinetics_conditions_meshtying_slave_side())
  {
    if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") !=
        static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux))
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_thermo_->set_condition_specific_sca_tra_parameters(
          *kinetics_slave_cond.second);
      // evaluate the condition
      thermo_->ScaTraField()->discretization()->evaluate_condition(
          condparams, strategyscatrastructures2i, "S2IKinetics", kinetics_slave_cond.first);
    }
  }

  // finalize thermo-structure matrix block
  switch (thermo_->ScaTraField()->MatrixType())
  {
    case Core::LinAlg::MatrixType::block_condition:
    {
      evaluate_matrix->Complete();

      auto evaluate_matrix_block =
          Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(evaluate_matrix);
      auto slavematrix_block =
          Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(slavematrix);

      // "slave side" from thermo and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssti_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->slave_slave_transformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            Core::Adapter::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        for (int iblock = 0; iblock < thermo_->ScaTraField()->BlockMaps()->NumMaps(); ++iblock)
        {
          auto evaluate_iblock = evaluate_matrix_block->Matrix(iblock, 0);
          auto slave_iblock = slavematrix_block->Matrix(iblock, 0);

          auto scatra_slave_block_mapi =
              Core::LinAlg::IntersectMap(*thermo_->ScaTraField()->BlockMaps()->Map(iblock),
                  *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

          Core::LinAlg::MatrixLogicalSplitAndTransform()(evaluate_iblock, *scatra_slave_block_mapi,
              *slave_map, 1.0, nullptr, &slave_slave_converter, slave_iblock, true, true);
        }
      }
      slavematrix->Complete();
      break;
    }

    case Core::LinAlg::MatrixType::sparse:
    {
      auto evaluate_matrix_sparse =
          Core::LinAlg::CastToConstSparseMatrixAndCheckSuccess(evaluate_matrix);
      auto slavematrix_sparse = Core::LinAlg::CastToSparseMatrixAndCheckSuccess(slavematrix);

      evaluate_matrix->Complete(
          *full_map_structure_, *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

      // "slave side" from thermo and from structure do not need to be the same nodes.
      // Linearization is evaluated on scatra slave side node --> Transformation needed
      for (const auto& meshtying : ssti_structure_meshtying_->MeshTyingHandlers())
      {
        auto slave_slave_transformation = meshtying->slave_slave_transformation();
        // converter between old slave dofs from input and actual slave dofs from current mesh tying
        // adapter
        auto slave_slave_converter =
            Core::Adapter::CouplingSlaveConverter(*slave_slave_transformation);

        // old slave dofs from input
        auto slave_map = slave_slave_transformation->SlaveDofMap();

        Core::LinAlg::MatrixLogicalSplitAndTransform()(*evaluate_matrix_sparse,
            *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap(), *slave_map, 1.0, nullptr,
            &slave_slave_converter, *slavematrix_sparse, true, true);
      }
      slavematrix->Complete(
          *full_map_structure_, *meshtying_strategy_thermo_->CouplingAdapter()->SlaveDofMap());

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
