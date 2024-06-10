/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

\level 2


*----------------------------------------------------------------------*/
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"

#include "4C_comm_utils_gid_vector.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_coupling_volmortar_shape.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_mortar_coupling3d_classes.hpp"
#include "4C_mortar_interface.hpp"
#include "4C_mortar_projector.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_boundary_calc.hpp"
#include "4C_scatra_ele_parameter_boundary.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i_elch.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::MeshtyingStrategyS2I::MeshtyingStrategyS2I(
    ScaTra::ScaTraTimIntImpl* scatratimint, const Teuchos::ParameterList& parameters)
    : MeshtyingStrategyBase(scatratimint),
      interfacemaps_(Teuchos::null),
      blockmaps_slave_(Teuchos::null),
      blockmaps_master_(Teuchos::null),
      icoup_(Teuchos::null),
      icoupmortar_(),
      imortarcells_(),
      imortarredistribution_(
          Core::UTILS::IntegralValue<Inpar::S2I::CouplingType>(parameters.sublist("S2I COUPLING"),
              "COUPLINGTYPE") == Inpar::S2I::coupling_mortar_standard and
          Teuchos::getIntegralValue<Inpar::Mortar::ParallelRedist>(
              Global::Problem::Instance()->mortar_coupling_params().sublist(
                  "PARALLEL REDISTRIBUTION"),
              "PARALLEL_REDIST") != Inpar::Mortar::ParallelRedist::redist_none),
      islavemap_(Teuchos::null),
      imastermap_(Teuchos::null),
      islavenodestomasterelements_(),
      islavenodesimpltypes_(),
      islavenodeslumpedareas_(),
      islavematrix_(Teuchos::null),
      imastermatrix_(Teuchos::null),
      imasterslavematrix_(Teuchos::null),
      couplingtype_(Core::UTILS::IntegralValue<Inpar::S2I::CouplingType>(
          parameters.sublist("S2I COUPLING"), "COUPLINGTYPE")),
      D_(Teuchos::null),
      M_(Teuchos::null),
      E_(Teuchos::null),
      P_(Teuchos::null),
      Q_(Teuchos::null),
      lm_(Teuchos::null),
      extendedmaps_(Teuchos::null),
      lmresidual_(Teuchos::null),
      lmincrement_(Teuchos::null),
      islavetomastercoltransform_(Teuchos::null),
      islavetomasterrowtransform_(Teuchos::null),
      islavetomasterrowcoltransform_(Teuchos::null),
      islaveresidual_(Teuchos::null),
      imasterresidual_(Teuchos::null),
      islavephidtnp_(Teuchos::null),
      imasterphidt_on_slave_side_np_(Teuchos::null),
      imasterphi_on_slave_side_np_(Teuchos::null),
      lmside_(Core::UTILS::IntegralValue<Inpar::S2I::InterfaceSides>(
          parameters.sublist("S2I COUPLING"), "LMSIDE")),
      matrixtype_(Teuchos::getIntegralValue<Core::LinAlg::MatrixType>(parameters, "MATRIXTYPE")),
      ntsprojtol_(parameters.sublist("S2I COUPLING").get<double>("NTSPROJTOL")),
      intlayergrowth_evaluation_(Core::UTILS::IntegralValue<Inpar::S2I::GrowthEvaluation>(
          parameters.sublist("S2I COUPLING"), "INTLAYERGROWTH_EVALUATION")),
      intlayergrowth_convtol_(
          parameters.sublist("S2I COUPLING").get<double>("INTLAYERGROWTH_CONVTOL")),
      intlayergrowth_itemax_(parameters.sublist("S2I COUPLING").get<int>("INTLAYERGROWTH_ITEMAX")),
      intlayergrowth_timestep_(
          parameters.sublist("S2I COUPLING").get<double>("INTLAYERGROWTH_TIMESTEP")),
      blockmapgrowth_(Teuchos::null),
      extendedblockmaps_(Teuchos::null),
      extendedsystemmatrix_(Teuchos::null),
      extendedsolver_(Teuchos::null),
      growthn_(Teuchos::null),
      growthnp_(Teuchos::null),
      growthdtn_(Teuchos::null),
      growthdtnp_(Teuchos::null),
      growthhist_(Teuchos::null),
      growthresidual_(Teuchos::null),
      growthincrement_(Teuchos::null),
      scatragrowthblock_(Teuchos::null),
      growthscatrablock_(Teuchos::null),
      growthgrowthblock_(Teuchos::null),
      equilibration_(Teuchos::null),
      output_interface_flux_(Core::UTILS::IntegralValue<bool>(
          parameters.sublist("S2I COUPLING"), "OUTPUT_INTERFACE_FLUX")),
      has_capacitive_contributions_(false),
      kinetics_conditions_meshtying_slaveside_(),
      slaveonly_(Core::UTILS::IntegralValue<bool>(parameters.sublist("S2I COUPLING"), "SLAVEONLY")),
      indepedent_setup_of_conditions_(Core::UTILS::IntegralValue<bool>(
          parameters.sublist("S2I COUPLING"), "MESHTYING_CONDITIONS_INDEPENDENT_SETUP"))
{
  // empty constructor
}  // ScaTra::meshtying_strategy_s2_i::meshtying_strategy_s2_i


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::CondenseMatAndRHS(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& residual, const bool calcinittimederiv) const
{
  switch (couplingtype_)
  {
    case Inpar::S2I::coupling_mortar_condensed_bubnov:
    {
      // extract global system matrix
      Teuchos::RCP<Core::LinAlg::SparseMatrix> sparsematrix = scatratimint_->SystemMatrix();
      if (sparsematrix == Teuchos::null) FOUR_C_THROW("System matrix is not a sparse matrix!");

      if (lmside_ == Inpar::S2I::side_slave)
      {
        // initialize temporary matrix for slave-side rows of global system matrix
        Core::LinAlg::SparseMatrix sparsematrixrowsslave(*interfacemaps_->Map(1), 81);

        // extract slave-side rows of global system matrix into temporary matrix
        ExtractMatrixRows(*sparsematrix, sparsematrixrowsslave, *interfacemaps_->Map(1));

        // finalize temporary matrix with slave-side rows of global system matrix
        sparsematrixrowsslave.Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(1));

        // zero out slave-side rows of global system matrix after having extracted them into
        // temporary matrix
        sparsematrix->Complete();
        sparsematrix->ApplyDirichlet(*interfacemaps_->Map(1), false);

        // apply scatra-scatra interface coupling
        if (not slaveonly_)
        {
          // replace slave-side rows of global system matrix by projected slave-side rows including
          // interface contributions
          sparsematrix->Add(*Core::LinAlg::MLMultiply(
                                *Q_, true, sparsematrixrowsslave, false, false, false, true),
              false, 1., 1.);
          // during calculation of initial time derivative, standard global system matrix is
          // replaced by global mass matrix, and hence interface contributions must not be included
          if (!calcinittimederiv) sparsematrix->Add(*islavematrix_, false, 1., 1.);
        }

        // apply standard meshtying
        else
        {
          sparsematrix->Add(*D_, false, 1., 1.);
          sparsematrix->Add(*M_, false, -1., 1.);
        }

        // add projected slave-side rows to master-side rows of global system matrix
        sparsematrix->Add(
            *Core::LinAlg::MLMultiply(*P_, true, sparsematrixrowsslave, false, false, false, true),
            false, 1., 1.);

        // extract slave-side entries of global residual vector
        Teuchos::RCP<Epetra_Vector> residualslave =
            interfacemaps_->ExtractVector(scatratimint_->Residual(), 1);

        // apply scatra-scatra interface coupling
        if (not slaveonly_)
        {
          // replace slave-side entries of global residual vector by projected slave-side entries
          // including interface contributions
          Epetra_Vector Q_residualslave(*interfacemaps_->Map(1));
          if (Q_->Multiply(true, *residualslave, Q_residualslave))
            FOUR_C_THROW("Matrix-vector multiplication failed!");
          interfacemaps_->InsertVector(Q_residualslave, 1, *scatratimint_->Residual());
          interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
        }

        // apply standard meshtying
        else
        {
          // zero out slave-side entries of global residual vector
          interfacemaps_->PutScalar(*scatratimint_->Residual(), 1, 0.);
        }

        // add projected slave-side entries to master-side entries of global residual vector
        Epetra_Vector P_residualslave(*interfacemaps_->Map(2));
        if (P_->Multiply(true, *residualslave, P_residualslave))
          FOUR_C_THROW("Matrix-vector multiplication failed!");
        interfacemaps_->AddVector(P_residualslave, 2, *scatratimint_->Residual());
      }

      else
      {
        // initialize temporary matrix for master-side rows of global system matrix
        Core::LinAlg::SparseMatrix sparsematrixrowsmaster(*interfacemaps_->Map(2), 81);

        // extract master-side rows of global system matrix into temporary matrix
        ExtractMatrixRows(*sparsematrix, sparsematrixrowsmaster, *interfacemaps_->Map(2));

        // finalize temporary matrix with master-side rows of global system matrix
        sparsematrixrowsmaster.Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(2));

        // zero out master-side rows of global system matrix after having extracted them into
        // temporary matrix and replace them by projected master-side rows including interface
        // contributions
        sparsematrix->Complete();
        sparsematrix->ApplyDirichlet(*interfacemaps_->Map(2), false);
        sparsematrix->Add(
            *Core::LinAlg::MLMultiply(*Q_, true, sparsematrixrowsmaster, false, false, false, true),
            false, 1., 1.);
        // during calculation of initial time derivative, standard global system matrix is replaced
        // by global mass matrix, and hence interface contributions must not be included
        if (!calcinittimederiv) sparsematrix->Add(*imastermatrix_, false, 1., 1.);

        // add projected master-side rows to slave-side rows of global system matrix
        sparsematrix->Add(
            *Core::LinAlg::MLMultiply(*P_, true, sparsematrixrowsmaster, false, false, false, true),
            false, 1., 1.);

        // extract master-side entries of global residual vector
        Teuchos::RCP<Epetra_Vector> residualmaster =
            interfacemaps_->ExtractVector(scatratimint_->Residual(), 2);

        // replace master-side entries of global residual vector by projected master-side entries
        // including interface contributions
        Epetra_Vector Q_residualmaster(*interfacemaps_->Map(2));
        if (Q_->Multiply(true, *residualmaster, Q_residualmaster))
          FOUR_C_THROW("Matrix-vector multiplication failed!");
        interfacemaps_->InsertVector(Q_residualmaster, 2, *scatratimint_->Residual());
        interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());

        // add projected master-side entries to slave-side entries of global residual vector
        Epetra_Vector P_residualmaster(*interfacemaps_->Map(1));
        if (P_->Multiply(true, *residualmaster, P_residualmaster))
          FOUR_C_THROW("Matrix-vector multiplication failed!");
        interfacemaps_->AddVector(P_residualmaster, 1, *scatratimint_->Residual());
      }

      break;
    }

    case Inpar::S2I::coupling_matching_nodes:
    case Inpar::S2I::coupling_mortar_standard:
    case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
    case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
    case Inpar::S2I::coupling_mortar_condensed_petrov:
    case Inpar::S2I::coupling_nts_standard:
    {
      // do nothing in these cases
      break;
    }

    default:
    {
      FOUR_C_THROW("Type of mortar meshtying for scatra-scatra interface coupling not recognized!");
      break;
    }
  }
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
const Epetra_Map& ScaTra::MeshtyingStrategyS2I::dof_row_map() const
{
  return extendedmaps_ != Teuchos::null ? *extendedmaps_->FullMap() : *scatratimint_->dof_row_map();
}


/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::EvaluateMeshtying()
{
  // time measurement: evaluate condition 'S2IMeshtying'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'S2IMeshtying'");

  switch (couplingtype_)
  {
    case Inpar::S2I::coupling_matching_nodes:
    {
      // create parameter list
      Teuchos::ParameterList condparams;

      // action for elements
      Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
          "action", ScaTra::BoundaryAction::calc_s2icoupling, condparams);

      // set global state vectors according to time-integration scheme
      scatratimint_->add_time_integration_specific_vectors();

      // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
      islavematrix_->Zero();
      if (not slaveonly_) imastermatrix_->Zero();
      islaveresidual_->PutScalar(0.);
      for (auto kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
      {
        if (kinetics_slave_cond.second->parameters().Get<int>("kinetic model") !=
                static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux) and
            kinetics_slave_cond.second->GType() != Core::Conditions::geometry_type_point)
        {
          // collect condition specific data and store to scatra boundary parameter class
          set_condition_specific_sca_tra_parameters(*kinetics_slave_cond.second);

          if (not slaveonly_)
          {
            scatratimint_->discretization()->evaluate_condition(condparams, islavematrix_,
                imastermatrix_, islaveresidual_, Teuchos::null, Teuchos::null, "S2IKinetics",
                kinetics_slave_cond.second->parameters().Get<int>("ConditionID"));
          }
          else
          {
            scatratimint_->discretization()->evaluate_condition(condparams, islavematrix_,
                Teuchos::null, islaveresidual_, Teuchos::null, Teuchos::null, "S2IKinetics",
                kinetics_slave_cond.second->parameters().Get<int>("ConditionID"));
          }
        }
      }

      // finalize interface matrices
      islavematrix_->Complete();
      if (not slaveonly_) imastermatrix_->Complete();

      // assemble global system matrix depending on matrix type
      switch (matrixtype_)
      {
        case Core::LinAlg::MatrixType::sparse:
        {
          // check matrix
          Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix = scatratimint_->SystemMatrix();
          FOUR_C_ASSERT(systemmatrix != Teuchos::null, "System matrix is not a sparse matrix!");

          // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
          systemmatrix->Add(*islavematrix_, false, 1., 1.);

          if (not slaveonly_)
          {
            // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomastercoltransform_)(imastermatrix_->RowMap(), imastermatrix_->ColMap(),
                *imastermatrix_, 1., Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix,
                true, true);

            // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global
            // system matrix
            (*islavetomasterrowtransform_)(*islavematrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true);

            // derive linearizations of master fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomasterrowcoltransform_)(*imastermatrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_),
                Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true, true);
          }

          // In case the interface linearizations and residuals are evaluated on slave side only,
          // we now apply a standard meshtying algorithm to condense out the slave-side degrees of
          // freedom.
          else if (!scatratimint_->discretization()->GetCondition("PointCoupling"))
          {
            // initialize temporary matrix for slave-side rows of system matrix
            Core::LinAlg::SparseMatrix systemmatrixrowsslave(*icoup_->SlaveDofMap(), 81);

            // extract slave-side rows of system matrix into temporary matrix
            ExtractMatrixRows(*systemmatrix, systemmatrixrowsslave, *icoup_->SlaveDofMap());

            // zero out slave-side rows of system matrix and put a one on the main diagonal
            systemmatrix->Complete();
            systemmatrix->ApplyDirichlet(*icoup_->SlaveDofMap(), true);
            systemmatrix->UnComplete();

            // loop over all slave-side rows of system matrix
            for (int slavedoflid = 0; slavedoflid < icoup_->SlaveDofMap()->NumMyElements();
                 ++slavedoflid)
            {
              // determine global ID of current matrix row
              const int slavedofgid = icoup_->SlaveDofMap()->GID(slavedoflid);
              if (slavedofgid < 0) FOUR_C_THROW("Couldn't find local ID %d in map!", slavedoflid);

              // determine global ID of associated master-side matrix column
              const int masterdofgid = icoup_->PermMasterDofMap()->GID(slavedoflid);
              if (masterdofgid < 0)
                FOUR_C_THROW("Couldn't find local ID %d in permuted map!", slavedoflid);

              // insert value -1. into intersection of slave-side row and master-side column in
              // system matrix this effectively forces the slave-side degree of freedom to assume
              // the same value as the master-side degree of freedom
              const double value(-1.);
              if (systemmatrix->EpetraMatrix()->InsertGlobalValues(
                      slavedofgid, 1, &value, &masterdofgid) < 0)
              {
                FOUR_C_THROW(
                    "Cannot insert value -1. into matrix row with global ID %d and matrix column "
                    "with global ID %d!",
                    slavedofgid, masterdofgid);
              }

              // insert zero into intersection of slave-side row and master-side column in temporary
              // matrix this prevents the system matrix from changing its graph when calling this
              // function again during the next Newton iteration
              const double zero(0.);
              if (systemmatrixrowsslave.EpetraMatrix()->InsertGlobalValues(
                      slavedofgid, 1, &zero, &masterdofgid) < 0)
              {
                FOUR_C_THROW(
                    "Cannot insert zero into matrix row with global ID %d and matrix column with "
                    "global ID %d!",
                    slavedofgid, masterdofgid);
              }
            }

            // finalize temporary matrix with slave-side rows of system matrix
            systemmatrixrowsslave.Complete(*scatratimint_->dof_row_map(), *icoup_->SlaveDofMap());

            // add slave-side rows of system matrix to corresponding master-side rows to finalize
            // matrix condensation of slave-side degrees of freedom
            (*islavetomasterrowtransform_)(systemmatrixrowsslave, 1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true);
          }
          break;
        }
        case Core::LinAlg::MatrixType::block_condition:
        case Core::LinAlg::MatrixType::block_condition_dof:
        {
          // check matrix
          Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blocksystemmatrix =
              scatratimint_->BlockSystemMatrix();
          FOUR_C_ASSERT(blocksystemmatrix != Teuchos::null, "System matrix is not a block matrix!");

          Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkss(
              islavematrix_->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                  *blockmaps_slave_, *blockmaps_slave_));
          blockkss->Complete();

          // assemble interface block matrix into global block system matrix
          blocksystemmatrix->Add(*blockkss, false, 1., 1.);

          if (not slaveonly_)
          {
            Teuchos::RCP<Core::LinAlg::SparseMatrix> ksm(
                Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->SlaveDofMap(), 81, false)));
            Teuchos::RCP<Core::LinAlg::SparseMatrix> kms(
                Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm(
                Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));

            // transform linearizations of slave fluxes w.r.t. master dofs
            (*islavetomastercoltransform_)(imastermatrix_->RowMap(), imastermatrix_->ColMap(),
                *imastermatrix_, 1., Core::Adapter::CouplingSlaveConverter(*icoup_), *ksm);
            ksm->Complete(*icoup_->MasterDofMap(), *icoup_->SlaveDofMap());

            // derive linearizations of master fluxes w.r.t. slave dofs
            (*islavetomasterrowtransform_)(
                *islavematrix_, -1., Core::Adapter::CouplingSlaveConverter(*icoup_), *kms);
            kms->Complete(*icoup_->SlaveDofMap(), *icoup_->MasterDofMap());

            // derive linearizations of master fluxes w.r.t. master dofs
            (*islavetomasterrowcoltransform_)(*imastermatrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_),
                Core::Adapter::CouplingSlaveConverter(*icoup_), *kmm);
            kmm->Complete();

            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockksm(
                ksm->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_slave_));
            blockksm->Complete();
            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkms(
                kms->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                    *blockmaps_slave_, *blockmaps_master_));
            blockkms->Complete();
            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkmm(
                kmm->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_master_));
            blockkmm->Complete();

            // assemble interface block matrices into global block system matrix
            blocksystemmatrix->Add(*blockksm, false, 1., 1.);
            blocksystemmatrix->Add(*blockkms, false, 1., 1.);
            blocksystemmatrix->Add(*blockkmm, false, 1., 1.);
          }

          // safety check
          else
          {
            FOUR_C_THROW(
                "Scatra-scatra interface coupling with evaluation of interface linearizations and "
                "residuals on slave side only is not yet available for block system matrices!");
          }

          break;
        }

        default:
        {
          FOUR_C_THROW(
              "Type of global system matrix for scatra-scatra interface coupling not recognized!");
          break;
        }
      }

      // assemble slave residuals into global residual vector
      interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());

      if (not slaveonly_)
      {
        // transform master residuals and assemble into global residual vector
        interfacemaps_->AddVector(
            icoup_->SlaveToMaster(islaveresidual_), 2, scatratimint_->Residual(), -1.);
      }
      // In case the interface linearizations and residuals are evaluated on slave side only,
      // we now apply a standard meshtying algorithm to condense out the slave-side degrees of
      // freedom.
      else if (!scatratimint_->discretization()->GetCondition("PointCoupling"))
      {
        // initialize temporary vector for slave-side entries of residual vector
        Teuchos::RCP<Epetra_Vector> residualslave =
            Teuchos::rcp(new Epetra_Vector(*icoup_->SlaveDofMap()));

        // loop over all slave-side entries of residual vector
        for (int slavedoflid = 0; slavedoflid < icoup_->SlaveDofMap()->NumMyElements();
             ++slavedoflid)
        {
          // determine global ID of current vector entry
          const int slavedofgid = icoup_->SlaveDofMap()->GID(slavedoflid);
          if (slavedofgid < 0) FOUR_C_THROW("Couldn't find local ID %d in map!", slavedoflid);

          // copy current vector entry into temporary vector
          if (residualslave->ReplaceGlobalValue(slavedofgid, 0,
                  (*scatratimint_->Residual())[scatratimint_->dof_row_map()->LID(slavedofgid)]))
            FOUR_C_THROW(
                "Cannot insert residual vector entry with global ID %d into temporary vector!",
                slavedofgid);

          // zero out current vector entry
          if (scatratimint_->Residual()->ReplaceGlobalValue(slavedofgid, 0, 0.))
            FOUR_C_THROW(
                "Cannot insert zero into residual vector entry with global ID %d!", slavedofgid);
        }

        // add slave-side entries of residual vector to corresponding master-side entries to
        // finalize vector condensation of slave-side degrees of freedom
        interfacemaps_->AddVector(
            icoup_->SlaveToMaster(residualslave), 2, scatratimint_->Residual());
      }

      if (has_capacitive_contributions_) evaluate_and_assemble_capacitive_contributions();

      break;
    }

    case Inpar::S2I::coupling_mortar_standard:
    case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
    case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
    case Inpar::S2I::coupling_mortar_condensed_petrov:
    case Inpar::S2I::coupling_mortar_condensed_bubnov:
    case Inpar::S2I::coupling_nts_standard:
    {
      // initialize auxiliary system matrix and vector for slave side
      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_slave or couplingtype_ == Inpar::S2I::coupling_nts_standard)
      {
        islavematrix_->Zero();
        islaveresidual_->PutScalar(0.);
      }

      // initialize auxiliary system matrix and vector for master side
      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_master or couplingtype_ == Inpar::S2I::coupling_nts_standard)
      {
        imastermatrix_->Zero();
        imasterresidual_->PutScalar(0.);
      }

      // loop over all scatra-scatra coupling interfaces
      for (auto& kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
      {
        // extract mortar interface discretization
        Core::FE::Discretization& idiscret =
            icoupmortar_[kinetics_slave_cond.first]->Interface()->Discret();

        // export global state vector to mortar interface
        Teuchos::RCP<Epetra_Vector> iphinp =
            Teuchos::rcp(new Epetra_Vector(*idiscret.DofColMap(), false));
        Core::LinAlg::Export(*scatratimint_->Phiafnp(), *iphinp);
        idiscret.set_state("iphinp", iphinp);

        // create parameter list for mortar integration cells
        Teuchos::ParameterList params;

        // add current condition to parameter list
        params.set<Core::Conditions::Condition*>("condition", kinetics_slave_cond.second);

        // collect condition specific data and store to scatra boundary parameter class
        set_condition_specific_sca_tra_parameters(*(kinetics_slave_cond.second));

        if (couplingtype_ != Inpar::S2I::coupling_nts_standard)
        {
          // set action
          params.set<int>("action", Inpar::S2I::evaluate_condition);

          // evaluate mortar integration cells at current interface
          evaluate_mortar_cells(idiscret, params, islavematrix_, Inpar::S2I::side_slave,
              Inpar::S2I::side_slave, islavematrix_, Inpar::S2I::side_slave,
              Inpar::S2I::side_master, imastermatrix_, Inpar::S2I::side_master,
              Inpar::S2I::side_slave, imastermatrix_, Inpar::S2I::side_master,
              Inpar::S2I::side_master, islaveresidual_, Inpar::S2I::side_slave, imasterresidual_,
              Inpar::S2I::side_master);
        }

        else
        {
          // set action
          params.set<int>("action", Inpar::S2I::evaluate_condition_nts);

          // evaluate note-to-segment coupling at current interface
          evaluate_nts(*islavenodestomasterelements_[kinetics_slave_cond.first],
              *islavenodeslumpedareas_[kinetics_slave_cond.first],
              *islavenodesimpltypes_[kinetics_slave_cond.first], idiscret, params, islavematrix_,
              Inpar::S2I::side_slave, Inpar::S2I::side_slave, islavematrix_, Inpar::S2I::side_slave,
              Inpar::S2I::side_master, imastermatrix_, Inpar::S2I::side_master,
              Inpar::S2I::side_slave, imastermatrix_, Inpar::S2I::side_master,
              Inpar::S2I::side_master, islaveresidual_, Inpar::S2I::side_slave, imasterresidual_,
              Inpar::S2I::side_master);
        }
      }

      // finalize auxiliary system matrix for slave side
      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_slave or couplingtype_ == Inpar::S2I::coupling_nts_standard)
        islavematrix_->Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(1));

      // finalize auxiliary system matrix and residual vector for master side
      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_master or couplingtype_ == Inpar::S2I::coupling_nts_standard)
      {
        imastermatrix_->Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(2));
        if (imasterresidual_->GlobalAssemble(Add, true))
          FOUR_C_THROW(
              "Assembly of auxiliary residual vector for master residuals not successful!");
      }

      // assemble global system of equations depending on matrix type
      switch (matrixtype_)
      {
        case Core::LinAlg::MatrixType::sparse:
        {
          // extract global system matrix from time integrator
          const Teuchos::RCP<Core::LinAlg::SparseMatrix> systemmatrix =
              scatratimint_->SystemMatrix();
          if (systemmatrix == Teuchos::null) FOUR_C_THROW("System matrix is not a sparse matrix!");

          // assemble interface contributions into global system of equations
          switch (couplingtype_)
          {
            case Inpar::S2I::coupling_mortar_standard:
            case Inpar::S2I::coupling_nts_standard:
            {
              const Teuchos::RCP<const Core::LinAlg::SparseMatrix> islavematrix =
                  not imortarredistribution_
                      ? islavematrix_
                      : Mortar::MatrixRowTransform(islavematrix_, islavemap_);
              const Teuchos::RCP<const Core::LinAlg::SparseMatrix> imastermatrix =
                  not imortarredistribution_
                      ? imastermatrix_
                      : Mortar::MatrixRowTransform(imastermatrix_, imastermap_);
              systemmatrix->Add(*islavematrix, false, 1., 1.);
              systemmatrix->Add(*imastermatrix, false, 1., 1.);
              interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
              interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());

              break;
            }

            case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
            case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
            {
              if (lmside_ == Inpar::S2I::side_slave)
              {
                // assemble slave-side interface contributions into global residual vector
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (D_->Multiply(true, *lm_, islaveresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual(), -1.);

                // assemble master-side interface contributions into global residual vector
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (M_->Multiply(true, *lm_, imasterresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual());

                // build constraint residual vector associated with Lagrange multiplier dofs
                Epetra_Vector ilmresidual(*islaveresidual_);
                if (ilmresidual.ReplaceMap(*extendedmaps_->Map(1)))
                  FOUR_C_THROW("Couldn't replace map!");
                if (lmresidual_->Update(1., ilmresidual, 0.)) FOUR_C_THROW("Vector update failed!");
                if (E_->Multiply(true, *lm_, ilmresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                if (lmresidual_->Update(1., ilmresidual, 1.)) FOUR_C_THROW("Vector update failed!");
              }
              else
              {
                // assemble slave-side interface contributions into global residual vector
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (M_->Multiply(true, *lm_, islaveresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual());

                // assemble master-side interface contributions into global residual vector
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (D_->Multiply(true, *lm_, imasterresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual(), -1.);

                // build constraint residual vector associated with Lagrange multiplier dofs
                Epetra_Vector ilmresidual(Copy, *imasterresidual_, 0);
                if (ilmresidual.ReplaceMap(*extendedmaps_->Map(1)))
                  FOUR_C_THROW("Couldn't replace map!");
                if (lmresidual_->Update(1., ilmresidual, 0.)) FOUR_C_THROW("Vector update failed!");
                if (E_->Multiply(true, *lm_, ilmresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                if (lmresidual_->Update(1., ilmresidual, 1.)) FOUR_C_THROW("Vector update failed!");
              }

              break;
            }

            case Inpar::S2I::coupling_mortar_condensed_petrov:
            {
              if (lmside_ == Inpar::S2I::side_slave)
              {
                systemmatrix->Add(*islavematrix_, false, 1., 1.);
                systemmatrix->Add(
                    *Core::LinAlg::MLMultiply(*P_, true, *islavematrix_, false, false, false, true),
                    false, -1., 1.);
                interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (P_->Multiply(true, *islaveresidual_, imasterresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual(), -1.);
              }
              else
              {
                systemmatrix->Add(*Core::LinAlg::MLMultiply(
                                      *P_, true, *imastermatrix_, false, false, false, true),
                    false, -1., 1.);
                systemmatrix->Add(*imastermatrix_, false, 1., 1.);
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (P_->Multiply(true, *imasterresidual_, islaveresidual))
                  FOUR_C_THROW("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual(), -1.);
                interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());
              }

              break;
            }

            case Inpar::S2I::coupling_mortar_condensed_bubnov:
            {
              // assemble interface contributions into global system of equations
              if (slaveonly_)
              {
                systemmatrix->Add(*islavematrix_, false, 1., 1.);
                interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
              }

              // during calculation of initial time derivative, condensation must not be performed
              // here, but after assembly of the modified global system of equations
              if (scatratimint_->Step() > 0)
                CondenseMatAndRHS(systemmatrix, scatratimint_->Residual());

              break;
            }

            default:
            {
              FOUR_C_THROW("Not yet implemented!");
              break;
            }
          }

          break;
        }

        case Core::LinAlg::MatrixType::block_condition:
        {
          // extract global system matrix from time integrator
          Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blocksystemmatrix =
              scatratimint_->BlockSystemMatrix();
          if (blocksystemmatrix == Teuchos::null)
            FOUR_C_THROW("System matrix is not a block matrix!");

          // assemble interface contributions into global system of equations
          switch (couplingtype_)
          {
            case Inpar::S2I::coupling_mortar_standard:
            {
              // split interface sparse matrices into block matrices
              const Teuchos::RCP<const Core::LinAlg::SparseMatrix> islavematrix =
                  not imortarredistribution_
                      ? islavematrix_
                      : Mortar::MatrixRowTransform(islavematrix_, islavemap_);
              Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockslavematrix(
                  islavematrix->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                      *scatratimint_->BlockMaps(), *blockmaps_slave_));
              blockslavematrix->Complete();
              const Teuchos::RCP<const Core::LinAlg::SparseMatrix> imastermatrix =
                  not imortarredistribution_
                      ? imastermatrix_
                      : Mortar::MatrixRowTransform(imastermatrix_, imastermap_);
              Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockmastermatrix(
                  imastermatrix->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                      *scatratimint_->BlockMaps(), *blockmaps_master_));
              blockmastermatrix->Complete();

              // assemble interface block matrices into global block system matrix
              blocksystemmatrix->Add(*blockslavematrix, false, 1., 1.);
              blocksystemmatrix->Add(*blockmastermatrix, false, 1., 1.);

              // assemble interface residual vectors into global residual vector
              interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
              interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());

              break;
            }

            default:
            {
              FOUR_C_THROW("Not yet implemented!");
              break;
            }
          }

          break;
        }

        default:
        {
          FOUR_C_THROW("Not yet implemented!");
          break;
        }
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Not yet implemented!");
      break;
    }
  }
  // extract boundary conditions for scatra-scatra interface layer growth
  std::vector<Core::Conditions::Condition*> s2icoupling_growth_conditions;
  scatratimint_->discretization()->GetCondition("S2IKineticsGrowth", s2icoupling_growth_conditions);

  // evaluate scatra-scatra interface layer growth
  if (s2icoupling_growth_conditions.size())
  {
    switch (couplingtype_)
    {
      case Inpar::S2I::coupling_matching_nodes:
      {
        // create parameter list for elements
        Teuchos::ParameterList conditionparams;

        // action for elements
        Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
            "action", ScaTra::BoundaryAction::calc_s2icoupling, conditionparams);

        // set global state vectors according to time-integration scheme
        scatratimint_->add_time_integration_specific_vectors();

        // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
        islavematrix_->Zero();
        imastermatrix_->Zero();
        islaveresidual_->PutScalar(0.);

        // collect condition specific data and store to scatra boundary parameter class
        set_condition_specific_sca_tra_parameters(*s2icoupling_growth_conditions[0]);
        // evaluate the condition
        scatratimint_->discretization()->evaluate_condition(conditionparams, islavematrix_,
            imastermatrix_, islaveresidual_, Teuchos::null, Teuchos::null, "S2IKineticsGrowth");

        // finalize interface matrices
        islavematrix_->Complete();
        imastermatrix_->Complete();

        // assemble interface matrices into global system matrix depending on matrix type
        switch (matrixtype_)
        {
          case Core::LinAlg::MatrixType::sparse:
          {
            // check matrix
            const Teuchos::RCP<Core::LinAlg::SparseMatrix>& systemmatrix =
                scatratimint_->SystemMatrix();
            if (systemmatrix == Teuchos::null)
              FOUR_C_THROW("System matrix is not a sparse matrix!");

            // We assume that the scatra-scatra interface layer growth is caused by master-side
            // fluxes to the interface, whereas there is no mass exchange between the interface
            // layer and the slave side of the interface. Hence, we only need to linearize the
            // master-side fluxes w.r.t. the slave-side and master-side degrees of freedom.

            // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global
            // system matrix
            Core::LinAlg::MatrixRowTransform()(*islavematrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true);

            // derive linearizations of master fluxes w.r.t. master dofs and assemble into global
            // system matrix
            Core::LinAlg::MatrixRowColTransform()(*imastermatrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_),
                Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true, true);

            break;
          }

          case Core::LinAlg::MatrixType::block_condition:
          case Core::LinAlg::MatrixType::block_condition_dof:
          {
            // check matrix
            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blocksystemmatrix =
                scatratimint_->BlockSystemMatrix();
            if (blocksystemmatrix == Teuchos::null)
              FOUR_C_THROW("System matrix is not a block matrix!");

            // We assume that the scatra-scatra interface layer growth is caused by master-side
            // fluxes to the interface, whereas there is no mass exchange between the interface
            // layer and the slave side of the interface. Hence, we only need to linearize the
            // master-side fluxes w.r.t. the slave-side and master-side degrees of freedom.

            // derive linearizations of master fluxes w.r.t. slave dofs
            Teuchos::RCP<Core::LinAlg::SparseMatrix> kms(
                Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            Core::LinAlg::MatrixRowTransform()(
                *islavematrix_, -1., Core::Adapter::CouplingSlaveConverter(*icoup_), *kms);
            kms->Complete(*icoup_->SlaveDofMap(), *icoup_->MasterDofMap());
            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkms(
                kms->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                    *blockmaps_slave_, *blockmaps_master_));
            blockkms->Complete();

            // derive linearizations of master fluxes w.r.t. master dofs
            Teuchos::RCP<Core::LinAlg::SparseMatrix> kmm(
                Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            Core::LinAlg::MatrixRowColTransform()(*imastermatrix_, -1.,
                Core::Adapter::CouplingSlaveConverter(*icoup_),
                Core::Adapter::CouplingSlaveConverter(*icoup_), *kmm);
            kmm->Complete();
            Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkmm(
                kmm->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_master_));
            blockkmm->Complete();

            // assemble interface block matrices into global block system matrix
            blocksystemmatrix->Add(*blockkms, false, 1., 1.);
            blocksystemmatrix->Add(*blockkmm, false, 1., 1.);

            break;
          }

          default:
          {
            FOUR_C_THROW(
                "Type of global system matrix for scatra-scatra interface coupling involving "
                "interface layer growth not recognized!");
            break;
          }
        }

        // As before, we only need to consider residual contributions from the master-side fluxes.

        // transform master residuals and assemble into global residual vector
        interfacemaps_->AddVector(
            icoup_->SlaveToMaster(islaveresidual_), 2, scatratimint_->Residual(), -1.);

        // compute additional linearizations and residuals in case of monolithic evaluation approach
        if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
        {
          // extract map associated with scalar transport degrees of freedom
          const Epetra_Map& dofrowmap_scatra = *scatratimint_->discretization()->dof_row_map();

          // extract map associated with scatra-scatra interface layer thicknesses
          const Epetra_Map& dofrowmap_growth = *scatratimint_->discretization()->dof_row_map(2);

          // extract ID of boundary condition for scatra-scatra interface layer growth
          // the corresponding boundary condition for scatra-scatra interface coupling is expected
          // to have the same ID
          const int condid = scatratimint_->discretization()
                                 ->GetCondition("S2IKineticsGrowth")
                                 ->parameters()
                                 .Get<int>("ConditionID");

          // set global state vectors according to time-integration scheme
          scatratimint_->add_time_integration_specific_vectors();

          // compute additional linearizations and residuals depending on type of scalar transport
          // system matrix
          switch (matrixtype_)
          {
            case Core::LinAlg::MatrixType::sparse:
            {
              // assemble off-diagonal scatra-growth block of global system matrix, containing
              // derivatives of discrete scatra residuals w.r.t. discrete scatra-scatra interface
              // layer thicknesses
              {
                // check matrix
                const Teuchos::RCP<Core::LinAlg::SparseMatrix> scatragrowthblock =
                    Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(scatragrowthblock_);
                if (scatragrowthblock == Teuchos::null)
                  FOUR_C_THROW("Matrix is not a sparse matrix!");

                // initialize matrix block
                scatragrowthblock->Zero();

                // initialize auxiliary matrix block for linearizations of slave fluxes w.r.t.
                // scatra-scatra interface layer thicknesses
                Teuchos::RCP<Core::LinAlg::SparseMatrix> islavematrix =
                    Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(icoup_)->SlaveDofMap(), 81));

                // initialize assembly strategy for auxiliary matrix block
                Core::FE::AssembleStrategy strategy(
                    0, 2, islavematrix, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
                    "action", ScaTra::BoundaryAction::calc_s2icoupling_scatragrowth, condparams);

                // evaluate off-diagonal linearizations arising from scatra-scatra interface
                // coupling
                for (auto [kinetics_slave_cond_id, kinetics_slave_cond] :
                    kinetics_conditions_meshtying_slaveside_)
                {
                  if (kinetics_slave_cond_id == condid)
                  {
                    // collect condition specific data and store to scatra boundary parameter class
                    set_condition_specific_sca_tra_parameters(*kinetics_slave_cond);

                    scatratimint_->discretization()->evaluate_condition(
                        condparams, strategy, "S2IKinetics", condid);
                  }
                }

                // finalize auxiliary matrix block
                islavematrix->Complete(dofrowmap_growth, dofrowmap_scatra);

                // assemble linearizations of slave fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses into global matrix
                // block
                scatragrowthblock->Add(*islavematrix, false, 1., 0.);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses and assemble into
                // global matrix block
                Core::LinAlg::MatrixRowTransform()(*islavematrix, -1.,
                    Core::Adapter::CouplingSlaveConverter(*icoup_), *scatragrowthblock, false);

                // zero out auxiliary matrix block for subsequent evaluation
                islavematrix->Zero();

                // evaluate off-diagonal linearizations arising from scatra-scatra interface layer
                // growth
                auto* s2i_coupling_growth_cond =
                    scatratimint_->discretization()->GetCondition("S2IKineticsGrowth");

                // collect condition specific data and store to scatra boundary parameter class
                set_condition_specific_sca_tra_parameters(*s2i_coupling_growth_cond);

                scatratimint_->discretization()->evaluate_condition(
                    condparams, strategy, "S2IKineticsGrowth", condid);

                // finalize auxiliary matrix block
                islavematrix->Complete(dofrowmap_growth, dofrowmap_scatra);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // layer growth w.r.t. scatra-scatra interface layer thicknesses and assemble into
                // global matrix block
                Core::LinAlg::MatrixRowTransform()(*islavematrix, -1.,
                    Core::Adapter::CouplingSlaveConverter(*icoup_), *scatragrowthblock, true);

                // finalize global matrix block
                scatragrowthblock->Complete(dofrowmap_growth, dofrowmap_scatra);

                // apply Dirichlet boundary conditions to global matrix block
                scatragrowthblock->ApplyDirichlet(*scatratimint_->DirichMaps()->CondMap(), false);
              }

              // assemble off-diagonal growth-scatra block of global system matrix, containing
              // derivatives of discrete scatra-scatra interface layer growth residuals w.r.t.
              // discrete scatra degrees of freedom
              {
                // check matrix
                const Teuchos::RCP<Core::LinAlg::SparseMatrix> growthscatrablock =
                    Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(growthscatrablock_);
                if (growthscatrablock == Teuchos::null)
                  FOUR_C_THROW("Matrix is not a sparse matrix!");

                // initialize matrix block
                growthscatrablock->Zero();

                // initialize auxiliary matrix blocks for linearizations of scatra-scatra interface
                // layer growth residuals w.r.t. slave-side and master-side scalar transport degrees
                // of freedom
                Teuchos::RCP<Core::LinAlg::SparseMatrix> islavematrix =
                    Teuchos::rcp(new Core::LinAlg::SparseMatrix(dofrowmap_growth, 81));
                Teuchos::RCP<Core::LinAlg::SparseMatrix> imastermatrix =
                    Teuchos::rcp(new Core::LinAlg::SparseMatrix(dofrowmap_growth, 81));

                // initialize assembly strategy for auxiliary matrix block
                Core::FE::AssembleStrategy strategy(
                    2, 0, islavematrix, imastermatrix, Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
                    "action", ScaTra::BoundaryAction::calc_s2icoupling_growthscatra, condparams);

                // evaluate off-diagonal linearizations
                scatratimint_->discretization()->evaluate_condition(
                    condparams, strategy, "S2IKineticsGrowth", condid);

                // finalize auxiliary matrix blocks
                islavematrix->Complete(dofrowmap_scatra, dofrowmap_growth);
                imastermatrix->Complete(dofrowmap_scatra, dofrowmap_growth);

                // assemble linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // slave-side scalar transport degrees of freedom into global matrix block
                growthscatrablock->Add(*islavematrix, false, 1., 0.);

                // derive linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // master-side scalar transport degrees of freedom and assemble into global matrix
                // block
                Core::LinAlg::MatrixColTransform()(imastermatrix->RowMap(), imastermatrix->ColMap(),
                    *imastermatrix, 1., Core::Adapter::CouplingSlaveConverter(*icoup_),
                    *growthscatrablock, true, true);

                // finalize global matrix block
                growthscatrablock->Complete(dofrowmap_scatra, dofrowmap_growth);
              }

              break;
            }

            case Core::LinAlg::MatrixType::block_condition:
            case Core::LinAlg::MatrixType::block_condition_dof:
            {
              // assemble off-diagonal scatra-growth block of global system matrix, containing
              // derivatives of discrete scatra residuals w.r.t. discrete scatra-scatra interface
              // layer thicknesses
              {
                // initialize auxiliary matrix block for linearizations of slave fluxes w.r.t.
                // scatra-scatra interface layer thicknesses
                const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockslavematrix =
                    Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<
                        Core::LinAlg::DefaultBlockMatrixStrategy>(
                        *blockmapgrowth_, *blockmaps_slave_, 81, false, true));

                // initialize assembly strategy for auxiliary matrix block
                Core::FE::AssembleStrategy strategy(0, 2, blockslavematrix, Teuchos::null,
                    Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
                    "action", ScaTra::BoundaryAction::calc_s2icoupling_scatragrowth, condparams);

                // evaluate off-diagonal linearizations arising from scatra-scatra interface
                // coupling
                scatratimint_->discretization()->evaluate_condition(
                    condparams, strategy, "S2IKinetics", condid);

                // finalize auxiliary matrix block
                blockslavematrix->Complete();

                // assemble linearizations of slave fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses into global matrix
                // block
                scatragrowthblock_->Add(*blockslavematrix, false, 1., 0.);

                // initialize auxiliary system matrix for linearizations of master fluxes associated
                // with scatra-scatra interface coupling w.r.t. scatra-scatra interface layer
                // thicknesses
                Core::LinAlg::SparseMatrix mastermatrix(*icoup_->MasterDofMap(), 27, false, true);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses
                for (int iblock = 0; iblock < blockmaps_slave_->NumMaps(); ++iblock)
                  Core::LinAlg::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.,
                      Core::Adapter::CouplingSlaveConverter(*icoup_), mastermatrix, true);

                // zero out auxiliary matrices for subsequent evaluation
                blockslavematrix->Zero();
                mastermatrix.Zero();

                // evaluate off-diagonal linearizations arising from scatra-scatra interface layer
                // growth
                scatratimint_->discretization()->evaluate_condition(
                    condparams, strategy, "S2IKineticsGrowth", condid);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // layer growth w.r.t. scatra-scatra interface layer thicknesses
                for (int iblock = 0; iblock < blockmaps_slave_->NumMaps(); ++iblock)
                  Core::LinAlg::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.,
                      Core::Adapter::CouplingSlaveConverter(*icoup_), mastermatrix, true);

                // finalize auxiliary system matrix
                mastermatrix.Complete(dofrowmap_growth, *icoup_->MasterDofMap());

                // split auxiliary system matrix and assemble into global matrix block
                const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockmastermatrix =
                    mastermatrix.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                        *blockmapgrowth_, *scatratimint_->BlockMaps());
                blockmastermatrix->Complete();
                scatragrowthblock_->Add(*blockmastermatrix, false, 1., 1.);

                // finalize global matrix block
                scatragrowthblock_->Complete();

                // apply Dirichlet boundary conditions to global matrix block
                scatragrowthblock_->ApplyDirichlet(*scatratimint_->DirichMaps()->CondMap(), false);
              }

              // assemble off-diagonal growth-scatra block of global system matrix, containing
              // derivatives of discrete scatra-scatra interface layer growth residuals w.r.t.
              // discrete scatra degrees of freedom
              {
                // initialize auxiliary matrix blocks for linearizations of scatra-scatra interface
                // layer growth residuals w.r.t. slave-side and master-side scalar transport degrees
                // of freedom
                const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockslavematrix =
                    Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<
                        Core::LinAlg::DefaultBlockMatrixStrategy>(
                        *blockmaps_slave_, *blockmapgrowth_, 81, false, true));
                const Teuchos::RCP<Core::LinAlg::SparseMatrix> imastermatrix =
                    Teuchos::rcp(new Core::LinAlg::SparseMatrix(dofrowmap_growth, 81));

                // initialize assembly strategy for auxiliary matrix blocks
                Core::FE::AssembleStrategy strategy(2, 0, blockslavematrix, imastermatrix,
                    Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
                    "action", ScaTra::BoundaryAction::calc_s2icoupling_growthscatra, condparams);

                // evaluate off-diagonal linearizations
                scatratimint_->discretization()->evaluate_condition(
                    condparams, strategy, "S2IKineticsGrowth", condid);

                // finalize auxiliary matrix blocks
                blockslavematrix->Complete();
                imastermatrix->Complete(dofrowmap_scatra, dofrowmap_growth);

                // assemble linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // slave-side scalar transport degrees of freedom into global matrix block
                growthscatrablock_->Add(*blockslavematrix, false, 1., 0.);

                // initialize temporary matrix
                Core::LinAlg::SparseMatrix kgm(dofrowmap_growth, 27, false, true);

                // derive linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // master-side scalar transport degrees of freedom
                Core::LinAlg::MatrixColTransform()(imastermatrix->RowMap(), imastermatrix->ColMap(),
                    *imastermatrix, 1., Core::Adapter::CouplingSlaveConverter(*icoup_), kgm);

                // finalize temporary matrix
                kgm.Complete(*icoup_->MasterDofMap(), dofrowmap_growth);

                // split temporary matrix and assemble into global matrix block
                const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> blockkgm(
                    kgm.Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
                        *scatratimint_->BlockMaps(), *blockmapgrowth_));
                blockkgm->Complete();
                growthscatrablock_->Add(*blockkgm, false, 1., 1.);

                // finalize global matrix block
                growthscatrablock_->Complete();
              }

              break;
            }

            default:
            {
              FOUR_C_THROW(
                  "Type of global system matrix for scatra-scatra interface coupling involving "
                  "interface layer growth not recognized!");
              break;
            }
          }  // type of scalar transport system matrix

          // assemble residual vector associated with scatra-scatra interface layer thicknesses and
          // main-diagonal growth-growth block of global system matrix, containing derivatives of
          // discrete scatra-scatra interface layer growth residuals w.r.t. discrete scatra-scatra
          // interface layer thicknesses
          {
            // initialize matrix block and corresponding residual vector
            growthgrowthblock_->Zero();
            growthresidual_->PutScalar(0.);

            // initialize assembly strategy for main-diagonal growth-growth block and
            Core::FE::AssembleStrategy strategy(2, 2, growthgrowthblock_, Teuchos::null,
                growthresidual_, Teuchos::null, Teuchos::null);

            // create parameter list for elements
            Teuchos::ParameterList condparams;

            // set action for elements
            Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
                "action", ScaTra::BoundaryAction::calc_s2icoupling_growthgrowth, condparams);

            // set history vector associated with discrete scatra-scatra interface layer thicknesses
            scatratimint_->discretization()->set_state(2, "growthhist", growthhist_);

            // evaluate main-diagonal linearizations and corresponding residuals
            scatratimint_->discretization()->evaluate_condition(
                condparams, strategy, "S2IKineticsGrowth", condid);

            // finalize global matrix block
            growthgrowthblock_->Complete();
          }
        }  // monolithic evaluation of scatra-scatra interface layer growth

        break;
      }

      default:
      {
        FOUR_C_THROW(
            "Evaluation of scatra-scatra interface layer growth only implemented for conforming "
            "interface discretizations!");
        break;
      }
    }
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_and_assemble_capacitive_contributions()
{
  // create parameter list for elements
  Teuchos::ParameterList capcondparas;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_s2icoupling_capacitance, capcondparas);

  // set global state vectors according to time-integration scheme
  scatratimint_->add_time_integration_specific_vectors();

  // zero out matrices and vectors
  islavematrix_->Zero();
  imasterslavematrix_->Zero();
  islaveresidual_->PutScalar(0.0);
  auto imasterresidual_on_slave_side = Teuchos::rcp(new Epetra_Vector(*interfacemaps_->Map(1)));
  imasterresidual_on_slave_side->PutScalar(0.0);

  // evaluate scatra-scatra interface coupling
  for (auto kinetics_slave_cond_cap : kinetics_conditions_meshtying_slaveside_)
  {
    if (kinetics_slave_cond_cap.second->parameters().Get<int>("kinetic model") ==
        static_cast<int>(Inpar::S2I::kinetics_butlervolmerreducedcapacitance))
    {
      // collect condition specific data and store to scatra boundary parameter class
      set_condition_specific_sca_tra_parameters(*kinetics_slave_cond_cap.second);

      scatratimint_->discretization()->evaluate_condition(capcondparas, islavematrix_,
          imasterslavematrix_, islaveresidual_, imasterresidual_on_slave_side, Teuchos::null,
          "S2IKinetics", kinetics_slave_cond_cap.second->parameters().Get<int>("ConditionID"));
    }
  }

  // finalize interface matrices
  islavematrix_->Complete();
  imasterslavematrix_->Complete();

  switch (matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
    {
      auto systemmatrix = scatratimint_->SystemMatrix();
      FOUR_C_ASSERT(systemmatrix != Teuchos::null, "System matrix is not a sparse matrix!");

      // assemble additional components of linearizations of slave fluxes due to capacitance
      // w.r.t. slave dofs into the global system matrix
      systemmatrix->Add(*islavematrix_, false, 1.0, 1.0);

      // assemble additional components of linearizations of slave fluxes due to capacitance
      // w.r.t. master dofs into the global system matrix
      Core::LinAlg::MatrixColTransform()(islavematrix_->RowMap(), islavematrix_->ColMap(),
          *islavematrix_, -1.0, Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true,
          true);

      // assemble additional components of linearizations of master fluxes due to capacitance
      // w.r.t. slave dofs into the global system matrix
      Core::LinAlg::MatrixRowTransform()(*imasterslavematrix_, 1.0,
          Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true);

      // assemble additional components of linearizations of master fluxes due to capacitance
      // w.r.t. master dofs into the global system matrix
      Core::LinAlg::MatrixRowColTransform()(*imasterslavematrix_, -1.0,
          Core::Adapter::CouplingSlaveConverter(*icoup_),
          Core::Adapter::CouplingSlaveConverter(*icoup_), *systemmatrix, true, true);
      break;
    }
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      // check matrix
      auto blocksystemmatrix = scatratimint_->BlockSystemMatrix();
      FOUR_C_ASSERT(blocksystemmatrix != Teuchos::null, "System matrix is not a block matrix!");

      // prepare linearizations of slave fluxes due to capacitance w.r.t. slave dofs
      auto blockkss = islavematrix_->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *blockmaps_slave_, *blockmaps_slave_);
      blockkss->Complete();

      // prepare linearizations of slave fluxes due to capacitance w.r.t. master dofs
      auto ksm = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->SlaveDofMap(), 81, false));
      Core::LinAlg::MatrixColTransform()(islavematrix_->RowMap(), islavematrix_->ColMap(),
          *islavematrix_, -1.0, Core::Adapter::CouplingSlaveConverter(*icoup_), *ksm);
      ksm->Complete(*icoup_->MasterDofMap(), *icoup_->SlaveDofMap());
      auto blockksm = ksm->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *blockmaps_master_, *blockmaps_slave_);
      blockksm->Complete();

      // prepare linearizations of master fluxes due to capacitance w.r.t. slave dofs
      auto kms = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false));
      Core::LinAlg::MatrixRowTransform()(
          *imasterslavematrix_, 1.0, Core::Adapter::CouplingSlaveConverter(*icoup_), *kms);
      kms->Complete(*icoup_->SlaveDofMap(), *icoup_->MasterDofMap());
      auto blockkms = kms->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *blockmaps_slave_, *blockmaps_master_);
      blockkms->Complete();

      // derive linearizations of master fluxes w.r.t. master dofs
      auto kmm = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*icoup_->MasterDofMap(), 81, false));
      Core::LinAlg::MatrixRowColTransform()(*imasterslavematrix_, -1.0,
          Core::Adapter::CouplingSlaveConverter(*icoup_),
          Core::Adapter::CouplingSlaveConverter(*icoup_), *kmm);
      kmm->Complete();
      auto blockkmm = kmm->Split<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *blockmaps_master_, *blockmaps_master_);
      blockkmm->Complete();

      // assemble interface block matrices into global block system matrix
      blocksystemmatrix->Add(*blockkss, false, 1.0, 1.0);
      blocksystemmatrix->Add(*blockksm, false, 1.0, 1.0);
      blocksystemmatrix->Add(*blockkms, false, 1.0, 1.0);
      blocksystemmatrix->Add(*blockkmm, false, 1.0, 1.0);

      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Type of global system matrix for scatra-scatra interface coupling not recognized!");
      break;
    }
  }

  // assemble slave residuals into global residual vector
  interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
  // transform master residuals and assemble into global residual vector
  interfacemaps_->AddVector(
      icoup_->SlaveToMaster(imasterresidual_on_slave_side), 2, scatratimint_->Residual(), 1.0);
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_mortar_cell(const Core::FE::Discretization& idiscret,
    Mortar::IntCell& cell, const Inpar::ScaTra::ImplType& impltype, Mortar::Element& slaveelement,
    Mortar::Element& masterelement, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& cellmatrix1, Core::LinAlg::SerialDenseMatrix& cellmatrix2,
    Core::LinAlg::SerialDenseMatrix& cellmatrix3, Core::LinAlg::SerialDenseMatrix& cellmatrix4,
    Core::LinAlg::SerialDenseVector& cellvector1,
    Core::LinAlg::SerialDenseVector& cellvector2) const
{
  // evaluate single mortar integration cell
  ScaTra::MortarCellFactory::mortar_cell_calc(
      impltype, slaveelement, masterelement, couplingtype_, lmside_, idiscret.Name())
      ->Evaluate(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_slave_node(const Core::FE::Discretization& idiscret,
    const Mortar::Node& slavenode, const double& lumpedarea,
    const Inpar::ScaTra::ImplType& impltype, Mortar::Element& slaveelement,
    Mortar::Element& masterelement, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& ntsmatrix1, Core::LinAlg::SerialDenseMatrix& ntsmatrix2,
    Core::LinAlg::SerialDenseMatrix& ntsmatrix3, Core::LinAlg::SerialDenseMatrix& ntsmatrix4,
    Core::LinAlg::SerialDenseVector& ntsvector1, Core::LinAlg::SerialDenseVector& ntsvector2) const
{
  // evaluate single slave-side node
  ScaTra::MortarCellFactory::mortar_cell_calc(
      impltype, slaveelement, masterelement, couplingtype_, lmside_, idiscret.Name())
      ->evaluate_nts(idiscret, slavenode, lumpedarea, slaveelement, masterelement, la_slave,
          la_master, params, ntsmatrix1, ntsmatrix2, ntsmatrix3, ntsmatrix4, ntsvector1,
          ntsvector2);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_mortar_element(const Core::FE::Discretization& idiscret,
    Mortar::Element& element, const Inpar::ScaTra::ImplType& impltype,
    Core::Elements::Element::LocationArray& la, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& elematrix1, Core::LinAlg::SerialDenseMatrix& elematrix2,
    Core::LinAlg::SerialDenseMatrix& elematrix3, Core::LinAlg::SerialDenseMatrix& elematrix4,
    Core::LinAlg::SerialDenseVector& elevector1, Core::LinAlg::SerialDenseVector& elevector2) const
{
  // evaluate single mortar element
  ScaTra::MortarCellFactory::mortar_cell_calc(
      impltype, element, element, couplingtype_, lmside_, idiscret.Name())
      ->evaluate_mortar_element(idiscret, element, la, params, elematrix1, elematrix2, elematrix3,
          elematrix4, elevector1, elevector2);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_mortar_cells(const Core::FE::Discretization& idiscret,
    const Teuchos::ParameterList& params,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix1,
    const Inpar::S2I::InterfaceSides matrix1_side_rows,
    const Inpar::S2I::InterfaceSides matrix1_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix2,
    const Inpar::S2I::InterfaceSides matrix2_side_rows,
    const Inpar::S2I::InterfaceSides matrix2_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix3,
    const Inpar::S2I::InterfaceSides matrix3_side_rows,
    const Inpar::S2I::InterfaceSides matrix3_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix4,
    const Inpar::S2I::InterfaceSides matrix4_side_rows,
    const Inpar::S2I::InterfaceSides matrix4_side_cols,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
    const Inpar::S2I::InterfaceSides vector1_side,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
    const Inpar::S2I::InterfaceSides vector2_side) const
{
  // instantiate assembly strategy for mortar integration cells
  ScaTra::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // evaluate mortar integration cells
  evaluate_mortar_cells(idiscret, params, strategy);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_mortar_cells(const Core::FE::Discretization& idiscret,
    const Teuchos::ParameterList& params, ScaTra::MortarCellAssemblyStrategy& strategy) const
{
  // extract scatra-scatra interface coupling condition from parameter list
  const Core::Conditions::Condition* const condition =
      params.get<Core::Conditions::Condition*>("condition");
  if (condition == nullptr)
    FOUR_C_THROW("Cannot access scatra-scatra interface coupling condition!");

  // extract mortar integration cells associated with current condition
  const std::vector<std::pair<Teuchos::RCP<Mortar::IntCell>, Inpar::ScaTra::ImplType>>& cells =
      imortarcells_.at(condition->parameters().Get<int>("ConditionID"));

  // loop over all mortar integration cells
  for (const auto& icell : cells)
  {
    // extract current mortar integration cell
    const Teuchos::RCP<Mortar::IntCell>& cell = icell.first;
    if (cell == Teuchos::null) FOUR_C_THROW("Invalid mortar integration cell!");

    // extract slave-side element associated with current cell
    auto* slaveelement = dynamic_cast<Mortar::Element*>(idiscret.gElement(cell->GetSlaveId()));
    if (!slaveelement)
      FOUR_C_THROW("Couldn't extract slave element from mortar interface discretization!");

    // extract master-side element associated with current cell
    auto* masterelement = dynamic_cast<Mortar::Element*>(idiscret.gElement(cell->GetMasterId()));
    if (!masterelement)
      FOUR_C_THROW("Couldn't extract master element from mortar interface discretization!");

    // safety check
    if (!slaveelement->IsSlave() or masterelement->IsSlave())
      FOUR_C_THROW("Something is wrong with the slave-master element pairing!");

    // construct slave-side and master-side location arrays
    Core::Elements::Element::LocationArray la_slave(idiscret.NumDofSets());
    slaveelement->LocationVector(idiscret, la_slave, false);
    Core::Elements::Element::LocationArray la_master(idiscret.NumDofSets());
    masterelement->LocationVector(idiscret, la_master, false);

    // initialize cell matrices and vectors
    strategy.init_cell_matrices_and_vectors(la_slave, la_master);

    // evaluate current cell
    evaluate_mortar_cell(idiscret, *cell, icell.second, *slaveelement, *masterelement, la_slave,
        la_master, params, strategy.CellMatrix1(), strategy.CellMatrix2(), strategy.CellMatrix3(),
        strategy.CellMatrix4(), strategy.CellVector1(), strategy.CellVector2());

    // assemble cell matrices and vectors into system matrices and vectors
    strategy.assemble_cell_matrices_and_vectors(la_slave, la_master, la_slave[0].lmowner_[0]);
  }
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_nts(const Epetra_IntVector& islavenodestomasterelements,
    const Epetra_Vector& islavenodeslumpedareas, const Epetra_IntVector& islavenodesimpltypes,
    const Core::FE::Discretization& idiscret, const Teuchos::ParameterList& params,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix1,
    const Inpar::S2I::InterfaceSides matrix1_side_rows,
    const Inpar::S2I::InterfaceSides matrix1_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix2,
    const Inpar::S2I::InterfaceSides matrix2_side_rows,
    const Inpar::S2I::InterfaceSides matrix2_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix3,
    const Inpar::S2I::InterfaceSides matrix3_side_rows,
    const Inpar::S2I::InterfaceSides matrix3_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix4,
    const Inpar::S2I::InterfaceSides matrix4_side_rows,
    const Inpar::S2I::InterfaceSides matrix4_side_cols,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
    const Inpar::S2I::InterfaceSides vector1_side,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
    const Inpar::S2I::InterfaceSides vector2_side) const
{
  // instantiate assembly strategy for node-to-segment coupling
  ScaTra::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // extract slave-side noderowmap
  const Epetra_BlockMap& noderowmap_slave = islavenodestomasterelements.Map();

  // loop over all slave-side nodes
  for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
  {
    // extract slave-side node
    auto* const slavenode =
        dynamic_cast<Mortar::Node* const>(idiscret.gNode(noderowmap_slave.GID(inode)));
    if (slavenode == nullptr) FOUR_C_THROW("Couldn't extract slave-side node from discretization!");

    // extract first slave-side element associated with current slave-side node
    auto* const slaveelement = dynamic_cast<Mortar::Element* const>(slavenode->Elements()[0]);
    if (!slaveelement) FOUR_C_THROW("Invalid slave-side mortar element!");

    // extract master-side element associated with current slave-side node
    auto* const masterelement =
        dynamic_cast<Mortar::Element* const>(idiscret.gElement(islavenodestomasterelements[inode]));
    if (!masterelement) FOUR_C_THROW("Invalid master-side mortar element!");

    // safety check
    if (!slaveelement->IsSlave() or masterelement->IsSlave())
      FOUR_C_THROW("Something is wrong with the slave-master element pairing!");

    // construct slave-side and master-side location arrays
    Core::Elements::Element::LocationArray la_slave(idiscret.NumDofSets());
    slaveelement->LocationVector(idiscret, la_slave, false);
    Core::Elements::Element::LocationArray la_master(idiscret.NumDofSets());
    masterelement->LocationVector(idiscret, la_master, false);

    // initialize cell matrices and vectors
    strategy.init_cell_matrices_and_vectors(la_slave, la_master);

    // evaluate current slave-side node
    evaluate_slave_node(idiscret, *slavenode, islavenodeslumpedareas[inode],
        (Inpar::ScaTra::ImplType)islavenodesimpltypes[inode], *slaveelement, *masterelement,
        la_slave, la_master, params, strategy.CellMatrix1(), strategy.CellMatrix2(),
        strategy.CellMatrix3(), strategy.CellMatrix4(), strategy.CellVector1(),
        strategy.CellVector2());

    // assemble cell matrices and vectors into system matrices and vectors
    strategy.assemble_cell_matrices_and_vectors(la_slave, la_master, slavenode->Owner());
  }
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::evaluate_mortar_elements(const Epetra_Map& ielecolmap,
    const Epetra_IntVector& ieleimpltypes, const Core::FE::Discretization& idiscret,
    const Teuchos::ParameterList& params,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix1,
    const Inpar::S2I::InterfaceSides matrix1_side_rows,
    const Inpar::S2I::InterfaceSides matrix1_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix2,
    const Inpar::S2I::InterfaceSides matrix2_side_rows,
    const Inpar::S2I::InterfaceSides matrix2_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix3,
    const Inpar::S2I::InterfaceSides matrix3_side_rows,
    const Inpar::S2I::InterfaceSides matrix3_side_cols,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix4,
    const Inpar::S2I::InterfaceSides matrix4_side_rows,
    const Inpar::S2I::InterfaceSides matrix4_side_cols,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,
    const Inpar::S2I::InterfaceSides vector1_side,
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,
    const Inpar::S2I::InterfaceSides vector2_side) const
{
  // instantiate assembly strategy for mortar elements
  ScaTra::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // loop over all mortar elements
  for (int ielement = 0; ielement < ielecolmap.NumMyElements(); ++ielement)
  {
    // extract current mortar element
    auto* const element =
        dynamic_cast<Mortar::Element* const>(idiscret.gElement(ielecolmap.GID(ielement)));
    if (!element) FOUR_C_THROW("Couldn't extract mortar element from mortar discretization!");

    // construct location array for current mortar element
    Core::Elements::Element::LocationArray la(idiscret.NumDofSets());
    element->LocationVector(idiscret, la, false);

    // initialize element matrices and vectors
    strategy.init_cell_matrices_and_vectors(
        la, la);  // second function argument only serves as dummy

    // evaluate current mortar element
    evaluate_mortar_element(idiscret, *element, (Inpar::ScaTra::ImplType)ieleimpltypes[ielement],
        la, params, strategy.CellMatrix1(), strategy.CellMatrix2(), strategy.CellMatrix3(),
        strategy.CellMatrix4(), strategy.CellVector1(), strategy.CellVector2());

    // assemble element matrices and vectors into system matrices and vectors
    strategy.assemble_cell_matrices_and_vectors(la, la, -1);
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
ScaTra::MortarCellInterface* ScaTra::MortarCellFactory::mortar_cell_calc(
    const Inpar::ScaTra::ImplType& impltype, const Mortar::Element& slaveelement,
    const Mortar::Element& masterelement, const Inpar::S2I::CouplingType& couplingtype,
    const Inpar::S2I::InterfaceSides& lmside, const std::string& disname)
{
  // extract number of slave-side degrees of freedom per node
  const int numdofpernode_slave = slaveelement.NumDofPerNode(*slaveelement.Nodes()[0]);

  switch (slaveelement.Shape())
  {
    case Core::FE::CellType::tri3:
    {
      return mortar_cell_calc<Core::FE::CellType::tri3>(
          impltype, masterelement, couplingtype, lmside, numdofpernode_slave, disname);
      break;
    }

    case Core::FE::CellType::quad4:
    {
      return mortar_cell_calc<Core::FE::CellType::quad4>(
          impltype, masterelement, couplingtype, lmside, numdofpernode_slave, disname);
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid slave-side discretization type!");
      break;
    }
  }

  return nullptr;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
template <Core::FE::CellType distypeS>
ScaTra::MortarCellInterface* ScaTra::MortarCellFactory::mortar_cell_calc(
    const Inpar::ScaTra::ImplType& impltype, const Mortar::Element& masterelement,
    const Inpar::S2I::CouplingType& couplingtype, const Inpar::S2I::InterfaceSides& lmside,
    const int& numdofpernode_slave, const std::string& disname)
{
  // extract number of master-side degrees of freedom per node
  const int numdofpernode_master = masterelement.NumDofPerNode(*masterelement.Nodes()[0]);

  switch (masterelement.Shape())
  {
    case Core::FE::CellType::tri3:
    {
      return mortar_cell_calc<distypeS, Core::FE::CellType::tri3>(
          impltype, couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case Core::FE::CellType::quad4:
    {
      return mortar_cell_calc<distypeS, Core::FE::CellType::quad4>(
          impltype, couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid master-side discretization type!");
      break;
    }
  }

  return nullptr;
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
ScaTra::MortarCellInterface* ScaTra::MortarCellFactory::mortar_cell_calc(
    const Inpar::ScaTra::ImplType& impltype, const Inpar::S2I::CouplingType& couplingtype,
    const Inpar::S2I::InterfaceSides& lmside, const int& numdofpernode_slave,
    const int& numdofpernode_master, const std::string& disname)
{
  // return instance of evaluation class for mortar integration cell depending on physical
  // implementation type
  switch (impltype)
  {
    case Inpar::ScaTra::impltype_std:
    {
      return ScaTra::MortarCellCalc<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case Inpar::ScaTra::impltype_elch_electrode:
    {
      return ScaTra::MortarCellCalcElch<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case Inpar::ScaTra::impltype_elch_electrode_thermo:
    {
      return ScaTra::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case Inpar::ScaTra::impltype_thermo_elch_electrode:
    {
      return ScaTra::MortarCellCalcSTIElch<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown physical implementation type of mortar integration cell!");
      break;
    }
  }

  return nullptr;
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::init_conv_check_strategy()
{
  if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_petrov or
      couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov)
  {
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyS2ILM(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
  else
    convcheckstrategy_ = Teuchos::rcp(new ScaTra::ConvCheckStrategyStd(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
}  // ScaTra::meshtying_strategy_s2_i::init_conv_check_strategy


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::setup_meshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<Core::Conditions::Condition*> s2imeshtying_conditions(0, nullptr);
  scatratimint_->discretization()->GetCondition("S2IMeshtying", s2imeshtying_conditions);
  std::vector<Core::Conditions::Condition*> s2ikinetics_conditions(0, nullptr);
  scatratimint_->discretization()->GetCondition("S2IKinetics", s2ikinetics_conditions);
  kinetics_conditions_meshtying_slaveside_.clear();
  master_conditions_.clear();
  runtime_csvwriter_.emplace(scatratimint_->discretization()->Comm().MyPID(),
      *scatratimint_->DiscWriter()->Output(), "kinetics_interface_flux");

  for (auto* s2imeshtying_cond : s2imeshtying_conditions)
  {
    for (auto* s2ikinetics_cond : s2ikinetics_conditions)
    {
      const int s2ikinetics_cond_id = s2ikinetics_cond->parameters().Get<int>("ConditionID");
      const int s2ikinetics_cond_interface_side =
          s2ikinetics_cond->parameters().Get<int>("interface side");

      if (s2ikinetics_cond_id < 0)
        FOUR_C_THROW("Invalid condition ID %i for S2IKinetics Condition!", s2ikinetics_cond_id);

      // only continue if ID's match
      if (s2imeshtying_cond->parameters().Get<int>("S2IKineticsID") != s2ikinetics_cond_id)
        continue;
      // only continue if sides match
      if (s2imeshtying_cond->parameters().Get<int>("interface side") !=
          s2ikinetics_cond_interface_side)
        continue;

      switch (s2ikinetics_cond_interface_side)
      {
        case Inpar::S2I::side_slave:
        {
          if (kinetics_conditions_meshtying_slaveside_.find(s2ikinetics_cond_id) ==
              kinetics_conditions_meshtying_slaveside_.end())
          {
            kinetics_conditions_meshtying_slaveside_.insert(
                std::make_pair(s2ikinetics_cond_id, s2ikinetics_cond));
          }
          else
          {
            FOUR_C_THROW(
                "Cannot have multiple slave-side scatra-scatra interface kinetics conditions with "
                "the same ID %i!",
                s2ikinetics_cond_id);
          }

          if (s2ikinetics_cond->parameters().Get<int>("kinetic model") ==
              static_cast<int>(Inpar::S2I::kinetics_butlervolmerreducedcapacitance))
          {
            has_capacitive_contributions_ = true;

            auto timeintscheme = Core::UTILS::IntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(
                *scatratimint_->ScatraParameterList(), "TIMEINTEGR");
            if (not(timeintscheme == Inpar::ScaTra::timeint_bdf2 or
                    timeintscheme == Inpar::ScaTra::timeint_one_step_theta))
            {
              FOUR_C_THROW(
                  "Solution of capacitive interface contributions, i.e. additional transient terms "
                  "is only implemented for OST and BDF2 time integration schemes.");
            }
          }
          std::stringstream condition_flux_name, condition_area_name;
          condition_flux_name << "flux_S2I_condition_" << s2ikinetics_cond_id;
          condition_area_name << "area_S2I_condition_" << s2ikinetics_cond_id;
          runtime_csvwriter_->register_data_vector(condition_flux_name.str(), 1, 16);
          runtime_csvwriter_->register_data_vector(condition_area_name.str(), 1, 16);

          break;
        }

        case Inpar::S2I::side_master:
        {
          if (master_conditions_.find(s2ikinetics_cond_id) == master_conditions_.end())
          {
            master_conditions_.insert(std::make_pair(s2ikinetics_cond_id, s2ikinetics_cond));
          }
          else
          {
            FOUR_C_THROW(
                "Cannot have multiple master-side scatra-scatra interface kinetics conditions with "
                "the same ID %i!",
                s2ikinetics_cond_id);
          }
          break;
        }

        default:
        {
          FOUR_C_THROW("Invalid scatra-scatra interface kinetics condition!");
          break;
        }
      }
    }
  }

  // determine type of mortar meshtying
  switch (couplingtype_)
  {
    // setup scatra-scatra interface coupling for interfaces with pairwise overlapping interface
    // nodes
    case Inpar::S2I::coupling_matching_nodes:
    {
      // overwrite IDs of master-side scatra-scatra interface coupling conditions with the value -1
      // to prevent them from being evaluated when calling evaluate_condition on the discretization
      // TODO: this is somewhat unclean, because changing the conditions, makes calling
      // setup_meshtying() twice invalid (which should not be necessary, but conceptually possible)
      for (auto& mastercondition : master_conditions_)
        mastercondition.second->parameters().Add("ConditionID", -1);

      if (scatratimint_->NumScal() < 1)
        FOUR_C_THROW("Number of transported scalars not correctly set!");

      // construct new (empty coupling adapter)
      icoup_ = Teuchos::rcp(new Core::Adapter::Coupling());

      int num_dof_per_condition = -1;

      // initialize int vectors for global ids of slave and master interface nodes
      if (indepedent_setup_of_conditions_)
      {
        std::vector<std::vector<int>> islavenodegidvec_cond;
        std::vector<std::vector<int>> imasternodegidvec_cond;

        for (auto& kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
        {
          std::vector<int> islavenodegidvec;
          std::vector<int> imasternodegidvec;

          const int kineticsID = kinetics_slave_cond.first;
          auto* kinetics_condition = kinetics_slave_cond.second;

          if (num_dof_per_condition == -1)
            num_dof_per_condition =
                scatratimint_->num_dof_per_node_in_condition(*kinetics_condition);
          else if (num_dof_per_condition !=
                   scatratimint_->num_dof_per_node_in_condition(*kinetics_condition))
            FOUR_C_THROW("all S2I conditions must have the same number of dof per node");

          if (kinetics_condition->parameters().Get<int>("kinetic model") !=
              static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux))
          {
            Core::Communication::AddOwnedNodeGIDFromList(*scatratimint_->discretization(),
                *kinetics_condition->GetNodes(), islavenodegidvec);

            auto mastercondition = master_conditions_.find(kineticsID);
            if (mastercondition == master_conditions_.end())
              FOUR_C_THROW("Could not find master condition");

            Core::Communication::AddOwnedNodeGIDFromList(*scatratimint_->discretization(),
                *mastercondition->second->GetNodes(), imasternodegidvec);

            islavenodegidvec_cond.push_back(islavenodegidvec);
            imasternodegidvec_cond.push_back(imasternodegidvec);
          }
        }

        icoup_->setup_coupling(*(scatratimint_->discretization()),
            *(scatratimint_->discretization()), imasternodegidvec_cond, islavenodegidvec_cond,
            num_dof_per_condition, true, 1.0e-8);
      }
      else
      {
        std::set<int> islavenodegidset;
        std::set<int> imasternodegidset;

        for (const auto& kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
        {
          const int kineticsID = kinetics_slave_cond.first;
          auto* kinetics_condition = kinetics_slave_cond.second;

          if (num_dof_per_condition == -1)
            num_dof_per_condition =
                scatratimint_->num_dof_per_node_in_condition(*kinetics_condition);
          else if (num_dof_per_condition !=
                   scatratimint_->num_dof_per_node_in_condition(*kinetics_condition))
            FOUR_C_THROW("all S2I conditions must have the same number of dof per node");

          if (kinetics_condition->parameters().Get<int>("kinetic model") !=
              static_cast<int>(Inpar::S2I::kinetics_nointerfaceflux))
          {
            Core::Communication::AddOwnedNodeGIDFromList(*scatratimint_->discretization(),
                *kinetics_condition->GetNodes(), islavenodegidset);

            auto mastercondition = master_conditions_.find(kineticsID);
            if (mastercondition == master_conditions_.end())
              FOUR_C_THROW("Could not find master condition");
            else
              Core::Communication::AddOwnedNodeGIDFromList(*scatratimint_->discretization(),
                  *mastercondition->second->GetNodes(), imasternodegidset);
          }
        }

        std::vector<int> islavenodegidvec(islavenodegidset.begin(), islavenodegidset.end());
        std::vector<int> imasternodegidvec(imasternodegidset.begin(), imasternodegidset.end());

        icoup_->setup_coupling(*(scatratimint_->discretization()),
            *(scatratimint_->discretization()), imasternodegidvec, islavenodegidvec,
            num_dof_per_condition, true, 1.0e-8);
      }

      // generate interior and interface maps
      auto ifullmap = Core::LinAlg::MergeMap(icoup_->SlaveDofMap(), icoup_->MasterDofMap());
      std::vector<Teuchos::RCP<const Epetra_Map>> imaps;
      imaps.emplace_back(
          Core::LinAlg::SplitMap(*(scatratimint_->discretization()->dof_row_map()), *ifullmap));
      imaps.emplace_back(icoup_->SlaveDofMap());
      imaps.emplace_back(icoup_->MasterDofMap());

      // initialize global map extractor
      interfacemaps_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
          *(scatratimint_->discretization()->dof_row_map()), imaps));
      interfacemaps_->check_for_valid_map_extractor();

      // initialize interface vector
      // Although the interface vector only contains the transformed master interface dofs, we still
      // initialize it with the full dof_row_map of the discretization to make it work for parallel
      // computations.
      islavephidtnp_ =
          Core::LinAlg::CreateVector(*(scatratimint_->discretization()->dof_row_map()), false);
      imasterphidt_on_slave_side_np_ =
          Core::LinAlg::CreateVector(*(scatratimint_->discretization()->dof_row_map()), false);
      imasterphi_on_slave_side_np_ =
          Core::LinAlg::CreateVector(*(scatratimint_->discretization()->dof_row_map()), false);

      // initialize auxiliary system matrices and associated transformation operators
      islavematrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(icoup_->SlaveDofMap()), 81));
      imastermatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(icoup_->SlaveDofMap()), 81));
      imasterslavematrix_ =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(icoup_->SlaveDofMap()), 81));
      islavetomasterrowtransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowTransform);
      if (not slaveonly_)
      {
        islavetomastercoltransform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
        islavetomasterrowcoltransform_ = Teuchos::rcp(new Core::LinAlg::MatrixRowColTransform);
      }

      // initialize auxiliary residual vector
      islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*(icoup_->SlaveDofMap())));

      break;
    }

    // setup scatra-scatra interface coupling for interfaces with non-overlapping interface nodes
    case Inpar::S2I::coupling_mortar_standard:
    case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
    case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
    case Inpar::S2I::coupling_mortar_condensed_petrov:
    case Inpar::S2I::coupling_mortar_condensed_bubnov:
    case Inpar::S2I::coupling_nts_standard:
    {
      // safety checks
      if (imortarredistribution_ and couplingtype_ != Inpar::S2I::coupling_mortar_standard)
      {
        FOUR_C_THROW(
            "Parallel redistribution only implemented for scatra-scatra interface coupling based "
            "on standard mortar approach!");
      }
      if (Core::UTILS::IntegralValue<Inpar::Mortar::MeshRelocation>(
              Global::Problem::Instance()->mortar_coupling_params(), "MESH_RELOCATION") !=
          Inpar::Mortar::relocation_none)
        FOUR_C_THROW("Mesh relocation not yet implemented for scatra-scatra interface coupling!");

      // initialize empty interface maps
      Teuchos::RCP<Epetra_Map> imastermap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->discretization()->Comm()));
      Teuchos::RCP<Epetra_Map> islavemap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->discretization()->Comm()));
      Teuchos::RCP<Epetra_Map> ifullmap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->discretization()->Comm()));
      if (imortarredistribution_)
      {
        imastermap_ = Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->discretization()->Comm()));
        islavemap_ = Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->discretization()->Comm()));
      }

      // loop over all slave-side scatra-scatra interface coupling conditions
      for (auto& kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
      {
        // extract condition ID
        const int condid = kinetics_slave_cond.first;

        // initialize maps for row nodes associated with current condition
        std::map<int, Core::Nodes::Node*> masternodes;
        std::map<int, Core::Nodes::Node*> slavenodes;

        // initialize maps for column nodes associated with current condition
        std::map<int, Core::Nodes::Node*> mastergnodes;
        std::map<int, Core::Nodes::Node*> slavegnodes;

        // initialize maps for elements associated with current condition
        std::map<int, Teuchos::RCP<Core::Elements::Element>> masterelements;
        std::map<int, Teuchos::RCP<Core::Elements::Element>> slaveelements;

        // extract current slave-side and associated master-side scatra-scatra interface coupling
        // conditions
        std::vector<Core::Conditions::Condition*> mastercondition(1, master_conditions_[condid]);
        std::vector<Core::Conditions::Condition*> slavecondition(1, kinetics_slave_cond.second);

        // fill maps
        Core::Conditions::FindConditionObjects(*scatratimint_->discretization(), masternodes,
            mastergnodes, masterelements, mastercondition);
        Core::Conditions::FindConditionObjects(*scatratimint_->discretization(), slavenodes,
            slavegnodes, slaveelements, slavecondition);

        // initialize mortar coupling adapter
        icoupmortar_[condid] =
            Teuchos::rcp(new Core::Adapter::CouplingMortar(Global::Problem::Instance()->NDim(),
                Global::Problem::Instance()->mortar_coupling_params(),
                Global::Problem::Instance()->contact_dynamic_params(),
                Global::Problem::Instance()->spatial_approximation_type()));
        Core::Adapter::CouplingMortar& icoupmortar = *icoupmortar_[condid];
        std::vector<int> coupleddof(scatratimint_->NumDofPerNode(), 1);
        icoupmortar.setup_interface(scatratimint_->discretization(),
            scatratimint_->discretization(), coupleddof, mastergnodes, slavegnodes, masterelements,
            slaveelements, scatratimint_->discretization()->Comm());

        // extract mortar interface
        Mortar::Interface& interface = *icoupmortar.Interface();

        // extract mortar discretization
        const Core::FE::Discretization& idiscret = interface.Discret();

        if (couplingtype_ != Inpar::S2I::coupling_nts_standard)
        {
          // provide each slave-side mortar element with material of corresponding parent element
          for (int iele = 0; iele < interface.SlaveColElements()->NumMyElements(); ++iele)
          {
            // determine global ID of current slave-side mortar element
            const int elegid = interface.SlaveColElements()->GID(iele);

            // add material
            idiscret.gElement(elegid)->SetMaterial(
                0, Mat::Factory(Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(
                       kinetics_slave_cond.second->Geometry()[elegid])
                                    ->parent_element()
                                    ->Material()
                                    ->Parameter()
                                    ->Id()));
          }

          // assign physical implementation type to each slave-side mortar element by copying the
          // physical implementation type of the corresponding parent volume element
          Epetra_IntVector impltypes_row(*interface.SlaveRowElements());
          for (int iele = 0; iele < interface.SlaveRowElements()->NumMyElements(); ++iele)
          {
            impltypes_row[iele] = dynamic_cast<const Discret::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<const Core::Elements::FaceElement>(
                    kinetics_slave_cond.second->Geometry()[interface.SlaveRowElements()->GID(iele)])
                    ->parent_element())
                                      ->ImplType();
          }

          // perform parallel redistribution if desired
          if (imortarredistribution_ and idiscret.Comm().NumProc() > 1)
          {
            interface.interface_params()
                .sublist("PARALLEL REDISTRIBUTION")
                .set<std::string>("PARALLEL_REDIST", Global::Problem::Instance()
                                                         ->mortar_coupling_params()
                                                         .sublist("PARALLEL REDISTRIBUTION")
                                                         .get<std::string>("PARALLEL_REDIST"));
            interface.Redistribute();
            interface.fill_complete(true);
            interface.print_parallel_distribution();
            interface.CreateSearchTree();
          }

          // generate mortar integration cells
          std::vector<Teuchos::RCP<Mortar::IntCell>> imortarcells(0, Teuchos::null);
          icoupmortar.EvaluateGeometry(imortarcells);

          // assign physical implementation type to each mortar integration cell by copying the
          // physical implementation type of the corresponding slave-side mortar element
          Epetra_IntVector impltypes_col(*interface.SlaveColElements());
          Core::LinAlg::Export(impltypes_row, impltypes_col);
          imortarcells_[condid].resize(imortarcells.size());
          for (unsigned icell = 0; icell < imortarcells.size(); ++icell)
          {
            imortarcells_[condid][icell] =
                std::pair<Teuchos::RCP<Mortar::IntCell>, Inpar::ScaTra::ImplType>(
                    imortarcells[icell], static_cast<Inpar::ScaTra::ImplType>(
                                             impltypes_col[interface.SlaveColElements()->LID(
                                                 imortarcells[icell]->GetSlaveId())]));
          }
        }

        else
        {
          // match slave-side and master-side elements at mortar interface
          switch (interface.SearchAlg())
          {
            case Inpar::Mortar::search_bfele:
            {
              interface.evaluate_search_brute_force(interface.SearchParam());
              break;
            }

            case Inpar::Mortar::search_binarytree:
            {
              interface.evaluate_search_binarytree();
              break;
            }

            default:
            {
              FOUR_C_THROW("Invalid search algorithm!");
              break;
            }
          }

          // evaluate normal vectors associated with slave-side nodes
          interface.evaluate_nodal_normals();

          // extract slave-side noderowmap
          const Epetra_Map& noderowmap_slave = *interface.SlaveRowNodes();

          // initialize vector for node-to-segment connectivity, i.e., for pairings between slave
          // nodes and master elements
          Teuchos::RCP<Epetra_IntVector>& islavenodestomasterelements =
              islavenodestomasterelements_[condid];
          islavenodestomasterelements = Teuchos::rcp(new Epetra_IntVector(noderowmap_slave, false));
          islavenodestomasterelements->PutValue(-1);

          // initialize vector for physical implementation types of slave-side nodes
          Teuchos::RCP<Epetra_IntVector>& islavenodesimpltypes = islavenodesimpltypes_[condid];
          islavenodesimpltypes = Teuchos::rcp(new Epetra_IntVector(noderowmap_slave, false));
          islavenodesimpltypes->PutValue(Inpar::ScaTra::impltype_undefined);

          // loop over all slave-side nodes
          for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
          {
            // extract slave-side node
            auto* const slavenode =
                dynamic_cast<Mortar::Node*>(idiscret.gNode(noderowmap_slave.GID(inode)));
            if (!slavenode)
              FOUR_C_THROW("Couldn't extract slave-side mortar node from mortar discretization!");

            // find associated master-side elements
            std::vector<Mortar::Element*> master_mortar_elements(0, nullptr);
            interface.FindMEles(*slavenode, master_mortar_elements);

            // loop over all master-side elements
            for (auto* master_mortar_ele : master_mortar_elements)
            {
              // extract master-side element
              // project slave-side node onto master-side element
              std::array<double, 2> coordinates_master;
              double dummy(0.);
              Mortar::Projector::Impl(*master_mortar_ele)
                  ->project_gauss_point_auxn3_d(slavenode->X().data(), slavenode->MoData().n(),
                      *master_mortar_ele, coordinates_master.data(), dummy);

              // check whether projected node lies inside master-side element
              if (master_mortar_ele->Shape() == Core::FE::CellType::quad4)
              {
                if (coordinates_master[0] < -1. - ntsprojtol_ or
                    coordinates_master[1] < -1. - ntsprojtol_ or
                    coordinates_master[0] > 1. + ntsprojtol_ or
                    coordinates_master[1] > 1. + ntsprojtol_)
                  // projected node lies outside master-side element
                  continue;
              }

              else if (master_mortar_ele->Shape() == Core::FE::CellType::tri3)
              {
                if (coordinates_master[0] < -ntsprojtol_ or coordinates_master[1] < -ntsprojtol_ or
                    coordinates_master[0] + coordinates_master[1] > 1. + 2 * ntsprojtol_)
                  // projected node lies outside master-side element
                  continue;
              }

              else
                FOUR_C_THROW("Invalid discretization type of master-side element!");

              // projected node lies inside master-side element
              (*islavenodestomasterelements)[inode] = master_mortar_ele->Id();
              break;
            }

            // safety check
            if ((*islavenodestomasterelements)[inode] == -1)
              FOUR_C_THROW("Couldn't match slave-side node with master-side element!");

            // determine physical implementation type of slave-side node based on first associated
            // element
            (*islavenodesimpltypes)[inode] = dynamic_cast<Discret::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(
                    kinetics_slave_cond.second->Geometry()[slavenode->Elements()[0]->Id()])
                    ->parent_element())
                                                 ->ImplType();
          }

          // extract slave-side elerowmap
          const Epetra_Map& elecolmap_slave = *interface.SlaveColElements();

          // initialize vector for physical implementation types of slave-side elements
          Epetra_IntVector islaveelementsimpltypes(elecolmap_slave, false);
          islaveelementsimpltypes.PutValue(Inpar::ScaTra::impltype_undefined);

          // loop over all slave-side elements
          for (int ielement = 0; ielement < elecolmap_slave.NumMyElements(); ++ielement)
          {
            // determine physical implementation type of current slave-side element
            islaveelementsimpltypes[ielement] = dynamic_cast<Discret::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<Core::Elements::FaceElement>(
                    kinetics_slave_cond.second->Geometry()[elecolmap_slave.GID(ielement)])
                    ->parent_element())
                                                    ->ImplType();
          }

          // create parameter list for slave-side elements
          Teuchos::ParameterList eleparams;

          // set action for slave-side elements
          eleparams.set<int>("action", Inpar::S2I::evaluate_nodal_area_fractions);

          // compute vector for lumped interface area fractions associated with slave-side nodes
          const Epetra_Map& dofrowmap_slave = *interface.SlaveRowDofs();
          Teuchos::RCP<Epetra_Vector> islavenodeslumpedareas_dofvector =
              Core::LinAlg::CreateVector(dofrowmap_slave);
          evaluate_mortar_elements(elecolmap_slave, islaveelementsimpltypes, idiscret, eleparams,
              Teuchos::null, Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
              Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
              Inpar::S2I::side_undefined, Inpar::S2I::side_undefined, Teuchos::null,
              Inpar::S2I::side_undefined, Inpar::S2I::side_undefined,
              islavenodeslumpedareas_dofvector, Inpar::S2I::side_slave, Teuchos::null,
              Inpar::S2I::side_undefined);

          // transform map of result vector
          Teuchos::RCP<Epetra_Vector>& islavenodeslumpedareas = islavenodeslumpedareas_[condid];
          islavenodeslumpedareas = Core::LinAlg::CreateVector(noderowmap_slave);
          for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
          {
            (*islavenodeslumpedareas)[inode] =
                (*islavenodeslumpedareas_dofvector)[dofrowmap_slave.LID(
                    idiscret.Dof(idiscret.gNode(noderowmap_slave.GID(inode)), 0))];
          }
        }

        // build interface maps
        imastermap = Core::LinAlg::MergeMap(imastermap, interface.MasterRowDofs(), false);
        islavemap = Core::LinAlg::MergeMap(islavemap, interface.SlaveRowDofs(), false);
        ifullmap = Core::LinAlg::MergeMap(ifullmap,
            Core::LinAlg::MergeMap(interface.MasterRowDofs(), interface.SlaveRowDofs(), false),
            false);
        if (imortarredistribution_)
        {
          imastermap_ = Core::LinAlg::MergeMap(imastermap_, icoupmortar.MasterDofMap(), false);
          islavemap_ = Core::LinAlg::MergeMap(islavemap_, icoupmortar.SlaveDofMap(), false);
        }
      }

      // generate interior and interface maps
      std::vector<Teuchos::RCP<const Epetra_Map>> imaps;
      imaps.emplace_back(
          Core::LinAlg::SplitMap(*(scatratimint_->discretization()->dof_row_map()), *ifullmap));
      imaps.emplace_back(islavemap);
      imaps.emplace_back(imastermap);

      // initialize global map extractor
      interfacemaps_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
          *(scatratimint_->discretization()->dof_row_map()), imaps));
      interfacemaps_->check_for_valid_map_extractor();

      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_slave or couplingtype_ == Inpar::S2I::coupling_nts_standard)
      {
        // initialize auxiliary system matrix for slave side
        islavematrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interfacemaps_->Map(1), 81));

        // initialize auxiliary residual vector for slave side
        islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*interfacemaps_->Map(1)));
      }

      if (couplingtype_ == Inpar::S2I::coupling_mortar_standard or
          lmside_ == Inpar::S2I::side_master or couplingtype_ == Inpar::S2I::coupling_nts_standard)
      {
        // initialize auxiliary system matrix for master side
        imastermatrix_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
            *interfacemaps_->Map(2), 81, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));

        // initialize auxiliary residual vector for master side
        imasterresidual_ = Teuchos::rcp(new Epetra_FEVector(*interfacemaps_->Map(2)));
      }

      switch (couplingtype_)
      {
        case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
        case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
        case Inpar::S2I::coupling_mortar_condensed_petrov:
        case Inpar::S2I::coupling_mortar_condensed_bubnov:
        {
          if (lmside_ == Inpar::S2I::side_slave)
          {
            D_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interfacemaps_->Map(1), 81));
            M_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interfacemaps_->Map(1), 81));
            if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov or
                couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
              E_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*interfacemaps_->Map(1), 81));
          }
          else
          {
            D_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
                *interfacemaps_->Map(2), 81, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
            M_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
                *interfacemaps_->Map(2), 81, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
            if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov or
                couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
              E_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(
                  *interfacemaps_->Map(2), 81, true, false, Core::LinAlg::SparseMatrix::FE_MATRIX));
          }

          // loop over all scatra-scatra coupling interfaces
          for (auto& kinetics_slave_cond : kinetics_conditions_meshtying_slaveside_)
          {
            // create parameter list for mortar integration cells
            Teuchos::ParameterList params;

            // add current condition to parameter list
            params.set<Core::Conditions::Condition*>("condition", kinetics_slave_cond.second);

            // set action
            params.set<int>("action", Inpar::S2I::evaluate_mortar_matrices);

            // evaluate mortar integration cells at current interface
            evaluate_mortar_cells(icoupmortar_[kinetics_slave_cond.first]->Interface()->Discret(),
                params, D_,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_slave
                                                  : Inpar::S2I::side_master,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_slave
                                                  : Inpar::S2I::side_master,
                M_,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_slave
                                                  : Inpar::S2I::side_master,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_master
                                                  : Inpar::S2I::side_slave,
                E_,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_slave
                                                  : Inpar::S2I::side_master,
                lmside_ == Inpar::S2I::side_slave ? Inpar::S2I::side_slave
                                                  : Inpar::S2I::side_master,
                Teuchos::null, Inpar::S2I::side_undefined, Inpar::S2I::side_undefined,
                Teuchos::null, Inpar::S2I::side_undefined, Teuchos::null,
                Inpar::S2I::side_undefined);
          }

          // finalize mortar matrices D, M, and E
          D_->Complete();
          if (lmside_ == Inpar::S2I::side_slave)
            M_->Complete(*interfacemaps_->Map(2), *interfacemaps_->Map(1));
          else
            M_->Complete(*interfacemaps_->Map(1), *interfacemaps_->Map(2));
          if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov or
              couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
            E_->Complete();

          switch (couplingtype_)
          {
            case Inpar::S2I::coupling_mortar_condensed_petrov:
            case Inpar::S2I::coupling_mortar_condensed_bubnov:
            {
              // set up mortar projector P
              Teuchos::RCP<Epetra_Vector> D_diag(Teuchos::null);
              if (lmside_ == Inpar::S2I::side_slave)
                D_diag = Core::LinAlg::CreateVector(*interfacemaps_->Map(1));
              else
                D_diag = Core::LinAlg::CreateVector(*interfacemaps_->Map(2));
              if (D_->ExtractDiagonalCopy(*D_diag))
                FOUR_C_THROW("Couldn't extract main diagonal from mortar matrix D!");
              if (D_diag->Reciprocal(*D_diag))
                FOUR_C_THROW("Couldn't invert main diagonal entries of mortar matrix D!");

              P_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*M_));
              if (P_->LeftScale(*D_diag)) FOUR_C_THROW("Setup of mortar projector P failed!");

              // free memory
              if (!slaveonly_)
              {
                D_ = Teuchos::null;
                M_ = Teuchos::null;
              }

              if (couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
              {
                // set up mortar projector Q
                Q_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*E_));
                if (Q_->LeftScale(*D_diag)) FOUR_C_THROW("Setup of mortar projector Q failed!");

                // free memory
                E_ = Teuchos::null;
              }

              break;
            }

            case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
            case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
            {
              // determine number of Lagrange multiplier dofs owned by each processor
              const Epetra_Comm& comm(scatratimint_->discretization()->Comm());
              const int numproc(comm.NumProc());
              const int mypid(comm.MyPID());
              std::vector<int> localnumlmdof(numproc, 0);
              std::vector<int> globalnumlmdof(numproc, 0);
              if (lmside_ == Inpar::S2I::side_slave)
                localnumlmdof[mypid] = interfacemaps_->Map(1)->NumMyElements();
              else
                localnumlmdof[mypid] = interfacemaps_->Map(2)->NumMyElements();
              comm.SumAll(localnumlmdof.data(), globalnumlmdof.data(), numproc);

              // for each processor, determine offset of minimum Lagrange multiplier dof GID w.r.t.
              // maximum standard dof GID
              int offset(0);
              for (int ipreviousproc = 0; ipreviousproc < mypid; ++ipreviousproc)
                offset += globalnumlmdof[ipreviousproc];

              // for each processor, determine Lagrange multiplier dof GIDs
              std::vector<int> lmdofgids(globalnumlmdof[mypid], 0);
              for (int lmdoflid = 0; lmdoflid < globalnumlmdof[mypid]; ++lmdoflid)
                lmdofgids[lmdoflid] =
                    scatratimint_->dof_row_map()->MaxAllGID() + 1 + offset + lmdoflid;

              // build Lagrange multiplier dofrowmap
              const Teuchos::RCP<Epetra_Map> lmdofrowmap = Teuchos::rcp(
                  new Epetra_Map(-1, (int)lmdofgids.size(), lmdofgids.data(), 0, comm));

              // initialize vectors associated with Lagrange multiplier dofs
              lm_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));
              lmresidual_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));
              lmincrement_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));

              // initialize extended map extractor
              Teuchos::RCP<Epetra_Map> extendedmap = Core::LinAlg::MergeMap(
                  *(scatratimint_->discretization()->dof_row_map()), *lmdofrowmap, false);
              extendedmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(
                  *extendedmap, lmdofrowmap, scatratimint_->discretization()->dof_row_map()));
              extendedmaps_->check_for_valid_map_extractor();

              // transform range map of mortar matrices D and M from slave-side dofrowmap to
              // Lagrange multiplier dofrowmap
              D_ = Mortar::MatrixRowTransformGIDs(D_, lmdofrowmap);
              M_ = Mortar::MatrixRowTransformGIDs(M_, lmdofrowmap);

              if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_petrov)
              {
                // transform domain map of mortar matrix D from slave-side dofrowmap to Lagrange
                // multiplier dofrowmap and store transformed matrix as mortar matrix E
                E_ = Mortar::MatrixColTransformGIDs(D_, lmdofrowmap);
              }
              else
              {
                // transform domain and range maps of mortar matrix E from slave-side dofrowmap to
                // Lagrange multiplier dofrowmap
                E_ = Mortar::MatrixRowColTransformGIDs(E_, lmdofrowmap, lmdofrowmap);
              }

              break;
            }

            default:
            {
              FOUR_C_THROW("Invalid type of mortar meshtying!");
              break;
            }
          }

          break;
        }

        default:
        {
          // do nothing
          break;
        }
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Type of mortar meshtying for scatra-scatra interface coupling not recognized!");
      break;
    }
  }

  // further initializations depending on type of global system matrix
  switch (matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
    {
      // nothing needs to be done in this case
      break;
    }
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      // safety check
      if (!scatratimint_->Solver()->Params().isSublist("AMGnxn Parameters"))
        FOUR_C_THROW(
            "Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize map extractors associated with blocks of global system matrix
      build_block_map_extractors();

      break;
    }
    default:
    {
      FOUR_C_THROW(
          "%i is not a valid 'ScaTra::MatrixType'. Set a valid 'ScaTra::MatrixType' in your input "
          "file!",
          static_cast<int>(matrixtype_));
      break;
    }
  }

  // extract boundary condition for scatra-scatra interface layer growth
  const Core::Conditions::Condition* const condition =
      scatratimint_->discretization()->GetCondition("S2IKineticsGrowth");

  // setup evaluation of scatra-scatra interface layer growth if applicable
  if (condition)
  {
    // perform setup depending on evaluation method
    switch (intlayergrowth_evaluation_)
    {
      case Inpar::S2I::growth_evaluation_monolithic:
      case Inpar::S2I::growth_evaluation_semi_implicit:
      {
        // extract map associated with scatra-scatra interface layer thicknesses
        const Teuchos::RCP<const Epetra_Map>& dofrowmap_growth = scatratimint_->dof_row_map(2);

        // initialize state vector of discrete scatra-scatra interface layer thicknesses at time n
        growthn_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));

        // additional initializations for monolithic solution approach
        if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
        {
          // initialize extended map extractor
          const Epetra_Map* const dofrowmap_scatra = scatratimint_->discretization()->dof_row_map();
          extendedmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(
              *Core::LinAlg::MergeMap(*dofrowmap_scatra, *dofrowmap_growth, false),
              scatratimint_->dof_row_map(2), dofrowmap_scatra));
          extendedmaps_->check_for_valid_map_extractor();

          // initialize additional state vectors
          growthnp_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthdtn_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthdtnp_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthhist_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthresidual_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthincrement_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));

          // initialize map extractors and global matrix blocks
          growthgrowthblock_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_growth, 81));
          switch (matrixtype_)
          {
            case Core::LinAlg::MatrixType::sparse:
            {
              // initialize extended map extractor associated with blocks of global system matrix
              extendedblockmaps_ = extendedmaps_;

              scatragrowthblock_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_scatra,
                  81));  // We actually don't really need the entire scalar transport dofrowmap
                         // here, but only a submap associated with all (slave-side and master-side)
                         // interfacial degrees of freedom. However, this will later cause an error
                         // in debug mode when assigning the scatra-growth matrix block to the
                         // global system matrix in the Solve() routine.
              growthscatrablock_ =
                  Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_growth, 81));

              break;
            }

            case Core::LinAlg::MatrixType::block_condition:
            case Core::LinAlg::MatrixType::block_condition_dof:
            {
              // initialize map extractor associated with all degrees of freedom for scatra-scatra
              // interface layer growth
              blockmapgrowth_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*dofrowmap_growth,
                  std::vector<Teuchos::RCP<const Epetra_Map>>(1, dofrowmap_growth)));
              blockmapgrowth_->check_for_valid_map_extractor();

              // initialize extended map extractor associated with blocks of global system matrix
              const unsigned nblockmaps = scatratimint_->BlockMaps()->NumMaps();
              std::vector<Teuchos::RCP<const Epetra_Map>> extendedblockmaps(
                  nblockmaps + 1, Teuchos::null);
              for (int iblockmap = 0; iblockmap < static_cast<int>(nblockmaps); ++iblockmap)
                extendedblockmaps[iblockmap] = scatratimint_->BlockMaps()->Map(iblockmap);
              extendedblockmaps[nblockmaps] = dofrowmap_growth;
              extendedblockmaps_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(
                  *extendedmaps_->FullMap(), extendedblockmaps));
              extendedblockmaps_->check_for_valid_map_extractor();

              scatragrowthblock_ = Teuchos::rcp(
                  new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                      *blockmapgrowth_, *scatratimint_->BlockMaps(), 81, false, true));
              growthscatrablock_ = Teuchos::rcp(
                  new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                      *scatratimint_->BlockMaps(), *blockmapgrowth_, 81, false, true));

              break;
            }

            default:
            {
              FOUR_C_THROW(
                  "Type of global system matrix for scatra-scatra interface coupling involving "
                  "interface layer growth not recognized!");
              break;
            }
          }

          // initialize extended system matrix including rows and columns associated with
          // scatra-scatra interface layer thickness variables
          extendedsystemmatrix_ = Teuchos::rcp(
              new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
                  *extendedblockmaps_, *extendedblockmaps_));
        }

        // loop over all boundary conditions for scatra-scatra interface coupling
        for (auto& icond : s2ikinetics_conditions)
        {
          // check whether current boundary condition is associated with boundary condition for
          // scatra-scatra interface layer growth
          if (icond->parameters().Get<int>("ConditionID") ==
              condition->parameters().Get<int>("ConditionID"))
            // copy conductivity parameter
            icond->parameters().Add(
                "conductivity", condition->parameters().Get<double>("conductivity"));
        }

        // extract initial scatra-scatra interface layer thickness from condition
        const double initthickness = condition->parameters().Get<double>("initial thickness");

        // extract nodal cloud from condition
        const std::vector<int>* nodegids = condition->GetNodes();

        // loop over all nodes
        for (int nodegid : *nodegids)
        {
          // extract global ID of current node
          // process only nodes stored by current processor
          if (scatratimint_->discretization()->HaveGlobalNode(nodegid))
          {
            // extract current node
            const Core::Nodes::Node* const node = scatratimint_->discretization()->gNode(nodegid);

            // process only nodes owned by current processor
            if (node->Owner() == scatratimint_->discretization()->Comm().MyPID())
            {
              // extract local ID of scatra-scatra interface layer thickness variable associated
              // with current node
              const int doflid_growth = scatratimint_->discretization()->dof_row_map(2)->LID(
                  scatratimint_->discretization()->Dof(2, node, 0));
              if (doflid_growth < 0)
              {
                FOUR_C_THROW(
                    "Couldn't extract local ID of scatra-scatra interface layer thickness "
                    "variable!");
              }

              // set initial value
              (*growthn_)[doflid_growth] = initthickness;
            }  // nodes owned by current processor
          }    // nodes stored by current processor
        }      // loop over all nodes

        // copy initial state
        if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
          growthnp_->Update(1., *growthn_, 0.);

        break;
      }

      default:
      {
        FOUR_C_THROW(
            "Unknown evaluation method for scatra-scatra interface coupling involving interface "
            "layer growth!");
        break;
      }
    }
  }

  // instantiate appropriate equilibration class
  auto equilibration_method =
      std::vector<Core::LinAlg::EquilibrationMethod>(1, scatratimint_->EquilibrationMethod());
  equilibration_ = Core::LinAlg::BuildEquilibration(matrixtype_, equilibration_method,
      (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic
              ? extendedmaps_->FullMap()
              : Teuchos::rcp(
                    new const Epetra_Map(*scatratimint_->discretization()->dof_row_map()))));
}  // ScaTra::meshtying_strategy_s2_i::setup_meshtying

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::compute_time_derivative() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
  {
    // compute inverse time factor 1./(theta*dt)
    const double timefac_inverse =
        1. / scatratimint_->scatra_time_parameter_list()->get<double>("time factor");

    // compute state vector of time derivatives of discrete scatra-scatra interface layer
    // thicknesses
    growthdtnp_->Update(timefac_inverse, *growthnp_, -timefac_inverse, *growthhist_, 0.);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::Update() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
  {
    // update state vectors
    growthn_->Update(1., *growthnp_, 0.);
    growthdtn_->Update(1., *growthdtnp_, 0.);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::set_element_general_parameters(
    Teuchos::ParameterList& parameters) const
{
  // add local Newton-Raphson convergence tolerance for scatra-scatra interface layer growth to
  // parameter list
  parameters.set<double>("intlayergrowth_convtol", intlayergrowth_convtol_);

  // add maximum number of local Newton-Raphson iterations for scatra-scatra interface layer growth
  // to parameter list
  parameters.set<unsigned>("intlayergrowth_itemax", intlayergrowth_itemax_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::set_condition_specific_sca_tra_parameters(
    Core::Conditions::Condition& s2icondition) const
{
  Teuchos::ParameterList conditionparams;

  // fill the parameter list
  write_s2_i_kinetics_specific_sca_tra_parameters_to_parameter_list(s2icondition, conditionparams);

  // call standard loop over elements
  scatratimint_->discretization()->Evaluate(
      conditionparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::
    write_s2_i_kinetics_specific_sca_tra_parameters_to_parameter_list(
        Core::Conditions::Condition& s2ikinetics_cond,
        Teuchos::ParameterList& s2icouplingparameters)
{
  // get kinetic model and condition type
  const int kineticmodel = s2ikinetics_cond.parameters().Get<int>("kinetic model");
  const Core::Conditions::ConditionType conditiontype = s2ikinetics_cond.Type();

  // set action, kinetic model, condition type and numscal
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::set_scatra_ele_boundary_parameter, s2icouplingparameters);
  s2icouplingparameters.set<int>("kinetic model", kineticmodel);
  s2icouplingparameters.set<Core::Conditions::ConditionType>("condition type", conditiontype);

  // set the condition type specific parameters
  switch (conditiontype)
  {
    case Core::Conditions::ConditionType::S2IKinetics:
    {
      // set the kinetic model specific parameters
      switch (kineticmodel)
      {
        case Inpar::S2I::kinetics_constperm:
        case Inpar::S2I::kinetics_linearperm:
        {
          s2icouplingparameters.set<int>(
              "numscal", s2ikinetics_cond.parameters().Get<int>("numscal"));
          s2icouplingparameters.set<const std::vector<double>*>("permeabilities",
              &s2ikinetics_cond.parameters().Get<std::vector<double>>("permeabilities"));
          s2icouplingparameters.set<int>(
              "is_pseudo_contact", s2ikinetics_cond.parameters().Get<int>("is_pseudo_contact"));
          break;
        }

        case Inpar::S2I::kinetics_constantinterfaceresistance:
        {
          s2icouplingparameters.set<double>(
              "resistance", s2ikinetics_cond.parameters().Get<double>("resistance"));
          s2icouplingparameters.set<const std::vector<int>*>(
              "onoff", &s2ikinetics_cond.parameters().Get<std::vector<int>>("onoff"));
          s2icouplingparameters.set<int>(
              "numelectrons", s2ikinetics_cond.parameters().Get<int>("e-"));
          s2icouplingparameters.set<int>(
              "is_pseudo_contact", s2ikinetics_cond.parameters().Get<int>("is_pseudo_contact"));
          break;
        }

        case Inpar::S2I::kinetics_nointerfaceflux:
        {
          // do nothing
          break;
        }

        case Inpar::S2I::kinetics_butlervolmer:
        case Inpar::S2I::kinetics_butlervolmerlinearized:
        case Inpar::S2I::kinetics_butlervolmerreduced:
        case Inpar::S2I::kinetics_butlervolmerreducedcapacitance:
        case Inpar::S2I::kinetics_butlervolmerreducedlinearized:
        case Inpar::S2I::kinetics_butlervolmerpeltier:
        case Inpar::S2I::kinetics_butlervolmerresistance:
        case Inpar::S2I::kinetics_butlervolmerreducedthermoresistance:
        case Inpar::S2I::kinetics_butlervolmerreducedresistance:
        {
          s2icouplingparameters.set<int>(
              "numscal", s2ikinetics_cond.parameters().Get<int>("numscal"));
          s2icouplingparameters.set<const std::vector<int>*>("stoichiometries",
              &s2ikinetics_cond.parameters().Get<std::vector<int>>("stoichiometries"));
          s2icouplingparameters.set<int>(
              "numelectrons", s2ikinetics_cond.parameters().Get<int>("e-"));
          s2icouplingparameters.set<double>(
              "k_r", s2ikinetics_cond.parameters().Get<double>("k_r"));
          s2icouplingparameters.set<double>(
              "alpha_a", s2ikinetics_cond.parameters().Get<double>("alpha_a"));
          s2icouplingparameters.set<double>(
              "alpha_c", s2ikinetics_cond.parameters().Get<double>("alpha_c"));
          s2icouplingparameters.set<int>(
              "is_pseudo_contact", s2ikinetics_cond.parameters().Get<int>("is_pseudo_contact"));

          if (kineticmodel == Inpar::S2I::kinetics_butlervolmerreducedcapacitance)
            s2icouplingparameters.set<double>(
                "capacitance", s2ikinetics_cond.parameters().Get<double>("capacitance"));

          if (kineticmodel == Inpar::S2I::kinetics_butlervolmerpeltier)
            s2icouplingparameters.set<double>(
                "peltier", s2ikinetics_cond.parameters().Get<double>("peltier"));

          if (kineticmodel == Inpar::S2I::kinetics_butlervolmerresistance or
              kineticmodel == Inpar::S2I::kinetics_butlervolmerreducedresistance)
          {
            s2icouplingparameters.set<double>(
                "resistance", s2ikinetics_cond.parameters().Get<double>("resistance"));
            s2icouplingparameters.set<double>("CONVTOL_IMPLBUTLERVOLMER",
                s2ikinetics_cond.parameters().Get<double>("CONVTOL_IMPLBUTLERVOLMER"));
            s2icouplingparameters.set<int>("ITEMAX_IMPLBUTLERVOLMER",
                s2ikinetics_cond.parameters().Get<int>("ITEMAX_IMPLBUTLERVOLMER"));
          }

          if (kineticmodel == Inpar::S2I::kinetics_butlervolmerreducedthermoresistance)
          {
            s2icouplingparameters.set<double>(
                "thermoperm", s2ikinetics_cond.parameters().Get<double>("thermoperm"));
            s2icouplingparameters.set<double>("molar_heat_capacity",
                s2ikinetics_cond.parameters().Get<double>("molar_heat_capacity"));
          }
          break;
        }

        default:
        {
          FOUR_C_THROW("Not implemented for this kinetic model: %i", kineticmodel);
          break;
        }
      }
      break;
    }

    case Core::Conditions::ConditionType::S2IKineticsGrowth:
    {
      // set the kinetic model specific parameters
      switch (kineticmodel)
      {
        case Inpar::S2I::growth_kinetics_butlervolmer:
        {
          s2icouplingparameters.set<int>(
              "numscal", s2ikinetics_cond.parameters().Get<int>("numscal"));
          s2icouplingparameters.set<const std::vector<int>*>("stoichiometries",
              &s2ikinetics_cond.parameters().Get<std::vector<int>>("stoichiometries"));
          s2icouplingparameters.set<int>(
              "numelectrons", s2ikinetics_cond.parameters().Get<int>("e-"));
          s2icouplingparameters.set<double>(
              "k_r", s2ikinetics_cond.parameters().Get<double>("k_r"));
          s2icouplingparameters.set<double>(
              "alpha_a", s2ikinetics_cond.parameters().Get<double>("alpha_a"));
          s2icouplingparameters.set<double>(
              "alpha_c", s2ikinetics_cond.parameters().Get<double>("alpha_c"));
          s2icouplingparameters.set<double>(
              "density", s2ikinetics_cond.parameters().Get<double>("density"));
          s2icouplingparameters.set<double>(
              "molar mass", s2ikinetics_cond.parameters().Get<double>("molar mass"));
          s2icouplingparameters.set<double>(
              "regpar", s2ikinetics_cond.parameters().Get<double>("regularization parameter"));
          s2icouplingparameters.set<int>(
              "regtype", s2ikinetics_cond.parameters().Get<int>("regularization type"));
          s2icouplingparameters.set<double>(
              "conductivity", s2ikinetics_cond.parameters().Get<double>("conductivity"));
          break;
        }

        default:
        {
          FOUR_C_THROW("Not implemented for this kinetic model: %i", kineticmodel);
          break;
        }
      }
      break;
    }

    default:
    {
      FOUR_C_THROW("Not implemented for this condition type: %i", conditiontype);
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::SetOldPartOfRHS() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
  {
    // compute factor dt*(1-theta)
    const double factor = scatratimint_->Dt() -
                          scatratimint_->scatra_time_parameter_list()->get<double>("time factor");

    // compute history vector
    growthhist_->Update(1., *growthn_, factor, *growthdtn_, 0.);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::write_restart() const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_semi_implicit)
  {
    // output state vector of discrete scatra-scatra interface layer thicknesses
    scatratimint_->DiscWriter()->WriteVector("growthn", growthn_);

    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
    {
      // output state vector of time derivatives of discrete scatra-scatra interface layer
      // thicknesses
      scatratimint_->DiscWriter()->WriteVector("growthdtn", growthdtn_);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input) const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_semi_implicit)
  {
    // initialize reader
    Teuchos::RCP<Core::IO::DiscretizationReader> reader(Teuchos::null);
    if (input == Teuchos::null)
      reader = Teuchos::rcp(new Core::IO::DiscretizationReader(
          scatratimint_->discretization(), Global::Problem::Instance()->InputControlFile(), step));
    else
      reader = Teuchos::rcp(
          new Core::IO::DiscretizationReader(scatratimint_->discretization(), input, step));

    // read state vector of discrete scatra-scatra interface layer thicknesses
    reader->ReadVector(growthn_, "growthn");

    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
    {
      // read state vector of time derivatives of discrete scatra-scatra interface layer thicknesses
      reader->ReadVector(growthdtn_, "growthdtn");

      // copy restart state
      growthnp_->Update(1., *growthn_, 0.);
      growthdtnp_->Update(1., *growthdtn_, 0.);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::FE::Discretization& ScaTra::MeshtyingStrategyS2I::mortar_discretization(
    const int& condid) const
{
  return icoupmortar_.at(condid)->Interface()->Discret();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::Output() const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_semi_implicit)
  {
    // extract relevant state vector of discrete scatra-scatra interface layer thicknesses based on
    // map of scatra-scatra interface layer thickness variables
    const Epetra_Vector& growth =
        intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic ? *growthnp_
                                                                               : *growthn_;

    // for proper output, initialize target state vector of discrete scatra-scatra interface layer
    // thicknesses based on map of row nodes
    const Teuchos::RCP<Epetra_Vector> intlayerthickness =
        Teuchos::rcp(new Epetra_Vector(*scatratimint_->discretization()->NodeRowMap(), true));

    // extract boundary condition for scatra-scatra interface layer growth
    const Core::Conditions::Condition* const condition =
        scatratimint_->discretization()->GetCondition("S2IKineticsGrowth");

    // extract nodal cloud from condition
    const std::vector<int>* nodegids = condition->GetNodes();

    // loop over all nodes
    for (int nodegid : *nodegids)
    {
      // extract global ID of current node
      // process only nodes stored by current processor
      if (scatratimint_->discretization()->HaveGlobalNode(nodegid))
      {
        // extract current node
        const Core::Nodes::Node* const node = scatratimint_->discretization()->gNode(nodegid);

        // process only nodes owned by current processor
        if (node->Owner() == scatratimint_->discretization()->Comm().MyPID())
        {
          // extract local ID of current node
          const int nodelid = scatratimint_->discretization()->NodeRowMap()->LID(nodegid);
          if (nodelid < 0) FOUR_C_THROW("Couldn't extract local node ID!");

          // extract local ID of scatra-scatra interface layer thickness variable associated with
          // current node
          const int doflid_growth = scatratimint_->discretization()->dof_row_map(2)->LID(
              scatratimint_->discretization()->Dof(2, node, 0));
          if (doflid_growth < 0)
            FOUR_C_THROW(
                "Couldn't extract local ID of scatra-scatra interface layer thickness variable!");

          // copy thickness variable into target state vector of discrete scatra-scatra interface
          // layer thicknesses
          (*intlayerthickness)[nodelid] = growth[doflid_growth];
        }  // nodes owned by current processor
      }    // nodes stored by current processor
    }      // loop over all nodes

    // output target state vector of discrete scatra-scatra interface layer thicknesses
    scatratimint_->DiscWriter()->WriteVector("intlayerthickness", intlayerthickness);
  }
  if (output_interface_flux_) OutputInterfaceFlux();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::OutputInterfaceFlux() const
{
  add_time_integration_specific_vectors();
  std::vector<Core::Conditions::Condition*> s2ikinetics_conditions(0, nullptr);
  scatratimint_->discretization()->GetCondition("S2IKinetics", s2ikinetics_conditions);

  std::map<std::string, std::vector<double>> output_data;

  for (auto* s2ikinetics_cond : s2ikinetics_conditions)
  {
    // only slave side has relevant information
    if (s2ikinetics_cond->parameters().Get<int>("interface side") ==
        static_cast<int>(Inpar::S2I::side_slave))
    {
      const int condition_id = s2ikinetics_cond->parameters().Get<int>("ConditionID");
      auto s2i_flux =
          Teuchos::rcp(new Core::LinAlg::SerialDenseVector(scatratimint_->NumDofPerNode()));
      auto boundaryint_vector = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(1));
      {
        Teuchos::ParameterList condparams;

        Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
            "action", ScaTra::BoundaryAction::calc_s2icoupling_flux, condparams);

        scatratimint_->discretization()->EvaluateScalars(
            condparams, s2i_flux, "S2IKinetics", condition_id);
      }

      {
        Teuchos::ParameterList condparams;

        // overwrite action in parameter list
        Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
            "action", ScaTra::BoundaryAction::calc_boundary_integral, condparams);

        // compute value of boundary integral
        scatratimint_->discretization()->EvaluateScalars(
            condparams, boundaryint_vector, "S2IKinetics", condition_id);
      }

      // extract value of boundary integral
      const double boundaryint = (*boundaryint_vector)(0);
      s2i_flux->scale(1.0 / boundaryint);

      std::stringstream condition_flux_name, condition_area_name;
      condition_flux_name << "flux_S2I_condition_" << condition_id;
      condition_area_name << "area_S2I_condition_" << condition_id;

      FOUR_C_ASSERT(
          runtime_csvwriter_.has_value(), "internal error: runtime csv writer not created.");
      output_data[condition_flux_name.str()] = {(*s2i_flux)[1]};
      output_data[condition_area_name.str()] = {boundaryint};
    }
  }
  runtime_csvwriter_->WriteDataToFile(scatratimint_->Time(), scatratimint_->Step(), output_data);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::explicit_predictor() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
    // predict state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    growthnp_->Update(scatratimint_->Dt(), *growthdtn_, 1.);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::ExtractMatrixRows(
    const Core::LinAlg::SparseMatrix& matrix,  //!< source matrix
    Core::LinAlg::SparseMatrix& rows,          //!< destination matrix
    const Epetra_Map& rowmap                   //!< map of matrix rows to be extracted
)
{
  // safety check
  if (rows.Filled())
    FOUR_C_THROW("Source matrix rows cannot be extracted into filled destination matrix!");

  // loop over all source matrix rows to be extracted
  for (int doflid = 0; doflid < rowmap.NumMyElements(); ++doflid)
  {
    // determine global ID of current matrix row
    const int dofgid = rowmap.GID(doflid);
    if (dofgid < 0) FOUR_C_THROW("Couldn't find local ID %d in map!", doflid);

    // extract current matrix row from source matrix
    const int length = matrix.EpetraMatrix()->NumGlobalEntries(dofgid);
    int numentries(0);
    std::vector<double> values(length, 0.);
    std::vector<int> indices(length, 0);
    if (matrix.EpetraMatrix()->ExtractGlobalRowCopy(
            dofgid, length, numentries, values.data(), indices.data()))
      FOUR_C_THROW("Cannot extract matrix row with global ID %d from source matrix!", dofgid);

    // copy current source matrix row into destination matrix
    if (rows.EpetraMatrix()->InsertGlobalValues(dofgid, numentries, values.data(), indices.data()) <
        0)
      FOUR_C_THROW("Cannot insert matrix row with global ID %d into destination matrix!", dofgid);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::add_time_integration_specific_vectors() const
{
  // only relevant for scatra-scatra interface coupling with pairwise coinciding interface nodes
  if (couplingtype_ == Inpar::S2I::coupling_matching_nodes)
  {
    // add state vector containing master-side scatra degrees of freedom to scatra discretization
    interfacemaps_->InsertVector(
        icoup_->MasterToSlave(interfacemaps_->ExtractVector(*(scatratimint_->Phiafnp()), 2)), 1,
        imasterphi_on_slave_side_np_);
    scatratimint_->discretization()->set_state("imasterphinp", imasterphi_on_slave_side_np_);

    if (has_capacitive_contributions_)
    {
      interfacemaps_->InsertVector(
          interfacemaps_->ExtractVector(*(scatratimint_->Phidtnp()), 1), 1, islavephidtnp_);
      scatratimint_->discretization()->set_state("islavephidtnp", islavephidtnp_);
      interfacemaps_->InsertVector(
          icoup_->MasterToSlave(interfacemaps_->ExtractVector(*(scatratimint_->Phidtnp()), 2)), 1,
          imasterphidt_on_slave_side_np_);
      scatratimint_->discretization()->set_state("imasterphidtnp", imasterphidt_on_slave_side_np_);
    }
  }

  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_semi_implicit)
  {
    // extract relevant state vector of discrete scatra-scatra interface layer thicknesses
    const Teuchos::RCP<Epetra_Vector>& growth =
        intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic ? growthnp_
                                                                               : growthn_;

    // set state vector
    scatratimint_->discretization()->set_state(2, "growth", growth);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::compute_time_step_size(double& dt)
{
  // not implemented for standard scalar transport
  if (intlayergrowth_timestep_ > 0.)
  {
    FOUR_C_THROW(
        "Adaptive time stepping for scatra-scatra interface layer growth not implemented for "
        "standard scalar transport!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  init_conv_check_strategy();

  // extract boundary conditions for scatra-scatra interface layer growth
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> conditions;
  scatratimint_->discretization()->GetCondition("S2IKineticsGrowth", conditions);

  // initialize scatra-scatra interface layer growth
  if (conditions.size())
  {
    // safety checks
    if (conditions.size() != 1)
    {
      FOUR_C_THROW(
          "Can't have more than one boundary condition for scatra-scatra interface layer growth at "
          "the moment!");
    }
    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_none)
    {
      FOUR_C_THROW(
          "Invalid flag for evaluation of scatra-scatra interface coupling involving interface "
          "layer growth!");
    }
    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic and
        scatratimint_->MethodName() != Inpar::ScaTra::timeint_one_step_theta)
    {
      FOUR_C_THROW(
          "Monolithic evaluation of scatra-scatra interface layer growth only implemented for "
          "one-step-theta time integration scheme at the moment!");
    }
    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_semi_implicit and
        conditions[0]->parameters().Get<int>("regularization type") !=
            Inpar::S2I::RegularizationType::regularization_none)
    {
      FOUR_C_THROW(
          "No regularization implemented for semi-implicit evaluation of scatra-scatra interface "
          "layer growth!");
    }
    if (couplingtype_ != Inpar::S2I::coupling_matching_nodes)
    {
      FOUR_C_THROW(
          "Evaluation of scatra-scatra interface layer growth only implemented for conforming "
          "interface discretizations!");
    }

    // provide scalar transport discretization with additional dofset for scatra-scatra interface
    // layer thickness
    const Teuchos::RCP<Epetra_IntVector> numdofpernode =
        Teuchos::rcp(new Epetra_IntVector(*scatratimint_->discretization()->NodeColMap()));
    for (int inode = 0; inode < scatratimint_->discretization()->NumMyColNodes(); ++inode)
    {
      // add one degree of freedom for scatra-scatra interface layer growth to current node if
      // applicable
      if (scatratimint_->discretization()->lColNode(inode)->GetCondition("S2IKineticsGrowth"))
        (*numdofpernode)[inode] = 1;
    }

    int number_dofsets = scatratimint_->GetMaxDofSetNumber();
    Teuchos::RCP<Core::DOFSets::DofSetInterface> dofset =
        Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
            numdofpernode, Teuchos::null, Teuchos::null, true));
    if (scatratimint_->discretization()->AddDofSet(dofset) != ++number_dofsets)
      FOUR_C_THROW("Scalar transport discretization exhibits invalid number of dofsets!");
    scatratimint_->set_number_of_dof_set_growth(number_dofsets);

    // initialize linear solver for monolithic scatra-scatra interface coupling involving interface
    // layer growth
    if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
    {
      const int extendedsolver = Global::Problem::Instance()
                                     ->scalar_transport_dynamic_params()
                                     .sublist("S2I COUPLING")
                                     .get<int>("INTLAYERGROWTH_LINEAR_SOLVER");
      if (extendedsolver < 1)
      {
        FOUR_C_THROW(
            "Invalid ID of linear solver for monolithic scatra-scatra interface coupling involving "
            "interface layer growth!");
      }
      extendedsolver_ = Teuchos::rcp(
          new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(extendedsolver),
              scatratimint_->discretization()->Comm()));
    }
  }  // initialize scatra-scatra interface layer growth

  // safety check
  else if (intlayergrowth_evaluation_ != Inpar::S2I::growth_evaluation_none)
  {
    FOUR_C_THROW(
        "Cannot evaluate scatra-scatra interface coupling involving interface layer growth without "
        "specifying a corresponding boundary condition!");
  }

  // safety checks associated with adaptive time stepping for scatra-scatra interface layer growth
  if (intlayergrowth_timestep_ > 0.)
  {
    if (not Core::UTILS::IntegralValue<bool>(
            *scatratimint_->ScatraParameterList(), "ADAPTIVE_TIMESTEPPING"))
    {
      FOUR_C_THROW(
          "Adaptive time stepping for scatra-scatra interface layer growth requires "
          "ADAPTIVE_TIMESTEPPING flag to be set!");
    }
    if (not scatratimint_->discretization()->GetCondition("S2IKineticsGrowth"))
    {
      FOUR_C_THROW(
          "Adaptive time stepping for scatra-scatra interface layer growth requires corresponding "
          "boundary condition!");
    }
    if (intlayergrowth_timestep_ >= scatratimint_->Dt())
    {
      FOUR_C_THROW(
          "Adaptive time stepping for scatra-scatra interface layer growth requires that the "
          "modified time step size is smaller than the original time step size!");
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::build_block_map_extractors()
{
  if (matrixtype_ == Core::LinAlg::MatrixType::block_condition or
      matrixtype_ == Core::LinAlg::MatrixType::block_condition_dof)
  {
    // initialize reduced interface map extractors associated with blocks of global system matrix
    const int nblocks = scatratimint_->BlockMaps()->NumMaps();
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps_slave(nblocks);
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps_master(nblocks);
    for (int iblock = 0; iblock < nblocks; ++iblock)
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> maps(2);
      maps[0] = scatratimint_->BlockMaps()->Map(iblock);
      maps[1] = not imortarredistribution_
                    ? interfacemaps_->Map(1)
                    : Teuchos::rcp_dynamic_cast<const Epetra_Map>(islavemap_);
      blockmaps_slave[iblock] = Core::LinAlg::MultiMapExtractor::IntersectMaps(maps);
      maps[1] = not imortarredistribution_
                    ? interfacemaps_->Map(2)
                    : Teuchos::rcp_dynamic_cast<const Epetra_Map>(imastermap_);
      blockmaps_master[iblock] = Core::LinAlg::MultiMapExtractor::IntersectMaps(maps);
    }
    blockmaps_slave_ =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*interfacemaps_->Map(1), blockmaps_slave));
    blockmaps_slave_->check_for_valid_map_extractor();
    blockmaps_master_ = Teuchos::rcp(
        new Core::LinAlg::MultiMapExtractor(*interfacemaps_->Map(2), blockmaps_master));
    blockmaps_master_->check_for_valid_map_extractor();
  }
}  // ScaTra::meshtying_strategy_s2_i::build_block_map_extractors

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::equip_extended_solver_with_null_space_info() const
{
  // consider extended linear solver for scatra-scatra interface layer growth if applicable
  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
  {
    // loop over blocks of scalar transport system matrix
    for (int iblock = 0; iblock < scatratimint_->BlockMaps()->NumMaps(); ++iblock)
    {
      // store number of current block as string, starting from 1
      std::ostringstream iblockstr;
      iblockstr << iblock + 1;

      // equip smoother for current matrix block with previously computed null space
      extendedsolver_->Params().sublist("Inverse" + iblockstr.str()) =
          scatratimint_->Solver()->Params().sublist("Inverse" + iblockstr.str());
    }
    // store number of matrix block associated with scatra-scatra interface layer growth as string
    std::stringstream iblockstr;
    iblockstr << scatratimint_->BlockMaps()->NumMaps() + 1;

    // equip smoother for extra matrix block with null space associated with all degrees of freedom
    // for scatra-scatra interface layer growth
    Teuchos::ParameterList& mllist =
        extendedsolver_->Params().sublist("Inverse" + iblockstr.str()).sublist("MueLu Parameters");
    mllist.set("PDE equations", 1);
    mllist.set("null space: dimension", 1);
    mllist.set("null space: type", "pre-computed");
    mllist.set("null space: add default vectors", false);

    const Teuchos::RCP<Epetra_MultiVector> nullspace =
        Teuchos::rcp(new Epetra_MultiVector(*(scatratimint_->dof_row_map(2)), 1, true));
    nullspace->PutScalar(1.0);

    mllist.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspace);
    mllist.set("null space: vectors", nullspace->Values());
    mllist.set("ML validate parameter list", false);
  }
}  // ScaTra::meshtying_strategy_s2_i::build_block_null_spaces

/*------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::Solve(const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& increment, const Teuchos::RCP<Epetra_Vector>& residual,
    const Teuchos::RCP<Epetra_Vector>& phinp, const int iteration,
    Core::LinAlg::SolverParams& solver_params) const
{
  switch (intlayergrowth_evaluation_)
  {
    // no or semi-implicit treatment of scatra-scatra interface layer growth
    case Inpar::S2I::growth_evaluation_none:
    case Inpar::S2I::growth_evaluation_semi_implicit:
    {
      switch (couplingtype_)
      {
        case Inpar::S2I::coupling_matching_nodes:
        case Inpar::S2I::coupling_mortar_standard:
        case Inpar::S2I::coupling_mortar_condensed_petrov:
        case Inpar::S2I::coupling_mortar_condensed_bubnov:
        case Inpar::S2I::coupling_nts_standard:
        {
          // equilibrate global system of equations if necessary
          equilibration_->EquilibrateSystem(systemmatrix, residual, scatratimint_->BlockMaps());

          // solve global system of equations
          solver_params.refactor = true;
          solver_params.reset = iteration == 1;
          solver->Solve(systemmatrix->EpetraOperator(), increment, residual, solver_params);

          // unequilibrate global increment vector if necessary
          equilibration_->unequilibrate_increment(increment);

          break;
        }

        case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
        case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
        {
          // check scalar transport system matrix
          Teuchos::RCP<Core::LinAlg::SparseMatrix> sparsematrix =
              Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix);
          if (sparsematrix == Teuchos::null) FOUR_C_THROW("System matrix is not a sparse matrix!");

          // assemble extended system matrix including rows and columns associated with Lagrange
          // multipliers
          Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>
              extendedsystemmatrix(*extendedmaps_, *extendedmaps_);
          extendedsystemmatrix.Assign(0, 0, Core::LinAlg::View, *sparsematrix);
          if (lmside_ == Inpar::S2I::side_slave)
          {
            extendedsystemmatrix.Matrix(0, 1).Add(*D_, true, 1., 0.);
            extendedsystemmatrix.Matrix(0, 1).Add(*M_, true, -1., 1.);
            extendedsystemmatrix.Matrix(1, 0).Add(
                *Mortar::MatrixRowTransformGIDs(islavematrix_, extendedmaps_->Map(1)), false, 1.,
                0.);
          }
          else
          {
            extendedsystemmatrix.Matrix(0, 1).Add(*M_, true, -1., 0.);
            extendedsystemmatrix.Matrix(0, 1).Add(*D_, true, 1., 1.);
            extendedsystemmatrix.Matrix(1, 0).Add(
                *Mortar::MatrixRowTransformGIDs(imastermatrix_, extendedmaps_->Map(1)), false, 1.,
                0.);
          }
          extendedsystemmatrix.Matrix(1, 1).Add(*E_, true, -1., 0.);
          extendedsystemmatrix.Complete();
          extendedsystemmatrix.Matrix(0, 1).ApplyDirichlet(
              *scatratimint_->DirichMaps()->CondMap(), false);

          Teuchos::RCP<Epetra_Vector> extendedresidual =
              Core::LinAlg::CreateVector(*extendedmaps_->FullMap());
          extendedmaps_->InsertVector(scatratimint_->Residual(), 0, extendedresidual);
          extendedmaps_->InsertVector(lmresidual_, 1, extendedresidual);

          Teuchos::RCP<Epetra_Vector> extendedincrement =
              Core::LinAlg::CreateVector(*extendedmaps_->FullMap());
          extendedmaps_->InsertVector(scatratimint_->Increment(), 0, extendedincrement);
          extendedmaps_->InsertVector(lmincrement_, 1, extendedincrement);

          // solve extended system of equations
          solver_params.refactor = true;
          solver_params.reset = iteration == 1;
          solver->Solve(extendedsystemmatrix.EpetraOperator(), extendedincrement, extendedresidual,
              solver_params);

          // store solution
          extendedmaps_->ExtractVector(extendedincrement, 0, increment);
          extendedmaps_->ExtractVector(extendedincrement, 1, lmincrement_);

          // update Lagrange multipliers
          lm_->Update(1., *lmincrement_, 1.);

          break;
        }

        default:
        {
          FOUR_C_THROW("Type of scatra-scatra interface coupling not recognized!");
          break;
        }
      }

      break;
    }

    // monolithic treatment of scatra-scatra interface layer growth
    case Inpar::S2I::growth_evaluation_monolithic:
    {
      switch (couplingtype_)
      {
        case Inpar::S2I::coupling_matching_nodes:
        {
          switch (matrixtype_)
          {
            case Core::LinAlg::MatrixType::sparse:
            {
              // check scalar transport system matrix
              const Teuchos::RCP<const Core::LinAlg::SparseMatrix> sparsematrix =
                  Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix);
              if (sparsematrix == Teuchos::null)
                FOUR_C_THROW("System matrix is not a sparse matrix!");

              // assemble extended system matrix including rows and columns associated with
              // scatra-scatra interface layer thickness variables
              extendedsystemmatrix_->Assign(0, 0, Core::LinAlg::View, *sparsematrix);
              extendedsystemmatrix_->Assign(0, 1, Core::LinAlg::View,
                  *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(scatragrowthblock_));
              extendedsystemmatrix_->Assign(1, 0, Core::LinAlg::View,
                  *Teuchos::rcp_dynamic_cast<const Core::LinAlg::SparseMatrix>(growthscatrablock_));
              extendedsystemmatrix_->Assign(1, 1, Core::LinAlg::View, *growthgrowthblock_);

              break;
            }

            case Core::LinAlg::MatrixType::block_condition:
            case Core::LinAlg::MatrixType::block_condition_dof:
            {
              // check scalar transport system matrix
              const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> blocksparsematrix =
                  Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(systemmatrix);
              if (blocksparsematrix == Teuchos::null)
                FOUR_C_THROW("System matrix is not a block sparse matrix!");

              // extract number of matrix row or column blocks associated with scalar transport
              // field
              const int nblockmaps = static_cast<int>(scatratimint_->BlockMaps()->NumMaps());

              // construct extended system matrix by assigning matrix blocks
              for (int iblock = 0; iblock < nblockmaps; ++iblock)
              {
                for (int jblock = 0; jblock < nblockmaps; ++jblock)
                  extendedsystemmatrix_->Assign(iblock, jblock, Core::LinAlg::View,
                      blocksparsematrix->Matrix(iblock, jblock));
                extendedsystemmatrix_->Assign(iblock, nblockmaps, Core::LinAlg::View,
                    Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                        scatragrowthblock_)
                        ->Matrix(iblock, 0));
                extendedsystemmatrix_->Assign(nblockmaps, iblock, Core::LinAlg::View,
                    Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(
                        growthscatrablock_)
                        ->Matrix(0, iblock));
              }
              extendedsystemmatrix_->Assign(
                  nblockmaps, nblockmaps, Core::LinAlg::View, *growthgrowthblock_);

              break;
            }

            default:
            {
              FOUR_C_THROW(
                  "Type of global system matrix for scatra-scatra interface coupling involving "
                  "interface layer growth not recognized!");
              break;
            }
          }

          // finalize extended system matrix
          extendedsystemmatrix_->Complete();

          // assemble extended residual vector
          Teuchos::RCP<Epetra_Vector> extendedresidual =
              Teuchos::rcp(new Epetra_Vector(*extendedmaps_->FullMap(), true));
          extendedmaps_->InsertVector(scatratimint_->Residual(), 0, extendedresidual);
          extendedmaps_->InsertVector(growthresidual_, 1, extendedresidual);

          // perform finite-difference check if desired
          if (scatratimint_->FDCheckType() == Inpar::ScaTra::fdcheck_global_extended)
            fd_check(*extendedsystemmatrix_, extendedresidual);

          // assemble extended increment vector
          Teuchos::RCP<Epetra_Vector> extendedincrement =
              Teuchos::rcp(new Epetra_Vector(*extendedmaps_->FullMap(), true));
          extendedmaps_->InsertVector(scatratimint_->Increment(), 0, extendedincrement);
          extendedmaps_->InsertVector(growthincrement_, 1, extendedincrement);

          // equilibrate global system of equations if necessary
          equilibration_->EquilibrateSystem(
              extendedsystemmatrix_, extendedresidual, extendedblockmaps_);

          // solve extended system of equations
          solver_params.refactor = true;
          solver_params.reset = iteration == 1;
          extendedsolver_->Solve(extendedsystemmatrix_->EpetraOperator(), extendedincrement,
              extendedresidual, solver_params);

          // unequilibrate global increment vector if necessary
          equilibration_->unequilibrate_increment(extendedincrement);

          // store solution
          extendedmaps_->ExtractVector(extendedincrement, 0, increment);
          extendedmaps_->ExtractVector(extendedincrement, 1, growthincrement_);

          // update state vector of discrete scatra-scatra interface layer thicknesses
          growthnp_->Update(1., *growthincrement_, 1.);

          break;
        }

        default:
        {
          FOUR_C_THROW("Type of scatra-scatra interface coupling not recognized!");
          break;
        }
      }  // switch(couplingtype_)

      break;
    }

    default:
    {
      FOUR_C_THROW(
          "Unknown evaluation method for scatra-scatra interface coupling involving interface "
          "layer growth!");
      break;
    }
  }
}  // ScaTra::meshtying_strategy_s2_i::Solve

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Core::LinAlg::Solver& ScaTra::MeshtyingStrategyS2I::Solver() const
{
  const Core::LinAlg::Solver* solver(nullptr);

  if (intlayergrowth_evaluation_ == Inpar::S2I::growth_evaluation_monolithic)
  {
    if (extendedsolver_ == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");
    solver = extendedsolver_.get();
  }

  else
  {
    if (scatratimint_->Solver() == Teuchos::null) FOUR_C_THROW("Invalid linear solver!");
    solver = scatratimint_->Solver().get();
  }

  return *solver;
}  // ScaTra::meshtying_strategy_s2_i::Solver()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MeshtyingStrategyS2I::fd_check(
    const Core::LinAlg::BlockSparseMatrixBase& extendedsystemmatrix,
    const Teuchos::RCP<Epetra_Vector>& extendedresidual) const
{
  // initial screen output
  if (scatratimint_->discretization()->Comm().MyPID() == 0)
  {
    std::cout << std::endl
              << "FINITE DIFFERENCE CHECK FOR EXTENDED SYSTEM MATRIX INVOLVING SCATRA-SCATRA "
                 "INTERFACE LAYER GROWTH"
              << std::endl;
  }

  // extract perturbation magnitude and relative tolerance
  const double fdcheckeps(scatratimint_->FDCheckEps());
  const double fdchecktol(scatratimint_->FDCheckTol());

  // create global state vector
  Epetra_Vector statenp(*extendedmaps_->FullMap(), true);
  extendedmaps_->InsertVector(*scatratimint_->Phinp(), 0, statenp);
  extendedmaps_->InsertVector(*growthnp_, 1, statenp);

  // make a copy of global state vector to undo perturbations later
  Epetra_Vector statenp_original(statenp);

  // make a copy of system matrix as Epetra_CrsMatrix
  Epetra_CrsMatrix sysmat_original =
      *Core::LinAlg::SparseMatrix(*extendedsystemmatrix.Merge()).EpetraMatrix();
  sysmat_original.FillComplete();

  // make a copy of system right-hand side vector
  Epetra_Vector rhs_original(*extendedresidual);

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over all columns of system matrix
  for (int colgid = 0; colgid <= sysmat_original.ColMap().MaxAllGID(); ++colgid)
  {
    // check whether current column index is a valid global column index and continue loop if not
    int collid(sysmat_original.ColMap().LID(colgid));
    int maxcollid(-1);
    scatratimint_->discretization()->Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill global state vector with original state variables
    statenp.Update(1., statenp_original, 0.);

    // impose perturbation
    if (statenp.Map().MyGID(colgid))
      if (statenp.SumIntoGlobalValue(colgid, 0, fdcheckeps))
        FOUR_C_THROW(
            "Perturbation could not be imposed on state vector for finite difference check!");
    scatratimint_->Phinp()->Update(1., *extendedmaps_->ExtractVector(statenp, 0), 0.);
    growthnp_->Update(1., *extendedmaps_->ExtractVector(statenp, 1), 0.);

    // calculate global right-hand side contributions based on perturbed state
    scatratimint_->assemble_mat_and_rhs();

    // assemble global residual vector
    extendedmaps_->InsertVector(scatratimint_->Residual(), 0, extendedresidual);
    extendedmaps_->InsertVector(growthresidual_, 1, extendedresidual);

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the
    // system matrix, the second comparison might yield good agreement in spite of the entries being
    // wrong!
    for (int rowlid = 0; rowlid < extendedmaps_->FullMap()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original.RowMap().GID(rowlid);
      if (rowgid < 0) FOUR_C_THROW("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original.NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original.ExtractMyRowCopy(rowlid, length, numentries, values.data(), indices.data());
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original.ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval =
          -(*extendedresidual)[rowlid] / fdcheckeps + rhs_original[rowlid] / fdcheckeps;

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if (abs(abserr1) > maxabserr) maxabserr = abs(abserr1);
      double relerr1(0.);
      if (abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if (abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if (abs(relerr1) > maxrelerr) maxrelerr = abs(relerr1);

      // evaluate first comparison
      if (abs(relerr1) > fdchecktol)
      {
        std::cout << std::setprecision(6);
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left = entry - rhs_original[rowlid] / fdcheckeps;

        // right-hand side in second comparison
        const double right = -(*extendedresidual)[rowlid] / fdcheckeps;

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if (abs(abserr2) > maxabserr) maxabserr = abs(abserr2);
        double relerr2(0.);
        if (abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if (abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if (abs(relerr2) > maxrelerr) maxrelerr = abs(relerr2);

        // evaluate second comparison
        if (abs(relerr2) > fdchecktol)
        {
          std::cout << std::setprecision(6);
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid
                    << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  scatratimint_->discretization()->Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  scatratimint_->discretization()->Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  scatratimint_->discretization()->Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (scatratimint_->discretization()->Comm().MyPID() == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      FOUR_C_THROW(
          "Finite difference check failed for extended system matrix involving scatra-scatra "
          "interface layer growth!");
    }
    else
    {
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
    }
  }

  // undo perturbations of state variables
  scatratimint_->Phinp()->Update(1., *extendedmaps_->ExtractVector(statenp_original, 0), 0.);
  growthnp_->Update(1., *extendedmaps_->ExtractVector(statenp_original, 1), 0.);

  // recompute system matrix and right-hand side vector based on original state variables
  scatratimint_->assemble_mat_and_rhs();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::MortarCellInterface::MortarCellInterface(const Inpar::S2I::CouplingType& couplingtype,
    const Inpar::S2I::InterfaceSides& lmside, const int& numdofpernode_slave,
    const int& numdofpernode_master)
    : lmside_(lmside),
      couplingtype_(couplingtype),
      numdofpernode_slave_(numdofpernode_slave),
      numdofpernode_master_(numdofpernode_master)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
ScaTra::MortarCellCalc<distypeS, distypeM>* ScaTra::MortarCellCalc<distypeS, distypeM>::Instance(
    const Inpar::S2I::CouplingType& couplingtype, const Inpar::S2I::InterfaceSides& lmside,
    const int& numdofpernode_slave, const int& numdofpernode_master, const std::string& disname)
{
  // static map assigning mortar discretization names to class instances
  static std::map<std::string, Core::UTILS::SingletonOwner<MortarCellCalc<distypeS, distypeM>>>
      owners;

  // add an owner for the given disname if not already present
  if (owners.find(disname) == owners.end())
  {
    owners.template emplace(
        disname, Core::UTILS::SingletonOwner<MortarCellCalc<distypeS, distypeM>>(
                     [&]()
                     {
                       return std::unique_ptr<MortarCellCalc<distypeS, distypeM>>(
                           new MortarCellCalc<distypeS, distypeM>(
                               couplingtype, lmside, numdofpernode_slave, numdofpernode_master));
                     }));
  }

  return owners.find(disname)->second.Instance(Core::UTILS::SingletonAction::create);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::Evaluate(const Core::FE::Discretization& idiscret,
    Mortar::IntCell& cell, Mortar::Element& slaveelement, Mortar::Element& masterelement,
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& cellmatrix1, Core::LinAlg::SerialDenseMatrix& cellmatrix2,
    Core::LinAlg::SerialDenseMatrix& cellmatrix3, Core::LinAlg::SerialDenseMatrix& cellmatrix4,
    Core::LinAlg::SerialDenseVector& cellvector1, Core::LinAlg::SerialDenseVector& cellvector2)
{
  // extract and evaluate action
  switch (Core::UTILS::GetAsEnum<Inpar::S2I::EvaluationActions>(params, "action"))
  {
    case Inpar::S2I::evaluate_mortar_matrices:
    {
      // evaluate mortar matrices
      evaluate_mortar_matrices(
          cell, slaveelement, masterelement, cellmatrix1, cellmatrix2, cellmatrix3);

      break;
    }

    case Inpar::S2I::evaluate_condition:
    {
      // evaluate and assemble interface linearizations and residuals
      evaluate_condition(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown action for mortar cell evaluation!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_nts(
    const Core::FE::Discretization& idiscret, const Mortar::Node& slavenode,
    const double& lumpedarea, Mortar::Element& slaveelement, Mortar::Element& masterelement,
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& ntsmatrix1, Core::LinAlg::SerialDenseMatrix& ntsmatrix2,
    Core::LinAlg::SerialDenseMatrix& ntsmatrix3, Core::LinAlg::SerialDenseMatrix& ntsmatrix4,
    Core::LinAlg::SerialDenseVector& ntsvector1, Core::LinAlg::SerialDenseVector& ntsvector2)
{
  // extract and evaluate action
  switch (Core::UTILS::GetAsEnum<Inpar::S2I::EvaluationActions>(params, "action"))
  {
    case Inpar::S2I::evaluate_condition_nts:
    {
      // extract condition from parameter list
      Core::Conditions::Condition* condition =
          params.get<Core::Conditions::Condition*>("condition");
      if (condition == nullptr)
        FOUR_C_THROW("Cannot access scatra-scatra interface coupling condition!");

      // extract nodal state variables associated with slave and master elements
      extract_node_values(idiscret, la_slave, la_master);

      // evaluate and assemble interface linearizations and residuals
      evaluate_condition_nts(*condition, slavenode, lumpedarea, slaveelement, masterelement,
          ephinp_slave_, ephinp_master_, ntsmatrix1, ntsmatrix2, ntsmatrix3, ntsmatrix4, ntsvector1,
          ntsvector2);

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown action for evaluation of node-to-segment coupling!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_mortar_element(
    const Core::FE::Discretization& idiscret, Mortar::Element& element,
    Core::Elements::Element::LocationArray& la, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& elematrix1, Core::LinAlg::SerialDenseMatrix& elematrix2,
    Core::LinAlg::SerialDenseMatrix& elematrix3, Core::LinAlg::SerialDenseMatrix& elematrix4,
    Core::LinAlg::SerialDenseVector& elevector1, Core::LinAlg::SerialDenseVector& elevector2)
{
  // extract and evaluate action
  switch (Core::UTILS::GetAsEnum<Inpar::S2I::EvaluationActions>(params, "action"))
  {
    case Inpar::S2I::evaluate_nodal_area_fractions:
    {
      // evaluate and assemble lumped interface area fractions associated with element nodes
      evaluate_nodal_area_fractions(element, elevector1);

      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown action for evaluation of mortar element!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
ScaTra::MortarCellCalc<distypeS, distypeM>::MortarCellCalc(
    const Inpar::S2I::CouplingType& couplingtype, const Inpar::S2I::InterfaceSides& lmside,
    const int& numdofpernode_slave, const int& numdofpernode_master)
    : MortarCellInterface(couplingtype, lmside, numdofpernode_slave, numdofpernode_master),
      scatraparamsboundary_(Discret::ELEMENTS::ScaTraEleParameterBoundary::Instance("scatra")),
      ephinp_slave_(numdofpernode_slave, Core::LinAlg::Matrix<nen_slave_, 1>(true)),
      ephinp_master_(numdofpernode_master, Core::LinAlg::Matrix<nen_master_, 1>(true)),
      funct_slave_(true),
      funct_master_(true),
      shape_lm_slave_(true),
      shape_lm_master_(true),
      test_lm_slave_(true),
      test_lm_master_(true)
{
  // safety check
  if (nsd_slave_ != 2 or nsd_master_ != 2)
  {
    FOUR_C_THROW(
        "Scatra-scatra interface coupling with non-matching interface discretization currently "
        "only implemented for two-dimensional interface manifolds!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::extract_node_values(
    const Core::FE::Discretization& idiscret, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master)
{
  // extract nodal state variables associated with mortar integration cell
  extract_node_values(ephinp_slave_, ephinp_master_, idiscret, la_slave, la_master);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::extract_node_values(
    Core::LinAlg::Matrix<nen_slave_, 1>& estate_slave, const Core::FE::Discretization& idiscret,
    Core::Elements::Element::LocationArray& la_slave, const std::string& statename,
    const int& nds) const
{
  // extract interface state vector from interface discretization
  const Teuchos::RCP<const Epetra_Vector> state = idiscret.GetState(nds, statename);
  if (state == Teuchos::null)
    FOUR_C_THROW(
        "Cannot extract state vector \"" + statename + "\" from interface discretization!");

  // extract nodal state variables associated with slave element
  Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_slave_, 1>>(
      *state, estate_slave, la_slave[nds].lm_);
}


/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::extract_node_values(
    std::vector<Core::LinAlg::Matrix<nen_slave_, 1>>& estate_slave,
    std::vector<Core::LinAlg::Matrix<nen_master_, 1>>& estate_master,
    const Core::FE::Discretization& idiscret, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const std::string& statename,
    const int& nds) const
{
  // extract interface state vector from interface discretization
  const Teuchos::RCP<const Epetra_Vector> state = idiscret.GetState(nds, statename);
  if (state == Teuchos::null)
    FOUR_C_THROW(
        "Cannot extract state vector \"" + statename + "\" from interface discretization!");

  // extract nodal state variables associated with slave and master elements
  Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_slave_, 1>>(
      *state, estate_slave, la_slave[nds].lm_);
  Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_master_, 1>>(
      *state, estate_master, la_master[nds].lm_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
double ScaTra::MortarCellCalc<distypeS, distypeM>::eval_shape_func_and_dom_int_fac_at_int_point(
    Mortar::Element& slaveelement, Mortar::Element& masterelement, Mortar::IntCell& cell,
    const Core::FE::IntPointsAndWeights<nsd_slave_>& intpoints, const int iquad)
{
  // reference coordinates of integration point
  std::array<double, nsd_slave_> coordinates_ref;
  for (int idim = 0; idim < nsd_slave_; ++idim)
    coordinates_ref[idim] = intpoints.IP().qxg[iquad][idim];

  // global coordinates of integration point
  std::array<double, nsd_slave_ + 1> coordinates_global;
  cell.LocalToGlobal(coordinates_ref.data(), coordinates_global.data(), 0);

  // project integration point onto slave and master elements
  std::array<double, nsd_slave_> coordinates_slave;
  std::array<double, nsd_master_> coordinates_master;
  double dummy(0.);
  Mortar::Projector::Impl(slaveelement)
      ->project_gauss_point_auxn3_d(
          coordinates_global.data(), cell.Auxn(), slaveelement, coordinates_slave.data(), dummy);
  Mortar::Projector::Impl(masterelement)
      ->project_gauss_point_auxn3_d(
          coordinates_global.data(), cell.Auxn(), masterelement, coordinates_master.data(), dummy);

  // evaluate shape functions at current integration point on slave and master elements
  Core::VolMortar::UTILS::shape_function<distypeS>(funct_slave_, coordinates_slave.data());
  Core::VolMortar::UTILS::shape_function<distypeM>(funct_master_, coordinates_master.data());
  switch (couplingtype_)
  {
    case Inpar::S2I::coupling_mortar_standard:
    {
      // there actually aren't any Lagrange multipliers, but we still need to set pseudo Lagrange
      // multiplier test functions equal to the standard shape and test functions for correct
      // evaluation of the scatra-scatra interface coupling conditions
      test_lm_slave_ = funct_slave_;
      test_lm_master_ = funct_master_;

      break;
    }

    case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
    case Inpar::S2I::coupling_mortar_condensed_petrov:
    {
      // dual Lagrange multiplier shape functions combined with standard Lagrange multiplier test
      // functions
      if (lmside_ == Inpar::S2I::side_slave)
      {
        Core::VolMortar::UTILS::dual_shape_function<distypeS>(
            shape_lm_slave_, coordinates_slave.data(), slaveelement);
        test_lm_slave_ = funct_slave_;
      }
      else
      {
        Core::VolMortar::UTILS::dual_shape_function<distypeM>(
            shape_lm_master_, coordinates_master.data(), masterelement);
        test_lm_master_ = funct_master_;
      }

      break;
    }

    case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
    case Inpar::S2I::coupling_mortar_condensed_bubnov:
    {
      // dual Lagrange multiplier shape functions combined with dual Lagrange multiplier test
      // functions
      if (lmside_ == Inpar::S2I::side_slave)
      {
        Core::VolMortar::UTILS::dual_shape_function<distypeS>(
            shape_lm_slave_, coordinates_slave.data(), slaveelement);
        test_lm_slave_ = shape_lm_slave_;
      }
      else
      {
        Core::VolMortar::UTILS::dual_shape_function<distypeM>(
            shape_lm_master_, coordinates_master.data(), masterelement);
        test_lm_master_ = shape_lm_master_;
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Not yet implemented!");
      break;
    }
  }

  // integration weight
  const double weight = intpoints.IP().qwgt[iquad];

  // Jacobian determinant
  const double jacobian = cell.Jacobian();

  // domain integration factor
  return jacobian * weight;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
double ScaTra::MortarCellCalc<distypeS, distypeM>::eval_shape_func_and_dom_int_fac_at_int_point(
    Mortar::Element& element, const Core::FE::IntPointsAndWeights<nsd_slave_>& intpoints,
    const int iquad)
{
  // extract global coordinates of element nodes
  Core::LinAlg::Matrix<nsd_slave_ + 1, nen_slave_> coordinates_nodes;
  Core::Geo::fillInitialPositionArray<distypeS, nsd_slave_ + 1,
      Core::LinAlg::Matrix<nsd_slave_ + 1, nen_slave_>>(&element, coordinates_nodes);

  // extract reference coordinates of integration point
  Core::LinAlg::Matrix<nsd_slave_, 1> coordinates_ref(intpoints.IP().qxg[iquad]);

  // evaluate slave-side shape functions and their first derivatives at integration point
  Core::LinAlg::Matrix<nsd_slave_, nen_slave_> deriv_slave;
  Core::FE::shape_function<distypeS>(coordinates_ref, funct_slave_);
  Core::FE::shape_function_deriv1<distypeS>(coordinates_ref, deriv_slave);

  // evaluate transposed Jacobian matrix at integration point
  Core::LinAlg::Matrix<nsd_slave_, nsd_slave_ + 1> jacobian;
  jacobian.MultiplyNT(deriv_slave, coordinates_nodes);

  // evaluate metric tensor at integration point
  Core::LinAlg::Matrix<nsd_slave_, nsd_slave_> metrictensor;
  metrictensor.MultiplyNT(jacobian, jacobian);

  // return domain integration factor, i.e., Jacobian determinant times integration weight, at
  // integration point
  return sqrt(metrictensor.Determinant()) * intpoints.IP().qwgt[iquad];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::eval_shape_func_at_slave_node(
    const Mortar::Node& slavenode, Mortar::Element& slaveelement, Mortar::Element& masterelement)
{
  // safety check
  if (couplingtype_ != Inpar::S2I::coupling_nts_standard)
    FOUR_C_THROW("This function should only be called when evaluating node-to-segment coupling!");

  // extract global ID of slave-side node
  const int& nodeid = slavenode.Id();

  // find out index of slave-side node w.r.t. slave-side element
  int index(-1);
  for (int inode = 0; inode < slaveelement.num_node(); ++inode)
  {
    if (nodeid == slaveelement.Nodes()[inode]->Id())
    {
      index = inode;
      break;
    }
  }
  if (index == -1)
    FOUR_C_THROW("Couldn't find out index of slave-side node w.r.t. slave-side element!");

  // set slave-side shape function array according to node position
  funct_slave_.Clear();
  funct_slave_(index) = 1.;

  // project slave-side node onto master-side element
  std::array<double, 2> coordinates_master;
  double dummy(0.);
  Mortar::Projector::Impl(masterelement)
      ->project_gauss_point_auxn3_d(slavenode.X().data(), slavenode.MoData().n(), masterelement,
          coordinates_master.data(), dummy);

  // evaluate master-side shape functions at projected node on master-side element
  Core::VolMortar::UTILS::shape_function<distypeM>(funct_master_, coordinates_master.data());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_mortar_matrices(Mortar::IntCell& cell,
    Mortar::Element& slaveelement, Mortar::Element& masterelement,
    Core::LinAlg::SerialDenseMatrix& D, Core::LinAlg::SerialDenseMatrix& M,
    Core::LinAlg::SerialDenseMatrix& E)
{
  // safety check
  if (numdofpernode_slave_ != numdofpernode_master_)
    FOUR_C_THROW("Must have same number of degrees of freedom per node on slave and master sides!");

  // determine quadrature rule
  const Core::FE::IntPointsAndWeights<2> intpoints(Core::FE::GaussRule2D::tri_7point);

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = eval_shape_func_and_dom_int_fac_at_int_point(
        slaveelement, masterelement, cell, intpoints, iquad);

    if (lmside_ == Inpar::S2I::side_slave)
    {
      // loop over all degrees of freedom per node
      for (int k = 0; k < numdofpernode_slave_; ++k)
      {
        for (int vi = 0; vi < nen_slave_; ++vi)
        {
          const int row_slave = vi * numdofpernode_slave_ + k;

          switch (couplingtype_)
          {
            case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
            case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
            case Inpar::S2I::coupling_mortar_condensed_petrov:
            case Inpar::S2I::coupling_mortar_condensed_bubnov:
            {
              D(row_slave, row_slave) += shape_lm_slave_(vi) * fac;

              if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov or
                  couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
              {
                for (int ui = 0; ui < nen_slave_; ++ui)
                  E(row_slave, ui * numdofpernode_slave_ + k) +=
                      shape_lm_slave_(vi) * test_lm_slave_(ui) * fac;
              }

              break;
            }

            default:
            {
              for (int ui = 0; ui < nen_slave_; ++ui)
                D(row_slave, ui * numdofpernode_slave_ + k) +=
                    shape_lm_slave_(vi) * funct_slave_(ui) * fac;

              break;
            }
          }

          for (int ui = 0; ui < nen_master_; ++ui)
            M(row_slave, ui * numdofpernode_master_ + k) +=
                shape_lm_slave_(vi) * funct_master_(ui) * fac;
        }
      }
    }

    else
    {
      // loop over all degrees of freedom per node
      for (int k = 0; k < numdofpernode_master_; ++k)
      {
        for (int vi = 0; vi < nen_master_; ++vi)
        {
          const int row_master = vi * numdofpernode_master_ + k;

          switch (couplingtype_)
          {
            case Inpar::S2I::coupling_mortar_saddlepoint_petrov:
            case Inpar::S2I::coupling_mortar_saddlepoint_bubnov:
            case Inpar::S2I::coupling_mortar_condensed_petrov:
            case Inpar::S2I::coupling_mortar_condensed_bubnov:
            {
              D(row_master, row_master) += shape_lm_master_(vi) * fac;

              if (couplingtype_ == Inpar::S2I::coupling_mortar_saddlepoint_bubnov or
                  couplingtype_ == Inpar::S2I::coupling_mortar_condensed_bubnov)
              {
                for (int ui = 0; ui < nen_master_; ++ui)
                  E(row_master, ui * numdofpernode_master_ + k) +=
                      shape_lm_master_(vi) * test_lm_master_(ui) * fac;
              }

              break;
            }

            default:
            {
              for (int ui = 0; ui < nen_master_; ++ui)
                D(row_master, ui * numdofpernode_master_ + k) +=
                    shape_lm_master_(vi) * funct_master_(ui) * fac;

              break;
            }
          }

          for (int ui = 0; ui < nen_slave_; ++ui)
            M(row_master, ui * numdofpernode_slave_ + k) +=
                shape_lm_master_(vi) * funct_slave_(ui) * fac;
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_condition(
    const Core::FE::Discretization& idiscret, Mortar::IntCell& cell, Mortar::Element& slaveelement,
    Mortar::Element& masterelement, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseMatrix& k_ss, Core::LinAlg::SerialDenseMatrix& k_sm,
    Core::LinAlg::SerialDenseMatrix& k_ms, Core::LinAlg::SerialDenseMatrix& k_mm,
    Core::LinAlg::SerialDenseVector& r_s, Core::LinAlg::SerialDenseVector& r_m)
{
  // extract nodal state variables associated with slave and master elements
  extract_node_values(idiscret, la_slave, la_master);

  // safety check
  if (numdofpernode_slave_ != 1 or numdofpernode_master_ != 1)
  {
    FOUR_C_THROW(
        "Invalid number of degrees of freedom per node! Code should theoretically work for more "
        "than one degree of freedom per node, but not yet tested!");
  }

  // always in contact
  const double pseudo_contact_fac = 1.0;

  // determine quadrature rule
  const Core::FE::IntPointsAndWeights<2> intpoints(Core::FE::GaussRule2D::tri_7point);

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = eval_shape_func_and_dom_int_fac_at_int_point(
        slaveelement, masterelement, cell, intpoints, iquad);

    // overall integration factors
    const double timefacfac =
        Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * fac;
    const double timefacrhsfac =
        Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) FOUR_C_THROW("Integration factor is negative!");

    Discret::ELEMENTS::ScaTraEleBoundaryCalc<
        distypeS>::template evaluate_s2_i_coupling_at_integration_point<distypeM>(ephinp_slave_,
        ephinp_master_, pseudo_contact_fac, funct_slave_, funct_master_, test_lm_slave_,
        test_lm_master_, numdofpernode_slave_, scatraparamsboundary_, timefacfac, timefacrhsfac,
        k_ss, k_sm, k_ms, k_mm, r_s, r_m);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_condition_nts(
    Core::Conditions::Condition& condition, const Mortar::Node& slavenode, const double& lumpedarea,
    Mortar::Element& slaveelement, Mortar::Element& masterelement,
    const std::vector<Core::LinAlg::Matrix<nen_slave_, 1>>& ephinp_slave,
    const std::vector<Core::LinAlg::Matrix<nen_master_, 1>>& ephinp_master,
    Core::LinAlg::SerialDenseMatrix& k_ss, Core::LinAlg::SerialDenseMatrix& k_sm,
    Core::LinAlg::SerialDenseMatrix& k_ms, Core::LinAlg::SerialDenseMatrix& k_mm,
    Core::LinAlg::SerialDenseVector& r_s, Core::LinAlg::SerialDenseVector& r_m)
{
  // safety check
  if (numdofpernode_slave_ != 1 or numdofpernode_master_ != 1)
  {
    FOUR_C_THROW(
        "Invalid number of degrees of freedom per node! Code should theoretically work for more "
        "than one degree of freedom per node, but not yet tested!");
  }

  // evaluate shape functions at position of slave-side node
  eval_shape_func_at_slave_node(slavenode, slaveelement, masterelement);

  // always in contact
  const double pseudo_contact_fac = 1.0;

  // overall integration factors
  const double timefacfac =
      Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * lumpedarea;
  const double timefacrhsfac =
      Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * lumpedarea;
  if (timefacfac < 0. or timefacrhsfac < 0.) FOUR_C_THROW("Integration factor is negative!");

  Discret::ELEMENTS::ScaTraEleBoundaryCalc<
      distypeS>::template evaluate_s2_i_coupling_at_integration_point<distypeM>(ephinp_slave,
      ephinp_master, pseudo_contact_fac, funct_slave_, funct_master_, funct_slave_, funct_master_,
      numdofpernode_slave_, scatraparamsboundary_, timefacfac, timefacrhsfac, k_ss, k_sm, k_ms,
      k_mm, r_s, r_m);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distypeS, Core::FE::CellType distypeM>
void ScaTra::MortarCellCalc<distypeS, distypeM>::evaluate_nodal_area_fractions(
    Mortar::Element& slaveelement, Core::LinAlg::SerialDenseVector& areafractions)
{
  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_slave_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distypeS>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = eval_shape_func_and_dom_int_fac_at_int_point(slaveelement, intpoints, iquad);

    // compute integrals of shape functions to obtain lumped interface area fractions associated
    // with element nodes
    for (int inode = 0; inode < nen_slave_; ++inode)
      areafractions[inode * numdofpernode_slave_] += funct_slave_(inode) * fac;
  }  // loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::MortarCellAssemblyStrategy::MortarCellAssemblyStrategy(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix1,
    const Inpar::S2I::InterfaceSides matrix1_side_rows,
    const Inpar::S2I::InterfaceSides matrix1_side_cols,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix2,
    const Inpar::S2I::InterfaceSides matrix2_side_rows,
    const Inpar::S2I::InterfaceSides matrix2_side_cols,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix3,
    const Inpar::S2I::InterfaceSides matrix3_side_rows,
    const Inpar::S2I::InterfaceSides matrix3_side_cols,
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix4,
    const Inpar::S2I::InterfaceSides matrix4_side_rows,
    const Inpar::S2I::InterfaceSides matrix4_side_cols,
    Teuchos::RCP<Epetra_MultiVector> systemvector1, const Inpar::S2I::InterfaceSides vector1_side,
    Teuchos::RCP<Epetra_MultiVector> systemvector2, const Inpar::S2I::InterfaceSides vector2_side,
    const int nds_rows, const int nds_cols)
    : matrix1_side_rows_(matrix1_side_rows),
      matrix1_side_cols_(matrix1_side_cols),
      matrix2_side_rows_(matrix2_side_rows),
      matrix2_side_cols_(matrix2_side_cols),
      matrix3_side_rows_(matrix3_side_rows),
      matrix3_side_cols_(matrix3_side_cols),
      matrix4_side_rows_(matrix4_side_rows),
      matrix4_side_cols_(matrix4_side_cols),
      systemmatrix1_(std::move(systemmatrix1)),
      systemmatrix2_(std::move(systemmatrix2)),
      systemmatrix3_(std::move(systemmatrix3)),
      systemmatrix4_(std::move(systemmatrix4)),
      systemvector1_(std::move(systemvector1)),
      systemvector2_(std::move(systemvector2)),
      vector1_side_(vector1_side),
      vector2_side_(vector2_side),
      nds_rows_(nds_rows),
      nds_cols_(nds_cols)
{
}


/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::assemble_cell_matrices_and_vectors(
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const int assembler_pid_master) const
{
  // assemble cell matrix 1 into system matrix 1
  if (AssembleMatrix1())
    assemble_cell_matrix(systemmatrix1_, cellmatrix1_, matrix1_side_rows_, matrix1_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 2 into system matrix 2
  if (AssembleMatrix2())
    assemble_cell_matrix(systemmatrix2_, cellmatrix2_, matrix2_side_rows_, matrix2_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 3 into system matrix 3
  if (AssembleMatrix3())
    assemble_cell_matrix(systemmatrix3_, cellmatrix3_, matrix3_side_rows_, matrix3_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 4 into system matrix 4
  if (AssembleMatrix4())
    assemble_cell_matrix(systemmatrix4_, cellmatrix4_, matrix4_side_rows_, matrix4_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell vector 1 into system vector 1
  if (AssembleVector1())
    assemble_cell_vector(
        systemvector1_, cellvector1_, vector1_side_, la_slave, la_master, assembler_pid_master);

  // assemble cell vector 2 into system vector 2
  if (AssembleVector2())
    assemble_cell_vector(
        systemvector2_, cellvector2_, vector2_side_, la_slave, la_master, assembler_pid_master);
}


/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::assemble_cell_matrix(
    const Teuchos::RCP<Core::LinAlg::SparseOperator>& systemmatrix,
    const Core::LinAlg::SerialDenseMatrix& cellmatrix, const Inpar::S2I::InterfaceSides side_rows,
    const Inpar::S2I::InterfaceSides side_cols, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const int assembler_pid_master) const
{
  // determine location array associated with matrix columns
  Core::Elements::Element::LocationArray& la_cols =
      side_cols == Inpar::S2I::side_slave ? la_slave : la_master;

  // assemble cell matrix into system matrix
  switch (side_rows)
  {
    case Inpar::S2I::side_slave:
    {
      systemmatrix->Assemble(-1, la_cols[nds_cols_].stride_, cellmatrix, la_slave[nds_rows_].lm_,
          la_slave[nds_rows_].lmowner_, la_cols[nds_cols_].lm_);
      break;
    }

    case Inpar::S2I::side_master:
    {
      Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(systemmatrix)
          ->FEAssemble(cellmatrix, la_master[nds_rows_].lm_,
              std::vector<int>(la_master[nds_rows_].lmowner_.size(), assembler_pid_master),
              la_cols[nds_cols_].lm_);
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid interface side!");
      break;
    }
  }
}


/*----------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::assemble_cell_vector(
    const Teuchos::RCP<Epetra_MultiVector>& systemvector,
    const Core::LinAlg::SerialDenseVector& cellvector, const Inpar::S2I::InterfaceSides side,
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master, const int assembler_pid_master) const
{
  // assemble cell vector into system vector
  switch (side)
  {
    case Inpar::S2I::side_slave:
    {
      if (systemvector->NumVectors() != 1)
        FOUR_C_THROW("Invalid number of vectors inside Epetra_MultiVector!");
      Core::LinAlg::Assemble(*(*systemvector)(nds_rows_), cellvector, la_slave[nds_rows_].lm_,
          la_slave[nds_rows_].lmowner_);

      break;
    }

    case Inpar::S2I::side_master:
    {
      if (assembler_pid_master == systemvector->Comm().MyPID())
      {
        if (Teuchos::rcp_dynamic_cast<Epetra_FEVector>(systemvector)
                ->SumIntoGlobalValues(static_cast<int>(la_master[nds_rows_].lm_.size()),
                    la_master[nds_rows_].lm_.data(), cellvector.values()))
          FOUR_C_THROW("Assembly into master-side system vector not successful!");
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid interface side!");
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::init_cell_matrices_and_vectors(
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master)
{
  // initialize system matrix 1
  if (AssembleMatrix1())
    init_cell_matrix(cellmatrix1_, matrix1_side_rows_, matrix1_side_cols_, la_slave, la_master);

  // initialize system matrix 2
  if (AssembleMatrix2())
    init_cell_matrix(cellmatrix2_, matrix2_side_rows_, matrix2_side_cols_, la_slave, la_master);

  // initialize system matrix 3
  if (AssembleMatrix3())
    init_cell_matrix(cellmatrix3_, matrix3_side_rows_, matrix3_side_cols_, la_slave, la_master);

  // initialize system matrix 4
  if (AssembleMatrix4())
    init_cell_matrix(cellmatrix4_, matrix4_side_rows_, matrix4_side_cols_, la_slave, la_master);

  // initialize system vector 1
  if (AssembleVector1()) init_cell_vector(cellvector1_, vector1_side_, la_slave, la_master);

  // initialize system vector 2
  if (AssembleVector2()) init_cell_vector(cellvector2_, vector2_side_, la_slave, la_master);
}


/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::init_cell_matrix(
    Core::LinAlg::SerialDenseMatrix& cellmatrix, const Inpar::S2I::InterfaceSides side_rows,
    const Inpar::S2I::InterfaceSides side_cols, Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master) const
{
  // determine number of matrix rows and number of matrix columns
  const int nrows = side_rows == Inpar::S2I::side_slave ? la_slave[nds_rows_].Size()
                                                        : la_master[nds_rows_].Size();
  const int ncols = side_cols == Inpar::S2I::side_slave ? la_slave[nds_cols_].Size()
                                                        : la_master[nds_cols_].Size();

  // reshape cell matrix if necessary
  if (cellmatrix.numRows() != nrows or cellmatrix.numCols() != ncols)
  {
    cellmatrix.shape(nrows, ncols);
  }

  // simply zero out otherwise
  else
    cellmatrix.putScalar(0.0);
}


/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void ScaTra::MortarCellAssemblyStrategy::init_cell_vector(
    Core::LinAlg::SerialDenseVector& cellvector, const Inpar::S2I::InterfaceSides side,
    Core::Elements::Element::LocationArray& la_slave,
    Core::Elements::Element::LocationArray& la_master) const
{
  // determine number of vector components
  const int ndofs =
      side == Inpar::S2I::side_slave ? la_slave[nds_rows_].Size() : la_master[nds_rows_].Size();

  // reshape cell vector if necessary
  if (cellvector.length() != ndofs)
  {
    cellvector.size(ndofs);

    // simply zero out otherwise
  }
  else
    cellvector.putScalar(0.0);
}


// forward declarations
template class ScaTra::MortarCellCalc<Core::FE::CellType::tri3, Core::FE::CellType::tri3>;
template class ScaTra::MortarCellCalc<Core::FE::CellType::tri3, Core::FE::CellType::quad4>;
template class ScaTra::MortarCellCalc<Core::FE::CellType::quad4, Core::FE::CellType::tri3>;
template class ScaTra::MortarCellCalc<Core::FE::CellType::quad4, Core::FE::CellType::quad4>;

FOUR_C_NAMESPACE_CLOSE
