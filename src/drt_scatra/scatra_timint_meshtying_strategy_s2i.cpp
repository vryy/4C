/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for standard scalar transport problems

\level 2

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_s2i.H"
#include "scatra_timint_implicit.H"
#include "scatra_timint_meshtying_strategy_s2i_elch.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"

#include "../drt_fluid/fluid_utils.H"
#include "../linalg/linalg_matrixtransform.H"

#include "../drt_geometry/position_array.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"

#include "../drt_mortar/mortar_coupling3d_classes.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_mortar/mortar_projector.H"
#include "../drt_mortar/mortar_utils.H"

#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_boundary_calc.H"
#include "../drt_scatra_ele/scatra_ele_calc_utils.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../drt_scatra_ele/scatra_ele_parameter_boundary.H"

#include "../drt_volmortar/volmortar_shape.H"

#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I(
    SCATRA::ScaTraTimIntImpl* scatratimint, const Teuchos::ParameterList& parameters)
    : MeshtyingStrategyBase(scatratimint),
      interfacemaps_(Teuchos::null),
      blockmaps_slave_(Teuchos::null),
      blockmaps_master_(Teuchos::null),
      icoup_(Teuchos::null),
      icoupmortar_(),
      imortarcells_(),
      imortarredistribution_(
          DRT::INPUT::IntegralValue<INPAR::S2I::CouplingType>(parameters.sublist("S2I COUPLING"),
              "COUPLINGTYPE") == INPAR::S2I::coupling_mortar_standard and
          DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
              DRT::Problem::Instance()->MortarCouplingParams().sublist("PARALLEL REDISTRIBUTION"),
              "PARALLEL_REDIST") != INPAR::MORTAR::parredist_none),
      islavemap_(Teuchos::null),
      imastermap_(Teuchos::null),
      islavenodestomasterelements_(),
      islavenodesimpltypes_(),
      islavenodeslumpedareas_(),
      islavematrix_(Teuchos::null),
      imastermatrix_(Teuchos::null),
      couplingtype_(DRT::INPUT::IntegralValue<INPAR::S2I::CouplingType>(
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
      imasterphinp_(Teuchos::null),
      invrowsums_(Teuchos::null),
      invcolsums_(Teuchos::null),
      lmside_(DRT::INPUT::IntegralValue<INPAR::S2I::InterfaceSides>(
          parameters.sublist("S2I COUPLING"), "LMSIDE")),
      matrixtype_(Teuchos::getIntegralValue<INPAR::SCATRA::MatrixType>(parameters, "MATRIXTYPE")),
      ntsprojtol_(parameters.sublist("S2I COUPLING").get<double>("NTSPROJTOL")),
      intlayergrowth_evaluation_(DRT::INPUT::IntegralValue<INPAR::S2I::GrowthEvaluation>(
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
      slaveconditions_(),
      slaveonly_(DRT::INPUT::IntegralValue<bool>(parameters.sublist("S2I COUPLING"), "SLAVEONLY"))
{
  // empty constructor
}  // SCATRA::MeshtyingStrategyS2I::MeshtyingStrategyS2I


/*-----------------------------------------------------------------------*
 | condense global system of equations                        fang 01/16 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::CondenseMatAndRHS(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
    const bool calcinittimederiv  //!< flag for calculation of initial time derivative
    ) const
{
  switch (couplingtype_)
  {
    case INPAR::S2I::coupling_mortar_condensed_bubnov:
    {
      // extract global system matrix
      Teuchos::RCP<LINALG::SparseMatrix> sparsematrix = scatratimint_->SystemMatrix();
      if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

      if (lmside_ == INPAR::S2I::side_slave)
      {
        // initialize temporary matrix for slave-side rows of global system matrix
        LINALG::SparseMatrix sparsematrixrowsslave(*interfacemaps_->Map(1), 81);

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
          sparsematrix->Add(
              *LINALG::MLMultiply(*Q_, true, sparsematrixrowsslave, false, false, false, true),
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
            *LINALG::MLMultiply(*P_, true, sparsematrixrowsslave, false, false, false, true), false,
            1., 1.);

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
            dserror("Matrix-vector multiplication failed!");
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
          dserror("Matrix-vector multiplication failed!");
        interfacemaps_->AddVector(P_residualslave, 2, *scatratimint_->Residual());
      }

      else
      {
        // initialize temporary matrix for master-side rows of global system matrix
        LINALG::SparseMatrix sparsematrixrowsmaster(*interfacemaps_->Map(2), 81);

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
            *LINALG::MLMultiply(*Q_, true, sparsematrixrowsmaster, false, false, false, true),
            false, 1., 1.);
        // during calculation of initial time derivative, standard global system matrix is replaced
        // by global mass matrix, and hence interface contributions must not be included
        if (!calcinittimederiv) sparsematrix->Add(*imastermatrix_, false, 1., 1.);

        // add projected master-side rows to slave-side rows of global system matrix
        sparsematrix->Add(
            *LINALG::MLMultiply(*P_, true, sparsematrixrowsmaster, false, false, false, true),
            false, 1., 1.);

        // extract master-side entries of global residual vector
        Teuchos::RCP<Epetra_Vector> residualmaster =
            interfacemaps_->ExtractVector(scatratimint_->Residual(), 2);

        // replace master-side entries of global residual vector by projected master-side entries
        // including interface contributions
        Epetra_Vector Q_residualmaster(*interfacemaps_->Map(2));
        if (Q_->Multiply(true, *residualmaster, Q_residualmaster))
          dserror("Matrix-vector multiplication failed!");
        interfacemaps_->InsertVector(Q_residualmaster, 2, *scatratimint_->Residual());
        interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());

        // add projected master-side entries to slave-side entries of global residual vector
        Epetra_Vector P_residualmaster(*interfacemaps_->Map(1));
        if (P_->Multiply(true, *residualmaster, P_residualmaster))
          dserror("Matrix-vector multiplication failed!");
        interfacemaps_->AddVector(P_residualmaster, 1, *scatratimint_->Residual());
      }

      break;
    }

    case INPAR::S2I::coupling_matching_nodes:
    case INPAR::S2I::coupling_mortar_standard:
    case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
    case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
    case INPAR::S2I::coupling_mortar_condensed_petrov:
    case INPAR::S2I::coupling_nts_standard:
    {
      // do nothing in these cases
      break;
    }

    default:
    {
      dserror("Type of mortar meshtying for scatra-scatra interface coupling not recognized!");
      break;
    }
  }

  return;
}


/*-----------------------------------------------------------------------*
 | return global map of degrees of freedom                    fang 02/18 |
 *-----------------------------------------------------------------------*/
const Epetra_Map& SCATRA::MeshtyingStrategyS2I::DofRowMap() const
{
  return extendedmaps_ != Teuchos::null ? *extendedmaps_->FullMap() : *scatratimint_->DofRowMap();
}


/*-----------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions       fang 10/14 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying()
{
  // time measurement: evaluate condition 'S2ICoupling'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'S2ICoupling'");

  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> conditions;
  scatratimint_->Discretization()->GetCondition("S2ICoupling", conditions);

  switch (couplingtype_)
  {
    case INPAR::S2I::coupling_matching_nodes:
    {
      // create parameter list
      Teuchos::ParameterList condparams;

      // action for elements
      condparams.set<int>("action", SCATRA::bd_calc_s2icoupling);

      // number of dofset associated with displacement dofs
      if (scatratimint_->IsALE()) condparams.set<int>("ndsdisp", scatratimint_->NdsDisp());

      // set global state vectors according to time-integration scheme
      scatratimint_->Discretization()->ClearState();
      scatratimint_->AddTimeIntegrationSpecificVectors();

      // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
      islavematrix_->Zero();
      if (not slaveonly_) imastermatrix_->Zero();
      islaveresidual_->PutScalar(0.);
      for (unsigned icondition = 0; icondition < conditions.size(); ++icondition)
      {
        if (conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
        {
          // collect condition specific data and store to scatra boundary parameter class
          SetConditionSpecificScaTraParameters(*conditions[icondition]);

          if (not slaveonly_)
            scatratimint_->Discretization()->EvaluateCondition(condparams, islavematrix_,
                imastermatrix_, islaveresidual_, Teuchos::null, Teuchos::null, "S2ICoupling",
                conditions[icondition]->GetInt("ConditionID"));
          else
            scatratimint_->Discretization()->EvaluateCondition(condparams, islavematrix_,
                Teuchos::null, islaveresidual_, Teuchos::null, Teuchos::null, "S2ICoupling",
                conditions[icondition]->GetInt("ConditionID"));
        }
      }
      scatratimint_->Discretization()->ClearState();

      // finalize interface matrices
      islavematrix_->Complete();
      if (not slaveonly_) imastermatrix_->Complete();

      // assemble global system matrix depending on matrix type
      switch (matrixtype_)
      {
        case INPAR::SCATRA::MatrixType::sparse:
        {
          // check matrix
          Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = scatratimint_->SystemMatrix();
          if (systemmatrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

          // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
          systemmatrix->Add(*islavematrix_, false, 1., 1.);

          if (not slaveonly_)
          {
            // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomastercoltransform_)(imastermatrix_->RowMap(), imastermatrix_->ColMap(),
                *imastermatrix_, 1., ADAPTER::CouplingSlaveConverter(*icoup_), *systemmatrix, true,
                true);

            // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global
            // system matrix
            (*islavetomasterrowtransform_)(
                *islavematrix_, -1., ADAPTER::CouplingSlaveConverter(*icoup_), *systemmatrix, true);

            // derive linearizations of master fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomasterrowcoltransform_)(*imastermatrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), ADAPTER::CouplingSlaveConverter(*icoup_),
                *systemmatrix, true, true);
          }

          // In case the interface linearizations and residuals are evaluated on slave side only,
          // we now apply a standard meshtying algorithm to condense out the slave-side degrees of
          // freedom.
          else if (!scatratimint_->Discretization()->GetCondition("PointCoupling"))
          {
            // initialize temporary matrix for slave-side rows of system matrix
            LINALG::SparseMatrix systemmatrixrowsslave(*icoup_->SlaveDofMap(), 81);

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
              if (slavedofgid < 0) dserror("Couldn't find local ID %d in map!", slavedoflid);

              // determine global ID of associated master-side matrix column
              const int masterdofgid = icoup_->PermMasterDofMap()->GID(slavedoflid);
              if (masterdofgid < 0)
                dserror("Couldn't find local ID %d in permuted map!", slavedoflid);

              // insert value -1. into intersection of slave-side row and master-side column in
              // system matrix this effectively forces the slave-side degree of freedom to assume
              // the same value as the master-side degree of freedom
              const double value(-1.);
              if (systemmatrix->EpetraMatrix()->InsertGlobalValues(
                      slavedofgid, 1, &value, &masterdofgid) < 0)
                dserror(
                    "Cannot insert value -1. into matrix row with global ID %d and matrix column "
                    "with global ID %d!",
                    slavedofgid, masterdofgid);

              // insert zero into intersection of slave-side row and master-side column in temporary
              // matrix this prevents the system matrix from changing its graph when calling this
              // function again during the next Newton iteration
              const double zero(0.);
              if (systemmatrixrowsslave.EpetraMatrix()->InsertGlobalValues(
                      slavedofgid, 1, &zero, &masterdofgid) < 0)
                dserror(
                    "Cannot insert zero into matrix row with global ID %d and matrix column with "
                    "global ID %d!",
                    slavedofgid, masterdofgid);
            }

            // finalize temporary matrix with slave-side rows of system matrix
            systemmatrixrowsslave.Complete(*scatratimint_->DofRowMap(), *icoup_->SlaveDofMap());

            // add slave-side rows of system matrix to corresponding master-side rows to finalize
            // matrix condensation of slave-side degrees of freedom
            (*islavetomasterrowtransform_)(systemmatrixrowsslave, 1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), *systemmatrix, true);
          }

          break;
        }

        case INPAR::SCATRA::MatrixType::block_geometry:
        {
          // check matrix
          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix =
              scatratimint_->BlockSystemMatrix();
          if (blocksystemmatrix == Teuchos::null) dserror("System matrix is not a block matrix!");

          // assemble linearizations of slave fluxes w.r.t. slave dofs into global system matrix
          blocksystemmatrix->Matrix(1, 1).Add(*islavematrix_, false, 1., 1.);

          if (not slaveonly_)
          {
            // transform linearizations of slave fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomastercoltransform_)(imastermatrix_->RowMap(), imastermatrix_->ColMap(),
                *imastermatrix_, 1., ADAPTER::CouplingSlaveConverter(*icoup_),
                blocksystemmatrix->Matrix(1, 2));

            // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global
            // system matrix
            (*islavetomasterrowtransform_)(*islavematrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), blocksystemmatrix->Matrix(2, 1));

            // derive linearizations of master fluxes w.r.t. master dofs and assemble into global
            // system matrix
            (*islavetomasterrowcoltransform_)(*imastermatrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), ADAPTER::CouplingSlaveConverter(*icoup_),
                blocksystemmatrix->Matrix(2, 2), true, true);
          }

          // safety check
          else
            dserror(
                "Scatra-scatra interface coupling with evaluation of interface linearizations and "
                "residuals on slave side only is not yet available for block system matrices!");

          break;
        }

        case INPAR::SCATRA::MatrixType::block_condition:
        case INPAR::SCATRA::MatrixType::block_condition_dof:
        {
          // check matrix
          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix =
              scatratimint_->BlockSystemMatrix();
          if (blocksystemmatrix == Teuchos::null) dserror("System matrix is not a block matrix!");

          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkss(
              islavematrix_->Split<LINALG::DefaultBlockMatrixStrategy>(
                  *blockmaps_slave_, *blockmaps_slave_));
          blockkss->Complete();

          // assemble interface block matrix into global block system matrix
          blocksystemmatrix->Add(*blockkss, false, 1., 1.);

          if (not slaveonly_)
          {
            Teuchos::RCP<LINALG::SparseMatrix> ksm(
                Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->SlaveDofMap(), 81, false)));
            Teuchos::RCP<LINALG::SparseMatrix> kms(
                Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            Teuchos::RCP<LINALG::SparseMatrix> kmm(
                Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));

            // transform linearizations of slave fluxes w.r.t. master dofs
            (*islavetomastercoltransform_)(imastermatrix_->RowMap(), imastermatrix_->ColMap(),
                *imastermatrix_, 1., ADAPTER::CouplingSlaveConverter(*icoup_), *ksm);
            ksm->Complete(*icoup_->MasterDofMap(), *icoup_->SlaveDofMap());

            // derive linearizations of master fluxes w.r.t. slave dofs
            (*islavetomasterrowtransform_)(
                *islavematrix_, -1., ADAPTER::CouplingSlaveConverter(*icoup_), *kms);
            kms->Complete(*icoup_->SlaveDofMap(), *icoup_->MasterDofMap());

            // derive linearizations of master fluxes w.r.t. master dofs
            (*islavetomasterrowcoltransform_)(*imastermatrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), ADAPTER::CouplingSlaveConverter(*icoup_),
                *kmm);
            kmm->Complete();

            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockksm(
                ksm->Split<LINALG::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_slave_));
            blockksm->Complete();
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkms(
                kms->Split<LINALG::DefaultBlockMatrixStrategy>(
                    *blockmaps_slave_, *blockmaps_master_));
            blockkms->Complete();
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkmm(
                kmm->Split<LINALG::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_master_));
            blockkmm->Complete();

            // assemble interface block matrices into global block system matrix
            blocksystemmatrix->Add(*blockksm, false, 1., 1.);
            blocksystemmatrix->Add(*blockkms, false, 1., 1.);
            blocksystemmatrix->Add(*blockkmm, false, 1., 1.);
          }

          // safety check
          else
            dserror(
                "Scatra-scatra interface coupling with evaluation of interface linearizations and "
                "residuals on slave side only is not yet available for block system matrices!");

          break;
        }

        default:
        {
          dserror(
              "Type of global system matrix for scatra-scatra interface coupling not recognized!");
          break;
        }
      }

      // assemble slave residuals into global residual vector
      interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());

      if (not slaveonly_)
        // transform master residuals and assemble into global residual vector
        interfacemaps_->AddVector(
            icoup_->SlaveToMaster(islaveresidual_), 2, scatratimint_->Residual(), -1.);

      // In case the interface linearizations and residuals are evaluated on slave side only,
      // we now apply a standard meshtying algorithm to condense out the slave-side degrees of
      // freedom.
      else if (!scatratimint_->Discretization()->GetCondition("PointCoupling"))
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
          if (slavedofgid < 0) dserror("Couldn't find local ID %d in map!", slavedoflid);

          // copy current vector entry into temporary vector
          if (residualslave->ReplaceGlobalValue(slavedofgid, 0,
                  (*scatratimint_->Residual())[scatratimint_->DofRowMap()->LID(slavedofgid)]))
            dserror("Cannot insert residual vector entry with global ID %d into temporary vector!",
                slavedofgid);

          // zero out current vector entry
          if (scatratimint_->Residual()->ReplaceGlobalValue(slavedofgid, 0, 0.))
            dserror(
                "Cannot insert zero into residual vector entry with global ID %d!", slavedofgid);
        }

        // add slave-side entries of residual vector to corresponding master-side entries to
        // finalize vector condensation of slave-side degrees of freedom
        interfacemaps_->AddVector(
            icoup_->SlaveToMaster(residualslave), 2, scatratimint_->Residual());
      }

      break;
    }

    case INPAR::S2I::coupling_mortar_standard:
    case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
    case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
    case INPAR::S2I::coupling_mortar_condensed_petrov:
    case INPAR::S2I::coupling_mortar_condensed_bubnov:
    case INPAR::S2I::coupling_nts_standard:
    {
      // initialize auxiliary system matrix and vector for slave side
      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_slave or couplingtype_ == INPAR::S2I::coupling_nts_standard)
      {
        islavematrix_->Zero();
        islaveresidual_->PutScalar(0.);
      }

      // initialize auxiliary system matrix and vector for master side
      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_master or couplingtype_ == INPAR::S2I::coupling_nts_standard)
      {
        imastermatrix_->Zero();
        imasterresidual_->PutScalar(0.);
      }

      // loop over all scatra-scatra coupling interfaces
      for (std::map<const int, DRT::Condition* const>::iterator islavecondition =
               slaveconditions_.begin();
           islavecondition != slaveconditions_.end(); ++islavecondition)
      {
        // extract mortar interface discretization
        DRT::Discretization& idiscret =
            icoupmortar_[islavecondition->first]->Interface()->Discret();

        // export global state vector to mortar interface
        Teuchos::RCP<Epetra_Vector> iphinp =
            Teuchos::rcp(new Epetra_Vector(*idiscret.DofColMap(), false));
        LINALG::Export(*scatratimint_->Phiafnp(), *iphinp);
        idiscret.SetState("iphinp", iphinp);

        // create parameter list for mortar integration cells
        Teuchos::ParameterList params;

        // add current condition to parameter list
        params.set<DRT::Condition*>("condition", islavecondition->second);

        // collect condition specific data and store to scatra boundary parameter class
        SetConditionSpecificScaTraParameters(*(islavecondition->second));

        if (couplingtype_ != INPAR::S2I::coupling_nts_standard)
        {
          // set action
          params.set<int>("action", INPAR::S2I::evaluate_condition);

          // evaluate mortar integration cells at current interface
          EvaluateMortarCells(idiscret, params, islavematrix_, INPAR::S2I::side_slave,
              INPAR::S2I::side_slave, islavematrix_, INPAR::S2I::side_slave,
              INPAR::S2I::side_master, imastermatrix_, INPAR::S2I::side_master,
              INPAR::S2I::side_slave, imastermatrix_, INPAR::S2I::side_master,
              INPAR::S2I::side_master, islaveresidual_, INPAR::S2I::side_slave, imasterresidual_,
              INPAR::S2I::side_master);
        }

        else
        {
          // set action
          params.set<int>("action", INPAR::S2I::evaluate_condition_nts);

          // evaluate note-to-segment coupling at current interface
          EvaluateNTS(*islavenodestomasterelements_[islavecondition->first],
              *islavenodeslumpedareas_[islavecondition->first],
              *islavenodesimpltypes_[islavecondition->first], idiscret, params, islavematrix_,
              INPAR::S2I::side_slave, INPAR::S2I::side_slave, islavematrix_, INPAR::S2I::side_slave,
              INPAR::S2I::side_master, imastermatrix_, INPAR::S2I::side_master,
              INPAR::S2I::side_slave, imastermatrix_, INPAR::S2I::side_master,
              INPAR::S2I::side_master, islaveresidual_, INPAR::S2I::side_slave, imasterresidual_,
              INPAR::S2I::side_master);
        }
      }

      // finalize auxiliary system matrix for slave side
      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_slave or couplingtype_ == INPAR::S2I::coupling_nts_standard)
        islavematrix_->Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(1));

      // finalize auxiliary system matrix and residual vector for master side
      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_master or couplingtype_ == INPAR::S2I::coupling_nts_standard)
      {
        imastermatrix_->Complete(*interfacemaps_->FullMap(), *interfacemaps_->Map(2));
        if (imasterresidual_->GlobalAssemble(Add, true))
          dserror("Assembly of auxiliary residual vector for master residuals not successful!");
      }

      // assemble global system of equations depending on matrix type
      switch (matrixtype_)
      {
        case INPAR::SCATRA::MatrixType::sparse:
        {
          // extract global system matrix from time integrator
          const Teuchos::RCP<LINALG::SparseMatrix> systemmatrix = scatratimint_->SystemMatrix();
          if (systemmatrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

          // assemble interface contributions into global system of equations
          switch (couplingtype_)
          {
            case INPAR::S2I::coupling_mortar_standard:
            case INPAR::S2I::coupling_nts_standard:
            {
              const Teuchos::RCP<const LINALG::SparseMatrix> islavematrix =
                  not imortarredistribution_
                      ? islavematrix_
                      : MORTAR::MatrixRowTransform(islavematrix_, islavemap_);
              const Teuchos::RCP<const LINALG::SparseMatrix> imastermatrix =
                  not imortarredistribution_
                      ? imastermatrix_
                      : MORTAR::MatrixRowTransform(imastermatrix_, imastermap_);
              systemmatrix->Add(*islavematrix, false, 1., 1.);
              systemmatrix->Add(*imastermatrix, false, 1., 1.);
              interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
              interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());

              break;
            }

            case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
            case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
            {
              if (lmside_ == INPAR::S2I::side_slave)
              {
                // assemble slave-side interface contributions into global residual vector
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (D_->Multiply(true, *lm_, islaveresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual(), -1.);

                // assemble master-side interface contributions into global residual vector
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (M_->Multiply(true, *lm_, imasterresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual());

                // build constraint residual vector associated with Lagrange multiplier dofs
                Epetra_Vector ilmresidual(*islaveresidual_);
                if (ilmresidual.ReplaceMap(*extendedmaps_->Map(1)))
                  dserror("Couldn't replace map!");
                if (lmresidual_->Update(1., ilmresidual, 0.)) dserror("Vector update failed!");
                if (E_->Multiply(true, *lm_, ilmresidual))
                  dserror("Matrix-vector multiplication failed!");
                if (lmresidual_->Update(1., ilmresidual, 1.)) dserror("Vector update failed!");
              }
              else
              {
                // assemble slave-side interface contributions into global residual vector
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (M_->Multiply(true, *lm_, islaveresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual());

                // assemble master-side interface contributions into global residual vector
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (D_->Multiply(true, *lm_, imasterresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual(), -1.);

                // build constraint residual vector associated with Lagrange multiplier dofs
                Epetra_Vector ilmresidual(Copy, *imasterresidual_, 0);
                if (ilmresidual.ReplaceMap(*extendedmaps_->Map(1)))
                  dserror("Couldn't replace map!");
                if (lmresidual_->Update(1., ilmresidual, 0.)) dserror("Vector update failed!");
                if (E_->Multiply(true, *lm_, ilmresidual))
                  dserror("Matrix-vector multiplication failed!");
                if (lmresidual_->Update(1., ilmresidual, 1.)) dserror("Vector update failed!");
              }

              break;
            }

            case INPAR::S2I::coupling_mortar_condensed_petrov:
            {
              if (lmside_ == INPAR::S2I::side_slave)
              {
                systemmatrix->Add(*islavematrix_, false, 1., 1.);
                systemmatrix->Add(
                    *LINALG::MLMultiply(*P_, true, *islavematrix_, false, false, false, true),
                    false, -1., 1.);
                interfacemaps_->AddVector(islaveresidual_, 1, scatratimint_->Residual());
                Epetra_Vector imasterresidual(*interfacemaps_->Map(2));
                if (P_->Multiply(true, *islaveresidual_, imasterresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(imasterresidual, 2, *scatratimint_->Residual(), -1.);
              }
              else
              {
                systemmatrix->Add(
                    *LINALG::MLMultiply(*P_, true, *imastermatrix_, false, false, false, true),
                    false, -1., 1.);
                systemmatrix->Add(*imastermatrix_, false, 1., 1.);
                Epetra_Vector islaveresidual(*interfacemaps_->Map(1));
                if (P_->Multiply(true, *imasterresidual_, islaveresidual))
                  dserror("Matrix-vector multiplication failed!");
                interfacemaps_->AddVector(islaveresidual, 1, *scatratimint_->Residual(), -1.);
                interfacemaps_->AddVector(*imasterresidual_, 2, *scatratimint_->Residual());
              }

              break;
            }

            case INPAR::S2I::coupling_mortar_condensed_bubnov:
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
              dserror("Not yet implemented!");
              break;
            }
          }

          break;
        }

        case INPAR::SCATRA::MatrixType::block_condition:
        {
          // extract global system matrix from time integrator
          Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix =
              scatratimint_->BlockSystemMatrix();
          if (blocksystemmatrix == Teuchos::null) dserror("System matrix is not a block matrix!");

          // assemble interface contributions into global system of equations
          switch (couplingtype_)
          {
            case INPAR::S2I::coupling_mortar_standard:
            {
              // split interface sparse matrices into block matrices
              const Teuchos::RCP<const LINALG::SparseMatrix> islavematrix =
                  not imortarredistribution_
                      ? islavematrix_
                      : MORTAR::MatrixRowTransform(islavematrix_, islavemap_);
              Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockslavematrix(
                  islavematrix->Split<LINALG::DefaultBlockMatrixStrategy>(
                      scatratimint_->BlockMaps(), *blockmaps_slave_));
              blockslavematrix->Complete();
              const Teuchos::RCP<const LINALG::SparseMatrix> imastermatrix =
                  not imortarredistribution_
                      ? imastermatrix_
                      : MORTAR::MatrixRowTransform(imastermatrix_, imastermap_);
              Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmastermatrix(
                  imastermatrix->Split<LINALG::DefaultBlockMatrixStrategy>(
                      scatratimint_->BlockMaps(), *blockmaps_master_));
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
              dserror("Not yet implemented!");
              break;
            }
          }

          break;
        }

        default:
        {
          dserror("Not yet implemented!");
          break;
        }
      }

      break;
    }

    default:
    {
      dserror("Not yet implemented!");
      break;
    }
  }

  // extract boundary conditions for scatra-scatra interface layer growth
  scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth", conditions);

  // evaluate scatra-scatra interface layer growth
  if (conditions.size())
  {
    switch (couplingtype_)
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // create parameter list for elements
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action", SCATRA::bd_calc_s2icoupling);

        // set global state vectors according to time-integration scheme
        scatratimint_->Discretization()->ClearState();
        scatratimint_->AddTimeIntegrationSpecificVectors();

        // evaluate scatra-scatra interface coupling at time t_{n+1} or t_{n+alpha_F}
        islavematrix_->Zero();
        imastermatrix_->Zero();
        islaveresidual_->PutScalar(0.);

        // collect condition specific data and store to scatra boundary parameter class
        SetConditionSpecificScaTraParameters(*conditions[0]);
        // evaluate the condition
        scatratimint_->Discretization()->EvaluateCondition(condparams, islavematrix_,
            imastermatrix_, islaveresidual_, Teuchos::null, Teuchos::null, "S2ICouplingGrowth");
        scatratimint_->Discretization()->ClearState();

        // finalize interface matrices
        islavematrix_->Complete();
        imastermatrix_->Complete();

        // assemble interface matrices into global system matrix depending on matrix type
        switch (matrixtype_)
        {
          case INPAR::SCATRA::MatrixType::sparse:
          {
            // check matrix
            const Teuchos::RCP<LINALG::SparseMatrix>& systemmatrix = scatratimint_->SystemMatrix();
            if (systemmatrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

            // We assume that the scatra-scatra interface layer growth is caused by master-side
            // fluxes to the interface, whereas there is no mass exchange between the interface
            // layer and the slave side of the interface. Hence, we only need to linearize the
            // master-side fluxes w.r.t. the slave-side and master-side degrees of freedom.

            // derive linearizations of master fluxes w.r.t. slave dofs and assemble into global
            // system matrix
            LINALG::MatrixRowTransform()(
                *islavematrix_, -1., ADAPTER::CouplingSlaveConverter(*icoup_), *systemmatrix, true);

            // derive linearizations of master fluxes w.r.t. master dofs and assemble into global
            // system matrix
            LINALG::MatrixRowColTransform()(*imastermatrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), ADAPTER::CouplingSlaveConverter(*icoup_),
                *systemmatrix, true, true);

            break;
          }

          case INPAR::SCATRA::MatrixType::block_condition:
          case INPAR::SCATRA::MatrixType::block_condition_dof:
          {
            // check matrix
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksystemmatrix =
                scatratimint_->BlockSystemMatrix();
            if (blocksystemmatrix == Teuchos::null) dserror("System matrix is not a block matrix!");

            // We assume that the scatra-scatra interface layer growth is caused by master-side
            // fluxes to the interface, whereas there is no mass exchange between the interface
            // layer and the slave side of the interface. Hence, we only need to linearize the
            // master-side fluxes w.r.t. the slave-side and master-side degrees of freedom.

            // derive linearizations of master fluxes w.r.t. slave dofs
            Teuchos::RCP<LINALG::SparseMatrix> kms(
                Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            LINALG::MatrixRowTransform()(
                *islavematrix_, -1., ADAPTER::CouplingSlaveConverter(*icoup_), *kms);
            kms->Complete(*icoup_->SlaveDofMap(), *icoup_->MasterDofMap());
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkms(
                kms->Split<LINALG::DefaultBlockMatrixStrategy>(
                    *blockmaps_slave_, *blockmaps_master_));
            blockkms->Complete();

            // derive linearizations of master fluxes w.r.t. master dofs
            Teuchos::RCP<LINALG::SparseMatrix> kmm(
                Teuchos::rcp(new LINALG::SparseMatrix(*icoup_->MasterDofMap(), 81, false)));
            LINALG::MatrixRowColTransform()(*imastermatrix_, -1.,
                ADAPTER::CouplingSlaveConverter(*icoup_), ADAPTER::CouplingSlaveConverter(*icoup_),
                *kmm);
            kmm->Complete();
            Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkmm(
                kmm->Split<LINALG::DefaultBlockMatrixStrategy>(
                    *blockmaps_master_, *blockmaps_master_));
            blockkmm->Complete();

            // assemble interface block matrices into global block system matrix
            blocksystemmatrix->Add(*blockkms, false, 1., 1.);
            blocksystemmatrix->Add(*blockkmm, false, 1., 1.);

            break;
          }

          default:
          {
            dserror(
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
        if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
        {
          // extract map associated with scalar transport degrees of freedom
          const Epetra_Map& dofrowmap_scatra = *scatratimint_->Discretization()->DofRowMap();

          // extract map associated with scatra-scatra interface layer thicknesses
          const Epetra_Map& dofrowmap_growth = *scatratimint_->Discretization()->DofRowMap(2);

          // extract ID of boundary condition for scatra-scatra interface layer growth
          // the corresponding boundary condition for scatra-scatra interface coupling is expected
          // to have the same ID
          const int condid = scatratimint_->Discretization()
                                 ->GetCondition("S2ICouplingGrowth")
                                 ->GetInt("ConditionID");

          // set global state vectors according to time-integration scheme
          scatratimint_->Discretization()->ClearState();
          scatratimint_->AddTimeIntegrationSpecificVectors();

          // compute additional linearizations and residuals depending on type of scalar transport
          // system matrix
          switch (matrixtype_)
          {
            case INPAR::SCATRA::MatrixType::sparse:
            {
              // assemble off-diagonal scatra-growth block of global system matrix, containing
              // derivatives of discrete scatra residuals w.r.t. discrete scatra-scatra interface
              // layer thicknesses
              {
                // check matrix
                const Teuchos::RCP<LINALG::SparseMatrix> scatragrowthblock =
                    Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(scatragrowthblock_);
                if (scatragrowthblock == Teuchos::null) dserror("Matrix is not a sparse matrix!");

                // initialize matrix block
                scatragrowthblock->Zero();

                // initialize auxiliary matrix block for linearizations of slave fluxes w.r.t.
                // scatra-scatra interface layer thicknesses
                Teuchos::RCP<LINALG::SparseMatrix> islavematrix =
                    Teuchos::rcp(new LINALG::SparseMatrix(*(icoup_)->SlaveDofMap(), 81));

                // initialize assembly strategy for auxiliary matrix block
                DRT::AssembleStrategy strategy(
                    0, 2, islavematrix, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_scatragrowth);

                // evaluate off-diagonal linearizations arising from scatra-scatra interface
                // coupling
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICoupling", condid);

                // finalize auxiliary matrix block
                islavematrix->Complete(dofrowmap_growth, dofrowmap_scatra);

                // assemble linearizations of slave fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses into global matrix
                // block
                scatragrowthblock->Add(*islavematrix, false, 1., 0.);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses and assemble into
                // global matrix block
                LINALG::MatrixRowTransform()(*islavematrix, -1.,
                    ADAPTER::CouplingSlaveConverter(*icoup_), *scatragrowthblock, false);

                // zero out auxiliary matrix block for subsequent evaluation
                islavematrix->Zero();

                // evaluate off-diagonal linearizations arising from scatra-scatra interface layer
                // growth
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICouplingGrowth", condid);

                // finalize auxiliary matrix block
                islavematrix->Complete(dofrowmap_growth, dofrowmap_scatra);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // layer growth w.r.t. scatra-scatra interface layer thicknesses and assemble into
                // global matrix block
                LINALG::MatrixRowTransform()(*islavematrix, -1.,
                    ADAPTER::CouplingSlaveConverter(*icoup_), *scatragrowthblock, true);

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
                const Teuchos::RCP<LINALG::SparseMatrix> growthscatrablock =
                    Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(growthscatrablock_);
                if (growthscatrablock == Teuchos::null) dserror("Matrix is not a sparse matrix!");

                // initialize matrix block
                growthscatrablock->Zero();

                // initialize auxiliary matrix blocks for linearizations of scatra-scatra interface
                // layer growth residuals w.r.t. slave-side and master-side scalar transport degrees
                // of freedom
                Teuchos::RCP<LINALG::SparseMatrix> islavematrix =
                    Teuchos::rcp(new LINALG::SparseMatrix(dofrowmap_growth, 81));
                Teuchos::RCP<LINALG::SparseMatrix> imastermatrix =
                    Teuchos::rcp(new LINALG::SparseMatrix(dofrowmap_growth, 81));

                // initialize assembly strategy for auxiliary matrix block
                DRT::AssembleStrategy strategy(
                    2, 0, islavematrix, imastermatrix, Teuchos::null, Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_growthscatra);

                // evaluate off-diagonal linearizations
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICouplingGrowth", condid);

                // finalize auxiliary matrix blocks
                islavematrix->Complete(dofrowmap_scatra, dofrowmap_growth);
                imastermatrix->Complete(dofrowmap_scatra, dofrowmap_growth);

                // assemble linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // slave-side scalar transport degrees of freedom into global matrix block
                growthscatrablock->Add(*islavematrix, false, 1., 0.);

                // derive linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // master-side scalar transport degrees of freedom and assemble into global matrix
                // block
                LINALG::MatrixColTransform()(imastermatrix->RowMap(), imastermatrix->ColMap(),
                    *imastermatrix, 1., ADAPTER::CouplingSlaveConverter(*icoup_),
                    *growthscatrablock, true, true);

                // finalize global matrix block
                growthscatrablock->Complete(dofrowmap_scatra, dofrowmap_growth);
              }

              break;
            }

            case INPAR::SCATRA::MatrixType::block_condition:
            case INPAR::SCATRA::MatrixType::block_condition_dof:
            {
              // assemble off-diagonal scatra-growth block of global system matrix, containing
              // derivatives of discrete scatra residuals w.r.t. discrete scatra-scatra interface
              // layer thicknesses
              {
                // initialize auxiliary matrix block for linearizations of slave fluxes w.r.t.
                // scatra-scatra interface layer thicknesses
                const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockslavematrix =
                    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                        *blockmapgrowth_, *blockmaps_slave_, 81, false, true));

                // initialize assembly strategy for auxiliary matrix block
                DRT::AssembleStrategy strategy(0, 2, blockslavematrix, Teuchos::null, Teuchos::null,
                    Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_scatragrowth);

                // evaluate off-diagonal linearizations arising from scatra-scatra interface
                // coupling
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICoupling", condid);

                // finalize auxiliary matrix block
                blockslavematrix->Complete();

                // assemble linearizations of slave fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses into global matrix
                // block
                scatragrowthblock_->Add(*blockslavematrix, false, 1., 0.);

                // initialize auxiliary system matrix for linearizations of master fluxes associated
                // with scatra-scatra interface coupling w.r.t. scatra-scatra interface layer
                // thicknesses
                LINALG::SparseMatrix mastermatrix(*icoup_->MasterDofMap(), 27, false, true);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // coupling w.r.t. scatra-scatra interface layer thicknesses
                for (int iblock = 0; iblock < blockmaps_slave_->NumMaps(); ++iblock)
                  LINALG::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.,
                      ADAPTER::CouplingSlaveConverter(*icoup_), mastermatrix, true);

                // zero out auxiliary matrices for subsequent evaluation
                blockslavematrix->Zero();
                mastermatrix.Zero();

                // evaluate off-diagonal linearizations arising from scatra-scatra interface layer
                // growth
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICouplingGrowth", condid);

                // derive linearizations of master fluxes associated with scatra-scatra interface
                // layer growth w.r.t. scatra-scatra interface layer thicknesses
                for (int iblock = 0; iblock < blockmaps_slave_->NumMaps(); ++iblock)
                  LINALG::MatrixRowTransform()(blockslavematrix->Matrix(iblock, 0), -1.,
                      ADAPTER::CouplingSlaveConverter(*icoup_), mastermatrix, true);

                // finalize auxiliary system matrix
                mastermatrix.Complete(dofrowmap_growth, *icoup_->MasterDofMap());

                // split auxiliary system matrix and assemble into global matrix block
                const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmastermatrix =
                    mastermatrix.Split<LINALG::DefaultBlockMatrixStrategy>(
                        *blockmapgrowth_, scatratimint_->BlockMaps());
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
                const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockslavematrix =
                    Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                        *blockmaps_slave_, *blockmapgrowth_, 81, false, true));
                const Teuchos::RCP<LINALG::SparseMatrix> imastermatrix =
                    Teuchos::rcp(new LINALG::SparseMatrix(dofrowmap_growth, 81));

                // initialize assembly strategy for auxiliary matrix blocks
                DRT::AssembleStrategy strategy(2, 0, blockslavematrix, imastermatrix, Teuchos::null,
                    Teuchos::null, Teuchos::null);

                // create parameter list for elements
                Teuchos::ParameterList condparams;

                // set action for elements
                condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_growthscatra);

                // evaluate off-diagonal linearizations
                scatratimint_->Discretization()->EvaluateCondition(
                    condparams, strategy, "S2ICouplingGrowth", condid);

                // finalize auxiliary matrix blocks
                blockslavematrix->Complete();
                imastermatrix->Complete(dofrowmap_scatra, dofrowmap_growth);

                // assemble linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // slave-side scalar transport degrees of freedom into global matrix block
                growthscatrablock_->Add(*blockslavematrix, false, 1., 0.);

                // initialize temporary matrix
                LINALG::SparseMatrix kgm(dofrowmap_growth, 27, false, true);

                // derive linearizations of scatra-scatra interface layer growth residuals w.r.t.
                // master-side scalar transport degrees of freedom
                LINALG::MatrixColTransform()(imastermatrix->RowMap(), imastermatrix->ColMap(),
                    *imastermatrix, 1., ADAPTER::CouplingSlaveConverter(*icoup_), kgm);

                // finalize temporary matrix
                kgm.Complete(*icoup_->MasterDofMap(), dofrowmap_growth);

                // split temporary matrix and assemble into global matrix block
                const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockkgm(
                    kgm.Split<LINALG::DefaultBlockMatrixStrategy>(
                        scatratimint_->BlockMaps(), *blockmapgrowth_));
                blockkgm->Complete();
                growthscatrablock_->Add(*blockkgm, false, 1., 1.);

                // finalize global matrix block
                growthscatrablock_->Complete();
              }

              break;
            }

            default:
            {
              dserror(
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
            DRT::AssembleStrategy strategy(2, 2, growthgrowthblock_, Teuchos::null, growthresidual_,
                Teuchos::null, Teuchos::null);

            // create parameter list for elements
            Teuchos::ParameterList condparams;

            // set action for elements
            condparams.set<int>("action", SCATRA::bd_calc_s2icoupling_growthgrowth);

            // set history vector associated with discrete scatra-scatra interface layer thicknesses
            scatratimint_->Discretization()->SetState(2, "growthhist", growthhist_);

            // evaluate main-diagonal linearizations and corresponding residuals
            scatratimint_->Discretization()->EvaluateCondition(
                condparams, strategy, "S2ICouplingGrowth", condid);

            // finalize global matrix block
            growthgrowthblock_->Complete();
          }

          // clear discretization
          scatratimint_->Discretization()->ClearState();
        }  // monolithic evaluation of scatra-scatra interface layer growth

        break;
      }

      default:
      {
        dserror(
            "Evaluation of scatra-scatra interface layer growth only implemented for conforming "
            "interface discretizations!");
        break;
      }
    }
  }

  return;
}  // SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying


/*--------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell                                   fang 01/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarCell(
    const DRT::Discretization& idiscret,  //!< interface discretization
    MORTAR::IntCell& cell,                //!< mortar integration cell
    const INPAR::SCATRA::ImplType&
        impltype,  //!< physical implementation type of mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
    Epetra_SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
    Epetra_SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
    Epetra_SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
    Epetra_SerialDenseVector& cellvector1,   //!< cell vector 1
    Epetra_SerialDenseVector& cellvector2    //!< cell vector 2
    ) const
{
  // evaluate single mortar integration cell
  SCATRA::MortarCellFactory::MortarCellCalc(
      impltype, slaveelement, masterelement, couplingtype_, lmside_, idiscret.Name())
      ->Evaluate(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate single slave-side node for node-to-segment coupling              fang 08/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateSlaveNode(
    const DRT::Discretization& idiscret,  //!< interface discretization
    const MORTAR::MortarNode& slavenode,  //!< slave-side node
    const double& lumpedarea,  //!< lumped interface area fraction associated with slave-side node
    const INPAR::SCATRA::ImplType&
        impltype,  //!< physical implementation type of mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& ntsmatrix1,    //!< node-to-segment matrix 1
    Epetra_SerialDenseMatrix& ntsmatrix2,    //!< node-to-segment matrix 2
    Epetra_SerialDenseMatrix& ntsmatrix3,    //!< node-to-segment matrix 3
    Epetra_SerialDenseMatrix& ntsmatrix4,    //!< node-to-segment matrix 4
    Epetra_SerialDenseVector& ntsvector1,    //!< node-to-segment vector 1
    Epetra_SerialDenseVector& ntsvector2     //!< node-to-segment vector 2
    ) const
{
  // evaluate single slave-side node
  SCATRA::MortarCellFactory::MortarCellCalc(
      impltype, slaveelement, masterelement, couplingtype_, lmside_, idiscret.Name())
      ->EvaluateNTS(idiscret, slavenode, lumpedarea, slaveelement, masterelement, la_slave,
          la_master, params, ntsmatrix1, ntsmatrix2, ntsmatrix3, ntsmatrix4, ntsvector1,
          ntsvector2);

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate single mortar element                                            fang 08/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarElement(
    const DRT::Discretization& idiscret,      //!< interface discretization
    MORTAR::MortarElement& element,           //!< mortar element
    const INPAR::SCATRA::ImplType& impltype,  //!< physical implementation type of mortar element
    DRT::Element::LocationArray& la,          //!< location array
    const Teuchos::ParameterList& params,     //!< parameter list
    Epetra_SerialDenseMatrix& elematrix1,     //!< element matrix 1
    Epetra_SerialDenseMatrix& elematrix2,     //!< element matrix 2
    Epetra_SerialDenseMatrix& elematrix3,     //!< element matrix 3
    Epetra_SerialDenseMatrix& elematrix4,     //!< element matrix 4
    Epetra_SerialDenseVector& elevector1,     //!< element vector 1
    Epetra_SerialDenseVector& elevector2      //!< element vector 2
    ) const
{
  // evaluate single mortar element
  SCATRA::MortarCellFactory::MortarCellCalc(
      impltype, element, element, couplingtype_, lmside_, idiscret.Name())
      ->EvaluateMortarElement(idiscret, element, la, params, elematrix1, elematrix2, elematrix3,
          elematrix4, elevector1, elevector2);

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate mortar integration cells                                         fang 05/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarCells(
    const DRT::Discretization& idiscret,  //!< interface discretization
    const Teuchos::ParameterList&
        params,  //!< parameter list for evaluation of mortar integration cells
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix1,  //!< system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_rows,  //!< interface side associated with rows of system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_cols,  //!< interface side associated with columns of system matrix 1
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix2,  //!< system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_rows,  //!< interface side associated with rows of system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_cols,  //!< interface side associated with columns of system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix3,  //!< system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_rows,  //!< interface side associated with rows of system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_cols,  //!< interface side associated with columns of system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix4,  //!< system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_rows,  //!< interface side associated with rows of system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_cols,  //!< interface side associated with columns of system matrix 4
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,  //!< system vector 1
    const INPAR::S2I::InterfaceSides
        vector1_side,  //!< interface side associated with system vector 1
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,  //!< system vector 2
    const INPAR::S2I::InterfaceSides
        vector2_side  //!< interface side associated with system vector 2
    ) const
{
  // instantiate assembly strategy for mortar integration cells
  SCATRA::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // evaluate mortar integration cells
  EvaluateMortarCells(idiscret, params, strategy);

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate mortar integration cells                                         fang 01/17 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarCells(
    const DRT::Discretization& idiscret,  //!< interface discretization
    const Teuchos::ParameterList&
        params,  //!< parameter list for evaluation of mortar integration cells
    SCATRA::MortarCellAssemblyStrategy&
        strategy  //!< assembly strategy for mortar integration cells
    ) const
{
  // extract scatra-scatra interface coupling condition from parameter list
  const DRT::Condition* const condition = params.get<DRT::Condition*>("condition");
  if (condition == nullptr) dserror("Cannot access scatra-scatra interface coupling condition!");

  // extract mortar integration cells associated with current condition
  const std::vector<std::pair<Teuchos::RCP<MORTAR::IntCell>, INPAR::SCATRA::ImplType>>& cells =
      imortarcells_.at(condition->GetInt("ConditionID"));

  // loop over all mortar integration cells
  for (unsigned icell = 0; icell < cells.size(); ++icell)
  {
    // extract current mortar integration cell
    const Teuchos::RCP<MORTAR::IntCell>& cell = cells[icell].first;
    if (cell == Teuchos::null) dserror("Invalid mortar integration cell!");

    // extract slave-side element associated with current cell
    MORTAR::MortarElement* slaveelement =
        dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(cell->GetSlaveId()));
    if (!slaveelement)
      dserror("Couldn't extract slave element from mortar interface discretization!");

    // extract master-side element associated with current cell
    MORTAR::MortarElement* masterelement =
        dynamic_cast<MORTAR::MortarElement*>(idiscret.gElement(cell->GetMasterId()));
    if (!masterelement)
      dserror("Couldn't extract master element from mortar interface discretization!");

    // safety check
    if (!slaveelement->IsSlave() or masterelement->IsSlave())
      dserror("Something is wrong with the slave-master element pairing!");

    // construct slave-side and master-side location arrays
    DRT::Element::LocationArray la_slave(idiscret.NumDofSets());
    slaveelement->LocationVector(idiscret, la_slave, false);
    DRT::Element::LocationArray la_master(idiscret.NumDofSets());
    masterelement->LocationVector(idiscret, la_master, false);

    // initialize cell matrices and vectors
    strategy.InitCellMatricesAndVectors(la_slave, la_master);

    // evaluate current cell
    EvaluateMortarCell(idiscret, *cell, cells[icell].second, *slaveelement, *masterelement,
        la_slave, la_master, params, strategy.CellMatrix1(), strategy.CellMatrix2(),
        strategy.CellMatrix3(), strategy.CellMatrix4(), strategy.CellVector1(),
        strategy.CellVector2());

    // assemble cell matrices and vectors into system matrices and vectors
    strategy.AssembleCellMatricesAndVectors(la_slave, la_master, la_slave[0].lmowner_[0]);
  }

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate node-to-segment coupling                                         fang 08/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateNTS(
    const Epetra_IntVector&
        islavenodestomasterelements,              //!< vector for node-to-segment connectivity
    const Epetra_Vector& islavenodeslumpedareas,  //!< vector for lumped interface area fractions
                                                  //!< associated with slave-side nodes
    const Epetra_IntVector&
        islavenodesimpltypes,  //!< vector for physical implementation types of slave-side nodes
    const DRT::Discretization& idiscret,  //!< interface discretization
    const Teuchos::ParameterList&
        params,  //!< parameter list for evaluation of mortar integration cells
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix1,  //!< system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_rows,  //!< interface side associated with rows of system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_cols,  //!< interface side associated with columns of system matrix 1
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix2,  //!< system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_rows,  //!< interface side associated with rows of system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_cols,  //!< interface side associated with columns of system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix3,  //!< system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_rows,  //!< interface side associated with rows of system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_cols,  //!< interface side associated with columns of system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix4,  //!< system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_rows,  //!< interface side associated with rows of system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_cols,  //!< interface side associated with columns of system matrix 4
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,  //!< system vector 1
    const INPAR::S2I::InterfaceSides
        vector1_side,  //!< interface side associated with system vector 1
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,  //!< system vector 2
    const INPAR::S2I::InterfaceSides
        vector2_side  //!< interface side associated with system vector 2
    ) const
{
  // instantiate assembly strategy for node-to-segment coupling
  SCATRA::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // extract slave-side noderowmap
  const Epetra_BlockMap& noderowmap_slave = islavenodestomasterelements.Map();

  // loop over all slave-side nodes
  for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
  {
    // extract slave-side node
    MORTAR::MortarNode* const slavenode =
        dynamic_cast<MORTAR::MortarNode* const>(idiscret.gNode(noderowmap_slave.GID(inode)));
    if (slavenode == nullptr) dserror("Couldn't extract slave-side node from discretization!");

    // extract first slave-side element associated with current slave-side node
    MORTAR::MortarElement* const slaveelement =
        dynamic_cast<MORTAR::MortarElement* const>(slavenode->Elements()[0]);
    if (!slaveelement) dserror("Invalid slave-side mortar element!");

    // extract master-side element associated with current slave-side node
    MORTAR::MortarElement* const masterelement = dynamic_cast<MORTAR::MortarElement* const>(
        idiscret.gElement(islavenodestomasterelements[inode]));
    if (!masterelement) dserror("Invalid master-side mortar element!");

    // safety check
    if (!slaveelement->IsSlave() or masterelement->IsSlave())
      dserror("Something is wrong with the slave-master element pairing!");

    // construct slave-side and master-side location arrays
    DRT::Element::LocationArray la_slave(idiscret.NumDofSets());
    slaveelement->LocationVector(idiscret, la_slave, false);
    DRT::Element::LocationArray la_master(idiscret.NumDofSets());
    masterelement->LocationVector(idiscret, la_master, false);

    // initialize cell matrices and vectors
    strategy.InitCellMatricesAndVectors(la_slave, la_master);

    // evaluate current slave-side node
    EvaluateSlaveNode(idiscret, *slavenode, islavenodeslumpedareas[inode],
        (INPAR::SCATRA::ImplType)islavenodesimpltypes[inode], *slaveelement, *masterelement,
        la_slave, la_master, params, strategy.CellMatrix1(), strategy.CellMatrix2(),
        strategy.CellMatrix3(), strategy.CellMatrix4(), strategy.CellVector1(),
        strategy.CellVector2());

    // assemble cell matrices and vectors into system matrices and vectors
    strategy.AssembleCellMatricesAndVectors(la_slave, la_master, slavenode->Owner());
  }

  return;
}


/*--------------------------------------------------------------------------------------*
 | evaluate mortar elements                                                  fang 08/16 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EvaluateMortarElements(
    const Epetra_Map& ielecolmap,  //!< column map of mortar elements
    const Epetra_IntVector&
        ieleimpltypes,  //!< vector for physical implementation types of mortar elements
    const DRT::Discretization& idiscret,  //!< interface discretization
    const Teuchos::ParameterList&
        params,  //!< parameter list for evaluation of mortar integration cells
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix1,  //!< system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_rows,  //!< interface side associated with rows of system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_cols,  //!< interface side associated with columns of system matrix 1
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix2,  //!< system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_rows,  //!< interface side associated with rows of system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_cols,  //!< interface side associated with columns of system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix3,  //!< system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_rows,  //!< interface side associated with rows of system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_cols,  //!< interface side associated with columns of system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix4,  //!< system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_rows,  //!< interface side associated with rows of system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_cols,  //!< interface side associated with columns of system matrix 4
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,  //!< system vector 1
    const INPAR::S2I::InterfaceSides
        vector1_side,  //!< interface side associated with system vector 1
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,  //!< system vector 2
    const INPAR::S2I::InterfaceSides
        vector2_side  //!< interface side associated with system vector 2
    ) const
{
  // instantiate assembly strategy for mortar elements
  SCATRA::MortarCellAssemblyStrategy strategy(systemmatrix1, matrix1_side_rows, matrix1_side_cols,
      systemmatrix2, matrix2_side_rows, matrix2_side_cols, systemmatrix3, matrix3_side_rows,
      matrix3_side_cols, systemmatrix4, matrix4_side_rows, matrix4_side_cols, systemvector1,
      vector1_side, systemvector2, vector2_side);

  // loop over all mortar elements
  for (int ielement = 0; ielement < ielecolmap.NumMyElements(); ++ielement)
  {
    // extract current mortar element
    MORTAR::MortarElement* const element =
        dynamic_cast<MORTAR::MortarElement* const>(idiscret.gElement(ielecolmap.GID(ielement)));
    if (!element) dserror("Couldn't extract mortar element from mortar discretization!");

    // construct location array for current mortar element
    DRT::Element::LocationArray la(idiscret.NumDofSets());
    element->LocationVector(idiscret, la, false);

    // initialize element matrices and vectors
    strategy.InitCellMatricesAndVectors(la, la);  // second function argument only serves as dummy

    // evaluate current mortar element
    EvaluateMortarElement(idiscret, *element, (INPAR::SCATRA::ImplType)ieleimpltypes[ielement], la,
        params, strategy.CellMatrix1(), strategy.CellMatrix2(), strategy.CellMatrix3(),
        strategy.CellMatrix4(), strategy.CellVector1(), strategy.CellVector2());

    // assemble element matrices and vectors into system matrices and vectors
    strategy.AssembleCellMatricesAndVectors(la, la, -1);
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------*
 | provide instance of mortar cell evaluation class of particular slave-side discretization type
 fang 01/16 |
 *------------------------------------------------------------------------------------------------------------*/
SCATRA::MortarCellInterface* SCATRA::MortarCellFactory::MortarCellCalc(
    const INPAR::SCATRA::ImplType&
        impltype,  //!< physical implementation type of mortar integration cell
    const MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    const MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const std::string& disname  //!< name of interface discretization
)
{
  // extract number of slave-side degrees of freedom per node
  const int numdofpernode_slave = slaveelement.NumDofPerNode(*slaveelement.Nodes()[0]);

  switch (slaveelement.Shape())
  {
    case DRT::Element::tri3:
    {
      return MortarCellCalc<DRT::Element::tri3>(
          impltype, masterelement, couplingtype, lmside, numdofpernode_slave, disname);
      break;
    }

    case DRT::Element::quad4:
    {
      return MortarCellCalc<DRT::Element::quad4>(
          impltype, masterelement, couplingtype, lmside, numdofpernode_slave, disname);
      break;
    }

    default:
    {
      dserror("Invalid slave-side discretization type!");
      break;
    }
  }

  return nullptr;
}


/*-----------------------------------------------------------------------------------------------------------------------------*
 | provide instance of mortar cell evaluation class of particular slave-side and master-side
 discretization types   fang 01/16 |
 *-----------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS>
SCATRA::MortarCellInterface* SCATRA::MortarCellFactory::MortarCellCalc(
    const INPAR::SCATRA::ImplType&
        impltype,  //!< physical implementation type of mortar integration cell
    const MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const std::string& disname       //!< name of interface discretization
)
{
  // extract number of master-side degrees of freedom per node
  const int numdofpernode_master = masterelement.NumDofPerNode(*masterelement.Nodes()[0]);

  switch (masterelement.Shape())
  {
    case DRT::Element::tri3:
    {
      return MortarCellCalc<distypeS, DRT::Element::tri3>(
          impltype, couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case DRT::Element::quad4:
    {
      return MortarCellCalc<distypeS, DRT::Element::quad4>(
          impltype, couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    default:
    {
      dserror("Invalid master-side discretization type!");
      break;
    }
  }

  return nullptr;
}


/*--------------------------------------------------------------------------------------*
 | provide specific instance of mortar cell evaluation class                 fang 01/16 |
 *--------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellInterface* SCATRA::MortarCellFactory::MortarCellCalc(
    const INPAR::SCATRA::ImplType&
        impltype,  //!< physical implementation type of mortar integration cell
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname        //!< name of interface discretization
)
{
  // return instance of evaluation class for mortar integration cell depending on physical
  // implementation type
  switch (impltype)
  {
    case INPAR::SCATRA::impltype_std:
    {
      return SCATRA::MortarCellCalc<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case INPAR::SCATRA::impltype_elch_electrode:
    {
      return SCATRA::MortarCellCalcElch<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    {
      return SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    case INPAR::SCATRA::impltype_thermo_elch_electrode:
    {
      return SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Instance(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master, disname);
      break;
    }

    default:
    {
      dserror("Unknown physical implementation type of mortar integration cell!");
      break;
    }
  }

  return nullptr;
}


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::InitConvCheckStrategy()
{
  if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_petrov or
      couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov)
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyS2ILM(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStd(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyS2I::InitConvCheckStrategy


/*----------------------------------------------------------------------*
 | perform setup of scatra-scatra interface coupling         fang 10/14 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::SetupMeshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> conditions(0, nullptr);
  scatratimint_->Discretization()->GetCondition("S2ICoupling", conditions);
  slaveconditions_.clear();
  std::map<const int, DRT::Condition* const> masterconditions;
  for (unsigned icondition = 0; icondition < conditions.size(); ++icondition)
  {
    DRT::Condition* const condition = conditions[icondition];
    const int condid = condition->GetInt("ConditionID");
    if (condid < 0) dserror("Invalid condition ID %i for S2ICoupling Condition!", condid);

    switch (condition->GetInt("interface side"))
    {
      case INPAR::S2I::side_slave:
      {
        if (slaveconditions_.find(condid) == slaveconditions_.end())
          slaveconditions_.insert(std::pair<const int, DRT::Condition* const>(condid, condition));
        else
          dserror(
              "Cannot have multiple slave-side scatra-scatra interface coupling conditions with "
              "the same ID %i!",
              condid);
        break;
      }

      case INPAR::S2I::side_master:
      {
        if (masterconditions.find(condid) == masterconditions.end())
          masterconditions.insert(std::pair<const int, DRT::Condition* const>(condid, condition));
        else
          dserror(
              "Cannot have multiple master-side scatra-scatra interface coupling conditions with "
              "the same ID %i!",
              condid);
        break;
      }

      default:
      {
        dserror("Invalid scatra-scatra interface coupling condition!");
        break;
      }
    }
  }

  // determine type of mortar meshtying
  switch (couplingtype_)
  {
    // setup scatra-scatra interface coupling for interfaces with pairwise overlapping interface
    // nodes
    case INPAR::S2I::coupling_matching_nodes:
    {
      // overwrite IDs of master-side scatra-scatra interface coupling conditions with the value -1
      // to prevent them from being evaluated when calling EvaluateCondition on the discretization
      // TODO: this is somewhat unclean, because changing the conditions, makes calling
      // SetupMeshtying() twice
      //       invalid (which should not be necessary, but conceptually possible)
      for (std::map<const int, DRT::Condition* const>::iterator imastercondition =
               masterconditions.begin();
           imastercondition != masterconditions.end(); ++imastercondition)
        imastercondition->second->Add("ConditionID", -1);

      // initialize int vectors for global ids of slave and master interface nodes
      std::vector<int> islavenodegidvec;
      std::vector<int> imasternodegidvec;

      // fill vectors
      for (std::map<const int, DRT::Condition* const>::iterator islavecondition =
               slaveconditions_.begin();
           islavecondition != slaveconditions_.end(); ++islavecondition)
      {
        const std::vector<int>* islavenodegids = islavecondition->second->Nodes();

        for (unsigned islavenode = 0; islavenode < islavenodegids->size(); ++islavenode)
        {
          const int islavenodegid = (*islavenodegids)[islavenode];

          // insert global id of current node into associated vector only if node is owned by
          // current processor need to make sure that node is stored on current processor, otherwise
          // cannot resolve "->Owner()"
          if (scatratimint_->Discretization()->HaveGlobalNode(islavenodegid) and
              scatratimint_->Discretization()->gNode(islavenodegid)->Owner() ==
                  scatratimint_->Discretization()->Comm().MyPID())
            islavenodegidvec.push_back(islavenodegid);
        }
      }
      for (std::map<const int, DRT::Condition* const>::iterator imastercondition =
               masterconditions.begin();
           imastercondition != masterconditions.end(); ++imastercondition)
      {
        const std::vector<int>* imasternodegids = imastercondition->second->Nodes();

        for (unsigned imasternode = 0; imasternode < imasternodegids->size(); ++imasternode)
        {
          const int imasternodegid = (*imasternodegids)[imasternode];

          // insert global id of current node into associated vector only if node is owned by
          // current processor need to make sure that node is stored on current processor, otherwise
          // cannot resolve "->Owner()"
          if (scatratimint_->Discretization()->HaveGlobalNode(imasternodegid) and
              scatratimint_->Discretization()->gNode(imasternodegid)->Owner() ==
                  scatratimint_->Discretization()->Comm().MyPID())
            imasternodegidvec.push_back(imasternodegid);
        }
      }

      // remove potential duplicates from vectors
      std::sort(islavenodegidvec.begin(), islavenodegidvec.end());
      islavenodegidvec.erase(
          unique(islavenodegidvec.begin(), islavenodegidvec.end()), islavenodegidvec.end());
      std::sort(imasternodegidvec.begin(), imasternodegidvec.end());
      imasternodegidvec.erase(
          unique(imasternodegidvec.begin(), imasternodegidvec.end()), imasternodegidvec.end());

      // initialize non-mortar coupling adapter
      if (scatratimint_->NumScal() < 1) dserror("Number of transported scalars not correctly set!");
      icoup_ = Teuchos::rcp(new ADAPTER::Coupling());
      icoup_->SetupCoupling(*(scatratimint_->Discretization()), *(scatratimint_->Discretization()),
          imasternodegidvec, islavenodegidvec, scatratimint_->NumDofPerNode(), true, 1.e-8);

      // generate interior and interface maps
      Teuchos::RCP<Epetra_Map> ifullmap =
          LINALG::MergeMap(icoup_->SlaveDofMap(), icoup_->MasterDofMap());
      std::vector<Teuchos::RCP<const Epetra_Map>> imaps;
      imaps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()), *ifullmap));
      imaps.push_back(icoup_->SlaveDofMap());
      imaps.push_back(icoup_->MasterDofMap());

      // initialize global map extractor
      interfacemaps_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()), imaps));
      interfacemaps_->CheckForValidMapExtractor();

      // initialize interface vector
      // Although the interface vector only contains the transformed master interface dofs, we still
      // initialize it with the full DofRowMap of the discretization to make it work for parallel
      // computations.
      imasterphinp_ = LINALG::CreateVector(*(scatratimint_->Discretization()->DofRowMap()), false);

      // initialize auxiliary system matrices and associated transformation operators
      islavematrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*(icoup_->SlaveDofMap()), 81));
      imastermatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*(icoup_->SlaveDofMap()), 81));
      islavetomasterrowtransform_ = Teuchos::rcp(new LINALG::MatrixRowTransform);
      if (not slaveonly_)
      {
        islavetomastercoltransform_ = Teuchos::rcp(new LINALG::MatrixColTransform);
        islavetomasterrowcoltransform_ = Teuchos::rcp(new LINALG::MatrixRowColTransform);
      }

      // initialize auxiliary residual vector
      islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*(icoup_->SlaveDofMap())));

      break;
    }

    // setup scatra-scatra interface coupling for interfaces with non-overlapping interface nodes
    case INPAR::S2I::coupling_mortar_standard:
    case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
    case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
    case INPAR::S2I::coupling_mortar_condensed_petrov:
    case INPAR::S2I::coupling_mortar_condensed_bubnov:
    case INPAR::S2I::coupling_nts_standard:
    {
      // safety checks
      if (imortarredistribution_ and couplingtype_ != INPAR::S2I::coupling_mortar_standard)
        dserror(
            "Parallel redistribution only implemented for scatra-scatra interface coupling based "
            "on standard mortar approach!");
      if (DRT::INPUT::IntegralValue<INPAR::MORTAR::MeshRelocation>(
              DRT::Problem::Instance()->MortarCouplingParams(), "MESH_RELOCATION") !=
          INPAR::MORTAR::relocation_none)
        dserror("Mesh relocation not yet implemented for scatra-scatra interface coupling!");

      // initialize empty interface maps
      Teuchos::RCP<Epetra_Map> imastermap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->Discretization()->Comm()));
      Teuchos::RCP<Epetra_Map> islavemap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->Discretization()->Comm()));
      Teuchos::RCP<Epetra_Map> ifullmap =
          Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->Discretization()->Comm()));
      if (imortarredistribution_)
      {
        imastermap_ = Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->Discretization()->Comm()));
        islavemap_ = Teuchos::rcp(new Epetra_Map(0, 0, scatratimint_->Discretization()->Comm()));
      }

      // loop over all slave-side scatra-scatra interface coupling conditions
      for (std::map<const int, DRT::Condition* const>::iterator islavecondition =
               slaveconditions_.begin();
           islavecondition != slaveconditions_.end(); ++islavecondition)
      {
        // extract condition ID
        const int condid = islavecondition->first;

        // initialize maps for row nodes associated with current condition
        std::map<int, DRT::Node*> masternodes;
        std::map<int, DRT::Node*> slavenodes;

        // initialize maps for column nodes associated with current condition
        std::map<int, DRT::Node*> mastergnodes;
        std::map<int, DRT::Node*> slavegnodes;

        // initialize maps for elements associated with current condition
        std::map<int, Teuchos::RCP<DRT::Element>> masterelements;
        std::map<int, Teuchos::RCP<DRT::Element>> slaveelements;

        // extract current slave-side and associated master-side scatra-scatra interface coupling
        // conditions
        std::vector<DRT::Condition*> mastercondition(1, masterconditions[condid]);
        std::vector<DRT::Condition*> slavecondition(1, islavecondition->second);

        // fill maps
        DRT::UTILS::FindConditionObjects(*scatratimint_->Discretization(), masternodes,
            mastergnodes, masterelements, mastercondition);
        DRT::UTILS::FindConditionObjects(*scatratimint_->Discretization(), slavenodes, slavegnodes,
            slaveelements, slavecondition);

        // initialize mortar coupling adapter
        icoupmortar_[condid] = Teuchos::rcp(new ADAPTER::CouplingMortar());
        ADAPTER::CouplingMortar& icoupmortar = *icoupmortar_[condid];
        std::vector<int> coupleddof(scatratimint_->NumDofPerNode(), 1);
        icoupmortar.SetupInterface(scatratimint_->Discretization(), scatratimint_->Discretization(),
            coupleddof, mastergnodes, slavegnodes, masterelements, slaveelements,
            scatratimint_->Discretization()->Comm());

        // extract mortar interface
        MORTAR::MortarInterface& interface = *icoupmortar.Interface();

        // extract mortar discretization
        const DRT::Discretization& idiscret = interface.Discret();

        if (couplingtype_ != INPAR::S2I::coupling_nts_standard)
        {
          // provide each slave-side mortar element with material of corresponding parent element
          for (int iele = 0; iele < interface.SlaveColElements()->NumMyElements(); ++iele)
          {
            // determine global ID of current slave-side mortar element
            const int elegid = interface.SlaveColElements()->GID(iele);

            // add material
            idiscret.gElement(elegid)->SetMaterial(Teuchos::rcp_dynamic_cast<DRT::FaceElement>(
                islavecondition->second->Geometry()[elegid])
                                                       ->ParentElement()
                                                       ->Material()
                                                       ->Parameter()
                                                       ->Id());
          }

          // assign physical implementation type to each slave-side mortar element by copying the
          // physical implementation type of the corresponding parent volume element
          Epetra_IntVector impltypes_row(*interface.SlaveRowElements());
          for (int iele = 0; iele < interface.SlaveRowElements()->NumMyElements(); ++iele)
            impltypes_row[iele] = dynamic_cast<const DRT::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<const DRT::FaceElement>(
                    islavecondition->second->Geometry()[interface.SlaveRowElements()->GID(iele)])
                    ->ParentElement())
                                      ->ImplType();

          // perform parallel redistribution if desired
          if (imortarredistribution_ and idiscret.Comm().NumProc() > 1)
          {
            interface.InterfaceParams()
                .sublist("PARALLEL REDISTRIBUTION")
                .set<std::string>("PARALLEL_REDIST", DRT::Problem::Instance()
                                                         ->MortarCouplingParams()
                                                         .sublist("PARALLEL REDISTRIBUTION")
                                                         .get<std::string>("PARALLEL_REDIST"));
            interface.Redistribute();
            interface.FillComplete(true);
            interface.PrintParallelDistribution();
            interface.CreateSearchTree();
          }

          // generate mortar integration cells
          std::vector<Teuchos::RCP<MORTAR::IntCell>> imortarcells(0, Teuchos::null);
          icoupmortar.EvaluateGeometry(imortarcells);

          // assign physical implementation type to each mortar integration cell by copying the
          // physical implementation type of the corresponding slave-side mortar element
          Epetra_IntVector impltypes_col(*interface.SlaveColElements());
          LINALG::Export(impltypes_row, impltypes_col);
          imortarcells_[condid].resize(imortarcells.size());
          for (unsigned icell = 0; icell < imortarcells.size(); ++icell)
            imortarcells_[condid][icell] =
                std::pair<Teuchos::RCP<MORTAR::IntCell>, INPAR::SCATRA::ImplType>(
                    imortarcells[icell], static_cast<INPAR::SCATRA::ImplType>(
                                             impltypes_col[interface.SlaveColElements()->LID(
                                                 imortarcells[icell]->GetSlaveId())]));
        }

        else
        {
          // match slave-side and master-side elements at mortar interface
          switch (interface.SearchAlg())
          {
            case INPAR::MORTAR::search_bfele:
            {
              interface.EvaluateSearchBruteForce(interface.SearchParam());
              break;
            }

            case INPAR::MORTAR::search_binarytree:
            {
              interface.EvaluateSearchBinarytree();
              break;
            }

            default:
            {
              dserror("Invalid search algorithm!");
              break;
            }
          }

          // evaluate normal vectors associated with slave-side nodes
          interface.EvaluateNodalNormals();

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
          islavenodesimpltypes->PutValue(INPAR::SCATRA::impltype_undefined);

          // loop over all slave-side nodes
          for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
          {
            // extract slave-side node
            MORTAR::MortarNode* const slavenode =
                dynamic_cast<MORTAR::MortarNode*>(idiscret.gNode(noderowmap_slave.GID(inode)));
            if (!slavenode)
              dserror("Couldn't extract slave-side mortar node from mortar discretization!");

            // find associated master-side elements
            std::vector<MORTAR::MortarElement*> masterelements(0, nullptr);
            interface.FindMEles(*slavenode, masterelements);

            // loop over all master-side elements
            for (unsigned imasterelement = 0; imasterelement < masterelements.size();
                 ++imasterelement)
            {
              // extract master-side element
              MORTAR::MortarElement* const masterelement = masterelements[imasterelement];

              // project slave-side node onto master-side element
              double coordinates_master[2] = {};
              double dummy(0.);
              MORTAR::MortarProjector::Impl(*masterelement)
                  ->ProjectGaussPointAuxn3D(slavenode->X(), slavenode->MoData().n(), *masterelement,
                      coordinates_master, dummy);

              // check whether projected node lies inside master-side element
              if (masterelement->Shape() == DRT::Element::quad4)
              {
                if (coordinates_master[0] < -1. - ntsprojtol_ or
                    coordinates_master[1] < -1. - ntsprojtol_ or
                    coordinates_master[0] > 1. + ntsprojtol_ or
                    coordinates_master[1] > 1. + ntsprojtol_)
                  // projected node lies outside master-side element
                  continue;
              }

              else if (masterelement->Shape() == DRT::Element::tri3)
              {
                if (coordinates_master[0] < -ntsprojtol_ or coordinates_master[1] < -ntsprojtol_ or
                    coordinates_master[0] + coordinates_master[1] > 1. + 2 * ntsprojtol_)
                  // projected node lies outside master-side element
                  continue;
              }

              else
                dserror("Invalid discretization type of master-side element!");

              // projected node lies inside master-side element
              (*islavenodestomasterelements)[inode] = masterelement->Id();
              break;
            }

            // safety check
            if ((*islavenodestomasterelements)[inode] == -1)
              dserror("Couldn't match slave-side node with master-side element!");

            // determine physical implementation type of slave-side node based on first associated
            // element
            (*islavenodesimpltypes)[inode] = dynamic_cast<DRT::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<DRT::FaceElement>(
                    islavecondition->second->Geometry()[slavenode->Elements()[0]->Id()])
                    ->ParentElement())
                                                 ->ImplType();
          }

          // extract slave-side elerowmap
          const Epetra_Map& elecolmap_slave = *interface.SlaveColElements();

          // initialize vector for physical implementation types of slave-side elements
          Epetra_IntVector islaveelementsimpltypes(elecolmap_slave, false);
          islaveelementsimpltypes.PutValue(INPAR::SCATRA::impltype_undefined);

          // loop over all slave-side elements
          for (int ielement = 0; ielement < elecolmap_slave.NumMyElements(); ++ielement)
            // determine physical implementation type of current slave-side element
            islaveelementsimpltypes[ielement] = dynamic_cast<DRT::ELEMENTS::Transport*>(
                Teuchos::rcp_dynamic_cast<DRT::FaceElement>(
                    islavecondition->second->Geometry()[elecolmap_slave.GID(ielement)])
                    ->ParentElement())
                                                    ->ImplType();

          // create parameter list for slave-side elements
          Teuchos::ParameterList eleparams;

          // set action for slave-side elements
          eleparams.set<int>("action", INPAR::S2I::evaluate_nodal_area_fractions);

          // compute vector for lumped interface area fractions associated with slave-side nodes
          const Epetra_Map& dofrowmap_slave = *interface.SlaveRowDofs();
          Teuchos::RCP<Epetra_Vector> islavenodeslumpedareas_dofvector =
              LINALG::CreateVector(dofrowmap_slave);
          EvaluateMortarElements(elecolmap_slave, islaveelementsimpltypes, idiscret, eleparams,
              Teuchos::null, INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
              INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
              INPAR::S2I::side_undefined, INPAR::S2I::side_undefined, Teuchos::null,
              INPAR::S2I::side_undefined, INPAR::S2I::side_undefined,
              islavenodeslumpedareas_dofvector, INPAR::S2I::side_slave, Teuchos::null,
              INPAR::S2I::side_undefined);

          // transform map of result vector
          Teuchos::RCP<Epetra_Vector>& islavenodeslumpedareas = islavenodeslumpedareas_[condid];
          islavenodeslumpedareas = LINALG::CreateVector(noderowmap_slave);
          for (int inode = 0; inode < noderowmap_slave.NumMyElements(); ++inode)
            (*islavenodeslumpedareas)[inode] =
                (*islavenodeslumpedareas_dofvector)[dofrowmap_slave.LID(
                    idiscret.Dof(idiscret.gNode(noderowmap_slave.GID(inode)), 0))];
        }

        // build interface maps
        imastermap = LINALG::MergeMap(imastermap, interface.MasterRowDofs(), false);
        islavemap = LINALG::MergeMap(islavemap, interface.SlaveRowDofs(), false);
        ifullmap = LINALG::MergeMap(ifullmap,
            LINALG::MergeMap(interface.MasterRowDofs(), interface.SlaveRowDofs(), false), false);
        if (imortarredistribution_)
        {
          imastermap_ = LINALG::MergeMap(imastermap_, icoupmortar.MasterDofMap(), false);
          islavemap_ = LINALG::MergeMap(islavemap_, icoupmortar.SlaveDofMap(), false);
        }
      }

      // generate interior and interface maps
      std::vector<Teuchos::RCP<const Epetra_Map>> imaps;
      imaps.push_back(LINALG::SplitMap(*(scatratimint_->Discretization()->DofRowMap()), *ifullmap));
      imaps.push_back(islavemap);
      imaps.push_back(imastermap);

      // initialize global map extractor
      interfacemaps_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*(scatratimint_->Discretization()->DofRowMap()), imaps));
      interfacemaps_->CheckForValidMapExtractor();

      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_slave or couplingtype_ == INPAR::S2I::coupling_nts_standard)
      {
        // initialize auxiliary system matrix for slave side
        islavematrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(1), 81));

        // initialize auxiliary residual vector for slave side
        islaveresidual_ = Teuchos::rcp(new Epetra_Vector(*interfacemaps_->Map(1)));
      }

      if (couplingtype_ == INPAR::S2I::coupling_mortar_standard or
          lmside_ == INPAR::S2I::side_master or couplingtype_ == INPAR::S2I::coupling_nts_standard)
      {
        // initialize auxiliary system matrix for master side
        imastermatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
            *interfacemaps_->Map(2), 81, true, false, LINALG::SparseMatrix::FE_MATRIX));

        // initialize auxiliary residual vector for master side
        imasterresidual_ = Teuchos::rcp(new Epetra_FEVector(*interfacemaps_->Map(2)));
      }

      switch (couplingtype_)
      {
        case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
        case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
        case INPAR::S2I::coupling_mortar_condensed_petrov:
        case INPAR::S2I::coupling_mortar_condensed_bubnov:
        {
          if (lmside_ == INPAR::S2I::side_slave)
          {
            D_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(1), 81));
            M_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(1), 81));
            if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov or
                couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
              E_ = Teuchos::rcp(new LINALG::SparseMatrix(*interfacemaps_->Map(1), 81));
          }
          else
          {
            D_ = Teuchos::rcp(new LINALG::SparseMatrix(
                *interfacemaps_->Map(2), 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
            M_ = Teuchos::rcp(new LINALG::SparseMatrix(
                *interfacemaps_->Map(2), 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
            if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov or
                couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
              E_ = Teuchos::rcp(new LINALG::SparseMatrix(
                  *interfacemaps_->Map(2), 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
          }

          // loop over all scatra-scatra coupling interfaces
          for (std::map<const int, DRT::Condition* const>::iterator islavecondition =
                   slaveconditions_.begin();
               islavecondition != slaveconditions_.end(); ++islavecondition)
          {
            // create parameter list for mortar integration cells
            Teuchos::ParameterList params;

            // add current condition to parameter list
            params.set<DRT::Condition*>("condition", islavecondition->second);

            // set action
            params.set<int>("action", INPAR::S2I::evaluate_mortar_matrices);

            // evaluate mortar integration cells at current interface
            EvaluateMortarCells(icoupmortar_[islavecondition->first]->Interface()->Discret(),
                params, D_,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_slave
                                                  : INPAR::S2I::side_master,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_slave
                                                  : INPAR::S2I::side_master,
                M_,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_slave
                                                  : INPAR::S2I::side_master,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_master
                                                  : INPAR::S2I::side_slave,
                E_,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_slave
                                                  : INPAR::S2I::side_master,
                lmside_ == INPAR::S2I::side_slave ? INPAR::S2I::side_slave
                                                  : INPAR::S2I::side_master,
                Teuchos::null, INPAR::S2I::side_undefined, INPAR::S2I::side_undefined,
                Teuchos::null, INPAR::S2I::side_undefined, Teuchos::null,
                INPAR::S2I::side_undefined);
          }

          // finalize mortar matrices D, M, and E
          D_->Complete();
          if (lmside_ == INPAR::S2I::side_slave)
            M_->Complete(*interfacemaps_->Map(2), *interfacemaps_->Map(1));
          else
            M_->Complete(*interfacemaps_->Map(1), *interfacemaps_->Map(2));
          if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov or
              couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
            E_->Complete();

          switch (couplingtype_)
          {
            case INPAR::S2I::coupling_mortar_condensed_petrov:
            case INPAR::S2I::coupling_mortar_condensed_bubnov:
            {
              // set up mortar projector P
              Teuchos::RCP<Epetra_Vector> D_diag(Teuchos::null);
              if (lmside_ == INPAR::S2I::side_slave)
                D_diag = LINALG::CreateVector(*interfacemaps_->Map(1));
              else
                D_diag = LINALG::CreateVector(*interfacemaps_->Map(2));
              if (D_->ExtractDiagonalCopy(*D_diag))
                dserror("Couldn't extract main diagonal from mortar matrix D!");
              if (D_diag->Reciprocal(*D_diag))
                dserror("Couldn't invert main diagonal entries of mortar matrix D!");
              ;
              P_ = Teuchos::rcp(new LINALG::SparseMatrix(*M_));
              if (P_->LeftScale(*D_diag)) dserror("Setup of mortar projector P failed!");

              // free memory
              if (!slaveonly_)
              {
                D_ = Teuchos::null;
                M_ = Teuchos::null;
              }

              if (couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
              {
                // set up mortar projector Q
                Q_ = Teuchos::rcp(new LINALG::SparseMatrix(*E_));
                if (Q_->LeftScale(*D_diag)) dserror("Setup of mortar projector Q failed!");

                // free memory
                E_ = Teuchos::null;
              }

              break;
            }

            case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
            case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
            {
              // determine number of Lagrange multiplier dofs owned by each processor
              const Epetra_Comm& comm(scatratimint_->Discretization()->Comm());
              const int numproc(comm.NumProc());
              const int mypid(comm.MyPID());
              std::vector<int> localnumlmdof(numproc, 0);
              std::vector<int> globalnumlmdof(numproc, 0);
              if (lmside_ == INPAR::S2I::side_slave)
                localnumlmdof[mypid] = interfacemaps_->Map(1)->NumMyElements();
              else
                localnumlmdof[mypid] = interfacemaps_->Map(2)->NumMyElements();
              comm.SumAll(&localnumlmdof[0], &globalnumlmdof[0], numproc);

              // for each processor, determine offset of minimum Lagrange multiplier dof GID w.r.t.
              // maximum standard dof GID
              int offset(0);
              for (int ipreviousproc = 0; ipreviousproc < mypid; ++ipreviousproc)
                offset += globalnumlmdof[ipreviousproc];

              // for each processor, determine Lagrange multiplier dof GIDs
              std::vector<int> lmdofgids(globalnumlmdof[mypid], 0);
              for (int lmdoflid = 0; lmdoflid < globalnumlmdof[mypid]; ++lmdoflid)
                lmdofgids[lmdoflid] =
                    scatratimint_->DofRowMap()->MaxAllGID() + 1 + offset + lmdoflid;

              // build Lagrange multiplier dofrowmap
              const Teuchos::RCP<Epetra_Map> lmdofrowmap =
                  Teuchos::rcp(new Epetra_Map(-1, (int)lmdofgids.size(), &lmdofgids[0], 0, comm));

              // initialize vectors associated with Lagrange multiplier dofs
              lm_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));
              lmresidual_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));
              lmincrement_ = Teuchos::rcp(new Epetra_Vector(*lmdofrowmap));

              // initialize extended map extractor
              Teuchos::RCP<Epetra_Map> extendedmap = LINALG::MergeMap(
                  *(scatratimint_->Discretization()->DofRowMap()), *lmdofrowmap, false);
              extendedmaps_ = Teuchos::rcp(new LINALG::MapExtractor(
                  *extendedmap, lmdofrowmap, scatratimint_->Discretization()->DofRowMap()));
              extendedmaps_->CheckForValidMapExtractor();

              // transform range map of mortar matrices D and M from slave-side dofrowmap to
              // Lagrange multiplier dofrowmap
              D_ = MORTAR::MatrixRowTransformGIDs(D_, lmdofrowmap);
              M_ = MORTAR::MatrixRowTransformGIDs(M_, lmdofrowmap);

              if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_petrov)
              {
                // transform domain map of mortar matrix D from slave-side dofrowmap to Lagrange
                // multiplier dofrowmap and store transformed matrix as mortar matrix E
                E_ = MORTAR::MatrixColTransformGIDs(D_, lmdofrowmap);
              }
              else
              {
                // transform domain and range maps of mortar matrix E from slave-side dofrowmap to
                // Lagrange multiplier dofrowmap
                E_ = MORTAR::MatrixRowColTransformGIDs(E_, lmdofrowmap, lmdofrowmap);
              }

              break;
            }

            default:
            {
              dserror("Invalid type of mortar meshtying!");
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
      dserror("Type of mortar meshtying for scatra-scatra interface coupling not recognized!");
      break;
    }
  }

  // further initializations depending on type of global system matrix
  switch (matrixtype_)
  {
    case INPAR::SCATRA::MatrixType::sparse:
    {
      // nothing needs to be done in this case
      break;
    }
    case INPAR::SCATRA::MatrixType::block_condition:
    case INPAR::SCATRA::MatrixType::block_condition_dof:
    case INPAR::SCATRA::MatrixType::block_geometry:
    {
      // safety check
      if (!scatratimint_->Solver()->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // initialize map extractors associated with blocks of global system matrix
      BuildBlockMapExtractors();

      break;
    }
    default:
    {
      dserror(
          "%i is not a valid 'SCATRA::MatrixType'. Set a valid 'SCATRA::MatrixType' in your input "
          "file!",
          static_cast<int>(matrixtype_));
      break;
    }
  }

  // extract boundary condition for scatra-scatra interface layer growth
  const DRT::Condition* const condition =
      scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth");

  // setup evaluation of scatra-scatra interface layer growth if applicable
  if (condition)
  {
    // perform setup depending on evaluation method
    switch (intlayergrowth_evaluation_)
    {
      case INPAR::S2I::growth_evaluation_monolithic:
      case INPAR::S2I::growth_evaluation_semi_implicit:
      {
        // extract map associated with scatra-scatra interface layer thicknesses
        const Teuchos::RCP<const Epetra_Map>& dofrowmap_growth = scatratimint_->DofRowMap(2);

        // initialize state vector of discrete scatra-scatra interface layer thicknesses at time n
        growthn_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));

        // additional initializations for monolithic solution approach
        if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
        {
          // initialize extended map extractor
          const Epetra_Map* const dofrowmap_scatra = scatratimint_->Discretization()->DofRowMap();
          extendedmaps_ = Teuchos::rcp(new LINALG::MapExtractor(
              *LINALG::MergeMap(*dofrowmap_scatra, *dofrowmap_growth, false),
              scatratimint_->DofRowMap(2), dofrowmap_scatra));
          extendedmaps_->CheckForValidMapExtractor();

          // initialize additional state vectors
          growthnp_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthdtn_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthdtnp_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthhist_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthresidual_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));
          growthincrement_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_growth, true));

          // initialize map extractors and global matrix blocks
          growthgrowthblock_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_growth, 81));
          switch (matrixtype_)
          {
            case INPAR::SCATRA::MatrixType::sparse:
            {
              // initialize extended map extractor associated with blocks of global system matrix
              extendedblockmaps_ = extendedmaps_;

              scatragrowthblock_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_scatra,
                  81));  // We actually don't really need the entire scalar transport dofrowmap
                         // here, but only a submap associated with all (slave-side and master-side)
                         // interfacial degrees of freedom. However, this will later cause an error
                         // in debug mode when assigning the scatra-growth matrix block to the
                         // global system matrix in the Solve() routine.
              growthscatrablock_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_growth, 81));

              break;
            }

            case INPAR::SCATRA::MatrixType::block_condition:
            case INPAR::SCATRA::MatrixType::block_condition_dof:
            {
              // initialize map extractor associated with all degrees of freedom for scatra-scatra
              // interface layer growth
              blockmapgrowth_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*dofrowmap_growth,
                  std::vector<Teuchos::RCP<const Epetra_Map>>(1, dofrowmap_growth)));
              blockmapgrowth_->CheckForValidMapExtractor();

              // initialize extended map extractor associated with blocks of global system matrix
              const unsigned nblockmaps = scatratimint_->BlockMaps().NumMaps();
              std::vector<Teuchos::RCP<const Epetra_Map>> extendedblockmaps(
                  nblockmaps + 1, Teuchos::null);
              for (unsigned iblockmap = 0; iblockmap < nblockmaps; ++iblockmap)
                extendedblockmaps[iblockmap] = scatratimint_->BlockMaps().Map(iblockmap);
              extendedblockmaps[nblockmaps] = dofrowmap_growth;
              extendedblockmaps_ = Teuchos::rcp(
                  new LINALG::MultiMapExtractor(*extendedmaps_->FullMap(), extendedblockmaps));
              extendedblockmaps_->CheckForValidMapExtractor();

              scatragrowthblock_ =
                  Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                      *blockmapgrowth_, scatratimint_->BlockMaps(), 81, false, true));
              growthscatrablock_ =
                  Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                      scatratimint_->BlockMaps(), *blockmapgrowth_, 81, false, true));

              break;
            }

            default:
            {
              dserror(
                  "Type of global system matrix for scatra-scatra interface coupling involving "
                  "interface layer growth not recognized!");
              break;
            }
          }

          // initialize extended system matrix including rows and columns associated with
          // scatra-scatra interface layer thickness variables
          extendedsystemmatrix_ =
              Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                  *extendedblockmaps_, *extendedblockmaps_));
        }

        // loop over all boundary conditions for scatra-scatra interface coupling
        for (unsigned icond = 0; icond < conditions.size(); ++icond)
        {
          // check whether current boundary condition is associated with boundary condition for
          // scatra-scatra interface layer growth
          if (conditions[icond]->GetInt("ConditionID") == condition->GetInt("ConditionID"))
            // copy conductivity parameter
            conditions[icond]->Add("conductivity", condition->GetDouble("conductivity"));
        }

        // extract initial scatra-scatra interface layer thickness from condition
        const double initthickness = condition->GetDouble("initial thickness");

        // extract nodal cloud from condition
        const std::vector<int>* nodegids = condition->Nodes();

        // loop over all nodes
        for (unsigned inode = 0; inode < nodegids->size(); ++inode)
        {
          // extract global ID of current node
          const int nodegid = (*nodegids)[inode];

          // process only nodes stored by current processor
          if (scatratimint_->Discretization()->HaveGlobalNode(nodegid))
          {
            // extract current node
            const DRT::Node* const node = scatratimint_->Discretization()->gNode(nodegid);

            // process only nodes owned by current processor
            if (node->Owner() == scatratimint_->Discretization()->Comm().MyPID())
            {
              // extract local ID of scatra-scatra interface layer thickness variable associated
              // with current node
              const int doflid_growth = scatratimint_->Discretization()->DofRowMap(2)->LID(
                  scatratimint_->Discretization()->Dof(2, node, 0));
              if (doflid_growth < 0)
                dserror(
                    "Couldn't extract local ID of scatra-scatra interface layer thickness "
                    "variable!");

              // set initial value
              (*growthn_)[doflid_growth] = initthickness;
            }  // nodes owned by current processor
          }    // nodes stored by current processor
        }      // loop over all nodes

        // copy initial state
        if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
          growthnp_->Update(1., *growthn_, 0.);

        break;
      }

      default:
      {
        dserror(
            "Unknown evaluation method for scatra-scatra interface coupling involving interface "
            "layer growth!");
        break;
      }
    }
  }

  // initialize vectors for row and column sums of global system matrix if necessary
  const Epetra_Map& dofrowmap =
      intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic
          ? *extendedmaps_->FullMap()
          : *scatratimint_->Discretization()->DofRowMap();
  if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_full or
      scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_maindiag or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
    invrowsums_ = Teuchos::rcp(new Epetra_Vector(dofrowmap, false));
  if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::columns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::columns_maindiag or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
    invcolsums_ = Teuchos::rcp(new Epetra_Vector(dofrowmap, false));
}  // SCATRA::MeshtyingStrategyS2I::SetupMeshtying


/*----------------------------------------------------------------------*
 | compute time derivatives of discrete state variables      fang 01/17 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeTimeDerivative() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
  {
    // compute inverse time factor 1./(theta*dt)
    const double timefac_inverse =
        1. / scatratimint_->ScatraTimeParameterList()->get<double>("time factor");

    // compute state vector of time derivatives of discrete scatra-scatra interface layer
    // thicknesses
    growthdtnp_->Update(timefac_inverse, *growthnp_, -timefac_inverse, *growthhist_, 0.);
  }

  return;
}


/*------------------------------------------------------------------------------------------*
 | update solution after convergence of the nonlinear Newton-Raphson iteration   fang 01/17 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::Update() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
  {
    // update state vectors
    growthn_->Update(1., *growthnp_, 0.);
    growthdtn_->Update(1., *growthdtnp_, 0.);
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 | set general parameters for element evaluation                            fang 01/17 |
 *-------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::SetElementGeneralParameters(
    Teuchos::ParameterList& parameters) const
{
  // add local Newton-Raphson convergence tolerance for scatra-scatra interface layer growth to
  // parameter list
  parameters.set<double>("intlayergrowth_convtol", intlayergrowth_convtol_);

  // add maximum number of local Newton-Raphson iterations for scatra-scatra interface layer growth
  // to parameter list
  parameters.set<unsigned>("intlayergrowth_itemax", intlayergrowth_itemax_);

  return;
}

/*----------------------------------------------------------------------*
 | set interface kinetic model specific parameters        civaner 08/19 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::SetConditionSpecificScaTraParameters(
    DRT::Condition& s2icondition) const
{
  Teuchos::ParameterList conditionparams;

  // get kinetic model and condition type
  const int kineticmodel = s2icondition.GetInt("kinetic model");
  const DRT::Condition::ConditionType conditiontype = s2icondition.Type();

  // set action, kinetic model, condition type and numscal
  conditionparams.set<int>("action", SCATRA::set_scatra_ele_boundary_parameter);
  conditionparams.set<int>("kinetic model", kineticmodel);
  conditionparams.set<DRT::Condition::ConditionType>("condition type", conditiontype);

  // set the condition type specific parameters
  switch (conditiontype)
  {
    case DRT::Condition::ConditionType::S2ICoupling:
    {
      // set the kinetic model specific parameters
      switch (kineticmodel)
      {
        case INPAR::S2I::kinetics_constperm:
        {
          conditionparams.set<int>("numscal", s2icondition.GetInt("numscal"));
          conditionparams.set<std::vector<double>*>(
              "permeabilities", s2icondition.GetMutable<std::vector<double>>("permeabilities"));
          break;
        }

        case INPAR::S2I::kinetics_constantinterfaceresistance:
        {
          conditionparams.set<double>("resistance", s2icondition.GetDouble("resistance"));
          break;
        }

        case INPAR::S2I::kinetics_nointerfaceflux:
        {
          // do nothing
          break;
        }

        case INPAR::S2I::kinetics_butlervolmer:
        case INPAR::S2I::kinetics_butlervolmerreduced:
        case INPAR::S2I::kinetics_butlervolmerpeltier:
        case INPAR::S2I::kinetics_butlervolmerresistance:
        case INPAR::S2I::kinetics_butlervolmerreducedwithresistance:
        {
          conditionparams.set<int>("numscal", s2icondition.GetInt("numscal"));
          conditionparams.set<std::vector<int>*>(
              "stoichiometries", s2icondition.GetMutable<std::vector<int>>("stoichiometries"));
          conditionparams.set<int>("numelectrons", s2icondition.GetInt("e-"));
          conditionparams.set<double>("k_r", s2icondition.GetDouble("k_r"));
          conditionparams.set<double>("alpha_a", s2icondition.GetDouble("alpha_a"));
          conditionparams.set<double>("alpha_c", s2icondition.GetDouble("alpha_c"));

          if (kineticmodel == INPAR::S2I::kinetics_butlervolmerpeltier)
            conditionparams.set<double>("peltier", s2icondition.GetDouble("peltier"));

          if (kineticmodel == INPAR::S2I::kinetics_butlervolmerresistance or
              kineticmodel == INPAR::S2I::kinetics_butlervolmerreducedwithresistance)
          {
            conditionparams.set<double>("resistance", s2icondition.GetDouble("resistance"));
            conditionparams.set<double>(
                "CONVTOL_IMPLBUTLERVOLMER", s2icondition.GetDouble("CONVTOL_IMPLBUTLERVOLMER"));
            conditionparams.set<double>(
                "ITEMAX_IMPLBUTLERVOLMER", s2icondition.GetDouble("ITEMAX_IMPLBUTLERVOLMER"));
          }
          break;
        }

        default:
        {
          dserror("Not implemented for this kinetic model: %i", kineticmodel);
          break;
        }
      }
      break;
    }

    case DRT::Condition::ConditionType::S2ICouplingGrowth:
    {
      // set the kinetic model specific parameters
      switch (kineticmodel)
      {
        case INPAR::S2I::growth_kinetics_butlervolmer:
        {
          conditionparams.set<int>("numscal", s2icondition.GetInt("numscal"));
          conditionparams.set<std::vector<int>*>(
              "stoichiometries", s2icondition.GetMutable<std::vector<int>>("stoichiometries"));
          conditionparams.set<int>("numelectrons", s2icondition.GetInt("e-"));
          conditionparams.set<double>("k_r", s2icondition.GetDouble("k_r"));
          conditionparams.set<double>("alpha_a", s2icondition.GetDouble("alpha_a"));
          conditionparams.set<double>("alpha_c", s2icondition.GetDouble("alpha_c"));
          conditionparams.set<double>("density", s2icondition.GetDouble("density"));
          conditionparams.set<double>("molar mass", s2icondition.GetDouble("molar mass"));
          conditionparams.set<double>("regpar", s2icondition.GetDouble("regularization parameter"));
          conditionparams.set<int>("regtype", s2icondition.GetInt("regularization type"));
          conditionparams.set<double>("conductivity", s2icondition.GetDouble("conductivity"));
          break;
        }

        default:
        {
          dserror("Not implemented for this kinetic model: %i", kineticmodel);
          break;
        }
      }
      break;
    }

    default:
    {
      dserror("Not implemented for this condition type: %i", conditiontype);
      break;
    }
  }

  // call standard loop over elements
  scatratimint_->Discretization()->Evaluate(
      conditionparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------------------------------------------------------------------------------*
 | compute history vector, i.e., the history part of the right-hand side vector with all
 contributions from the previous time step   fang 01/17 |
 *----------------------------------------------------------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::SetOldPartOfRHS() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
  {
    // compute factor dt*(1-theta)
    const double factor =
        scatratimint_->Dt() - scatratimint_->ScatraTimeParameterList()->get<double>("time factor");

    // compute history vector
    growthhist_->Update(1., *growthn_, factor, *growthdtn_, 0.);
  }

  return;
}


/*-----------------------------------------------------------------------*
 | output restart information                                 fang 01/17 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::OutputRestart() const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // output state vector of discrete scatra-scatra interface layer thicknesses
    scatratimint_->DiscWriter()->WriteVector("growthn", growthn_);

    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
      // output state vector of time derivatives of discrete scatra-scatra interface layer
      // thicknesses
      scatratimint_->DiscWriter()->WriteVector("growthdtn", growthdtn_);
  }

  return;
}


/*-----------------------------------------------------------------------*
 | read restart data                                          fang 01/17 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ReadRestart(const int step,  //!< restart step
    Teuchos::RCP<IO::InputControl> input                        //!< control file manager
    ) const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // initialize reader
    Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
    if (input == Teuchos::null)
      reader = Teuchos::rcp(new IO::DiscretizationReader(scatratimint_->Discretization(), step));
    else
      reader =
          Teuchos::rcp(new IO::DiscretizationReader(scatratimint_->Discretization(), input, step));

    // read state vector of discrete scatra-scatra interface layer thicknesses
    reader->ReadVector(growthn_, "growthn");

    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
    {
      // read state vector of time derivatives of discrete scatra-scatra interface layer thicknesses
      reader->ReadVector(growthdtn_, "growthdtn");

      // copy restart state
      growthnp_->Update(1., *growthn_, 0.);
      growthdtnp_->Update(1., *growthdtn_, 0.);
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------------*
 | return mortar interface discretization associated with particular condition ID   fang 01/17 |
 *---------------------------------------------------------------------------------------------*/
DRT::Discretization& SCATRA::MeshtyingStrategyS2I::MortarDiscretization(const int& condid) const
{
  return icoupmortar_.at(condid)->Interface()->Discret();
}


/*-----------------------------------------------------------------------*
 | output solution for post-processing                        fang 01/17 |
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::Output() const
{
  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // extract relevant state vector of discrete scatra-scatra interface layer thicknesses based on
    // map of scatra-scatra interface layer thickness variables
    const Epetra_Vector& growth =
        intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic ? *growthnp_
                                                                               : *growthn_;

    // for proper output, initialize target state vector of discrete scatra-scatra interface layer
    // thicknesses based on map of row nodes
    const Teuchos::RCP<Epetra_Vector> intlayerthickness =
        Teuchos::rcp(new Epetra_Vector(*scatratimint_->Discretization()->NodeRowMap(), true));

    // extract boundary condition for scatra-scatra interface layer growth
    const DRT::Condition* const condition =
        scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth");

    // extract nodal cloud from condition
    const std::vector<int>* nodegids = condition->Nodes();

    // loop over all nodes
    for (unsigned inode = 0; inode < nodegids->size(); ++inode)
    {
      // extract global ID of current node
      const int nodegid((*nodegids)[inode]);

      // process only nodes stored by current processor
      if (scatratimint_->Discretization()->HaveGlobalNode(nodegid))
      {
        // extract current node
        const DRT::Node* const node = scatratimint_->Discretization()->gNode(nodegid);

        // process only nodes owned by current processor
        if (node->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        {
          // extract local ID of current node
          const int nodelid = scatratimint_->Discretization()->NodeRowMap()->LID(nodegid);
          if (nodelid < 0) dserror("Couldn't extract local node ID!");

          // extract local ID of scatra-scatra interface layer thickness variable associated with
          // current node
          const int doflid_growth = scatratimint_->Discretization()->DofRowMap(2)->LID(
              scatratimint_->Discretization()->Dof(2, node, 0));
          if (doflid_growth < 0)
            dserror(
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

  return;
}


/*---------------------------------------------------------------------------------------------------*
 | explicit predictor step to obtain better starting value for Newton-Raphson iteration   fang 01/17
 |
 *---------------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ExplicitPredictor() const
{
  // only relevant for monolithic evaluation of scatra-scatra interface layer growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
    // predict state vector of discrete scatra-scatra interface layer thicknesses at time n+1
    growthnp_->Update(scatratimint_->Dt(), *growthdtn_, 1.);

  return;
}


/*----------------------------------------------------------------------------*
 | extract selected rows from a sparse matrix                      fang 02/17 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ExtractMatrixRows(
    const LINALG::SparseMatrix& matrix,  //!< source matrix
    LINALG::SparseMatrix& rows,          //!< destination matrix
    const Epetra_Map& rowmap             //!< map of matrix rows to be extracted
)
{
  // safety check
  if (rows.Filled())
    dserror("Source matrix rows cannot be extracted into filled destination matrix!");

  // loop over all source matrix rows to be extracted
  for (int doflid = 0; doflid < rowmap.NumMyElements(); ++doflid)
  {
    // determine global ID of current matrix row
    const int dofgid = rowmap.GID(doflid);
    if (dofgid < 0) dserror("Couldn't find local ID %d in map!", doflid);

    // extract current matrix row from source matrix
    const int length = matrix.EpetraMatrix()->NumGlobalEntries(dofgid);
    int numentries(0);
    std::vector<double> values(length, 0.);
    std::vector<int> indices(length, 0);
    if (matrix.EpetraMatrix()->ExtractGlobalRowCopy(
            dofgid, length, numentries, &values[0], &indices[0]))
      dserror("Cannot extract matrix row with global ID %d from source matrix!", dofgid);

    // copy current source matrix row into destination matrix
    if (rows.EpetraMatrix()->InsertGlobalValues(dofgid, numentries, &values[0], &indices[0]) < 0)
      dserror("Cannot insert matrix row with global ID %d into destination matrix!", dofgid);
  }

  return;
}


/*---------------------------------------------------------------------------*
 | provide global state vectors for element evaluation            fang 01/17 |
 *---------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::AddTimeIntegrationSpecificVectors() const
{
  // only relevant for scatra-scatra interface coupling with pairwise coinciding interface nodes
  if (couplingtype_ == INPAR::S2I::coupling_matching_nodes)
  {
    // add state vector containing master-side scatra degrees of freedom to scatra discretization
    interfacemaps_->InsertVector(
        icoup_->MasterToSlave(interfacemaps_->ExtractVector(*(scatratimint_->Phiafnp()), 2)), 1,
        imasterphinp_);
    scatratimint_->Discretization()->SetState("imasterphinp", imasterphinp_);
  }

  // only relevant for monolithic or semi-implicit evaluation of scatra-scatra interface layer
  // growth
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic or
      intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // extract relevant state vector of discrete scatra-scatra interface layer thicknesses
    const Teuchos::RCP<Epetra_Vector>& growth =
        intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic ? growthnp_
                                                                               : growthn_;

    // set state vector
    scatratimint_->Discretization()->SetState(2, "growth", growth);
  }

  return;
}


/*---------------------------------------------------------------------------*
 | compute time step size                                         fang 02/18 |
 *---------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeTimeStepSize(double& dt)
{
  // not implemented for standard scalar transport
  if (intlayergrowth_timestep_ > 0.)
    dserror(
        "Adaptive time stepping for scatra-scatra interface layer growth not implemented for "
        "standard scalar transport!");

  return;
}


/*-------------------------------------------------------------------------------*
 | perform initialization of scatra-scatra interface coupling         fang 10/14 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::InitMeshtying()
{
  // instantiate strategy for Newton-Raphson convergence check
  InitConvCheckStrategy();

  // extract boundary conditions for scatra-scatra interface layer growth
  std::vector<Teuchos::RCP<DRT::Condition>> conditions;
  scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth", conditions);

  // initialize scatra-scatra interface layer growth
  if (conditions.size())
  {
    // safety checks
    if (conditions.size() != 1)
      dserror(
          "Can't have more than one boundary condition for scatra-scatra interface layer growth at "
          "the moment!");
    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_none)
      dserror(
          "Invalid flag for evaluation of scatra-scatra interface coupling involving interface "
          "layer growth!");
    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic and
        scatratimint_->MethodName() != INPAR::SCATRA::timeint_one_step_theta)
      dserror(
          "Monolithic evaluation of scatra-scatra interface layer growth only implemented for "
          "one-step-theta time integration scheme at the moment!");
    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit and
        conditions[0]->GetInt("regularization type") !=
            INPAR::S2I::RegularizationType::regularization_none)
      dserror(
          "No regularization implemented for semi-implicit evaluation of scatra-scatra interface "
          "layer growth!");
    if (couplingtype_ != INPAR::S2I::coupling_matching_nodes)
      dserror(
          "Evaluation of scatra-scatra interface layer growth only implemented for conforming "
          "interface discretizations!");

    // provide scalar transport discretization with additional dofset for scatra-scatra interface
    // layer thickness
    const Teuchos::RCP<Epetra_IntVector> numdofpernode =
        Teuchos::rcp(new Epetra_IntVector(*scatratimint_->Discretization()->NodeColMap()));
    for (int inode = 0; inode < scatratimint_->Discretization()->NumMyColNodes(); ++inode)
      // add one degree of freedom for scatra-scatra interface layer growth to current node if
      // applicable
      if (scatratimint_->Discretization()->lColNode(inode)->GetCondition("S2ICouplingGrowth"))
        (*numdofpernode)[inode] = 1;
    Teuchos::RCP<DRT::DofSetInterface> dofset = Teuchos::rcp(
        new DRT::DofSetPredefinedDoFNumber(numdofpernode, Teuchos::null, Teuchos::null, true));
    if (scatratimint_->Discretization()->AddDofSet(dofset) != 2)
      dserror("Scalar transport discretization exhibits invalid number of dofsets!");

    // initialize linear solver for monolithic scatra-scatra interface coupling involving interface
    // layer growth
    if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
    {
      const int extendedsolver = DRT::Problem::Instance()
                                     ->ScalarTransportDynamicParams()
                                     .sublist("S2I COUPLING")
                                     .get<int>("INTLAYERGROWTH_LINEAR_SOLVER");
      if (extendedsolver < 1)
        dserror(
            "Invalid ID of linear solver for monolithic scatra-scatra interface coupling involving "
            "interface layer growth!");
      extendedsolver_ =
          Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(extendedsolver),
              scatratimint_->Discretization()->Comm(),
              DRT::Problem::Instance()->ErrorFile()->Handle()));
    }
  }  // initialize scatra-scatra interface layer growth

  // safety check
  else if (intlayergrowth_evaluation_ != INPAR::S2I::growth_evaluation_none)
    dserror(
        "Cannot evaluate scatra-scatra interface coupling involving interface layer growth without "
        "specifying a corresponding boundary condition!");

  // safety checks associated with adaptive time stepping for scatra-scatra interface layer growth
  if (intlayergrowth_timestep_ > 0.)
  {
    if (not DRT::INPUT::IntegralValue<bool>(
            *scatratimint_->ScatraParameterList(), "ADAPTIVE_TIMESTEPPING"))
      dserror(
          "Adaptive time stepping for scatra-scatra interface layer growth requires "
          "ADAPTIVE_TIMESTEPPING flag to be set!");
    if (not scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth"))
      dserror(
          "Adaptive time stepping for scatra-scatra interface layer growth requires corresponding "
          "boundary condition!");
    if (not(intlayergrowth_timestep_ < scatratimint_->Dt()))
      dserror(
          "Adaptive time stepping for scatra-scatra interface layer growth requires that the "
          "modified time step size is smaller than the original time step size!");
  }

  return;
}


/*----------------------------------------------------------------------------------*
 | build map extractors associated with blocks of global system matrix   fang 07/15 |
 *----------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::BuildBlockMapExtractors()
{
  if (matrixtype_ == INPAR::SCATRA::MatrixType::block_condition or
      matrixtype_ == INPAR::SCATRA::MatrixType::block_condition_dof)
  {
    // initialize reduced interface map extractors associated with blocks of global system matrix
    const int nblocks = scatratimint_->BlockMaps().NumMaps();
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps_slave(nblocks);
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps_master(nblocks);
    for (int iblock = 0; iblock < nblocks; ++iblock)
    {
      std::vector<Teuchos::RCP<const Epetra_Map>> maps(2);
      maps[0] = scatratimint_->BlockMaps().Map(iblock);
      maps[1] = not imortarredistribution_
                    ? interfacemaps_->Map(1)
                    : Teuchos::rcp_dynamic_cast<const Epetra_Map>(islavemap_);
      blockmaps_slave[iblock] = LINALG::MultiMapExtractor::IntersectMaps(maps);
      maps[1] = not imortarredistribution_
                    ? interfacemaps_->Map(2)
                    : Teuchos::rcp_dynamic_cast<const Epetra_Map>(imastermap_);
      blockmaps_master[iblock] = LINALG::MultiMapExtractor::IntersectMaps(maps);
    }
    blockmaps_slave_ =
        Teuchos::rcp(new LINALG::MultiMapExtractor(*interfacemaps_->Map(1), blockmaps_slave));
    blockmaps_slave_->CheckForValidMapExtractor();
    blockmaps_master_ =
        Teuchos::rcp(new LINALG::MultiMapExtractor(*interfacemaps_->Map(2), blockmaps_master));
    blockmaps_master_->CheckForValidMapExtractor();
  }
}  // SCATRA::MeshtyingStrategyS2I::BuildBlockMapExtractors

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquipExtendedSolverWithNullSpaceInfo() const
{
  // consider extended linear solver for scatra-scatra interface layer growth if applicable
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
  {
    // loop over blocks of scalar transport system matrix
    for (int iblock = 0; iblock < scatratimint_->BlockMaps().NumMaps(); ++iblock)
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
    iblockstr << scatratimint_->BlockMaps().NumMaps() + 1;

    // equip smoother for extra matrix block with null space associated with all degrees of freedom
    // for scatra-scatra interface layer growth
    Teuchos::ParameterList& mllist =
        extendedsolver_->Params().sublist("Inverse" + iblockstr.str()).sublist("MueLu Parameters");
    mllist.set("PDE equations", 1);
    mllist.set("null space: dimension", 1);
    mllist.set("null space: type", "pre-computed");
    mllist.set("null space: add default vectors", false);
    const Teuchos::RCP<std::vector<double>> ns =
        Teuchos::rcp(new std::vector<double>(scatratimint_->DofRowMap(2)->NumMyElements(), 1.));
    mllist.set<Teuchos::RCP<std::vector<double>>>("nullspace", ns);
    mllist.set("null space: vectors", &((*ns)[0]));
    mllist.set<bool>("ML validate parameter list", false);
  }
}  // SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces

/*------------------------------------------------------------------------------------*
 | solve linear system of equations for scatra-scatra interface coupling   fang 12/14 |
 *------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::Solve(const Teuchos::RCP<LINALG::Solver>& solver,  //!< solver
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Teuchos::RCP<Epetra_Vector>& increment,              //!< increment vector
    const Teuchos::RCP<Epetra_Vector>& residual,               //!< residual vector
    const Teuchos::RCP<Epetra_Vector>& phinp,                  //!< state vector at time n+1
    const int& iteration,  //!< number of current Newton-Raphson iteration
    const Teuchos::RCP<LINALG::KrylovProjector>& projector  //!< Krylov projector
    ) const
{
  switch (intlayergrowth_evaluation_)
  {
    // no or semi-implicit treatment of scatra-scatra interface layer growth
    case INPAR::S2I::growth_evaluation_none:
    case INPAR::S2I::growth_evaluation_semi_implicit:
    {
      switch (couplingtype_)
      {
        case INPAR::S2I::coupling_matching_nodes:
        case INPAR::S2I::coupling_mortar_standard:
        case INPAR::S2I::coupling_mortar_condensed_petrov:
        case INPAR::S2I::coupling_mortar_condensed_bubnov:
        case INPAR::S2I::coupling_nts_standard:
        {
          // equilibrate global system of equations if necessary
          EquilibrateSystem(systemmatrix, residual, scatratimint_->BlockMaps());

          // solve global system of equations
          solver->Solve(
              systemmatrix->EpetraOperator(), increment, residual, true, iteration == 1, projector);

          // unequilibrate global increment vector if necessary
          UnequilibrateIncrement(increment);

          break;
        }

        case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
        case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
        {
          // check scalar transport system matrix
          Teuchos::RCP<LINALG::SparseMatrix> sparsematrix =
              Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
          if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

          // assemble extended system matrix including rows and columns associated with Lagrange
          // multipliers
          LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> extendedsystemmatrix(
              *extendedmaps_, *extendedmaps_);
          extendedsystemmatrix.Assign(0, 0, LINALG::View, *sparsematrix);
          if (lmside_ == INPAR::S2I::side_slave)
          {
            extendedsystemmatrix.Matrix(0, 1).Add(*D_, true, 1., 0.);
            extendedsystemmatrix.Matrix(0, 1).Add(*M_, true, -1., 1.);
            extendedsystemmatrix.Matrix(1, 0).Add(
                *MORTAR::MatrixRowTransformGIDs(islavematrix_, extendedmaps_->Map(1)), false, 1.,
                0.);
          }
          else
          {
            extendedsystemmatrix.Matrix(0, 1).Add(*M_, true, -1., 0.);
            extendedsystemmatrix.Matrix(0, 1).Add(*D_, true, 1., 1.);
            extendedsystemmatrix.Matrix(1, 0).Add(
                *MORTAR::MatrixRowTransformGIDs(imastermatrix_, extendedmaps_->Map(1)), false, 1.,
                0.);
          }
          extendedsystemmatrix.Matrix(1, 1).Add(*E_, true, -1., 0.);
          extendedsystemmatrix.Complete();
          extendedsystemmatrix.Matrix(0, 1).ApplyDirichlet(
              *scatratimint_->DirichMaps()->CondMap(), false);

          Teuchos::RCP<Epetra_Vector> extendedresidual =
              LINALG::CreateVector(*extendedmaps_->FullMap());
          extendedmaps_->InsertVector(scatratimint_->Residual(), 0, extendedresidual);
          extendedmaps_->InsertVector(lmresidual_, 1, extendedresidual);

          Teuchos::RCP<Epetra_Vector> extendedincrement =
              LINALG::CreateVector(*extendedmaps_->FullMap());
          extendedmaps_->InsertVector(scatratimint_->Increment(), 0, extendedincrement);
          extendedmaps_->InsertVector(lmincrement_, 1, extendedincrement);

          // solve extended system of equations
          solver->Solve(extendedsystemmatrix.EpetraOperator(), extendedincrement, extendedresidual,
              true, iteration == 1, projector);

          // store solution
          extendedmaps_->ExtractVector(extendedincrement, 0, increment);
          extendedmaps_->ExtractVector(extendedincrement, 1, lmincrement_);

          // update Lagrange multipliers
          lm_->Update(1., *lmincrement_, 1.);

          break;
        }

        default:
        {
          dserror("Type of scatra-scatra interface coupling not recognized!");
          break;
        }
      }

      break;
    }

    // monolithic treatment of scatra-scatra interface layer growth
    case INPAR::S2I::growth_evaluation_monolithic:
    {
      switch (couplingtype_)
      {
        case INPAR::S2I::coupling_matching_nodes:
        {
          switch (matrixtype_)
          {
            case INPAR::SCATRA::MatrixType::sparse:
            {
              // check scalar transport system matrix
              const Teuchos::RCP<const LINALG::SparseMatrix> sparsematrix =
                  Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
              if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

              // assemble extended system matrix including rows and columns associated with
              // scatra-scatra interface layer thickness variables
              extendedsystemmatrix_->Assign(0, 0, LINALG::View, *sparsematrix);
              extendedsystemmatrix_->Assign(0, 1, LINALG::View,
                  *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(scatragrowthblock_));
              extendedsystemmatrix_->Assign(1, 0, LINALG::View,
                  *Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(growthscatrablock_));
              extendedsystemmatrix_->Assign(1, 1, LINALG::View, *growthgrowthblock_);

              break;
            }

            case INPAR::SCATRA::MatrixType::block_condition:
            case INPAR::SCATRA::MatrixType::block_condition_dof:
            {
              // check scalar transport system matrix
              const Teuchos::RCP<const LINALG::BlockSparseMatrixBase> blocksparsematrix =
                  Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
              if (blocksparsematrix == Teuchos::null)
                dserror("System matrix is not a block sparse matrix!");

              // extract number of matrix row or column blocks associated with scalar transport
              // field
              const unsigned nblockmaps = scatratimint_->BlockMaps().NumMaps();

              // construct extended system matrix by assigning matrix blocks
              for (unsigned iblock = 0; iblock < nblockmaps; ++iblock)
              {
                for (unsigned jblock = 0; jblock < nblockmaps; ++jblock)
                  extendedsystemmatrix_->Assign(
                      iblock, jblock, LINALG::View, blocksparsematrix->Matrix(iblock, jblock));
                extendedsystemmatrix_->Assign(iblock, nblockmaps, LINALG::View,
                    Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                        scatragrowthblock_)
                        ->Matrix(iblock, 0));
                extendedsystemmatrix_->Assign(nblockmaps, iblock, LINALG::View,
                    Teuchos::rcp_dynamic_cast<const LINALG::BlockSparseMatrixBase>(
                        growthscatrablock_)
                        ->Matrix(0, iblock));
              }
              extendedsystemmatrix_->Assign(
                  nblockmaps, nblockmaps, LINALG::View, *growthgrowthblock_);

              break;
            }

            default:
            {
              dserror(
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
          if (scatratimint_->FDCheckType() == INPAR::SCATRA::fdcheck_global_extended)
            FDCheck(*extendedsystemmatrix_, extendedresidual);

          // assemble extended increment vector
          Teuchos::RCP<Epetra_Vector> extendedincrement =
              Teuchos::rcp(new Epetra_Vector(*extendedmaps_->FullMap(), true));
          extendedmaps_->InsertVector(scatratimint_->Increment(), 0, extendedincrement);
          extendedmaps_->InsertVector(growthincrement_, 1, extendedincrement);

          // equilibrate global system of equations if necessary
          EquilibrateSystem(extendedsystemmatrix_, extendedresidual, *extendedblockmaps_);

          // solve extended system of equations
          extendedsolver_->Solve(extendedsystemmatrix_->EpetraOperator(), extendedincrement,
              extendedresidual, true, iteration == 1, projector);

          // unequilibrate global increment vector if necessary
          UnequilibrateIncrement(extendedincrement);

          // store solution
          extendedmaps_->ExtractVector(extendedincrement, 0, increment);
          extendedmaps_->ExtractVector(extendedincrement, 1, growthincrement_);

          // update state vector of discrete scatra-scatra interface layer thicknesses
          growthnp_->Update(1., *growthincrement_, 1.);

          break;
        }

        default:
        {
          dserror("Type of scatra-scatra interface coupling not recognized!");
          break;
        }
      }  // switch(couplingtype_)

      break;
    }

    default:
    {
      dserror(
          "Unknown evaluation method for scatra-scatra interface coupling involving interface "
          "layer growth!");
      break;
    }
  }

  return;
}  // SCATRA::MeshtyingStrategyS2I::Solve


/*-------------------------------------------------------------------------*
 | return linear solver for global system of linear equations   fang 01/18 |
 *-------------------------------------------------------------------------*/
const LINALG::Solver& SCATRA::MeshtyingStrategyS2I::Solver() const
{
  const LINALG::Solver* solver(nullptr);

  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_monolithic)
  {
    if (extendedsolver_ == Teuchos::null) dserror("Invalid linear solver!");
    solver = extendedsolver_.get();
  }

  else
  {
    if (scatratimint_->Solver() == Teuchos::null) dserror("Invalid linear solver!");
    solver = scatratimint_->Solver().get();
  }

  return *solver;
}  // SCATRA::MeshtyingStrategyS2I::Solver()


/*----------------------------------------------------------------------*
 | equilibrate global system of equations if necessary       fang 05/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateSystem(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& residual, const LINALG::MultiMapExtractor& blockmaps) const
{
  if (scatratimint_->EquilibrationMethod() != INPAR::SCATRA::EquilibrationMethod::none)
  {
    // perform equilibration depending on type of global system matrix
    switch (matrixtype_)
    {
      case INPAR::SCATRA::MatrixType::sparse:
      {
        // check matrix
        Teuchos::RCP<LINALG::SparseMatrix> sparsematrix =
            Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix);
        if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

        // perform row equilibration
        if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rows_maindiag or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
        {
          // compute inverse row sums of global system matrix
          ComputeInvRowSums(*sparsematrix, invrowsums_);

          // perform row equilibration of global system matrix
          EquilibrateMatrixRows(*sparsematrix, invrowsums_);
        }

        // perform column equilibration
        if (scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::columns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::columns_maindiag or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
        {
          // compute inverse column sums of global system matrix
          ComputeInvColSums(*sparsematrix, invcolsums_);

          // perform column equilibration of global system matrix
          EquilibrateMatrixColumns(*sparsematrix, invcolsums_);
        }

        break;
      }

      case INPAR::SCATRA::MatrixType::block_condition:
      case INPAR::SCATRA::MatrixType::block_condition_dof:
      case INPAR::SCATRA::MatrixType::block_geometry:
      {
        // check matrix
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocksparsematrix =
            Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
        if (blocksparsematrix == Teuchos::null)
          dserror("System matrix is not a block sparse matrix!");

        // perform row equilibration
        if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rows_maindiag or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
        {
          for (int i = 0; i < blocksparsematrix->Rows(); ++i)
          {
            // initialize vector for inverse row sums
            Teuchos::RCP<Epetra_Vector> invrowsums(
                Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).RowMap())));

            // compute inverse row sums of current main diagonal matrix block
            if (scatratimint_->EquilibrationMethod() ==
                    INPAR::SCATRA::EquilibrationMethod::rows_maindiag or
                scatratimint_->EquilibrationMethod() ==
                    INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
              ComputeInvRowSums(blocksparsematrix->Matrix(i, i), invrowsums);

            // compute inverse row sums of current row block of global system matrix
            else
            {
              // loop over all column blocks of global system matrix
              for (int j = 0; j < blocksparsematrix->Cols(); ++j)
              {
                // extract current block of global system matrix
                const LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

                // loop over all rows of current matrix block
                for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
                {
                  // determine length of current matrix row
                  const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

                  if (length > 0)
                  {
                    // extract current matrix row from matrix block
                    int numentries(0);
                    std::vector<double> values(length, 0.);
                    if (matrix.EpetraMatrix()->ExtractMyRowCopy(
                            irow, length, numentries, &values[0]))
                      dserror(
                          "Cannot extract matrix row with local ID %d from matrix block!", irow);

                    // compute and store current row sum
                    double rowsum(0.);
                    for (int ientry = 0; ientry < numentries; ++ientry)
                      rowsum += std::abs(values[ientry]);
                    (*invrowsums)[irow] += rowsum;
                  }
                }
              }

              // invert row sums
              if (invrowsums->Reciprocal(*invrowsums)) dserror("Vector could not be inverted!");

              // take square root of inverse row sums if matrix is scaled from left and right
              if (scatratimint_->EquilibrationMethod() ==
                      INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
                  scatratimint_->EquilibrationMethod() ==
                      INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
                for (int i = 0; i < invrowsums->MyLength(); ++i)
                  (*invrowsums)[i] = sqrt((*invrowsums)[i]);
            }

            // perform row equilibration of matrix blocks in current row block of global system
            // matrix
            for (int j = 0; j < blocksparsematrix->Cols(); ++j)
              EquilibrateMatrixRows(blocksparsematrix->Matrix(i, j), invrowsums);

            // insert inverse row sums of current main diagonal matrix block into global vector
            blockmaps.InsertVector(invrowsums, i, invrowsums_);
          }
        }

        // perform column equilibration
        if (scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::columns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::columns_maindiag or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
            scatratimint_->EquilibrationMethod() ==
                INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
        {
          for (int j = 0; j < blocksparsematrix->Cols(); ++j)
          {
            // initialize vector for inverse column sums
            Teuchos::RCP<Epetra_Vector> invcolsums(
                Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(j, j).DomainMap())));

            // compute inverse column sums of current main diagonal matrix block
            if (scatratimint_->EquilibrationMethod() ==
                    INPAR::SCATRA::EquilibrationMethod::columns_maindiag or
                scatratimint_->EquilibrationMethod() ==
                    INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
              ComputeInvColSums(blocksparsematrix->Matrix(j, j), invcolsums);

            // compute inverse column sums of current column block of global system matrix
            else
            {
              // loop over all row blocks of global system matrix
              for (int i = 0; i < blocksparsematrix->Rows(); ++i)
              {
                // extract current block of global system matrix
                const LINALG::SparseMatrix& matrix = blocksparsematrix->Matrix(i, j);

                // loop over all rows of current matrix block
                for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
                {
                  // determine length of current matrix row
                  const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

                  if (length > 0)
                  {
                    // extract current matrix row from matrix block
                    int numentries(0);
                    std::vector<double> values(length, 0.);
                    std::vector<int> indices(length, 0);
                    if (matrix.EpetraMatrix()->ExtractMyRowCopy(
                            irow, length, numentries, &values[0], &indices[0]))
                      dserror(
                          "Cannot extract matrix row with local ID %d from matrix block!", irow);

                    // add entries of current matrix row to column sums
                    for (int ientry = 0; ientry < numentries; ++ientry)
                      invcolsums->SumIntoGlobalValue(
                          matrix.ColMap().GID(indices[ientry]), 0, std::abs(values[ientry]));
                  }
                }
              }

              // invert column sums
              if (invcolsums->Reciprocal(*invcolsums)) dserror("Vector could not be inverted!");

              // take square root of inverse column sums if matrix is scaled from left and right
              if (scatratimint_->EquilibrationMethod() ==
                      INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
                  scatratimint_->EquilibrationMethod() ==
                      INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
                for (int i = 0; i < invcolsums->MyLength(); ++i)
                  (*invcolsums)[i] = sqrt((*invcolsums)[i]);
            }

            // perform column equilibration of matrix blocks in current column block of global
            // system matrix
            for (int i = 0; i < blocksparsematrix->Rows(); ++i)
              EquilibrateMatrixColumns(blocksparsematrix->Matrix(i, j), invcolsums);

            // insert inverse column sums of current main diagonal matrix block into global vector
            blockmaps.InsertVector(invcolsums, j, invcolsums_);
          }
        }

        break;
      }

      default:
      {
        dserror(
            "Equilibration of global system of equations for scatra-scatra interface coupling is "
            "not implemented for chosen type of global system matrix!");
        break;
      }
    }

    // perform equilibration of global residual vector
    if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_full or
        scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::rows_maindiag or
        scatratimint_->EquilibrationMethod() ==
            INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
        scatratimint_->EquilibrationMethod() ==
            INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
      if (residual->Multiply(1., *invrowsums_, *residual, 0.))
        dserror("Equilibration of global residual vector failed!");
  }

}  // SCATRA::MeshtyingStrategyS2I::EquilibrateSystem


/*-------------------------------------------------------------------------------------------------------------------------------------*
 | finite difference check for extended system matrix involving scatra-scatra interface layer growth
 (for debugging only)   fang 01/17 |
 *-------------------------------------------------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::FDCheck(
    const LINALG::BlockSparseMatrixBase& extendedsystemmatrix,  //!< global system matrix
    const Teuchos::RCP<Epetra_Vector>& extendedresidual         //!< global residual vector
    ) const
{
  // initial screen output
  if (scatratimint_->Discretization()->Comm().MyPID() == 0)
    std::cout << std::endl
              << "FINITE DIFFERENCE CHECK FOR EXTENDED SYSTEM MATRIX INVOLVING SCATRA-SCATRA "
                 "INTERFACE LAYER GROWTH"
              << std::endl;

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
      *LINALG::SparseMatrix(*extendedsystemmatrix.Merge()).EpetraMatrix();
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
    scatratimint_->Discretization()->Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill global state vector with original state variables
    statenp.Update(1., statenp_original, 0.);

    // impose perturbation
    if (statenp.Map().MyGID(colgid))
      if (statenp.SumIntoGlobalValue(colgid, 0, fdcheckeps))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    scatratimint_->Phinp()->Update(1., *extendedmaps_->ExtractVector(statenp, 0), 0.);
    growthnp_->Update(1., *extendedmaps_->ExtractVector(statenp, 1), 0.);

    // calculate global right-hand side contributions based on perturbed state
    scatratimint_->AssembleMatAndRHS();

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
      if (rowgid < 0) dserror("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original.NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original.ExtractMyRowCopy(rowlid, length, numentries, &values[0], &indices[0]);
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
  scatratimint_->Discretization()->Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  scatratimint_->Discretization()->Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  scatratimint_->Discretization()->Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (scatratimint_->Discretization()->Comm().MyPID() == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      dserror(
          "Finite difference check failed for extended system matrix involving scatra-scatra "
          "interface layer growth!");
    }
    else
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatratimint_->Phinp()->Update(1., *extendedmaps_->ExtractVector(statenp_original, 0), 0.);
  growthnp_->Update(1., *extendedmaps_->ExtractVector(statenp_original, 1), 0.);

  // recompute system matrix and right-hand side vector based on original state variables
  scatratimint_->AssembleMatAndRHS();

  return;
}


/*----------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix row entries   fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeInvRowSums(
    const LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< inverse sums of absolute values of row entries in matrix
    ) const
{
  // compute inverse row sums of matrix
  if (matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");

  // take square root of inverse row sums if matrix is scaled from left and right
  if (scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invrowsums->MyLength(); ++i) (*invrowsums)[i] = sqrt((*invrowsums)[i]);
}  // SCATRA::MeshtyingStrategyS2I::ComputeInvRowSums


/*-------------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix column entries   fang 06/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::ComputeInvColSums(
    const LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invcolsums  //!< inverse sums of absolute values of column entries in matrix
    ) const
{
  // compute inverse column sums of matrix
  if (matrix.EpetraMatrix()->InvColSums(*invcolsums))
    dserror("Inverse column sums of matrix could not be successfully computed!");

  // take square root of inverse column sums if matrix is scaled from left and right
  if (scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
    for (int i = 0; i < invcolsums->MyLength(); ++i) (*invcolsums)[i] = sqrt((*invcolsums)[i]);
}  // SCATRA::MeshtyingStrategyS2I::ComputeInvColSums


/*----------------------------------------------------------------------*
 | equilibrate matrix rows                                   fang 06/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixRows(LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< sums of absolute values of row entries in matrix
    ) const
{
  if (matrix.LeftScale(*invrowsums)) dserror("Row equilibration of matrix failed!");

  return;
}  // SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixRows


/*----------------------------------------------------------------------*
 | equilibrate matrix columns                                fang 06/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixColumns(
    LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invcolsums  //!< sums of absolute values of column entries in matrix
    ) const
{
  if (matrix.RightScale(*invcolsums)) dserror("Column equilibration of matrix failed!");

  return;
}  // SCATRA::MeshtyingStrategyS2I::EquilibrateMatrixColumns


/*----------------------------------------------------------------------*
 | unequilibrate global increment vector if necessary        fang 05/15 |
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2I::UnequilibrateIncrement(
    const Teuchos::RCP<Epetra_Vector>& increment  //!< increment vector
    ) const
{
  // unequilibrate global increment vector if necessary
  if (scatratimint_->EquilibrationMethod() == INPAR::SCATRA::EquilibrationMethod::columns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::columns_maindiag or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_full or
      scatratimint_->EquilibrationMethod() ==
          INPAR::SCATRA::EquilibrationMethod::rowsandcolumns_maindiag)
  {
    if (increment->Multiply(1., *invcolsums_, *increment, 0.))
      dserror("Unequilibration of global increment vector failed!");
  }
}  // SCATRA::MeshtyingStrategyS2I::UnequilibrateIncrement


/*--------------------------------------------------------------------------------------*
 | protected constructor for singletons                                      fang 01/16 |
 *--------------------------------------------------------------------------------------*/
SCATRA::MortarCellInterface::MortarCellInterface(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    : lmside_(lmside),
      couplingtype_(couplingtype),
      numdofpernode_slave_(numdofpernode_slave),
      numdofpernode_master_(numdofpernode_master)
{
  return;
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalc<distypeS, distypeM>* SCATRA::MortarCellCalc<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname,       //!< name of mortar discretization
    const MortarCellCalc* delete_me   //!< pointer to instance to be deleted
)
{
  // static map assigning mortar discretization names to class instances
  static std::map<std::string, MortarCellCalc<distypeS, distypeM>*> instances;

  // create new instance or return existing one
  if (!delete_me)
  {
    // create new instance if not yet available
    if (instances.find(disname) == instances.end())
      instances[disname] = new MortarCellCalc<distypeS, distypeM>(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master);
  }

  // delete existing instance
  else
  {
    // loop over all existing instances
    for (typename std::map<std::string, MortarCellCalc<distypeS, distypeM>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
    {
      // check whether current instance should be deleted
      if (i->second == delete_me)
      {
        // delete current instance
        delete i->second;

        // remove deleted instance from map
        instances.erase(i);

        // return null pointer
        return nullptr;
      }
    }

    // catch internal error
    dserror("Instance to be deleted couldn't be found in static map!");
  }

  // return existing instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::Done()
{
  // delete singleton
  Instance(INPAR::S2I::coupling_undefined, INPAR::S2I::side_undefined, 0, 0, "", this);

  return;
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization
 types   fang 01/16 |
 *--------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::Evaluate(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
    Epetra_SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
    Epetra_SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
    Epetra_SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
    Epetra_SerialDenseVector& cellvector1,   //!< cell vector 1
    Epetra_SerialDenseVector& cellvector2    //!< cell vector 2
)
{
  // extract and evaluate action
  switch (DRT::INPUT::get<INPAR::S2I::EvaluationActions>(params, "action"))
  {
    case INPAR::S2I::evaluate_mortar_matrices:
    {
      // evaluate mortar matrices
      EvaluateMortarMatrices(
          cell, slaveelement, masterelement, cellmatrix1, cellmatrix2, cellmatrix3);

      break;
    }

    case INPAR::S2I::evaluate_condition:
    {
      // evaluate and assemble interface linearizations and residuals
      EvaluateCondition(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);

      break;
    }

    default:
    {
      dserror("Unknown action for mortar cell evaluation!");
      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate single slave-side node for node-to-segment coupling   fang 08/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateNTS(
    const DRT::Discretization& idiscret,  //!< interface discretization
    const MORTAR::MortarNode& slavenode,  //!< slave-side node
    const double& lumpedarea,  //!< lumped interface area fraction associated with slave-side node
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& ntsmatrix1,    //!< node-to-segment matrix 1
    Epetra_SerialDenseMatrix& ntsmatrix2,    //!< node-to-segment matrix 2
    Epetra_SerialDenseMatrix& ntsmatrix3,    //!< node-to-segment matrix 3
    Epetra_SerialDenseMatrix& ntsmatrix4,    //!< node-to-segment matrix 4
    Epetra_SerialDenseVector& ntsvector1,    //!< node-to-segment vector 1
    Epetra_SerialDenseVector& ntsvector2     //!< node-to-segment vector 2
)
{
  // extract and evaluate action
  switch (DRT::INPUT::get<INPAR::S2I::EvaluationActions>(params, "action"))
  {
    case INPAR::S2I::evaluate_condition_nts:
    {
      // extract condition from parameter list
      DRT::Condition* condition = params.get<DRT::Condition*>("condition");
      if (condition == nullptr)
        dserror("Cannot access scatra-scatra interface coupling condition!");

      // extract nodal state variables associated with slave and master elements
      ExtractNodeValues(idiscret, la_slave, la_master);

      // evaluate and assemble interface linearizations and residuals
      EvaluateConditionNTS(*condition, slavenode, lumpedarea, slaveelement, masterelement,
          ephinp_slave_, ephinp_master_, ntsmatrix1, ntsmatrix2, ntsmatrix3, ntsmatrix4, ntsvector1,
          ntsvector2);

      break;
    }

    default:
    {
      dserror("Unknown action for evaluation of node-to-segment coupling!");
      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate single mortar element                                 fang 08/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateMortarElement(
    const DRT::Discretization& idiscret,   //!< interface discretization
    MORTAR::MortarElement& element,        //!< mortar element
    DRT::Element::LocationArray& la,       //!< location array
    const Teuchos::ParameterList& params,  //!< parameter list
    Epetra_SerialDenseMatrix& elematrix1,  //!< element matrix 1
    Epetra_SerialDenseMatrix& elematrix2,  //!< element matrix 2
    Epetra_SerialDenseMatrix& elematrix3,  //!< element matrix 3
    Epetra_SerialDenseMatrix& elematrix4,  //!< element matrix 4
    Epetra_SerialDenseVector& elevector1,  //!< element vector 1
    Epetra_SerialDenseVector& elevector2   //!< element vector 2
)
{
  // extract and evaluate action
  switch (DRT::INPUT::get<INPAR::S2I::EvaluationActions>(params, "action"))
  {
    case INPAR::S2I::evaluate_nodal_area_fractions:
    {
      // evaluate and assemble lumped interface area fractions associated with element nodes
      EvaluateNodalAreaFractions(element, elevector1);

      break;
    }

    default:
    {
      dserror("Unknown action for evaluation of mortar element!");
      break;
    }
  }

  return;
}


/*--------------------------------------------------------------------------------------*
 | protected constructor for singletons                                      fang 01/16 |
 *--------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalc<distypeS, distypeM>::MortarCellCalc(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    : MortarCellInterface(couplingtype, lmside, numdofpernode_slave, numdofpernode_master),
      scatraparamsboundary_(DRT::ELEMENTS::ScaTraEleParameterBoundary::Instance("scatra")),
      ephinp_slave_(numdofpernode_slave, LINALG::Matrix<nen_slave_, 1>(true)),
      ephinp_master_(numdofpernode_master, LINALG::Matrix<nen_master_, 1>(true)),
      funct_slave_(true),
      funct_master_(true),
      shape_lm_slave_(true),
      shape_lm_master_(true),
      test_lm_slave_(true),
      test_lm_master_(true)
{
  // safety check
  if (nsd_slave_ != 2 or nsd_master_ != 2)
    dserror(
        "Scatra-scatra interface coupling with non-matching interface discretization currently "
        "only implemented for two-dimensional interface manifolds!");

  return;
}


/*------------------------------------------------------------------------------------*
 | extract nodal state variables associated with mortar integration cell   fang 01/17 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::ExtractNodeValues(
    const DRT::Discretization& idiscret,    //!< interface discretization
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
)
{
  // extract nodal state variables associated with mortar integration cell
  ExtractNodeValues(ephinp_slave_, ephinp_master_, idiscret, la_slave, la_master);

  return;
}


/*--------------------------------------------------------------------------*
 | extract nodal state variables associated with slave element   fang 01/17 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::ExtractNodeValues(
    LINALG::Matrix<nen_slave_, 1>& estate_slave,  //!< state variables at slave-side nodes
    const DRT::Discretization& idiscret,          //!< interface discretization
    DRT::Element::LocationArray& la_slave,        //!< slave-side location array
    const std::string& statename,                 //!< name of relevant state
    const int& nds                                //!< number of relevant dofset
    ) const
{
  // extract interface state vector from interface discretization
  const Teuchos::RCP<const Epetra_Vector> state = idiscret.GetState(nds, statename);
  if (state == Teuchos::null)
    dserror("Cannot extract state vector \"" + statename + "\" from interface discretization!");

  // extract nodal state variables associated with slave element
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_slave_, 1>>(
      *state, estate_slave, la_slave[nds].lm_);

  return;
}


/*--------------------------------------------------------------------------------------*
 | extract nodal state variables associated with slave and master elements   fang 01/16 |
 *--------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::ExtractNodeValues(
    std::vector<LINALG::Matrix<nen_slave_, 1>>&
        estate_slave,  //!< state variables at slave-side nodes
    std::vector<LINALG::Matrix<nen_master_, 1>>&
        estate_master,                       //!< state variables at master-side nodes
    const DRT::Discretization& idiscret,     //!< interface discretization
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const std::string& statename,            //!< name of relevant state
    const int& nds                           //!< number of relevant dofset
    ) const
{
  // extract interface state vector from interface discretization
  const Teuchos::RCP<const Epetra_Vector> state = idiscret.GetState(nds, statename);
  if (state == Teuchos::null)
    dserror("Cannot extract state vector \"" + statename + "\" from interface discretization!");

  // extract nodal state variables associated with slave and master elements
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_slave_, 1>>(
      *state, estate_slave, la_slave[nds].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_master_, 1>>(
      *state, estate_master, la_master[nds].lm_);

  return;
}


/*--------------------------------------------------------------------------------------------------------------------------*
 | evaluate slave-side and master-side shape functions and domain integration factor at cell
 integration point   fang 01/16 |
 *--------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
double SCATRA::MortarCellCalc<distypeS, distypeM>::EvalShapeFuncAndDomIntFacAtIntPoint(
    MORTAR::MortarElement& slaveelement,                           //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,                          //!< master-side mortar element
    MORTAR::IntCell& cell,                                         //!< mortar integration cell
    const DRT::UTILS::IntPointsAndWeights<nsd_slave_>& intpoints,  //!< quadrature rule
    const int iquad                                                //!< ID of integration point
)
{
  // reference coordinates of integration point
  double coordinates_ref[nsd_slave_] = {};
  for (int idim = 0; idim < nsd_slave_; ++idim)
    coordinates_ref[idim] = intpoints.IP().qxg[iquad][idim];

  // global coordinates of integration point
  double coordinates_global[nsd_slave_ + 1] = {};
  cell.LocalToGlobal(coordinates_ref, coordinates_global, 0);

  // project integration point onto slave and master elements
  double coordinates_slave[nsd_slave_] = {};
  double coordinates_master[nsd_master_] = {};
  double dummy(0.);
  MORTAR::MortarProjector::Impl(slaveelement)
      ->ProjectGaussPointAuxn3D(
          coordinates_global, cell.Auxn(), slaveelement, coordinates_slave, dummy);
  MORTAR::MortarProjector::Impl(masterelement)
      ->ProjectGaussPointAuxn3D(
          coordinates_global, cell.Auxn(), masterelement, coordinates_master, dummy);

  // evaluate shape functions at current integration point on slave and master elements
  VOLMORTAR::UTILS::shape_function<distypeS>(funct_slave_, coordinates_slave);
  VOLMORTAR::UTILS::shape_function<distypeM>(funct_master_, coordinates_master);
  switch (couplingtype_)
  {
    case INPAR::S2I::coupling_mortar_standard:
    {
      // there actually aren't any Lagrange multipliers, but we still need to set pseudo Lagrange
      // multiplier test functions equal to the standard shape and test functions for correct
      // evaluation of the scatra-scatra interface coupling conditions
      test_lm_slave_ = funct_slave_;
      test_lm_master_ = funct_master_;

      break;
    }

    case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
    case INPAR::S2I::coupling_mortar_condensed_petrov:
    {
      // dual Lagrange multiplier shape functions combined with standard Lagrange multiplier test
      // functions
      if (lmside_ == INPAR::S2I::side_slave)
      {
        VOLMORTAR::UTILS::dual_shape_function<distypeS>(
            shape_lm_slave_, coordinates_slave, slaveelement);
        test_lm_slave_ = funct_slave_;
      }
      else
      {
        VOLMORTAR::UTILS::dual_shape_function<distypeM>(
            shape_lm_master_, coordinates_master, masterelement);
        test_lm_master_ = funct_master_;
      }

      break;
    }

    case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
    case INPAR::S2I::coupling_mortar_condensed_bubnov:
    {
      // dual Lagrange multiplier shape functions combined with dual Lagrange multiplier test
      // functions
      if (lmside_ == INPAR::S2I::side_slave)
      {
        VOLMORTAR::UTILS::dual_shape_function<distypeS>(
            shape_lm_slave_, coordinates_slave, slaveelement);
        test_lm_slave_ = shape_lm_slave_;
      }
      else
      {
        VOLMORTAR::UTILS::dual_shape_function<distypeM>(
            shape_lm_master_, coordinates_master, masterelement);
        test_lm_master_ = shape_lm_master_;
      }

      break;
    }

    default:
    {
      dserror("Not yet implemented!");
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


/*-------------------------------------------------------------------------------------------------------------*
 | evaluate slave-side shape functions and domain integration factor at element integration point
 fang 08/16 |
 *-------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
double SCATRA::MortarCellCalc<distypeS, distypeM>::EvalShapeFuncAndDomIntFacAtIntPoint(
    MORTAR::MortarElement& element,                                //!< mortar element
    const DRT::UTILS::IntPointsAndWeights<nsd_slave_>& intpoints,  //!< quadrature rule
    const int iquad                                                //!< ID of integration point
)
{
  // extract global coordinates of element nodes
  LINALG::Matrix<nsd_slave_ + 1, nen_slave_> coordinates_nodes;
  GEO::fillInitialPositionArray<distypeS, nsd_slave_ + 1,
      LINALG::Matrix<nsd_slave_ + 1, nen_slave_>>(&element, coordinates_nodes);

  // extract reference coordinates of integration point
  LINALG::Matrix<nsd_slave_, 1> coordinates_ref(intpoints.IP().qxg[iquad]);

  // evaluate slave-side shape functions and their first derivatives at integration point
  LINALG::Matrix<nsd_slave_, nen_slave_> deriv_slave;
  DRT::UTILS::shape_function<distypeS>(coordinates_ref, funct_slave_);
  DRT::UTILS::shape_function_deriv1<distypeS>(coordinates_ref, deriv_slave);

  // evaluate transposed Jacobian matrix at integration point
  LINALG::Matrix<nsd_slave_, nsd_slave_ + 1> jacobian;
  jacobian.MultiplyNT(deriv_slave, coordinates_nodes);

  // evaluate metric tensor at integration point
  LINALG::Matrix<nsd_slave_, nsd_slave_> metrictensor;
  metrictensor.MultiplyNT(jacobian, jacobian);

  // return domain integration factor, i.e., Jacobian determinant times integration weight, at
  // integration point
  return sqrt(metrictensor.Determinant()) * intpoints.IP().qwgt[iquad];
}


/*------------------------------------------------------------------------------------------*
 | evaluate shape functions at position of slave-side node                       fang 08/16 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvalShapeFuncAtSlaveNode(
    const MORTAR::MortarNode& slavenode,  //!< slave-side node
    MORTAR::MortarElement& slaveelement,  //!< slave-side element
    MORTAR::MortarElement& masterelement  //!< master-side element
)
{
  // safety check
  if (couplingtype_ != INPAR::S2I::coupling_nts_standard)
    dserror("This function should only be called when evaluating node-to-segment coupling!");

  // extract global ID of slave-side node
  const int& nodeid = slavenode.Id();

  // find out index of slave-side node w.r.t. slave-side element
  int index(-1);
  for (int inode = 0; inode < slaveelement.NumNode(); ++inode)
    if (nodeid == slaveelement.Nodes()[inode]->Id())
    {
      index = inode;
      break;
    }
  if (index == -1) dserror("Couldn't find out index of slave-side node w.r.t. slave-side element!");

  // set slave-side shape function array according to node position
  funct_slave_.Clear();
  funct_slave_(index) = 1.;

  // project slave-side node onto master-side element
  double coordinates_master[2] = {};
  double dummy(0.);
  MORTAR::MortarProjector::Impl(masterelement)
      ->ProjectGaussPointAuxn3D(
          slavenode.X(), slavenode.MoData().n(), masterelement, coordinates_master, dummy);

  // evaluate master-side shape functions at projected node on master-side element
  VOLMORTAR::UTILS::shape_function<distypeM>(funct_master_, coordinates_master);

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate mortar matrices                                       fang 01/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateMortarMatrices(
    MORTAR::IntCell& cell,                 //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,   //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,  //!< master-side mortar element
    Epetra_SerialDenseMatrix& D,           //!< mortar matrix D
    Epetra_SerialDenseMatrix& M,           //!< mortar matrix M
    Epetra_SerialDenseMatrix& E            //!< mortar matrix E
)
{
  // safety check
  if (numdofpernode_slave_ != numdofpernode_master_)
    dserror("Must have same number of degrees of freedom per node on slave and master sides!");

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac =
        EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, masterelement, cell, intpoints, iquad);

    if (lmside_ == INPAR::S2I::side_slave)
    {
      // loop over all degrees of freedom per node
      for (int k = 0; k < numdofpernode_slave_; ++k)
      {
        for (int vi = 0; vi < nen_slave_; ++vi)
        {
          const int row_slave = vi * numdofpernode_slave_ + k;

          switch (couplingtype_)
          {
            case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
            case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
            case INPAR::S2I::coupling_mortar_condensed_petrov:
            case INPAR::S2I::coupling_mortar_condensed_bubnov:
            {
              D(row_slave, row_slave) += shape_lm_slave_(vi) * fac;

              if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov or
                  couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
                for (int ui = 0; ui < nen_slave_; ++ui)
                  E(row_slave, ui * numdofpernode_slave_ + k) +=
                      shape_lm_slave_(vi) * test_lm_slave_(ui) * fac;

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
            case INPAR::S2I::coupling_mortar_saddlepoint_petrov:
            case INPAR::S2I::coupling_mortar_saddlepoint_bubnov:
            case INPAR::S2I::coupling_mortar_condensed_petrov:
            case INPAR::S2I::coupling_mortar_condensed_bubnov:
            {
              D(row_master, row_master) += shape_lm_master_(vi) * fac;

              if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov or
                  couplingtype_ == INPAR::S2I::coupling_mortar_condensed_bubnov)
                for (int ui = 0; ui < nen_master_; ++ui)
                  E(row_master, ui * numdofpernode_master_ + k) +=
                      shape_lm_master_(vi) * test_lm_master_(ui) * fac;

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

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateCondition(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&
        k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
    Epetra_SerialDenseVector& r_m   //!< master-side residual vector
)
{
  // extract nodal state variables associated with slave and master elements
  ExtractNodeValues(idiscret, la_slave, la_master);

  // safety check
  if (numdofpernode_slave_ != 1 or numdofpernode_master_ != 1)
    dserror(
        "Invalid number of degrees of freedom per node! Code should theoretically work for more "
        "than one degree of freedom per node, but not yet tested!");

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  const int kineticmodel = scatraparamsboundary_->KineticModel();
  const std::vector<double>* permeabilities = scatraparamsboundary_->Permeabilities();

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac =
        EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, masterelement, cell, intpoints, iquad);

    // overall integration factors
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * fac;
    const double timefacrhsfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalc<distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<
        distypeM>(ephinp_slave_, ephinp_master_, funct_slave_, funct_master_, test_lm_slave_,
        test_lm_master_, numdofpernode_slave_, kineticmodel, permeabilities, timefacfac,
        timefacrhsfac, k_ss, k_sm, k_ms, k_mm, r_s, r_m);
  }

  return;
}


/*--------------------------------------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals for node-to-segment coupling   fang
 08/16 |
 *--------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateConditionNTS(
    DRT::Condition& condition,            //!< scatra-scatra interface coupling condition
    const MORTAR::MortarNode& slavenode,  //!< slave-side node
    const double& lumpedarea,  //!< lumped interface area fraction associated with slave-side node
    MORTAR::MortarElement& slaveelement,   //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,  //!< master-side mortar element
    const std::vector<LINALG::Matrix<nen_slave_, 1>>&
        ephinp_slave,  //!< state variables at slave-side nodes
    const std::vector<LINALG::Matrix<nen_master_, 1>>&
        ephinp_master,  //!< state variables at master-side nodes
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&
        k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
    Epetra_SerialDenseVector& r_m   //!< master-side residual vector
)
{
  // safety check
  if (numdofpernode_slave_ != 1 or numdofpernode_master_ != 1)
    dserror(
        "Invalid number of degrees of freedom per node! Code should theoretically work for more "
        "than one degree of freedom per node, but not yet tested!");

  // evaluate shape functions at position of slave-side node
  EvalShapeFuncAtSlaveNode(slavenode, slaveelement, masterelement);

  const int kineticmodel = scatraparamsboundary_->KineticModel();
  const std::vector<double>* permeabilities = scatraparamsboundary_->Permeabilities();

  // overall integration factors
  const double timefacfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * lumpedarea;
  const double timefacrhsfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * lumpedarea;
  if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

  DRT::ELEMENTS::ScaTraEleBoundaryCalc<distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<
      distypeM>(ephinp_slave, ephinp_master, funct_slave_, funct_master_, funct_slave_,
      funct_master_, numdofpernode_slave_, kineticmodel, permeabilities, timefacfac, timefacrhsfac,
      k_ss, k_sm, k_ms, k_mm, r_s, r_m);

  return;
}


/*-------------------------------------------------------------------------------------------------------------*
 | evaluate and assemble lumped interface area fractions associated with slave-side element nodes
 fang 08/16 |
 *-------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalc<distypeS, distypeM>::EvaluateNodalAreaFractions(
    MORTAR::MortarElement& slaveelement,  //!< slave-side mortar element
    Epetra_SerialDenseVector&
        areafractions  //!< lumped interface area fractions associated with slave-side element nodes
)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_slave_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distypeS>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, intpoints, iquad);

    // compute integrals of shape functions to obtain lumped interface area fractions associated
    // with element nodes
    for (int inode = 0; inode < nen_slave_; ++inode)
      areafractions[inode * numdofpernode_slave_] += funct_slave_(inode) * fac;
  }  // loop over integration points

  return;
}


/*---------------------------------------------------------------------------*
 | constructor                                                    fang 05/16 |
 *---------------------------------------------------------------------------*/
SCATRA::MortarCellAssemblyStrategy::MortarCellAssemblyStrategy(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix1,  //!< system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_rows,  //!< interface side associated with rows of system matrix 1
    const INPAR::S2I::InterfaceSides
        matrix1_side_cols,  //!< interface side associated with columns of system matrix 1
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix2,  //!< system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_rows,  //!< interface side associated with rows of system matrix 2
    const INPAR::S2I::InterfaceSides
        matrix2_side_cols,  //!< interface side associated with columns of system matrix 2
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix3,  //!< system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_rows,  //!< interface side associated with rows of system matrix 3
    const INPAR::S2I::InterfaceSides
        matrix3_side_cols,  //!< interface side associated with columns of system matrix 3
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix4,  //!< system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_rows,  //!< interface side associated with rows of system matrix 4
    const INPAR::S2I::InterfaceSides
        matrix4_side_cols,  //!< interface side associated with columns of system matrix 4
    const Teuchos::RCP<Epetra_MultiVector>& systemvector1,  //!< system vector 1
    const INPAR::S2I::InterfaceSides
        vector1_side,  //!< interface side associated with system vector 1
    const Teuchos::RCP<Epetra_MultiVector>& systemvector2,  //!< system vector 2
    const INPAR::S2I::InterfaceSides
        vector2_side,         //!< interface side associated with system vector 2
    const unsigned nds_rows,  //!< number of dofset associated with matrix rows
    const unsigned nds_cols   //!< number of dofset associated with matrix columns
    )
    : matrix1_side_rows_(matrix1_side_rows),
      matrix1_side_cols_(matrix1_side_cols),
      matrix2_side_rows_(matrix2_side_rows),
      matrix2_side_cols_(matrix2_side_cols),
      matrix3_side_rows_(matrix3_side_rows),
      matrix3_side_cols_(matrix3_side_cols),
      matrix4_side_rows_(matrix4_side_rows),
      matrix4_side_cols_(matrix4_side_cols),
      systemmatrix1_(systemmatrix1),
      systemmatrix2_(systemmatrix2),
      systemmatrix3_(systemmatrix3),
      systemmatrix4_(systemmatrix4),
      systemvector1_(systemvector1),
      systemvector2_(systemvector2),
      vector1_side_(vector1_side),
      vector2_side_(vector2_side),
      nds_rows_(nds_rows),
      nds_cols_(nds_cols)
{
  return;
}


/*----------------------------------------------------------------------------------*
 | assemble cell matrices and vectors into system matrices and vectors   fang 05/16 |
 *----------------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::AssembleCellMatricesAndVectors(
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const int
        assembler_pid_master  //!< ID of processor performing master-side matrix and vector assembly
    ) const
{
  // assemble cell matrix 1 into system matrix 1
  if (AssembleMatrix1())
    AssembleCellMatrix(systemmatrix1_, cellmatrix1_, matrix1_side_rows_, matrix1_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 2 into system matrix 2
  if (AssembleMatrix2())
    AssembleCellMatrix(systemmatrix2_, cellmatrix2_, matrix2_side_rows_, matrix2_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 3 into system matrix 3
  if (AssembleMatrix3())
    AssembleCellMatrix(systemmatrix3_, cellmatrix3_, matrix3_side_rows_, matrix3_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell matrix 4 into system matrix 4
  if (AssembleMatrix4())
    AssembleCellMatrix(systemmatrix4_, cellmatrix4_, matrix4_side_rows_, matrix4_side_cols_,
        la_slave, la_master, assembler_pid_master);

  // assemble cell vector 1 into system vector 1
  if (AssembleVector1())
    AssembleCellVector(
        systemvector1_, cellvector1_, vector1_side_, la_slave, la_master, assembler_pid_master);

  // assemble cell vector 2 into system vector 2
  if (AssembleVector2())
    AssembleCellVector(
        systemvector2_, cellvector2_, vector2_side_, la_slave, la_master, assembler_pid_master);

  return;
}


/*----------------------------------------------------------------------------------*
 | assemble cell matrix into system matrix                               fang 05/16 |
 *----------------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::AssembleCellMatrix(
    const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,  //!< system matrix
    const Epetra_SerialDenseMatrix& cellmatrix,                //!< cell matrix
    const INPAR::S2I::InterfaceSides side_rows,  //!< interface side associated with matrix rows
    const INPAR::S2I::InterfaceSides side_cols,  //!< interface side associated with matrix columns
    DRT::Element::LocationArray& la_slave,       //!< slave-side location array
    DRT::Element::LocationArray& la_master,      //!< master-side location array
    const int assembler_pid_master  //!< ID of processor performing master-side matrix assembly
    ) const
{
  // determine location array associated with matrix columns
  DRT::Element::LocationArray& la_cols = side_cols == INPAR::S2I::side_slave ? la_slave : la_master;

  // assemble cell matrix into system matrix
  switch (side_rows)
  {
    case INPAR::S2I::side_slave:
    {
      systemmatrix->Assemble(-1, la_cols[nds_cols_].stride_, cellmatrix, la_slave[nds_rows_].lm_,
          la_slave[nds_rows_].lmowner_, la_cols[nds_cols_].lm_);
      break;
    }

    case INPAR::S2I::side_master:
    {
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix)
          ->FEAssemble(cellmatrix, la_master[nds_rows_].lm_,
              std::vector<int>(la_master[nds_rows_].lmowner_.size(), assembler_pid_master),
              la_cols[nds_cols_].lm_);
      break;
    }

    default:
    {
      dserror("Invalid interface side!");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------------------*
 | assemble cell vector into system vector                               fang 05/16 |
 *----------------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::AssembleCellVector(
    const Teuchos::RCP<Epetra_MultiVector>& systemvector,  //!< system vector
    const Epetra_SerialDenseVector& cellvector,            //!< cell vector
    const INPAR::S2I::InterfaceSides
        side,  //!< interface side associated with system and cell vectors
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const int assembler_pid_master  //!< ID of processor performing master-side vector assembly
    ) const
{
  // assemble cell vector into system vector
  switch (side)
  {
    case INPAR::S2I::side_slave:
    {
      if (systemvector->NumVectors() != 1)
        dserror("Invalid number of vectors inside Epetra_MultiVector!");
      LINALG::Assemble(*(*systemvector)(nds_rows_), cellvector, la_slave[nds_rows_].lm_,
          la_slave[nds_rows_].lmowner_);

      break;
    }

    case INPAR::S2I::side_master:
    {
      if (assembler_pid_master == systemvector->Comm().MyPID())
        if (Teuchos::rcp_dynamic_cast<Epetra_FEVector>(systemvector)
                ->SumIntoGlobalValues(
                    la_master[nds_rows_].lm_.size(), &la_master[nds_rows_].lm_[0], cellvector.A()))
          dserror("Assembly into master-side system vector not successful!");

      break;
    }

    default:
    {
      dserror("Invalid interface side!");
      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | initialize cell matrices and vectors                           fang 05/16 |
 *---------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::InitCellMatricesAndVectors(
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
)
{
  // initialize system matrix 1
  if (AssembleMatrix1())
    InitCellMatrix(cellmatrix1_, matrix1_side_rows_, matrix1_side_cols_, la_slave, la_master);

  // initialize system matrix 2
  if (AssembleMatrix2())
    InitCellMatrix(cellmatrix2_, matrix2_side_rows_, matrix2_side_cols_, la_slave, la_master);

  // initialize system matrix 3
  if (AssembleMatrix3())
    InitCellMatrix(cellmatrix3_, matrix3_side_rows_, matrix3_side_cols_, la_slave, la_master);

  // initialize system matrix 4
  if (AssembleMatrix4())
    InitCellMatrix(cellmatrix4_, matrix4_side_rows_, matrix4_side_cols_, la_slave, la_master);

  // initialize system vector 1
  if (AssembleVector1()) InitCellVector(cellvector1_, vector1_side_, la_slave, la_master);

  // initialize system vector 2
  if (AssembleVector2()) InitCellVector(cellvector2_, vector2_side_, la_slave, la_master);

  return;
}


/*---------------------------------------------------------------------------*
 | initialize cell matrix                                         fang 05/16 |
 *---------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::InitCellMatrix(
    Epetra_SerialDenseMatrix& cellmatrix,  //!< cell matrix
    const INPAR::S2I::InterfaceSides
        side_rows,  //!< interface side associated with rows of cell matrix
    const INPAR::S2I::InterfaceSides
        side_cols,  //!< interface side associated with columns of cell matrix
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
    ) const
{
  // determine number of matrix rows and number of matrix columns
  const int nrows = side_rows == INPAR::S2I::side_slave ? la_slave[nds_rows_].Size()
                                                        : la_master[nds_rows_].Size();
  const int ncols = side_cols == INPAR::S2I::side_slave ? la_slave[nds_cols_].Size()
                                                        : la_master[nds_cols_].Size();

  // reshape cell matrix if necessary
  if (cellmatrix.M() != nrows or cellmatrix.N() != ncols) cellmatrix.Shape(nrows, ncols);

  // simply zero out otherwise
  else
    memset(cellmatrix.A(), 0., nrows * ncols * sizeof(double));

  return;
}


/*---------------------------------------------------------------------------*
 | initialize cell vector                                         fang 05/16 |
 *---------------------------------------------------------------------------*/
void SCATRA::MortarCellAssemblyStrategy::InitCellVector(
    Epetra_SerialDenseVector& cellvector,   //!< cell vector
    const INPAR::S2I::InterfaceSides side,  //!< interface side associated with cell vector
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
    ) const
{
  // determine number of vector components
  const int ndofs =
      side == INPAR::S2I::side_slave ? la_slave[nds_rows_].Size() : la_master[nds_rows_].Size();

  // reshape cell vector if necessary
  if (cellvector.Length() != ndofs) cellvector.Size(ndofs);

  // simply zero out otherwise
  else
    memset(cellvector.Values(), 0., ndofs * sizeof(double));

  return;
}


// forward declarations
template class SCATRA::MortarCellCalc<DRT::Element::tri3, DRT::Element::tri3>;
template class SCATRA::MortarCellCalc<DRT::Element::tri3, DRT::Element::quad4>;
template class SCATRA::MortarCellCalc<DRT::Element::quad4, DRT::Element::tri3>;
template class SCATRA::MortarCellCalc<DRT::Element::quad4, DRT::Element::quad4>;
