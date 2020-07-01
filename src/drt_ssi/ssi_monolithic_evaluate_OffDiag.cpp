/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of off-diagonal blocks for monolithic SSI
\level 2


 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_evaluate_OffDiag.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_matrixtransform.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraStructureOffDiagCoupling::ScatraStructureOffDiagCoupling(
    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapsscatra,
    Teuchos::RCP<const LINALG::MultiMapExtractor> blockmapstructure,
    const Teuchos::RCP<const Epetra_Map> full_map_scatra,
    const Teuchos::RCP<const Epetra_Map> full_map_structure,
    Teuchos::RCP<ADAPTER::Coupling> icoup_structure,
    const Teuchos::RCP<const Epetra_Map> interface_map_scatra,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtying_strategy_s2i,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
    : blockmapsscatra_(blockmapsscatra),
      blockmapstructure_(blockmapstructure),
      full_map_scatra_(full_map_scatra),
      full_map_structure_(full_map_structure),
      icoup_structure_(icoup_structure),
      interface_map_scatra_(interface_map_scatra),
      meshtying_strategy_s2i_(meshtying_strategy_s2i),
      scatra_(scatra),
      structure_(structure)
{
  return;
}


/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureDomain(
    Teuchos::RCP<LINALG::SparseOperator>& scatrastructureblock)
{
  // initialize scatra-structure matrix block
  scatrastructureblock->Zero();

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
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

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
  scatra_->ScaTraField()->Discretization()->Evaluate(eleparams, strategyscatrastructure);

  // finalize scatra-structure matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      scatrastructureblock->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      scatrastructureblock->Complete(*full_map_structure_, *full_map_scatra_);
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

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrastructureinterface)
{
  scatrastructureinterface->Zero();

  // slave and master matrix for evaluation of conditions
  Teuchos::RCP<LINALG::SparseOperator> slavematrix(Teuchos::null);
  Teuchos::RCP<LINALG::SparseOperator> mastermatrix(Teuchos::null);
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, meshtying_strategy_s2i_->BlockMapsSlave(), 81, false, true));
      mastermatrix = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, meshtying_strategy_s2i_->BlockMapsMaster(), 81, false, true));
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
  EvaluateScatraStructureInterfaceSlaveSide(slavematrix);

  // copy interface contributions from slave side to master side
  CopySlaveToMasterScatraStructureInterface(slavematrix, mastermatrix);

  // add contributions from slave side and master side
  scatrastructureinterface->Add(*slavematrix, false, 1.0, 1.0);
  scatrastructureinterface->Add(*mastermatrix, false, 1.0, 1.0);

  // finalize scatra-structure matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      scatrastructureinterface->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // finalize auxiliary system matrix
      scatrastructureinterface->Complete(*full_map_structure_, *interface_map_scatra_);
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

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::EvaluateOffDiagBlockStructureScatraDomain(
    Teuchos::RCP<LINALG::SparseOperator>& structurescatradomain) const
{
  // initialize structure-scatra matrix block
  structurescatradomain->Zero();

  // create parameter list for element evaluation and fill it
  Teuchos::ParameterList eleparams;
  // set action
  eleparams.set("action", "calc_struct_stiffscalar");
  // set time
  eleparams.set<double>("total time", structure_->Time());
  // set numscatradofspernode
  eleparams.set<int>("numscatradofspernode", scatra_->ScaTraField()->NumDofPerNode());

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

  // finalize structure-scatra matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      structurescatradomain->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      structurescatradomain->Complete(*full_map_scatra_, *full_map_structure_);
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

/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureOffDiagCoupling::CopySlaveToMasterScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator> slavematrix,
    Teuchos::RCP<LINALG::SparseOperator>& mastermatrix)
{
  mastermatrix->Zero();
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      const int numberscatrablocks = scatra_->ScaTraField()->BlockMaps().NumMaps();

      // cast scatrastructureinterfaceslaveside
      // cast master and slave matrix
      auto blockslavematrix = Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(slavematrix);
      auto blockmastermatrix =
          Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(mastermatrix);

      // initialize auxiliary system matrix for linearizations of master-side scatra fluxes w.r.t.
      // master-side structural dofs
      LINALG::SparseMatrix mastermatrixsparse(
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(), 27, false, true);

      // derive linearizations of master-side scatra fluxes w.r.t. master-side structural dofs and
      // assemble into auxiliary system matrix
      for (int iblock = 0; iblock < numberscatrablocks; ++iblock)
        LINALG::MatrixRowColTransform()(blockslavematrix->Matrix(iblock, 0), -1.0,
            ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()),
            ADAPTER::CouplingSlaveConverter(*icoup_structure_), mastermatrixsparse, true, true);

      // finalize auxiliary system matrix
      mastermatrixsparse.Complete(*icoup_structure_->MasterDofMap(),
          *meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap());

      // split auxiliary system matrix and assemble into scatra-structure matrix block
      blockmastermatrix = mastermatrixsparse.Split<LINALG::DefaultBlockMatrixStrategy>(
          *blockmapstructure_, *blockmapsscatra_);
      mastermatrix->Complete();

      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // cast scatrastructureinterfaceslaveside
      auto sparseslavematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(slavematrix);
      auto sparsemastermatrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mastermatrix);

      // copy slave side values to master side and scale with minus 1. Insert into
      // scatrastructureinterface_sparse
      LINALG::MatrixRowColTransform()(*sparseslavematrix, -1.0,
          ADAPTER::CouplingSlaveConverter(*meshtying_strategy_s2i_->CouplingAdapter()),
          ADAPTER::CouplingSlaveConverter(*icoup_structure_), *sparsemastermatrix, true, true);

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

  return;
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

  // number of dofset associated with displacement-related dofs on scalar transport discretization
  condparams.set<int>("ndsdisp", 1);

  // remove state vectors from scalar transport discretization
  scatra_->ScaTraField()->Discretization()->ClearState();

  // add state vectors to scalar transport discretization
  scatra_->ScaTraField()->AddTimeIntegrationSpecificVectors();

  // create strategy for assembly of auxiliary system matrix
  DRT::AssembleStrategy strategyscatrastructures2i(
      0,  // row assembly based on number of dofset associated with scalar transport dofs on
          // scalar transport discretization
      1,  // column assembly based on number of dofset associated with structural dofs on
          // structural discretization
      slavematrix,    // auxiliary system matrix
      Teuchos::null,  // no additional matrices of vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  std::vector<DRT::Condition*> conditions;
  scatra_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (const auto& condition : conditions)
    if (condition->GetInt("interface side") == INPAR::S2I::side_slave)
    {
      // collect condition specific data and store to scatra boundary parameter class
      meshtying_strategy_s2i_->SetConditionSpecificScaTraParameters(*condition);
      // evaluate the condition now
      scatra_->ScaTraField()->Discretization()->EvaluateCondition(
          condparams, strategyscatrastructures2i, "S2ICoupling", condition->GetInt("ConditionID"));
    }

  // finalize scatra-structure matrix block
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      slavematrix->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // finalize auxiliary system matrix
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
  scatra_->ScaTraField()->Discretization()->ClearState();

  return;
}
