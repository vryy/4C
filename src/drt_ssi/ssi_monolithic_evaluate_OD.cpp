/*----------------------------------------------------------------------*/
/*! \file
\brief Evaluation of OD blocks for monolithic SSI
\level 2

\maintainer Stephan Sinzig

 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_evaluate_OD.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::ScatraStructureODCoupling::ScatraStructureODCoupling(
    const Teuchos::RCP<const Epetra_Map>& full_map_scatra,
    const Teuchos::RCP<const Epetra_Map>& full_map_structure,
    Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I>& meshtying_strategy_s2i,
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>& scatra,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper>& structure)
    : full_map_scatra_(full_map_scatra),
      full_map_structure_(full_map_structure),
      meshtying_strategy_s2i_(meshtying_strategy_s2i),
      scatra_(scatra),
      structure_(structure)
{
  return;
}


/*-----------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------*/
void SSI::ScatraStructureODCoupling::EvaluateODBlockScatraStructureDomain(
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
void SSI::ScatraStructureODCoupling::EvaluateODBlockScatraStructureInterface(
    Teuchos::RCP<LINALG::SparseOperator>& scatrastructureinterface)
{
  scatrastructureinterface->Zero();
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
      scatrastructureinterface,  // auxiliary system matrix
      Teuchos::null,             // no additional matrices of vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate scatra-scatra interface coupling
  std::vector<DRT::Condition*> conditions;
  scatra_->ScaTraField()->Discretization()->GetCondition("S2ICoupling", conditions);
  for (auto const& condition : conditions)
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
      scatrastructureinterface->Complete();
      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // finalize auxiliary system matrix
      scatrastructureinterface->Complete(*full_map_structure_, *full_map_scatra_);
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
void SSI::ScatraStructureODCoupling::EvaluateODBlockStructureScatraDomain(
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
