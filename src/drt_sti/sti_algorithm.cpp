/*----------------------------------------------------------------------*/
/*! \file

\brief general coupling algorithm for scatra-thermo interaction

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#include "sti_algorithm.H"

#include <Epetra_Time.h>

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 04/15 |
 *--------------------------------------------------------------------------------*/
STI::Algorithm::Algorithm(const Epetra_Comm& comm,  //! communicator
    const Teuchos::ParameterList& stidyn,           //! parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&
        scatradyn,  //! scalar transport parameter list for scatra and thermo fields
    const Teuchos::ParameterList& solverparams_scatra,  //! solver parameter list for scatra field
    const Teuchos::ParameterList& solverparams_thermo   //! solver parameter list for thermo field
    )
    :  // instantiate base class
      AlgorithmBase(comm, scatradyn),

      scatra_(Teuchos::null),
      thermo_(Teuchos::null),
      strategyscatra_(Teuchos::null),
      strategythermo_(Teuchos::null),
      fieldparameters_(Teuchos::rcp(new Teuchos::ParameterList(scatradyn))),
      iter_(0),
      itermax_(0),
      itertol_(0.),
      stiparameters_(Teuchos::rcp(new Teuchos::ParameterList(stidyn))),

      // initialize timer for Newton-Raphson iteration
      timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
  // check input parameters for scatra and thermo fields
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*fieldparameters_, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_zero)
    dserror("Scatra-thermo interaction with convection not yet implemented!");

  // initialize scatra time integrator
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  scatra->Init(*fieldparameters_, *fieldparameters_, solverparams_scatra);
  scatra->Setup();
  scatra_ = scatra->ScaTraField();

  // modify field parameters for thermo field
  ModifyFieldParametersForThermoField();

  // initialize thermo time integrator
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> thermo =
      Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());
  thermo->Init(*fieldparameters_, *fieldparameters_, solverparams_thermo, "thermo");
  thermo->Setup();
  thermo_ = thermo->ScaTraField();

  // check maps from scatra and thermo discretizations
  if (scatra_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Scatra discretization does not have any degrees of freedom!");
  if (thermo_->Discretization()->DofRowMap()->NumGlobalElements() == 0)
    dserror("Thermo discretization does not have any degrees of freedom!");

  // additional safety check
  if (thermo_->NumScal() != 1) dserror("Thermo field must involve exactly one transported scalar!");

  // perform initializations associated with scatra-scatra interface coupling
  if (scatra_->S2ICoupling())
  {
    // safety check
    if (!thermo_->S2ICoupling())
      dserror(
          "Can't evaluate scatra-scatra interface coupling in scatra field, but not in thermo "
          "field!");

    // extract meshtying strategies for scatra-scatra interface coupling from scatra and thermo time
    // integrators
    strategyscatra_ = Teuchos::rcp_dynamic_cast<SCATRA::MeshtyingStrategyS2I>(scatra_->Strategy());
    strategythermo_ = Teuchos::rcp_dynamic_cast<SCATRA::MeshtyingStrategyS2I>(thermo_->Strategy());

    // perform initializations depending on type of meshtying method
    switch (strategyscatra_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // safety check
        if (strategythermo_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
          dserror(
              "Must have matching nodes at scatra-scatra coupling interfaces in both the scatra "
              "and the thermo fields!");

        break;
      }

      case INPAR::S2I::coupling_mortar_standard:
      {
        // safety check
        if (strategythermo_->CouplingType() != INPAR::S2I::coupling_mortar_condensed_bubnov)
          dserror("Invalid type of scatra-scatra interface coupling for thermo field!");

        // extract scatra-scatra interface coupling conditions
        std::vector<DRT::Condition*> conditions;
        scatra_->Discretization()->GetCondition("S2ICoupling", conditions);

        // loop over all conditions
        for (unsigned icondition = 0; icondition < conditions.size(); ++icondition)
          // consider conditions for slave side only
          if (conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
          {
            // extract ID of current condition
            const int condid = conditions[icondition]->GetInt("ConditionID");
            if (condid < 0) dserror("Invalid condition ID!");

            // extract mortar discretizations associated with current condition
            DRT::Discretization& scatradis = strategyscatra_->MortarDiscretization(condid);
            DRT::Discretization& thermodis = strategythermo_->MortarDiscretization(condid);

            // exchange dofsets between discretizations
            scatradis.AddDofSet(thermodis.GetDofSetProxy());
            thermodis.AddDofSet(scatradis.GetDofSetProxy());
          }

        break;
      }

      default:
      {
        dserror("Invalid type of scatra-scatra interface coupling!");
        break;
      }
    }
  }

  return;
}  // STI::Algorithm::Algorithm


/*----------------------------------------------------------------------*
 | modify field parameters for thermo field                  fang 06/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ModifyFieldParametersForThermoField()
{
  // extract parameters for initial temperature field from parameter list for scatra-thermo
  // interaction and overwrite corresponding parameters in parameter list for thermo field
  if (!fieldparameters_->isParameter("INITIALFIELD") or
      !fieldparameters_->isParameter("INITFUNCNO"))
    dserror(
        "Initial field parameters not properly set in input file section SCALAR TRANSPORT "
        "DYNAMIC!");
  if (!stiparameters_->isParameter("THERMO_INITIALFIELD") or
      !stiparameters_->isParameter("THERMO_INITFUNCNO"))
    dserror(
        "Initial field parameters not properly set in input file section SCALAR TRANSPORT "
        "DYNAMIC!");
  fieldparameters_->set<std::string>(
      "INITIALFIELD", stiparameters_->get<std::string>("THERMO_INITIALFIELD"));
  fieldparameters_->set<int>("INITFUNCNO", stiparameters_->get<int>("THERMO_INITFUNCNO"));

  // perform additional manipulations associated with scatra-scatra interface coupling
  if (scatra_->S2ICoupling())
  {
    // set flag for matrix type associated with thermo field
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("MATRIXTYPE", "sparse");

    // set flag in thermo meshtying strategy for evaluation of interface linearizations and
    // residuals on slave side only
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("SLAVEONLY", "Yes");

    // adapt type of meshtying method for thermo field
    if (fieldparameters_->sublist("S2I COUPLING").get<std::string>("COUPLINGTYPE") ==
        "StandardMortar")
      fieldparameters_->sublist("S2I COUPLING")
          .set<std::string>("COUPLINGTYPE", "CondensedMortar_Bubnov");
    else if (fieldparameters_->sublist("S2I COUPLING").get<std::string>("COUPLINGTYPE") !=
             "MatchingNodes")
      dserror("Invalid type of scatra-scatra interface coupling!");

    // make sure that interface side underlying Lagrange multiplier definition is slave side
    fieldparameters_->sublist("S2I COUPLING").set<std::string>("LMSIDE", "slave");
  }

  return;
}  // STI::Algorithm::ModifyFieldParametersForThermoField()


/*----------------------------------------------------------------------*
 | output solution to screen and files                       fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::Output()
{
  // output scatra field
  scatra_->Output();

  // output thermo field
  thermo_->Output();

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // provide scatra and thermo fields with velocities
  scatra_->SetVelocityField(1);
  thermo_->SetVelocityField(1);

  // pass thermo degrees of freedom to scatra discretization for preparation of first time step
  // (calculation of initial time derivatives etc.)
  if (Step() == 1) TransferThermoToScatra(thermo_->Phiafnp());

  // prepare time step for scatra field
  scatra_->PrepareTimeStep();

  // pass scatra degrees of freedom to thermo discretization for preparation of first time step
  // (calculation of initial time derivatives etc.)
  if (Step() == 1) TransferScatraToThermo(scatra_->Phiafnp());

  // prepare time step for thermo field
  thermo_->PrepareTimeStep();

  return;
}  // STI::Algorithm::PrepareTimeStep()


/*----------------------------------------------------------------------*
 | read restart data                                         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::ReadRestart(int step  //! time step for restart
)
{
  // read scatra and thermo restart variables
  scatra_->ReadRestart(step);
  thermo_->ReadRestart(step);

  // set time and time step
  SetTimeStep(scatra_->Time(), step);

  return;
}  // STI::Algorithm::ReadRestart


/*----------------------------------------------------------------------*
 | time loop                                                 fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::Algorithm::TimeLoop()
{
  // output initial solution to screen and files
  if (Step() == 0) Output();

  // time loop
  while (NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    double time = timer_->WallTime();

    // evaluate time step
    Solve();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->WallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if applicable
    if (DRT::INPUT::IntegralValue<int>(*fieldparameters_, "OUTPUTNONLINSOLVERSTATS"))
      scatra_->OutputNonlinSolverStats(iter_, dtnonlinsolve, Step(), Comm());

    // update scatra and thermo fields
    Update();

    // output solution to screen and files
    Output();
  }  // while(NotFinished())

  return;
}  // STI::Algorithm::TimeLoop()


/*-------------------------------------------------------------------------------------*
 | pass scatra degrees of freedom to thermo discretization                  fang 09/17 |
 *-------------------------------------------------------------------------------------*/
void STI::Algorithm::TransferScatraToThermo(
    const Teuchos::RCP<const Epetra_Vector> scatra  //!< scatra state vector
    ) const
{
  // pass scatra degrees of freedom to thermo discretization
  thermo_->Discretization()->SetState(2, "scatra", scatra);

  // transfer state vector for evaluation of scatra-scatra interface coupling
  if (thermo_->S2ICoupling())
  {
    switch (strategythermo_->CouplingType())
    {
      case INPAR::S2I::coupling_matching_nodes:
      {
        // pass master-side scatra degrees of freedom to thermo discretization
        const Teuchos::RCP<Epetra_Vector> imasterphinp =
            LINALG::CreateVector(*scatra_->Discretization()->DofRowMap(), true);
        strategyscatra_->InterfaceMaps()->InsertVector(
            strategyscatra_->CouplingAdapter()->MasterToSlave(
                strategyscatra_->InterfaceMaps()->ExtractVector(*scatra, 2)),
            1, imasterphinp);
        thermo_->Discretization()->SetState(2, "imasterscatra", imasterphinp);

        break;
      }

      case INPAR::S2I::coupling_mortar_condensed_bubnov:
      {
        // extract scatra-scatra interface coupling conditions
        std::vector<DRT::Condition*> conditions;
        thermo_->Discretization()->GetCondition("S2ICoupling", conditions);

        // loop over all conditions
        for (unsigned icondition = 0; icondition < conditions.size(); ++icondition)
          // consider conditions for slave side only
          if (conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
          {
            // extract ID of current condition
            const int condid = conditions[icondition]->GetInt("ConditionID");
            if (condid < 0) dserror("Invalid condition ID!");

            // extract mortar discretization associated with current condition
            DRT::Discretization& thermodis = strategythermo_->MortarDiscretization(condid);

            // pass interfacial scatra degrees of freedom to thermo discretization
            const Teuchos::RCP<Epetra_Vector> iscatra =
                Teuchos::rcp(new Epetra_Vector(*thermodis.DofRowMap(1)));
            LINALG::Export(*scatra, *iscatra);
            thermodis.SetState(1, "scatra", iscatra);
          }

        break;
      }

      default:
      {
        dserror("You must be kidding me...");
        break;
      }
    }
  }

  return;
}  // STI::Algorithm::TransferScatraToThermo()


/*-------------------------------------------------------------------------------------*
 | pass thermo degrees of freedom to scatra discretization                  fang 09/17 |
 *-------------------------------------------------------------------------------------*/
void STI::Algorithm::TransferThermoToScatra(
    const Teuchos::RCP<const Epetra_Vector> thermo  //!< thermo state vector
    ) const
{
  // pass thermo degrees of freedom to scatra discretization
  scatra_->Discretization()->SetState(2, "thermo", thermo);

  // transfer state vector for evaluation of scatra-scatra interface coupling
  if (scatra_->S2ICoupling() and
      strategyscatra_->CouplingType() == INPAR::S2I::coupling_mortar_standard)
  {
    // extract scatra-scatra interface coupling conditions
    std::vector<DRT::Condition*> conditions;
    scatra_->Discretization()->GetCondition("S2ICoupling", conditions);

    // loop over all conditions
    for (unsigned icondition = 0; icondition < conditions.size(); ++icondition)
      // consider conditions for slave side only
      if (conditions[icondition]->GetInt("interface side") == INPAR::S2I::side_slave)
      {
        // extract ID of current condition
        const int condid = conditions[icondition]->GetInt("ConditionID");
        if (condid < 0) dserror("Invalid condition ID!");

        // extract mortar discretization associated with current condition
        DRT::Discretization& scatradis = strategyscatra_->MortarDiscretization(condid);

        // pass interfacial thermo degrees of freedom to scatra discretization
        const Teuchos::RCP<Epetra_Vector> ithermo =
            Teuchos::rcp(new Epetra_Vector(*scatradis.DofRowMap(1)));
        LINALG::Export(*thermo, *ithermo);
        scatradis.SetState(1, "thermo", ithermo);
      }
  }

  return;
}  // STI::Algorithm::TransferThermoToScatra()


/*-------------------------------------------------------------------------*
 | update scatra and thermo fields after time step evaluation   fang 04/15 |
 *-------------------------------------------------------------------------*/
void STI::Algorithm::Update()
{
  // update scatra field
  scatra_->Update();

  // compare scatra field to analytical solution if applicable
  scatra_->EvaluateErrorComparedToAnalyticalSol();

  // update thermo field
  thermo_->Update();

  // compare thermo field to analytical solution if applicable
  thermo_->EvaluateErrorComparedToAnalyticalSol();

  return;
}  // STI::Algorithm::Update()
