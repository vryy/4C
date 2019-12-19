/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "poromultiphase_base.H"

#include "../drt_porofluidmultiphase/porofluidmultiphase_utils.H"

#include "../drt_inpar/inpar_porofluidmultiphase.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

// new structural time integration
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_factory.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16  |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseBase::PoroMultiPhaseBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      structure_(Teuchos::null),
      fluid_(Teuchos::null),
      struct_zeros_(Teuchos::null),
      solve_structure_(true),
      artery_coupl_(DRT::INPUT::IntegralValue<int>(globaltimeparams, "ARTERY_COUPLING"))
{
}

/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::Init(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
    const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // access the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);

  // build structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase =
      ADAPTER::STR::BuildStructureAlgorithm(structparams);
  adapterbase->Init(globaltimeparams, const_cast<Teuchos::ParameterList&>(structparams), structdis);
  adapterbase->Setup();
  structure_ = Teuchos::rcp_dynamic_cast<::ADAPTER::Structure>(adapterbase->StructureField());

  // initialize zero vector for convenience
  struct_zeros_ = LINALG::CreateVector(*structure_->DofRowMap(), true);
  // do we also solve the structure, this is helpful in case of fluid-scatra coupling without mesh
  // deformation
  solve_structure_ = DRT::INPUT::IntegralValue<int>(algoparams, "SOLVE_STRUCTURE");

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = fluidparams.get<int>("LINEAR_SOLVER");


  // -------------------------------------------------------------------
  // access the fluid discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis(fluid_disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!fluiddis->Filled()) fluiddis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = fluiddis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
      DRT::INPUT::IntegralValue<INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(
          fluidparams, "TIMEINTEGR");

  // build poro fluid time integrator
  Teuchos::RCP<ADAPTER::PoroFluidMultiphase> porofluid =
      POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(timintscheme, fluiddis, linsolvernumber,
          globaltimeparams, fluidparams, DRT::Problem::Instance()->ErrorFile()->Handle(), output);

  // wrap it
  fluid_ = Teuchos::rcp(new ADAPTER::PoroFluidMultiphaseWrapper(porofluid));
  // initialize it
  fluid_->Init(isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra, nearbyelepairs);

  // done.
  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::ReadRestart(int restart)
{
  if (restart)
  {
    // read restart data for structure field (will set time and step internally)
    structure_->ReadRestart(restart);

    // read restart data for fluid field (will set time and step internally)
    fluid_->ReadRestart(restart);

    // reset time and step for the global algorithm
    SetTimeStep(structure_->TimeOld(), restart);
  }


  return;
}

/*----------------------------------------------------------------------*
 | time loop                                            kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::Timeloop()
{
  // prepare the loop
  PrepareTimeLoop();

  // time loop
  while (NotFinished())
  {
    PrepareTimeStep();

    TimeStep();

    UpdateAndOutput();
  }

  return;
}

/*----------------------------------------------------------------------*
 | prepare the time loop                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::PrepareTimeLoop()
{
  // initial output
  if (solve_structure_)
  {
    StructureField()->PrepareOutput();
    StructureField()->Output();
    SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  }
  else
  {
    // Inform user that structure field has been disabled
    PrintStructureDisabledInfo();
    // just set displacements and velocities to zero
    SetStructSolution(struct_zeros_, struct_zeros_);
  }
  FluidField()->PrepareTimeLoop();

  return;
}

/*----------------------------------------------------------------------*
 | prepare one time step                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::PrepareTimeStep()
{
  IncrementTimeAndStep();

  StructureField()->Discretization()->SetState(1, "porofluid", FluidField()->Phinp());

  if (solve_structure_)
  {
    // NOTE: the predictor of the structure is called in here
    StructureField()->PrepareTimeStep();
    SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  }
  else
    SetStructSolution(struct_zeros_, struct_zeros_);

  FluidField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | Test the results of all subproblems                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::CreateFieldTest()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(structure_->CreateFieldTest());
  problem->AddFieldTest(fluid_->CreateFieldTest());
}

/*------------------------------------------------------------------------*
 | communicate the solution of the structure to the fluid    vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetStructSolution(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*------------------------------------------------------------------------*
 | communicate the structure velocity  to the fluid           vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetVelocityFields(Teuchos::RCP<const Epetra_Vector> vel)
{
  fluid_->SetVelocityField(vel);
}

/*------------------------------------------------------------------------*
 | communicate the scatra solution to the fluid             vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetScatraSolution(
    unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars)
{
  fluid_->SetScatraSolution(nds, scalars);
}


/*------------------------------------------------------------------------*
 | communicate the structure displacement to the fluid        vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetMeshDisp(Teuchos::RCP<const Epetra_Vector> disp)
{
  fluid_->ApplyMeshMovement(disp);
}

/*----------------------------------------------------------------------*
 | update fields and output results                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::UpdateAndOutput()
{
  // prepare the output
  StructureField()->PrepareOutput();

  // update single fields
  StructureField()->Update();
  FluidField()->Update();

  // evaluate error if desired
  FluidField()->EvaluateErrorComparedToAnalyticalSol();

  // set structure on fluid (necessary for possible domain integrals)
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // output single fields
  StructureField()->Output();
  FluidField()->Output();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of structure field           vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::StructDofRowMap() const
{
  return structure_->DofRowMap();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of fluid field           vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::FluidDofRowMap() const
{
  return fluid_->DofRowMap();
}

/*------------------------------------------------------------------------*
 | coupled artery-porofluid system matrix                kremheller 05/18 |
 *------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase>
POROMULTIPHASE::PoroMultiPhaseBase::ArteryPorofluidSysmat() const
{
  return fluid_->ArteryPorofluidSysmat();
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of artery field         kremheller 05/18 |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseBase::ArteryDofRowMap() const
{
  return fluid_->ArteryDofRowMap();
}

/*------------------------------------------------------------------------*
 | return structure displacements                             vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::StructDispnp() const
{
  return structure_->Dispnp();
}

/*------------------------------------------------------------------------*
 | return structure velocities                               vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::StructVelnp() const
{
  return structure_->Velnp();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> POROMULTIPHASE::PoroMultiPhaseBase::FluidFlux() const
{
  return fluid_->Flux();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidPhinp() const
{
  return fluid_->Phinp();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidSaturation() const
{
  return fluid_->Saturation();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::FluidPressure() const
{
  return fluid_->Pressure();
}

/*------------------------------------------------------------------------*
 | return fluid Flux                                         vuong 08/16  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROMULTIPHASE::PoroMultiPhaseBase::SolidPressure() const
{
  return fluid_->SolidPressure();
}

/*----------------------------------------------------------------------*
 | inform user that structure is not solved            kremheller 04/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::PrintStructureDisabledInfo()
{
  // print out Info
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
    std::cout << " INFO:    STRUCTURE FIELD IS NOT SOLVED   \n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
  }
}
