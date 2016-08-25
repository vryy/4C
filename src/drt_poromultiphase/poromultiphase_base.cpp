/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_base.cpp

 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "poromultiphase_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/ad_porofluidmultiphase.H"
#include "../drt_adapter/ad_str_wrapper.H"

#include "../linalg/linalg_utils.H"



/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16  |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseBase::PoroMultiPhaseBase(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    AlgorithmBase(comm, globaltimeparams),
    structure_(Teuchos::null),
    fluid_(Teuchos::null),
    zeros_(Teuchos::null)
{

}

/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure)
{
  // access the global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // Create the two uncoupled subproblems.
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis(struct_disname);

  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
          globaltimeparams,
          structparams,
          structdis));

  // build structural time integrator
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::Structure>(structure->StructureField());

  // initialize zero vector for convenience
  zeros_ = LINALG::CreateVector(*structure_->DofRowMap(), true);

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = fluidparams.get<int>("LINEAR_SOLVER");

  // build poro fluid time integrator
  fluid_ = Teuchos::rcp(new ADAPTER::PoroFluidMultiphase());
  // initialize it
  fluid_->Init(
      globaltimeparams,
      fluidparams,
      linsolvernumber,
      fluid_disname,
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure);

  //done.
  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::ReadRestart( int restart )
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
    Teuchos::RCP<const Epetra_Vector> disp,
    Teuchos::RCP<const Epetra_Vector> vel )
{
  SetMeshDisp(disp);
  SetVelocityFields(vel);
}

/*------------------------------------------------------------------------*
 | communicate the structure velocity  to the fluid           vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetVelocityFields(
    Teuchos::RCP<const Epetra_Vector> vel)
{
  fluid_->SetVelocityField(vel);
}

/*------------------------------------------------------------------------*
 | communicate the structure displacement to the fluid        vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetMeshDisp(
    Teuchos::RCP<const Epetra_Vector> disp )
{
  fluid_->ApplyMeshMovement(disp);
}

/*------------------------------------------------------------------------*
 | communicate the pressure of the fluid to the structure    vuong 08/16  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseBase::SetSolidPressure(
    Teuchos::RCP<const Epetra_Vector> pressure )
{
  const int nds_solidpressure = fluid_->GetDofSetNumberOfSolidPressure();
  structure_->Discretization()->SetState(nds_solidpressure,"solid_pressure",pressure);
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
