/*----------------------------------------------------------------------*/
/*!
 \file ad_porofluidmultiphase.cpp

 \brief adapter for multiphase porous fluid problem

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/



#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_porofluidmultiphase.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// general time integration schemes
#include "../drt_porofluidmultiphase/porofluidmultiphase_timint_implicit.H"
#include "../drt_porofluidmultiphase/porofluidmultiphase_timint_ost.H"
#include "ad_porofluidmultiphase.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::PoroFluidMultiphase::PoroFluidMultiphase(
    )
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::Init(
    const Teuchos::ParameterList&   globalparams,         ///< parameter list of global control algorithm
    const Teuchos::ParameterList&   porofluiddyn,         ///< parameter list for multiphase porous flow  subproblem
    const int                       linsolvernumber,      ///< number of linear solver
    const std::string&              disname,              ///< name of multiphase porous flow discretization
    const bool                      isale ,                ///< ALE flag
    const int nds_disp,
    const int nds_vel,
    const int nds_solidpressure
    )
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
    DRT::INPUT::IntegralValue<INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(porofluiddyn,"TIMEINTEGR");

  switch(timintscheme)
  {
  case INPAR::POROFLUIDMULTIPHASE::timeint_one_step_theta:
  {
    // create instance of time integration class (call the constructor)
    porofluid_ =
        Teuchos::rcp(new POROFLUIDMULTIPHASE::TimIntOneStepTheta(
            actdis,
            linsolvernumber,
            globalparams,
            porofluiddyn,
            DRT::Problem::Instance()->ErrorFile()->Handle(),
            output));
    break;
  }
  default:
    dserror("Unknown time-integration scheme for multiphase poro fluid problem");
    break;
  }

  // initialize algorithm for specific time-integration scheme
  porofluid_->Init(
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::PoroFluidMultiphase::CreateFieldTest()
{
  return porofluid_->CreateFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::PoroFluidMultiphase::DofRowMap(unsigned nds) const
{
  return porofluid_->DofRowMap(nds);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::ReadRestart(int restart)
{
  return porofluid_->ReadRestart(restart);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> dispnp //!< displacement vector
    )
{
  porofluid_->ApplyMeshMovement(dispnp);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::SetVelocityField(
    Teuchos::RCP<const Epetra_Vector>   vel                       //!< velocity vector
   )
{
  porofluid_->SetVelocityField(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::PoroFluidMultiphase::Phinp() const
{
  return porofluid_->Phinp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::PoroFluidMultiphase::SolidPressure() const
{
  return porofluid_->SolidPressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::PoroFluidMultiphase::Pressure() const
{
  return porofluid_->Pressure();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::PoroFluidMultiphase::Saturation() const
{
  return porofluid_->Saturation();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> ADAPTER::PoroFluidMultiphase::Flux() const
{
  return porofluid_->Flux();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::PrepareTimeStep()
{
  porofluid_->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::Output()
{
  porofluid_->Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::Update()
{
  porofluid_->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::EvaluateErrorComparedToAnalyticalSol()
{
  porofluid_->EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::Solve()
{
  porofluid_->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::PoroFluidMultiphase::PrepareTimeLoop()
{
  porofluid_->PrepareTimeLoop();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::PoroFluidMultiphase::GetDofSetNumberOfSolidPressure()
{
  return porofluid_->GetDofSetNumberOfSolidPressure();
}
