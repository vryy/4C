/*------------------------------------------------------------------------------------------------*/
/*!
\file adapter_fluid_combust.cpp

\brief Fluid field adapter

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_combust.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_combust/combust_fluidresulttest.H"
#include "../drt_combust/combust_interface.H"

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 08/08 |
 |
 | It initializes its associated time integration scheme fluid_.
 *------------------------------------------------------------------------------------------------*/
ADAPTER::FluidCombust::FluidCombust(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<ParameterList> params,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : fluid_(dis, *solver, *params, *output), // calls COMBUST::CombustFluidImplicitTimeInt()
    fluiddis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  // std::cout << "Proc " << fluiddis_->Comm().MyPID() << "ADAPTER::FluidCombust constructor done \n" << endl;
}

void ADAPTER::FluidCombust::SetInitialFlowField(const INPAR::FLUID::InitialField initfield, const int startfuncno)
{
  // This function is called from the Fluid Base Algorithm, but has no effect.
  // The fluid flow field is initialized from the COMBUST::Algorithm.
  return;
}

void ADAPTER::FluidCombust::SetInitialFlowField(const INPAR::COMBUST::InitialField initfield, const int initfuncno)
{
  return fluid_.SetInitialFlowField(initfield, initfuncno);
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::TrueResidual()
{
  return fluid_.TrueResidual();
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::Veln()
{
  return fluid_.Veln();
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::Velnp()
{
  return fluid_.Velnp();
}

/*------------------------------------------------------------------------------------------------*
 | return history vector                                                          rasthofer 01/10 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::Hist()
{
  return fluid_.Hist();
}

/*------------------------------------------------------------------------------------------------*
 |                                                                                    bauer 04/10 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::TimeLoop()
{
  fluid_.TimeLoop();
}

/*------------------------------------------------------------------------------------------------*
 |                                                                                    henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();
}

void ADAPTER::FluidCombust::ClearTimeInt()
{
	fluid_.ClearTimeInt();
}

/*------------------------------------------------------------------------------------------------*
 | Wozu ist diese Abfrage nötig?                                                      henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  dserror("Thou shalt not call this function!");
  if (stepinc!=Teuchos::null)
  {
    fluid_.Evaluate(stepinc);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}

/*------------------------------------------------------------------------------------------------*
 | In der Zeitintegration gibt es keine Update Funktion, das passiert hier!          henke 10/08  |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Update()
{
  fluid_.TimeUpdate();
}

/*------------------------------------------------------------------------------------------------*
 | Was davon brauche ich? henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Output()
{
  fluid_.Output();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::NonlinearSolve()
{

  // std::cout << "Proc " << fluiddis_->Comm().MyPID()<< "FluidCombust::NonlinearSolve()" << endl;
  fluid_.NonlinearSolve();
}

/*----------------------------------------------------------------------*
 |                                                             vg 04/11 |
 *----------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Predictor()
{
  fluid_.Predictor();
}


/*----------------------------------------------------------------------*
 |                                                             vg 04/11 |
 *----------------------------------------------------------------------*/
void ADAPTER::FluidCombust::MultiCorrector()
{
  fluid_.MultiCorrector();
}


/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidCombust::VelocityRowMap()
{
  return fluid_.VelocityRowMap();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidCombust::PressureRowMap()
{
  return fluid_.PressureRowMap();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::ResidualScaling() const
{
  dserror("Thou shalt not call this function!");
  return 0.0;
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::TimeScaling() const
{
  return 1./fluid_.Dt();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::Time() const
{
  return fluid_.Time();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
int ADAPTER::FluidCombust::Step() const
{
  return fluid_.Step();
}

double ADAPTER::FluidCombust::Dt() const
{
  return fluid_.Dt();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
int ADAPTER::FluidCombust::Itemax() const
{
  return fluid_.Itemax();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::SetItemax(int itemax)
{
  dserror("Thou shalt not call this function!");
  //fluid_.SetItemax(itemax);
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::LiftDrag()
{
  fluid_.LiftDrag();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidCombust::ExtractInterfaceForces()
{
  dserror("Don't use this function, I don't know what it does!");
  return Teuchos::null;
}

/*------------------------------------------------------------------------------------------------*
 | henke 07/09 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidCombust::ExtractInterfaceVeln()
{
  // get convection velocity vector (including pressure) for transfer to scalar transport field
  return fluid_.ConVelnp();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  dserror("Don't use this function, I don't know what it does!");
  return;
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidCombust::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::CombustFluidResultTest(fluid_));
}

/*------------------------------------------------------------------------------------------------*
 | Zum Exportieren der Konvektionsgeschw. Warum die Bezeichnung velpres? Druck?       henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}

void ADAPTER::FluidCombust::ImportDiscretization(Teuchos::RCP<DRT::Discretization> importdis)
{
  /* This function is not needed, because the ImportInterface is used instead. Besides, it does not work
   * anymore, since the InterfaceHandle receives a FlameFront and the discretizations. henke 01/09 */
  dserror("Thou shalt not call this function!");

  // import level set discretization from combustion algorithm
  gfuncdis_ = importdis;

  // construct a combustion interface handle (couldn't we also construct a flame front here!) and
  // thus automatically process the flame front
//  interfacehandle = rcp(new COMBUST::InterfaceHandleCombust(fluiddis_,gfuncdis_));

  // pass geometrical information aboout flame front to fluid time integration scheme
//  fluid_.IncorporateInterface(interfacehandle);
}

void ADAPTER::FluidCombust::ImportInterface(
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle,
    const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle_old)
{
  // pass geometrical information about flame front to fluid time integration scheme
  fluid_.IncorporateInterface(interfacehandle,interfacehandle_old);
}

void ADAPTER::FluidCombust::ImportFlameFront(const Teuchos::RCP<COMBUST::FlameFront>& flamefront)
{
  // pass geometrical information aboout flame front to fluid time integration scheme
  fluid_.StoreFlameFront(flamefront);
}

#endif  // #ifdef CCADISCRET
