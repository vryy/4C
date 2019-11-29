/*----------------------------------------------------------------------*/
/*! \file
\file ad_fld_fluidbeam_immersed.cpp

\brief Fluid field adapter for immersed fluids (beam)

\level 3

\maintainer Nora Hagmeyer
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "ad_fld_fluidbeam_immersed.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBeamImmersed::FluidBeamImmersed(
    const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluid_ = Teuchos::rcp(new FluidBaseAlgorithm(
                            prbdyn, DRT::Problem::Instance()->FluidDynamicParams(), "fluid", false))
               ->FluidField();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidBeamImmersed::Discretization()
{
  return FluidField()->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& ADAPTER::FluidBeamImmersed::Interface() const
{
  return fluid_->Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBeamImmersed::PrepareTimeStep() { FluidField()->PrepareTimeStep(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBeamImmersed::Update() { FluidField()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBeamImmersed::Output() { FluidField()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidBeamImmersed::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  return FluidField()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBeamImmersed::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidBeamImmersed::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  dserror("RelaxationSolve not yet implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidBeamImmersed::ExtractInterfaceForces()
{
  return FluidField()->ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidBeamImmersed::ExtractInterfaceVelnp()
{
  return FluidField()->ExtractInterfaceVelnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidBeamImmersed::ExtractInterfaceVeln()
{
  return FluidField()->ExtractInterfaceVeln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidBeamImmersed::IntegrateInterfaceShape()
{
  // Actually we do not need this here, because this will be handled in the coupling.
  return FluidField()->IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidBeamImmersed::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}
