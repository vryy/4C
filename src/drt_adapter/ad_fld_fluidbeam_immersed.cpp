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
#include "ad_fld_fbi_wrapper.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBeamImmersed::FluidBeamImmersed(
    const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluidadapter_ = Teuchos::rcp(new FluidBaseAlgorithm(prbdyn,
                                   DRT::Problem::Instance()->FluidDynamicParams(), "fluid", false))
                      ->FluidField();
  // make sure
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true) == Teuchos::null)
    dserror("Failed to create the correct underlying fluid adapter");

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
  return fluidadapter_->Interface();
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FluidBeamImmersed::SetCouplingContributions(Teuchos::RCP<LINALG::SparseMatrix> matrix)
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true)
      ->SetCouplingContributions(matrix);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FluidBeamImmersed::ApplyInterfaceValues(
    Teuchos::RCP<Epetra_Vector> iforce, Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->AddContributionToExternalLoads(iforce);
}
