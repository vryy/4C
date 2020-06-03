/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed fluids (beam)

\level 3

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fbi_wrapper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "ad_fld_fbi_movingboundary.H"
#include "ad_fld_base_algorithm.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FBIFluidMB::FBIFluidMB(const Teuchos::ParameterList& prbdyn, std::string condname)
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
Teuchos::RCP<DRT::Discretization> ADAPTER::FBIFluidMB::Discretization()
{
  return FluidField()->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& ADAPTER::FBIFluidMB::Interface() const
{
  return fluidadapter_->Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::PrepareTimeStep() { FluidField()->PrepareTimeStep(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::Update() { FluidField()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::Output() { FluidField()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FBIFluidMB::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  return FluidField()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  dserror("RelaxationSolve not yet implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::ExtractInterfaceForces()
{
  return FluidField()->ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::ExtractInterfaceVelnp()
{
  return FluidField()->ExtractInterfaceVelnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::ExtractInterfaceVeln()
{
  return FluidField()->ExtractInterfaceVeln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::IntegrateInterfaceShape()
{
  // Actually we do not need this here, because this will be handled in the coupling.
  return FluidField()->IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FBIFluidMB::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FBIFluidMB::SetCouplingContributions(Teuchos::RCP<const LINALG::SparseMatrix> matrix)
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true)
      ->SetCouplingContributions(matrix);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::ApplyInterfaceValues(
    Teuchos::RCP<Epetra_Vector> iforce, Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->AddContributionToExternalLoads(iforce);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Get velocity at timestep n+1
Teuchos::RCP<const Epetra_Vector> ADAPTER::FBIFluidMB::Velnp() { return FluidField()->Velnp(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::FBIFluidMB::Itemax() const { return fluidadapter_->Itemax(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// set the maximum number of iterations for the fluid field
void ADAPTER::FBIFluidMB::SetItemax(int itemax) { FluidField()->SetItemax(itemax); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::ResetExternalForces()
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true)->ResetExternalForces();
}
