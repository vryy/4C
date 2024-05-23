/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed fluids (beam)

\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fbi_movingboundary.hpp"

#include "4C_adapter_fld_base_algorithm.hpp"
#include "4C_adapter_fld_fbi_wrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_linalg_sparseoperator.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FBIFluidMB::FBIFluidMB(const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluidadapter_ =
      Teuchos::rcp(new FluidBaseAlgorithm(
                       prbdyn, GLOBAL::Problem::Instance()->FluidDynamicParams(), "fluid", false))
          ->FluidField();
  // make sure
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true) == Teuchos::null)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");

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
  FOUR_C_THROW("RelaxationSolve not yet implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::extract_interface_forces()
{
  return FluidField()->extract_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::extract_interface_velnp()
{
  return FluidField()->extract_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::extract_interface_veln()
{
  return FluidField()->extract_interface_veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIFluidMB::integrate_interface_shape()
{
  // Actually we do not need this here, because this will be handled in the coupling.
  return FluidField()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> ADAPTER::FBIFluidMB::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FBIFluidMB::set_coupling_contributions(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> matrix)
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true)
      ->set_coupling_contributions(matrix);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIFluidMB::apply_interface_values(
    Teuchos::RCP<Epetra_Vector> iforce, Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->add_contribution_to_external_loads(iforce);
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
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const FLD::Meshtying> ADAPTER::FBIFluidMB::GetMeshtying()
{
  return Teuchos::rcp_dynamic_cast<ADAPTER::FluidFBI>(FluidField(), true)->GetMeshtying();
}

FOUR_C_NAMESPACE_CLOSE
