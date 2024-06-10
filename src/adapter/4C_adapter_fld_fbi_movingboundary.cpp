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
Adapter::FBIFluidMB::FBIFluidMB(const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluidadapter_ =
      Teuchos::rcp(new FluidBaseAlgorithm(
                       prbdyn, Global::Problem::Instance()->FluidDynamicParams(), "fluid", false))
          ->fluid_field();
  // make sure
  if (Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true) == Teuchos::null)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> Adapter::FBIFluidMB::discretization()
{
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& Adapter::FBIFluidMB::Interface() const
{
  return fluidadapter_->Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::Update() { fluid_field()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::Output() { fluid_field()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::FBIFluidMB::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::nonlinear_solve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  fluid_field()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  FOUR_C_THROW("RelaxationSolve not yet implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::extract_interface_velnp()
{
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::integrate_interface_shape()
{
  // Actually we do not need this here, because this will be handled in the coupling.
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> Adapter::FBIFluidMB::CreateFieldTest()
{
  return fluid_field()->CreateFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FBIFluidMB::set_coupling_contributions(
    Teuchos::RCP<const Core::LinAlg::SparseOperator> matrix)
{
  Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true)
      ->set_coupling_contributions(matrix);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::apply_interface_values(
    Teuchos::RCP<Epetra_Vector> iforce, Teuchos::RCP<Epetra_Vector> ivel)
{
  fluid_field()->add_contribution_to_external_loads(iforce);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Get velocity at timestep n+1
Teuchos::RCP<const Epetra_Vector> Adapter::FBIFluidMB::Velnp() { return fluid_field()->Velnp(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::FBIFluidMB::Itemax() const { return fluidadapter_->Itemax(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// set the maximum number of iterations for the fluid field
void Adapter::FBIFluidMB::SetItemax(int itemax) { fluid_field()->SetItemax(itemax); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::ResetExternalForces()
{
  Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true)->ResetExternalForces();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const FLD::Meshtying> Adapter::FBIFluidMB::GetMeshtying()
{
  return Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true)->GetMeshtying();
}

FOUR_C_NAMESPACE_CLOSE
