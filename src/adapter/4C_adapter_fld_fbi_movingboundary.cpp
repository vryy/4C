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
                       prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid", false))
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
Teuchos::RCP<FLD::UTILS::MapExtractor> const& Adapter::FBIFluidMB::interface() const
{
  return fluidadapter_->interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::update() { fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::output() { fluid_field()->statistics_and_output(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::FBIFluidMB::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::nonlinear_solve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIFluidMB::relaxation_solve(
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
Teuchos::RCP<Core::UTILS::ResultTest> Adapter::FBIFluidMB::create_field_test()
{
  return fluid_field()->create_field_test();
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
Teuchos::RCP<const Epetra_Vector> Adapter::FBIFluidMB::velnp() { return fluid_field()->velnp(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int Adapter::FBIFluidMB::itemax() const { return fluidadapter_->itemax(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// set the maximum number of iterations for the fluid field
void Adapter::FBIFluidMB::set_itemax(int itemax) { fluid_field()->set_itemax(itemax); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::reset_external_forces()
{
  Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true)->reset_external_forces();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const FLD::Meshtying> Adapter::FBIFluidMB::get_meshtying()
{
  return Teuchos::rcp_dynamic_cast<Adapter::FluidFBI>(fluid_field(), true)->get_meshtying();
}

FOUR_C_NAMESPACE_CLOSE
