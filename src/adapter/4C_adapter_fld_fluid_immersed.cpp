/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for immersed fluids

\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_immersed.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidImmersed::FluidImmersed(const Teuchos::ParameterList& prbdyn, std::string condname)
{
  fluid_ = Teuchos::rcp(new FluidBaseAlgorithm(
      prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid", false));
  fluidadapter_ = fluid_->fluid_field();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> Adapter::FluidImmersed::discretization()
{
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& Adapter::FluidImmersed::interface() const
{
  return fluidadapter_->interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::update() { fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::output() { fluid_field()->statistics_and_output(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::FluidImmersed::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::nonlinear_solve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidImmersed::relaxation_solve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->relaxation_solve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidImmersed::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidImmersed::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidImmersed::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FluidImmersed::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> Adapter::FluidImmersed::create_field_test()
{
  return fluid_field()->create_field_test();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::add_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  fluid_field()->add_dirich_cond(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidImmersed::remove_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  fluid_field()->remove_dirich_cond(maptoremove);
}

FOUR_C_NAMESPACE_CLOSE
