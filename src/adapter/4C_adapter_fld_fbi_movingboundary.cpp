// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
  fluidadapter_ = std::make_shared<FluidBaseAlgorithm>(
      prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid", false)
                      ->fluid_field();
  // make sure
  if (std::dynamic_pointer_cast<Adapter::FluidFBI>(fluid_field()) == nullptr)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::FE::Discretization> Adapter::FBIFluidMB::discretization()
{
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<FLD::Utils::MapExtractor> const& Adapter::FBIFluidMB::interface() const
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
void Adapter::FBIFluidMB::nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
    std::shared_ptr<Core::LinAlg::Vector<double>> ivel)
{
  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::relaxation_solve(
    std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt)
{
  FOUR_C_THROW("RelaxationSolve not yet implemented");
  return nullptr;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::extract_interface_forces()
{
  return fluid_field()->extract_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::extract_interface_velnp()
{
  return fluid_field()->extract_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::extract_interface_veln()
{
  return fluid_field()->extract_interface_veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::integrate_interface_shape()
{
  // Actually we do not need this here, because this will be handled in the coupling.
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::Utils::ResultTest> Adapter::FBIFluidMB::create_field_test()
{
  return fluid_field()->create_field_test();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FBIFluidMB::set_coupling_contributions(
    std::shared_ptr<const Core::LinAlg::SparseOperator> matrix)
{
  std::dynamic_pointer_cast<Adapter::FluidFBI>(fluid_field())->set_coupling_contributions(matrix);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIFluidMB::apply_interface_values(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce,
    std::shared_ptr<Core::LinAlg::Vector<double>> ivel)
{
  fluid_field()->add_contribution_to_external_loads(iforce);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Get velocity at timestep n+1
std::shared_ptr<const Core::LinAlg::Vector<double>> Adapter::FBIFluidMB::velnp()
{
  return fluid_field()->velnp();
}
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
  std::dynamic_pointer_cast<Adapter::FluidFBI>(fluid_field())->reset_external_forces();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const FLD::Meshtying> Adapter::FBIFluidMB::get_meshtying()
{
  return std::dynamic_pointer_cast<Adapter::FluidFBI>(fluid_field())->get_meshtying();
}

FOUR_C_NAMESPACE_CLOSE
