/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to output at runtime for the fluid field

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_fluid_discretization_runtime_output_params.hpp"

#include "baci_inpar_fluid.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeOutputParams::Init(
    const Teuchos::ParameterList& IO_fluid_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values

  output_velocity_state_ = CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "VELOCITY");
  output_pressure_state_ = CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "PRESSURE");
  output_acceleration_state_ =
      CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "ACCELERATION");
  output_displacement_state_ =
      CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "DISPLACEMENT");
  output_gridvelocity_state_ =
      CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "GRIDVELOCITY");
  output_element_owner_ = CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "ELEMENT_OWNER");
  output_element_gid_ = CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "ELEMENT_GID");
  output_node_gid_ = CORE::UTILS::IntegralValue<bool>(IO_fluid_paramslist, "NODE_GID");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeOutputParams::Setup()
{
  if (not IsInit()) FOUR_C_THROW("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeOutputParams::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) FOUR_C_THROW("Call Init() and Setup() first!");
}

FOUR_C_NAMESPACE_CLOSE
