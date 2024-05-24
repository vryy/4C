/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to output at runtime for the fluid field

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_fluid_discretization_runtime_output_params.hpp"

#include "4C_inpar_fluid.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

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
  if (not is_init()) FOUR_C_THROW("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeOutputParams::check_init_setup() const
{
  if (not is_init() or not is_setup()) FOUR_C_THROW("Call Init() and Setup() first!");
}

FOUR_C_NAMESPACE_CLOSE
