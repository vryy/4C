/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTU output at runtime for the fluid field

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "fluid_discretization_runtime_vtu_output_params.H"

#include "dserror.H"

#include "inpar_parameterlist_utils.H"
#include "inpar_fluid.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeVtuOutputParams::Init(
    const Teuchos::ParameterList& IO_vtk_fluid_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values

  output_velocity_state_ = DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "VELOCITY");
  output_pressure_state_ = DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "PRESSURE");
  output_acceleration_state_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "ACCELERATION");
  output_displacement_state_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "DISPLACEMENT");
  output_gridvelocity_state_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "GRIDVELOCITY");
  output_element_owner_ = DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "ELEMENT_OWNER");
  output_element_gid_ = DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "ELEMENT_GID");
  output_node_gid_ = DRT::INPUT::IntegralValue<bool>(IO_vtk_fluid_paramslist, "NODE_GID");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeVtuOutputParams::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidRuntimeVtuOutputParams::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
