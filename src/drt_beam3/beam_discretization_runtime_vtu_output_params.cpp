/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_discretization_runtime_vtu_output_params.cpp

\brief input parameters related to VTU output at runtime for beams

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "beam_discretization_runtime_vtu_output_params.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::BeamRuntimeVtuOutputParams::BeamRuntimeVtuOutputParams()
    : isinit_(false),
      issetup_(false),
      output_displacement_state_(false),
      use_absolute_positions_visualizationpoint_coordinates_(true),
      write_internal_energy_element_(false),
      write_kinetic_energy_element_(false),
      write_triads_visualizationpoints_(false),
      write_material_crosssection_strains_gausspoints_(false),
      write_material_crosssection_stresses_gausspoints_(false),
      write_spatial_crosssection_stresses_gausspoints_(false),
      write_filament_condition_(false),
      write_orientation_parameter_(false),
      write_rve_crosssection_forces_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::Init(
    const Teuchos::ParameterList& IO_vtk_structure_beams_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values

  output_displacement_state_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_beams_paramslist, "DISPLACEMENT");

  use_absolute_positions_visualizationpoint_coordinates_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "USE_ABSOLUTE_POSITIONS");

  write_internal_energy_element_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "INTERNAL_ENERGY_ELEMENT");

  write_kinetic_energy_element_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "KINETIC_ENERGY_ELEMENT");

  write_triads_visualizationpoints_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "TRIAD_VISUALIZATIONPOINT");

  write_material_crosssection_strains_gausspoints_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_beams_paramslist, "STRAINS_GAUSSPOINT");

  write_material_crosssection_stresses_gausspoints_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "MATERIAL_FORCES_GAUSSPOINT");

  write_spatial_crosssection_stresses_gausspoints_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "SPATIAL_FORCES_GAUSSPOINT");

  write_orientation_parameter_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "ORIENTATION_PARAMETER");

  write_rve_crosssection_forces_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_beams_paramslist, "RVE_CROSSSECTION_FORCES");


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
