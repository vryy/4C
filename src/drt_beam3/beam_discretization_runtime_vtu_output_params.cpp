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
      write_triads_visualizationpoints_(false),
      write_material_crosssection_strains_gausspoints_(false),
      write_material_crosssection_stresses_gausspoints_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::Init(
    const Teuchos::ParameterList& IO_vtk_structure_beams_paramslist )
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values

  write_triads_visualizationpoints_ =
      (bool) DRT::INPUT::IntegralValue<int>(
          IO_vtk_structure_beams_paramslist, "TRIAD_VISUALIZATIONPOINT");

  write_material_crosssection_strains_gausspoints_ =
      (bool) DRT::INPUT::IntegralValue<int>(
          IO_vtk_structure_beams_paramslist, "STRAINS_GAUSSPOINT");

  write_material_crosssection_stresses_gausspoints_ =
      (bool) DRT::INPUT::IntegralValue<int>(
          IO_vtk_structure_beams_paramslist, "STRESSES_GAUSSPOINT");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::Setup()
{
  if ( not IsInit() )
    dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::BeamRuntimeVtuOutputParams::CheckInitSetup() const
{
  if ( not IsInit() or not IsSetup() )
    dserror("Call Init() and Setup() first!");
}
