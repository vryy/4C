/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTU output at runtime for beams

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "str_discretization_runtime_vtu_output_params.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::StructureRuntimeVtuOutputParams::StructureRuntimeVtuOutputParams()
    : isinit_(false),
      issetup_(false),
      output_displacement_state_(false),
      output_element_owner_(false),
      output_element_gid_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::StructureRuntimeVtuOutputParams::Init(
    const Teuchos::ParameterList& IO_vtk_structure_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values

  output_displacement_state_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_structure_paramslist, "DISPLACEMENT");

  output_element_owner_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_structure_paramslist, "ELEMENT_OWNER");

  output_element_gid_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_structure_paramslist, "ELEMENT_GID");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::StructureRuntimeVtuOutputParams::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::StructureRuntimeVtuOutputParams::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
