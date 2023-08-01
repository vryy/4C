/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTU output at runtime for beams

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_structure_new_discretization_runtime_vtu_output_params.H"

#include "baci_utils_exceptions.H"

#include "baci_inpar_parameterlist_utils.H"
#include "baci_inpar_structure.H"
#include "baci_lib_globalproblem.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DRT::ELEMENTS::StructureRuntimeVtuOutputParams::StructureRuntimeVtuOutputParams()
    : isinit_(false),
      issetup_(false),
      output_displacement_state_(false),
      output_velocity_state_(false),
      output_element_owner_(false),
      output_element_gid_(false),
      output_element_ghosting_(false),
      output_node_gid_(false),
      output_stress_strain_(false),
      gauss_point_data_output_type_(INPAR::STR::GaussPointDataOutputType::none)
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
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "DISPLACEMENT");
  output_velocity_state_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "VELOCITY");
  output_element_owner_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "ELEMENT_OWNER");
  output_element_gid_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "ELEMENT_GID");
  output_element_ghosting_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "ELEMENT_GHOSTING");
  output_node_gid_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "NODE_GID");
  output_stress_strain_ =
      DRT::INPUT::IntegralValue<bool>(IO_vtk_structure_structure_paramslist, "STRESS_STRAIN");
  gauss_point_data_output_type_ = Teuchos::getIntegralValue<INPAR::STR::GaussPointDataOutputType>(
      IO_vtk_structure_structure_paramslist, "GAUSS_POINT_DATA_OUTPUT_TYPE");

  if (output_stress_strain_)
  {
    // If stress / strain data should be output, check that the relevant parameters in the --IO
    // section are set.
    const Teuchos::ParameterList& io_parameter_list = DRT::Problem::Instance()->IOParams();
    INPAR::STR::StressType io_stress =
        DRT::INPUT::IntegralValue<INPAR::STR::StressType>(io_parameter_list, "STRUCT_STRESS");
    INPAR::STR::StrainType io_strain =
        DRT::INPUT::IntegralValue<INPAR::STR::StrainType>(io_parameter_list, "STRUCT_STRAIN");
    if (io_stress == INPAR::STR::stress_none and io_strain == INPAR::STR::strain_none)
      dserror(
          "If stress / strain runtime output is required, one or two of the flags STRUCT_STRAIN / "
          "STRUCT_STRESS in the --IO section has to be set.");
  }

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
