/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTK output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_structure_new_timint_basedataio_runtime_vtk_output.H"

#include "baci_beam3_discretization_runtime_output_params.H"
#include "baci_global_data.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_structure_new_discretization_runtime_vtu_output_params.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeOutput::Init(
    const Teuchos::ParameterList& IO_vtk_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // Set general output parameters
  visualization_parameters_ = IO::VisualizationParametersFactory(
      DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"));

  // initialize the parameter values
  output_interval_steps_ = IO_vtk_structure_paramslist.get<int>("INTERVAL_STEPS");
  output_step_offset_ = IO_vtk_structure_paramslist.get<int>("STEP_OFFSET");
  output_every_iteration_ =
      (bool)INPUT::IntegralValue<int>(IO_vtk_structure_paramslist, "EVERY_ITERATION");

  // Overwrite non default values in the visualization parameters
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  // check for output of structure discretization which is to be handled by an own writer object
  output_structure_ = (bool)INPUT::IntegralValue<int>(
      IO_vtk_structure_paramslist.sublist("STRUCTURE"), "OUTPUT_STRUCTURE");

  // create and initialize parameter container object for structure specific runtime vtk output
  if (output_structure_)
  {
    params_runtime_vtu_output_structure_ =
        Teuchos::rcp(new DRT::ELEMENTS::StructureRuntimeVtuOutputParams());

    params_runtime_vtu_output_structure_->Init(IO_vtk_structure_paramslist.sublist("STRUCTURE"));
    params_runtime_vtu_output_structure_->Setup();
  }


  // check for special beam output which is to be handled by an own writer object
  output_beams_ =
      (bool)INPUT::IntegralValue<int>(IO_vtk_structure_paramslist.sublist("BEAMS"), "OUTPUT_BEAMS");

  // create and initialize parameter container object for beam specific runtime vtk output
  if (output_beams_)
  {
    params_runtime_vtu_output_beams_ = Teuchos::rcp(new DRT::ELEMENTS::BeamRuntimeOutputParams());

    params_runtime_vtu_output_beams_->Init(IO_vtk_structure_paramslist.sublist("BEAMS"));
    params_runtime_vtu_output_beams_->Setup();
  }


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeOutput::Setup()
{
  dsassert(IsInit(), "Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeOutput::CheckInitSetup() const
{
  dsassert(IsInit() and IsSetup(), "Call Init() and Setup() first!");
}

BACI_NAMESPACE_CLOSE
