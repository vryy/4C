/*-----------------------------------------------------------------------------------------------*/
/*!

\brief input parameters related to VTK output at runtime for the structural (time) integration

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "str_timint_basedataio_runtime_vtk_output.H"
#include "str_discretization_runtime_vtu_output_params.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_beam3/beam_discretization_runtime_vtu_output_params.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
STR::TIMINT::ParamsRuntimeVtkOutput::ParamsRuntimeVtkOutput()
    : isinit_(false),
      issetup_(false),
      output_data_format_(INPAR::IO_RUNTIME_VTK::vague),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_structure_(false),
      output_beams_(false),
      params_runtime_vtu_output_structure_(Teuchos::null),
      params_runtime_vtu_output_beams_(Teuchos::null)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtkOutput::Init(
    const Teuchos::ParameterList& IO_vtk_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values
  output_data_format_ = DRT::INPUT::IntegralValue<INPAR::IO_RUNTIME_VTK::OutputDataFormat>(
      IO_vtk_structure_paramslist, "OUTPUT_DATA_FORMAT");

  output_interval_steps_ = IO_vtk_structure_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtk_structure_paramslist, "EVERY_ITERATION");


  // check for output of structure discretization which is to be handled by an own writer object
  output_structure_ = (bool)DRT::INPUT::IntegralValue<int>(
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
  output_beams_ = (bool)DRT::INPUT::IntegralValue<int>(
      IO_vtk_structure_paramslist.sublist("BEAMS"), "OUTPUT_BEAMS");

  // create and initialize parameter container object for beam specific runtime vtk output
  if (output_beams_)
  {
    params_runtime_vtu_output_beams_ =
        Teuchos::rcp(new DRT::ELEMENTS::BeamRuntimeVtuOutputParams());

    params_runtime_vtu_output_beams_->Init(IO_vtk_structure_paramslist.sublist("BEAMS"));
    params_runtime_vtu_output_beams_->Setup();
  }


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtkOutput::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtkOutput::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
