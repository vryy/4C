/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_structure_new_timint_basedataio_runtime_vtp_output.H"

#include "baci_beam3_discretization_runtime_output_params.H"
#include "baci_global_data.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::Init(
    const Teuchos::ParameterList& IO_vtp_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values
  visualization_parameters_ = IO::VisualizationParametersFactory(
      GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"));

  output_interval_steps_ = IO_vtp_structure_paramslist.get<int>("INTERVAL_STEPS");

  output_step_offset_ = IO_vtp_structure_paramslist.get<int>("STEP_OFFSET");

  output_every_iteration_ =
      (bool)INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "EVERY_ITERATION");

  // Overwrite non default values in the visualization parameters
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  if (output_every_iteration_)
    dserror("Every iteration output not implemented for structure vtp output!");

  /*  output_displacement_state_ =
        (bool) INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "DISPLACEMENT");*/

  output_owner_ = (bool)INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "OWNER");

  output_orientationandlength_ =
      (bool)INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "ORIENTATIONANDLENGTH");

  output_numberofbonds_ =
      (bool)INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "NUMBEROFBONDS");

  output_linkingforce_ =
      (bool)INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "LINKINGFORCE");


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::Setup()
{
  dsassert(IsInit(), "Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::CheckInitSetup() const
{
  dsassert(IsInit() and IsSetup(), "Call Init() and Setup() first!");
}

BACI_NAMESPACE_CLOSE
