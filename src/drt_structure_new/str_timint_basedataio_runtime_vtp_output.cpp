/*-----------------------------------------------------------------------------------------------*/
/*!

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "str_timint_basedataio_runtime_vtp_output.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include "../drt_beam3/beam_discretization_runtime_vtu_output_params.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
STR::TIMINT::ParamsRuntimeVtpOutput::ParamsRuntimeVtpOutput()
    : isinit_(false),
      issetup_(false),
      output_data_format_(INPAR::IO_RUNTIME_VTP_STRUCTURE::vague),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_owner_(false),
      output_orientationandlength_(false),
      output_numberofbonds_(false),
      output_linkingforce_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::Init(
    const Teuchos::ParameterList& IO_vtp_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize the parameter values
  output_data_format_ =
      DRT::INPUT::IntegralValue<INPAR::IO_RUNTIME_VTP_STRUCTURE::OutputDataFormat>(
          IO_vtp_structure_paramslist, "OUTPUT_DATA_FORMAT");

  output_interval_steps_ = IO_vtp_structure_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "EVERY_ITERATION");

  /*  output_displacement_state_ =
        (bool) DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "DISPLACEMENT");*/

  if (output_every_iteration_) dserror("not implemented yet!");

  output_owner_ = (bool)DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "OWNER");

  output_orientationandlength_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "ORIENTATIONANDLENGTH");

  output_numberofbonds_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "NUMBEROFBONDS");

  output_linkingforce_ =
      (bool)DRT::INPUT::IntegralValue<int>(IO_vtp_structure_paramslist, "LINKINGFORCE");


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
