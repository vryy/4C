/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTP output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_structure_new_timint_basedataio_runtime_vtp_output.hpp"

#include "baci_beam3_discretization_runtime_output_params.hpp"
#include "baci_global_data.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsRuntimeVtpOutput::Init(
    const Teuchos::ParameterList& IO_vtp_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  output_interval_steps_ = IO_vtp_structure_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool)CORE::UTILS::IntegralValue<int>(IO_vtp_structure_paramslist, "EVERY_ITERATION");

  if (output_every_iteration_)
    dserror("Every iteration output not implemented for structure vtp output!");

  output_owner_ = (bool)CORE::UTILS::IntegralValue<int>(IO_vtp_structure_paramslist, "OWNER");

  output_orientationandlength_ =
      (bool)CORE::UTILS::IntegralValue<int>(IO_vtp_structure_paramslist, "ORIENTATIONANDLENGTH");

  output_numberofbonds_ =
      (bool)CORE::UTILS::IntegralValue<int>(IO_vtp_structure_paramslist, "NUMBEROFBONDS");

  output_linkingforce_ =
      (bool)CORE::UTILS::IntegralValue<int>(IO_vtp_structure_paramslist, "LINKINGFORCE");


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

FOUR_C_NAMESPACE_CLOSE
