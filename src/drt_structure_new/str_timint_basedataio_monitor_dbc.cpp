/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related monitoring reaction forces for the structural (time) integration

\level 3

\maintainer Jonas Eichinger
*/
/*-----------------------------------------------------------------------------------------------*/

#include "str_timint_basedataio_monitor_dbc.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
STR::TIMINT::ParamsMonitorDBC::ParamsMonitorDBC()
    : isinit_(false),
      issetup_(false),
      output_interval_steps_(-1),
      of_precision_(-1),
      os_precision_(-1),
      file_type_("none"),
      write_header_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsMonitorDBC::Init(
    const Teuchos::ParameterList& IO_monitor_dbc_structure_paramslist)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // output intervall in steps
  output_interval_steps_ = IO_monitor_dbc_structure_paramslist.get<int>("INTERVAL_STEPS");

  // file precision
  of_precision_ = IO_monitor_dbc_structure_paramslist.get<int>("PRECISION_FILE");

  // screen precision
  os_precision_ = IO_monitor_dbc_structure_paramslist.get<int>("PRECISION_SCREEN");

  // file type
  file_type_ = IO_monitor_dbc_structure_paramslist.get<std::string>("FILE_TYPE");
  ;

  // write header in csv file
  write_header_ =
      DRT::INPUT::IntegralValue<int>(IO_monitor_dbc_structure_paramslist, "WRITE_HEADER");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsMonitorDBC::Setup()
{
  if (not IsInit()) dserror("Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TIMINT::ParamsMonitorDBC::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
