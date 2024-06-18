/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief input parameters related to VTK output at runtime for the structural (time) integration

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"

#include "4C_beam3_discretization_runtime_output_params.hpp"
#include "4C_global_data.hpp"
#include "4C_structure_new_discretization_runtime_output_params.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TimeInt::ParamsRuntimeOutput::Init(
    const Teuchos::ParameterList& IO_vtk_structure_paramslist)
{
  // We have to call setup() after Init()
  issetup_ = false;

  // initialize the parameter values
  output_interval_steps_ = IO_vtk_structure_paramslist.get<int>("INTERVAL_STEPS");
  output_step_offset_ = IO_vtk_structure_paramslist.get<int>("STEP_OFFSET");
  output_every_iteration_ =
      (bool)Core::UTILS::IntegralValue<int>(IO_vtk_structure_paramslist, "EVERY_ITERATION");

  // check for output of structure discretization which is to be handled by an own writer object
  output_structure_ = (bool)Core::UTILS::IntegralValue<int>(
      IO_vtk_structure_paramslist.sublist("STRUCTURE"), "OUTPUT_STRUCTURE");

  // create and initialize parameter container object for structure specific runtime output
  if (output_structure_)
  {
    params_runtime_output_structure_ =
        Teuchos::rcp(new Discret::ELEMENTS::StructureRuntimeOutputParams());

    params_runtime_output_structure_->Init(IO_vtk_structure_paramslist.sublist("STRUCTURE"));
    params_runtime_output_structure_->setup();
  }


  // check for special beam output which is to be handled by an own writer object
  output_beams_ = (bool)Core::UTILS::IntegralValue<int>(
      IO_vtk_structure_paramslist.sublist("BEAMS"), "OUTPUT_BEAMS");

  // create and initialize parameter container object for beam specific runtime output
  if (output_beams_)
  {
    params_runtime_output_beams_ = Teuchos::rcp(new Discret::ELEMENTS::BeamRuntimeOutputParams());

    params_runtime_output_beams_->Init(IO_vtk_structure_paramslist.sublist("BEAMS"));
    params_runtime_output_beams_->setup();
  }


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TimeInt::ParamsRuntimeOutput::setup()
{
  FOUR_C_ASSERT(is_init(), "Init() has not been called, yet!");

  // Nothing to do here at the moment

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void STR::TimeInt::ParamsRuntimeOutput::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and setup() first!");
}

FOUR_C_NAMESPACE_CLOSE
