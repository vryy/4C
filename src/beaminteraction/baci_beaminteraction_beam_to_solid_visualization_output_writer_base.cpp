/*----------------------------------------------------------------------*/
/*! \file

\brief Base object that stores all relevant data for beam to solid output

\level 3

*/


#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_base.H"

#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.H"
#include "baci_io_control.H"
#include "baci_lib_globalproblem.H"
#include "baci_utils_exceptions.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::BeamToSolidVisualizationOutputWriterBase(
    const std::string& base_output_name,
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
    double restart_time)
    : base_output_name_(base_output_name),
      visualization_output_params_(visualization_output_params),
      restart_time_(restart_time)
{
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::AddVisualizationWriter(
    const std::string& writer_name, const std::string& writer_name_key)
{
  const auto& it = visualization_writers_.find(writer_name_key);
  if (it != visualization_writers_.end())
  {
    dserror("The output writer key '%s' you want to add already exists.", writer_name_key.c_str());
    return Teuchos::null;
  }
  else
  {
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> new_writer =
        Teuchos::rcp<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>(
            new BEAMINTERACTION::BeamToSolidOutputWriterVisualization(
                base_output_name_ + "-" + writer_name, visualization_output_params_,
                restart_time_));
    visualization_writers_[writer_name_key] = new_writer;
    return new_writer;
  }
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::AddVisualizationWriter(
    const std::string& writer_name)
{
  return AddVisualizationWriter(writer_name, writer_name);
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::GetVisualizationWriter(
    const std::string& writer_name)
{
  const auto& it = visualization_writers_.find(writer_name);
  if (it != visualization_writers_.end())
    return it->second;
  else
    return Teuchos::null;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::Write(
    const unsigned int timestep_number, const double time)
{
  for (auto& it : visualization_writers_) it.second->Write(timestep_number, time);
}
