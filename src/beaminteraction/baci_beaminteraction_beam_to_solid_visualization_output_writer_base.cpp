/*----------------------------------------------------------------------*/
/*! \file

\brief Base object that stores all relevant data for beam to solid output

\level 3

*/


#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"

#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "baci_utils_exceptions.hpp"

#include <utility>

BACI_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::BeamToSolidVisualizationOutputWriterBase(
    const std::string& base_output_name,
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
    IO::VisualizationParameters visualization_params)
    : base_output_name_(base_output_name),
      visualization_output_params_(visualization_output_params),
      visualization_params_(std::move(visualization_params))
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
  }
  else
  {
    auto new_writer = Teuchos::rcp<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>(
        new BEAMINTERACTION::BeamToSolidOutputWriterVisualization(
            base_output_name_ + "-" + writer_name, visualization_params_,
            visualization_output_params_));
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

BACI_NAMESPACE_CLOSE
