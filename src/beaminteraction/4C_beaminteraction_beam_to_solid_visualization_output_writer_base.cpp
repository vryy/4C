/*----------------------------------------------------------------------*/
/*! \file

\brief Base object that stores all relevant data for beam to solid output

\level 3

*/


#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"

#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::BeamToSolidVisualizationOutputWriterBase(
    const std::string& base_output_name,
    Teuchos::RCP<const Solid::TimeInt::ParamsRuntimeOutput> visualization_output_params,
    Core::IO::VisualizationParameters visualization_params)
    : base_output_name_(base_output_name),
      visualization_output_params_(visualization_output_params),
      visualization_params_(std::move(visualization_params))
{
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::add_visualization_writer(
    const std::string& writer_name, const std::string& writer_name_key)
{
  const auto& it = visualization_writers_.find(writer_name_key);
  if (it != visualization_writers_.end())
  {
    FOUR_C_THROW(
        "The output writer key '%s' you want to add already exists.", writer_name_key.c_str());
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
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::add_visualization_writer(
    const std::string& writer_name)
{
  return add_visualization_writer(writer_name, writer_name);
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::get_visualization_writer(
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
void BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase::write(
    const unsigned int timestep_number, const double time)
{
  for (auto& it : visualization_writers_) it.second->write(timestep_number, time);
}

FOUR_C_NAMESPACE_CLOSE
