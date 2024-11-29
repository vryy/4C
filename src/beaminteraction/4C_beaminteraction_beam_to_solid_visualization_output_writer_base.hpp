// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_BASE_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_BASE_HPP


#include "4C_config.hpp"

#include "4C_io_visualization_parameters.hpp"

#include <map>
#include <memory>
#include <string>

FOUR_C_NAMESPACE_OPEN


// Forward declarations.
namespace BeamInteraction
{
  class BeamToSolidOutputWriterVisualization;
}
namespace Solid::TimeInt
{
  class ParamsRuntimeOutput;
}


namespace BeamInteraction
{
  /**
   * \brief A class that stores and manages the output for a visualization in ParaView. This object
   * contains multiple visualization writers, which can be used to store the actual output data in.
   * The way this object is designed it can be passed to a GetVisualization function and all
   * visualization output can be written stored through this object.
   */
  class BeamToSolidVisualizationOutputWriterBase
  {
   public:
    /**
     * \brief Empty constructor, set the class variables.
     * @param base_output_name (in) Base name for the created output files.
     * @param visualization_params (in) visualization parameters
     */
    BeamToSolidVisualizationOutputWriterBase(const std::string& base_output_name,
        Core::IO::VisualizationParameters visualization_params);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVisualizationOutputWriterBase() = default;

    /**
     * \brief Create a new visualization writer in this object.
     * @param writer_name (in) Internal name for the new writer (this one will be used for actual
     * created output files).
     * @param writer_name_key (in) Key for the new writer in the writer map.
     * @return RCP to the newly created writer.
     */
    std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> add_visualization_writer(
        const std::string& writer_name, const std::string& writer_name_key);

    /**
     * \brief Create a new visualization writer in this object.
     * @param writer_name (in) Name of the new writer. If the name already exists, throw an error.
     * @return RCP to the newly created writer.
     */
    std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> add_visualization_writer(
        const std::string& writer_name);

    /**
     * \brief Return the RCP to one of the visualization writers in this object.
     * @param writer_name (in) Name of the writer. If the name does not exist an null pointer will
     * be returned.
     * @return RCP to the writer.
     */
    std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization> get_visualization_writer(
        const std::string& writer_name);

    /**
     * \brief Write all visualization writers to disc. After writing them, the data in the objects
     * will be deleted.
     */
    void write(const unsigned int timestep_number, const double time);

   private:
    //! Base name of the output files create from this object.
    std::string base_output_name_;

    //! Map of the sub output writers.
    std::map<std::string, std::shared_ptr<BeamInteraction::BeamToSolidOutputWriterVisualization>>
        visualization_writers_;

    //! visualization parameters
    const Core::IO::VisualizationParameters visualization_params_;
  };

}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
