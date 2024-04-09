/*----------------------------------------------------------------------*/
/*! \file

\brief Base object that stores all relevant data for beam to solid output

\level 3

*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_BASE_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_BASE_HPP


#include "baci_config.hpp"

#include "baci_io_visualization_parameters.hpp"

#include <Teuchos_RCP.hpp>

#include <map>
#include <string>

BACI_NAMESPACE_OPEN


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidOutputWriterVisualization;
}
namespace STR::TIMINT
{
  class ParamsRuntimeOutput;
}


namespace BEAMINTERACTION
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
     * @param visualization_output_params (in) RCP to the global visualization parameter list.
     * @param visualization_params (in) visualization parameters
     */
    BeamToSolidVisualizationOutputWriterBase(const std::string& base_output_name,
        Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
        IO::VisualizationParameters visualization_params);

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
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> AddVisualizationWriter(
        const std::string& writer_name, const std::string& writer_name_key);

    /**
     * \brief Create a new visualization writer in this object.
     * @param writer_name (in) Name of the new writer. If the name already exists, throw an error.
     * @return RCP to the newly created writer.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> AddVisualizationWriter(
        const std::string& writer_name);

    /**
     * \brief Return the RCP to one of the visualization writers in this object.
     * @param writer_name (in) Name of the writer. If the name does not exist an null pointer will
     * be returned.
     * @return RCP to the writer.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> GetVisualizationWriter(
        const std::string& writer_name);

    /**
     * \brief Write all visualization writers to disc. After writing them, the data in the objects
     * will be deleted.
     */
    void Write(const unsigned int timestep_number, const double time);

   private:
    //! Base name of the output files create from this object.
    std::string base_output_name_;

    //! Map of the sub output writers.
    std::map<std::string, Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>>
        visualization_writers_;

    //! Pointer to the global visualization input file options.
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params_;

    //! visualization parameters
    const IO::VisualizationParameters visualization_params_;
  };

}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
