/*----------------------------------------------------------------------*/
/*! \file

\brief Base object that stores all relevant data for beam to solid output

\level 3

*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VTU_OUTPUT_WRITER_BASE_H
#define BEAM_TO_SOLID_VTU_OUTPUT_WRITER_BASE_H


#include <string>
#include <map>
#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidVtuOutputWriterVisualization;
}
namespace STR
{
  namespace TIMINT
  {
    class ParamsRuntimeVtkOutput;
  }
}  // namespace STR


namespace BEAMINTERACTION
{
  /**
   * \brief A class that stores and manages the output for a visualization in ParaView. This object
   * contains multiple visualization writers, which can be used to store the actual output data in.
   * The way this object is designed it can be passed to a GetVisualization function and all
   * visualization output can be written stored through this object.
   */
  class BeamToSolidVtuOutputWriterBase
  {
   public:
    /**
     * \brief Empty constructor, set the class variables.
     * @param base_output_name (in) Base name for the created output files.
     * @param vtk_params (in) RCP to the global vtk parameter list.
     * @param restart_time (in) Time at creation of this object.
     */
    BeamToSolidVtuOutputWriterBase(const std::string& base_output_name,
        Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params, double restart_time);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVtuOutputWriterBase(){};

    /**
     * \brief Create a new visualization writer in this object.
     * @param writer_name (in) Internal name for the new writer (this one will be used for actual
     * created output files).
     * @param writer_name_key (in) Key for the new writer in the writer map.
     * @return RCP to the newly created writer.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> AddVisualizationWriter(
        const std::string& writer_name, const std::string& writer_name_key);

    /**
     * \brief Create a new visualization writer in this object.
     * @param writer_name (in) Name of the new writer. If the name already exists, throw an error.
     * @return RCP to the newly created writer.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> AddVisualizationWriter(
        const std::string& writer_name);

    /**
     * \brief Return the RCP to one of the visualization writers in this object.
     * @param writer_name (in) Name of the writer. If the name does not exist an null pointer will
     * be returned.
     * @return RCP to the writer.
     */
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> GetVisualizationWriter(
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
    std::map<std::string, Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>>
        visualization_writers_;

    //! Pointer to the global vtk input file options.
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params_;

    //! Time at which a restart simulation is started.
    const double restart_time_;
  };

}  // namespace BEAMINTERACTION

#endif
