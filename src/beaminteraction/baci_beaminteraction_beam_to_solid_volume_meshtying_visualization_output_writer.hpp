/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid volume meshtying output creation.

\level 3

*/
// End doxygen header.


#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_VISUALIZATION_OUTPUT_WRITER_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VOLUME_MESHTYING_VISUALIZATION_OUTPUT_WRITER_HPP


#include "baci_config.hpp"

#include "baci_io_visualization_parameters.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingVisualizationOutputParams;
  class BeamToSolidVisualizationOutputWriterBase;
  namespace SUBMODELEVALUATOR
  {
    class BeamContact;
  }
}  // namespace BEAMINTERACTION
namespace STR::TIMINT
{
  class ParamsRuntimeOutput;
}


namespace BEAMINTERACTION
{
  /**
   * \brief This class manages and creates all visualization output for beam to solid volume
   * meshtying.
   */
  class BeamToSolidVolumeMeshtyingVisualizationOutputWriter
  {
   public:
    /**
     * \brief Constructor.
     */
    explicit BeamToSolidVolumeMeshtyingVisualizationOutputWriter(
        IO::VisualizationParameters visualization_params);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingVisualizationOutputWriter() = default;

    /**
     * \brief Initialize the object.
     */
    void Init();

    /**
     * \brief Setup the output writer base and the desired field data.
     * @param visualization_output_params (in) RCP to parameter container for global visualization
     * output options.
     * @param output_params_ptr (in) RCP to parameter container for beam to solid output.
     */
    void Setup(Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
        Teuchos::RCP<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputParams>
            output_params_ptr);

    /**
     * \brief Setup time step output creation, and call WriteOutputData.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator. This is a raw
     * pointer since this function is called from within the sub model evaluator, which does not
     * (and probably can not) have a RCP to itself.
     */
    void WriteOutputRuntime(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact) const;

    /**
     * \brief Setup post iteration output creation, and call WriteOutputData.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator. This is a raw
     * pointer since this function is called from within the sub model evaluator, which does not
     * (and probably can not) have a RCP to itself.
     * @param i_iteration (in) current number of iteration.
     */
    void WriteOutputRuntimeIteration(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_iteration) const;

   private:
    /**
     * \brief Gather all output data after for beam to solid volume mesh tying and write the files
     * to disc.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator.
     * @param i_step (in) Number of this visualization step (does not have to be continuous, e.g. in
     * iteration visualization).
     * @param time (in) Scalar time value for this visualization step.
     */
    void WriteOutputBeamToSolidVolumeMeshTying(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_step,
        double time) const;

    /**
     * \brief Checks the init and setup status.
     */
    void CheckInitSetup() const;

    /**
     * \brief Checks the init status.
     */
    void CheckInit() const;

   private:
    //! Flag if object is initialized.
    bool isinit_;

    //! Flag if object is set up.
    bool issetup_;

    //! Parameter container for output.
    Teuchos::RCP<const BeamToSolidVolumeMeshtyingVisualizationOutputParams> output_params_ptr_;

    //! Pointer to the output writer, which handles the actual output data for this object.
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase> output_writer_base_ptr_;

    //! visualization parameters
    const IO::VisualizationParameters visualization_params_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
