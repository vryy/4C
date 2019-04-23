/*!
\file beam_to_solid_volume_meshtying_vtk_output_writer.H

\brief Object to handle beam to solid volume meshtying output creation.

\level 3

\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_VTK_OUTPUT_WRITER_H
#define BEAM_TO_SOLID_VOLUME_MESHTYING_VTK_OUTPUT_WRITER_H


#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace BEAMINTERACTION
{
  class BeamToSolidVolumeMeshtyingVtkOutputParams;
  class BeamToSolidVtuOutputWriterBase;
  namespace SUBMODELEVALUATOR
  {
    class BeamContact;
  }
}  // namespace BEAMINTERACTION
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
   * \brief This class manages and creates all visualization output for beam to solid volume
   * meshtying.
   */
  class BeamToSolidVolumeMeshtyingVtkOutputWriter
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToSolidVolumeMeshtyingVtkOutputWriter();

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVolumeMeshtyingVtkOutputWriter(){};

    /**
     * \brief Initialize the object.
     */
    void Init();

    /**
     * \brief Setup the output writer base and the desired field data.
     * @param vtk_params (in) RCP to parameter container for global vtk options.
     * @param output_params_ptr (in) RCP to parameter container for beam to solid output.
     * @param restart_time (in) This is needed in case the simulation is a restart.
     */
    void Setup(Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params,
        Teuchos::RCP<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams>
            output_params_ptr,
        double restart_time);

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
    Teuchos::RCP<const BeamToSolidVolumeMeshtyingVtkOutputParams> output_params_ptr_;

    //! Pointer to the output writer, which handles the actual output data for this object.
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterBase> output_writer_base_ptr_;
  };

}  // namespace BEAMINTERACTION

#endif
