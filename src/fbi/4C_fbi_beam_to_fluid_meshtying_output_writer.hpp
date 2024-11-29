// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_OUTPUT_WRITER_HPP
#define FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_OUTPUT_WRITER_HPP


#include "4C_config.hpp"

#include "4C_io_visualization_parameters.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


// Forward declarations.
namespace Adapter
{
  class FBIConstraintenforcer;
}
namespace BeamInteraction
{
  class BeamToSolidVisualizationOutputWriterBase;
}
namespace FBI
{
  class BeamToFluidMeshtyingVtkOutputParams;
}
namespace Solid::TimeInt
{
  class ParamsRuntimeOutput;
}

namespace BeamInteraction
{
  /**
   * \brief This class manages and creates all visualization output for beam to solid volume
   * meshtying.
   */
  class BeamToFluidMeshtyingVtkOutputWriter
  {
   public:
    /**
     * \brief Constructor
     */
    BeamToFluidMeshtyingVtkOutputWriter(
        const Core::IO::VisualizationParameters& visualization_params,
        std::shared_ptr<const FBI::BeamToFluidMeshtyingVtkOutputParams> output_params_ptr);

    /**
     * \brief empty Destructor
     */
    virtual ~BeamToFluidMeshtyingVtkOutputWriter() = default;

    /**
     * \brief Setup time step output creation, and call WriteOutputData.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator. This is a raw
     * pointer since this function is called from within the sub model evaluator, which does not
     * (and probably can not) have a RCP to itself.
     */
    void write_output_runtime(
        Adapter::FBIConstraintenforcer& couplingenforcer, int i_step, double time) const;

    /**
     * \brief Setup post iteration output creation, and call WriteOutputData.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator. This is a raw
     * pointer since this function is called from within the sub model evaluator, which does not
     * (and probably can not) have a RCP to itself.
     * @param i_iteration (in) current number of iteration.
     */
    void write_output_runtime_iteration(
        const std::shared_ptr<Adapter::FBIConstraintenforcer>& couplingenforcer, int i_iteration,
        int i_step, double time) const;

   private:
    /**
     * \brief Gather all output data after for beam to solid volume mesh tying and write the files
     * to disc.
     * @param beam_contact (in) Pointer to the beam contact sub model evaluator.
     * @param i_step (in) Number of this visualization step (does not have to be continuous, e.g. in
     * iteration visualization).
     * @param time (in) Scalar time value for this visualization step.
     */
    void write_output_beam_to_fluid_mesh_tying(
        Adapter::FBIConstraintenforcer& couplingenforcer, int i_step, double time) const;

   private:
    //! Parameter container for output.
    std::shared_ptr<const FBI::BeamToFluidMeshtyingVtkOutputParams> output_params_ptr_;

    //! Pointer to the output writer, which handles the actual output data for this object.
    std::shared_ptr<BeamInteraction::BeamToSolidVisualizationOutputWriterBase>
        output_writer_base_ptr_;

    //! visualization parameters
    Core::IO::VisualizationParameters visualization_params_;
  };

}  // namespace BeamInteraction

FOUR_C_NAMESPACE_CLOSE

#endif
