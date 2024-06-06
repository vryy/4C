/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particles in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ENGINE_RUNTIME_VTP_WRITER_HPP
#define FOUR_C_PARTICLE_ENGINE_RUNTIME_VTP_WRITER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_io_visualization_manager.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::IO
{
  class DiscretizationReader;
}  // namespace Core::IO

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief particle runtime vtp writer class
   *
   * A class that writes visualization output for particles in vtk/vtu format at runtime. For each
   * particle type (and particle status) a RuntimeVtuWriter is initialized that writes the output to
   * a separate name.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  class ParticleRuntimeVtpWriter final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] comm communicator
     */
    explicit ParticleRuntimeVtpWriter(const Epetra_Comm& comm);

    /*!
     * \brief init particle runtime vtp writer
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] particlecontainerbundle particle container bundle
     */
    void Init(const ParticleContainerBundleShrdPtr particlecontainerbundle);

    /*!
     * \brief setup particle runtime vtp writer
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] write_ghosted_particles flag for output of ghosted particles
     */
    void Setup(bool write_ghosted_particles);

    /*!
     * \brief read restart of runtime vtp writer
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] reader discretization reader
     */
    void read_restart(const std::shared_ptr<Core::IO::DiscretizationReader> reader);

    /*!
     * \brief set positions and states of particles
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void set_particle_positions_and_states();

    /*!
     * \brief Write visualization files to disk
     */
    void WriteToDisk(const double visualization_time, const int visualization_step);

   private:
    //! communicator
    const Epetra_Comm& comm_;

    //! particle container bundle
    ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! collection of vtu writer objects indexed by particle type enum and particle status enum
    std::vector<std::vector<std::shared_ptr<Core::IO::VisualizationManager>>>
        runtime_visualization_managers_;

    //! setup time of runtime vtp writer
    double setuptime_;

    /*!
     * \brief black list of particle state enums
     *
     * All particle states on the black list are ignored when writing output. The states on the
     * black list are hard coded and typically contain states not relevant for output.
     *
     * \note Reconsider the concept of a black list with hard coded particle state enums. Think
     *       about defining all particle states for the output via a list in the input file.
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    std::set<ParticleState> blackliststates_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
