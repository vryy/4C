/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particles in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RUNTIME_VTP_WRITER_H
#define PARTICLE_RUNTIME_VTP_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_container.H"
#include "particle_container_bundle.H"

#include "Epetra_Comm.h"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
class RuntimeVtpWriter;

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief particle runtime vtp writer class
   *
   * A class that writes visualization output for particles in vtk/vtp format at runtime. For each
   * particle type (and particle status) a RuntimeVtpWriter is initialized that writes the output to
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

    //! destructor
    ~ParticleRuntimeVtpWriter() = default;

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
     * \param[in] write_binary_output     flag for output in binary format
     * \param[in] write_ghosted_particles flag for output of ghosted particles
     */
    void Setup(bool write_binary_output, bool write_ghosted_particles);

    /*!
     * \brief read restart of runtime vtp writer
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief reset current simulation time and time step number
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] time     output time
     * \param[in] timestep output step
     */
    void ResetTimeAndTimeStep(double time, unsigned int timestep);

    /*!
     * \brief set positions and states of particles
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void SetParticlePositionsAndStates();

    /*!
     * \brief write all required vtp files to filesystem
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void WriteFiles();

    /*!
     * \brief write a vtp collection file to filesystem
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void WriteCollectionFileOfAllWrittenFiles();

   private:
    //! communicator
    const Epetra_Comm& comm_;

    //! particle container bundle
    ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! collection of vtp writer objects indexed by particle type enum and particle status enum
    std::vector<std::vector<std::shared_ptr<RuntimeVtpWriter>>> runtime_vtpwriters_;

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
    std::set<StateEnum> blackliststates_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
