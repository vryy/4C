/*---------------------------------------------------------------------------*/
/*!
\file particle_runtime_vtp_writer.H

\brief write visualization output for particles in vtk/vtp format at runtime

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RUNTIME_VTP_WRITER_H
#define PARTICLE_RUNTIME_VTP_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_container.H"
#include "particle_container_bundle.H"

#include "Epetra_Comm.h"
#include "../drt_io/io_control.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
class RuntimeVtpWriter;

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleRuntimeVtpWriter
  {
   public:
    //! constructor
    explicit ParticleRuntimeVtpWriter(const Epetra_Comm& comm);

    //! destructor
    ~ParticleRuntimeVtpWriter();

    //! init particle runtime vtp writer
    void Init(const ParticleContainerBundleShrdPtr particlecontainerbundle);

    //! setup particle runtime vtp writer
    void Setup(bool write_binary_output, bool write_ghosted_particles);

    //! write restart of runtime vtp writer
    void WriteRestart(const int step, const double time) const;

    //! read restart of runtime vtp writer
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! reset current simulation time and time step number
    void ResetTimeAndTimeStep(double time, unsigned int timestep);

    //! set positions and states of particles
    void SetParticlePositionsAndStates();

    //! write all required vtp files to filesystem
    void WriteFiles();

    //! write a vtp collection file to filesystem
    void WriteCollectionFileOfAllWrittenFiles();

   private:
    //! communication
    const Epetra_Comm& comm_;

    //! particle container bundle
    ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! holds all vtp writer objects for each particle type and status
    std::map<TypeEnum, std::map<StatusEnum, std::shared_ptr<RuntimeVtpWriter>>> runtime_vtpwriters_;

    //! setup time of runtime vtp writer
    double setuptime_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
