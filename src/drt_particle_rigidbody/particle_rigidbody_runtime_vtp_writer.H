/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for rigid bodies in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_RUNTIME_VTP_WRITER_H
#define PARTICLE_RIGIDBODY_RUNTIME_VTP_WRITER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "Epetra_Comm.h"
#include <vector>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
class RuntimeVtpWriter;

namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLERIGIDBODY
{
  class RigidBodyDataState;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  /*!
   * \brief rigid body runtime vtp writer class
   *
   * A class that writes visualization output for rigid bodies in vtk/vtp format at runtime.
   *
   * \author Sebastian Fuchs \date 09/2020
   */
  class RigidBodyRuntimeVtpWriter final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] comm communicator
     */
    explicit RigidBodyRuntimeVtpWriter(const Epetra_Comm& comm);

    //! destructor
    ~RigidBodyRuntimeVtpWriter() = default;

    /*!
     * \brief init rigid body runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] rigidbodydatastate rigid body data state container
     */
    void Init(const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate);

    /*!
     * \brief setup rigid body runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] write_binary_output flag for output in binary format
     */
    void Setup(bool write_binary_output);

    /*!
     * \brief read restart of runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief reset current simulation time and time step number
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] time     output time
     * \param[in] timestep output step
     */
    void ResetTimeAndTimeStep(double time, unsigned int timestep);

    /*!
     * \brief set positions and states of rigid bodies
     *
     * Set positions and states of rigid bodies owned by this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] ownedrigidbodies owned rigid bodies by this processor
     */
    void SetRigidBodyPositionsAndStates(const std::vector<int>& ownedrigidbodies);

    /*!
     * \brief write all required vtp files to filesystem
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void WriteFiles();

    /*!
     * \brief write a vtp collection file to filesystem
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void WriteCollectionFileOfAllWrittenFiles();

   private:
    //! communicator
    const Epetra_Comm& comm_;

    //! setup time of runtime vtp writer
    double setuptime_;

    //! field name
    const std::string fieldname_;

    //! rigid body data state container
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate_;

    //! vtp writer object
    std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter_;
  };

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
