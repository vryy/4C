/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for rigid bodies in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_RIGIDBODY_RUNTIME_VTP_WRITER_HPP
#define BACI_PARTICLE_RIGIDBODY_RUNTIME_VTP_WRITER_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_io_visualization_manager.hpp"

#include <Epetra_Comm.h>

#include <memory>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}  // namespace IO

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

    /*!
     * \brief init rigid body runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] rigidbodydatastate rigid body data state container
     */
    void Init(const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate);

    /*!
     * \brief read restart of runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

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
     * \brief Write the visualization files to disk
     */
    void WriteToDisk(const double time, const unsigned int timestep_number);

   private:
    //! communicator
    const Epetra_Comm& comm_;

    //! setup time of runtime vtp writer
    double setuptime_;

    //! rigid body data state container
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate_;

    //! visualization manager
    std::shared_ptr<IO::VisualizationManager> visualization_manager_;
  };

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
