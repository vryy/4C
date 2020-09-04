/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_H
#define PARTICLE_RIGIDBODY_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_rigidbody_interface.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class UniqueGlobalIdHandler;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  /*!
   * \brief rigid body handler for particle problem
   *
   * \author Sebastian Fuchs \date 08/2020
   */
  class RigidBodyHandler final : public RigidBodyHandlerInterface
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit RigidBodyHandler(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~RigidBodyHandler() override;

    /*!
     * \brief init rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Init();

    /*!
     * \brief setup rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief write restart of rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void WriteRestart() const;

    /*!
     * \brief read restart of rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

   private:
    /*!
     * \brief init rigid body unique global identifier handler
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void InitRigidBodyUniqueGlobalIdHandler();

    //! communicator
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! rigid body unique global identifier handler
    std::unique_ptr<PARTICLEENGINE::UniqueGlobalIdHandler> rigidbodyuniqueglobalidhandler_;
  };
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
