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

#include "../drt_particle_engine/particle_typedefs.H"

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

namespace PARTICLERIGIDBODY
{
  class RigidBodyDataState;
  class RigidBodyRuntimeVtpWriter;
  class RigidBodyAffiliationPairs;
}  // namespace PARTICLERIGIDBODY

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

    /*!
     * \brief insert rigid body handler dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    /*!
     * \brief write rigid body runtime output
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] step output step
     * \param[in] time output time
     */
    void WriteRigidBodyRuntimeOutput(const int step, const double time) const;

    /*!
     * \brief set unique global ids for all rigid bodies
     *
     * Currently, the rigid body color of each rigid particle defined via the input file in the
     * particle section defines the affiliation of a rigid particle to a rigid body. First, the
     * affiliation pair data is initialized. Then, the rigid body unique global identifier handler
     * is initialized. This requires knowledge of all global ids currently assigned to rigid bodies
     * over all processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetUniqueGlobalIdsForAllRigidBodies();

    /*!
     * \brief allocate rigid body states
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void AllocateRigidBodyStates();

    /*!
     * \brief distribute rigid body
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void DistributeRigidBody();

    /*!
     * \brief communicate rigid body
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void CommunicateRigidBody();

    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> GetRigidBodyDataState() const override
    {
      return rigidbodydatastate_;
    }

    const std::vector<int>& GetOwnedRigidBodies() const override { return ownedrigidbodies_; }

   private:
    /*!
     * \brief init rigid body unique global identifier handler
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void InitRigidBodyUniqueGlobalIdHandler();

    /*!
     * \brief init rigid body data state container
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void InitRigidBodyDataState();

    /*!
     * \brief init rigid body runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void InitRigidBodyVtpWriter();

    /*!
     * \brief init affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void InitAffiliationPairHandler();

    /*!
     * \brief setup rigid body runtime vtp writer
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetupRigidBodyVtpWriter();

    /*!
     * \brief get packed rigid body state data
     *
     * Get the packed required states of all rigid bodies owned by this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[out] buffer buffer of packed rigid body state data
     */
    void GetPackedRigidBodyStates(std::vector<char>& buffer) const;

    /*!
     * \brief extract packed rigid body state data
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] buffer buffer of packed rigid body state data
     */
    void ExtractPackedRigidBodyStates(std::vector<char>& buffer);

    /*!
     * \brief update rigid body ownership
     *
     * Update the ownership of rigid bodies on all processors. First, the global ids of rigid bodies
     * hosted, i.e., owned and non-owned, by this processor are determined using the affiliation of
     * rigid particles to rigid bodies. Then, based on the maximum number of rigid particles per
     * processor the owning processor of all rigid bodies is determined. Finally, each processor
     * stores the global ids of rigid bodies owned.
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void UpdateRigidBodyOwnership();

    /*!
     * \brief relate owned rigid bodies to all hosting processors
     *
     * Relate the global ids of rigid bodies owned by this processor to the corresponding hosting
     * processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void RelateOwnedRigidBodiesToHostingProcs();

    /*!
     * \brief communicate rigid body states
     *
     * Communicate required states of all rigid bodies previously owned by this processor to the new
     * owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] previouslyownedrigidbodies rigid bodies previously owned by this processor
     */
    void CommunicateRigidBodyStates(std::vector<int>& previouslyownedrigidbodies);

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

    //! rigid body data state container
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> rigidbodydatastate_;

    //! rigid body runtime vtp writer
    std::unique_ptr<PARTICLERIGIDBODY::RigidBodyRuntimeVtpWriter> rigidbodyvtpwriter_;

    //! affiliation pair handler
    std::unique_ptr<PARTICLERIGIDBODY::RigidBodyAffiliationPairs> affiliationpairs_;

    //! rigid bodies owned by this processor
    std::vector<int> ownedrigidbodies_;

    //! rigid bodies hosted (owned and non-owned) by this processor
    std::vector<int> hostedrigidbodies_;

    //! owner of all rigid bodies
    std::vector<int> ownerofrigidbodies_;

    //! owned rigid bodies related to hosting processors (without this processor)
    std::vector<std::vector<int>> ownedrigidbodiestohostingprocs_;
  };
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
