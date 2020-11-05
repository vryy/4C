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
#include <unordered_map>

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
     * \brief set initial affiliation pair data
     *
     * The rigid body color of each rigid particle given via the input file in the particle section
     * defines the affiliation of a rigid particle to a rigid body. Upon this information the
     * affiliation pair data is initialized.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetInitialAffiliationPairData();

    /*!
     * \brief set unique global ids for all rigid bodies
     *
     * The rigid body unique global identifier handler is initialized. This requires knowledge of
     * all global ids currently assigned to rigid bodies over all processors.
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
     * \brief set initial states
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetInitialStates();

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

    /*!
     * \brief compute mass quantities of rigid bodies
     *
     * Compute the mass quantities, i.e., mass, center of gravity, and mass moment of inertia, of
     * rigid bodies over all processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ComputeRigidBodyMassQuantities();

    /*!
     * \brief clear forces and torques
     *
     * Clear the forces and torques acting on rigid bodies and on the corresponding rigid particles
     * affiliated with the rigid bodies.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearForcesAndTorques();

    /*!
     * \brief add gravity acceleration
     *
     * Add the gravity acceleration acting on rigid bodies.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] gravity gravity acceleration
     */
    void AddGravityAcceleration(std::vector<double>& gravity);

    /*!
     * \brief compute accelerations of rigid bodies
     *
     * First compute the force and the torque acting on rigid bodies over all processors. Then, the
     * accelerations are computed from force and torque using the mass quantities on the owning
     * processor. Finally, the corresponding accelerations of rigid particles affiliated with the
     * rigid bodies are set.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ComputeAccelerations();

    /*!
     * \brief update positions with given time increment
     *
     * Update the positions of rigid bodies with a given time increment and set the corresponding
     * positions of rigid particles affiliated with the rigid bodies.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    void UpdatePositions(const double timeincrement) override;

    /*!
     * \brief update velocities with given time increment
     *
     * Update the velocities of rigid bodies with a given time increment and set the corresponding
     * velocities of rigid particles affiliated with the rigid bodies.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    void UpdateVelocities(const double timeincrement) override;

    /*!
     * \brief clear accelerations
     *
     * Clear the accelerations of rigid bodies. However, note that the accelerations of rigid
     * particles remain untouched!
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearAccelerations() override;

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

    /*!
     * \brief clear partial mass quantities of rigid bodies
     *
     * Clear partial mass quantities of hosted rigid bodies on this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearPartialMassQuantities();

    /*!
     * \brief compute partial mass quantities of rigid bodies
     *
     * Compute the partial mass quantities of hosted rigid bodies on respective hosting processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ComputePartialMassQuantities();

    /*!
     * \brief gather partial mass quantities of rigid bodies
     *
     * Gather the partial mass quantities for owned rigid bodies from respective hosting processors
     * on owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[out] gatheredpartialmass     gathered partial mass
     * \param[out] gatheredpartialinertia  gathered partial inertia
     * \param[out] gatheredpartialposition gathered partial position
     */
    void GatherPartialMassQuantities(
        std::unordered_map<int, std::vector<double>>& gatheredpartialmass,
        std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialinertia,
        std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialposition);

    /*!
     * \brief compute full mass quantities of rigid bodies
     *
     * Compute the full mass quantities of owned rigid bodies on owning processor using all gathered
     * partial mass quantities.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] gatheredpartialmass     gathered partial mass
     * \param[in] gatheredpartialinertia  gathered partial inertia
     * \param[in] gatheredpartialposition gathered partial position
     */
    void ComputeFullMassQuantities(
        std::unordered_map<int, std::vector<double>>& gatheredpartialmass,
        std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialinertia,
        std::unordered_map<int, std::vector<std::vector<double>>>& gatheredpartialposition);

    /*!
     * \brief clear force and torque acting on rigid bodies
     *
     * Clear the force and torque of hosted rigid bodies on this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearRigidBodyForceAndTorque();

    /*!
     * \brief clear force acting on rigid particles
     *
     * Clear the force of rigid particles owned on this processor.
     */
    void ClearRigidParticleForce();

    /*!
     * \brief compute partial force and torque acting on rigid bodies
     *
     * Compute the partial force and torque acting on hosted rigid bodies on respective hosting
     * processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ComputePartialForceAndTorque();

    /*!
     * \brief gather partial and compute full force and torque acting on rigid bodies
     *
     * Gather the partial force and torque acting on owned rigid bodies from respective hosting
     * processors on owning processor and compute the full force and torque.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void GatherPartialAndComputeFullForceAndTorque();

    /*!
     * \brief compute accelerations of rigid bodies from force and torque
     *
     * Compute the accelerations of owned rigid bodies from force and torque using mass quantities
     * on the owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ComputeAccelerationsFromForceAndTorque();

    /*!
     * \brief clear orientation of rigid bodies
     *
     * Clear the orientation of rigid bodies owned on this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearRigidBodyOrientation();

    /*!
     * \brief update positions of rigid bodies with given time increment
     *
     * Update the positions of rigid bodies with a given time increment on the owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    void UpdateRigidBodyPositions(const double timeincrement);

    /*!
     * \brief update velocities of rigid bodies with given time increment
     *
     * Update the velocities of rigid bodies with a given time increment on the owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    void UpdateRigidBodyVelocities(const double timeincrement);

    /*!
     * \brief clear accelerations of rigid bodies
     *
     * Clear the accelerations of rigid bodies on the owning processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void ClearRigidBodyAccelerations();

    /*!
     * \brief broadcast positions of rigid bodies
     *
     * Broadcast positions of rigid bodies from the owning processor to respective local processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void BroadcastRigidBodyPositions();

    /*!
     * \brief broadcast velocities of rigid bodies
     *
     * Broadcast velocities of rigid bodies from the owning processor to respective local
     * processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void BroadcastRigidBodyVelocities();

    /*!
     * \brief broadcast accelerations of rigid bodies
     *
     * Broadcast accelerations of rigid bodies from the owning processor to respective local
     * processors.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void BroadcastRigidBodyAccelerations();

    /*!
     * \brief set relative position of rigid particles in reference frame
     *
     * Set the relative position of rigid particles in the reference frame owned on this processor.
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void SetRigidParticleReferenceRelativePosition();

    /*!
     * \brief update relative position of rigid particles
     *
     * Update the relative position of rigid particles in the rotating frame owned on this
     * processor.
     */
    void UpdateRigidParticleRelativePosition();

    /*!
     * \brief set position of rigid particles
     *
     * Set the position of rigid particles owned on this processor consistent to the affiliated
     * rigid body.
     */
    void SetRigidParticlePosition();

    /*!
     * \brief set velocities of rigid particles
     *
     * Set the velocities of rigid particles owned on this processor consistent to the affiliated
     * rigid body.
     */
    void SetRigidParticleVelocities();

    /*!
     * \brief set accelerations of rigid particles
     *
     * Set the accelerations of rigid particles owned on this processor consistent to the affiliated
     * rigid body.
     */
    void SetRigidParticleAccelerations();

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
