/*---------------------------------------------------------------------------*/
/*!
\file particle_algorithm.H

\brief algorithm to control particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ALGORITHM_H
#define PARTICLE_ALGORITHM_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"

#include "../drt_particle_engine/particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class TimInt;
  class GravityHandler;
  class ResultTest;
}  // namespace PARTICLEALGORITHM

namespace PARTICLEINTERACTION
{
  class ParticleInteractionBase;
}

namespace PARTICLEENGINE
{
  class ParticleEngine;
  class ParticleObject;
}  // namespace PARTICLEENGINE

namespace DRT
{
  class ResultTest;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class ParticleAlgorithm : public ADAPTER::AlgorithmBase
  {
   public:
    //! constructor
    explicit ParticleAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~ParticleAlgorithm();

    //! init particle algorithm
    void Init(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles);

    //! setup particle algorithm
    void Setup();

    //! read restart information for given time step
    void ReadRestart(const int restartstep);

    //! time loop for particle problem
    void Timeloop();

    //! prepare time step
    void PrepareTimeStep(bool print_header = true);

    //! integrate particle time step
    void Integrate();

    //! output particle time step
    void Output() const;

    //! create result test
    std::shared_ptr<DRT::ResultTest> CreateResultTest();

   private:
    //! init particle engine
    void InitParticleEngine();

    //! init particle time integration
    void InitParticleTimeIntegration();

    //! init particle interaction handler
    void InitParticleInteraction();

    //! init particle gravity handler
    void InitParticleGravity();

    //! generate initial particles
    void GenerateInitialParticles();

    //! determine all particle types
    void DetermineParticleTypes();

    //! determine particle states of all particle types
    void DetermineParticleStatesOfParticleTypes();

    //! setup initial particles
    void SetupInitialParticles();

    //! setup initial states
    void SetupInitialStates();

    //! check particle types for modified states
    bool HaveModifiedStates();

    //! update connectivity
    void UpdateConnectivity();

    //! check particle transfer including bin size safety checks
    bool CheckParticleTransfer();

    //! check load balance
    bool CheckLoadBalancing();

    //! set initial conditions
    void SetInitialConditions();

    //! set current time
    void SetCurrentTime();

    //! set current step size
    void SetCurrentStepSize();

    //! set gravity acceleration
    void SetGravityAcceleration();

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngine> particleengine_;

    //! particle time integration
    std::unique_ptr<PARTICLEALGORITHM::TimInt> particletimint_;

    //! particle interaction
    std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionBase> particleinteraction_;

    //! particle gravity handler
    std::unique_ptr<PARTICLEALGORITHM::GravityHandler> particlegravity_;

    //! map of particle types and corresponding states
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes_;

    //! vector of initial or generated particles to distribute
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr> particlestodistribute_;

    //! number of particles on this processor after last load balance
    int numparticlesafterlastloadbalance_;

    //! transfer particles to new bins every time step
    bool transferevery_;

    //! write results interval
    const int writeresultsevery_;

    //! write restart interval
    const int writerestartevery_;

    //! result control flag
    bool writeresultsthisstep_;

    //! restart control flag
    bool writerestartthisstep_;

    //! simulation is restarted
    bool isrestarted_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
