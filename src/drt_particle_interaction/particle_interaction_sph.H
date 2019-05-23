/*---------------------------------------------------------------------------*/
/*!

\brief smoothed particle hydrodynamics (SPH) interaction handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_H
#define PARTICLE_INTERACTION_SPH_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_base.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
  class SPHDensityBase;
  class SPHPressure;
  class SPHTemperature;
  class SPHMomentum;
  class SPHSurfaceTensionBase;
  class SPHBoundaryParticleBase;
  class SPHPhaseChangeBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class ParticleInteractionSPH : public ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionSPH(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~ParticleInteractionSPH();

    //! init particle interaction handler
    virtual void Init();

    //! setup particle interaction handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface);

    //! write restart of particle interaction handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of particle interaction handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert interaction dependent states of all particle types
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! set initial states
    virtual void SetInitialStates() override;

    //! evaluate particle interactions
    virtual void EvaluateInteractions() override;

    //! maximum interaction distance (on this processor)
    virtual double MaxInteractionDistance() const override;

    //! distribute interaction history
    virtual void DistributeInteractionHistory() const override;

    //! communicate interaction history
    virtual void CommunicateInteractionHistory() const override;

    //! set current time
    virtual void SetCurrentTime(const double currenttime);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize);

   private:
    //! init kernel handler
    virtual void InitKernelHandler();

    //! init equation of state bundle
    virtual void InitEquationOfStateBundle();

    //! init neighbor pair handler
    virtual void InitNeighborPairHandler();

    //! init density handler
    virtual void InitDensityHandler();

    //! init pressure handler
    virtual void InitPressureHandler();

    //! init temperature handler
    virtual void InitTemperatureHandler();

    //! init momentum handler
    virtual void InitMomentumHandler();

    //! init surface tension handler
    virtual void InitSurfaceTensionHandler();

    //! init boundary particle handler
    virtual void InitBoundaryParticleHandler();

    //! init phase change handler
    virtual void InitPhaseChangeHandler();

    //! compute initial consistent particle volume
    virtual double ComputeConsistentParticleVolume() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! density handler
    std::unique_ptr<PARTICLEINTERACTION::SPHDensityBase> density_;

    //! pressure handler
    std::unique_ptr<PARTICLEINTERACTION::SPHPressure> pressure_;

    //! temperature handler
    std::unique_ptr<PARTICLEINTERACTION::SPHTemperature> temperature_;

    //! momentum handler
    std::unique_ptr<PARTICLEINTERACTION::SPHMomentum> momentum_;

    //! surface tension handler
    std::unique_ptr<PARTICLEINTERACTION::SPHSurfaceTensionBase> surfacetension_;

    //! boundary particle handler
    std::unique_ptr<PARTICLEINTERACTION::SPHBoundaryParticleBase> boundaryparticle_;

    //! phase change handler
    std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeBase> phasechange_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
