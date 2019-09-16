/*---------------------------------------------------------------------------*/
/*! \file
\brief smoothed particle hydrodynamics (SPH) interaction handler

\level 3

\maintainer  Sebastian Fuchs
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
  class ParticleInteractionSPH final : public ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionSPH(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~ParticleInteractionSPH() override;

    //! init particle interaction handler
    void Init() override;

    //! setup particle interaction handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface) override;

    //! write restart of particle interaction handler
    void WriteRestart(const int step, const double time) const override;

    //! read restart of particle interaction handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader) override;

    //! insert interaction dependent states of all particle types
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! set initial states
    void SetInitialStates() override;

    //! evaluate particle interactions
    void EvaluateInteractions() override;

    //! maximum interaction distance (on this processor)
    double MaxInteractionDistance() const override;

    //! distribute interaction history
    void DistributeInteractionHistory() const override;

    //! communicate interaction history
    void CommunicateInteractionHistory() const override;

    //! set current time
    void SetCurrentTime(const double currenttime);

    //! set current step size
    void SetCurrentStepSize(const double currentstepsize);

   private:
    //! init kernel handler
    void InitKernelHandler();

    //! init equation of state bundle
    void InitEquationOfStateBundle();

    //! init neighbor pair handler
    void InitNeighborPairHandler();

    //! init density handler
    void InitDensityHandler();

    //! init pressure handler
    void InitPressureHandler();

    //! init temperature handler
    void InitTemperatureHandler();

    //! init momentum handler
    void InitMomentumHandler();

    //! init surface tension handler
    void InitSurfaceTensionHandler();

    //! init boundary particle handler
    void InitBoundaryParticleHandler();

    //! init phase change handler
    void InitPhaseChangeHandler();

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
