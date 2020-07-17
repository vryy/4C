/*---------------------------------------------------------------------------*/
/*! \file
\brief smoothed particle hydrodynamics (SPH) interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_H
#define PARTICLE_INTERACTION_SPH_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_base.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
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
  class SPHOpenBoundaryBase;
  class SPHVirtualWallParticle;
  class SPHPhaseChangeBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  /*!
   * \brief smoothed particle hydrodynamics (SPH) interaction
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class ParticleInteractionSPH final : public ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionSPH(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
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

    //! prepare time step
    void PrepareTimeStep() override;

    //! evaluate particle interactions
    void EvaluateInteractions() override;

    //! post evaluate time step
    void PostEvaluateTimeStep() override;

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

    //! init dirichlet open boundary handler
    void InitDirichletOpenBoundaryHandler();

    //! init neumann open boundary handler
    void InitNeumannOpenBoundaryHandler();

    //! init virtual wall particle handler
    void InitVirtualWallParticleHandler();

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

    //! dirichlet open boundary handler
    std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase> dirichletopenboundary_;

    //! neumann open boundary handler
    std::unique_ptr<PARTICLEINTERACTION::SPHOpenBoundaryBase> neumannopenboundary_;

    //! virtual wall particle handler
    std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle_;

    //! phase change handler
    std::unique_ptr<PARTICLEINTERACTION::SPHPhaseChangeBase> phasechange_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
