/*---------------------------------------------------------------------------*/
/*! \file
\brief density handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_DENSITY_H
#define PARTICLE_INTERACTION_SPH_DENSITY_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

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
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
  class MaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
  class SPHDensityCorrectionBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHDensityBase
  {
   public:
    //! constructor
    explicit SPHDensityBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHDensityBase(){};

    //! init density handler
    virtual void Init();

    //! setup density handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of density handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of density handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize) final;

    //! insert density evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const = 0;

    //! compute density field
    virtual void ComputeDensity() const = 0;

   protected:
    //! evaluate sum of weighted mass and colorfield
    virtual void SumWeightedMassAndColorfield() const final;

    //! clear density sum and colorfield states
    virtual void ClearDensitySumAndColorfieldStates() const final;

    //! sum weighted mass/colorfield (self contribution)
    virtual void SumWeightedMassAndColorfieldSelfContribution() const final;

    //! sum weighted mass/colorfield (particle contribution)
    virtual void SumWeightedMassAndColorfieldParticleContribution() const final;

    //! evaluate continuity equation
    virtual void ContinuityEquation() const final;

    //! clear density dot state
    virtual void ClearDensityDotState() const final;

    //! continuity equation (particle contribution)
    virtual void ContinuityEquationParticleContribution() const final;

    //! set density sum to density field
    virtual void SetDensitySum() const final;

    //! add time step scaled density dot to density field
    virtual void AddTimeStepScaledDensityDot() const final;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! density of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh densitytorefresh_;

    //! time step size
    double dt_;

    //! compute colorfield
    bool computecolorfield_;

    //! apply transport velocity formulation
    bool applytransportvelocity_;
  };

  class SPHDensitySummation final : public SPHDensityBase
  {
   public:
    //! constructor
    explicit SPHDensitySummation(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHDensitySummation() override{};

    //! insert density evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const override;

    //! compute density field
    void ComputeDensity() const override;
  };

  class SPHDensityIntegration final : public SPHDensityBase
  {
   public:
    //! constructor
    explicit SPHDensityIntegration(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHDensityIntegration() override{};

    //! insert density evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const override;

    //! compute density field
    void ComputeDensity() const override;
  };

  class SPHDensityPredictCorrect final : public SPHDensityBase
  {
   public:
    //! constructor
    explicit SPHDensityPredictCorrect(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHDensityPredictCorrect() override;

    //! init density handler
    void Init() override;

    //! setup density handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! write restart of density handler
    void WriteRestart(const int step, const double time) const override;

    //! read restart of density handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader) override;

    //! insert density evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const override;

    //! compute density field
    void ComputeDensity() const override;

   private:
    //! init density correction handler
    void InitDensityCorrectionHandler();

    //! correct density of interior/surface particles
    void CorrectDensity() const;

    //! density correction handler
    std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionBase> densitycorrection_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
