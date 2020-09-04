/*---------------------------------------------------------------------------*/
/*! \file
\brief density handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
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

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
  class MaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
  class SPHVirtualWallParticle;
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
    virtual ~SPHDensityBase() = default;

    //! init density handler
    virtual void Init();

    //! setup density handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize) final;

    //! insert density evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const = 0;

    //! compute density field
    virtual void ComputeDensity() const = 0;

   protected:
    //! evaluate sum of weighted mass
    virtual void SumWeightedMass() const final;

    //! clear density sum state
    virtual void ClearDensitySumState() const final;

    //! sum weighted mass (self contribution)
    virtual void SumWeightedMassSelfContribution() const final;

    //! sum weighted mass (particle contribution)
    virtual void SumWeightedMassParticleContribution() const final;

    //! sum weighted mass (particle-wall contribution)
    virtual void SumWeightedMassParticleWallContribution() const final;

    //! evaluate sum of colorfield
    virtual void SumColorfield() const final;

    //! clear colorfield state
    virtual void ClearColorfieldState() const final;

    //! sum colorfield (self contribution)
    virtual void SumColorfieldSelfContribution() const final;

    //! sum colorfield (particle contribution)
    virtual void SumColorfieldParticleContribution() const final;

    //! sum colorfield (particle-wall contribution)
    virtual void SumColorfieldParticleWallContribution() const final;

    //! evaluate continuity equation
    virtual void ContinuityEquation() const final;

    //! clear density dot state
    virtual void ClearDensityDotState() const final;

    //! continuity equation (particle contribution)
    virtual void ContinuityEquationParticleContribution() const final;

    //! continuity equation (particle-wall contribution)
    virtual void ContinuityEquationParticleWallContribution() const final;

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

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! virtual wall particle handler
    std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle_;

    //! density of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh densitytorefresh_;

    //! set of fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> fluidtypes_;

    //! time step size
    double dt_;
  };

  class SPHDensitySummation final : public SPHDensityBase
  {
   public:
    //! constructor
    explicit SPHDensitySummation(const Teuchos::ParameterList& params);

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

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~SPHDensityPredictCorrect() override;

    //! init density handler
    void Init() override;

    //! setup density handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle)
        override;

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
