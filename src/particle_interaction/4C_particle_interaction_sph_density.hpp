/*---------------------------------------------------------------------------*/
/*! \file
\brief density handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_DENSITY_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_DENSITY_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_particle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
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
    virtual void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const = 0;

    //! compute density field
    virtual void ComputeDensity() const = 0;

   protected:
    //! evaluate sum of weighted mass
    virtual void SumWeightedMass() const final;

    //! clear density sum state
    virtual void clear_density_sum_state() const final;

    //! sum weighted mass (self contribution)
    virtual void sum_weighted_mass_self_contribution() const final;

    //! sum weighted mass (particle contribution)
    virtual void sum_weighted_mass_particle_contribution() const final;

    //! sum weighted mass (particle-wall contribution)
    virtual void sum_weighted_mass_particle_wall_contribution() const final;

    //! evaluate sum of colorfield
    virtual void SumColorfield() const final;

    //! clear colorfield state
    virtual void clear_colorfield_state() const final;

    //! sum colorfield (self contribution)
    virtual void sum_colorfield_self_contribution() const final;

    //! sum colorfield (particle contribution)
    virtual void sum_colorfield_particle_contribution() const final;

    //! sum colorfield (particle-wall contribution)
    virtual void sum_colorfield_particle_wall_contribution() const final;

    //! evaluate continuity equation
    virtual void ContinuityEquation() const final;

    //! clear density dot state
    virtual void clear_density_dot_state() const final;

    //! continuity equation (particle contribution)
    virtual void continuity_equation_particle_contribution() const final;

    //! continuity equation (particle-wall contribution)
    virtual void continuity_equation_particle_wall_contribution() const final;

    //! set density sum to density field
    virtual void SetDensitySum() const final;

    //! add time step scaled density dot to density field
    virtual void add_time_step_scaled_density_dot() const final;

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
    void insert_particle_states_of_particle_types(
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
    void insert_particle_states_of_particle_types(
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
    void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const override;

    //! compute density field
    void ComputeDensity() const override;

   private:
    //! init density correction handler
    void init_density_correction_handler();

    //! correct density of interior/surface particles
    void CorrectDensity() const;

    //! density correction handler
    std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionBase> densitycorrection_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
