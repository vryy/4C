/*---------------------------------------------------------------------------*/
/*! \file
\brief evaporation induced heat loss handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_HEATLOSS_EVAPORATION_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_HEATLOSS_EVAPORATION_HPP

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

namespace ParticleInteraction
{
  class MaterialHandler;
}

namespace Mat
{
  namespace PAR
  {
    class ParticleMaterialThermo;
  }
}  // namespace Mat

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
{
  class SPHHeatLossEvaporation
  {
   public:
    //! constructor
    explicit SPHHeatLossEvaporation(const Teuchos::ParameterList& params);

    //! init evaporation induced heat loss handler
    void init();

    //! setup evaporation induced heat loss handler
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<ParticleInteraction::MaterialHandler> particlematerial);

    //! evaluate evaporation induced heat loss
    void evaluate_evaporation_induced_heat_loss() const;

   protected:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle material handler
    std::shared_ptr<ParticleInteraction::MaterialHandler> particlematerial_;

    //! pointer to thermo material of particle types
    std::vector<const Mat::PAR::ParticleMaterialThermo*> thermomaterial_;

    //! evaporating phase
    PARTICLEENGINE::TypeEnum evaporatingphase_;

    //! boiling temperature in recoil pressure formula
    double recoilboilingtemp_;

    //! pressure factor in recoil pressure formula
    double recoil_pfac_;

    //! temperature factor in recoil pressure formula
    double recoil_tfac_;

    //! latent heat in heat loss formula
    double latentheat_;

    //! enthalpy reference temperature in heat loss formula
    double enthalpyreftemp_;

    //! pressure factor in heat loss formula
    double heatloss_pfac_;

    //! temperature factor in heat loss formula
    double heatloss_tfac_;
  };

}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
