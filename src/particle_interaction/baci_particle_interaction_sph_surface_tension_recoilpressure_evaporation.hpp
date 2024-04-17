/*---------------------------------------------------------------------------*/
/*! \file
\brief evaporation induced recoil pressure handler for smoothed particle hydrodynamics (SPH)
       interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_RECOILPRESSURE_EVAPORATION_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_RECOILPRESSURE_EVAPORATION_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_inpar_particle.hpp"
#include "baci_particle_engine_enums.hpp"
#include "baci_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHRecoilPressureEvaporation
  {
   public:
    //! constructor
    explicit SPHRecoilPressureEvaporation(const Teuchos::ParameterList& params);

    //! init evaporation induced recoil pressure handler
    void Init();

    //! setup evaporation induced recoil pressure handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! compute evaporation induced recoil pressure contribution
    void ComputeRecoilPressureContribution() const;

   protected:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! evaporating phase
    PARTICLEENGINE::TypeEnum evaporatingphase_;

    //! boiling temperature in recoil pressure formula
    double recoilboilingtemp_;

    //! pressure factor in recoil pressure formula
    double recoil_pfac_;

    //! temperature factor in recoil pressure formula
    double recoil_tfac_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
