/*---------------------------------------------------------------------------*/
/*! \file
\brief pressure handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_PRESSURE_H
#define PARTICLE_INTERACTION_SPH_PRESSURE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class SPHEquationOfStateBundle;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHPressure final
  {
   public:
    //! constructor
    explicit SPHPressure();

    //! destructor
    ~SPHPressure() = default;

    //! init pressure handler
    void Init();

    //! setup pressure handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle);

    //! compute pressure using equation of state and density
    void ComputePressure() const;

   protected:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! pressure of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh pressuretorefresh_;

    //! set of fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> fluidtypes_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
