/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_pressure.H

\brief pressure handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_PRESSURE_H
#define PARTICLE_INTERACTION_SPH_PRESSURE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
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
  class MaterialHandler;
  class SPHEquationOfStateBundle;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHPressure
  {
   public:
    //! constructor
    explicit SPHPressure();

    //! virtual destructor
    virtual ~SPHPressure(){};

    //! init pressure handler
    virtual void Init();

    //! setup pressure handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle);

    //! write restart of pressure handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of pressure handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! compute pressure using equation of state and density
    virtual void ComputePressure() const;

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
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
