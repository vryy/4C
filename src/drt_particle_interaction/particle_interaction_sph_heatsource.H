/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_heatsource.H

\brief heat source handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_HEATSOURCE_H
#define PARTICLE_INTERACTION_SPH_HEATSOURCE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 02/2019 |
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
}

namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialThermo;
  }
}  // namespace MAT

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHHeatSourceBase
  {
   public:
    //! constructor
    explicit SPHHeatSourceBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHHeatSourceBase(){};

    //! init heat source handler
    virtual void Init();

    //! setup heat source handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial);

    //! write restart of heat source handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of heat source handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! evaluate heat source
    virtual void EvaluateHeatSource(const double& evaltime) const = 0;

   protected:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! pointer to thermo material of particle types
    std::vector<const MAT::PAR::ParticleMaterialThermo*> thermomaterial_;

    //! heat source function number
    const int heatsourcefctnumber_;
  };

  class SPHHeatSourceVolume : public SPHHeatSourceBase
  {
   public:
    //! constructor
    explicit SPHHeatSourceVolume(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHHeatSourceVolume(){};

    //! evaluate heat source
    virtual void EvaluateHeatSource(const double& evaltime) const override;
  };

  class SPHHeatSourceSurface : public SPHHeatSourceBase
  {
   public:
    //! constructor
    explicit SPHHeatSourceSurface(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHHeatSourceSurface(){};

    //! evaluate heat source
    virtual void EvaluateHeatSource(const double& evaltime) const override;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
