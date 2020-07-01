/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_TEMPERATURE_H
#define PARTICLE_INTERACTION_SPH_TEMPERATURE_H

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
  class MaterialHandler;
  class SPHNeighborPairs;
  class SPHHeatSourceBase;
}  // namespace PARTICLEINTERACTION

namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialThermo;
  }
}  // namespace MAT

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHTemperature final
  {
   public:
    //! constructor
    explicit SPHTemperature(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHTemperature();

    //! init temperature handler
    void Init();

    //! setup temperature handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of temperature handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of temperature handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current time
    void SetCurrentTime(const double currenttime);

    //! set current step size
    void SetCurrentStepSize(const double currentstepsize);

    //! insert temperature evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! compute temperature field using energy equation
    void ComputeTemperature() const;

   private:
    //! init heat source handler
    void InitHeatSourceHandler();

    //! evaluate energy equation
    void EnergyEquation() const;

    //! evaluate temperature gradient
    void TemperatureGradient() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! heat source handler
    std::unique_ptr<PARTICLEINTERACTION::SPHHeatSourceBase> heatsource_;

    //! temperature of ghosted particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh temptorefresh_;

    //! current time
    double time_;

    //! time step size
    double dt_;

    //! evaluate temperature gradient
    bool temperaturegradient_;

    //! pointer to thermo material of particle types
    std::vector<const MAT::PAR::ParticleMaterialThermo*> thermomaterial_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
