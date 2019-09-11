/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature boundary condition handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                 meier 09/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_TEMPERATURE_BC_H
#define PARTICLE_TEMPERATURE_BC_H

/*---------------------------------------------------------------------------*
 | headers                                                     meier 09/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                        meier 09/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                           meier 09/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class TemperatureBoundaryConditionHandler
  {
   public:
    //! constructor
    explicit TemperatureBoundaryConditionHandler(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~TemperatureBoundaryConditionHandler(){};

    //! init temperature boundary condition handler
    void Init();

    //! setup temperature boundary condition handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! write restart of temperature boundary condition handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of temperature oundary condition handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to particle type set
    const std::set<PARTICLEENGINE::TypeEnum>& GetParticleTypesSubjectedToTemperatureBCSet() const
    {
      return typessubjectedtotemperaturebc_;
    };

    //! insert tempbc dependent states of all particle types
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! set particle reference position
    void SetParticleReferencePosition() const;

    //! evaluate temperature boundary condition
    void EvaluateTemperatureBoundaryCondition(const double& evaltime) const;

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! relating particle types to function ids of temperature boundary conditions
    std::map<PARTICLEENGINE::TypeEnum, int> temperaturebctypetofunctid_;

    //! set of particle types subjected to temperature boundary conditions
    std::set<PARTICLEENGINE::TypeEnum> typessubjectedtotemperaturebc_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
