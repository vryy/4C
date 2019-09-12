/*---------------------------------------------------------------------------*/
/*! \file
\brief dirichlet boundary condition handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_DIRICHLET_BC_H
#define PARTICLE_DIRICHLET_BC_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2018 |
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
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class DirichletBoundaryConditionHandler
  {
   public:
    //! constructor
    explicit DirichletBoundaryConditionHandler(const Teuchos::ParameterList& params);

    //! destructor
    ~DirichletBoundaryConditionHandler(){};

    //! init dirichlet boundary condition handler
    void Init();

    //! setup dirichlet boundary condition handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! write restart of dirichlet boundary condition handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of dirichlet oundary condition handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to particle type set
    const std::set<PARTICLEENGINE::TypeEnum>& GetParticleTypesSubjectedToDirichletBCSet() const
    {
      return typessubjectedtodirichletbc_;
    };

    //! insert dbc dependent states of all particle types
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! set particle reference position
    void SetParticleReferencePosition() const;

    //! evaluate dirichlet boundary condition
    void EvaluateDirichletBoundaryCondition(
        const double& evaltime, const bool evalpos, const bool evalvel, const bool evalacc) const;

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! relating particle types to function ids of dirichlet boundary conditions
    std::map<PARTICLEENGINE::TypeEnum, int> dirichletbctypetofunctid_;

    //! set of particle types subjected to dirichlet boundary conditions
    std::set<PARTICLEENGINE::TypeEnum> typessubjectedtodirichletbc_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
