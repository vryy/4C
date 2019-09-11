/*---------------------------------------------------------------------------*/
/*! \file
\brief initial field handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INITIAL_FIELD_H
#define PARTICLE_INITIAL_FIELD_H

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

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class InitialFieldHandler
  {
   public:
    //! constructor
    explicit InitialFieldHandler(const Teuchos::ParameterList& params);

    //! destructor
    ~InitialFieldHandler(){};

    //! init initial field handler
    void Init();

    //! setup initial field handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! set initial fields
    void SetInitialFields();

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! relating particle types to function ids
    std::map<PARTICLEENGINE::StateEnum, std::map<PARTICLEENGINE::TypeEnum, int>>
        statetotypetofunctidmap_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
