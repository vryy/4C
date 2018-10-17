/*---------------------------------------------------------------------------*/
/*!
\file particle_input_generator.H

\brief particle input generator for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INPUT_GENERATOR_H
#define PARTICLE_INPUT_GENERATOR_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleObject;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class InputGenerator
  {
   public:
    //! constructor
    explicit InputGenerator(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~InputGenerator(){};

    //! init input generator
    virtual void Init();

    //! generate particles
    void GenerateParticles(
        std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const;

   protected:
    //! add generated particle
    void AddGeneratedParticle(const std::vector<double>& position,
        const PARTICLEENGINE::TypeEnum particletype,
        std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
