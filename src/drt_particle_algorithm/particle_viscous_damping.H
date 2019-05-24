/*---------------------------------------------------------------------------*/
/*!
\brief viscous damping handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_VISCOUS_DAMPING_H
#define PARTICLE_VISCOUS_DAMPING_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

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
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class ViscousDampingHandler
  {
   public:
    //! constructor
    explicit ViscousDampingHandler(const double viscdampfac);

    //! virtual destructor
    ~ViscousDampingHandler(){};

    //! init viscous damping handler
    void Init();

    //! setup viscous damping handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! write restart of viscous damping handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of viscous damping handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! apply viscous damping contribution
    void ApplyViscousDamping();

   private:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! viscous damping factor
    const double viscdampfac_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
