/*---------------------------------------------------------------------------*/
/*! \file
\brief viscous damping handler for particle simulations

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_VISCOUS_DAMPING_H
#define PARTICLE_VISCOUS_DAMPING_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

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
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief viscous damping handler for particle simulations
   *
   * \author Sebastian Fuchs \date 02/2019
   */
  class ViscousDampingHandler
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \param[in] viscdampfac viscous damping factor
     */
    explicit ViscousDampingHandler(const double viscdampfac);

    //! destructor
    ~ViscousDampingHandler(){};

    /*!
     * \brief init viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     */
    void Init();

    /*!
     * \brief setup viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief write restart of viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief apply viscous damping contribution
     *
     * \author Sebastian Fuchs \date 02/2019
     */
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
