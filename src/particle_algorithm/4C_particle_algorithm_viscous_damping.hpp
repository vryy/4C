/*---------------------------------------------------------------------------*/
/*! \file
\brief viscous damping handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_VISCOUS_DAMPING_HPP
#define FOUR_C_PARTICLE_ALGORITHM_VISCOUS_DAMPING_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
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
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

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
FOUR_C_NAMESPACE_CLOSE

#endif
