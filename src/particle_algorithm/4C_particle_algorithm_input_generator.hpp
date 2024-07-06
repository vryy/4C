/*---------------------------------------------------------------------------*/
/*! \file
\brief particle input generator for particle simulations
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_INPUT_GENERATOR_HPP
#define FOUR_C_PARTICLE_ALGORITHM_INPUT_GENERATOR_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleObject;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief particle input generator
   *
   * \author Sebastian Fuchs \date 07/2018
   */
  class InputGenerator
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit InputGenerator(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief init input generator
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void init();

    /*!
     * \brief generate particles
     *
     * Generate initial particles in addition to particles read in from input file. Add all
     * generated particles to the vector using the function add_generated_particle().
     *
     * \note Either generate particles only on one processor or be sure that particles are not
     *       generated twice (or even more) on different processors.
     *
     * \note There is no need to set a global id at this stage. The unique global id handler sets
     *       the global id of all particles read in from the input file and generated herein later
     *       on.
     *
     * \note Think about reserving (not resizing!) the vector storing the generated particles in
     *       advance if you have a rough estimate of the total number of particles being generated.
     *       Keep in mind that particles read in from the input files are also stored in that
     *       vector.
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[out] particlesgenerated particle objects generated
     */
    void generate_particles(
        std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const;

   protected:
    /*!
     * \brief add generated particle
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[in]  position           position of particle
     * \param[in]  particletype       particle type enum
     * \param[out] particlesgenerated particle objects generated
     */
    void add_generated_particle(const std::vector<double>& position,
        const PARTICLEENGINE::TypeEnum particletype,
        std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& particlesgenerated) const;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
