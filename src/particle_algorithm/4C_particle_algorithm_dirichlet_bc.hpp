/*---------------------------------------------------------------------------*/
/*! \file
\brief dirichlet boundary condition handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_DIRICHLET_BC_HPP
#define FOUR_C_PARTICLE_ALGORITHM_DIRICHLET_BC_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

#include <Teuchos_ParameterList.hpp>

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
   * \brief dirichlet boundary condition handler for particle simulations
   *
   * \author Sebastian Fuchs \date 07/2018
   */
  class DirichletBoundaryConditionHandler
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] params particle simulation parameter list
     */
    explicit DirichletBoundaryConditionHandler(const Teuchos::ParameterList& params);

    /*!
     * \brief init dirichlet boundary condition handler
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void init();

    /*!
     * \brief setup dirichlet boundary condition handler
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief get reference to set of particle types subjected to dirichlet boundary conditions
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \return set of particle types subjected to dirichlet boundary conditions
     */
    const std::set<PARTICLEENGINE::TypeEnum>& get_particle_types_subjected_to_dirichlet_bc_set()
        const
    {
      return typessubjectedtodirichletbc_;
    };

    /*!
     * \brief insert dirichlet boundary condition dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    /*!
     * \brief set particle reference position
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void set_particle_reference_position() const;

    /*!
     * \brief evaluate dirichlet boundary condition
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] evaltime evaluation time
     * \param[in] evalpos  flag to indicate evaluation of position
     * \param[in] evalvel  flag to indicate evaluation of velocity
     * \param[in] evalacc  flag to indicate evaluation of acceleration
     */
    void evaluate_dirichlet_boundary_condition(
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
FOUR_C_NAMESPACE_CLOSE

#endif
