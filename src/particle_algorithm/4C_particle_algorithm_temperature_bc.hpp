/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature boundary condition handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_TEMPERATURE_BC_HPP
#define FOUR_C_PARTICLE_ALGORITHM_TEMPERATURE_BC_HPP

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
   * \brief temperature boundary condition handler for particle simulations
   *
   * \author Sebastian Fuchs \date 09/2018
   */
  class TemperatureBoundaryConditionHandler
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[in] params particle simulation parameter list
     */
    explicit TemperatureBoundaryConditionHandler(const Teuchos::ParameterList& params);

    /*!
     * \brief init temperature boundary condition handler
     *
     * \author Sebastian Fuchs \date 09/2018
     */
    void Init();

    /*!
     * \brief setup temperature boundary condition handler
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief get reference to set of particle types subjected to temperature boundary conditions
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \return set of particle types subjected to temperature boundary conditions
     */
    const std::set<PARTICLEENGINE::TypeEnum>& get_particle_types_subjected_to_temperature_bc_set()
        const
    {
      return typessubjectedtotemperaturebc_;
    };

    /*!
     * \brief insert temperature boundary condition dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    /*!
     * \brief set particle reference position
     *
     * \author Sebastian Fuchs \date 09/2018
     */
    void set_particle_reference_position() const;

    /*!
     * \brief evaluate temperature boundary condition
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[in] evaltime evaluation time
     */
    void evaluate_temperature_boundary_condition(const double& evaltime) const;

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
FOUR_C_NAMESPACE_CLOSE

#endif
