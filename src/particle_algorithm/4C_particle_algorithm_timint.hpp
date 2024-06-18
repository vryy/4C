/*---------------------------------------------------------------------------*/
/*! \file
\brief time integration for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_TIMINT_HPP
#define FOUR_C_PARTICLE_ALGORITHM_TIMINT_HPP

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
namespace PARTICLEALGORITHM
{
  class DirichletBoundaryConditionHandler;
  class TemperatureBoundaryConditionHandler;
}  // namespace PARTICLEALGORITHM

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace ParticleRigidBody
{
  class RigidBodyHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief time integration base
   *
   * \author Sebastian Fuchs \date 04/2018
   */
  class TimInt
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] params particle simulation parameter list
     */
    explicit TimInt(const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    virtual ~TimInt();

    /*!
     * \brief init particle time integration
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    virtual void init();

    /*!
     * \brief time integration scheme specific initialization routine
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    virtual void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<ParticleRigidBody::RigidBodyHandlerInterface>
            particlerigidbodyinterface);

    /*!
     * \brief insert integration dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    /*!
     * \brief time integration scheme specific initialization routine
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    virtual void SetInitialStates();

    /*!
     * \brief set current time
     *
     * \author Sebastian Fuchs \date 08/2018
     *
     * \param[in] currenttime current time
     */
    virtual void set_current_time(const double currenttime) final;

    /*!
     * \brief time integration scheme specific pre-interaction routine
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    virtual void pre_interaction_routine() = 0;

    /*!
     * \brief time integration scheme specific post-interaction routine
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    virtual void post_interaction_routine() = 0;

   private:
    /*!
     * \brief init dirichlet boundary condition handler
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void init_dirichlet_boundary_condition();

    /*!
     * \brief init temperature boundary condition handler
     *
     * \author Sebastian Fuchs \date 09/2018
     */
    void init_temperature_boundary_condition();

    /*!
     * \brief add initial random noise to particle position
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void add_initial_random_noise_to_position();

   protected:
    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! interface to rigid body handler
    std::shared_ptr<ParticleRigidBody::RigidBodyHandlerInterface> particlerigidbodyinterface_;

    //! dirichlet boundary condition handler
    std::unique_ptr<PARTICLEALGORITHM::DirichletBoundaryConditionHandler>
        dirichletboundarycondition_;

    //! temperature boundary condition handler
    std::unique_ptr<PARTICLEALGORITHM::TemperatureBoundaryConditionHandler>
        temperatureboundarycondition_;

    //! set of particle types to integrate in time
    std::set<PARTICLEENGINE::TypeEnum> typestointegrate_;

    //! current time
    double time_;

    //! time step size
    double dt_;
  };

  /*!
   * \brief semi-implicit Euler time integration scheme
   *
   * semi-implicit Euler time integration scheme of first order accuracy (also denoted as
   * semi-explicit Euler or symplectic Euler scheme)
   *
   * \$f v_{n+1} = v_{n} + dt * a_{n}   \$f with \$f a_{n} = a( r_{n-1}, v_{n-1} ) \$f
   * \$f x_{n+1} = x_{n} + dt * v_{n+1} \$f
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class TimIntSemiImplicitEuler : public TimInt
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] params particle simulation parameter list
     */
    TimIntSemiImplicitEuler(const Teuchos::ParameterList& params);

    /*!
     * \brief time integration scheme specific initialization routine
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<ParticleRigidBody::RigidBodyHandlerInterface>
            particlerigidbodyinterface) override;

    /*!
     * \brief time integration scheme specific pre-interaction routine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void pre_interaction_routine() override;

    /*!
     * \brief time integration scheme specific post-interaction routine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void post_interaction_routine() override;
  };

  /*!
   * \brief explicit velocity Verlet time integration scheme
   *
   * explicit velocity Verlet time integration scheme of second order accuracy (also denoted as
   * leapfrog scheme in kick-drift-kick form)
   *
   * \$f v_{n+1/2} = v_{n}     + dt/2 * a_{n}     \$f with \$f a_{n}   = a( r_{n}, v_{n-1/2} )   \$f
   * \$f x_{n+1}   = x_{n}     + dt   * v_{n+1/2} \$f
   * \$f v_{n+1}   = v_{n+1/2} + dt/2 * a_{n+1}   \$f with \$f a_{n+1} = a( r_{n+1}, v_{n+1/2} ) \$f
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class TimIntVelocityVerlet : public TimInt
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] params particle simulation parameter list
     */
    TimIntVelocityVerlet(const Teuchos::ParameterList& params);

    /*!
     * \brief time integration scheme specific initialization routine
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void SetInitialStates() override;

    /*!
     * \brief time integration scheme specific pre-interaction routine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void pre_interaction_routine() override;

    /*!
     * \brief time integration scheme specific post-interaction routine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void post_interaction_routine() override;

   private:
    //! half time step size
    double dthalf_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
