/*---------------------------------------------------------------------------*/
/*! \file
\brief algorithm to control particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_HPP
#define FOUR_C_PARTICLE_ALGORITHM_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_particle_engine.hpp"
#include "4C_particle_engine_typedefs.hpp"
#include "4C_particle_wall.hpp"
#include "4C_utils_result_test.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class TimInt;
  class GravityHandler;
  class ViscousDampingHandler;
}  // namespace PARTICLEALGORITHM

namespace PARTICLEINTERACTION
{
  class ParticleInteractionBase;
}

namespace PARTICLEENGINE
{
  class ParticleObject;
}

namespace PARTICLERIGIDBODY
{
  class RigidBodyHandler;
}

namespace DRT
{
  class ResultTest;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief algorithm to control particle simulations
   *
   * The particle algorithm controls all major steps of particle simulations. Besides initialization
   * and setup of the problem this mainly concerns the control of the time loop (with the
   * integration of the particle step) consisting of the following steps:
   * - update and maintain valid connectivity (parallel distribution, load transfer, ghosting,
   *   neighbor pair relations)
   * - evaluation of particle accelerations from gravity acceleration, interactions, and viscous
   *   damping
   * - (runtime) output of fields
   * - write restart information
   *
   * \author Sebastian Fuchs \date 04/2018
   */

  class ParticleAlgorithm final : public ADAPTER::AlgorithmBase
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit ParticleAlgorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~ParticleAlgorithm() override;

    /*!
     * \brief init particle algorithm
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] initialparticles particle objects read in from restart
     */
    void Init(std::vector<PARTICLEENGINE::ParticleObjShrdPtr>& initialparticles);

    /*!
     * \brief setup particle algorithm
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void Setup();

    /*!
     * \brief read restart information for given time step
     *
     * Read the restart information in the same order as it is written.
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] restartstep restart step
     */
    void read_restart(const int restartstep) override;

    /*!
     * \brief time loop for particle problem
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void Timeloop();

    /*!
     * \brief prepare time step
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] print_header flag to control output of time step header
     */
    void prepare_time_step(bool print_header = true);

    /*!
     * \brief pre evaluate time step
     *
     * \author Sebastian Fuchs \date 11/2020
     */
    void PreEvaluateTimeStep();

    /*!
     * \brief integrate time step
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void IntegrateTimeStep();

    /*!
     * \brief post evaluate time step
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void post_evaluate_time_step();

    /*!
     * \brief write output
     *
     * \author Sebastian Fuchs \date 07/2020
     */
    void WriteOutput() const;

    /*!
     * \brief write restart information
     *
     * The order of writing the restart information in the particle framework is important: First
     * restart information is written by the bin discretization writer, afterwards restart
     * information is written by the wall discretization writer. This order is important, otherwise
     * the restart fails.
     *
     * \author Sebastian Fuchs \date 07/2020
     */
    void WriteRestart() const;

    /*!
     * \brief create particle field specific result test objects
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \return particle field specific result test objects
     */
    std::vector<std::shared_ptr<CORE::UTILS::ResultTest>> CreateResultTests();

    /*!
     * \brief get interface to particle engine
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \return interface to particle engine
     */
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> get_particle_engine_interface() const
    {
      return particleengine_;
    }

    /*!
     * \brief get interface to particle wall handler
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \return interface to particle wall handler
     */
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> get_particle_wall_handler_interface() const
    {
      return particlewall_;
    }

   private:
    //! \name init and setup methods
    //! @{

    /*!
     * \brief init particle engine
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void init_particle_engine();

    /*!
     * \brief init particle wall handler
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void init_particle_wall();

    /*!
     * \brief init rigid body handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void init_particle_rigid_body();

    /*!
     * \brief init particle time integration
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void init_particle_time_integration();

    /*!
     * \brief init particle interaction handler
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void init_particle_interaction();

    /*!
     * \brief init particle gravity handler
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void init_particle_gravity();

    /*!
     * \brief init viscous damping handler
     *
     * \author Sebastian Fuchs \date 02/2019
     */
    void init_viscous_damping();

    /*!
     * \brief generate initial particles
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void generate_initial_particles();

    /*!
     * \brief determine all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void determine_particle_types();

    /*!
     * \brief determine particle states of all particle types
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void determine_particle_states_of_particle_types();

    /*!
     * \brief setup initial particles
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void setup_initial_particles();

    /*!
     * \brief setup initial rigid bodies
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void setup_initial_rigid_bodies();

    /*!
     * \brief setup initial states
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void setup_initial_states();

    //! @}

    /*!
     * \brief update and maintain valid connectivity in particle problems
     *
     * The main control routine of particle problems to maintain valid particle connectivity:
     * - parallel distribution of particles
     * - transfer of particles to new bins
     * - ghosting of particles
     * - neighbor pair relations
     *
     * This method includes several safety checks that are part of the called sub-methods.
     *
     * \author Sebastian Fuchs \date 06/2018
     */
    void update_connectivity();

    /*!
     * \brief check load transfer
     *
     * Check if a load transfer is needed based on the following criteria:
     * - maximum position increment of particles or walls
     * - invalid particle connectivity
     * - invalid particle neighbor pair relation
     * - invalid particle wall neighbor pair relation
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load transfer is needed
     */
    bool check_load_transfer_needed();

    /*!
     * \brief check maximum position increment
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load transfer is needed due to maximum position increment
     */
    bool check_max_position_increment();

    /*!
     * \brief get maximum particle position increment since last transfer
     *
     * \author Sebastian Fuchs \date 06/2018
     *
     * \param allprocmaxpositionincrement maximum particle position increment
     */
    void get_max_particle_position_increment(double& allprocmaxpositionincrement);

    /*!
     * \brief transfer load between processors
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void transfer_load_between_procs();

    /*!
     * \brief check load redistribution
     *
     * Check if a load redistribution is needed based on the difference in the absolute number of
     * particles since the last redistribution.
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \return flag indicating load redistribution is needed
     */
    bool check_load_redistribution_needed();

    /*!
     * \brief distribute load among processors
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void distribute_load_among_procs();

    /*!
     * \brief build potential neighbor relation
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void build_potential_neighbor_relation();

    /*!
     * \brief set initial conditions
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void set_initial_conditions();

    /*!
     * \brief set current time
     *
     * \author Sebastian Fuchs \date 08/2018
     */
    void set_current_time();

    /*!
     * \brief set current step size
     *
     * \author Sebastian Fuchs \date 08/2018
     */
    void set_current_step_size();

    /*!
     * \brief set current write result flag
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    void set_current_write_result_flag();

    /*!
     * \brief evaluate time step
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void evaluate_time_step();

    /*!
     * \brief set gravity acceleration
     *
     * \author Sebastian Fuchs \date 06/2018
     */
    void set_gravity_acceleration();

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngine> particleengine_;

    //! particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerBase> particlewall_;

    //! rigid body handler
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandler> particlerigidbody_;

    //! particle time integration
    std::unique_ptr<PARTICLEALGORITHM::TimInt> particletimint_;

    //! particle interaction
    std::unique_ptr<PARTICLEINTERACTION::ParticleInteractionBase> particleinteraction_;

    //! particle gravity handler
    std::unique_ptr<PARTICLEALGORITHM::GravityHandler> particlegravity_;

    //! particle viscous damping handler
    std::unique_ptr<PARTICLEALGORITHM::ViscousDampingHandler> viscousdamping_;

    //! map of particle types and corresponding states
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes_;

    //! vector of initial or generated particles to distribute
    std::vector<PARTICLEENGINE::ParticleObjShrdPtr> particlestodistribute_;

    //! number of particles on this processor after last load balance
    int numparticlesafterlastloadbalance_;

    //! transfer particles to new bins every time step
    bool transferevery_;

    //! write results interval
    const int writeresultsevery_;

    //! write restart interval
    const int writerestartevery_;

    //! result control flag
    bool writeresultsthisstep_;

    //! restart control flag
    bool writerestartthisstep_;

    //! simulation is restarted
    bool isrestarted_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
