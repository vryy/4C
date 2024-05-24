/*---------------------------------------------------------------------------*/
/*! \file
\brief base particle interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_BASE_HPP
#define FOUR_C_PARTICLE_INTERACTION_BASE_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class InteractionWriter;
}  // namespace PARTICLEINTERACTION

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  /*!
   * \brief base particle interaction
   *
   * \author Sebastian Fuchs \date 05/2018
   */
  class ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~ParticleInteractionBase() = default;

    //! init particle interaction handler
    virtual void Init();

    //! setup particle interaction handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! write restart of particle interaction handler
    virtual void WriteRestart() const;

    //! read restart of particle interaction handler
    virtual void read_restart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert interaction dependent states of all particle types
    virtual void insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) = 0;

    //! set initial states
    virtual void SetInitialStates() = 0;

    //! pre evaluate time step
    virtual void PreEvaluateTimeStep() = 0;

    //! evaluate particle interactions
    virtual void evaluate_interactions() = 0;

    //! post evaluate time step
    virtual void post_evaluate_time_step(
        std::vector<PARTICLEENGINE::ParticleTypeToType>& particlesfromphasetophase) = 0;

    //! check particle interaction distance concerning bin size
    virtual void check_particle_interaction_distance_concerning_bin_size() const final;

    //! maximum interaction distance (on this processor)
    virtual double max_interaction_distance() const = 0;

    //! distribute interaction history
    virtual void distribute_interaction_history() const = 0;

    //! communicate interaction history
    virtual void communicate_interaction_history() const = 0;

    //! set current time
    virtual void set_current_time(const double currenttime);

    //! set current step size
    virtual void set_current_step_size(const double currentstepsize);

    //! set current write result flag
    virtual void set_current_write_result_flag(bool writeresultsthisstep);

    //! set gravity
    virtual void SetGravity(std::vector<double>& gravity) final;

    //! write interaction runtime output
    virtual void write_interaction_runtime_output(const int step, const double time);

   private:
    //! init particle material handler
    virtual void init_particle_material_handler();

    //! init particle interaction writer
    virtual void init_particle_interaction_writer();

   protected:
    //! maximum particle radius (on this processor)
    virtual double MaxParticleRadius() const;

    //! communication
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! particle interaction writer
    std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter_;

    //! current time
    double time_;

    //! time step size
    double dt_;

    //! current gravity vector
    std::vector<double> gravity_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
