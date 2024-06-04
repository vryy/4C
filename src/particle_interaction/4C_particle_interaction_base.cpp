/*---------------------------------------------------------------------------*/
/*! \file
\brief base particle interaction handler
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_base.hpp"

#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_runtime_writer.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionBase::ParticleInteractionBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), myrank_(comm.MyPID()), params_(params), time_(0.0), dt_(0.0)
{
  // empty constructor
}

void PARTICLEINTERACTION::ParticleInteractionBase::Init()
{
  // init particle material handler
  init_particle_material_handler();

  // init particle interaction writer
  init_particle_interaction_writer();
}

void PARTICLEINTERACTION::ParticleInteractionBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->get_particle_container_bundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // setup particle material handler
  particlematerial_->Setup();

  // setup particle interaction writer
  particleinteractionwriter_->Setup();

  // init vector
  gravity_.resize(3, 0.0);
}

void PARTICLEINTERACTION::ParticleInteractionBase::write_restart() const
{
  // nothing to do
}

void PARTICLEINTERACTION::ParticleInteractionBase::read_restart(
    const std::shared_ptr<CORE::IO::DiscretizationReader> reader)
{
  // read restart of particle interaction writer
  particleinteractionwriter_->read_restart(reader);
}

void PARTICLEINTERACTION::ParticleInteractionBase::
    check_particle_interaction_distance_concerning_bin_size() const
{
  // get maximum particle interaction distance
  double allprocmaxinteractiondistance = 0.0;
  double maxinteractiondistance = max_interaction_distance();
  comm_.MaxAll(&maxinteractiondistance, &allprocmaxinteractiondistance, 1);

  // bin size safety check
  if (allprocmaxinteractiondistance > particleengineinterface_->MinBinSize())
    FOUR_C_THROW("the particle interaction distance is larger than the minimal bin size (%f > %f)!",
        allprocmaxinteractiondistance, particleengineinterface_->MinBinSize());

  // periodic length safety check
  if (particleengineinterface_->have_periodic_boundary_conditions())
  {
    // loop over all spatial directions
    for (int dim = 0; dim < 3; ++dim)
    {
      // check for periodic boundary condition in current spatial direction
      if (not particleengineinterface_->have_periodic_boundary_conditions_in_spatial_direction(dim))
        continue;

      // check periodic length in current spatial direction
      if ((2.0 * allprocmaxinteractiondistance) >
          particleengineinterface_->length_of_binning_domain_in_a_spatial_direction(dim))
        FOUR_C_THROW(
            "particles are not allowed to interact directly and across the periodic boundary!");
    }
  }
}

void PARTICLEINTERACTION::ParticleInteractionBase::set_current_time(const double currenttime)
{
  time_ = currenttime;
}

void PARTICLEINTERACTION::ParticleInteractionBase::set_current_step_size(
    const double currentstepsize)
{
  dt_ = currentstepsize;
}

void PARTICLEINTERACTION::ParticleInteractionBase::set_current_write_result_flag(
    bool writeresultsthisstep)
{
  // set current write result flag in particle interaction writer
  particleinteractionwriter_->set_current_write_result_flag(writeresultsthisstep);
}

void PARTICLEINTERACTION::ParticleInteractionBase::SetGravity(std::vector<double>& gravity)
{
  gravity_ = gravity;
}

void PARTICLEINTERACTION::ParticleInteractionBase::write_interaction_runtime_output(
    const int step, const double time)
{
  // write particle interaction runtime output
  particleinteractionwriter_->write_particle_interaction_runtime_output(step, time);
}

void PARTICLEINTERACTION::ParticleInteractionBase::init_particle_material_handler()
{
  // create particle material handler
  particlematerial_ = std::make_shared<PARTICLEINTERACTION::MaterialHandler>(params_);

  // init particle material handler
  particlematerial_->Init();
}

void PARTICLEINTERACTION::ParticleInteractionBase::init_particle_interaction_writer()
{
  // create particle interaction writer
  particleinteractionwriter_ =
      std::make_shared<PARTICLEINTERACTION::InteractionWriter>(comm_, params_);

  // init particle interaction writer
  particleinteractionwriter_->Init();
}

double PARTICLEINTERACTION::ParticleInteractionBase::max_particle_radius() const
{
  // init value of maximum radius
  double maxrad = 0.0;

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->get_specific_container(type_i, PARTICLEENGINE::Owned);

    // get maximum stored value of state
    double currmaxrad = container->GetMaxValueOfState(PARTICLEENGINE::Radius);

    // compare to current maximum
    maxrad = std::max(maxrad, currmaxrad);
  }

  return maxrad;
}

FOUR_C_NAMESPACE_CLOSE
