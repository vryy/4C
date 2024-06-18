/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particles in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_engine_runtime_vtp_writer.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEENGINE::ParticleRuntimeVtpWriter::ParticleRuntimeVtpWriter(const Epetra_Comm& comm)
    : comm_(comm), setuptime_(0.0)
{
  // empty constructor
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::Init(
    const ParticleContainerBundleShrdPtr particlecontainerbundle)
{
  // set particle container bundle
  particlecontainerbundle_ = particlecontainerbundle;

  // insert specific particle states in black list
  blackliststates_.insert({DensitySum, DensityDot});
  blackliststates_.insert(TemperatureDot);
  blackliststates_.insert({LastTransferPosition, ReferencePosition});
  blackliststates_.insert({ModifiedVelocity, ModifiedAcceleration});
  blackliststates_.insert({InterfaceNormal, Curvature, WallColorfield, WallInterfaceNormal});
  blackliststates_.insert({LastIterPosition, LastIterVelocity, LastIterAcceleration,
      LastIterAngularVelocity, LastIterAngularAcceleration, LastIterModifiedAcceleration,
      LastIterDensity, LastIterTemperature});
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::setup(bool write_ghosted_particles)
{
  // determine size of vector indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  runtime_visualization_managers_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // allocate memory for vtp writer objects of owned and ghosted states
    (runtime_visualization_managers_[type]).resize(2);

    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      if (status == Ghosted and write_ghosted_particles == false) continue;

      std::ostringstream fieldname;
      fieldname << "particle-" << EnumToTypeName(type) << "-" << EnumToStatusName(status);

      // construct visualiation manager object for current particle type and status
      (runtime_visualization_managers_[type])[status] =
          std::make_shared<Core::IO::VisualizationManager>(
              Core::IO::VisualizationParametersFactory(
                  Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                  *Global::Problem::Instance()->OutputControlFile(), setuptime_),
              comm_, fieldname.str());
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::read_restart(
    const std::shared_ptr<Core::IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->read_double("time");
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::set_particle_positions_and_states()
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_visualization_managers_[type])[status]) continue;

      // get container of current particle type and status
      ParticleContainer* container = particlecontainerbundle_->get_specific_container(type, status);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // get particle states stored in container
      const std::set<ParticleState>& states = container->GetStoredStates();

#ifdef FOUR_C_ENABLE_ASSERTIONS
      // safety check
      if (not container->HaveStoredState(Position))
        FOUR_C_THROW("particle state '%s' not found!", EnumToStateName(Position).c_str());
#endif

      // iterate over particle states
      for (const auto& state : states)
      {
        // get particle state dimension
        int statedim = container->GetStateDim(state);

        // get name of particle state
        std::string statename = EnumToStateName(state);

        // get pointer to particle state
        const double* state_ptr =
            (particlestored > 0) ? container->GetPtrToState(state, 0) : nullptr;

        if (state == Position)
        {
          // get and prepare storage for position data
          std::vector<double>& positiondata = (runtime_visualization_managers_[type])[status]
                                                  ->get_visualization_data()
                                                  .GetPointCoordinates();
          positiondata.clear();
          positiondata.reserve(statedim * particlestored);

          // copy particle position data
          for (int i = 0; i < (statedim * particlestored); ++i)
            positiondata.push_back(state_ptr[i]);

#ifdef FOUR_C_ENABLE_ASSERTIONS
          // safety check
          if (static_cast<int>(positiondata.size()) != statedim * particlestored)
            FOUR_C_THROW("ParticleRuntimeVtpWriter expected %d coordinate values, but got %d!",
                statedim * particlestored, positiondata.size());
#endif
        }
        else if (not blackliststates_.count(state))
        {
          // prepare particle state data
          std::vector<double> statedata;
          statedata.reserve(statedim * particlestored);

          // copy particle state data
          for (int i = 0; i < (statedim * particlestored); ++i) statedata.push_back(state_ptr[i]);

          // append particle state data to vtp writer
          (runtime_visualization_managers_[type])[status]
              ->get_visualization_data()
              .SetPointDataVector<double>(statename, statedata, statedim);
        }
      }

      // get pointer to global id of particles
      const int* globalids = (particlestored > 0) ? container->GetPtrToGlobalID(0) : nullptr;

      // prepare particle global id data
      std::vector<double> globaliddata;
      globaliddata.reserve(particlestored);

      // copy particle global id data
      for (int i = 0; i < particlestored; ++i) globaliddata.push_back(globalids[i]);

      // append global id of particles to vtp writer
      (runtime_visualization_managers_[type])[status]
          ->get_visualization_data()
          .SetPointDataVector<double>("globalid", globaliddata, 1);

      // set particle owner data
      std::vector<double> ownerdata(particlestored, comm_.MyPID());

      // append owner of particles to vtp writer
      (runtime_visualization_managers_[type])[status]
          ->get_visualization_data()
          .SetPointDataVector<double>("owner", ownerdata, 1);
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::WriteToDisk(
    const double visualization_time, const int visualization_step)
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_visualization_managers_[type])[status]) continue;

      (runtime_visualization_managers_[type])[status]->WriteToDisk(
          visualization_time, visualization_step);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
