/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particles in vtk/vtp format at runtime
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_runtime_vtp_writer.H"

#include "io.H"
#include "io_control.H"
#include "runtime_vtp_writer.H"

#include "globalproblem.H"
#include "dserror.H"

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

void PARTICLEENGINE::ParticleRuntimeVtpWriter::Setup(
    bool write_binary_output, bool write_ghosted_particles)
{
  // determine size of vector indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  runtime_vtpwriters_.resize(typevectorsize);

  // determine path of output directory
  const std::string outputfilename(DRT::Problem::Instance()->OutputControlFile()->FileName());
  size_t pos = outputfilename.find_last_of("/");
  if (pos == outputfilename.npos)
    pos = 0ul;
  else
    pos++;
  const std::string output_directory_path(outputfilename.substr(0ul, pos));

  // we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int max_number_timesteps_to_be_written = 1.0e+6;

  // construct and initialize all vtp writer objects
  std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter;

  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // allocate memory for vtp writer objects of owned and ghosted states
    (runtime_vtpwriters_[type]).resize(2);

    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      if (status == Ghosted and write_ghosted_particles == false) continue;

      // construct vtp writer object for current particle type and status
      runtime_vtpwriter = std::make_shared<RuntimeVtpWriter>();

      std::ostringstream fieldname;
      fieldname << "particle-" << EnumToTypeName(type) << "-" << EnumToStatusName(status);

      // initialize the writer object
      runtime_vtpwriter->Initialize(comm_.MyPID(), comm_.NumProc(),
          max_number_timesteps_to_be_written, output_directory_path,
          DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(), fieldname.str(),
          DRT::Problem::Instance()->OutputControlFile()->RestartName(), setuptime_,
          write_binary_output);

      // insert into data structure holding all vtp writer objects for each particle type and status
      (runtime_vtpwriters_[type])[status] = runtime_vtpwriter;
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::ResetTimeAndTimeStep(
    double time, unsigned int timestep)
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_vtpwriters_[type])[status]) continue;

      // particle field name
      std::ostringstream particlefieldname;
      particlefieldname << "particle-" << EnumToTypeName(type) << "-" << EnumToStatusName(status);

      (runtime_vtpwriters_[type])[status]->SetupForNewTimeStepAndGeometry(
          time, timestep, particlefieldname.str());
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::SetParticlePositionsAndStates()
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_vtpwriters_[type])[status]) continue;

      // get container of current particle type and status
      ParticleContainer* container = particlecontainerbundle_->GetSpecificContainer(type, status);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // get particle states stored in container
      const std::set<ParticleState>& states = container->GetStoredStates();

#ifdef DEBUG
      // safety check
      if (not container->HaveStoredState(Position))
        dserror("particle state '%s' not found!", EnumToStateName(Position).c_str());
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
          std::vector<double>& positiondata =
              (runtime_vtpwriters_[type])[status]->GetMutablePointCoordinateVector();
          positiondata.clear();
          positiondata.reserve(statedim * particlestored);

          // copy particle position data
          for (int i = 0; i < (statedim * particlestored); ++i)
            positiondata.push_back(state_ptr[i]);

#ifdef DEBUG
          // safety check
          if (static_cast<int>(positiondata.size()) != statedim * particlestored)
            dserror("ParticleRuntimeVtpWriter expected %d coordinate values, but got %d!",
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
          (runtime_vtpwriters_[type])[status]->AppendVisualizationPointDataVector(
              statedata, statedim, statename);
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
      (runtime_vtpwriters_[type])[status]->AppendVisualizationPointDataVector(
          globaliddata, 1, "globalid");

      // set particle owner data
      std::vector<double> ownerdata(particlestored, comm_.MyPID());

      // append owner of particles to vtp writer
      (runtime_vtpwriters_[type])[status]->AppendVisualizationPointDataVector(
          ownerdata, 1, "owner");
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::WriteFiles()
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_vtpwriters_[type])[status]) continue;

      (runtime_vtpwriters_[type])[status]->WriteFiles();
    }
  }
}

void PARTICLEENGINE::ParticleRuntimeVtpWriter::WriteCollectionFileOfAllWrittenFiles()
{
  // iterate over particle types
  for (const auto& type : particlecontainerbundle_->GetParticleTypes())
  {
    // iterate over particle statuses
    for (const auto& status : {Owned, Ghosted})
    {
      // check for runtime vtp writer for current particle type and status
      if (not(runtime_vtpwriters_[type])[status]) continue;

      // particle field name
      std::ostringstream particlefieldname;
      particlefieldname << "particle-" << EnumToTypeName(type) << "-" << EnumToStatusName(status);

      (runtime_vtpwriters_[type])[status]->WriteCollectionFileOfAllWrittenFiles(
          DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" +
          particlefieldname.str());
    }
  }
}
