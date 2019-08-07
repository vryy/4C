/*---------------------------------------------------------------------------*/
/*!
\brief write visualization output for particle interaction in vtk/vtp format at runtime

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_runtime_vtp_writer.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/runtime_vtp_writer.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::InteractionWriter::InteractionWriter(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), params_(params), setuptime_(0.0), writeresultsthisstep_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::InteractionWriter::~InteractionWriter()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init interaction writer                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup interaction writer                                   sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of interaction writer                        sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of interaction writer                         sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

/*---------------------------------------------------------------------------*
 | register specific runtime vtp writer                       sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::RegisterSpecificRuntimeVtpWriter(
    const std::string& fieldname)
{
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

  // get data format for written numeric data
  bool write_binary_output = (DRT::INPUT::IntegralValue<INPAR::PARTICLE::OutputDataFormat>(
                                  params_, "OUTPUT_DATA_FORMAT") == INPAR::PARTICLE::binary);

  // construct the writer object
  std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter = std::make_shared<RuntimeVtpWriter>();

  // initialize the writer object
  runtime_vtpwriter->Initialize(comm_.MyPID(), comm_.NumProc(), max_number_timesteps_to_be_written,
      output_directory_path, DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      fieldname, DRT::Problem::Instance()->OutputControlFile()->RestartName(), setuptime_,
      write_binary_output);

  // set the writer object
  runtime_vtpwriters_[fieldname] = runtime_vtpwriter;
}

/*---------------------------------------------------------------------------*
 | write particle interaction runtime output                  sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::InteractionWriter::WriteParticleInteractionRuntimeOutput(
    const int step, const double time) const
{
  // iterate over writer objects
  for (auto& writerIt : runtime_vtpwriters_)
  {
    const std::string& fieldname = writerIt.first;
    std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter = writerIt.second;

    // reset time and time step of the writer object
    runtime_vtpwriter->SetupForNewTimeStepAndGeometry(time, step, fieldname);

    // position and states preset in particle interaction evaluation

    // finalize everything and write all required files to filesystem
    runtime_vtpwriter->WriteFiles();
    runtime_vtpwriter->WriteCollectionFileOfAllWrittenFiles(
        DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" + fieldname);
  }
}
