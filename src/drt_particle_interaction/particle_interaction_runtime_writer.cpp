/*---------------------------------------------------------------------------*/
/*! \file
\brief write output for particle interaction at runtime
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_runtime_writer.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/runtime_vtp_writer.H"
#include "../drt_io/runtime_csv_writer.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::InteractionWriter::InteractionWriter(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), params_(params), setuptime_(0.0), writeresultsthisstep_(true)
{
  // empty constructor
}

void PARTICLEINTERACTION::InteractionWriter::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::InteractionWriter::Setup()
{
  // nothing to do
}

void PARTICLEINTERACTION::InteractionWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

void PARTICLEINTERACTION::InteractionWriter::RegisterSpecificRuntimeVtpWriter(
    const std::string& fieldname)
{
  // safety check
  if (runtime_vtpwriters_.count(fieldname))
    dserror("a runtime vtp writer for field '%s' is already stored!", fieldname.c_str());

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

  // construct and init the vtp writer object
  std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter = std::make_shared<RuntimeVtpWriter>();

  runtime_vtpwriter->Initialize(comm_.MyPID(), comm_.NumProc(), max_number_timesteps_to_be_written,
      output_directory_path, DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(),
      fieldname, DRT::Problem::Instance()->OutputControlFile()->RestartName(), setuptime_,
      write_binary_output);

  // set the vtp writer object
  runtime_vtpwriters_[fieldname] = runtime_vtpwriter;
}

void PARTICLEINTERACTION::InteractionWriter::RegisterSpecificRuntimeCsvWriter(
    const std::string& fieldname)
{
  // safety check
  if (runtime_csvwriters_.count(fieldname))
    dserror("a runtime csv writer for field '%s' is already stored!", fieldname.c_str());

  // construct and init the csv writer object
  std::shared_ptr<RuntimeCsvWriter> runtime_csvwriter =
      std::make_shared<RuntimeCsvWriter>(comm_.MyPID());

  runtime_csvwriter->Init(fieldname);

  // set the csv writer object
  runtime_csvwriters_[fieldname] = runtime_csvwriter;
}

void PARTICLEINTERACTION::InteractionWriter::WriteParticleInteractionRuntimeOutput(
    const int step, const double time) const
{
  // iterate over vtp writer objects
  for (auto& writerIt : runtime_vtpwriters_)
  {
    const std::string& fieldname = writerIt.first;
    std::shared_ptr<RuntimeVtpWriter> runtime_vtpwriter = writerIt.second;

    // reset time and time step of the writer object
    runtime_vtpwriter->SetupForNewTimeStepAndGeometry(time, step, fieldname);

    // data to be written preset in particle interaction evaluation

    // finalize everything and write all required files to filesystem
    runtime_vtpwriter->WriteFiles();
    runtime_vtpwriter->WriteCollectionFileOfAllWrittenFiles(
        DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" + fieldname);
  }

  // iterate over csv writer objects
  for (auto& writerIt : runtime_csvwriters_)
  {
    std::shared_ptr<RuntimeCsvWriter> runtime_csvwriter = writerIt.second;

    // reset time and time step of the writer object
    runtime_csvwriter->ResetTimeAndTimeStep(time, step);

    // data to be written preset in particle interaction evaluation

    // write file to filesystem
    runtime_csvwriter->WriteFile();
  }
}
