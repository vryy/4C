/*---------------------------------------------------------------------------*/
/*! \file
\brief write output for particle interaction at runtime
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_runtime_writer.H"

#include "baci_inpar_particle.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_runtime_csv_writer.H"
#include "baci_io_visualization_manager.H"
#include "baci_lib_globalproblem.H"

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

void PARTICLEINTERACTION::InteractionWriter::RegisterSpecificRuntimeVtuWriter(
    const std::string& fieldname)
{
  // safety check
  if (runtime_visualization_managers_.count(fieldname))
    dserror("a runtime vtu writer for field '%s' is already stored!", fieldname.c_str());

  // construct and init the vtp writer object
  std::shared_ptr<IO::VisualizationManager> runtime_visualization_manager =
      std::make_shared<IO::VisualizationManager>(
          IO::VisualizationParametersFactory(
              DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT")),
          comm_, fieldname);

  // set the vtp writer object
  runtime_visualization_managers_[fieldname] = runtime_visualization_manager;
}

void PARTICLEINTERACTION::InteractionWriter::RegisterSpecificRuntimeCsvWriter(
    const std::string& fieldname)
{
  // safety check
  if (runtime_csvwriters_.count(fieldname))
    dserror("a runtime csv writer for field '%s' is already stored!", fieldname.c_str());

  // set the csv writer object
  runtime_csvwriters_[fieldname] = std::make_shared<IO::RuntimeCsvWriter>(comm_.MyPID(), fieldname);
}

void PARTICLEINTERACTION::InteractionWriter::WriteParticleInteractionRuntimeOutput(
    const int step, const double time) const
{
  // iterate over vtp writer objects
  for (auto& writerIt : runtime_visualization_managers_)
  {
    std::shared_ptr<IO::VisualizationManager> runtime_visualization_manager = writerIt.second;

    // data to be written preset in particle interaction evaluation

    // finalize everything and write all required files to filesystem
    runtime_visualization_manager->WriteToDisk(time, step);
  }

  // iterate over csv writer objects
  for (auto& writerIt : runtime_csvwriters_)
  {
    std::shared_ptr<IO::RuntimeCsvWriter> runtime_csvwriter = writerIt.second;

    // reset time and time step of the writer object
    runtime_csvwriter->ResetTimeAndTimeStep(time, step);

    // data to be written preset in particle interaction evaluation

    // write file to filesystem
    runtime_csvwriter->WriteFile();
  }
}
