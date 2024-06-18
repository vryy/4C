/*---------------------------------------------------------------------------*/
/*! \file
\brief write output for particle interaction at runtime
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_runtime_writer.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_particle.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_io_visualization_manager.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::InteractionWriter::InteractionWriter(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), setuptime_(0.0), writeresultsthisstep_(true)
{
  // empty constructor
}

void ParticleInteraction::InteractionWriter::Init()
{
  // nothing to do
}

void ParticleInteraction::InteractionWriter::setup()
{
  // nothing to do
}

void ParticleInteraction::InteractionWriter::read_restart(
    const std::shared_ptr<Core::IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->read_double("time");
}

void ParticleInteraction::InteractionWriter::register_specific_runtime_output_writer(
    const std::string& fieldname)
{
  // safety check
  if (runtime_visualization_managers_.count(fieldname))
    FOUR_C_THROW("a runtime output writer for field '%s' is already stored!", fieldname.c_str());

  // construct and init the output writer object
  std::shared_ptr<Core::IO::VisualizationManager> runtime_visualization_manager =
      std::make_shared<Core::IO::VisualizationManager>(
          Core::IO::VisualizationParametersFactory(
              Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
              *Global::Problem::Instance()->OutputControlFile(), setuptime_),
          comm_, fieldname);

  // set the output writer object
  runtime_visualization_managers_[fieldname] = runtime_visualization_manager;
}

void ParticleInteraction::InteractionWriter::register_specific_runtime_csv_writer(
    const std::string& fieldname)
{
  // safety check
  if (runtime_csvwriters_.count(fieldname))
    FOUR_C_THROW("a runtime csv writer for field '%s' is already stored!", fieldname.c_str());

  // set the csv writer object
  runtime_csvwriters_[fieldname] = std::make_shared<Core::IO::RuntimeCsvWriter>(
      comm_.MyPID(), *Global::Problem::Instance()->OutputControlFile(), fieldname);
}

void ParticleInteraction::InteractionWriter::write_particle_interaction_runtime_output(
    const int step, const double time) const
{
  // iterate over output writer objects
  for (auto& writerIt : runtime_visualization_managers_)
  {
    std::shared_ptr<Core::IO::VisualizationManager> runtime_visualization_manager = writerIt.second;

    // data to be written preset in particle interaction evaluation

    // finalize everything and write all required files to filesystem
    runtime_visualization_manager->WriteToDisk(time, step);
  }

  // iterate over csv writer objects
  for (auto& writerIt : runtime_csvwriters_)
  {
    std::shared_ptr<Core::IO::RuntimeCsvWriter> runtime_csvwriter = writerIt.second;

    // reset time and time step of the writer object
    runtime_csvwriter->reset_time_and_time_step(time, step);

    // data to be written preset in particle interaction evaluation

    // write file to filesystem
    runtime_csvwriter->write_collected_data_to_file();
  }
}

FOUR_C_NAMESPACE_CLOSE
