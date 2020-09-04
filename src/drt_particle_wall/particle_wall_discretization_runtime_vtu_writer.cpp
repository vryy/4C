/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_wall_discretization_runtime_vtu_writer.H"

#include "particle_wall_datastate.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/discretization_runtime_vtu_writer.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter()
    : setuptime_(0.0)
{
  // empty constructor
}

PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::~WallDiscretizationRuntimeVtuWriter() = default;

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::Init(
    const Teuchos::RCP<DRT::Discretization> walldiscretization,
    const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate)
{
  // set wall discretization
  walldiscretization_ = walldiscretization;

  // set wall data state container
  walldatastate_ = walldatastate;

  // construct the writer object
  runtime_vtuwriter_ =
      std::unique_ptr<DiscretizationRuntimeVtuWriter>(new DiscretizationRuntimeVtuWriter());
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::Setup(bool write_binary_output)
{
  // we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int max_number_timesteps_to_be_written = 1.0e+6;

  // initialize the writer object
  runtime_vtuwriter_->Initialize(
      walldiscretization_, max_number_timesteps_to_be_written, setuptime_, write_binary_output);
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WriteWallDiscretizationRuntimeOutput(
    const int step, const double time) const
{
  // reset time and time step of the writer object
  runtime_vtuwriter_->ResetTimeAndTimeStep(time, step);

  // node displacements
  {
    if (walldatastate_->GetDispCol() != Teuchos::null)
      runtime_vtuwriter_->AppendDofBasedResultDataVector(
          walldatastate_->GetRefMutableDispCol(), 3, 0, "disp");
  }

  // node owner
  {
    Teuchos::RCP<Epetra_Vector> nodeowner =
        Teuchos::rcp(new Epetra_Vector(*walldiscretization_->NodeColMap(), true));
    for (int inode = 0; inode < walldiscretization_->NumMyColNodes(); ++inode)
    {
      const DRT::Node* node = walldiscretization_->lColNode(inode);
      (*nodeowner)[inode] = node->Owner();
    }
    runtime_vtuwriter_->AppendNodeBasedResultDataVector(nodeowner, 1, "owner");
  }

  // element owner
  {
    runtime_vtuwriter_->AppendElementOwner("owner");
  }

  // element id
  {
    Teuchos::RCP<Epetra_Vector> eleid =
        Teuchos::rcp(new Epetra_Vector(*walldiscretization_->ElementRowMap(), true));
    for (int iele = 0; iele < walldiscretization_->NumMyRowElements(); ++iele)
    {
      const DRT::Element* ele = walldiscretization_->lRowElement(iele);
      (*eleid)[iele] = ele->Id();
    }
    runtime_vtuwriter_->AppendElementBasedResultDataVector(eleid, 1, "id");
  }

  // finalize everything and write all required files to filesystem
  runtime_vtuwriter_->WriteFiles();
  runtime_vtuwriter_->WriteCollectionFileOfAllWrittenFiles();
}
