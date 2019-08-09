/*---------------------------------------------------------------------------*/
/*!
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_wall_discretization_runtime_vtu_writer.H"

#include "particle_wall_datastate.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/discretization_runtime_vtu_writer.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter()
    : setuptime_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::~WallDiscretizationRuntimeVtuWriter()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init wall discretization runtime vtu writer                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | setup wall discretization runtime vtu writer               sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::Setup(bool write_binary_output)
{
  // we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  const unsigned int max_number_timesteps_to_be_written = 1.0e+6;

  // initialize the writer object
  runtime_vtuwriter_->Initialize(
      walldiscretization_, max_number_timesteps_to_be_written, setuptime_, write_binary_output);
}

/*---------------------------------------------------------------------------*
 | write restart of wall discretization runtime vtu writer    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of wall discretization runtime vtu writer     sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // get restart time
  setuptime_ = reader->ReadDouble("time");
}

/*---------------------------------------------------------------------------*
 | write wall discretization runtime output                   sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
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
    Teuchos::RCP<Epetra_Vector> eleowner =
        Teuchos::rcp(new Epetra_Vector(*walldiscretization_->ElementColMap(), true));
    for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
    {
      const DRT::Element* ele = walldiscretization_->lColElement(iele);
      (*eleowner)[iele] = ele->Owner();
    }
    runtime_vtuwriter_->AppendElementBasedResultDataVector(eleowner, 1, "owner");
  }

  // element id
  {
    Teuchos::RCP<Epetra_Vector> eleid =
        Teuchos::rcp(new Epetra_Vector(*walldiscretization_->ElementColMap(), true));
    for (int iele = 0; iele < walldiscretization_->NumMyColElements(); ++iele)
    {
      const DRT::Element* ele = walldiscretization_->lColElement(iele);
      (*eleid)[iele] = ele->Id();
    }
    runtime_vtuwriter_->AppendElementBasedResultDataVector(eleid, 1, "id");
  }

  // finalize everything and write all required files to filesystem
  runtime_vtuwriter_->WriteFiles();
  runtime_vtuwriter_->WriteCollectionFileOfAllWrittenFiles();
}
