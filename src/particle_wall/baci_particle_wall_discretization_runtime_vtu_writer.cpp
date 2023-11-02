/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_wall_discretization_runtime_vtu_writer.H"

#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_discretization_visualization_writer_mesh.H"
#include "baci_io_visualization_parameters.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_particle_wall_datastate.H"

#include <memory>


/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter()
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
  runtime_vtuwriter_ = std::make_unique<IO::DiscretizationVisualizationWriterMesh>(
      walldiscretization_, IO::VisualizationParametersFactory(
                               DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT")));
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WriteWallDiscretizationRuntimeOutput(
    const int step, const double time) const
{
  // reset the writer object
  runtime_vtuwriter_->Reset();

  // node displacements
  {
    if (walldatastate_->GetDispCol() != Teuchos::null)
      runtime_vtuwriter_->AppendDofBasedResultDataVector(
          walldatastate_->GetRefDispCol(), 3, 0, "disp");
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
  runtime_vtuwriter_->WriteToDisk(time, step);
}
