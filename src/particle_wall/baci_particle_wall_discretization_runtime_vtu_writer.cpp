/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_wall_discretization_runtime_vtu_writer.hpp"

#include "baci_global_data.hpp"
#include "baci_io.hpp"
#include "baci_io_control.hpp"
#include "baci_io_discretization_visualization_writer_mesh.hpp"
#include "baci_io_visualization_parameters.hpp"
#include "baci_lib_discret.hpp"
#include "baci_particle_wall_datastate.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter(
    const Teuchos::RCP<DRT::Discretization> walldiscretization,
    const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate, const double restart_time)
    : walldiscretization_(walldiscretization), walldatastate_(walldatastate)
{
  // construct the writer object
  runtime_vtuwriter_ = std::make_unique<IO::DiscretizationVisualizationWriterMesh>(
      walldiscretization, IO::VisualizationParametersFactory(
                              GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                              *GLOBAL::Problem::Instance()->OutputControlFile(), restart_time));
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

FOUR_C_NAMESPACE_CLOSE
