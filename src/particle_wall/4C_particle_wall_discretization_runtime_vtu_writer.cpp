/*---------------------------------------------------------------------------*/
/*! \file
\brief write visualization output for particle wall discretization in vtk/vtu format at runtime
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_wall_discretization_runtime_vtu_writer.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_discretization_visualization_writer_mesh.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_particle_wall_datastate.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::WallDiscretizationRuntimeVtuWriter(
    const Teuchos::RCP<Discret::Discretization> walldiscretization,
    const std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate, const double restart_time)
    : walldiscretization_(walldiscretization), walldatastate_(walldatastate)
{
  // construct the writer object
  runtime_vtuwriter_ = std::make_unique<Core::IO::DiscretizationVisualizationWriterMesh>(
      walldiscretization, Core::IO::VisualizationParametersFactory(
                              Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
                              *Global::Problem::Instance()->OutputControlFile(), restart_time));
}

void PARTICLEWALL::WallDiscretizationRuntimeVtuWriter::write_wall_discretization_runtime_output(
    const int step, const double time) const
{
  // reset the writer object
  runtime_vtuwriter_->Reset();

  // node displacements
  {
    if (walldatastate_->GetDispCol() != Teuchos::null)
      runtime_vtuwriter_->append_dof_based_result_data_vector(
          walldatastate_->GetRefDispCol(), 3, 0, "disp");
  }

  // node owner
  {
    Teuchos::RCP<Epetra_Vector> nodeowner =
        Teuchos::rcp(new Epetra_Vector(*walldiscretization_->NodeColMap(), true));
    for (int inode = 0; inode < walldiscretization_->NumMyColNodes(); ++inode)
    {
      const Core::Nodes::Node* node = walldiscretization_->lColNode(inode);
      (*nodeowner)[inode] = node->Owner();
    }
    runtime_vtuwriter_->append_node_based_result_data_vector(nodeowner, 1, "owner");
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
      const Core::Elements::Element* ele = walldiscretization_->lRowElement(iele);
      (*eleid)[iele] = ele->Id();
    }
    runtime_vtuwriter_->append_element_based_result_data_vector(eleid, 1, "id");
  }

  // finalize everything and write all required files to filesystem
  runtime_vtuwriter_->WriteToDisk(time, step);
}

FOUR_C_NAMESPACE_CLOSE
