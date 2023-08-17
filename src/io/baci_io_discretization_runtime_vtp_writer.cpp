/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization in vtk/vtp format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "baci_io_discretization_runtime_vtp_writer.H"

#include "baci_io_control.H"
#include "baci_io_visualization_manager.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_node.H"
#include "baci_utils_exceptions.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DiscretizationRuntimeVtpWriter::DiscretizationRuntimeVtpWriter(
    const Teuchos::RCP<const DRT::Discretization>& discretization,
    IO::VisualizationParameters parameters)
    : discretization_(discretization),
      visualization_manager_(Teuchos::rcp(
          new IO::VisualizationManager(parameters, discretization->Comm(), discretization->Name())))
{
}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DiscretizationRuntimeVtpWriter::SetGeometryFromParticleDiscretization()
{
  // Note: this will only work for 'particle' discretizations that have nodes but no elements


  // Todo assume 3D for now
  const unsigned int num_spatial_dimensions = 3;

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  unsigned int num_row_nodes = (unsigned int)discretization_->NumMyRowNodes();

  // get and prepare storage for point coordinate values
  std::vector<double>& point_coordinates =
      visualization_manager_->GetVisualizationDataMutable().GetPointCoordinatesMutable();
  point_coordinates.clear();
  point_coordinates.reserve(num_spatial_dimensions * num_row_nodes);


  // loop over my nodes and collect the geometry/grid data, i.e. reference positions of nodes
  for (unsigned int inode = 0; inode < num_row_nodes; ++inode)
  {
    const DRT::Node* node = discretization_->lRowNode(inode);

    for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
      point_coordinates.push_back(node->X()[idim]);
  }


  // safety check
  if (point_coordinates.size() != num_spatial_dimensions * num_row_nodes)
  {
    dserror("DiscretizationRuntimeVtpWriter expected %d coordinate values, but got %d",
        num_spatial_dimensions * num_row_nodes, point_coordinates.size());
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DiscretizationRuntimeVtpWriter::AppendDofBasedResultDataVector(
    const Teuchos::RCP<Epetra_Vector>& result_data_dofbased, unsigned int result_num_dofs_per_node,
    const std::string& resultname)
{
  /* the idea is to transform the given data to a 'point data vector' and append it to the
   * collected solution data vectors by calling AppendVisualizationPointDataVector() */

  std::vector<double> vtp_point_result_data;
  vtp_point_result_data.reserve(result_data_dofbased->MyLength());

  for (int lid = 0; lid < result_data_dofbased->MyLength(); ++lid)
    vtp_point_result_data.push_back((*result_data_dofbased)[lid]);

  visualization_manager_->GetVisualizationDataMutable().SetPointDataVector<double>(
      resultname, vtp_point_result_data, result_num_dofs_per_node);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DiscretizationRuntimeVtpWriter::AppendNodeBasedResultDataVector(
    const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
    unsigned int result_num_components_per_node, const std::string& resultname)
{
  /*  the idea is to transform the given data to a 'point data vector' and append it to the
   *  collected solution data vectors by calling AppendVisualizationPointDataVector()
   */

  // count number of nodes for each processor
  const unsigned int num_row_nodes = (unsigned int)result_data_nodebased->Map().NumMyElements();

  // safety check
  if ((unsigned int)result_data_nodebased->NumVectors() != result_num_components_per_node)
    dserror(
        "DiscretizationRuntimeVtpWriter: expected Epetra_MultiVector with %d columns but got %d",
        result_num_components_per_node, result_data_nodebased->NumVectors());


  std::vector<double> vtp_point_result_data;
  vtp_point_result_data.reserve(result_num_components_per_node * num_row_nodes);

  for (unsigned int lid = 0; lid < num_row_nodes; ++lid)
  {
    for (unsigned int idf = 0; idf < result_num_components_per_node; ++idf)
    {
      Epetra_Vector* column = (*result_data_nodebased)(idf);
      vtp_point_result_data.push_back((*column)[lid]);
    }
  }

  visualization_manager_->GetVisualizationDataMutable().SetPointDataVector<double>(
      resultname, vtp_point_result_data, result_num_components_per_node);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DiscretizationRuntimeVtpWriter::WriteToDisk(
    const double visualziation_time, const int visualization_step)
{
  visualization_manager_->WriteToDisk(visualziation_time, visualization_step);
}
