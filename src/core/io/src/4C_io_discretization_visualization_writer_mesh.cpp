/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization, i.e., write the mesh and results on the mesh
to disk

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "4C_io_discretization_visualization_writer_mesh.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_io.hpp"
#include "4C_io_element_vtk_cell_type_register.hpp"
#include "4C_io_visualization_manager.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_FEVector.h>

#include <utility>

FOUR_C_NAMESPACE_OPEN


namespace Core::IO
{
  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  DiscretizationVisualizationWriterMesh::DiscretizationVisualizationWriterMesh(
      const Teuchos::RCP<const Core::FE::Discretization>& discretization,
      VisualizationParameters parameters,
      std::function<bool(const Core::Elements::Element* element)> element_filter)
      : discretization_(discretization),
        visualization_manager_(Teuchos::rcp(new VisualizationManager(
            std::move(parameters), discretization->get_comm(), discretization->name()))),
        element_filter_(std::move(element_filter))
  {
    set_geometry_from_discretization();
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::set_geometry_from_discretization()
  {
    // Todo assume 3D for now
    const unsigned int num_spatial_dimensions = 3;

    // count number of nodes; output is completely independent of the number of processors involved
    unsigned int num_row_elements = discretization_->num_my_row_elements();
    unsigned int num_nodes = 0;
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      num_nodes += ele->num_node();
    }

    // do not need to store connectivity indices here because we create a
    // contiguous array by the order in which we fill the coordinates (otherwise
    // need to adjust order of filling in the coordinates).
    auto& visualization_data = visualization_manager_->get_visualization_data();

    std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
    point_coordinates.clear();
    point_coordinates.reserve(num_spatial_dimensions * num_nodes);

    std::vector<uint8_t>& cell_types = visualization_data.get_cell_types();
    cell_types.clear();
    cell_types.reserve(num_row_elements);

    std::vector<int32_t>& cell_offsets = visualization_data.get_cell_offsets();
    cell_offsets.clear();
    cell_offsets.reserve(num_row_elements);


    // loop over my elements and collect the geometry/grid data
    unsigned int pointcounter = 0;
    unsigned int num_skipped_eles = 0;

    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (element_filter_(ele))
      {
        pointcounter +=
            ele->append_visualization_geometry(*discretization_, cell_types, point_coordinates);
        cell_offsets.push_back(pointcounter);
      }
      else
        ++num_skipped_eles;
    }


    // safety checks
    if (point_coordinates.size() != num_spatial_dimensions * pointcounter)
    {
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh expected %d coordinate values, but got %d",
          num_spatial_dimensions * pointcounter, point_coordinates.size());
    }

    if (cell_types.size() != num_row_elements - num_skipped_eles)
    {
      FOUR_C_THROW("DiscretizationVisualizationWriterMesh expected %d cell type values, but got %d",
          num_row_elements, cell_types.size());
    }

    if (cell_offsets.size() != num_row_elements - num_skipped_eles)
    {
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh expected %d cell offset values, but got %d",
          num_row_elements, cell_offsets.size());
    }

    // store node row and col maps (needed to check for changed parallel distribution)
    noderowmap_last_geometry_set_ = Teuchos::rcp(new Epetra_Map(*discretization_->node_row_map()));
    nodecolmap_last_geometry_set_ = Teuchos::rcp(new Epetra_Map(*discretization_->node_col_map()));
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::reset()
  {
    // check if parallel distribution of discretization changed
    int map_changed =
        ((not noderowmap_last_geometry_set_->SameAs(*discretization_->node_row_map())) or
            (not nodecolmap_last_geometry_set_->SameAs(*discretization_->node_col_map())));
    int map_changed_allproc(0);
    discretization_->get_comm().MaxAll(&map_changed, &map_changed_allproc, 1);

    // reset geometry of visualization writer
    if (map_changed_allproc) set_geometry_from_discretization();
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_dof_based_result_data_vector(
      const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int result_num_dofs_per_node, unsigned int read_result_data_from_dofindex,
      const std::string& resultname)
  {
    /* the idea is to transform the given data to a 'point data vector' and append it to the
     * collected solution data vectors by calling AppendVisualizationPointDataVector() */

    // safety checks
    if (!discretization_->dof_col_map()->SameAs(result_data_dofbased->Map()))
    {
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh: Received DofBasedResult's map does not match the "
          "discretization's dof col map.");
    }

    // count number of nodes for this visualization
    unsigned int num_nodes = 0;
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      num_nodes += ele->num_node();
    }

    std::vector<double> point_result_data;
    point_result_data.reserve(result_num_dofs_per_node * num_nodes);

    unsigned int pointcounter = 0;

    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (element_filter_(ele))
      {
        pointcounter += ele->append_visualization_dof_based_result_data_vector(*discretization_,
            result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
            point_result_data);
      }
    }

    // sanity check
    if (point_result_data.size() != result_num_dofs_per_node * pointcounter)
    {
      FOUR_C_THROW("DiscretizationVisualizationWriterMesh expected %d result values, but got %d",
          result_num_dofs_per_node * pointcounter, point_result_data.size());
    }

    visualization_manager_->get_visualization_data().set_point_data_vector(
        resultname, point_result_data, result_num_dofs_per_node);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_node_based_result_data_vector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
      unsigned int result_num_components_per_node, const std::string& resultname)
  {
    /* the idea is to transform the given data to a 'point data vector' and append it to the
     * collected solution data vectors by calling AppendVisualizationPointDataVector() */

    // safety checks
    if ((unsigned int)result_data_nodebased->NumVectors() != result_num_components_per_node)
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh: expected Epetra_MultiVector with %d columns but "
          "got %d",
          result_num_components_per_node, result_data_nodebased->NumVectors());

    if (!discretization_->node_col_map()->SameAs(result_data_nodebased->Map()))
    {
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh: Received NodeBasedResult's map does not match "
          "the "
          "discretization's node col map.");
    }

    // count number of nodes
    unsigned int num_nodes = 0;
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      num_nodes += ele->num_node();
    }

    std::vector<double> point_result_data;
    point_result_data.reserve(result_num_components_per_node * num_nodes);

    unsigned int pointcounter = 0;

    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (!element_filter_(ele)) continue;

      const std::vector<int>& numbering =
          Core::IO::GetVtkCellTypeFromFourCElementShapeType(ele->shape()).second;

      for (unsigned int inode = 0; inode < (unsigned int)ele->num_node(); ++inode)
      {
        const Core::Nodes::Node* node = ele->nodes()[numbering[inode]];

        const int lid = node->lid();

        for (unsigned int icpn = 0; icpn < result_num_components_per_node; ++icpn)
        {
          Epetra_Vector* column = (*result_data_nodebased)(icpn);

          if (lid > -1)
            point_result_data.push_back((*column)[lid]);
          else
            FOUR_C_THROW("received illegal node local id: %d", lid);
        }
      }

      pointcounter += ele->num_node();
    }

    // sanity check
    if (point_result_data.size() != result_num_components_per_node * pointcounter)
    {
      FOUR_C_THROW("DiscretizationVisualizationWriterMesh expected %d result values, but got %d",
          result_num_components_per_node * pointcounter, point_result_data.size());
    }

    visualization_manager_->get_visualization_data().set_point_data_vector<double>(
        resultname, point_result_data, result_num_components_per_node);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_element_based_result_data_vector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_elementbased,
      unsigned int result_num_components_per_element, const std::string& resultname)
  {
    /* the idea is to transform the given data to a 'cell data vector' and append it to the
     *  collected solution data vectors by calling AppendVisualizationCellDataVector() */

    // safety check
    if ((unsigned int)result_data_elementbased->NumVectors() != result_num_components_per_element)
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh: expected Epetra_MultiVector with %d columns but "
          "got %d",
          result_num_components_per_element, result_data_elementbased->NumVectors());

    if (!discretization_->element_row_map()->SameAs(result_data_elementbased->Map()))
    {
      FOUR_C_THROW(
          "DiscretizationVisualizationWriterMesh: Received ElementBasedResult's map does not match "
          "the "
          "discretization's element row map.");
    }

    // count number of elements for each processor
    unsigned int num_row_elements = (unsigned int)discretization_->num_my_row_elements();

    std::vector<double> cell_result_data;
    cell_result_data.reserve(result_num_components_per_element * num_row_elements);

    unsigned int cellcounter = 0;

    for (unsigned int iele = 0; iele < num_row_elements; ++iele)
    {
      const Core::Elements::Element* ele = discretization_->l_row_element(iele);

      if (!element_filter_(ele)) continue;

      for (unsigned int icpe = 0; icpe < result_num_components_per_element; ++icpe)
      {
        Epetra_Vector* column = (*result_data_elementbased)(icpe);

        cell_result_data.push_back((*column)[iele]);
      }

      ++cellcounter;
    }

    // sanity check
    if (cell_result_data.size() != result_num_components_per_element * cellcounter)
    {
      FOUR_C_THROW("DiscretizationVisualizationWriterMesh expected %d result values, but got %d",
          result_num_components_per_element * cellcounter, cell_result_data.size());
    }

    visualization_manager_->get_visualization_data().set_cell_data_vector(
        resultname, cell_result_data, result_num_components_per_element);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_element_owner(const std::string resultname)
  {
    // Vector with element owner for elements in the row map.
    std::vector<double> owner_of_row_elements;
    owner_of_row_elements.reserve(discretization_->num_my_row_elements());

    const int my_pid = discretization_->get_comm().MyPID();
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (element_filter_(ele)) owner_of_row_elements.push_back(my_pid);
    }

    // Pass data to the output writer.
    visualization_manager_->get_visualization_data().set_cell_data_vector(
        resultname, owner_of_row_elements, 1);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_element_gid(const std::string& resultname)
  {
    // Vector with element IDs for elements in the row map.
    std::vector<double> gid_of_row_elements;
    gid_of_row_elements.reserve(discretization_->num_my_row_elements());

    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (element_filter_(ele)) gid_of_row_elements.push_back(ele->id());
    }

    // Pass data to the output writer.
    visualization_manager_->get_visualization_data().set_cell_data_vector(
        resultname, gid_of_row_elements, 1);
  }


  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_element_ghosting_information()
  {
    Core::IO::append_element_ghosting_information(
        *discretization_, *visualization_manager_, element_filter_);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::append_node_gid(const std::string& resultname)
  {
    // count number of nodes; output is completely independent of the number of processors involved
    int num_nodes = 0;
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      num_nodes += ele->num_node();
    }

    // Setup the vector with the GIDs of the nodes.
    std::vector<double> gid_of_nodes;
    gid_of_nodes.reserve(num_nodes);

    // Loop over each element and add the node GIDs.
    for (const Core::Elements::Element* ele : discretization_->my_row_element_range())
    {
      if (!element_filter_(ele)) continue;

      // Add the node GIDs.
      const std::vector<int>& numbering =
          Core::IO::GetVtkCellTypeFromFourCElementShapeType(ele->shape()).second;
      const Core::Nodes::Node* const* nodes = ele->nodes();
      for (int inode = 0; inode < ele->num_node(); ++inode)
        gid_of_nodes.push_back(nodes[numbering[inode]]->id());
    }

    visualization_manager_->get_visualization_data().set_point_data_vector<double>(
        resultname, gid_of_nodes, 1);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void DiscretizationVisualizationWriterMesh::write_to_disk(
      const double visualization_time, const int visualization_step)
  {
    visualization_manager_->write_to_disk(visualization_time, visualization_step);
  }

  /*-----------------------------------------------------------------------------------------------*
   *-----------------------------------------------------------------------------------------------*/
  void append_element_ghosting_information(const Core::FE::Discretization& discretization,
      VisualizationManager& visualization_manager,
      const std::function<bool(const Core::Elements::Element* ele)>& element_predicate)
  {
    // Set up a multivector which will be populated with all ghosting informations.
    const Epetra_Comm& comm = discretization.element_col_map()->Comm();
    const int n_proc = comm.NumProc();
    const int my_proc = comm.MyPID();

    // Create Vectors to store the ghosting information.
    Epetra_FEVector ghosting_information(*discretization.element_row_map(), n_proc);

    // Get elements ghosted by this rank.
    std::vector<int> my_ghost_elements;
    my_ghost_elements.clear();
    int count = 0;
    for (const Core::Elements::Element* ele : discretization.my_col_element_range())
    {
      if (element_predicate(ele))
      {
        count++;
        if (ele->owner() != my_proc) my_ghost_elements.push_back(ele->id());
      }
    }

    // Add to the multi vector.
    std::vector<double> values(my_ghost_elements.size(), 1.0);
    ghosting_information.SumIntoGlobalValues(
        my_ghost_elements.size(), my_ghost_elements.data(), values.data(), my_proc);

    // Assemble over all processors.
    ghosting_information.GlobalAssemble();

    // Output the ghosting data of the elements owned by this proc.
    std::vector<double> ghosted_elements;
    ghosted_elements.reserve(count * n_proc);
    for (const Core::Elements::Element* ele : discretization.my_row_element_range())
    {
      if (element_predicate(ele))
      {
        const int local_row = ghosting_information.Map().LID(ele->id());
        if (local_row == -1) FOUR_C_THROW("The element has to exist in the row map.");
        for (int i_proc = 0; i_proc < n_proc; i_proc++)
          ghosted_elements.push_back(ghosting_information[i_proc][local_row]);
      }
    }

    visualization_manager.get_visualization_data().set_cell_data_vector(
        "element_ghosting", ghosted_elements, n_proc);
  }
}  // namespace Core::IO
FOUR_C_NAMESPACE_CLOSE
