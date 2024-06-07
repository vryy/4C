/*---------------------------------------------------------------------*/
/*! \file

\brief Functions to append visualization for runtime output

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_IO_ELEMENT_APPEND_VISUALIZATION_HPP
#define FOUR_C_IO_ELEMENT_APPEND_VISUALIZATION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_io_element_vtk_cell_type_register.hpp"
#include "4C_nurbs_discret_nurbs_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * \brief Add the element geometry visualization for elements that use Lagrange shape functions
   *
   * @param ele (in) Element
   * @param discret (in) discretization
   * @param cell_types (in/out) cell type data vector
   * @param point_coordinates (in/out) point coordinates for the representation of this element
   * @return Number of added points
   */
  unsigned int AppendVisualizationGeometryLagrangeEle(const Core::Elements::Element& ele,
      const Discret::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    const unsigned int num_spatial_dimensions = 3;
    auto vtk_cell_info = GetVtkCellTypeFromFourCElementShapeType(ele.Shape());
    const std::vector<int>& numbering = vtk_cell_info.second;

    // Add the cell type to the output.
    cell_types.push_back(vtk_cell_info.first);

    // Add each node to the output.
    const Core::Nodes::Node* const* nodes = ele.Nodes();

    for (int inode = 0; inode < ele.num_node(); ++inode)
      for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        point_coordinates.push_back(nodes[numbering[inode]]->X()[idim]);

    // Return the number of added points.
    return ele.num_node();
  }

  /**
   * \brief Helper function to return the proper VTK cell type for nurbs elements.
   */
  template <Core::FE::CellType celltype>
  auto GetVTKCellTypeForNURBSElements()
  {
    switch (celltype)
    {
      case Core::FE::CellType::nurbs9:
        return GetVtkCellTypeFromFourCElementShapeType(Core::FE::CellType::quad9);
        break;
      case Core::FE::CellType::nurbs27:
        return GetVtkCellTypeFromFourCElementShapeType(Core::FE::CellType::hex27);
        break;
      default:
        FOUR_C_THROW("The VTK cell type for the NURBS element %s is not implemented",
            Core::FE::CellTypeToString(celltype).c_str());
        break;
    }
  }

  /**
   * \brief Helper function to evaluate the NURBS interpolation inside the element.
   */
  template <unsigned int n_points, unsigned int n_dim_nurbs, unsigned int n_dim = n_dim_nurbs>
  Core::LinAlg::Matrix<n_dim, 1, double> EvalNurbsInterpolation(
      const Core::LinAlg::Matrix<n_points * n_dim, 1, double>& controlpoint_data,
      const Core::LinAlg::Matrix<n_dim_nurbs, 1, double>& xi,
      const Core::LinAlg::Matrix<n_points, 1, double>& weights,
      const std::vector<Core::LinAlg::SerialDenseVector>& knots, const Core::FE::CellType& distype)
  {
    Core::LinAlg::Matrix<n_dim, 1, double> point_result;

    // Get the shape functions.
    Core::LinAlg::Matrix<n_points, 1, double> N;

    if (n_dim_nurbs == 3)
      Core::FE::Nurbs::nurbs_get_3D_funct(N, xi, knots, weights, distype);
    else if (n_dim_nurbs == 2)
      Core::FE::Nurbs::nurbs_get_2D_funct(N, xi, knots, weights, distype);
    else
      FOUR_C_THROW("Unable to compute the shape functions for this nurbs element case");

    for (unsigned int i_node_nurbs = 0; i_node_nurbs < n_points; i_node_nurbs++)
    {
      for (unsigned int i_dim = 0; i_dim < n_dim; i_dim++)
        point_result(i_dim) += N(i_node_nurbs) * controlpoint_data(i_node_nurbs * 3 + i_dim);
    }

    return point_result;
  }

  /**
   * \brief Add the element geometry visualization for elements that use NURBS shape functions
   *
   * @param ele (in) Element
   * @param discret (in) discretization
   * @param cell_types (in/out) cell type data vector
   * @param point_coordinates (in/out) point coordinates for the representation of this element
   * @return Number of added points
   */
  template <Core::FE::CellType celltype>
  unsigned int AppendVisualizationGeometryNURBS(const Core::Elements::Element& ele,
      const Discret::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    constexpr int number_of_output_points = Core::FE::num_nodes<celltype>;
    constexpr int dim_nurbs = Core::FE::dim<celltype>;
    constexpr int dim_output = 3;

    // Get the vtk cell information
    const auto vtk_cell_info = GetVTKCellTypeForNURBSElements<celltype>();
    const std::vector<int>& numbering = vtk_cell_info.second;

    // Add the cell type to the output.
    cell_types.push_back(vtk_cell_info.first);

    // Create the vertices for the visualization.
    {
      // Get the knots and weights of this element.
      Core::LinAlg::Matrix<number_of_output_points, 1, double> weights(true);
      std::vector<Core::LinAlg::SerialDenseVector> knots(true);
      const bool zero_size =
          Discret::Nurbs::GetMyNurbsKnotsAndWeights(discret, &ele, knots, weights);
      if (zero_size) FOUR_C_THROW("GetMyNurbsKnotsAndWeights has to return a non zero size.");

      // Get the position of the control points in the reference configuration.
      Core::LinAlg::Matrix<number_of_output_points * dim_output, 1, double> pos_controlpoints;
      for (unsigned int i_controlpoint = 0; i_controlpoint < (unsigned int)ele.num_node();
           ++i_controlpoint)
      {
        const Core::Nodes::Node* controlpoint = ele.Nodes()[i_controlpoint];
        for (int i_dim = 0; i_dim < dim_output; ++i_dim)
          pos_controlpoints(dim_output * i_controlpoint + i_dim) = controlpoint->X()[i_dim];
      }

      Core::LinAlg::Matrix<dim_output, 1, double> point_result;
      Core::LinAlg::Matrix<dim_nurbs, 1, double> xi;
      for (unsigned int i_node_nurbs = 0; i_node_nurbs < number_of_output_points; i_node_nurbs++)
      {
        for (unsigned int i_dim_nurbs = 0; i_dim_nurbs < dim_nurbs; i_dim_nurbs++)
        {
          switch (celltype)
          {
            case Core::FE::CellType::nurbs9:
              xi(i_dim_nurbs) =
                  Core::FE::eleNodeNumbering_quad9_nodes_reference[numbering[i_node_nurbs]]
                                                                  [i_dim_nurbs];
              break;
            case Core::FE::CellType::nurbs27:
              xi(i_dim_nurbs) =
                  Core::FE::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_nurbs]]
                                                                  [i_dim_nurbs];
              break;
            default:
              FOUR_C_THROW("The node numbering for the nurbs element shape %s is not implemented",
                  Core::FE::CellTypeToString(ele.Shape()).c_str());
              break;
          }
        }

        // Get the position at the parameter coordinate.
        point_result = EvalNurbsInterpolation<number_of_output_points, dim_nurbs, dim_output>(
            pos_controlpoints, xi, weights, knots, ele.Shape());

        for (unsigned int i_dim = 0; i_dim < dim_output; i_dim++)
          point_coordinates.push_back(point_result(i_dim));
      }
    }

    return number_of_output_points;
  }

  /**
   * \brief Helper function to append the coordinates of vertices of NURBS elements.
   */
  unsigned int AppendVisualizationGeometryNURBSEle(const Core::Elements::Element& ele,
      const Discret::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    using implemented_celltypes =
        Core::FE::CelltypeSequence<Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs27>;
    return Core::FE::CellTypeSwitch<implemented_celltypes>(ele.Shape(),
        [&](auto celltype_t)
        {
          return AppendVisualizationGeometryNURBS<celltype_t()>(
              ele, discret, cell_types, point_coordinates);
        });
  }

  /**
   * \brief Add dof based results to point data vector for elements
   * that use Lagrange shape functions.
   *
   * @param ele (in) Element
   * @param discret (in) discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param result_num_dofs_per_node (in/out) Number of scalar values per point.
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  unsigned int AppendVisualizationDofBasedResultDataVectorLagrangeEle(
      const Core::Elements::Element& ele, const Discret::Discretization& discret,
      const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
      std::vector<double>& vtu_point_result_data)
  {
    const std::vector<int>& numbering = GetVtkCellTypeFromFourCElementShapeType(ele.Shape()).second;

    for (unsigned int inode = 0; inode < (unsigned int)ele.num_node(); ++inode)
    {
      std::vector<int> nodedofs;
      nodedofs.clear();

      // local storage position of desired dof gid in dofset number 0
      discret.Dof((unsigned)0, ele.Nodes()[numbering[inode]], nodedofs);

      // adjust resultdofs according to elements dof
      if (nodedofs.size() < result_num_dofs_per_node)
      {
        result_num_dofs_per_node = nodedofs.size();
      }

      for (unsigned int idof = 0; idof < result_num_dofs_per_node; ++idof)
      {
        const int lid =
            result_data_dofbased->Map().LID(nodedofs[idof + read_result_data_from_dofindex]);

        if (lid > -1)
          vtu_point_result_data.push_back((*result_data_dofbased)[lid]);
        else
          FOUR_C_THROW("received illegal dof local id: %d", lid);
      }
    }

    return ele.num_node();
  }

  /**
   * \brief Add dof based results to point data vector for 3-dimensional elements
   * that use second-order NURBS shape functions.
   *
   * @param ele (in) Element
   * @param discret (in) discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  template <Core::FE::CellType celltype, unsigned int result_num_dofs_per_node>
  unsigned int AppendVisualizationDofBasedResultDataVectorNURBS(const Core::Elements::Element& ele,
      const Discret::Discretization& discret,
      const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      const unsigned int read_result_data_from_dofindex, std::vector<double>& vtu_point_result_data)
  {
    if (read_result_data_from_dofindex != 0)
      FOUR_C_THROW("Nurbs output is only implemented for read_result_data_from_dofindex == 0");

    constexpr int number_of_output_points = Core::FE::num_nodes<celltype>;
    constexpr int dim_nurbs = Core::FE::dim<celltype>;

    // Get the vtk cell information
    const auto vtk_cell_info = GetVTKCellTypeForNURBSElements<celltype>();
    const std::vector<int>& numbering = vtk_cell_info.second;

    // Add the data at the nodes of the nurbs visualization.
    {
      // Get the knots and weights for this element.
      Core::LinAlg::Matrix<number_of_output_points, 1, double> weights(true);
      std::vector<Core::LinAlg::SerialDenseVector> knots(true);
      const bool zero_size =
          Discret::Nurbs::GetMyNurbsKnotsAndWeights(discret, &ele, knots, weights);
      if (zero_size) FOUR_C_THROW("GetMyNurbsKnotsAndWeights has to return a non zero size.");

      // Get the element result vector.
      Core::LinAlg::Matrix<number_of_output_points * result_num_dofs_per_node, 1, double>
          dof_result;
      std::vector<double> eledisp;
      std::vector<int> lm, lmowner, lmstride;
      ele.LocationVector(discret, lm, lmowner, lmstride);
      Core::FE::ExtractMyValues(*result_data_dofbased, eledisp, lm);
      dof_result.SetView(eledisp.data());

      // Loop over the nodes of the nurbs element.
      Core::LinAlg::Matrix<result_num_dofs_per_node, 1, double> point_result;
      Core::LinAlg::Matrix<dim_nurbs, 1, double> xi;
      for (unsigned int i_node_nurbs = 0; i_node_nurbs < number_of_output_points; i_node_nurbs++)
      {
        for (unsigned int i = 0; i < dim_nurbs; i++)
        {
          switch (celltype)
          {
            case Core::FE::CellType::nurbs9:
              xi(i) = Core::FE::eleNodeNumbering_quad9_nodes_reference[numbering[i_node_nurbs]][i];
              break;
            case Core::FE::CellType::nurbs27:
              xi(i) = Core::FE::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_nurbs]][i];
              break;
            default:
              FOUR_C_THROW("The node numbering for the nurbs element shape %s is not implemented",
                  Core::FE::CellTypeToString(ele.Shape()).c_str());
          }
        }

        // Get the field value at the parameter coordinate
        point_result =
            EvalNurbsInterpolation<number_of_output_points, dim_nurbs, result_num_dofs_per_node>(
                dof_result, xi, weights, knots, ele.Shape());

        for (unsigned int i_dim = 0; i_dim < result_num_dofs_per_node; i_dim++)
          vtu_point_result_data.push_back(point_result(i_dim));
      }
    }

    return number_of_output_points;
  }

  /**
   * \brief Add dof based results to point data vector for elements
   * that use NURBS shape functions.
   *
   * @param ele (in) Element
   * @param discret (in) discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param result_num_dofs_per_node (in/out) Number of scalar values per point.
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  unsigned int AppendVisualizationDofBasedResultDataVectorNURBSEle(
      const Core::Elements::Element& ele, const Discret::Discretization& discret,
      const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
      std::vector<double>& vtu_point_result_data)
  {
    using implemented_celltypes =
        Core::FE::CelltypeSequence<Core::FE::CellType::nurbs9, Core::FE::CellType::nurbs27>;
    return Core::FE::CellTypeSwitch<implemented_celltypes>(ele.Shape(),
        [&](auto celltype_t)
        {
          switch (result_num_dofs_per_node)
          {
            case 1:
              return AppendVisualizationDofBasedResultDataVectorNURBS<celltype_t(), 1>(ele, discret,
                  result_data_dofbased, read_result_data_from_dofindex, vtu_point_result_data);
            case 2:
              return AppendVisualizationDofBasedResultDataVectorNURBS<celltype_t(), 2>(ele, discret,
                  result_data_dofbased, read_result_data_from_dofindex, vtu_point_result_data);
            case 3:
              return AppendVisualizationDofBasedResultDataVectorNURBS<celltype_t(), 3>(ele, discret,
                  result_data_dofbased, read_result_data_from_dofindex, vtu_point_result_data);
            default:
              FOUR_C_THROW("The case of result_num_dofs_per_node = %d is not implemented",
                  result_num_dofs_per_node);
          }
        });
  }

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
