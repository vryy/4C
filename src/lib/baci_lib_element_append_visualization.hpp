/*---------------------------------------------------------------------*/
/*! \file

\brief Functions to append visualization for runtime output

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LIB_ELEMENT_APPEND_VISUALIZATION_HPP
#define BACI_LIB_ELEMENT_APPEND_VISUALIZATION_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "baci_lib_element_vtk_cell_type_register.hpp"
#include "baci_lib_utils.hpp"
#include "baci_nurbs_discret_nurbs_utils.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /**
   * \brief Add the element geometry visualization for elements that use Lagrange shape functions
   *
   * @param ele (in) Element
   * @param discret (in) Discretization
   * @param cell_types (in/out) cell type data vector
   * @param point_coordinates (in/out) point coordinates for the representation of this element
   * @return Number of added points
   */
  unsigned int AppendVisualizationGeometryLagrangeEle(const DRT::Element& ele,
      const DRT::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    const unsigned int num_spatial_dimensions = 3;
    auto vtk_cell_info = DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele.Shape());
    const std::vector<int>& numbering = vtk_cell_info.second;

    // Add the cell type to the output.
    cell_types.push_back(vtk_cell_info.first);

    // Add each node to the output.
    const DRT::Node* const* nodes = ele.Nodes();

    for (int inode = 0; inode < ele.NumNode(); ++inode)
      for (unsigned int idim = 0; idim < num_spatial_dimensions; ++idim)
        point_coordinates.push_back(nodes[numbering[inode]]->X()[idim]);

    // Return the number of added points.
    return ele.NumNode();
  }

  /**
   * \brief Helper function to evaluate the NURBS interpolation inside the element.
   */
  template <unsigned int n_points, unsigned int n_dof>
  CORE::LINALG::Matrix<n_dof, 1, double> EvalNurbs3DInterpolation(
      const CORE::LINALG::Matrix<n_points * n_dof, 1, double>& controlpoint_data,
      const CORE::LINALG::Matrix<n_dof, 1, double>& xi,
      const CORE::LINALG::Matrix<n_points, 1, double>& weights,
      const std::vector<CORE::LINALG::SerialDenseVector>& knots, const CORE::FE::CellType& distype)
  {
    CORE::LINALG::Matrix<n_dof, 1, double> point_result;

    // Get the shape functions.
    CORE::LINALG::Matrix<n_points, 1, double> N;
    CORE::FE::NURBS::nurbs_get_3D_funct(N, xi, knots, weights, distype);

    for (unsigned int i_node_nurbs = 0; i_node_nurbs < n_points; i_node_nurbs++)
    {
      for (unsigned int i_dim = 0; i_dim < n_dof; i_dim++)
        point_result(i_dim) += N(i_node_nurbs) * controlpoint_data(i_node_nurbs * 3 + i_dim);
    }

    return point_result;
  }

  /**
   * \brief Add the element geometry visualization for elements that use NURBS shape functions
   *
   * @param ele (in) Element
   * @param discret (in) Discretization
   * @param cell_types (in/out) cell type data vector
   * @param point_coordinates (in/out) point coordinates for the representation of this element
   * @return Number of added points
   */
  unsigned int AppendVisualizationGeometryNURBS27(const DRT::Element& ele,
      const DRT::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    const int number_of_output_points = 27;

    // For vtk output of nurbs27, we use the vtk cell of a hex27 element
    const auto vtk_cell_info =
        DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(CORE::FE::CellType::hex27);
    const std::vector<int>& numbering = vtk_cell_info.second;

    // Add the cell type to the output.
    cell_types.push_back(vtk_cell_info.first);

    // Create the vertices for the visualization of a nurbs27 element.
    {
      // Get the knots and weights of this element.
      CORE::LINALG::Matrix<27, 1, double> weights(true);
      std::vector<CORE::LINALG::SerialDenseVector> knots(true);
      const bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discret, &ele, knots, weights);
      if (zero_size) dserror("GetMyNurbsKnotsAndWeights has to return a non zero size.");

      // Get the position of the control points in the reference configuration.
      CORE::LINALG::Matrix<27 * 3, 1, double> pos_controlpoints;
      for (unsigned int i_controlpoint = 0; i_controlpoint < (unsigned int)ele.NumNode();
           ++i_controlpoint)
      {
        const DRT::Node* controlpoint = ele.Nodes()[i_controlpoint];
        for (int i_dim = 0; i_dim < 3; ++i_dim)
          pos_controlpoints(3 * i_controlpoint + i_dim) = controlpoint->X()[i_dim];
      }

      // Loop over the control points of the nurbs27 element.
      CORE::LINALG::Matrix<3, 1, double> point_result;
      CORE::LINALG::Matrix<3, 1, double> xi;
      for (unsigned int i_node_nurbs = 0; i_node_nurbs < number_of_output_points; i_node_nurbs++)
      {
        for (unsigned int i = 0; i < 3; i++)
          xi(i) = CORE::FE::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_nurbs]][i];

        // Get the reference position at the parameter coordinate.
        point_result = EvalNurbs3DInterpolation(pos_controlpoints, xi, weights, knots, ele.Shape());

        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          point_coordinates.push_back(point_result(i_dim));
      }
    }

    return number_of_output_points;
  }

  /**
   * \brief Helper function to append the coordinates of vertices of NURBS elements.
   */
  unsigned int AppendVisualizationGeometryNURBSEle(const DRT::Element& ele,
      const DRT::Discretization& discret, std::vector<uint8_t>& cell_types,
      std::vector<double>& point_coordinates)
  {
    switch (ele.Shape())
    {
      case CORE::FE::CellType::nurbs27:
        return AppendVisualizationGeometryNURBS27(ele, discret, cell_types, point_coordinates);
      default:
        dserror("The visualization for the nurbs element shape %s is not implemented",
            CORE::FE::CellTypeToString(ele.Shape()).c_str());
    }
  }

  /**
   * \brief Add dof based results to point data vector for elements
   * that use Lagrange shape functions.
   *
   * @param ele (in) Element
   * @param discret (in) Discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param result_num_dofs_per_node (in/out) Number of scalar values per point.
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  unsigned int AppendVisualizationDofBasedResultDataVectorLagrangeEle(const DRT::Element& ele,
      const DRT::Discretization& discret, const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
      std::vector<double>& vtu_point_result_data)
  {
    const std::vector<int>& numbering =
        DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele.Shape()).second;

    for (unsigned int inode = 0; inode < (unsigned int)ele.NumNode(); ++inode)
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
          dserror("received illegal dof local id: %d", lid);
      }
    }

    return ele.NumNode();
  }

  /**
   * \brief Add dof based results to point data vector for 3-dimensional elements
   * that use second-order NURBS shape functions.
   *
   * @param ele (in) Element
   * @param discret (in) Discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param result_num_dofs_per_node (in/out) Number of scalar values per point.
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  unsigned int AppendVisualizationDofBasedResultDataVectorNURBS27(const DRT::Element& ele,
      const DRT::Discretization& discret, const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
      std::vector<double>& vtu_point_result_data)
  {
    if (read_result_data_from_dofindex != 0)
      dserror("Nurbs output is only implemented for read_result_data_from_dofindex == 0");

    if (result_num_dofs_per_node != 3)
      dserror(
          "The nurbs elements can only output nodal data with dimension 3, e.g., displacements");

    const int number_of_output_points = 27;

    // For vtk output of nurbs27, we use the vtk cell of a hex27 element
    const std::vector<int>& numbering =
        DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(CORE::FE::CellType::hex27).second;

    // Add the data at the nodes of the nurbs27 visualization.
    {
      // Get the knots and weights for this element.
      CORE::LINALG::Matrix<27, 1, double> weights(true);
      std::vector<CORE::LINALG::SerialDenseVector> knots(true);
      const bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discret, &ele, knots, weights);
      if (zero_size) dserror("GetMyNurbsKnotsAndWeights has to return a non zero size.");

      // Get the element result vector.
      CORE::LINALG::Matrix<27 * 3, 1, double> dof_result;
      std::vector<double> eledisp;
      std::vector<int> lm, lmowner, lmstride;
      ele.LocationVector(discret, lm, lmowner, lmstride);
      DRT::UTILS::ExtractMyValues(*result_data_dofbased, eledisp, lm);
      dof_result.SetView(eledisp.data());

      // Loop over the nodes of the nurbs27 element.
      CORE::LINALG::Matrix<3, 1, double> point_result;
      CORE::LINALG::Matrix<3, 1, double> xi;
      for (unsigned int i_node_nurbs = 0; i_node_nurbs < number_of_output_points; i_node_nurbs++)
      {
        for (unsigned int i = 0; i < 3; i++)
          xi(i) = CORE::FE::eleNodeNumbering_hex27_nodes_reference[numbering[i_node_nurbs]][i];

        // Get the reference position at the parameter coordinate.
        point_result = EvalNurbs3DInterpolation(dof_result, xi, weights, knots, ele.Shape());

        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
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
   * @param discret (in) Discretization
   * @param result_data_dofbased (in) Global vector with results
   * @param result_num_dofs_per_node (in/out) Number of scalar values per point.
   * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
   * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
   * @param vtu_point_result_data (in/out) Result data vector.
   * @return Number of points added by this element.
   */
  unsigned int AppendVisualizationDofBasedResultDataVectorNURBSEle(const DRT::Element& ele,
      const DRT::Discretization& discret, const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int& result_num_dofs_per_node, const unsigned int read_result_data_from_dofindex,
      std::vector<double>& vtu_point_result_data)
  {
    switch (ele.Shape())
    {
      case CORE::FE::CellType::nurbs27:
        return AppendVisualizationDofBasedResultDataVectorNURBS27(ele, discret,
            result_data_dofbased, result_num_dofs_per_node, read_result_data_from_dofindex,
            vtu_point_result_data);
      default:
        dserror("The visualization for the nurbs element shape %s is not implemented",
            CORE::FE::CellTypeToString(ele.Shape()).c_str());
    }
  }

}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif
