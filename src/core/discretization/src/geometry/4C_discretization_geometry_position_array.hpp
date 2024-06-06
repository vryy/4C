/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations

\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_POSITION_ARRAY_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_POSITION_ARRAY_HPP


#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_discretization_fem_general_node.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note xyze is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  template <class M>
  void fillInitialPositionArray(const Core::Elements::Element* const ele, M& xyze)
  {
    const int numnode = ele->num_node();

    const Core::Nodes::Node* const* nodes = ele->Nodes();
    FOUR_C_ASSERT(nodes != nullptr,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    for (int inode = 0; inode < numnode; inode++)
    {
      const auto& x = nodes[inode]->X();
      xyze(0, inode) = x[0];
      xyze(1, inode) = x[1];
      xyze(2, inode) = x[2];
    }
  }


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  template <Core::FE::CellType distype, class M>
  void fillInitialPositionArray(const Core::Elements::Element* const ele, M& xyze)
  {
    FOUR_C_ASSERT(distype == ele->Shape(), "mismatch in distype");
    const int numnode = Core::FE::num_nodes<distype>;

    const Core::Nodes::Node* const* nodes = ele->Nodes();
    FOUR_C_ASSERT(nodes != nullptr,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    for (int inode = 0; inode < numnode; inode++)
    {
      const auto& x = nodes[inode]->X();
      xyze(0, inode) = x[0];
      xyze(1, inode) = x[1];
      xyze(2, inode) = x[2];
    }
  }


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note array is defined as (dim,numnode) with user-specified number of space dimensions
   *       that are of interest for the element application calling this method
   *
   * \return Array with 1, 2 or 3 dimensional position of all element nodes in the coordinate system
   * of the nodes
   */
  template <Core::FE::CellType distype, int dim, class M>
  void fillInitialPositionArray(const Core::Elements::Element* const ele, M& xyze)
  {
    FOUR_C_ASSERT(distype == ele->Shape(), "mismatch in distype");
    const int numnode = Core::FE::num_nodes<distype>;

    const Core::Nodes::Node* const* nodes = ele->Nodes();
    FOUR_C_ASSERT(nodes != nullptr,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    FOUR_C_ASSERT((dim > 0) && (dim < 4), "Illegal number of space dimensions");

    for (int inode = 0; inode < numnode; inode++)
    {
      const double* x = nodes[inode]->X().data();
      // copy the values in the current column
      std::copy(x, x + dim, &xyze(0, inode));
      // fill the remaining entries of the column with zeros, if the given matrix has
      // the wrong row dimension (this is primarily for safety reasons)
      std::fill(&xyze(0, inode) + dim, &xyze(0, inode) + xyze.numRows(), 0.0);
    }
  }


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  void InitialPositionArray(
      Core::LinAlg::SerialDenseMatrix& xyze, const Core::Elements::Element* const ele);


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  Core::LinAlg::SerialDenseMatrix InitialPositionArray(const Core::Elements::Element* const
          ele  ///< pointer to element, whose nodes we evaluate for their position
  );


  /*!
   * get current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  Core::LinAlg::SerialDenseMatrix getCurrentNodalPositions(
      const Core::Elements::Element* const ele,  ///< element with nodal pointers
      const std::map<int, Core::LinAlg::Matrix<3, 1>>&
          currentcutterpositions  ///< current positions of all cutter nodes
  );


  /*!
   * get current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  Core::LinAlg::SerialDenseMatrix getCurrentNodalPositions(
      const Teuchos::RCP<const Core::Elements::Element> ele,  ///< pointer on element
      const std::map<int, Core::LinAlg::Matrix<3, 1>>&
          currentpositions  ///< current positions of all cutter nodes
  );

}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
