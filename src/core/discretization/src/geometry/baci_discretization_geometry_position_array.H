/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations

\level 2

*----------------------------------------------------------------------*/
#ifndef BACI_DISCRETIZATION_GEOMETRY_POSITION_ARRAY_H
#define BACI_DISCRETIZATION_GEOMETRY_POSITION_ARRAY_H


#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_lib_node.H"
#include "baci_linalg_serialdensematrix.H"


namespace CORE::GEO
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
  void fillInitialPositionArray(const ::DRT::Element* const ele, M& xyze)
  {
    const int numnode = ele->NumNode();

    const ::DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != nullptr,
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
  template <CORE::FE::CellType distype, class M>
  void fillInitialPositionArray(const ::DRT::Element* const ele, M& xyze)
  {
    dsassert(distype == ele->Shape(), "mismatch in distype");
    const int numnode = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    const ::DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != nullptr,
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
  template <CORE::FE::CellType distype, int dim, class M>
  void fillInitialPositionArray(const ::DRT::Element* const ele, M& xyze)
  {
    dsassert(distype == ele->Shape(), "mismatch in distype");
    const int numnode = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    const ::DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != nullptr,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    dsassert((dim > 0) && (dim < 4), "Illegal number of space dimensions");

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
  void InitialPositionArray(CORE::LINALG::SerialDenseMatrix& xyze, const ::DRT::Element* const ele);


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   */
  CORE::LINALG::SerialDenseMatrix InitialPositionArray(const ::DRT::Element* const
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
  CORE::LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const ::DRT::Element* const ele,  ///< element with nodal pointers
      const std::map<int, CORE::LINALG::Matrix<3, 1>>&
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
  CORE::LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const Teuchos::RCP<const ::DRT::Element> ele,  ///< pointer on element
      const std::map<int, CORE::LINALG::Matrix<3, 1>>&
          currentpositions  ///< current positions of all cutter nodes
  );

}  // namespace CORE::GEO

#endif
