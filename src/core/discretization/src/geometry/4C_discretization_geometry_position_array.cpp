/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations


\level 2

*----------------------------------------------------------------------*/


#include "4C_discretization_geometry_position_array.hpp"

#include "4C_discretization_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

CORE::LINALG::SerialDenseMatrix CORE::GEO::InitialPositionArray(
    const CORE::Elements::Element* const ele)
{
  const int numnode = ele->num_node();
  CORE::LINALG::SerialDenseMatrix xyze(3, numnode);
  const CORE::Nodes::Node* const* nodes = ele->Nodes();
  if (nodes == nullptr)
  {
    FOUR_C_THROW("element has no nodal pointers, so getting a position array doesn't make sense!");
  }
  for (int inode = 0; inode < numnode; inode++)
  {
    const double* x = nodes[inode]->X().data();
    std::copy(x, x + 3, &xyze(0, inode));
  }
  return xyze;
}


void CORE::GEO::InitialPositionArray(
    CORE::LINALG::SerialDenseMatrix& xyze, const CORE::Elements::Element* const ele)
{
  const int numnode = ele->num_node();
  xyze.shape(3, numnode);
  const CORE::Nodes::Node* const* nodes = ele->Nodes();
  if (nodes == nullptr)
  {
    FOUR_C_THROW("element has no nodal pointers, so getting a position array doesn't make sense!");
  }
  for (int inode = 0; inode < numnode; inode++)
  {
    const double* x = nodes[inode]->X().data();
    std::copy(x, x + 3, &xyze(0, inode));
  }
}


CORE::LINALG::SerialDenseMatrix CORE::GEO::getCurrentNodalPositions(
    const CORE::Elements::Element* const ele,  ///< element with nodal pointers
    const std::map<int, CORE::LINALG::Matrix<3, 1>>&
        currentcutterpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->num_node();
  CORE::LINALG::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const CORE::LINALG::Matrix<3, 1>& x = currentcutterpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}


CORE::LINALG::SerialDenseMatrix CORE::GEO::getCurrentNodalPositions(
    const Teuchos::RCP<const CORE::Elements::Element> ele,  ///< pointer on element
    const std::map<int, CORE::LINALG::Matrix<3, 1>>&
        currentpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->num_node();
  CORE::LINALG::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const CORE::LINALG::Matrix<3, 1>& x = currentpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}

FOUR_C_NAMESPACE_CLOSE
