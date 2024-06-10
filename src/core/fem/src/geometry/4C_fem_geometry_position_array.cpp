/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations


\level 2

*----------------------------------------------------------------------*/


#include "4C_fem_geometry_position_array.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN

Core::LinAlg::SerialDenseMatrix Core::Geo::InitialPositionArray(
    const Core::Elements::Element* const ele)
{
  const int numnode = ele->num_node();
  Core::LinAlg::SerialDenseMatrix xyze(3, numnode);
  const Core::Nodes::Node* const* nodes = ele->Nodes();
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


void Core::Geo::InitialPositionArray(
    Core::LinAlg::SerialDenseMatrix& xyze, const Core::Elements::Element* const ele)
{
  const int numnode = ele->num_node();
  xyze.shape(3, numnode);
  const Core::Nodes::Node* const* nodes = ele->Nodes();
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


Core::LinAlg::SerialDenseMatrix Core::Geo::getCurrentNodalPositions(
    const Core::Elements::Element* const ele,  ///< element with nodal pointers
    const std::map<int, Core::LinAlg::Matrix<3, 1>>&
        currentcutterpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->num_node();
  Core::LinAlg::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const Core::LinAlg::Matrix<3, 1>& x = currentcutterpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}


Core::LinAlg::SerialDenseMatrix Core::Geo::getCurrentNodalPositions(
    const Teuchos::RCP<const Core::Elements::Element> ele,  ///< pointer on element
    const std::map<int, Core::LinAlg::Matrix<3, 1>>&
        currentpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->num_node();
  Core::LinAlg::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const Core::LinAlg::Matrix<3, 1>& x = currentpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}

FOUR_C_NAMESPACE_CLOSE
