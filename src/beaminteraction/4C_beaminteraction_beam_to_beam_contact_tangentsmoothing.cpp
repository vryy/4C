/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of functions in order to calculate a smooth tangent field

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"

#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_node.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::B3CNeighbor::B3CNeighbor(const DRT::Element* left_neighbor,
    const DRT::Element* right_neighbor, int connecting_node_left, int connecting_node_right)
    : left_neighbor_(left_neighbor),
      right_neighbor_(right_neighbor),
      connecting_node_left_(connecting_node_left),
      connecting_node_right_(connecting_node_right)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> BEAMINTERACTION::B3TANGENTSMOOTHING::DetermineNeigbors(
    const DRT::Element* element1)
{
  const DRT::Element* left_neighbor = nullptr;
  const DRT::Element* right_neighbor = nullptr;
  int connecting_node_left = 0;
  int connecting_node_right = 0;

  // number of nodes element1
  int numnode = element1->num_node();

  // n_right is the local node-ID of the elements right node (at xi = 1) whereas the elements left
  // node (at xi = -1) allways has the local ID 1
  int n_right;
  if (numnode == 2)
    n_right = 1;
  else
    n_right = numnode - 2;

  int globalneighborId = 0;
  int globalnodeId = 0;

  // local node 1 of element1 --> xi(element1)=-1 --> left neighbor
  globalnodeId = *element1->NodeIds();

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->Nodes())).NumElement() > 2)
  {
    FOUR_C_THROW(
        "ERROR: The implemented smoothing routine does not work for more than 2 adjacent Elements "
        "per node");
  }

  // loop over all elements adjacent to the left node of element 1
  for (int y = 0; y < (**(element1->Nodes())).NumElement(); ++y)
  {
    globalneighborId = (*((**(element1->Nodes())).Elements() + y))->Id();

    if (globalneighborId != element1->Id())
    {
      left_neighbor = (*((**(element1->Nodes())).Elements() + y));

      // if node 1 of the left neighbor is the connecting node:
      if (*((*((**(element1->Nodes())).Elements() + y))->NodeIds()) == globalnodeId)
      {
        connecting_node_left = 0;
      }
      // otherwise node n_right of the left neighbor is the connecting node
      else
      {
        connecting_node_left = n_right;
      }
    }
  }

  // ocal node n_right of element1 --> xi(element1)=1 --> right neighbor
  globalnodeId = *(element1->NodeIds() + n_right);

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->Nodes() + n_right)).NumElement() > 2)
  {
    FOUR_C_THROW(
        "ERROR: The implemented smoothing routine does not work for more than 2 adjacent Elements "
        "per node");
  }

  // loop over all elements adjacent to the right node of element 1
  for (int y = 0; y < (**(element1->Nodes() + n_right)).NumElement(); ++y)
  {
    globalneighborId = (*((**(element1->Nodes() + n_right)).Elements() + y))->Id();

    if (globalneighborId != element1->Id())
    {
      right_neighbor = (*((**(element1->Nodes() + n_right)).Elements() + y));

      // if node 1 of the right neighbor is the connecting node:
      if (*((*((**(element1->Nodes() + n_right)).Elements() + y))->NodeIds()) == globalnodeId)
      {
        connecting_node_right = 0;
      }

      // otherwise node n_right of the right neighbor is the connecting node
      else
      {
        connecting_node_right = n_right;
      }
    }
  }

  return Teuchos::rcp(new BEAMINTERACTION::B3CNeighbor(
      left_neighbor, right_neighbor, connecting_node_left, connecting_node_right));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int BEAMINTERACTION::B3TANGENTSMOOTHING::GetBoundaryNode(const int nnode)
{
  if (nnode == 2)
    return 1;
  else
    return nnode - 2;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double BEAMINTERACTION::B3TANGENTSMOOTHING::GetEleLength(
    const CORE::LINALG::SerialDenseMatrix& elepos, const int nright)
{
  double length = 0.0;
  for (int i = 0; i < 3; i++)
    length += (elepos(i, nright) - elepos(i, 0)) * (elepos(i, nright) - elepos(i, 0));

  return sqrt(length);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix BEAMINTERACTION::B3TANGENTSMOOTHING::GetNodalDerivatives(
    const int node, const int nnode, const double length, const CORE::FE::CellType distype)
{
  CORE::LINALG::SerialDenseMatrix deriv1(1, nnode);

  if (node == nnode)
    CORE::FE::shape_function_1D_deriv1(deriv1, -1.0 + 2.0 / (nnode - 1), distype);
  else
  {
    if (node == 1)
      CORE::FE::shape_function_1D_deriv1(deriv1, -1.0, distype);
    else
      CORE::FE::shape_function_1D_deriv1(deriv1, -1.0 + node * 2.0 / (nnode - 1), distype);
  }

  for (int i = 0; i < nnode; i++) deriv1(0, i) = 2.0 * deriv1(0, i) / length;

  return deriv1;
}

FOUR_C_NAMESPACE_CLOSE
