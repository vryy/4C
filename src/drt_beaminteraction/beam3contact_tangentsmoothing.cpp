/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of functions in order to calculate a smooth tangent field

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3contact_tangentsmoothing.H"

#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

/*----------------------------------------------------------------------*
 |  Constructor of class B3CNeighbor                         meier 04/11|
 *----------------------------------------------------------------------*/
CONTACT::B3CNeighbor::B3CNeighbor(const DRT::Element* left_neighbor,
    const DRT::Element* right_neighbor, int connecting_node_left, int connecting_node_right)
    : left_neighbor_(left_neighbor),
      right_neighbor_(right_neighbor),
      connecting_node_left_(connecting_node_left),
      connecting_node_right_(connecting_node_right)
{
  return;
}
/*----------------------------------------------------------------------*
 |  end: Constructor of class B3CNeighbor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Determine Neighbor Elements                              meier 04/11|
 *----------------------------------------------------------------------*/
Teuchos::RCP<CONTACT::B3CNeighbor> CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(
    const DRT::Element* element1)
{
  const DRT::Element* left_neighbor = NULL;
  const DRT::Element* right_neighbor = NULL;
  int connecting_node_left = 0;
  int connecting_node_right = 0;

  // number of nodes element1
  int numnode = element1->NumNode();

  // n_right is the local node-ID of the elements right node (at xi = 1) whereas the elements left
  // node (at xi = -1) allways has the local ID 1 For documentation of the node numbering see also
  // the file beam3.H
  int n_right = 0;
  if (numnode == 2)
  {
    n_right = 1;
  }
  else
  {
    n_right = numnode - 2;
  }


  int globalneighborId = 0;
  int globalnodeId = 0;

  //*******************local node 1 of element1 --> xi(element1)=-1 --> left
  // neighbor******************
  globalnodeId = *element1->NodeIds();

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->Nodes())).NumElement() > 2)
  {
    dserror(
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
  //*************************************************************************************************

  //*******************************local node n_right of element1 --> xi(element1)=1 --> right
  // neighbor*******************
  globalnodeId = *(element1->NodeIds() + n_right);

  // only one neighbor element on each side of the considered element is allowed
  if ((**(element1->Nodes() + n_right)).NumElement() > 2)
  {
    dserror(
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
  //********************************************************************************************************

  return Teuchos::rcp(new CONTACT::B3CNeighbor(
      left_neighbor, right_neighbor, connecting_node_left, connecting_node_right));
}
/*----------------------------------------------------------------------*
 |  end: Determine Neighbor Elements
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  compute right boundary node of an element                meier 05/11|
 *----------------------------------------------------------------------*/
void CONTACT::B3TANGENTSMOOTHING::GetBoundaryNode(int& nright, int nnode)
{
  if (nnode == 2)
  {
    nright = 1;
  }
  else
  {
    nright = nnode - 2;
  }
}
/*----------------------------------------------------------------------*
 |  end: compute right boundary node of an element           meier 05/11|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  compute element length                                   meier 05/11|
 *----------------------------------------------------------------------*/
void CONTACT::B3TANGENTSMOOTHING::GetEleLength(
    Epetra_SerialDenseMatrix& elepos, int& nright, double& length)
{
  length = 0;
  for (int i = 0; i < 3; i++)
  {
    length += (elepos(i, nright) - elepos(i, 0)) * (elepos(i, nright) - elepos(i, 0));
  }
  length = sqrt(length);
}

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 05/11|
 *----------------------------------------------------------------------*/
void CONTACT::B3TANGENTSMOOTHING::GetNodalDerivatives(Epetra_SerialDenseMatrix& deriv1, int node,
    const int nnode, double length, const DRT::Element::DiscretizationType distype)
{
  if (node == nnode)
  {
    DRT::UTILS::shape_function_1D_deriv1(deriv1, -1.0 + 2.0 / (nnode - 1), distype);
  }
  else
  {
    if (node == 1)
    {
      DRT::UTILS::shape_function_1D_deriv1(deriv1, -1.0, distype);
    }
    else
    {
      DRT::UTILS::shape_function_1D_deriv1(deriv1, -1.0 + node * 2.0 / (nnode - 1), distype);
    }
  }


  for (int i = 0; i < nnode; i++)
  {
    deriv1(0, i) = 2 * deriv1(0, i) / length;
  }



  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/
