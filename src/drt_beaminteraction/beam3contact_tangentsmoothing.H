/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of functions in order to calculate a smooth tangent field

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef BEAM3CONTACT_TANGENTSMOOTHING_H
#define BEAM3CONTACT_TANGENTSMOOTHING_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

namespace CONTACT
{
  // Small class that stores neighbor element information necessary for nodal tangent smoothing
  class B3CNeighbor
  {
   public:
    B3CNeighbor(const DRT::Element* left_neighbor, const DRT::Element* right_neighbor,
        int connecting_node_left, int connecting_node_right);  // standard constructor

    virtual ~B3CNeighbor(){};  // destructor

    const DRT::Element* GetLeftNeighbor() { return left_neighbor_; };
    const DRT::Element* GetRightNeighbor() { return right_neighbor_; };
    int GetLeftConNode() { return connecting_node_left_; };
    int GetRightConNode() { return connecting_node_right_; };


   private:
    const DRT::Element* left_neighbor_;   // pointer to the Element on the left side (eta=-1)
    const DRT::Element* right_neighbor_;  // pointer to the Element on the right side (eta=1)
    int connecting_node_left_;   // local node-ID of the connecting node of the left neighbor
    int connecting_node_right_;  // local node-ID of the connecting node of the right neighbor
  };                             // class B3CNeighbor

  namespace B3TANGENTSMOOTHING
  {
    // determine the neighbour elements of an element
    Teuchos::RCP<B3CNeighbor> DetermineNeigbors(const DRT::Element* element1);

    void GetBoundaryNode(int& nright, int nnode);

    void GetEleLength(Epetra_SerialDenseMatrix& elepos, int& nright, double& length);

    void GetNodalDerivatives(Epetra_SerialDenseMatrix& deriv1, int node, const int nnode,
        double length, const DRT::Element::DiscretizationType distype);

    /*----------------------------------------------------------------------*
     |  Determine nodal tangents                                 meier 05/11|
     *----------------------------------------------------------------------*/
    template <int numnodes>
    LINALG::Matrix<3 * numnodes, 1> CalculateNodalTangents(
        std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        LINALG::Matrix<3 * numnodes, 1> ele1pos_linalg, DRT::Element* ele1,
        Teuchos::RCP<B3CNeighbor> neighbors1)
    {
      LINALG::Matrix<3 * numnodes, 1> node_tangent1(true);  // = vector of element1's node tangents

      // So far, the next lines are only implemented and for linear Reissner beam elements -->
      // dserror() otherwise
      if (numnodes > 2)
        dserror("Tangent smoothing only verified for LINEAR Reissner beam elements!!!");

      // determine boundary node of element1
      const int numnode1 = 2;
      int n_right1 = 0;
      GetBoundaryNode(n_right1, numnode1);

      // vectors for shape functions and their derivatives
      Epetra_SerialDenseMatrix deriv1(
          1, numnode1);  // = matrix for derivatives of shape functions of element1
      Epetra_SerialDenseMatrix ele1pos(3, numnode1);

      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < numnode1; j++)
        {
          ele1pos(i, j) = ele1pos_linalg(3 * j + i);
        }
      }

      double length_ele1 = 0;
      GetEleLength(ele1pos, n_right1, length_ele1);

      const DRT::Element::DiscretizationType distype1 = ele1->Shape();

      //******************************contribution of element
      // 1*****************************************
      for (int node = 1; node < numnode1 + 1; node++)
      {
        // calculate nodal derivatives
        GetNodalDerivatives(deriv1, node, numnode1, length_ele1, distype1);

        for (int k = 0; k < 3; k++)
        {
          for (int j = 0; j < numnode1; j++)
          {
            node_tangent1(k + 3 * (node - 1)) += deriv1(0, j) * ele1pos(k, j);
          }
        }
      }
      //******************************end: contribution of element
      // 1*****************************************

      //****************************** contribution of left neighbor of element
      // 1*****************************************
      if (neighbors1->GetLeftNeighbor() != NULL)
      {
        DRT::Element::DiscretizationType distype_ele1l = neighbors1->GetLeftNeighbor()->Shape();
        int numnode_ele1l = neighbors1->GetLeftNeighbor()->NumNode();
        Epetra_SerialDenseMatrix deriv_neighbors_ele1l(
            1, numnode_ele1l);  // =matrix for derivatives of shape functions for neighbors
        int node_ele1l = neighbors1->GetLeftConNode() + 1;

        Epetra_SerialDenseMatrix temppos(3, numnode_ele1l);

        for (int j = 0; j < numnode_ele1l; j++)
        {
          int tempGID = (neighbors1->GetLeftNeighbor()->NodeIds())[j];
          LINALG::Matrix<3, 1> tempposvector = currentpositions[tempGID];
          for (int i = 0; i < 3; i++) temppos(i, j) = tempposvector(i);
        }

        int n_boundary = 0;
        GetBoundaryNode(n_boundary, numnode_ele1l);
        double length_ele1l = 0;
        GetEleLength(temppos, n_boundary, length_ele1l);

        GetNodalDerivatives(
            deriv_neighbors_ele1l, node_ele1l, numnode_ele1l, length_ele1l, distype_ele1l);

        for (int j = 0; j < numnode_ele1l; j++)
        {
          int tempGID = (neighbors1->GetLeftNeighbor()->NodeIds())[j];
          LINALG::Matrix<3, 1> temppos = currentpositions[tempGID];

          for (int k = 0; k < 3; k++)
          {
            node_tangent1(k) += deriv_neighbors_ele1l(0, j) * temppos(k);
          }
        }

        for (int k = 0; k < 3; k++)
        {
          node_tangent1(k) = 0.5 * node_tangent1(k);
        }
      }
      //******************************end: contribution of left neighbor of element
      // 1*****************************************

      //******************************contribution of right neighbor of element
      // 1*****************************************
      if (neighbors1->GetRightNeighbor() != NULL)
      {
        DRT::Element::DiscretizationType distype_ele1r = neighbors1->GetRightNeighbor()->Shape();
        int numnode_ele1r = neighbors1->GetRightNeighbor()->NumNode();
        Epetra_SerialDenseMatrix deriv_neighbors_ele1r(
            1, numnode_ele1r);  // =matrix for derivatives of shape functions for neighbors
        int node_ele1r = neighbors1->GetRightConNode() + 1;

        Epetra_SerialDenseMatrix temppos(3, numnode_ele1r);

        for (int j = 0; j < numnode_ele1r; j++)
        {
          int tempGID = (neighbors1->GetRightNeighbor()->NodeIds())[j];
          LINALG::Matrix<3, 1> tempposvector = currentpositions[tempGID];
          for (int i = 0; i < 3; i++) temppos(i, j) = tempposvector(i);
        }

        int n_boundary = 0;
        GetBoundaryNode(n_boundary, numnode_ele1r);
        double length_ele1r = 0;
        GetEleLength(temppos, n_boundary, length_ele1r);

        GetNodalDerivatives(
            deriv_neighbors_ele1r, node_ele1r, numnode_ele1r, length_ele1r, distype_ele1r);

        for (int j = 0; j < numnode_ele1r; j++)
        {
          int tempGID = (neighbors1->GetRightNeighbor()->NodeIds())[j];
          LINALG::Matrix<3, 1> temppos = currentpositions[tempGID];

          for (int k = 0; k < 3; k++)
          {
            node_tangent1(k + 3 * n_right1) += deriv_neighbors_ele1r(0, j) * temppos(k);
          }
        }

        for (int k = 0; k < 3; k++)
        {
          node_tangent1(k + 3 * n_right1) = 0.5 * node_tangent1(k + 3 * n_right1);
        }
      }
      //******************************end: contribution of right neighbor of element
      // 1*****************************************

      return node_tangent1;
      //***************** Uncomment this output to check the tangent
      // calculation*************************
      //             std::cout << "nodal tangents Element " << element1_->Id() << " : " <<
      //             ele1tangent_ << std::endl;
    }
    /*----------------------------------------------------------------------*
     |  end: Determine nodal tangents
     *----------------------------------------------------------------------*/

    template <int numnodes, int numnodalvalues>
    void ComputeTangentsAndDerivs(LINALG::Matrix<3, 1, TYPE>& t, LINALG::Matrix<3, 1, TYPE>& t_xi,
        const LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth1,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        const LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xi)
    {
      t.Clear();
      t_xi.Clear();

      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3 * numnodes; j++)
        {
          t(i) += N(i, j) * nodaltangentssmooth1(j);
          t_xi(i) += N_xi(i, j) * nodaltangentssmooth1(j);
        }
      }

      return;
    };

  }  // namespace B3TANGENTSMOOTHING

}  // namespace CONTACT

#endif  // #ifndef BEAM3CONTACT_TANGENTSMOOTHING_H
