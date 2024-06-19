/*----------------------------------------------------------------------------*/
/*! \file

\brief A set of functions in order to calculate a smooth tangent field

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_TANGENTSMOOTHING_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_BEAM_CONTACT_TANGENTSMOOTHING_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace BEAMINTERACTION
{
  //! \brief Small class that stores neighbor element information necessary for nodal tangent
  //! smoothing
  class B3CNeighbor
  {
   public:
    B3CNeighbor(const Core::Elements::Element* left_neighbor,
        const Core::Elements::Element* right_neighbor, int connecting_node_left,
        int connecting_node_right);

    const Core::Elements::Element* GetLeftNeighbor() { return left_neighbor_; };
    const Core::Elements::Element* GetRightNeighbor() { return right_neighbor_; };
    int GetLeftConNode() { return connecting_node_left_; };
    int GetRightConNode() { return connecting_node_right_; };

   private:
    const Core::Elements::Element*
        left_neighbor_;  // pointer to the Element on the left side (eta=-1)
    const Core::Elements::Element*
        right_neighbor_;         // pointer to the Element on the right side (eta=1)
    int connecting_node_left_;   // local node-ID of the connecting node of the left neighbor
    int connecting_node_right_;  // local node-ID of the connecting node of the right neighbor
  };                             // class B3CNeighbor

  namespace Beam3TangentSmoothing
  {
    //! \brief Determine the neighbour elements of an element
    Teuchos::RCP<B3CNeighbor> DetermineNeigbors(const Core::Elements::Element* element1);

    //! \brief Get boundary node
    int GetBoundaryNode(const int nnode);

    //! \brief Get the element length based on element positions
    double GetEleLength(const Core::LinAlg::SerialDenseMatrix& elepos, const int nright);

    //! \brief Get the nodal derivative matrix
    Core::LinAlg::SerialDenseMatrix GetNodalDerivatives(
        const int node, const int nnode, double length, const Core::FE::CellType distype);

    //! \brief Get the nodal tangent matrix
    template <int numnodes>
    Core::LinAlg::Matrix<3 * numnodes, 1> CalculateNodalTangents(
        std::map<int, Core::LinAlg::Matrix<3, 1>>& currentpositions,
        Core::LinAlg::Matrix<3 * numnodes, 1> ele1pos_linalg, Core::Elements::Element* ele1,
        Teuchos::RCP<B3CNeighbor> neighbors1)
    {
      Core::LinAlg::Matrix<3 * numnodes, 1> node_tangent1(
          true);  // = vector of element1's node tangents

      // So far, the next lines are only implemented and for linear Reissner beam elements -->
      // FOUR_C_THROW() otherwise
      if (numnodes > 2)
        FOUR_C_THROW("Tangent smoothing only verified for LINEAR Reissner beam elements!!!");

      // determine boundary node of element1
      const int numnode1 = 2;
      const int n_right1 = GetBoundaryNode(numnode1);

      // vectors for shape functions and their derivatives
      Core::LinAlg::SerialDenseMatrix ele1pos(3, numnode1);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < numnode1; j++) ele1pos(i, j) = ele1pos_linalg(3 * j + i);

      const double length_ele1 = GetEleLength(ele1pos, n_right1);
      const Core::FE::CellType distype1 = ele1->Shape();

      // contribution of element
      for (int node = 1; node < numnode1 + 1; node++)
      {
        // calculate nodal derivatives
        Core::LinAlg::SerialDenseMatrix deriv1 =
            GetNodalDerivatives(node, numnode1, length_ele1, distype1);

        for (int k = 0; k < 3; k++)
          for (int j = 0; j < numnode1; j++)
            node_tangent1(k + 3 * (node - 1)) += deriv1(0, j) * ele1pos(k, j);
      }

      // contribution of left neighbor of element
      if (neighbors1->GetLeftNeighbor() != nullptr)
      {
        const Core::FE::CellType distype_ele1l = neighbors1->GetLeftNeighbor()->Shape();
        const int numnode_ele1l = neighbors1->GetLeftNeighbor()->num_node();
        const int node_ele1l = neighbors1->GetLeftConNode() + 1;

        Core::LinAlg::SerialDenseMatrix temppos(3, numnode_ele1l);

        for (int j = 0; j < numnode_ele1l; j++)
        {
          const int tempGID = (neighbors1->GetLeftNeighbor()->NodeIds())[j];
          Core::LinAlg::Matrix<3, 1> tempposvector = currentpositions[tempGID];
          for (int i = 0; i < 3; i++) temppos(i, j) = tempposvector(i);
        }

        const int n_boundary = GetBoundaryNode(numnode_ele1l);
        const double length_ele1l = GetEleLength(temppos, n_boundary);

        Core::LinAlg::SerialDenseMatrix deriv_neighbors_ele1l =
            GetNodalDerivatives(node_ele1l, numnode_ele1l, length_ele1l, distype_ele1l);

        for (int j = 0; j < numnode_ele1l; j++)
        {
          const int tempGID = (neighbors1->GetLeftNeighbor()->NodeIds())[j];
          Core::LinAlg::Matrix<3, 1> tempposvector = currentpositions[tempGID];

          for (int k = 0; k < 3; k++)
            node_tangent1(k) += deriv_neighbors_ele1l(0, j) * tempposvector(k);
        }

        for (int k = 0; k < 3; k++) node_tangent1(k) = 0.5 * node_tangent1(k);
      }

      // contribution of right neighbor of element
      if (neighbors1->GetRightNeighbor() != nullptr)
      {
        const Core::FE::CellType distype_ele1r = neighbors1->GetRightNeighbor()->Shape();
        const int numnode_ele1r = neighbors1->GetRightNeighbor()->num_node();
        const int node_ele1r = neighbors1->GetRightConNode() + 1;

        Core::LinAlg::SerialDenseMatrix temppos(3, numnode_ele1r);

        for (int j = 0; j < numnode_ele1r; j++)
        {
          const int tempGID = (neighbors1->GetRightNeighbor()->NodeIds())[j];
          Core::LinAlg::Matrix<3, 1> tempposvector = currentpositions[tempGID];
          for (int i = 0; i < 3; i++) temppos(i, j) = tempposvector(i);
        }

        const int n_boundary = GetBoundaryNode(numnode_ele1r);
        const double length_ele1r = GetEleLength(temppos, n_boundary);

        Core::LinAlg::SerialDenseMatrix deriv_neighbors_ele1r =
            GetNodalDerivatives(node_ele1r, numnode_ele1r, length_ele1r, distype_ele1r);

        for (int j = 0; j < numnode_ele1r; j++)
        {
          const int tempGID = (neighbors1->GetRightNeighbor()->NodeIds())[j];
          Core::LinAlg::Matrix<3, 1> tempposvector = currentpositions[tempGID];

          for (int k = 0; k < 3; k++)
            node_tangent1(k + 3 * n_right1) += deriv_neighbors_ele1r(0, j) * tempposvector(k);
        }

        for (int k = 0; k < 3; k++)
          node_tangent1(k + 3 * n_right1) = 0.5 * node_tangent1(k + 3 * n_right1);
      }

      return node_tangent1;
    }

    //! \brief Get the tangent and derivative matrix
    template <int numnodes, int numnodalvalues>
    void ComputeTangentsAndDerivs(Core::LinAlg::Matrix<3, 1, TYPE>& t,
        Core::LinAlg::Matrix<3, 1, TYPE>& t_xi,
        const Core::LinAlg::Matrix<3 * numnodes, 1> nodaltangentssmooth1,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        const Core::LinAlg::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xi)
    {
      t.clear();
      t_xi.clear();

      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3 * numnodes; j++)
        {
          t(i) += N(i, j) * nodaltangentssmooth1(j);
          t_xi(i) += N_xi(i, j) * nodaltangentssmooth1(j);
        }
      }
    }

  }  // namespace Beam3TangentSmoothing
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
