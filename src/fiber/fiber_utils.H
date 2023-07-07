/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration of utility functions for fiber interpolation

\level 3
*/
/*----------------------------------------------------------------------*/

#ifndef FIBER_UTILS_H
#define FIBER_UTILS_H

#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "inpar_structure.H"

namespace DRT
{
  namespace FIBER
  {
    class NodalFiberHolder;

    namespace UTILS
    {
      /*!
       * \brief Projects the fibers and angles to the Gauss points with the shape functions
       *
       * Fibers are projected with the shape functions and then normalized to unit length. A
       * RAD-AXI-CIR coordinate system is orthogonalized with preserving the CIR direction over the
       * AXI direction.
       *
       * \tparam distype distype type of the element used
       * \param nodes list of nodes of the element
       * \param shapefcts list of shape functions evaluated at the gauss points
       * \param gpFiberHolder output of the projected fibers and angles
       */
      template <DRT::Element::DiscretizationType distype>
      void ProjectFibersToGaussPoints(DRT::Node** nodes,
          const std::vector<LINALG::Matrix<
              CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>>& shapefcts,
          FIBER::NodalFiberHolder& gpFiberHolder);

      /*!
       * @brief Projects an std::array<double, n> from the nodes with the shape functions into the
       * inner of the element
       *
       * @tparam distype type of the element used
       * @tparam dim dimension of the elemenet
       * @param quantity array of arrays to be projected with the shape functions. The first index
       * are the points, the second index the dimension of the quantity
       * @param shapefcts Shape functions
       * @param quantityProjected Output quantity projected with the shape functions
       */
      template <DRT::Element::DiscretizationType distype, std::size_t dim>
      void ProjectQuantityWithShapeFunctions(
          const std::array<std::array<double, dim>,
              CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
              quantity,
          const std::vector<LINALG::Matrix<
              CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>>& shapefcts,
          std::vector<LINALG::Matrix<dim, 1>>& quantityProjected);

      /*!
       * @brief Projectsd a double from the nodes with the shapefunctions into the innder of the
       * element
       *
       * @tparam distype type of the element used
       * @param quantity array of doubles to be projected with the shape functions.
       * @param shapefcts Shape functions
       * @param quantityProjected Output quantity projected with the shape functions
       */
      template <DRT::Element::DiscretizationType distype>
      void ProjectQuantityWithShapeFunctions(
          const std::array<double,
              CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
              quantity,
          const std::vector<LINALG::Matrix<
              CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>>& shapefcts,
          std::vector<double>& quantityProjected);

      /*!
       * @brief Setup cardiac fibers. If the fiber direction fiber1 is directly given, it is
       * returned. Otherwise, the fiber will be calculated with the RAD-AXI-CIR coordinate system
       * and the helix and transversal angles
       *
       * @tparam dim space dimension of the problem
       * @param fibers Fiber data projected
       * @param f Cardiac fiber direction setup from the fiber or coordinate system data
       */
      template <std::size_t dim>
      void SetupCardiacFibers(
          const NodalFiberHolder& fibers, std::vector<LINALG::Matrix<dim, 1>>& f);

      /*!
       * \brief Check whether all nodes of the element have fibers.
       *
       * \tparam distype Discretization type
       * \param nodes Pointer to the nodes of the element
       * \return true All nodes have fibers
       * \return false At least one one does not have a fiber
       */
      template <DRT::Element::DiscretizationType distype>
      bool HaveNodalFibers(DRT::Node** nodes);
    }  // namespace UTILS
  }    // namespace FIBER
}  // namespace DRT

#endif
