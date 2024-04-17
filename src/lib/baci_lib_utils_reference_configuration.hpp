/*---------------------------------------------------------------------*/
/*! \file

 \brief Utility functions to calculate properties in the reference configuration

 \level 2


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_UTILS_REFERENCE_CONFIGURATION_HPP
#define FOUR_C_LIB_UTILS_REFERENCE_CONFIGURATION_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace UTILS
  {
    /**
     * \brief Calculate the global position in the reference configuration for a given element at a
     * given position in parameter space \f$ \vec{\xi} \f$
     *
     *  \param [in]   element: element for which position shall be calculated
     *  \param [in]        xi: position in parameter space
     *  \param [in,out] coord: global position in reference configuration for element at parameter
     * space position xi
     *
     *  \author cschmidt \date 11/18 */
    template <int probdim, CORE::FE::CellType distype>
    static void LocalToGlobalPositionAtXiRefConfig(const DRT::Element* element,
        const CORE::LINALG::Matrix<CORE::FE::dim<distype>, 1>& xi,
        CORE::LINALG::Matrix<probdim, 1>& coord)
    {
      static CORE::LINALG::Matrix<CORE::FE::num_nodes<distype>, 1> funct(true);
      static CORE::LINALG::Matrix<probdim, CORE::FE::num_nodes<distype>> nodecoords(true);

      CORE::FE::shape_function<distype>(xi, funct);

      const DRT::Node* const* nodes = element->Nodes();
      const int nodedim = nodes[0]->Dim();

      if (!nodes)
      {
        dserror("ERROR: Did not get nodes of element!");
      }
      if (probdim != nodedim)
      {
        dserror(
            "Problem dimension: %i and dimension of nodes: %i does not match!", probdim, nodedim);
      }

      for (int i = 0; i < CORE::FE::num_nodes<distype>; ++i)
      {
        for (int j = 0; j < nodedim; ++j)
        {
          nodecoords(j, i) = nodes[i]->X()[j];
        }
      }

      coord.Multiply(1.0, nodecoords, funct, 0.0);
    }

    /**
     * \brief Calculate the normal in the reference configuration for a given element at a given
     * position in parameter space \f$ \vec{\xi} \f$
     *
     *  \param [in]    element: element for which position shall be calculated
     *  \param [in]         xi: position in parameter space
     *  \param [in,out] normal: normal in reference configuration for element at parameter space
     * position xi
     *
     *  \author cschmidt \date 11/18 */
    template <CORE::FE::CellType distype>
    static void ComputeUnitNormalAtXiRefConfig(const DRT::Element* element,
        const CORE::LINALG::Matrix<CORE::FE::dim<distype>, 1>& xi,
        CORE::LINALG::Matrix<3, 1>& normal)
    {
      static CORE::LINALG::Matrix<3, CORE::FE::dim<distype>> gxieta(true);
      static CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>> deriv(true);
      static CORE::LINALG::Matrix<3, CORE::FE::num_nodes<distype>> nodecoords(true);

      CORE::FE::shape_function_deriv1<distype>(xi, deriv);

      const DRT::Node* const* nodes = element->Nodes();
      const int nodedim = nodes[0]->Dim();

      if (!nodes)
      {
        dserror("ERROR: Did not get nodes of element!");
      }
      if (nodedim != 3)
      {
        dserror("ERROR: Only implemented for 3D cases so far!");
      }

      for (int i = 0; i < CORE::FE::num_nodes<distype>; ++i)
      {
        for (int j = 0; j < nodedim; ++j)
        {
          nodecoords(j, i) = nodes[i]->X()[j];
        }
      }

      gxieta.MultiplyNT(1.0, nodecoords, deriv, 0.0);
      static CORE::LINALG::Matrix<3, 1> gxi(true);
      static CORE::LINALG::Matrix<3, 1> geta(true);
      static CORE::LINALG::Matrix<2, 1> first(true);
      static CORE::LINALG::Matrix<2, 1> second(true);
      first(0, 0) = 1.0;
      second(1, 0) = 1.0;
      gxi.Multiply(1.0, gxieta, first, 0.0);
      geta.Multiply(1.0, gxieta, second, 0.0);

      // clear, calculate and scale normal
      normal.Clear();
      normal.CrossProduct(gxi, geta);
      const double normnormal = normal.Norm2();
      normal.Scale(1.0 / normnormal);
    }
  }  // namespace UTILS
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
