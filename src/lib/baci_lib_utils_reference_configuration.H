/*---------------------------------------------------------------------*/
/*! \file

 \brief Utility functions to calculate properties in the reference configuration

 \level 2


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_LIB_UTILS_REFERENCE_CONFIGURATION_H
#define BACI_LIB_UTILS_REFERENCE_CONFIGURATION_H

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_lib_node.H"
#include "baci_linalg_fixedsizematrix.H"

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
    template <int probdim, DRT::Element::DiscretizationType distype>
    static void LocalToGlobalPositionAtXiRefConfig(const DRT::Element* element,
        const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim, 1>& xi,
        CORE::LINALG::Matrix<probdim, 1>& coord)
    {
      static CORE::LINALG::Matrix<
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement, 1>
          funct(true);
      static CORE::LINALG::Matrix<probdim,
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
          nodecoords(true);

      CORE::DRT::UTILS::shape_function<distype>(xi, funct);

      const DRT::Node* const* nodes = element->Nodes();
      const int nodedim = nodes[0]->Dim();

      if (!nodes)
      {
        dserror("ERROR: Did not get nodes of element!");
      }
      if (probdim != nodedim)
      {
        dserror("Problem dimension and dimension of nodes does not match!");
      }

      for (int i = 0; i < CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
      {
        for (int j = 0; j < nodedim; ++j)
        {
          nodecoords(j, i) = nodes[i]->X()[j];
        }
      }

      coord.Multiply(1.0, nodecoords, funct, 0.0);

      return;
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
    template <DRT::Element::DiscretizationType distype>
    static void ComputeUnitNormalAtXiRefConfig(const DRT::Element* element,
        const CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim, 1>& xi,
        CORE::LINALG::Matrix<3, 1>& normal)
    {
      static CORE::LINALG::Matrix<3, CORE::DRT::UTILS::DisTypeToDim<distype>::dim> gxieta(true);
      static CORE::LINALG::Matrix<CORE::DRT::UTILS::DisTypeToDim<distype>::dim,
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
          deriv(true);
      static CORE::LINALG::Matrix<3,
          CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement>
          nodecoords(true);

      CORE::DRT::UTILS::shape_function_deriv1<distype>(xi, deriv);

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

      for (int i = 0; i < CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
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

      return;
    }
  }  // namespace UTILS
}  // namespace DRT

#endif  // LIB_UTILS_REFERENCE_CONFIGURATION_H
