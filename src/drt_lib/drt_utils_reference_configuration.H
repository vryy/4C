/*!----------------------------------------------------------------------
 \file drt_utils_reference_configuration.H

 \brief Utility functions to calculate properties in the reference configuration

 \level 2

 \maintainer Christoph Schmidt
 http://www.lnm.mw.tum.de
 089 - 289-15251

 *----------------------------------------------------------------------*/

#ifndef DRT_UTILS_REFERENCE_CONFIGURATION
#define DRT_UTILS_REFERENCE_CONFIGURATION

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_node.H"

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
        const LINALG::Matrix<DisTypeToDim<distype>::dim, 1>& xi, LINALG::Matrix<probdim, 1>& coord)
    {
      static LINALG::Matrix<DisTypeToNumNodePerEle<distype>::numNodePerElement, 1> funct(true);
      static LINALG::Matrix<probdim, DisTypeToNumNodePerEle<distype>::numNodePerElement> nodecoords(
          true);

      shape_function<distype>(xi, funct);

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

      for (int i = 0; i < DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
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
        const LINALG::Matrix<DisTypeToDim<distype>::dim, 1>& xi, LINALG::Matrix<3, 1>& normal)
    {
      static LINALG::Matrix<3, DisTypeToDim<distype>::dim> gxieta(true);
      static LINALG::Matrix<DisTypeToDim<distype>::dim,
          DisTypeToNumNodePerEle<distype>::numNodePerElement>
          deriv(true);
      static LINALG::Matrix<3, DisTypeToNumNodePerEle<distype>::numNodePerElement> nodecoords(true);

      shape_function_deriv1<distype>(xi, deriv);

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

      for (int i = 0; i < DisTypeToNumNodePerEle<distype>::numNodePerElement; ++i)
      {
        for (int j = 0; j < nodedim; ++j)
        {
          nodecoords(j, i) = nodes[i]->X()[j];
        }
      }

      gxieta.MultiplyNT(1.0, nodecoords, deriv, 0.0);
      static LINALG::Matrix<3, 1> gxi(true);
      static LINALG::Matrix<3, 1> geta(true);
      static LINALG::Matrix<2, 1> first(true);
      static LINALG::Matrix<2, 1> second(true);
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

#endif  // #ifdef DRT_UTILS_REFERENCE_CONFIGURATION
