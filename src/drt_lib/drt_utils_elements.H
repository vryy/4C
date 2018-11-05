/*----------------------------------------------------------------------------*/
/*!
\file drt_utils_elements.H

\brief Collection of general element utility functions

\maintainer Michael Hiermeier

\date Sep 27, 2018

\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_DRT_UTILS_ELEMENTS_H_
#define SRC_DRT_LIB_DRT_UTILS_ELEMENTS_H_

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

namespace DRT
{
  namespace UTILS
  {
    /// \brief get the minimal Jacobian determinant value calculated at the node positions
    /**
     *  \param xcurr  Current nodal positions
     *
     *  \return minimal value.
     *
     *  \author hiermeier \date 09/18 */
    template <DRT::Element::DiscretizationType type, unsigned numnode, unsigned numdim>
    double GetMinimalJacDeterminantAtNodes(const LINALG::Matrix<numdim, numnode>& xcurr)
    {
      // check input
      static_assert(
          numnode == DisTypeToNumNodePerEle<type>::numNodePerElement, "Wrong matrix dimension.");
      static_assert(numdim == DisTypeToDim<type>::dim, "Wrong matrix dimension.");

      LINALG::Matrix<numdim, numnode> deriv_at_c(false);
      LINALG::Matrix<numdim, numdim> jac_at_c(false);

      // parametric coordinates of the HEX8 corners
      static const LINALG::SerialDenseMatrix rst =
          DRT::UTILS::getEleNodeNumbering_nodes_paramspace(type);
      double min_detJ = std::numeric_limits<double>::max();

      for (unsigned c = 0; c < numnode; ++c)
      {
        const LINALG::Matrix<numdim, 1> rst_c(&rst(0, c), true);

        DRT::UTILS::shape_function_deriv1<type>(rst_c, deriv_at_c);
        jac_at_c.MultiplyNT(deriv_at_c, xcurr);

        const double detJ_at_c = jac_at_c.Determinant();
        if (detJ_at_c < min_detJ) min_detJ = detJ_at_c;
      }

      return min_detJ;
    }
  }  // namespace UTILS
}  // namespace DRT


#endif /* SRC_DRT_LIB_DRT_UTILS_ELEMENTS_H_ */
