/*---------------------------------------------------------------------*/
/*! \file

\brief Collection of general element utility functions


\level 1

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_ELEMENTS_JACOBIAN_HPP
#define FOUR_C_FEM_GENERAL_ELEMENTS_JACOBIAN_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  /// \brief get the minimal Jacobian determinant value calculated at the node positions
  /**
   *  \param xcurr  Current nodal positions
   *
   *  \return minimal value.
   *
   *  \author hiermeier \date 09/18 */
  template <Core::FE::CellType type, unsigned numnode, unsigned numdim>
  double GetMinimalJacDeterminantAtNodes(const Core::LinAlg::Matrix<numdim, numnode>& xcurr)
  {
    // check input
    static_assert(numnode == Core::FE::num_nodes<type>, "Wrong matrix dimension.");
    static_assert(numdim == Core::FE::dim<type>, "Wrong matrix dimension.");

    Core::LinAlg::Matrix<numdim, numnode> deriv_at_c(false);
    Core::LinAlg::Matrix<numdim, numdim> jac_at_c(false);

    // parametric coordinates of the HEX8 corners
    static const Core::LinAlg::SerialDenseMatrix rst =
        Core::FE::getEleNodeNumbering_nodes_paramspace(type);
    double min_detJ = std::numeric_limits<double>::max();

    for (unsigned c = 0; c < numnode; ++c)
    {
      const Core::LinAlg::Matrix<numdim, 1> rst_c(&rst(0, c), true);

      Core::FE::shape_function_deriv1<type>(rst_c, deriv_at_c);
      jac_at_c.multiply_nt(deriv_at_c, xcurr);

      const double detJ_at_c = jac_at_c.determinant();
      if (detJ_at_c < min_detJ) min_detJ = detJ_at_c;
    }

    return min_detJ;
  }
}  // namespace Core::Elements

FOUR_C_NAMESPACE_CLOSE

#endif
