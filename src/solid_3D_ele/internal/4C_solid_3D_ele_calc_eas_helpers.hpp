/*! \file

\brief Declaration of routines for calculation of solid element with EAS element technology

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_EAS_HELPERS_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_EAS_HELPERS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype>
  inline static constexpr int num_nodes = Core::FE::num_nodes<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dim = Core::FE::dim<celltype>;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_str = num_dim<celltype>*(num_dim<celltype> + 1) / 2;

  template <Core::FE::CellType celltype>
  inline static constexpr int num_dof_per_ele = num_nodes<celltype>* num_dim<celltype>;

  template <Core::FE::CellType celltype>
  struct CentroidTransformation
  {
    // transformation matrix T0^{-T}, which maps the matrix M from parameter space to the material
    // configuration see Andelfinger et al., EAS-elements, 1993, doi: 10.1002/nme.1620360805
    Core::LinAlg::Matrix<num_str<celltype>, num_str<celltype>> T0invT_;

    // Jacobi determinant evaluated at the element centroid
    double detJ0_;
  };
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
