// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_GEOMETRY_INTEGRATIONCELL_COORDTRAFO_HPP
#define FOUR_C_FEM_GEOMETRY_INTEGRATIONCELL_COORDTRAFO_HPP


#include "4C_config.hpp"

#include "4C_fem_geometry_element_coordtrafo.hpp"
#include "4C_fem_geometry_integrationcell.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::Geo
{
  //! map position from eta^boundary to xi^domain space
  inline void map_eta_b_to_xi_d(const Core::Geo::BoundaryIntCell& cell,
      const Core::LinAlg::Matrix<2, 1>& pos_eta_boundary,
      Core::LinAlg::Matrix<3, 1>& pos_xsi_domain)
  {
    // get cell node coordinates in xi_domain
    const Core::LinAlg::SerialDenseMatrix& xyze_cell(cell.cell_nodal_pos_xi_domain());
    Core::Geo::element_to_current_coordinates(
        cell.shape(), xyze_cell, pos_eta_boundary, pos_xsi_domain);
    return;
  }
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
