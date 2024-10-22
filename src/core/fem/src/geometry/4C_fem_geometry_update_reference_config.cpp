// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_geometry_update_reference_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Geo::update_reference_config_with_disp(
    const Core::FE::Discretization& dis, const Core::LinAlg::Vector<double>& disp)
{
  // Export row-displacments to col-displacements
  Core::LinAlg::Vector<double> coldisp(*dis.dof_col_map());
  Core::LinAlg::export_to(disp, coldisp);

  for (const auto& mynode : dis.my_col_node_range())
  {
    const unsigned int ndim = mynode->n_dim();

#ifdef FOUR_C_ENABLE_ASSERTIONS
    FOUR_C_ASSERT(static_cast<int>(ndim * dis.node_row_map()->NumGlobalElements()) ==
                      disp.Map().NumGlobalElements(),
        "Number of space dimensions does not fit to displacement vector.");

    for (int disp_lid = 0; disp_lid < disp.Map().NumMyElements(); ++disp_lid)
    {
      const int disp_gid = disp.Map().GID(disp_lid);
      FOUR_C_ASSERT(
          dis.dof_row_map()->LID(disp_gid) >= 0, "Displacement dofs not part of dof_row_map()");
    }
#endif

    const auto globaldofs = dis.dof(0, mynode);

    std::vector<double> nvector(ndim, 0.0);

    for (unsigned int i = 0; i < ndim; ++i)
    {
      const int gid = globaldofs[0] + static_cast<int>(i);
      const int lid = coldisp.Map().LID(gid);

      FOUR_C_ASSERT(lid >= 0, "Proc %d: Cannot find gid=%d in Core::LinAlg::Vector<double>",
          coldisp.Comm().MyPID(), globaldofs[i]);

      nvector[i] = (coldisp)[lid];
    }

    mynode->change_pos(nvector);
  }
}

FOUR_C_NAMESPACE_CLOSE
