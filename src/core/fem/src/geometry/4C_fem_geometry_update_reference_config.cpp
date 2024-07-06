/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace Discret

\level 1


*/
/*---------------------------------------------------------------------*/


#include "4C_fem_geometry_update_reference_config.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::Geo::update_reference_config_with_disp(
    Teuchos::RCP<const Core::FE::Discretization> dis, Teuchos::RCP<const Epetra_Vector> disp)
{
  // Export row-displacments to col-displacements
  auto coldisp = Teuchos::rcp(new Epetra_Vector(*dis->dof_col_map()));
  Core::LinAlg::Export(*disp, *coldisp);

  for (const auto& mynode : dis->my_col_node_range())
  {
    const unsigned int ndim = mynode->n_dim();

#ifdef DEBUG
    FOUR_C_ASSERT(ndim * dis->NodeRowMap()->NumGlobalElements() == disp->Map().NumGlobalElements(),
        "Number of space dimensions does not fit to displacement vector.");

    for (int disp_lid = 0; disp_lid < disp->Map().NumMyElements(); ++disp_lid)
    {
      const int disp_gid = disp->Map().GID(disp_lid);
      FOUR_C_ASSERT(
          dis->dof_row_map()->LID(disp_gid) >= 0, "Displacement dofs not part of dof_row_map()");
    }
#endif

    const auto globaldofs = dis->dof(0, mynode);

    std::vector<double> nvector(ndim, 0.0);

    for (unsigned int i = 0; i < ndim; ++i)
    {
      const int gid = globaldofs[0] + static_cast<int>(i);
      const int lid = coldisp->Map().LID(gid);

      FOUR_C_ASSERT(lid >= 0, "Proc %d: Cannot find gid=%d in Epetra_Vector",
          coldisp->Comm().MyPID(), globaldofs[i]);

      nvector[i] = (*coldisp)[lid];
    }

    mynode->change_pos(nvector);
  }
}

FOUR_C_NAMESPACE_CLOSE
