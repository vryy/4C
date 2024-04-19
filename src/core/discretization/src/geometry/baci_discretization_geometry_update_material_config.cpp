/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 1


*/
/*---------------------------------------------------------------------*/


#include "baci_discretization_geometry_update_material_config.hpp"

#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::GEO::UpdateMaterialConfigWithDispVector(
    Teuchos::RCP<const DRT::Discretization> dis, Teuchos::RCP<const Epetra_Vector> disp)
{
  // Export row-displacments to col-displacements
  auto coldisp = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap()));
  CORE::LINALG::Export(*disp, *coldisp);

  for (const auto& mynode : dis->MyColNodeRange())
  {
    const unsigned int ndim = mynode->Dim();

#ifdef DEBUG
    FOUR_C_ASSERT(ndim * dis->NodeRowMap()->NumGlobalElements() == disp->Map().NumGlobalElements(),
        "Number of space dimensions does not fit to displacement vector.");

    for (int disp_lid = 0; disp_lid < disp->Map().NumMyElements(); ++disp_lid)
    {
      const int disp_gid = disp->Map().GID(disp_lid);
      FOUR_C_ASSERT(
          dis->DofRowMap()->LID(disp_gid) >= 0, "Displacement dofs not part of DofRowMap()");
    }
#endif

    const auto globaldofs = dis->Dof(0, mynode);

    std::vector<double> nvector(ndim, 0.0);

    for (unsigned int i = 0; i < ndim; ++i)
    {
      const int gid = globaldofs[0] + static_cast<int>(i);
      const int lid = coldisp->Map().LID(gid);

      FOUR_C_ASSERT(lid >= 0, "Proc %d: Cannot find gid=%d in Epetra_Vector",
          coldisp->Comm().MyPID(), globaldofs[i]);

      nvector[i] = (*coldisp)[lid];
    }

    mynode->ChangePos(nvector);
  }
}

FOUR_C_NAMESPACE_CLOSE
