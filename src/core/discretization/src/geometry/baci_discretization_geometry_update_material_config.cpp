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

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::GEO::UpdateMaterialConfigWithDispVector(
    Teuchos::RCP<const DRT::Discretization> dis, Teuchos::RCP<const Epetra_Vector> disp)
{
  // Create Vector which holds all col-displacments of processor
  auto coldisp = Teuchos::rcp(new Epetra_Vector(*dis->DofColMap()));

  // Export row-displacments to col-displacements
  CORE::LINALG::Export(*disp, *coldisp);

  // loop over all nodes
  for (const auto& mynode : dis->MyColNodeRange())
  {
    std::vector<int> globaldofs = dis->Dof(0, mynode);
    // Since ChangePos requires vector of length 3, we use hardcoded length here
    // init with zero just to be sure
    std::vector<double> nvector(3, 0.0);

    for (unsigned int i = 0; i < globaldofs.size(); ++i)
    {
      const int lid = coldisp->Map().LID(globaldofs[i]);

      if (lid < 0)
        dserror(
            "Proc %d: Cannot find gid=%d in Epetra_Vector", coldisp->Comm().MyPID(), globaldofs[i]);
      nvector[i] += (*coldisp)[lid];
    }
    // changepos takes vector with length = 3
    mynode->ChangePos(nvector);
  }
}

BACI_NAMESPACE_CLOSE
