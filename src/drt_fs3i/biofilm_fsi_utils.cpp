/*----------------------------------------------------------------------*/
/*! \file

\brief utils for biofilm fs3i

\level 3


 *----------------------------------------------------------------------*/

#include "biofilm_fsi_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BIOFILM::UTILS::ScatraChangeConfig(Teuchos::RCP<DRT::Discretization> scatradis,
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Vector> disp)
{
  const int numnode = (scatradis->NodeColMap())->NumMyElements();

  // Create Vector which holds all col-displacments of processor
  Teuchos::RCP<Epetra_Vector> coldisp = Teuchos::rcp(new Epetra_Vector(*(dis->DofColMap())));

  // Export row-displacments to col-displacements
  LINALG::Export(*disp, *coldisp);


  const Epetra_Vector& gvector = *coldisp;

  // loop over all nodes
  for (int index = 0; index < numnode; ++index)
  {
    // get current node
    int gid = (scatradis->NodeColMap())->GID(index);
    DRT::Node* mynode = scatradis->gNode(gid);

    // get local fluid/structure node with the same local node id
    DRT::Node* lnode = dis->lColNode(index);

    // get degrees of freedom associated with this fluid/structure node
    std::vector<int> nodedofs = dis->Dof(0, lnode);

    // Since ChangePos requires vector of length 3, we use hardcoded length here
    // init with zero just to be sure
    std::vector<double> nvector(3, 0.0);

    // determine number of space dimensions
    const int numdim = DRT::Problem::Instance()->NDim();

    for (int i = 0; i < numdim; ++i)
    {
      const int lid = gvector.Map().LID(nodedofs[i]);

      if (lid < 0)
        dserror(
            "Proc %d: Cannot find gid=%d in Epetra_Vector", gvector.Comm().MyPID(), nodedofs[i]);
      nvector[i] += gvector[lid];
    }

    mynode->ChangePos(nvector);
  }

  return;
}
