/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of calculation routines related to Nullspaces

\level 0

*----------------------------------------------------------------------*/

#include "linalg_nullspace.H"
#include "linalg_utils_sparse_algebra_manipulation.H"
#include "linalg_utils_sparse_algebra_assemble.H"
#include "linalg_utils_sparse_algebra_print.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> LINALG::NULLSPACE::ComputeNullSpace(const DRT::Discretization& dis,
    const int numdf, const int dimns, const Teuchos::RCP<Epetra_Map> dofmap)
{
  if (dimns > 10) dserror("Nullspace size only up to 10 supported!");

  Teuchos::RCP<Epetra_MultiVector> nullspace =
      Teuchos::rcp(new Epetra_MultiVector(*dofmap, dimns, true));

  if (dimns == 1 && numdf == 1)
  {
    // compute nullspace for simple case: vector of ones
    nullspace->PutScalar(1.0);
  }
  else
  {
    // for rigid body rotations compute nodal center of the discretization
    double x0send[3] = {0.0, 0.0, 0.0};
    for (int i = 0; i < dis.NumMyRowNodes(); ++i)
      for (int j = 0; j < 3; ++j) x0send[j] += dis.lRowNode(i)->X()[j];
    double x0[3];
    dis.Comm().SumAll(x0send, x0, 3);
    for (int i = 0; i < 3; ++i) x0[i] /= dis.NumGlobalNodes();

    // assembly process of the nodalNullspace into the actual nullspace
    for (int i = 0; i < dis.NumMyRowNodes(); ++i)
    {
      DRT::Node* actnode = dis.lRowNode(i);
      std::vector<int> dofs = dis.Dof(0, actnode);
      const int localLength = dofs.size();

      // check if degrees of freedom are zero
      if (localLength == 0) continue;

      // check if dof is exisiting as index
      if (dofmap->LID(dofs[0]) == -1) continue;

      // check size of degrees of freedom
      if (localLength != numdf)
      {
        std::cout << "Warning: At local node " << i << " : nullspace degrees of freedom ( " << numdf
                  << " ) "
                  << "and rowmap degrees of freedom ( " << localLength << " ) are not consistent"
                  << std::endl;
      }

      // Here we check the first element type of the node. One node can be owned by several elements
      // we restrict the routine, that a node is only owned by elements with the same physics
      if (actnode->NumElement() > 1)
        for (int i = 0; i < actnode->NumElement() - 1; i++)
          // if element types are different, check nullspace dimension and dofs
          if (actnode->Elements()[i + 1]->ElementType() != actnode->Elements()[i]->ElementType())
          {
            int numdof1, dimnsp1, nv1, np1;
            actnode->Elements()[i]->ElementType().NodalBlockInformation(
                actnode->Elements()[i], numdof1, dimnsp1, nv1, np1);
            int numdof2, dimnsp2, nv2, np2;
            actnode->Elements()[i + 1]->ElementType().NodalBlockInformation(
                actnode->Elements()[i + 1], numdof2, dimnsp2, nv2, np2);

            if (numdof1 != numdof2 || dimnsp1 != dimnsp2)
              dserror("Node is owned by different element types, nullspace calculation aborted!");
          }

      Epetra_SerialDenseMatrix nodalNullspace =
          actnode->Elements()[0]->ElementType().ComputeNullSpace(*actnode, x0, localLength, dimns);

      for (int dim = 0; dim < dimns; ++dim)
      {
        double** arrayOfPointers;
        nullspace->ExtractView(&arrayOfPointers);
        double* data = arrayOfPointers[dim];
        Teuchos::ArrayRCP<double> dataVector(data, dofmap->LID(dofs[0]), localLength, false);

        for (int j = 0; j < localLength; ++j)
        {
          const int lid = dofmap->LID(dofs[j]);
          dataVector[lid] = nodalNullspace(j, dim);
        }
      }
    }
  }

  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::NULLSPACE::FixNullSpace(std::string field, const Epetra_Map& oldmap,
    const Epetra_Map& newmap, Teuchos::ParameterList& solveparams)
{
  if (!oldmap.Comm().MyPID()) printf("Fixing %s Nullspace\n", field.c_str());

  // there is no ML or MueLu list, do nothing
  if (!solveparams.isSublist("ML Parameters") && !solveparams.isSublist("MueLu Parameters")) return;

  // find the ML or MueLu list
  Teuchos::ParameterList* params_ptr = NULL;
  if (solveparams.isSublist("ML Parameters"))
    params_ptr = &(solveparams.sublist("ML Parameters"));
  else if (solveparams.isSublist("MueLu Parameters"))
    params_ptr = &(solveparams.sublist("MueLu Parameters"));
  else
    return;
  Teuchos::ParameterList& params = *params_ptr;

  const int ndim = params.get("null space: dimension", -1);
  if (ndim == -1) dserror("List does not contain nullspace dimension");

  Teuchos::RCP<Epetra_MultiVector> nullspace =
      params.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);
  if (nullspace == Teuchos::null) dserror("List does not contain nullspace");

  const int nullspaceLength = nullspace->MyLength();
  const int newmapLength = newmap.NumMyElements();

  if (nullspaceLength == newmapLength) return;
  if (nullspaceLength != oldmap.NumMyElements())
    dserror("Nullspace map of length %d does not match old map length of %d", nullspaceLength,
        oldmap.NumMyElements());
  if (newmapLength > nullspaceLength)
    dserror("New problem size larger than old - full rebuild of nullspace neccessary");

  Teuchos::RCP<Epetra_MultiVector> nullspaceNew =
      Teuchos::rcp(new Epetra_MultiVector(newmap, ndim, true));

  for (int i = 0; i < ndim; i++)
  {
    Epetra_Vector* nullspaceData = (*nullspace)(i);
    Epetra_Vector* nullspaceDataNew = (*nullspaceNew)(i);
    const size_t myLength = nullspaceDataNew->MyLength();
    for (size_t j = 0; j < myLength; j++)
    {
      int gid = newmap.GID(j);
      int olid = oldmap.LID(gid);
      if (olid == -1) continue;
      (*nullspaceDataNew)[j] = (*nullspaceData)[olid];
    }
  }

  params.set<Teuchos::RCP<Epetra_MultiVector>>("nullspace", nullspaceNew);
  params.set("null space: vectors", nullspaceNew->Values());
}
