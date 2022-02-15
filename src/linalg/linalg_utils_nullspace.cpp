/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for namespace DRT

\level 0


*/
/*---------------------------------------------------------------------*/

#include "linalg_utils_nullspace.H"
#include "../drt_s8/shell8.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeStructure3DNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  /* the rigid body modes for structures are:

        xtrans   ytrans  ztrans   xrot       yrot       zrot
        mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
    -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0

  valid element types: ale3, so_hex8, so_hex20, so_hex27, sosh8, so_tet4,
                       so_tet10, so_weg6, sodisp, so_shw6, truss3, torsion3
  */

  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const double* x = actnode->X();
    std::vector<int> dofs = dis.Dof(0, actnode);  // use current dofset

    if (dofs.size() != 3)
      dserror(
          "The computation of the solid nullspace in three dimensions requires three DOFs"
          "per solid node, however the current node carries %d DOFs.",
          dofs.size());

    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find dof");
      switch (j)  // j is degree of freedom
      {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
          break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
          break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
          break;
        default:
          dserror("Only dofs 0 - 5 supported");
          break;
      }  // switch (j)
    }    // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  }      // for (int i=0; i<NumMyRowNodes(); ++i)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeStructure2DNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  /* the rigid body modes for structures are:

        xtrans   ytrans   zrot
        mode[0]  mode[1]  mode[3]
      ----------------------------
  x   |    1       0       -y+y0
  y   |    0       1       x-x0

  valid element types: wall1, ale2, torsion2

   */

  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const double* x = actnode->X();
    std::vector<int> dofs = dis.Dof(0, actnode);

    if (dofs.size() != 2)
      dserror(
          "The computation of the solid nullspace in two dimensions requires two DOFs"
          "per solid node, however the current node carries %d DOFs.",
          dofs.size());

    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find dof");
      switch (j)  // j is degree of freedom
      {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = -x[1] + x0[1];
          break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = x[0] - x0[0];
          break;
        default:
          dserror("Only dofs 0 - 1 supported");
          break;
      }  // switch (j)
    }    // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  }      // for (int i=0; i<NumMyRowNodes(); ++i)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeShell3DNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  /* the rigid body modes for structures are:
      xtrans   ytrans  ztrans   xrot       yrot       zrot
      mode[0]  mode[1] mode[2]  mode[3]    mode[4]    mode[5]
    -----------------------------------------------------------
  x   |    1       0       0       0          z-z0      -y+y0
  y   |    0       1       0      -z+z0       0          x-x0
  z   |    0       0       1       y-y0      -x+x0       0
  dx  |    0       0       0       0          a3        -a2
  dy  |    0       0       0      -a3         0          a1
  dz  |    0       0       0       a2        -a1         0

  element types: shell8

   */

  const Epetra_Map* rowmap = dis.DofRowMap(0);
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  Epetra_SerialDenseMatrix dir;

  dir.Shape(dis.NumMyRowNodes(), 3);
  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    DRT::ELEMENTS::Shell8* s8 = dynamic_cast<DRT::ELEMENTS::Shell8*>(actnode->Elements()[0]);
    if (!s8) dserror("Cannot cast to Shell8");
    int j;
    for (j = 0; j < s8->NumNode(); ++j)
      if (s8->Nodes()[j]->Id() == actnode->Id()) break;
    if (j == s8->NumNode()) dserror("Can't find matching node - weird!");
    double h2 = (*s8->GetThickness())[j] / 2.0;
    // get director
    const Epetra_SerialDenseMatrix* a3ref = s8->GetDirectors();
    dir(i, 0) = (*a3ref)(0, j) * h2;
    dir(i, 1) = (*a3ref)(1, j) * h2;
    dir(i, 2) = (*a3ref)(2, j) * h2;
  }

  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    const double* x = actnode->X();
    std::vector<int> dofs = dis.Dof(0, actnode);  // use current dofset

    if (dofs.size() != 6)
      dserror(
          "The computation of the shell nullspace in three dimensions requires six DOFs"
          "per solid node, however the current node carries %d DOFs.",
          dofs.size());

    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find dof");
      switch (j)  // j is degree of freedom
      {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = x[2] - x0[2];
          mode[5][lid] = -x[1] + x0[1];
          break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -x[2] + x0[2];
          mode[4][lid] = 0.0;
          mode[5][lid] = x[0] - x0[0];
          break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = x[1] - x0[1];
          mode[4][lid] = -x[0] + x0[0];
          mode[5][lid] = 0.0;
          break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          mode[4][lid] = dir(i, 2);
          mode[5][lid] = -dir(i, 1);
          break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = -dir(i, 2);
          mode[4][lid] = 0.0;
          mode[5][lid] = dir(i, 0);
          break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = dir(i, 1);
          mode[4][lid] = -dir(i, 0);
          mode[5][lid] = 0.0;
          break;
        default:
          dserror("Only dofs 0 - 5 supported");
          break;
      }  // switch (j)
    }    // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  }      // for (int i=0; i<NumMyRowNodes(); ++i)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeXFluidDNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  /* the rigid body modes for fluids are:

        xtrans   ytrans  ztrans   pressure
        mode[0]  mode[1] mode[2]  mode[3]
    ----------------------------------------
  x   |    1       0       0       0
  y   |    0       1       0       0
  z   |    0       0       1       0
  p   |    0       0       0       1

  valid element types: fluid3, xfluid3

  */

  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[6];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    std::vector<int> dofs = dis.Dof(actnode);
    for (unsigned j = 0; j < dofs.size(); ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find dof");
      switch (j)  // j is degree of freedom
      {
        case 0:
          mode[0][lid] = 1.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        case 1:
          mode[0][lid] = 0.0;
          mode[1][lid] = 1.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        case 2:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 1.0;
          mode[3][lid] = 0.0;
          break;
        case 3:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 1.0;
          break;
        case 4:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        case 5:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        case 6:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        case 7:
          mode[0][lid] = 0.0;
          mode[1][lid] = 0.0;
          mode[2][lid] = 0.0;
          mode[3][lid] = 0.0;
          break;
        default:
          dserror("Only dofs 0 - 7 supported");
          break;
      }  // switch (j)
    }    // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  }      // for (int i=0; i<NumMyRowNodes(); ++i)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::ComputeFluidDNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  const Epetra_Map* rowmap = dis.DofRowMap();
  const int lrows = rowmap->NumMyElements();
  double* mode[10];
  for (int i = 0; i < dimns; ++i) mode[i] = &(ns[i * lrows]);

  for (int i = 0; i < dis.NumMyRowNodes(); ++i)
  {
    DRT::Node* actnode = dis.lRowNode(i);
    std::vector<int> dofs = dis.Dof(0, actnode);
    const unsigned int ndof = dofs.size();

    if (numdf > 10) dserror("Cannot define more than 10 modes");
    for (unsigned j = 0; j < ndof; ++j)
    {
      const int dof = dofs[j];
      const int lid = rowmap->LID(dof);
      if (lid < 0) dserror("Cannot find dof");

      for (unsigned k = 0; k < ndof; ++k)
      {
        if (k % numdf == j % numdf)
          mode[k % numdf][lid] = 1.0;
        else
          mode[k % numdf][lid] = 0.0;
      }
    }  // for (int j=0; j<actnode->Dof().NumDof(); ++j)
  }    // for (int i=0; i<NumMyRowNodes(); ++i)
}