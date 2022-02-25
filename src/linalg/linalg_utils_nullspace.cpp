/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for the nullspace calculation at
       node level

\level 0


*/
/*---------------------------------------------------------------------*/

#include "linalg_utils_nullspace.H"
#include "../drt_s8/shell8.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix LINALG::ComputeSolid3DNullSpace(
    const DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
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

  const double* x = node.X();

  if (numdof != 3)
    dserror(
        "The computation of the solid nullspace in three dimensions requires three DOFs"
        "per solid node, however the current node carries %d DOFs.",
        numdof);

  if (dimnsp != 6)
    dserror(
        "The computation of the solid nullspace in three dimensions requires six nullspace"
        "vectors per node, however the current node carries %d vectors.",
        numdof);

  // this is new old code!
  Epetra_SerialDenseMatrix nullspace = Epetra_SerialDenseMatrix(numdof, dimnsp);
  // x-modes
  nullspace(0, 0) = 1.0;
  nullspace(0, 1) = 0.0;
  nullspace(0, 2) = 0.0;
  nullspace(0, 3) = 0.0;
  nullspace(0, 4) = x[2] - x0[2];
  ;
  nullspace(0, 5) = -x[1] + x0[1];
  // y-modes
  nullspace(1, 0) = 0.0;
  nullspace(1, 1) = 1.0;
  nullspace(1, 2) = 0.0;
  nullspace(1, 3) = -x[2] + x0[2];
  nullspace(1, 4) = 0.0;
  nullspace(1, 5) = x[0] - x0[0];
  // z-modes
  nullspace(2, 0) = 0.0;
  nullspace(2, 1) = 0.0;
  nullspace(2, 2) = 1.0;
  nullspace(2, 3) = x[1] - x0[1];
  nullspace(2, 4) = -x[0] + x0[0];
  nullspace(2, 5) = 0.0;

  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix LINALG::ComputeSolid2DNullSpace(
    const DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  /* the rigid body modes for structures are:

        xtrans   ytrans   zrot
        mode[0]  mode[1]  mode[3]
      ----------------------------
  x   |    1       0       -y+y0
  y   |    0       1       x-x0

  valid element types: wall1, ale2, torsion2

   */

  const double* x = node.X();

  if (numdof != 2)
    dserror(
        "The computation of the solid nullspace in two dimensions requires two DOFs"
        "per solid node, however the current node carries %d DOFs.",
        numdof);

  if (dimnsp != 3)
    dserror(
        "The computation of the solid nullspace in two dimensions requires three nullspace"
        "vectors per node, however the current node carries %d vectors.",
        numdof);

  // this is the new code!
  Epetra_SerialDenseMatrix nullspace = Epetra_SerialDenseMatrix(numdof, dimnsp);
  // x-modes
  nullspace(0, 0) = 1.0;
  nullspace(0, 1) = 0.0;
  nullspace(0, 2) = -x[1] + x0[1];
  // y-modes
  nullspace(1, 0) = 0.0;
  nullspace(1, 1) = 1.0;
  nullspace(1, 2) = x[0] - x0[0];

  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix LINALG::ComputeShell3DNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
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

  Epetra_SerialDenseMatrix dir(1, 3);

  DRT::ELEMENTS::Shell8* s8 = dynamic_cast<DRT::ELEMENTS::Shell8*>(node.Elements()[0]);
  if (!s8) dserror("Cannot cast to Shell8");
  int j;
  for (j = 0; j < s8->NumNode(); ++j)
    if (s8->Nodes()[j]->Id() == node.Id()) break;
  if (j == s8->NumNode()) dserror("Can't find matching node - weird!");
  double h2 = (*s8->GetThickness())[j] / 2.0;

  // get director
  const Epetra_SerialDenseMatrix* a3ref = s8->GetDirectors();
  dir(0, 0) = (*a3ref)(0, j) * h2;
  dir(0, 1) = (*a3ref)(1, j) * h2;
  dir(0, 2) = (*a3ref)(2, j) * h2;

  const double* x = node.X();

  if (numdof != 6)
    dserror(
        "The computation of the shell nullspace in three dimensions requires six DOFs"
        "per node, however the current node carries %d DOFs.",
        numdof);

  if (dimnsp != 6)
    dserror(
        "The computation of the shell nullspace in three dimensions requires six nullspace"
        "vectors per node, however the current node carries %d vectors.",
        numdof);


  // this is the new code!
  Epetra_SerialDenseMatrix nullspace = Epetra_SerialDenseMatrix(numdof, dimnsp);
  // x-modes
  nullspace(0, 0) = 1.0;
  nullspace(0, 1) = 0.0;
  nullspace(0, 2) = 0.0;
  nullspace(0, 3) = 0.0;
  nullspace(0, 4) = x[2] - x0[2];
  nullspace(0, 5) = -x[1] + x0[1];
  // y-modes
  nullspace(1, 0) = 0.0;
  nullspace(1, 1) = 1.0;
  nullspace(1, 2) = 0.0;
  nullspace(1, 3) = -x[2] + x0[2];
  nullspace(1, 4) = 0.0;
  nullspace(1, 5) = x[0] - x0[0];
  // z-modes
  nullspace(2, 0) = 0.0;
  nullspace(2, 1) = 0.0;
  nullspace(2, 2) = 1.0;
  nullspace(2, 3) = x[1] - x0[1];
  nullspace(2, 4) = -x[0] + x0[0];
  nullspace(2, 5) = 0.0;
  // dx-modes
  nullspace(3, 0) = 0.0;
  nullspace(3, 1) = 0.0;
  nullspace(3, 2) = 0.0;
  nullspace(3, 3) = 0.0;
  nullspace(3, 4) = dir(0, 2);
  nullspace(3, 5) = -dir(0, 1);
  // dy-modes
  nullspace(4, 0) = 0.0;
  nullspace(4, 1) = 0.0;
  nullspace(4, 2) = 0.0;
  nullspace(4, 3) = -dir(0, 2);
  nullspace(4, 4) = 0.0;
  nullspace(4, 5) = dir(0, 0);
  // dz-modes
  nullspace(5, 0) = 0.0;
  nullspace(5, 1) = 0.0;
  nullspace(5, 2) = 0.0;
  nullspace(5, 3) = dir(0, 1);
  nullspace(5, 4) = -dir(0, 0);
  nullspace(5, 5) = 0.0;

  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix LINALG::ComputeFluidDNullSpace(
    const DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
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

  if (numdof > 10) dserror("Cannot define more than 10 degrees of freedom!");

  Epetra_SerialDenseMatrix nullspace = Epetra_SerialDenseMatrix(numdof, dimnsp);
  for (int i = 0; i < numdof; i++)
  {
    for (int j = 0; j < dimnsp; j++)
    {
      if (i == j)
        nullspace(i, j) = 1.0;
      else
        nullspace(i, j) = 0.0;
    }
  }

  return nullspace;
}
