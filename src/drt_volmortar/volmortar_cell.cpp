/*!----------------------------------------------------------------------
\file volmortar_cell.cpp

\level 1

<pre>
\maintainer Alexander Popp
</pre>

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "volmortar_cell.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

/*---------------------------------------------------------------------*
 | constructor                                             farah 01/14 |
 *---------------------------------------------------------------------*/
VOLMORTAR::Cell::Cell(int id, int nvertices, const Epetra_SerialDenseMatrix& coords,
    const DRT::Element::DiscretizationType& shape)
    : id_(id), coords_(coords), shape_(shape)
{
  if (shape_ == DRT::Element::tet4)
  {
    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1    1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **             [ z_1  z_2  z_3  z_4 ]
    */

    LINALG::Matrix<4, 4> jac;
    for (int i = 0; i < 4; i++) jac(0, i) = 1;
    for (int row = 0; row < 3; row++)
      for (int col = 0; col < 4; col++) jac(row + 1, col) = coords_(row, col);

    // volume of the element
    vol_ = jac.Determinant() / 6.0;
    if (vol_ <= 0.0) dserror("Element volume %10.5e <= 0.0", vol_);
  }
  else
    vol_ = 0.0;

  // std::cout << "SHAPE=     "<< shape_ << std::endl;
  if (shape_ != DRT::Element::tet4) dserror("wrong shape");

  return;
}

/*---------------------------------------------------------------------*
 | calculate jacobian for hex elements                     farah 04/14 |
 *---------------------------------------------------------------------*/
double VOLMORTAR::Cell::CalcJac(const double* xi)
{
  double jac = 0.0;

  LINALG::Matrix<3, 8> derivs;
  const double r = xi[0];
  const double s = xi[1];
  const double t = xi[2];

  DRT::UTILS::shape_function_3D_deriv1(derivs, r, s, t, DRT::Element::hex8);


  LINALG::Matrix<8, 3> xrefe;
  for (int i = 0; i < 8; ++i)
  {
    xrefe(i, 0) = coords_(0, i);
    xrefe(i, 1) = coords_(1, i);
    xrefe(i, 2) = coords_(2, i);
  }

  LINALG::Matrix<3, 3> invJ;
  invJ.Clear();

  invJ.Multiply(derivs, xrefe);
  jac = invJ.Invert();
  if (jac <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0", jac);

  return jac;
}
/*---------------------------------------------------------------------*
 | mapping from parameter space to global space            farah 01/14 |
 *---------------------------------------------------------------------*/
void VOLMORTAR::Cell::LocalToGlobal(double* local, double* global)
{
  if (shape_ == DRT::Element::tet4)
  {
    // check input
    if (!local) dserror("ERROR: LocalToGlobal called with xi=NULL");
    if (!global) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

    static const int n = 4;
    static const int ndim = 3;

    for (int i = 0; i < ndim; ++i) global[i] = 0.0;

    LINALG::Matrix<n, 1> val;
    DRT::UTILS::shape_function_3D(val, local[0], local[1], local[2], shape_);

    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < ndim; ++j)
      {
        // use shape function values for interpolation
        global[j] += val(i) * coords_(j, i);
      }
    }
  }
  else if (shape_ == DRT::Element::hex8)
  {
    // check input
    if (!local) dserror("ERROR: LocalToGlobal called with xi=NULL");
    if (!global) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

    static const int n = 8;
    static const int ndim = 3;

    for (int i = 0; i < ndim; ++i) global[i] = 0.0;

    LINALG::Matrix<n, 1> val;
    DRT::UTILS::shape_function_3D(val, local[0], local[1], local[2], shape_);

    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < ndim; ++j)
      {
        // use shape function values for interpolation
        global[j] += val(i) * coords_(j, i);
      }
    }
  }
  else
    dserror("ERROR: Shape of integration cell not supported!");

  return;
}

/*---------------------------------------------------------------------*
 | output                                                  farah 03/14 |
 *---------------------------------------------------------------------*/
void VOLMORTAR::Cell::Print()
{
  for (int i = 0; i < 4; ++i)
    std::cout << "coords= " << coords_(0, i) << " " << coords_(1, i) << " " << coords_(2, i)
              << std::endl;

  return;
}
