/*----------------------------------------------------------------------*/
/*! \file

\level 1


*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "4C_coupling_volmortar_cell.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_linalg_serialdensematrix.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | constructor                                             farah 01/14 |
 *---------------------------------------------------------------------*/
Core::VolMortar::Cell::Cell(int id, int nvertices, const Core::LinAlg::SerialDenseMatrix& coords,
    const Core::FE::CellType& shape)
    : id_(id), coords_(coords), shape_(shape)
{
  if (shape_ == Core::FE::CellType::tet4)
  {
    /* get the matrix of the coordinates of edges needed to compute the volume,
    ** which is used here as detJ in the quadrature rule.
    ** ("Jacobian matrix") for the quadrarture rule:
    **             [  1    1    1    1  ]
    ** jac_coord = [ x_1  x_2  x_3  x_4 ]
    **             [ y_1  y_2  y_3  y_4 ]
    **             [ z_1  z_2  z_3  z_4 ]
    */

    Core::LinAlg::Matrix<4, 4> jac;
    for (int i = 0; i < 4; i++) jac(0, i) = 1;
    for (int row = 0; row < 3; row++)
      for (int col = 0; col < 4; col++) jac(row + 1, col) = coords_(row, col);

    // volume of the element
    vol_ = jac.Determinant() / 6.0;
    if (vol_ <= 0.0) FOUR_C_THROW("Element volume %10.5e <= 0.0", vol_);
  }
  else
    vol_ = 0.0;

  // std::cout << "SHAPE=     "<< shape_ << std::endl;
  if (shape_ != Core::FE::CellType::tet4) FOUR_C_THROW("wrong shape");

  return;
}

/*---------------------------------------------------------------------*
 | calculate jacobian for hex elements                     farah 04/14 |
 *---------------------------------------------------------------------*/
double Core::VolMortar::Cell::CalcJac(const double* xi)
{
  double jac = 0.0;

  Core::LinAlg::Matrix<3, 8> derivs;
  const double r = xi[0];
  const double s = xi[1];
  const double t = xi[2];

  Core::FE::shape_function_3D_deriv1(derivs, r, s, t, Core::FE::CellType::hex8);


  Core::LinAlg::Matrix<8, 3> xrefe;
  for (int i = 0; i < 8; ++i)
  {
    xrefe(i, 0) = coords_(0, i);
    xrefe(i, 1) = coords_(1, i);
    xrefe(i, 2) = coords_(2, i);
  }

  Core::LinAlg::Matrix<3, 3> invJ;
  invJ.clear();

  invJ.Multiply(derivs, xrefe);
  jac = invJ.Invert();
  if (jac <= 0.0) FOUR_C_THROW("Element Jacobian mapping %10.5e <= 0.0", jac);

  return jac;
}
/*---------------------------------------------------------------------*
 | mapping from parameter space to global space            farah 01/14 |
 *---------------------------------------------------------------------*/
void Core::VolMortar::Cell::LocalToGlobal(double* local, double* global)
{
  if (shape_ == Core::FE::CellType::tet4)
  {
    // check input
    if (!local) FOUR_C_THROW("ERROR: LocalToGlobal called with xi=nullptr");
    if (!global) FOUR_C_THROW("ERROR: LocalToGlobal called with globcoord=nullptr");

    static const int n = 4;
    static const int ndim = 3;

    for (int i = 0; i < ndim; ++i) global[i] = 0.0;

    Core::LinAlg::Matrix<n, 1> val;
    Core::FE::shape_function_3D(val, local[0], local[1], local[2], shape_);

    for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < ndim; ++j)
      {
        // use shape function values for interpolation
        global[j] += val(i) * coords_(j, i);
      }
    }
  }
  else if (shape_ == Core::FE::CellType::hex8)
  {
    // check input
    if (!local) FOUR_C_THROW("ERROR: LocalToGlobal called with xi=nullptr");
    if (!global) FOUR_C_THROW("ERROR: LocalToGlobal called with globcoord=nullptr");

    static const int n = 8;
    static const int ndim = 3;

    for (int i = 0; i < ndim; ++i) global[i] = 0.0;

    Core::LinAlg::Matrix<n, 1> val;
    Core::FE::shape_function_3D(val, local[0], local[1], local[2], shape_);

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
    FOUR_C_THROW("ERROR: Shape of integration cell not supported!");

  return;
}

/*---------------------------------------------------------------------*
 | output                                                  farah 03/14 |
 *---------------------------------------------------------------------*/
void Core::VolMortar::Cell::Print()
{
  for (int i = 0; i < 4; ++i)
    std::cout << "coords= " << coords_(0, i) << " " << coords_(1, i) << " " << coords_(2, i)
              << std::endl;

  return;
}

FOUR_C_NAMESPACE_CLOSE
