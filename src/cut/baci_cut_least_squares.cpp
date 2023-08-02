/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3


*----------------------------------------------------------------------*/

#include "baci_cut_least_squares.H"

#include <cmath>
#include <iostream>

// solve the rectangular system with linear least squares
CORE::LINALG::SerialDenseVector CORE::GEO::CUT::LeastSquares::linear_least_square()
{
  CORE::LINALG::SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());
  CORE::LINALG::SerialDenseVector rhs(matri_[0].size());
  sqr = get_square_matrix(rhs);
  unknown_.size(matri_[0].size());

  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_GPweights;
  solve_for_GPweights.setMatrix(Teuchos::rcpFromRef(sqr));
  solve_for_GPweights.setVectors(Teuchos::rcpFromRef(unknown_), Teuchos::rcpFromRef(rhs));
  solve_for_GPweights.factorWithEquilibration(true);
  int err2 = solve_for_GPweights.factor();
  int err = solve_for_GPweights.solve();
  if ((err != 0) && (err2 != 0))
    dserror(
        "Computation of Gauss weights failed, Ill"
        "conditioned matrix in least square");

  return unknown_;
}

// premultiplying the matrix with its transpose to get the square matrix
// the source terms also get multiplied
CORE::LINALG::SerialDenseMatrix CORE::GEO::CUT::LeastSquares::get_square_matrix(
    CORE::LINALG::SerialDenseVector &rhs)
{
  CORE::LINALG::SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());

  for (unsigned i = 0; i < matri_[0].size(); i++)
  {
    for (unsigned j = 0; j < matri_[0].size(); j++)
    {
      sqr(j, i) = 0.0;

      // it is sqr(j,i) because the Epetra elements are column ordered first
      for (unsigned k = 0; k < matri_.size(); k++) sqr(j, i) += matri_[k][i] * matri_[k][j];
    }
  }

  for (unsigned i = 0; i < matri_[0].size(); i++)
  {
    rhs(i) = 0.0;
    for (unsigned j = 0; j < matri_.size(); j++) rhs(i) += matri_[j][i] * sourc_(j);
  }

  return sqr;
}
