/*---------------------------------------------------------------------*/
/*! \file

\brief Implementation of least squares by Sudhakar for Moment-fitting

\level 3


*----------------------------------------------------------------------*/

#include "baci_cut_least_squares.H"
#include <iostream>
#include <cmath>

// solve the rectangular system with linear least squares
Epetra_SerialDenseVector CORE::GEO::CUT::LeastSquares::linear_least_square()
{
  Epetra_SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());
  Epetra_SerialDenseVector rhs(matri_[0].size());
  sqr = get_square_matrix(rhs);
  unknown_.Size(matri_[0].size());

  Epetra_SerialDenseSolver solve_for_GPweights;
  solve_for_GPweights.SetMatrix(sqr);
  solve_for_GPweights.SetVectors(unknown_, rhs);
  solve_for_GPweights.FactorWithEquilibration(true);
  int err2 = solve_for_GPweights.Factor();
  int err = solve_for_GPweights.Solve();
  if ((err != 0) && (err2 != 0))
    dserror(
        "Computation of Gauss weights failed, Ill"
        "conditioned matrix in least square");


  /*  Epetra_SerialDenseMatrix matt(sqr.size(),sqr.size());
    Epetra_SerialDenseVector unn(sqr.size());
    Epetra_SerialDenseVector rrr(sqr.size());
    for(unsigned i=0;i<sqr.size();i++)
    {
      for(unsigned j=0;j<sqr.size();j++)
        matt(j,i) = sqr[i][j];
      unn(i) = 0.0;
      rrr(i) = rhs[i];
    }
    unknown_ = ConjugateGradient(sqr, rhs);*/

  return unknown_;
}

// premultiplying the matrix with its transpose to get the square matrix
// the source terms also get multiplied
Epetra_SerialDenseMatrix CORE::GEO::CUT::LeastSquares::get_square_matrix(
    Epetra_SerialDenseVector &rhs)
{
  Epetra_SerialDenseMatrix sqr(matri_[0].size(), matri_[0].size());

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
