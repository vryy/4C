/*----------------------------------------------------------------------*/
/*! \file

\brief A class that wraps Teuchos::SerialDenseMatrix

\level 0

*----------------------------------------------------------------------*/

#include "baci_linalg_serialdensematrix.H"


/*----------------------------------------------------------------------*
 |  B = alpha*A + beta*B                                                |
 *----------------------------------------------------------------------*/
void CORE::LINALG::Update(double alpha, const CORE::LINALG::SerialDenseMatrix& A, double beta,
    CORE::LINALG::SerialDenseMatrix& B)
{
  B.scale(beta);
  CORE::LINALG::SerialDenseMatrix Acopy(A);
  Acopy.scale(alpha);
  B += Acopy;
}

/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 | (uses long double to boost accuracy). Do not use for n > 4.          |
 *----------------------------------------------------------------------*/
long double CORE::LINALG::Det_long(const CORE::LINALG::SerialDenseMatrix& A)
{
  if (A.numCols() == 1)
  {
    return A(0, 0);
  }
  else if (A.numCols() == 2)
  {
    long double out_det;
    out_det = ((long double)(A(0, 0)) * (long double)(A(1, 1))) -
              ((long double)(A(0, 1)) * (long double)(A(1, 0)));
    return out_det;
  }
  else if (A.numCols() > 2)
  {
    long double out_det = 0;
    int sign = 1;
    for (int i_col = 0; i_col < A.numCols(); i_col++)
    {
      SerialDenseMatrix temp_matrix(A.numCols() - 1, A.numCols() - 1);
      for (int c_col = 0; c_col < i_col; c_col++)
      {
        for (int row = 1; row < A.numCols(); row++) temp_matrix(row - 1, c_col) = A(row, c_col);
      }
      for (int c_col = i_col + 1; c_col < A.numCols(); c_col++)
      {
        for (int row = 1; row < A.numCols(); row++) temp_matrix(row - 1, c_col - 1) = A(row, c_col);
      }
      out_det = out_det + ((long double)(sign) * (long double)(A(0, i_col)) *
                              (long double)(Det_long(temp_matrix)));
      sign *= -1;
    }
    return out_det;
  }
  else
    return 0;
}

/*----------------------------------------------------------------------*
 |   Set matrix components to zero                                      |
 |   this(0:Length) = 0.0                                               |
 *----------------------------------------------------------------------*/
void CORE::LINALG::Zero(CORE::LINALG::SerialDenseMatrix& mat, int Length)
{
  int cnt = 0;
  for (int j = 0; j < mat.numCols(); ++j)
  {
    for (int i = 0; i < mat.numRows(); ++i)
    {
      if (cnt++ < Length)
        mat(i, j) = 0.0;
      else
        break;
    }
  }
}

/*----------------------------------------------------------------------*
 |   Fill the SerialDenseMatrix column-wise with vector data            |
 *----------------------------------------------------------------------*/
void CORE::LINALG::copy(const double* vec, CORE::LINALG::SerialDenseMatrix::Base& mat)
{
  int cnt = 0;
  for (int j = 0; j < mat.numCols(); ++j)
  {
    for (int i = 0; i < mat.numRows(); ++i)
    {
      mat(i, j) = vec[cnt++];
    }
  }
}
