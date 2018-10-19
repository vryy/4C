/*----------------------------------------------------------------------*/
/*!
\file FADmatrix_utils.H
\brief A collection of helper methods and classes using automatic differentiation with FAD;

\level 3

<pre>
\maintainer Fabian Braeu
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef FADMATRIX_UTILS_H_
#define FADMATRIX_UTILS_H_

#include "../linalg/linalg_utils.H"
#include <Sacado.hpp>
typedef Sacado::Fad::DFad<double> FAD;

namespace LINALG
{
  /// Serial dense matrix with FAD (automatic differentiation data type) numbers and templated
  /// dimensions
  template <unsigned int rows, unsigned int cols>
  class FADMatrix : public TMatrix<FAD, rows, cols>
  {
   public:
    /// Constructor
    explicit FADMatrix(bool setzero = true) : TMatrix<FAD, rows, cols>(setzero) {}

    /// Constructor
    explicit FADMatrix(double* d, bool view = false) : TMatrix<FAD, rows, cols>(d, view) {}

    /// Constructor
    explicit FADMatrix(const double* d, bool view = false) : TMatrix<FAD, rows, cols>(d, view) {}

    /// Constructor
    explicit FADMatrix(Epetra_SerialDenseMatrix& d, bool view = false)
        : TMatrix<FAD, rows, cols>(d, view)
    {
    }

    /// Constructor
    FADMatrix(FADMatrix<rows, cols>& source, bool view) : TMatrix<FAD, rows, cols>(source, view) {}

    /// Copy constructor
    FADMatrix(const FADMatrix<rows, cols>& source) : TMatrix<FAD, rows, cols>(source) {}

    /// = operator
    inline FADMatrix<rows, cols>& operator=(FADMatrix<rows, cols> const& other)
    {
      TMatrix<FAD, rows, cols>::operator=(other);
      return *this;
    }

    /// = operator
    inline FADMatrix<rows, cols>& operator=(TMatrix<FAD, rows, cols> const& other)
    {
      TMatrix<FAD, rows, cols>::operator=(other);
      return *this;
    }

    /// = operator
    inline FADMatrix<rows, cols>& operator=(TMatrix<double, rows, cols> const& other)
    {
      for (unsigned i = 0; i < rows; ++i)
        for (unsigned j = 0; j < cols; ++j) (*this)(i, j) = other(i, j);

      return *this;
    };

    /// Set all components of the matrix as independent variable
    inline void diff(const int pos,  ///< appending array of derivatives starts from this position
        const int n)                 ///< total length of derivative array
    {
      if (rows != cols) dserror("diff does only work for quadratic matrices");
      for (unsigned i = 0; i < rows; ++i) (*this)(i, i).diff(pos + i, n);
      (*this)(0, 1).diff(pos + 2 + 1, n);
      (*this)(1, 2).diff(pos + 2 + 2, n);
      (*this)(0, 2).diff(pos + 2 + 3, n);
      (*this)(1, 0).diff(pos + 2 + 4, n);
      (*this)(2, 1).diff(pos + 2 + 5, n);
      (*this)(2, 0).diff(pos + 2 + 6, n);
    };

    /// Convert FADMatrix to TMatrix<double,rows,cols>
    inline Matrix<rows, cols> ConverttoDouble() const
    {
      LINALG::Matrix<rows, cols> tmp(true);
      for (unsigned i = 0; i < rows; ++i)
        for (unsigned j = 0; j < cols; ++j) tmp(i, j) = (*this)(i, j).val();

      return tmp;
    };
  };
}  // namespace LINALG


/// Save first derivatives in a 3x3 double matrix
template <typename T>
inline void FirstDerivToMatrix(FAD const& r_fad,  ///< FAD function
    LINALG::TMatrix<T, 3, 3>& out)                ///< First derivatives
{
  for (int i = 0; i < 3; ++i) out(i, i) = r_fad.dx(i);
  out(0, 1) = r_fad.dx(3);
  out(1, 2) = r_fad.dx(4);
  out(0, 2) = r_fad.dx(5);
  out(1, 0) = r_fad.dx(6);
  out(2, 1) = r_fad.dx(7);
  out(2, 0) = r_fad.dx(8);
};

#endif  // FADMATRIX_UTILS_H_
